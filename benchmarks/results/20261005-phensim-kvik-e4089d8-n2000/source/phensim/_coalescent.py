"""A self-contained, Numba-JIT coalescent-with-recombination simulator.

This is a compact, from-scratch stand-in for the piece of ``msprime`` that
:mod:`phensim.genotypes` actually uses: simulate the ancestry of a sample under
the coalescent *with recombination* (Hudson's algorithm) for a single
constant-size population, drop infinite-sites (binary) mutations on the
resulting ancestral recombination graph, and return a genotype/dosage matrix
with recombination-driven LD structure.

Why reimplement it?

* **No C-extension dependency.** ``msprime`` (and ``tskit``/``GSL``) is a
  compiled dependency; this backend is pure Python + Numba, so structured-LD
  simulation works with ``phensim[fast]`` and no coalescent extension.
* **Genotypes directly.** ``msprime`` builds a full tree-sequence table
  collection and then materialises the genotype matrix in a separate pass. Here
  the whole pipeline -- ancestry, mutations, densification -- is JIT-compiled and
  writes dosages in one sweep. Runtime and memory relative to msprime vary
  by workload and environment; sibling benchmark drivers measure them.

The algorithm (all times in generations):

1. **Ancestry (Hudson).** Start with ``2*n`` haploid lineages, each carrying the
   whole segment ``[0, L)`` and one descendant sample. Going back in time, draw
   the next event from competing exponentials:

   - coalescence, total rate ``C(k, 2) / (2*Ne)`` (pick two lineages, merge);
   - recombination, total rate ``r * sum_i links_i`` where ``links_i`` is the
     span between a lineage's leftmost and rightmost ancestral base (pick a
     lineage in proportion to its links via a Fenwick tree, split at a uniform
     breakpoint).

   Each segment tracks the number of samples it is ancestral to (``nsamp``);
   this is additive at coalescence and preserved at recombination, so a segment
   reaching ``nsamp == 2*n`` has found its most-recent common ancestor for that
   interval and is dropped -- which is exactly what makes the process terminate.
   Coalescences record tree-sequence edges ``(parent, child, left, right)`` and
   node times.

2. **Mutations (infinite sites).** For each edge, draw
   ``Poisson(mut_rate * branch_length_generations * span_bp)`` mutations at
   uniform positions; each sits on the branch above ``child``.

3. **Densification.** Sweep marginal trees left-to-right (tskit-style edge
   insertion/removal), and for each mutation add 1 to the dosage of every sample
   below its node. Consecutive haplotypes are paired into diploids.

Tests check linked-segment/recombination-weight invariants and compare selected
site-count, diversity and LD summaries with msprime. These bounded checks do not
establish full statistical equivalence; this implementation also has its own RNG.
"""

from __future__ import annotations

import numpy as np

from phensim._numba import HAVE_NUMBA, _jit


# --------------------------------------------------------------------------- #
# Fenwick (binary-indexed) tree over per-lineage recombination "links".
# Used to sample a recombining lineage in proportion to its link length in
# O(log k) and to keep the running total up to date under swap-remove.
# --------------------------------------------------------------------------- #
@_jit
def _fw_update(fw, i, delta):
    i += 1
    n = fw.shape[0]
    while i < n:
        fw[i] += delta
        i += i & (-i)


@_jit
def _fw_find(fw, value):
    """Smallest 0-based index whose inclusive prefix sum exceeds ``value``."""
    n = fw.shape[0]
    pos = 0
    logn = 0
    j = 1
    while j < n:
        j <<= 1
        logn += 1
    bit = 1 << logn
    while bit > 0:
        nxt = pos + bit
        if nxt < n and fw[nxt] <= value:
            pos = nxt
            value -= fw[nxt]
        bit >>= 1
    return pos  # 0-based lineage slot


@_jit
def _slot_set(fw, slot_link, slot, newval):
    """Set lineage ``slot``'s link length to ``newval``; return the delta so the
    caller can keep the running ``total_links`` in step with the Fenwick tree."""
    delta = newval - slot_link[slot]
    slot_link[slot] = newval
    _fw_update(fw, slot, delta)
    return delta


# --------------------------------------------------------------------------- #
# Hudson coalescent-with-recombination. Operates on pre-allocated buffers so it
# can be JIT-compiled with no dynamic allocation; the Python wrapper below sizes
# the buffers and grows/retries on the (rare) overflow.
# --------------------------------------------------------------------------- #
@_jit
def _hudson(n_samples, L, rec_rate, Ne, seed,
            seg_left, seg_right, seg_node, seg_nsamp, seg_prev, seg_next,
            slot_head, slot_tail, slot_link, fw,
            edge_left, edge_right, edge_parent, edge_child, node_time):
    """Simulate the ARG. Returns ``(status, num_nodes, num_edges)``.

    ``status`` is 0 on success, or 1/2/3 if the segment / edge / node buffer
    overflowed (the wrapper then retries with larger buffers).
    """
    np.random.seed(seed)

    max_seg = seg_left.shape[0]
    max_edge = edge_left.shape[0]
    max_node = node_time.shape[0]
    cap = slot_head.shape[0]

    # Segment free list (stack of unused segment slots).
    free_top = max_seg
    free_stack = np.empty(max_seg, dtype=np.int64)
    for i in range(max_seg):
        free_stack[i] = max_seg - 1 - i

    # Sample nodes 0..n_samples-1 exist at time 0.
    for s in range(n_samples):
        node_time[s] = 0.0
    num_nodes = n_samples
    num_edges = 0

    # Initial lineages: one per sample, full segment [0, L), nsamp = 1.
    num_lineages = n_samples
    total_links = 0.0
    for s in range(n_samples):
        sidx = free_stack[free_top - 1]
        free_top -= 1
        seg_left[sidx] = 0.0
        seg_right[sidx] = L
        seg_node[sidx] = s
        seg_nsamp[sidx] = 1
        seg_prev[sidx] = -1
        seg_next[sidx] = -1
        slot_head[s] = sidx
        slot_tail[s] = sidx
        link = L  # right(tail) - left(head)
        slot_link[s] = link
        _fw_update(fw, s, link)
        total_links += link

    t = 0.0
    inv2Ne = 1.0 / (2.0 * Ne)

    while num_lineages > 1:
        k = num_lineages
        rate_coal = (k * (k - 1) * 0.5) * inv2Ne
        rate_rec = rec_rate * total_links
        rate_tot = rate_coal + rate_rec
        if rate_tot <= 0.0:
            break
        t += np.random.exponential(1.0 / rate_tot)

        if np.random.random() * rate_tot < rate_coal:
            # ---- Coalescence: pick two distinct lineage slots. ----
            i = np.random.randint(0, k)
            j = np.random.randint(0, k - 1)
            if j >= i:
                j += 1
            head_a = slot_head[i]
            head_b = slot_head[j]

            # Remove the two slots (larger index first) via swap-remove: move the
            # current last slot into the vacated one, then shrink.
            hi = i if i > j else j
            lo = j if i > j else i
            last = num_lineages - 1
            if hi != last:
                slot_head[hi] = slot_head[last]
                slot_tail[hi] = slot_tail[last]
                total_links += _slot_set(fw, slot_link, hi, slot_link[last])
            total_links += _slot_set(fw, slot_link, last, 0.0)
            num_lineages -= 1
            last = num_lineages - 1
            if lo != last:
                slot_head[lo] = slot_head[last]
                slot_tail[lo] = slot_tail[last]
                total_links += _slot_set(fw, slot_link, lo, slot_link[last])
            total_links += _slot_set(fw, slot_link, last, 0.0)
            num_lineages -= 1

            # Merge the two lineages into (at most) one parent lineage.
            x = head_a
            y = head_b
            new_node = -1
            zhead = -1
            ztail = -1
            while x != -1 or y != -1:
                alpha = -1
                if x == -1 or y == -1:
                    if x != -1:
                        alpha = x
                        x = seg_next[x]
                    else:
                        alpha = y
                        y = seg_next[y]
                    seg_next[alpha] = -1
                else:
                    if seg_left[x] > seg_left[y]:
                        tmp = x
                        x = y
                        y = tmp
                    if seg_right[x] <= seg_left[y]:
                        alpha = x
                        x = seg_next[x]
                        seg_next[alpha] = -1
                    elif seg_left[x] < seg_left[y]:
                        if free_top <= 0:
                            return 1, num_nodes, num_edges
                        aidx = free_stack[free_top - 1]
                        free_top -= 1
                        seg_left[aidx] = seg_left[x]
                        seg_right[aidx] = seg_left[y]
                        seg_node[aidx] = seg_node[x]
                        seg_nsamp[aidx] = seg_nsamp[x]
                        seg_next[aidx] = -1
                        seg_left[x] = seg_left[y]
                        alpha = aidx
                    else:
                        # Overlap [ov_left, r): the two lineages coalesce here.
                        ov_left = seg_left[x]  # == seg_left[y]
                        r = seg_right[x] if seg_right[x] < seg_right[y] else seg_right[y]
                        if new_node == -1:
                            if num_nodes >= max_node:
                                return 3, num_nodes, num_edges
                            new_node = num_nodes
                            node_time[new_node] = t
                            num_nodes += 1
                        if num_edges + 2 > max_edge:
                            return 2, num_nodes, num_edges
                        edge_parent[num_edges] = new_node
                        edge_child[num_edges] = seg_node[x]
                        edge_left[num_edges] = ov_left
                        edge_right[num_edges] = r
                        num_edges += 1
                        edge_parent[num_edges] = new_node
                        edge_child[num_edges] = seg_node[y]
                        edge_left[num_edges] = ov_left
                        edge_right[num_edges] = r
                        num_edges += 1
                        merged = seg_nsamp[x] + seg_nsamp[y]
                        # Advance/trim x and y past the overlap.
                        if seg_right[x] == r:
                            xn = seg_next[x]
                            free_stack[free_top] = x
                            free_top += 1
                            x = xn
                        else:
                            seg_left[x] = r
                        if seg_right[y] == r:
                            yn = seg_next[y]
                            free_stack[free_top] = y
                            free_top += 1
                            y = yn
                        else:
                            seg_left[y] = r
                        if merged == n_samples:
                            alpha = -1  # MRCA reached for this interval: drop
                        else:
                            if free_top <= 0:
                                return 1, num_nodes, num_edges
                            aidx = free_stack[free_top - 1]
                            free_top -= 1
                            seg_left[aidx] = ov_left
                            seg_right[aidx] = r
                            seg_node[aidx] = new_node
                            seg_nsamp[aidx] = merged
                            seg_next[aidx] = -1
                            alpha = aidx
                if alpha != -1:
                    if zhead == -1:
                        zhead = alpha
                        seg_prev[alpha] = -1
                        seg_next[alpha] = -1
                        ztail = alpha
                    else:
                        seg_next[ztail] = alpha
                        seg_prev[alpha] = ztail
                        seg_next[alpha] = -1
                        ztail = alpha

            if zhead != -1:
                # Add the merged lineage back as a new slot.
                if num_lineages >= cap:
                    return 1, num_nodes, num_edges
                link = seg_right[ztail] - seg_left[zhead]
                slot = num_lineages
                slot_head[slot] = zhead
                slot_tail[slot] = ztail
                total_links += _slot_set(fw, slot_link, slot, link)
                num_lineages += 1

        else:
            # ---- Recombination: pick a lineage in proportion to its links. ----
            u = np.random.random() * total_links
            slot = _fw_find(fw, u)
            head = slot_head[slot]
            # breakpoint b uniform in (left(head), right(tail)); clamp to the
            # representable interior so an extreme draw can never sit on or
            # past an endpoint (an endpoint draw used to corrupt the linked
            # lists through a -1 segment index).
            lo_pos = seg_left[head]
            hi_pos = seg_right[slot_tail[slot]]
            b = lo_pos + np.random.random() * slot_link[slot]
            b = min(max(b, np.nextafter(lo_pos, hi_pos)),
                    np.nextafter(hi_pos, lo_pos))
            if not lo_pos < b < hi_pos:
                continue  # no representable interior point

            # Walk to the split point.
            yseg = head
            while yseg != -1 and seg_right[yseg] <= b:
                yseg = seg_next[yseg]
            # yseg is the first segment with right > b; it exists because the
            # clamp keeps b strictly below right(tail).
            if seg_left[yseg] < b:
                # Split inside segment yseg.
                if free_top <= 0:
                    return 1, num_nodes, num_edges
                nidx = free_stack[free_top - 1]
                free_top -= 1
                seg_left[nidx] = b
                seg_right[nidx] = seg_right[yseg]
                seg_node[nidx] = seg_node[yseg]
                seg_nsamp[nidx] = seg_nsamp[yseg]
                seg_next[nidx] = seg_next[yseg]
                seg_prev[nidx] = -1
                if seg_next[yseg] != -1:
                    seg_prev[seg_next[yseg]] = nidx
                seg_right[yseg] = b
                seg_next[yseg] = -1
                left_tail = yseg
                right_head = nidx
            else:
                # Split in the gap just before yseg.
                left_tail = seg_prev[yseg]
                right_head = yseg
                seg_next[left_tail] = -1
                seg_prev[right_head] = -1

            old_tail = slot_tail[slot]
            # Splitting the tail creates a new rightmost segment. The old
            # segment now ends at b and belongs to the left lineage.
            if old_tail == left_tail:
                old_tail = right_head
            # Left lineage stays in this slot (head unchanged, tail = left_tail).
            new_left_link = seg_right[left_tail] - seg_left[head]
            total_links += _slot_set(fw, slot_link, slot, new_left_link)
            slot_tail[slot] = left_tail

            # Right lineage becomes a new slot.
            if num_lineages >= cap:
                return 1, num_nodes, num_edges
            right_link = seg_right[old_tail] - seg_left[right_head]
            ns = num_lineages
            slot_head[ns] = right_head
            slot_tail[ns] = old_tail
            slot_link[ns] = 0.0
            total_links += _slot_set(fw, slot_link, ns, right_link)
            num_lineages += 1

    return 0, num_nodes, num_edges


# --------------------------------------------------------------------------- #
# Infinite-sites mutations: one Poisson draw per edge.
# --------------------------------------------------------------------------- #
@_jit
def _draw_mutations(edge_child, edge_parent, edge_left, edge_right, node_time,
                    num_edges, mut_rate, seed, mut_pos, mut_node):
    """Fill ``mut_pos`` / ``mut_node`` with mutations. Returns ``(status, nmut)``
    (``status == 1`` if the mutation buffer overflowed)."""
    np.random.seed(seed)
    cap = mut_pos.shape[0]
    nmut = 0
    for e in range(num_edges):
        c = edge_child[e]
        p = edge_parent[e]
        blen = node_time[p] - node_time[c]
        span = edge_right[e] - edge_left[e]
        lam = mut_rate * blen * span
        if lam <= 0.0:
            continue
        cnt = np.random.poisson(lam)
        for _ in range(cnt):
            if nmut >= cap:
                return 1, nmut
            mut_pos[nmut] = edge_left[e] + np.random.random() * span
            mut_node[nmut] = c
            nmut += 1
    return 0, nmut


# --------------------------------------------------------------------------- #
# Densify: sweep marginal trees left-to-right and, for each mutation, add 1 to
# the dosage of every sample below its node. Consecutive haplotypes are paired
# into diploids (sample s -> individual s // 2), matching msprime's layout.
# --------------------------------------------------------------------------- #
@_jit
def _build_dosages(edge_child, edge_parent, edge_left, edge_right, num_edges,
                   ins_order, rem_order, mut_pos, mut_node, nmut, L,
                   n_samples, parent, left_child, right_child, left_sib,
                   right_sib, stack, dosage):
    for i in range(parent.shape[0]):
        parent[i] = -1
        left_child[i] = -1
        right_child[i] = -1
        left_sib[i] = -1
        right_sib[i] = -1

    j = 0   # next edge to insert (ins_order)
    kk = 0  # next edge to remove (rem_order)
    mptr = 0
    x = 0.0
    while x < L:
        # Remove edges whose right endpoint is the current coordinate.
        while kk < num_edges and edge_right[rem_order[kk]] == x:
            e = rem_order[kk]
            c = edge_child[e]
            p = edge_parent[e]
            ls = left_sib[c]
            rs = right_sib[c]
            if ls == -1:
                left_child[p] = rs
            else:
                right_sib[ls] = rs
            if rs == -1:
                right_child[p] = ls
            else:
                left_sib[rs] = ls
            parent[c] = -1
            left_sib[c] = -1
            right_sib[c] = -1
            kk += 1
        # Insert edges whose left endpoint is the current coordinate.
        while j < num_edges and edge_left[ins_order[j]] == x:
            e = ins_order[j]
            c = edge_child[e]
            p = edge_parent[e]
            parent[c] = p
            rc = right_child[p]
            if rc == -1:
                left_child[p] = c
                right_child[p] = c
                left_sib[c] = -1
                right_sib[c] = -1
            else:
                right_sib[rc] = c
                left_sib[c] = rc
                right_sib[c] = -1
                right_child[p] = c
            j += 1
        # Right boundary of the current marginal tree.
        right = L
        if j < num_edges and edge_left[ins_order[j]] < right:
            right = edge_left[ins_order[j]]
        if kk < num_edges and edge_right[rem_order[kk]] < right:
            right = edge_right[rem_order[kk]]
        # Score every mutation falling in [x, right) against this tree.
        while mptr < nmut and mut_pos[mptr] < right:
            node = mut_node[mptr]
            sp = 1
            stack[0] = node
            while sp > 0:
                sp -= 1
                uu = stack[sp]
                if uu < n_samples:
                    dosage[uu // 2, mptr] += 1
                else:
                    v = left_child[uu]
                    while v != -1:
                        stack[sp] = v
                        sp += 1
                        v = right_sib[v]
            mptr += 1
        x = right


# --------------------------------------------------------------------------- #
# Python wrapper: size buffers, run the ARG, mutate, densify. Grows and retries
# on the (rare) buffer overflow, so callers never see partial results.
# --------------------------------------------------------------------------- #
def _call_random_kernel(kernel, *args):
    """Invoke a seeding random kernel without leaking RNG state.

    Under Numba the compiled kernels draw from Numba's private RNG, so the
    call is direct. The pure-Python fallback kernels instead seed and draw
    from NumPy's *global* RNG; save and restore its state around the call
    so the caller's stream is unperturbed, including when the kernel
    raises.
    """
    if HAVE_NUMBA:
        return kernel(*args)
    state = np.random.get_state()
    try:
        return kernel(*args)
    finally:
        np.random.set_state(state)


def simulate_dosages(n, seq_len, *, recomb_rate=1e-8, mut_rate=1e-8, Ne=10000,
                     seed=None):
    """Coalescent-with-recombination diploid dosages via the Numba backend.

    Returns ``(G, pos, af)`` where ``G`` is int8 ``(n, n_sites)`` diploid dosages
    (0/1/2), ``pos`` the physical site positions in bp (ascending), and ``af``
    the derived-allele frequency per site. No MAF filtering is applied here.

    ``seed`` is masked to 31 bits for the Numba backend, so it must be a
    non-negative integer below ``2**31`` for reproducibility (values differing
    only in bit 31 collide). ``seed=None`` here means the fixed seed 0, **not** a
    fresh random draw -- the ``simulate.py`` wrappers resolve ``None`` to a random
    seed before calling this; direct callers who want independent replicates must
    pass distinct seeds themselves.
    """
    if (isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer))
            or n < 1):
        raise ValueError("n must be a positive integer")
    if not np.isfinite(seq_len) or seq_len < 1:
        raise ValueError("seq_len must be finite and at least 1")
    for name, rate in (("recomb_rate", recomb_rate), ("mut_rate", mut_rate)):
        if not np.isfinite(rate) or rate < 0:
            raise ValueError(f"{name} must be finite and nonnegative")
    if not np.isfinite(Ne) or Ne <= 0:
        raise ValueError("Ne must be positive and finite")
    n = int(n)
    L = float(int(seq_len))
    n_samples = 2 * n              # diploid -> 2n haploid lineages
    seed = 0 if seed is None else int(seed) & 0x7FFFFFFF

    # Expected recombination events ~ rho * ln(n_samples); size buffers with a
    # generous margin and grow on overflow.
    rho = 4.0 * Ne * recomb_rate * L
    exp_events = n_samples + rho * (1.0 + np.log(n_samples + 1.0))
    seg_cap = int(8 * n_samples + 12 * exp_events + 256)
    node_cap = int(2 * n_samples + 4 * exp_events + 256)
    edge_cap = int(4 * n_samples + 12 * exp_events + 256)
    slot_cap = 1
    while slot_cap < 2 * n_samples + int(2 * exp_events) + 8:
        slot_cap <<= 1  # Fenwick tree wants a power-of-two-friendly capacity

    for _attempt in range(8):
        seg_left = np.empty(seg_cap, dtype=np.float64)
        seg_right = np.empty(seg_cap, dtype=np.float64)
        seg_node = np.empty(seg_cap, dtype=np.int64)
        seg_nsamp = np.empty(seg_cap, dtype=np.int64)
        seg_prev = np.empty(seg_cap, dtype=np.int64)
        seg_next = np.empty(seg_cap, dtype=np.int64)
        slot_head = np.empty(slot_cap, dtype=np.int64)
        slot_tail = np.empty(slot_cap, dtype=np.int64)
        slot_link = np.zeros(slot_cap, dtype=np.float64)
        fw = np.zeros(slot_cap + 1, dtype=np.float64)
        edge_left = np.empty(edge_cap, dtype=np.float64)
        edge_right = np.empty(edge_cap, dtype=np.float64)
        edge_parent = np.empty(edge_cap, dtype=np.int64)
        edge_child = np.empty(edge_cap, dtype=np.int64)
        node_time = np.empty(node_cap, dtype=np.float64)

        status, num_nodes, num_edges = _call_random_kernel(
            _hudson,
            n_samples, L, float(recomb_rate), float(Ne), seed,
            seg_left, seg_right, seg_node, seg_nsamp, seg_prev, seg_next,
            slot_head, slot_tail, slot_link, fw,
            edge_left, edge_right, edge_parent, edge_child, node_time)
        if status == 0:
            break
        if status == 1:
            seg_cap *= 2
            slot_cap <<= 1
        elif status == 2:
            edge_cap *= 2
        else:  # status == 3
            node_cap *= 2
    else:
        raise RuntimeError("coalescent buffers overflowed after repeated growth")

    edge_child = edge_child[:num_edges]
    edge_parent = edge_parent[:num_edges]
    edge_left = edge_left[:num_edges]
    edge_right = edge_right[:num_edges]
    node_time = node_time[:num_nodes]

    # Mutations: expected count sizes the buffer; grow on overflow.
    total_branch = float(np.sum((node_time[edge_parent] - node_time[edge_child])
                                * (edge_right - edge_left)))
    exp_mut = mut_rate * total_branch
    mut_cap = int(exp_mut * 1.5 + 4096)
    for _attempt in range(8):
        mut_pos = np.empty(mut_cap, dtype=np.float64)
        mut_node = np.empty(mut_cap, dtype=np.int64)
        mstatus, nmut = _call_random_kernel(
            _draw_mutations,
            edge_child, edge_parent, edge_left, edge_right, node_time,
            num_edges, float(mut_rate), (seed * 2654435761) & 0x7FFFFFFF,
            mut_pos, mut_node)
        if mstatus == 0:
            break
        mut_cap *= 2
    else:
        raise RuntimeError("mutation buffer overflowed after repeated growth")

    mut_pos = mut_pos[:nmut]
    mut_node = mut_node[:nmut]
    order = np.argsort(mut_pos, kind="mergesort")
    mut_pos = np.ascontiguousarray(mut_pos[order])
    mut_node = np.ascontiguousarray(mut_node[order])

    # Edge orderings for the incremental marginal-tree sweep: insert by (left,
    # parent-time asc), remove by (right, parent-time desc).
    ptime = node_time[edge_parent]
    ins_order = np.lexsort((ptime, edge_left)).astype(np.int64)
    rem_order = np.lexsort((-ptime, edge_right)).astype(np.int64)

    dosage = np.zeros((n, nmut), dtype=np.int8)
    parent = np.empty(num_nodes, dtype=np.int64)
    left_child = np.empty(num_nodes, dtype=np.int64)
    right_child = np.empty(num_nodes, dtype=np.int64)
    left_sib = np.empty(num_nodes, dtype=np.int64)
    right_sib = np.empty(num_nodes, dtype=np.int64)
    stack = np.empty(num_nodes, dtype=np.int64)

    _build_dosages(edge_child, edge_parent, edge_left, edge_right, num_edges,
                   ins_order, rem_order, mut_pos, mut_node, nmut, L,
                   n_samples, parent, left_child, right_child, left_sib,
                   right_sib, stack, dosage)

    af = dosage.sum(axis=0) / (2.0 * n)
    return dosage, mut_pos, af
