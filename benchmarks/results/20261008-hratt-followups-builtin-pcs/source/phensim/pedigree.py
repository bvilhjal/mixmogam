"""Pedigree simulation and pedigree kinship.

Trio-column pedigrees (``ids``, ``father``, ``mother``), generation-
coherent birth times, the additive relationship matrix ``A`` from the
recursive tabular method (Henderson 1976), and an O(n)-storage Mendelian
draw of genetic values with covariance ``A``. Extracted from ltpred's
population-register machinery; the family-liability (LT-FH) layer built
on top of these stays in ltpred.
"""

from __future__ import annotations

import warnings
from collections import OrderedDict
from typing import Optional, Sequence, Tuple, Union

import numpy as np

__all__ = [
    "simulate_pedigree",
    "pedigree_birth_times",
    "kinship_from_pedigree",
    "mendelian_draw",
]


def simulate_pedigree(
    n_founder_pairs: int = 150,
    gens: int = 3,
    remarry: float = 0.10,
    seed: Union[int, np.random.Generator, None] = None,
) -> Tuple[list, list, list]:
    """Simulate a multi-generation pedigree as trio columns.

    Returns ``(ids, father, mother)``. Founders have ``None`` parents;
    ``gens`` generations of descendants follow ``n_founder_pairs``
    unrelated founder couples (``gens + 1`` generations in all). Each
    couple has 2-4 children, and with probability ``remarry`` the father
    has one further child with a new founder partner, which is what
    produces half-siblings. People sharing any recorded parent are never
    paired. Same RNG call sequence as ltpred's ``simulate_pedigree``
    given the same generator, so seeds reproduce the same pedigrees.
    """
    if (isinstance(n_founder_pairs, (bool, np.bool_))
            or not isinstance(n_founder_pairs, (int, np.integer))
            or n_founder_pairs < 1):
        raise ValueError("n_founder_pairs must be a positive integer")
    if (isinstance(gens, (bool, np.bool_))
            or not isinstance(gens, (int, np.integer)) or gens < 0):
        raise ValueError("gens must be a nonnegative integer")
    try:
        remarry = float(remarry)
    except (TypeError, ValueError):
        raise ValueError("remarry must be a finite value in [0, 1]") from None
    if not 0.0 <= remarry <= 1.0:
        raise ValueError("remarry must be a finite value in [0, 1]")
    rng = np.random.default_rng(seed)
    ids, father, mother = [], [], []
    index = {}  # id -> row; keeps parent lookup O(1) across generations

    def add(f, m):
        """Append person ``p<k>`` with father ``f``, mother ``m``; return the id."""
        pid = f"p{len(ids)}"
        index[pid] = len(ids)
        ids.append(pid)
        father.append(f)
        mother.append(m)
        return pid

    couples = [(add(None, None), add(None, None)) for _ in range(n_founder_pairs)]
    prev_children = []
    for g in range(gens):
        if g > 0:
            pool = list(prev_children)
            rng.shuffle(pool)
            couples = []
            i = 0
            while i + 1 < len(pool):
                a, b = pool[i], pool[i + 1]
                i += 2
                # avoid mating recorded siblings (same recorded parent)
                fa, ma = father[index[a]], mother[index[a]]
                fb, mb = father[index[b]], mother[index[b]]
                if fa is not None and (fa in (fb, mb) or ma in (fb, mb)):
                    continue
                couples.append((a, b))
        next_children = []
        for fa, mo in couples:
            for _ in range(int(rng.integers(2, 5))):
                next_children.append(add(fa, mo))
            if rng.uniform() < remarry:  # second union -> half-sibs
                mate = add(None, None)
                next_children.append(add(fa, mate))
        prev_children = next_children
    return ids, father, mother


def pedigree_birth_times(
    ids: Sequence,
    father: Sequence,
    mother: Sequence,
    *,
    base_birth_year: float = 1920.0,
    generation_years: float = 30.0,
) -> np.ndarray:
    """Assign generation-coherent calendar birth times to a pedigree.

    Co-parents are placed in the same generation (a union-find over
    couples) and every child one generation later than its recorded
    parents, so birth times are ``base_birth_year + generation_years *
    generation``. This discrete-generation model can reject valid acyclic
    pedigrees with generation-skipping matings (e.g. uncle and niece).
    Raises ``ValueError`` for an ancestry cycle or incompatible generation
    constraints; kinship and Mendelian draws do not require this placement.
    Deterministic for a given pedigree (no randomness). Parent values
    that are not listed ids warn once per call (``UserWarning``) and are
    treated as unknown founders.
    """
    father, mother = list(father), list(mother)
    ids, pos, sire, dam, parent_children, unresolved = _parent_links(ids, father, mother)
    _warn_unresolved(unresolved)
    n = len(ids)
    if not np.isfinite(base_birth_year) or not np.isfinite(generation_years) or generation_years <= 0:
        raise ValueError("birth year must be finite and generation_years positive and finite")

    # Check ancestry before merging co-parents: a cycle created by that
    # merge can instead be an acyclic pedigree with overlapping generations.
    remaining = [int(sire[i] != -1) + int(dam[i] != -1) for i in range(n)]
    ready = [i for i, count in enumerate(remaining) if count == 0]
    for i in ready:
        for child in parent_children[i]:
            remaining[child] -= 1
            if remaining[child] == 0:
                ready.append(child)
    if len(ready) != n:
        raise ValueError("pedigree has an ancestry cycle (an individual is its own ancestor)")
    placement_error = (
        "pedigree cannot be placed in discrete generations: co-parents must "
        "share a generation and every child must be exactly one generation "
        "later; generation-skipping matings are unsupported"
    )
    representative = list(range(n))

    def find(i):
        """Union-find root of ``i``, halving the path."""
        while representative[i] != i:
            representative[i] = representative[representative[i]]
            i = representative[i]
        return i

    def union(i, j):
        """Merge the co-parent classes of ``i`` and ``j``."""
        left, right = find(i), find(j)
        if left != right:
            representative[right] = left

    for fa, mo in zip(father, mother):
        if fa in pos and mo in pos:
            union(pos[fa], pos[mo])

    component = np.array([find(i) for i in range(n)], dtype=np.intp)
    members = {root: [] for root in set(component.tolist())}
    for i, root in enumerate(component):
        members[root].append(i)
    children = {root: set() for root in members}
    indegree = {root: 0 for root in members}
    for child, (fa, mo) in enumerate(zip(father, mother)):
        child_root = int(component[child])
        for parent_id in (fa, mo):
            if parent_id not in pos:
                continue
            parent_root = int(component[pos[parent_id]])
            if child_root == parent_root:
                raise ValueError(placement_error)
            if child_root not in children[parent_root]:
                children[parent_root].add(child_root)
                indegree[child_root] += 1

    frontier = [root for root, count in indegree.items() if count == 0]
    generation = {root: 0 for root in frontier}
    visited = 0
    while frontier:
        root = frontier.pop()
        visited += 1
        for child_root in children[root]:
            generation[child_root] = max(
                generation.get(child_root, 0), generation[root] + 1)
            indegree[child_root] -= 1
            if indegree[child_root] == 0:
                frontier.append(child_root)
    if visited != len(members):
        raise ValueError(placement_error)
    for parent_root, child_roots in children.items():
        if any(generation[c] != generation[parent_root] + 1 for c in child_roots):
            raise ValueError(placement_error)
    return np.array([
        base_birth_year + generation_years * generation[int(root)]
        for root in component
    ])


# --------------------------------------------------------------------------- #
# Parent indexing and the dense relationship matrix A
# --------------------------------------------------------------------------- #
def _missing_id(p) -> bool:
    """``None``, ``""`` or NaN: never a usable id."""
    if p is None or (isinstance(p, str) and p == ""):
        return True
    return isinstance(p, (float, np.floating)) and bool(np.isnan(p))


def _missing_marker(p) -> bool:
    """A parent value meaning "unknown" (also ``0``/``"0"``) unless listed."""
    return _missing_id(p) or p == 0 or p == "0"


def _parent_links(ids, father, mother):
    """Index a pedigree: ``(ids, index, sire, dam, children, n_unresolved)``.

    Missing-id markers and parents not among ``ids`` map to ``-1``
    (unknown founder); the latter are counted in ``n_unresolved``.
    Raises on a length mismatch, duplicate ids, a missing id (``None``,
    ``""``, NaN), and a person recorded as their own parent. ``0``/``"0"``
    are valid ids, and resolve as parents when listed in ``ids``.
    """
    ids, father, mother = list(ids), list(father), list(mother)
    n = len(ids)
    if not (len(father) == len(mother) == n):
        raise ValueError("ids, father and mother must share length")
    if any(_missing_id(pid) for pid in ids):
        raise ValueError("ids must not contain missing values")
    if len(set(ids)) != n:
        raise ValueError("ids must be unique")
    index = {pid: i for i, pid in enumerate(ids)}
    unresolved = 0

    def _idx(p):
        """Index of parent ``p``, ``-1`` if missing or unlisted (counted)."""
        nonlocal unresolved
        if p in index:
            return index[p]
        if _missing_marker(p):
            return -1
        unresolved += 1
        return -1

    sire = [_idx(p) for p in father]
    dam = [_idx(p) for p in mother]
    children = [[] for _ in range(n)]
    for i in range(n):
        if sire[i] == i or dam[i] == i:
            raise ValueError(f"individual {ids[i]!r} is its own parent")
        if sire[i] != -1:
            children[sire[i]].append(i)
        if dam[i] != -1:
            children[dam[i]].append(i)
    return ids, index, sire, dam, children, unresolved


def _warn_unresolved(unresolved: int) -> None:
    """Warn once that unlisted parents are treated as unknown founders."""
    if unresolved:
        warnings.warn(
            f"{unresolved} unlisted parent reference(s) treated as unknown "
            "founders", UserWarning, stacklevel=3)


def _kinship_A(sire, dam, children=None) -> np.ndarray:
    """Dense ``A`` from integer parent indices (``-1`` for unknown founders).

    Order members by Kahn's algorithm (parents first, O(n + edges)), then
    for each member ``k`` in that order fill row ``k`` against all earlier
    members at once, ``A[k, :k] = (A[s_k, :k] + A[d_k, :k]) / 2``, mirror
    it, and set ``A[k, k] = 1 + A[s_k, d_k] / 2``. An all-zero sentinel
    row stands in for an unknown parent. O(n**2) time and memory; the
    result is exactly symmetric. Raises on a cycle.
    """
    n = len(sire)
    if children is None:
        children = [[] for _ in range(n)]
        for i in range(n):
            if sire[i] != -1:
                children[sire[i]].append(i)
            if dam[i] != -1:
                children[dam[i]].append(i)

    n_parents = [int(sire[i] != -1) + int(dam[i] != -1) for i in range(n)]
    ready = [i for i in range(n) if n_parents[i] == 0]
    order = []
    head = 0
    while head < len(ready):
        i = ready[head]
        head += 1
        order.append(i)
        for child in children[i]:
            n_parents[child] -= 1
            if n_parents[child] == 0:
                ready.append(child)
    if len(order) < n:
        raise ValueError("pedigree has a cycle (an individual is its own ancestor)")

    # Row/column ``n`` is an all-zero stand-in for an unknown parent.
    pos = np.empty(n, dtype=np.intp)
    pos[order] = np.arange(n)
    s = [pos[sire[i]] if sire[i] != -1 else n for i in order]
    d = [pos[dam[i]] if dam[i] != -1 else n for i in order]
    A = np.zeros((n + 1, n + 1), dtype=np.float64)
    for k in range(n):
        row = 0.5 * (A[s[k], :k] + A[d[k], :k])
        A[k, :k] = row
        A[:k, k] = row
        A[k, k] = 1.0 + 0.5 * A[s[k], d[k]]
    return A[np.ix_(pos, pos)]


def kinship_from_pedigree(ids: Sequence, father: Sequence,
                          mother: Sequence) -> np.ndarray:
    """Additive relationship matrix ``A`` (= 2x kinship) from a pedigree.

    ``ids`` is a sequence of unique ids; ``father``/``mother`` the
    same-length parent columns. A parent that is not itself one of
    ``ids`` (``None``, ``0``, ``""``, NaN, or any unlisted value) is an
    unknown founder. Returns the ``(n, n)`` matrix in ``ids`` order:
    ``A[i, i] = 1 + F_i`` (``F_i`` the inbreeding coefficient) and
    ``A[i, j] = 2 x kinship(i, j)`` -- 0.5 for parent-offspring and full
    sibs, 0.25 for grandparent/half-sib, 0.125 for first cousins.
    Handles inbreeding and any pedigree depth; raises on duplicate ids, a
    self-parent, or a cycle. A parent value that is not a missing marker
    and is not among ``ids`` still resolves to an unknown founder, but
    warns once per call (``UserWarning``) since it usually signals a typo.
    """
    ids, _index, sire, dam, children, unresolved = _parent_links(
        ids, father, mother)
    _warn_unresolved(unresolved)
    return _kinship_A(sire, dam, children)


# --------------------------------------------------------------------------- #
# O(n)-storage Mendelian draw of genetic values with covariance A
# --------------------------------------------------------------------------- #
class _SelectedKinship:
    """Selected entries of A = twice kinship; the graph must not change.

    Algorithm K with an explicit continuation stack (no Python
    recursion) and an LRU memo bounded by ``cache_max_entries``. For a
    requested pair ``(i, j)``: swap so ``rank[i] >= rank[j]``; ``A_ii =
    1`` when either parent is unknown else ``1 + A(s_i, d_i)/2``;
    ``A_ij = (A(s_i, j) + A(d_i, j))/2``; unknown parents contribute 0.
    Ranks come from a Kahn topological order of the whole graph fixed at
    construction (which also rejects a cycle). Pairs whose larger rank
    strictly decreases, so evaluation terminates; entries are dyadic
    rationals, so for realistic depths the floats are exact and agree
    bit for bit with the dense :func:`kinship_from_pedigree` fill.
    """

    def __init__(self, sire, dam, children, cache_max_entries=100_000):
        if isinstance(cache_max_entries, (bool, np.bool_)) or not isinstance(
                cache_max_entries, (int, np.integer)):
            raise TypeError("kinship_cache_size must be an integer")
        if cache_max_entries < 0:
            raise ValueError("kinship_cache_size must be nonnegative")
        self.cache_max_entries = int(cache_max_entries)
        self._cache = OrderedDict()
        self._sire, self._dam = sire, dam

        n = len(sire)
        remaining = [int(s != -1) + int(d != -1) for s, d in zip(sire, dam)]
        ready = [i for i, count in enumerate(remaining) if count == 0]
        self._rank = np.empty(n, dtype=np.intp)
        head = 0
        while head < len(ready):
            i = ready[head]
            self._rank[i] = head
            head += 1
            for child in children[i]:
                remaining[child] -= 1
                if remaining[child] == 0:
                    ready.append(child)
        if head != n:
            raise ValueError("pedigree has a cycle (an individual is its own ancestor)")

    def _remember(self, key, value):
        """Insert as most recently used; evict the LRU entry beyond the cap."""
        if self.cache_max_entries:
            self._cache[key] = value
            self._cache.move_to_end(key)
            if len(self._cache) > self.cache_max_entries:
                self._cache.popitem(last=False)

    def _pair(self, i, j):
        """``A_ij``. Frame stages: 0 evaluate, 1 sire term done, 2 combine
        ``(first + value) / 2``, 3 combine ``1 + value / 2``."""
        stack = [(i, j, 0, 0.0)]
        value = 0.0
        while stack:
            i, j, stage, first = stack.pop()
            if stage == 0:
                if i == -1 or j == -1:
                    value = 0.0
                    continue
                if self._rank[i] < self._rank[j]:
                    i, j = j, i
                key = (i, j)
                cached = self._cache.get(key)
                if cached is not None:
                    self._cache.move_to_end(key)
                    value = cached
                    continue
                s, d = self._sire[i], self._dam[i]
                if i == j:
                    if s == -1 or d == -1:
                        value = 1.0
                        self._remember(key, value)
                    else:
                        stack.append((i, j, 3, 0.0))
                        stack.append((s, d, 0, 0.0))
                else:
                    stack.append((i, j, 1, 0.0))
                    stack.append((s, j, 0, 0.0))
            elif stage == 1:
                stack.append((i, j, 2, value))
                stack.append((self._dam[i], j, 0, 0.0))
            else:
                value = 1.0 + 0.5 * value if stage == 3 else 0.5 * (first + value)
                self._remember((i, j), value)
        return value


def mendelian_draw(
    ids: Sequence,
    father: Sequence,
    mother: Sequence,
    innovations: Optional[np.ndarray] = None,
    seed: Union[int, np.random.Generator, None] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Draw genetic values with covariance ``A`` without materializing A.

    Visit people parents-first and set ``a_i = sum_known a_p / 2 +
    sqrt(w_i) * z_i`` with ``z`` the innovations (drawn from ``seed``
    when not supplied) and ``w_i = 1 - sum_known A_pp / 4``. With both
    parents known this is the Mendelian sampling variance
    ``(1 - (F_s + F_d) / 2) / 2``; a founder has ``w_i = 1``, one known
    parent ``(3 - F_p) / 4``. Then ``a ~ N(0, A)`` exactly, including
    inbreeding on the diagonal. Diagonals ``A_ii`` come from the selected
    -pair recursion. Storage is O(n) plus a bounded cache. Returns
    ``(a, diag(A))``. Unlisted parent references warn once per call
    (``UserWarning``) and are treated as unknown founders.
    """
    ids, _index, sire, dam, children, unresolved = _parent_links(
        ids, father, mother)
    _warn_unresolved(unresolved)
    n = len(ids)
    if innovations is None:
        innovations = np.random.default_rng(seed).standard_normal(n)
    innovations = np.asarray(innovations, dtype=float)
    if innovations.shape != (n,) or not np.isfinite(innovations).all():
        raise ValueError("innovations must be a finite vector with one entry per individual")
    kinship = _SelectedKinship(sire, dam, children)
    order = np.empty(n, dtype=np.intp)
    order[kinship._rank] = np.arange(n)
    diagonal = np.empty(n)
    genetic = np.empty(n)
    for i in order:
        diagonal[i] = kinship._pair(int(i), int(i))
        inherited, variance = 0.0, 1.0
        for parent in (sire[i], dam[i]):
            if parent != -1:
                inherited += 0.5 * genetic[parent]
                variance -= 0.25 * diagonal[parent]
        genetic[i] = inherited + np.sqrt(max(0.0, variance)) * innovations[i]
    return genetic, diagonal
