from pathlib import Path
from types import SimpleNamespace
import datetime as dt
import shutil
import sys

root = Path('/Users/au507860/REPOS/mixmogam')
sys.path.insert(0, str(root / 'benchmarks'))
import scale_study as study
from _harness import digest, environment_manifest, power_state, save_json, snapshot, write_csv
import phensim

out = root / 'benchmarks/results/20261029-cache-extension-dev8'
out.mkdir(parents=True, exist_ok=False)
panel = dict(n=5000, m=800000, nchr=20, kind='blocks', block_size=100, seed=20261039)
phenotype = dict(h2=0.5, n_qtl=10, n_background=40000, seed=20261038,
                 confounding=0.0, pcs=0)
power_state()
manifest = {'started_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
            'extension': 'One panel beyond the default cache ladder',
            'panel': panel, 'phenotype': phenotype, 'max_load': 5,
            'driver': str(Path(__file__).resolve()),
            'driver_sha256': digest(Path(__file__)),
            'package_driver_sha256': digest(study.DRIVER),
            'environment': environment_manifest()}
save_json(out / 'manifest.json', manifest)
shutil.copy2(Path(__file__), out / 'run_cache_extension.py')
provenance = {'mixmogam': snapshot(root, out / 'src', extra_files=(study.DRIVER, study.PLAN)),
              'phensim': study.snapshot_package(Path(phensim.__file__).resolve().parent,
                                                out / 'src' / 'phensim')}
save_json(out / 'src' / 'snapshot.json', provenance)
study.warmup(out)
jobs = []
for name, budget in (('default', 1e9), ('large', 2e9)):
    config = study.configuration(name, cache_bytes=budget)
    jobs.append({'arm': 'C2-extension', 'label': f'P3e__{name}',
                 'method': 'hratt-he', 'config': config, 'panel': panel,
                 'panel_name': 'P3e', 'phenotype': phenotype})
rows = study.run_jobs(jobs, out, SimpleNamespace(max_load=5), out / '.jit-cache')
write_csv(out / 'measurements.csv', rows)
comparison = study.agreement_table(out / 'C2-extension',
                                  ['P3e__default', 'P3e__large'], 'P3e__default')
save_json(out / 'agreement.json', comparison)
manifest['finished_utc'] = dt.datetime.now(dt.timezone.utc).isoformat()
manifest['measurements'] = len(rows)
save_json(out / 'manifest.json', manifest)
print('wrote', out, flush=True)
