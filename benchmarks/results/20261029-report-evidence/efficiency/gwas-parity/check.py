"""Seeded public-HRATT parity check against the frozen pre-change source."""
import json
import os
from pathlib import Path
import warnings
import numpy as np
from mixmogam import twostep as ts
from mixmogam.genotypes import Genotypes

out=Path(os.environ['PARITY_OUT']); out.mkdir(parents=True,exist_ok=True)
records=[]
for strong,binary,weighted,packed in [(False,True,False,False),(False,True,True,True),
                                      (False,False,True,True),(True,False,False,False),
                                      (True,True,True,True)]:
    rng=np.random.default_rng(611+strong)
    n,m=180,1024
    f=rng.uniform(.15,.45,m)
    if strong:
        pop=np.repeat(np.arange(3),n//3)
        direction=rng.choice((-1,1),m)
        af=np.clip(f[None,:]+(pop[:,None]-1)*.28*direction[None,:],.02,.98)
    else:
        af=np.broadcast_to(f,(n,m))
    calls=rng.binomial(2,af).astype(np.int8)
    z=(calls-calls.mean(0))/np.maximum(calls.std(0),1e-12)
    genetic=z@rng.normal(size=m)
    y=np.sqrt(.7)*genetic/genetic.std()+rng.normal(0,np.sqrt(.3),n)
    if binary:
        y=(y>np.median(y)).astype(float)
    calls[rng.random(calls.shape)<.01]=-1
    gt=Genotypes(np.asfortranarray(calls),chromosome=np.repeat(np.arange(2),m//2),
                 position=np.arange(m),packed=packed)
    weights=np.exp(rng.normal(0,.3,n)) if weighted else None
    X=rng.normal(size=(n,2))
    options=dict(heritability_method='he',alphas=(-.25,),he_probes=16,random_state=3,
                 grid=[(.1,.3),(.5,.5)],vb_max_iter=300,vb_tol=1e-5,
                 structure_snps=128,n_calibration=8,spa_threshold=.5,
                 block=256,loco_folds=5,trait='binary' if binary else 'quantitative',
                 sample_weights=weights)
    q=3; tri=128*129//2*4
    if not strong:
        options['gram_cache_bytes']=8*(tri+2*128*q*8)+3*(5*tri+6*128*q*8)+1
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        res=ts.hratt(y,gt,X,**options)
    key=f'{int(strong)}-{int(binary)}-{int(weighted)}-{int(packed)}'
    np.savez(out/f'{key}.npz',p=res.p,f_stat=res.f_stat,beta=res.beta,se=res.se,
             af=res.af,cv_mse=res.extra['cv_mse'])
    rec={'case':key,'strong':res.extra['structure']['strong'],'refit':res.extra['loco_refit'],
         'cv_converged':res.extra['cv_converged'],'loco_converged':res.extra['loco_converged'],
         'n_spa':res.extra.get('n_spa'),'cv_best':res.extra['cv_best'],
         'iterations':res.extra['loco_iterations'],'h2':res.extra['h2'],
         'warnings':[str(w.message) for w in caught]}
    records.append(rec); print(rec,flush=True)
(out/'diagnostics.json').write_text(json.dumps(records,indent=2,default=str)+'\n')
