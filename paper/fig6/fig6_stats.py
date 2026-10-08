"""Numbers quoted in the main text for Fig. 6I: detection of layer-fixed variants and the true VAF
detected in half of cases (VAF50), per mean depth. Usage: python3 fig6_stats.py DATA_DIR
Output with plantsoma_obs 1.1.0 (2026-10-08): depth 20: L1 0.377 L3 0.644 VAF50 0.0855; depth 60: L1 0.780
L3 0.930 VAF50 0.032; depth 150: L1 0.961 VAF50 0.0126."""
import sys, glob, numpy as np, pandas as pd
D=sys.argv[1]
c={'L1':0.10,'L2':0.70,'L3':0.20}
for dp in (20,60,150):
    d=pd.concat([pd.read_csv(f) for f in glob.glob(f'{D}/detection_rep*_d{dp}_bg0.csv.gz')])
    fx={k: d[(d.layer==k)&np.isclose(d.true_vaf,0.5*v)].detected.mean() for k,v in c.items()}
    b=np.logspace(-4,0,41); d['bin']=pd.cut(d.true_vaf,b)
    g=d.groupby('bin',observed=True).agg(v=('true_vaf','median'),p=('detected','mean'))
    g=g[g.v<0.3]
    i=np.argmax(g.p.values>=0.5); v50=np.exp(np.interp(0.5,g.p.values[i-1:i+1],np.log(g.v.values[i-1:i+1])))
    print(dp, 'fixed det', {k:round(x,3) for k,x in fx.items()}, 'VAF50', round(v50,4), 'overall', round(d.detected.mean(),3))
