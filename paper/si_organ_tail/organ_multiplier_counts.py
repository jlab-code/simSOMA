"""SI Sec. 7.6 numbers: variants per organ by VAF class for organ_mu_multiplier g = 1 and 5
(topology 01, m=4, rho=0, C=64, P_b=4, P_o=8, O=1e4, 1000 sampled cells, 5 replicates).
Run from the repository root. Output (2026-10-08): g=1: 933.94 low, 62.06 intermediate, 175.17 fixed;
g=5: 4682.26, 74.29, 175.17."""
import sys, numpy as np
sys.path.insert(0,'simSOMA_corefunc')
import topology_io, pipeline_wrapper, self_renewal
topo = topology_io.load_topology_auto(__import__("pathlib").Path("simSOMA_inputs/topologies/topology_01_multibranch_20organs_order1to5.json"),
        mapping={"unit":"years","rate":5.0,"mode":"deterministic"})
for g in (1.0,5.0):
    low=[];mid=[];fix=[]
    for s in range(5):
        sr=self_renewal.SelfRenewalParams(m=4,rho=0.0,mu_unit=1.0,kappa_sr=5.0)
        r=pipeline_wrapper.run_pipeline(topo,sr,sam_boundary_cells=64,branch_precursor_number=4,organ_precursor_number=8,
            organ_total_cells=10000,sequenced_cells=1000,seed=500+s,organ_mu_multiplier=g)
        for oe in r["organ_events"]:
            v=np.array(list(oe["allele_counts_by_mutation"].values()))/1000
            low.append((v<0.05).sum()); mid.append(((v>=0.05)&(v<1)).sum()); fix.append((v>=1).sum())
    print(g, np.mean(low), np.mean(mid), np.mean(fix))
