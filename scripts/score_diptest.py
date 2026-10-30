"""Hartigan dip test on AmphiBac Bd inhibition score ((1 - Proportional.Growth.Bd) x 100, CFS assays).
Inputs: tables/audit_per_sequence.csv (score, Func, Reference) and tables/genus_assignments.csv (genus_full).
Requires: pip install diptest scipy statsmodels pandas numpy"""
import pandas as pd, numpy as np, diptest
from scipy.stats import gaussian_kde
from statsmodels.stats.multitest import multipletests
A = pd.read_csv('tables/audit_per_sequence.csv', index_col=0)
A['genus'] = pd.read_csv('tables/genus_assignments.csv').set_index('FastaID').genus_full.reindex(A.index)
s = A[A.score.notna()]
rng = np.random.default_rng(1)
def dip(x):   # +/-0.05 jitter breaks ties (scores reported to ~0.1); 2,000 bootstrap reps
    x = np.asarray(x, float) + rng.uniform(-0.05, 0.05, len(x)); return diptest.diptest(x, boot_pval=True, n_boot=2000)
print('pooled', dip(s.score), 'trimmed -150..150', dip(s.score[s.score.between(-150, 150)]))
rows = []
for grp, col, minn in [('study', 'Reference', 30), ('genus', 'genus', 40)]:
    for k, g in s.groupby(col):
        if len(g) >= minn: d, p = dip(g.score); rows.append(dict(grouping=grp, group=k, n=len(g), dip=d, p=p))
O = pd.DataFrame(rows); O['q_BH'] = O.groupby('grouping').p.transform(lambda p: multipletests(p, method='fdr_bh')[1])
print(O.sort_values('p').head(10))
