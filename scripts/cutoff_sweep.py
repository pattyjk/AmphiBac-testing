"""Inhibitory-cutoff sweep: does moving the 80 cutoff (or using 50, the single-dose IC50 analogue) improve leave-study-out prediction?
Consolidated from the interactive session. Inputs: amphibac_2025.1_testbundle.tar.gz (db/), amphibac_window_fragments.tar.gz (frags/V4_254.fasta),
genus_assignments.csv. vsearch 2.32.0."""
import subprocess, numpy as np, pandas as pd
from sklearn.metrics import cohen_kappa_score, balanced_accuracy_score
from scipy.stats import spearmanr
def readfa(f):
    d = {}; k = None
    for l in open(f):
        l = l.strip()
        if l.startswith('>'): k = l[1:].split()[0]; d[k] = []
        elif k: d[k].append(l)
    return {k: ''.join(v) for k, v in d.items()}
au = pd.read_csv('db/audit_per_sequence.csv'); full = readfa('db/FullDatabase_2025.1.fasta'); v4 = readfa('frags/V4_254.fasta')
sc = au[au.score.notna() & au['Assay.Method'].str.contains('CFS', case=False, na=False)].set_index('FastaID'); sc = sc[sc.index.isin(full)]
g = pd.read_csv('genus_assignments.csv').set_index('FastaID')
S = sc[['score', 'Reference']].join(g[['genus_full']]).rename(columns={'Reference': 'study'})
with open('refs.fa', 'w') as h:
    for i in sc.index: h.write(f'>{i}\n{full[i]}\n')
cuts = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95]; rows = []
for L in [120, 150, 254]:
    with open(f'q{L}.fa', 'w') as h:
        for i in sc.index:
            if i in v4 and (len(v4[i][:L]) == L or L == 254): h.write(f'>{i}\n{v4[i][:L]}\n')
    subprocess.run(f'vsearch --usearch_global q{L}.fa --db refs.fa --id 0.99 --strand both --qmask none --minseqlength 1 --maxaccepts 0 --maxrejects 0 '
                   f'--threads 12 --userout h{L}.tsv --userfields query+target+id+alnlen+qlo+qhi+ql+mism+gaps --quiet', shell=True, check=True)
    h = pd.read_csv(f'h{L}.tsv', sep='\t', header=None, names='q t id alnlen qlo qhi ql mism gaps'.split())
    h = h[(h.q != h.t) & (h.qlo == 1) & (h.qhi == h.ql)]                      # full query coverage
    h = h.merge(S[['study', 'score']].add_prefix('q_'), left_on='q', right_index=True).merge(S[['study', 'score']].add_prefix('t_'), left_on='t', right_index=True)
    for scope, hh in [('other studies', h[h.q_study != h.t_study]), ('same study', h[h.q_study == h.t_study])]:
        agg = hh.groupby('q').agg(mean_s=('t_score', 'mean'), scores=('t_score', list)); P = S.loc[agg.index].join(agg)
        if scope == 'other studies':   # genus majority from other studies; common panel
            P = P.dropna(subset=['genus_full'])
            P['g_scores'] = [S[(S.genus_full == r.genus_full) & (S.study != r.study)].score.tolist() for r in P.itertuples()]
            P = P[P.g_scores.str.len() > 0]
        print(L, scope, 'n', len(P), 'rho hits', spearmanr(P.score, P.mean_s)[0])
        for c in cuts:
            y = (P.score > c).values
            pred = {'R4 majority of hits': P.scores.apply(lambda l: np.mean(np.array(l) > c) > 0.5).values}
            if scope == 'other studies':
                pred.update({'R2 any inhibitory hit': P.scores.apply(lambda l: any(v > c for v in l)).values,
                             'Mean score of hits': (P.mean_s > c).values,
                             'Genus majority': P.g_scores.apply(lambda l: np.mean(np.array(l) > c) > 0.5).values})
            for k, p in pred.items():
                rows.append(dict(L=L, scope=scope, cutoff=c, rule=k, prevalence=y.mean(), kappa=cohen_kappa_score(y, p), bal_acc=balanced_accuracy_score(y, p), n=len(y)))
pd.DataFrame(rows).to_csv('cutoff_sweep_kappa.csv', index=False)
# 95% CIs in cutoff_sweep_results.xlsx: 300 study-cluster bootstraps (resample studies with replacement; percentile).
