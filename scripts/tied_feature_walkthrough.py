"""Trace one tied SERDP feature (6a40d31c25c63fc295917d2a6e3c007b) through reference-order permutations.
Consolidated from the interactive session. Inputs: amphibac_2025.1_testbundle.tar.gz (db/), serdp_inputs.tar.gz,
AmphiBac_order_permutation_results.xlsx (sheets serdp_per_feature, serdp_per_sample). vsearch 2.32.0."""
import subprocess, numpy as np, pandas as pd, os, scipy.stats as st
from concurrent.futures import ThreadPoolExecutor
fid = '6a40d31c25c63fc295917d2a6e3c007b'
x = pd.ExcelFile('AmphiBac_order_permutation_results.xlsx')
pf = pd.read_excel(x, 'serdp_per_feature').set_index('feature_id'); seq = pf.loc[fid, 'sequence']
open('q.fa', 'w').write(f'>{fid}\n{seq}\n')
db = {}; k = None
for l in open('db/FullDatabase_2025.1.fasta'):
    l = l.strip()
    if l.startswith('>'): k = l[1:].split()[0]; db[k] = []
    else: db[k].append(l)
db = {k: ''.join(v) for k, v in db.items()}; ids = list(db)
# 1. exhaustive hits >=99%
subprocess.run('vsearch --usearch_global q.fa --db db/FullDatabase_2025.1.fasta --id 0.97 --strand both --qmask none --minseqlength 1 '
               '--maxaccepts 0 --maxrejects 0 --userout ex.tsv --userfields query+target+id+alnlen+mism+gaps+qcov+tcov+ql+tl+qstrand --quiet', shell=True, check=True)
ex = pd.read_csv('ex.tsv', sep='\t', header=None, names='query target id alnlen mism gaps qcov tcov ql tl strand'.split()); ex = ex[ex.id >= 99]
# 2. 1,001 real runs with the q2-vsearch closed-reference settings (seed 0 = distributed order)
def run(seed):
    order = ids if seed == 0 else list(np.random.default_rng(seed).permutation(ids)); f = f'db{seed}.fa'
    with open(f, 'w') as h:
        for i in order: h.write(f'>{i}\n{db[i]}\n')
    o = subprocess.run(f'vsearch --usearch_global q.fa --db {f} --id 0.99 --strand both --qmask none --minseqlength 1 --threads 1 '
                       '--userout - --userfields target+id --quiet', shell=True, capture_output=True, text=True).stdout.split()
    os.remove(f); return seed, (o[0] if o else None)
with ThreadPoolExecutor(10) as e: R = pd.DataFrame(list(e.map(run, range(1001))), columns=['seed', 'target'])
R.to_csv('tied_feature_real_runs.csv', index=False)
# 3. unique shared 8-mers (either strand) per tied reference
def rc(s): return s[::-1].translate(str.maketrans('ACGTN', 'TGCAN'))
def km(s): return {s[i:i+8] for i in range(len(s)-7)}
qk, qrk = km(seq), km(rc(seq))
md = pd.read_csv('db/AmphibBac_Metadata_2025.1.tsv', sep='\t', low_memory=False).set_index('FastaID')
au = pd.read_csv('db/audit_per_sequence.csv').set_index('FastaID')
T = ex.set_index('target')[['id', 'mism', 'tl', 'strand']].join(md[['Host.Species', 'Country', 'State', 'Reference', 'Isolate.ID',
     'Proportional.Growth.Bd', 'AntiFungal.Function.with.AP.function', 'Assay.Method']]).join(au[['score', 'in_Inhibitory', 'length']])
T['kmer_shared'] = [max(len(qk & km(db[h])), len(qrk & km(db[h]))) for h in T.index]
T['wins_of_1001'] = [int((R.target == h).sum()) for h in T.index]; T.to_csv('tied_feature_references.csv')
# 4. per-frog effect and Bd-load slope with this one call flipped (all other features at the default-order call)
ps = pd.read_excel(x, 'serdp_per_sample').set_index('SampleID')
tab = pd.read_csv('serdp_table_rare2900.csv', index_col=0); ps = ps.join((tab.loc[fid] / 2900).rename('feat_frac'))
x_inh = ps.tophit_inh_prop_default_order.values; x_no = x_inh - ps.feat_frac.values   # default order assigns AmphiBac_5440 (Inhibitory)
inf = (ps.Inf == 1).values; ly = ps.logBd.values[inf]
gs = pd.get_dummies(ps.stratum.values[inf]).values.astype(float); H = gs @ np.linalg.pinv(gs); ry = ly - H @ ly; dfree = inf.sum() - gs.shape[1] - 1
def beta_se(xv):
    xm = xv[inf]; z = (xm - xm.mean()) / xm.std(ddof=1); rx = z - H @ z; b = (rx * ry).sum() / (rx ** 2).sum()
    se = np.sqrt(((ry - b * rx) ** 2).sum() / dfree / (rx ** 2).sum()); return b, se, 2 * st.t.sf(abs(b / se), dfree)
print('called Inhibitory', beta_se(x_inh), 'called NoEffect', beta_se(x_no))
C = ps[ps.feat_frac > 0].assign(xno=lambda d: d.tophit_inh_prop_default_order - d.feat_frac)
C[['Species', 'Inf', 'logBd', 'feat_frac', 'xno', 'tophit_inh_prop_default_order', 'tophit_inh_prop_q025', 'tophit_inh_prop_q975']].to_csv('tied_feature_frogs.csv')
# figure fig_tied_feature_walkthrough.png drawn from T, R, C and the two slopes (matplotlib; 4 panels)
