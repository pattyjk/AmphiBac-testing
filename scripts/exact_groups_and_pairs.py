"""Exact-sequence aggregation and exact-match pairs in AmphiBac 2025.1 (summary section 8j).
Consolidated from the interactive session code; not re-run end to end as a single file.
Inputs: amphibac_2025.1_testbundle.tar.gz (db/), amphibac_window_fragments.tar.gz (frags/V4_254.fasta, frags/V3V4_429.fasta), genus_assignments.csv."""
import re, numpy as np, pandas as pd
from itertools import combinations
def readfa(f):
    d = {}; k = None
    for l in open(f):
        l = l.strip()
        if l.startswith('>'): k = l[1:].split()[0]; d[k] = []
        elif k: d[k].append(l)
    return {k: ''.join(v).upper() for k, v in d.items()}
IUPAC = {'A':'A','C':'C','G':'G','T':'T','R':'[AG]','Y':'[CT]','S':'[GC]','W':'[AT]','K':'[GT]','M':'[AC]','B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]','N':'[ACGT]'}
rx = lambda p: ''.join(IUPAC[c] for c in p)
rc = lambda s: s[::-1].translate(str.maketrans('ACGTRYKMSWBDHVN', 'TGCAYRMKSWVHDBN'))
# V3 = between 341F (CCTACGGGNGGCWGCAG) and 518R (ATTACCGCGGCTGCTGG) sites, either orientation, insert 100-250 nt
pf = re.compile(rx('CCTACGGGNGGCWGCAG')); pr = re.compile(rx(rc('ATTACCGCGGCTGCTGG')))
full = readfa('db/FullDatabase_2025.1.fasta'); v3 = {}
for k, s in full.items():
    for q in (s, rc(s)):
        a = pf.search(q)
        if a:
            b = pr.search(q, a.end())
            if b and 100 <= b.start() - a.end() <= 250: v3[k] = q[a.end():b.start()]; break
md = pd.read_csv('db/AmphibBac_Metadata_2025.1.tsv', sep='\t', low_memory=False).set_index('FastaID')
au = pd.read_csv('db/audit_per_sequence.csv').set_index('FastaID')
M = md[['Isolate.ID','Reference','Host.Species','Country','Latitude','Longitude','AntiFungal.Function.with.AP.function']].rename(columns={'AntiFungal.Function.with.AP.function':'Bd_label'}).join(au[['score']])
M['host'] = M['Host.Species'].astype(str).str.replace('.', ' ', regex=False).str.strip()
lat = pd.to_numeric(M.Latitude, errors='coerce'); lon = pd.to_numeric(M.Longitude, errors='coerce')
ok = lat.between(-90, 90) & lon.between(-180, 180) & ~((lat == 0) & (lon == 0)); M['lat'] = lat.where(ok); M['lon'] = lon.where(ok)
def hav(a1, o1, a2, o2):
    a1, o1, a2, o2 = map(np.radians, (a1, o1, a2, o2))
    return 2 * 6371.0088 * np.arcsin(np.sqrt(np.sin((a2-a1)/2)**2 + np.cos(a1)*np.cos(a2)*np.sin((o2-o1)/2)**2))
T = ['Inhibitory', 'NoEffect', 'Facilitating']
regions = {'V3': v3, 'V4': readfa('frags/V4_254.fasta'), 'V3V4': readfa('frags/V3V4_429.fasta')}
for r, s in regions.items():
    df = pd.DataFrame({'FastaID': list(s), 'seq': list(s.values())}).set_index('FastaID').join(M)
    vc = df.seq.value_counts(); df['group'] = df.seq.map({q: f'{r}_g{i+1:05d}' for i, q in enumerate(vc.index)})
    # group consensus: single isolate / unanimous / majority (>50% of tested) / conflict (no majority); see AmphiBac_exact_V4_V3V4_groups.xlsx
    if r == 'V3V4': continue
    rows = [(g, a, b) for g, d in df.groupby('group') if len(d) > 1 for a, b in combinations(d.index, 2)]
    P = pd.DataFrame(rows, columns=['group', 'FastaID_1', 'FastaID_2'])
    for i in (1, 2):
        P = P.join(df.loc[P[f'FastaID_{i}'], ['Reference','host','lat','lon','Bd_label','score']].reset_index(drop=True).add_suffix(f'_{i}'))
    P['distance_km'] = hav(P.lat_1, P.lon_1, P.lat_2, P.lon_2)
    P['same_host_species'] = P.host_1 == P.host_2; P['same_study'] = P.Reference_1 == P.Reference_2
    P['both_tested'] = P.Bd_label_1.isin(T) & P.Bd_label_2.isin(T)
    P['same_Bd_label'] = np.where(P.both_tested, P.Bd_label_1 == P.Bd_label_2, np.nan)
    P['score_abs_diff'] = (P.score_1 - P.score_2).abs()
    P.to_csv(f'pairs_{r}.csv.gz', index=False)
    # summaries: group-weighted % same label = mean over groups of within-group mean; distance bins 0, <1, 1-10, 10-100, 100-1k, 1k-5k, >5k km;
    # chance baseline = sum of squared label frequencies among tested isolates (44%). Score gap = median over groups of within-group median |diff|.
