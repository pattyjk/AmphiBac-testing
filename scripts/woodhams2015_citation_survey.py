"""Survey of works citing Woodhams et al. 2015 (Ecology 96:595, doi 10.1890/14-1837.1; OpenAlex W2167232079):
how many used the database predictively (matched their own sequences to it). Consolidated from the interactive session.
Steps: (1) OpenAlex filter=cites:W2167232079 (api_key required); (2) Europe PMC DOI->PMCID lookup + fullTextXML;
(3) remaining DOIs fetched via Unpaywall/S2/PMC/publisher OA (session full-text tool); PDFs parsed with pypdf;
(4) extract citing paragraphs (XML xref to the reference whose text matches 14-1837/'Antifungal isolates database') and
methods sections; (5) LLM classification (Sonnet-class, forced tool schema) into predictive_amplicon / predictive_isolate /
predictive_other / data_use_nonpredictive / background_citation / unclear, with matching method and threshold;
(6) second pass on methods sections; 25 weak calls re-run on whole text; texts <10k chars treated as abstract-only;
(7) drop PeerJ peer-review reports and preprint/article duplicates (normalised title)."""
import os, re, requests, pandas as pd
key = os.environ['OPENALEX_API_KEY']; rows = []; cur = '*'
while cur:
    r = requests.get('https://api.openalex.org/works', params={'filter': 'cites:W2167232079', 'per-page': 200, 'cursor': cur, 'api_key': key}).json()
    for w in r['results']:
        ab = w.get('abstract_inverted_index') or {}
        abs_ = ' '.join(t for t, _ in sorted(((t, p) for t, ps in ab.items() for p in ps), key=lambda z: z[1]))
        ids = w.get('ids', {})
        rows.append(dict(oa=w['id'], doi=(w.get('doi') or '').replace('https://doi.org/', ''), title=w['title'], year=w['publication_year'],
                         type=w['type'], pmcid=(ids.get('pmcid') or '').split('/')[-1], abstract=abs_))
    cur = r['meta'].get('next_cursor')
pd.DataFrame(rows).to_csv('citers_raw.csv', index=False)
def method_family(method, thr):
    s = f'{method} {thr}'.lower()
    if re.search(r'genus|family|taxonom|rdp|sintax|classifier', s): return 'taxonomy-based (genus/family names)'
    if re.search(r'phylogen|character mapping', s): return 'phylogenetic'
    if '100' in s or 'exact' in s: return '100% identity'
    if re.search(r'99', s): return '99% identity'
    if re.search(r'9[78]', s): return '97-98% identity'
    if re.search(r'95|e-value', s): return '<97% / e-value only'
    return 'not stated'
# Classification prompts/schema: see woodhams2015_citers_classified.csv columns (use_category, matching_method, identity_threshold, evidence_quote).
