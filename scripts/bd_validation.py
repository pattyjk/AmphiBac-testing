"""Round 2 - Bd validation of AmphiBac community metrics in the bullfrog nutrient-enrichment experiment.
Single-threaded (set OMP/OPENBLAS threads=1 before running)."""
import os, sys, warnings
for v in ["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS"]: os.environ[v]="1"
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")
MAP, PS, ASV, OUT = sys.argv[1:5]
os.makedirs(OUT, exist_ok=True)

m = pd.read_csv(MAP, sep="\t")
ps = pd.read_csv(PS)
with open(ASV) as f: first = f.readline()
asv = pd.read_csv(ASV, sep="\t", skiprows=1 if first.startswith("#") else 0, index_col=0)

# ---------------- metadata audit ----------------
notes = {}
notes["n_map"] = len(m); notes["frogs_per_pond_map"] = m.Pond.value_counts().sort_index().to_dict()
notes["map_missing_16S"] = sorted(set(m.SampleID) - set(asv.columns))
reads = asv.sum()
notes["16S_low_read"] = {k: int(v) for k, v in reads[reads < 1000].items()}
notes["n_analysed"] = len(ps)
notes["logBd_is_log10"] = bool(np.allclose(np.log10(m.Bd_load), m.LogBd, atol=1e-4))
g = m.groupby("Type").Bd_load.agg(["min", "max", "count"]); notes["Bd_by_type"] = g.to_dict()
notes["nitrate_monotone_in_file_order"] = bool(m.NitrateuM.is_monotonic_increasing)
notes["nitrate_unique_per_pond"] = m.groupby("Pond").NitrateuM.nunique().to_dict()
notes["nitrate_vs_LogBd_spearman"] = float(stats.spearmanr(m.NitrateuM, m.LogBd)[0])
notes["nitrate_sorted_equals_column"] = bool((np.sort(m.NitrateuM.values) == m.NitrateuM.values).all())
# High-N frogs with a Low-N-like nitrate value -> proof of misalignment
notes["highN_rows_with_nitrate<5"] = m.loc[(m.Type == "High N") & (m.NitrateuM < 5), "SampleID"].tolist()
notes["lowN_rows_with_nitrate>=4.8"] = m.loc[(m.Type == "Low N") & (m.NitrateuM >= 4.8), "SampleID"].tolist()

# ---------------- join ----------------
METRICS = {
 "Full_Proportion": "Cultured share (Full DB, 99%)",
 "Full_Richness": "Full DB richness (99%)",
 "Inhibitory_Proportion": "Inhibitory proportion (99%)",
 "Inhibitory@0.97_Proportion": "Inhibitory proportion (97%)",
 "Inhibitory@1.00_Proportion": "Inhibitory proportion (100%)",
 "Inhibitory_Richness": "Inhibitory richness (99%)",
 "Inhibitory_strict_Proportion": "InhibitoryStrict proportion (99%)",
 "Inhibitory_strict_Richness": "InhibitoryStrict richness (99%)",
 "Facilitating_Proportion": "Facilitating proportion (99%)",
 "Facilitating_Richness": "Facilitating richness (99%)",
 "Majority_Proportion": "Majority-vote inhibitory proportion (99%)",
 "Consensus75_Proportion": ">=75% consensus inhibitory proportion (99%)",
 "Exh_AnyInhibitory_Proportion": "Any-hit inhibitory proportion (exhaustive, 99%)",
 "ScoreWeighted_mean": "Score-weighted mean CFS score",
 "Inhibitory_over_Full": "Inhibitory / Full ratio (null-adjusted for cultured share)",
}
d = ps[["SampleID"] + list(METRICS)].merge(m[["SampleID", "Type", "Pond", "LogBd"]], on="SampleID", how="left")
assert d.LogBd.notna().all() and len(d) == 40
d["HighN"] = (d.Type == "High N").astype(int)
d.to_csv(f"{OUT}/bf_metrics_with_Bd.csv", index=False)
notes["analysed_per_pond"] = d.Pond.value_counts().sort_index().to_dict()
notes["analysed_per_type"] = d.Type.value_counts().to_dict()

def z(x): return (x - x.mean()) / x.std(ddof=1) if x.std(ddof=1) > 0 else x * 0

rows = []
pm = d.groupby("Pond").agg({**{k: "mean" for k in METRICS}, "LogBd": "mean", "HighN": "first"}).reset_index()
pm["n_frogs"] = d.groupby("Pond").size().values
for k, lab in METRICS.items():
    r = dict(metric=k, label=lab)
    # (1) naive Spearman n=40
    r["naive_rho"], r["naive_p"] = stats.spearmanr(d[k], d.LogBd); r["naive_n"] = len(d)
    # (2) within-treatment Spearman
    for t, tag in [("Low N", "lowN"), ("High N", "highN")]:
        s = d[d.Type == t]; r[f"{tag}_rho"], r[f"{tag}_p"] = stats.spearmanr(s[k], s.LogBd); r[f"{tag}_n"] = len(s)
    # (2b) within-pond (pond-centred) Spearman
    xc = d[k] - d.groupby("Pond")[k].transform("mean"); yc = d.LogBd - d.groupby("Pond").LogBd.transform("mean")
    r["withinpond_rho"], r["withinpond_p"] = stats.spearmanr(xc, yc)
    # (3) mixed model LogBd ~ z(metric) + Type + (1|Pond), ML fit
    dd = d.assign(x=z(d[k]), fp=z(d.Full_Proportion))
    try:
        f = smf.mixedlm("LogBd ~ x + HighN", dd, groups=dd.Pond).fit(reml=True, method=["lbfgs", "powell"])
        ci = f.conf_int().loc["x"]
        r.update(mm_beta_perSD=f.params["x"], mm_lo=ci[0], mm_hi=ci[1], mm_p=f.pvalues["x"],
                 mm_pond_var=float(f.cov_re.iloc[0, 0]), mm_resid_var=f.scale)
    except Exception as e: r["mm_err"] = str(e)
    # (4) mixed model adding cultured share as covariate (does 'inhibitory' add beyond cultured share?)
    if k != "Full_Proportion":
        try:
            f2 = smf.mixedlm("LogBd ~ x + fp + HighN", dd, groups=dd.Pond).fit(reml=True, method=["lbfgs", "powell"])
            ci = f2.conf_int().loc["x"]
            r.update(mmFP_beta_perSD=f2.params["x"], mmFP_lo=ci[0], mmFP_hi=ci[1], mmFP_p=f2.pvalues["x"])
        except Exception as e: r["mmFP_err"] = str(e)
    # (3b) pond-mean OLS LogBd ~ metric + Type, n=10 (unweighted; pond is the experimental unit)
    pp = pm.assign(x=z(pm[k]), fp=z(pm.Full_Proportion))
    o = smf.ols("LogBd ~ x + HighN", pp).fit(); ci = o.conf_int().loc["x"]
    r.update(pond_beta_perSD=o.params["x"], pond_lo=ci[0], pond_hi=ci[1], pond_p=o.pvalues["x"], pond_n=len(pp), pond_df=o.df_resid,
             pond_partial_r=np.sign(o.tvalues["x"]) * np.sqrt(o.tvalues["x"]**2 / (o.tvalues["x"]**2 + o.df_resid)))
    ow = smf.wls("LogBd ~ x + HighN", pp, weights=pp.n_frogs).fit()
    r.update(pondW_beta_perSD=ow.params["x"], pondW_p=ow.pvalues["x"])
    if k != "Full_Proportion":
        o2 = smf.ols("LogBd ~ x + fp + HighN", pp).fit(); r.update(pondFP_beta_perSD=o2.params["x"], pondFP_p=o2.pvalues["x"], pondFP_df=o2.df_resid)
    # treatment association of the metric itself (is the metric just tracking N treatment?)
    r["metric_vs_type_MWU_p"] = stats.mannwhitneyu(d.loc[d.HighN == 1, k], d.loc[d.HighN == 0, k]).pvalue
    r["metric_vs_type_pond_t_p"] = stats.ttest_ind(pm.loc[pm.HighN == 1, k], pm.loc[pm.HighN == 0, k]).pvalue
    rows.append(r)
res = pd.DataFrame(rows)
# (5) BH within each test family across the 15 metrics
for col in ["naive_p", "lowN_p", "highN_p", "withinpond_p", "mm_p", "mmFP_p", "pond_p", "pondW_p", "pondFP_p", "metric_vs_type_MWU_p"]:
    ok = res[col].notna(); res.loc[ok, col.replace("_p", "_q")] = multipletests(res.loc[ok, col], method="fdr_bh")[1]
res.to_csv(f"{OUT}/bd_validation_tests.csv", index=False)
pm.to_csv(f"{OUT}/bf_pond_means.csv", index=False)

# ---------------- power ----------------
# pond-level partial correlation (n=10, covariate Type -> df=7); exact via noncentral t
def power_partial(rho, n, k_cov, alpha, nsim=4000, seed=7):
    """Monte-Carlo power of the t-test on a regression slope with random normal x, k_cov
    binary/normal nuisance covariates (balanced 0/1 treatment first), bivariate-normal partial correlation rho."""
    rng = np.random.default_rng(seed)
    df = n - 2 - k_cov; tc = stats.t.ppf(1 - alpha / 2, df)
    cov = [np.r_[np.zeros(n // 2), np.ones(n - n // 2)]] + [rng.normal(size=n) for _ in range(k_cov - 1)]
    C = np.c_[np.ones(n), *cov] if k_cov else np.ones((n, 1))
    x = rng.normal(size=(nsim, n)); e = rng.normal(size=(nsim, n))
    y = rho / np.sqrt(1 - rho**2) * x + e
    P = np.eye(n) - C @ np.linalg.pinv(C)          # residualise on intercept+covariates
    xr = x @ P; yr = y @ P
    b = (xr * yr).sum(1) / (xr * xr).sum(1); res = yr - b[:, None] * xr
    se = np.sqrt((res**2).sum(1) / df / (xr * xr).sum(1))
    return float(np.mean(np.abs(b / se) > tc))
def mde(n, k_cov, alpha, target=0.8):
    for g in np.arange(0.05, 0.995, 0.005):
        if power_partial(g, n, k_cov, alpha) >= target: return round(float(g), 3)
    return np.nan
nm = len(METRICS)
pw = []
for lab, n, kc in [("pond means, n=10, +Type", 10, 1), ("pond means, n=10, +Type+cultured share", 10, 2),
                   ("frogs Low N, n=%d (ignores pond clustering)" % notes["analysed_per_type"]["Low N"], notes["analysed_per_type"]["Low N"], 0),
                   ("frogs High N, n=%d (ignores pond clustering)" % notes["analysed_per_type"]["High N"], notes["analysed_per_type"]["High N"], 0),
                   ("frogs all, n=40, +Type (ignores pond clustering)", 40, 1)]:
    pw.append(dict(design=lab, n=n, covariates=kc, mde_r_alpha05=mde(n, kc, 0.05), mde_r_bonf=mde(n, kc, 0.05 / nm),
                   power_at_r0p3=power_partial(0.3, n, kc, 0.05), power_at_r0p5=power_partial(0.5, n, kc, 0.05)))
pw = pd.DataFrame(pw)
# Monte-Carlo check of pond-level analytic power at r=0.8
rng = np.random.default_rng(1); hits = 0; NS = 4000
for _ in range(NS):
    t = np.r_[np.zeros(5), np.ones(5)]; x = rng.normal(size=10); y = 0.8 / np.sqrt(1 - .64) * x + rng.normal(size=10) + 2 * t
    X = np.c_[np.ones(10), x, t]; b, *_ = np.linalg.lstsq(X, y, rcond=None); e = y - X @ b; s2 = e @ e / 7
    se = np.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1]); hits += abs(b[1] / se) > stats.t.ppf(.975, 7)
notes["MC_power_pond_r0.8"] = hits / NS
notes["function_power_pond_r0.8"] = float(power_partial(0.8, 10, 1, 0.05))
# Residual SD of pond-mean LogBd within treatment -> MDE in log10 Bd per SD of metric
resid_sd = smf.ols("LogBd ~ HighN", pm).fit().resid.std(ddof=2)
notes["pond_LogBd_within_treatment_SD"] = float(resid_sd)
pw.to_csv(f"{OUT}/bd_power.csv", index=False)
import json; json.dump(notes, open(f"{OUT}/metadata_audit.json", "w"), indent=1, default=str)
print(json.dumps(notes, default=str, indent=0)[:3000])
print(pw.to_string())
