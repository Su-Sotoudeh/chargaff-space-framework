"""
scoring.py
-----------
Shared, fast scoring for formula and ML long tables.

Long-table schema (both families):
    RowID, UnitCode, Group2, Method, Base, Repeat, Fold, Predicted, Observed, Split
  UnitCode = integer id of (Group2, Organ_ID): the INDEPENDENT unit.
  Repeat   = repeat index of the organism-level CV (formula methods: always 0).

Metrics per (Group2, Method, Base)  [or (Method, Base) when pooled]:
    Pearson_r, CCC (Lin), MAE, RMSE, R2,
    Within_org_r   : r after removing each organism's own mean
                     (can the method predict chromosome-to-chromosome
                     differences inside one organism?)
    Between_org_r  : r on organism means, one point per organism
    + 95% CIs by bootstrapping ORGANISMS (not rows) for Pearson_r, CCC,
      MAE, Within_org_r.
For ML with several CV repeats every metric is computed per repeat and then
averaged (Pearson_r_sd = SD across repeats).

Everything is computed from per-organism sufficient statistics, so the
bootstrap is a single matrix product (fast even for 28,000-row groups).
"""
import numpy as np
import pandas as pd

_K = ["n", "sx", "sy", "sxx", "syy", "sxy", "sabs", "ssq",
      "wc", "xm", "ym", "xm2", "ym2", "xmym", "wxy", "wxx", "wyy"]
_IDX = {k: i for i, k in enumerate(_K)}
_NK = len(_K)
_VAR_EPS = 1e-10
CI_METRICS = ["Pearson_r", "CCC", "MAE", "Within_org_r"]
ALL_METRICS = ["Pearson_r", "CCC", "MAE", "RMSE", "R2", "Within_org_r", "Between_org_r"]


def unit_stats(long, keys):
    """Per (keys, Repeat, UnitCode) sums. Values are shifted by the combo's
    observed mean first (r, CCC, MAE, RMSE, R2 are shift-invariant) to avoid
    floating-point cancellation."""
    keys = list(keys)
    t = long[keys + ["Repeat", "UnitCode"]].copy()
    obs = long["Observed"].to_numpy(float)
    prd = long["Predicted"].to_numpy(float)
    c0 = pd.Series(obs).groupby([long[k].to_numpy() for k in keys]).transform("mean").to_numpy()
    x, y = obs - c0, prd - c0
    t["n"] = 1.0
    t["sx"], t["sy"] = x, y
    t["sxx"], t["syy"], t["sxy"] = x * x, y * y, x * y
    t["sabs"], t["ssq"] = np.abs(x - y), (x - y) ** 2
    g = t.groupby(keys + ["Repeat", "UnitCode"], observed=True, sort=True).sum().reset_index()
    n = g["n"]
    g["wc"] = 1.0
    g["xm"], g["ym"] = g["sx"] / n, g["sy"] / n
    g["xm2"], g["ym2"], g["xmym"] = g["xm"] ** 2, g["ym"] ** 2, g["xm"] * g["ym"]
    g["wxy"] = g["sxy"] - g["sx"] * g["sy"] / n
    g["wxx"] = g["sxx"] - g["sx"] ** 2 / n
    g["wyy"] = g["syy"] - g["sy"] ** 2 / n
    return g


def _metrics_from_sums(S):
    c = lambda k: S[..., _IDX[k]]
    n = c("n")
    mx, my = c("sx") / n, c("sy") / n
    vx, vy = c("sxx") / n - mx ** 2, c("syy") / n - my ** 2
    cov = c("sxy") / n - mx * my
    with np.errstate(all="ignore"):
        ok = (vx > _VAR_EPS) & (vy > _VAR_EPS)
        r = np.where(ok, cov / np.sqrt(vx * vy), np.nan)
        ccc = np.where(ok, 2 * cov / (vx + vy + (mx - my) ** 2), np.nan)
        mae = c("sabs") / n
        rmse = np.sqrt(c("ssq") / n)
        sst = c("sxx") - c("sx") ** 2 / n
        r2 = np.where(sst > _VAR_EPS * n, 1 - c("ssq") / sst, np.nan)
        wxx, wyy = c("wxx"), c("wyy")
        wr = np.where((wxx > _VAR_EPS) & (wyy > _VAR_EPS), c("wxy") / np.sqrt(wxx * wyy), np.nan)
        M = c("wc")
        mxm, mym = c("xm") / M, c("ym") / M
        vxm, vym = c("xm2") / M - mxm ** 2, c("ym2") / M - mym ** 2
        cvm = c("xmym") / M - mxm * mym
        br = np.where((vxm > _VAR_EPS) & (vym > _VAR_EPS), cvm / np.sqrt(vxm * vym), np.nan)
    out = dict(Pearson_r=r, CCC=ccc, MAE=mae, RMSE=rmse, R2=r2, Within_org_r=wr, Between_org_r=br)
    for k in ("Pearson_r", "CCC", "Within_org_r", "Between_org_r"):
        out[k] = np.clip(out[k], -1, 1)
    return out


def _tensor(sub):
    reps = np.unique(sub["Repeat"].to_numpy())
    units = np.unique(sub["UnitCode"].to_numpy())
    ri = np.searchsorted(reps, sub["Repeat"].to_numpy())
    ui = np.searchsorted(units, sub["UnitCode"].to_numpy())
    T = np.zeros((len(reps), len(units), _NK))
    T[ri, ui, :] = sub[_K].to_numpy(float)
    return T


def _nanmean(a, axis):
    with np.errstate(all="ignore"):
        valid = np.isfinite(a)
        s = np.where(valid, a, 0.0).sum(axis=axis)
        c = valid.sum(axis=axis)
        return np.where(c > 0, s / np.maximum(c, 1), np.nan)


def _bootstrap(T, n_boot, rng, chunk=200):
    R, U, K = T.shape
    T2 = T.transpose(1, 0, 2).reshape(U, R * K)
    p = np.full(U, 1.0 / U)
    acc = {m: [] for m in CI_METRICS}
    done = 0
    while done < n_boot:
        b = min(chunk, n_boot - done)
        C = rng.multinomial(U, p, size=b).astype(float)
        S = (C @ T2).reshape(b, R, K)
        m = _metrics_from_sums(S)
        for name in CI_METRICS:
            acc[name].append(_nanmean(m[name], axis=1))
        done += b
    out = {}
    for name in CI_METRICS:
        v = np.concatenate(acc[name])
        v = v[np.isfinite(v)]
        out[name] = tuple(np.percentile(v, [2.5, 97.5])) if len(v) >= 20 else (np.nan, np.nan)
    return out


def score_from_stats(stats, keys, n_boot=1000, seed=0, min_n=3, min_units_ci=10):
    keys = list(keys)
    rng = np.random.default_rng(seed)
    rows = []
    for key, sub in stats.groupby(keys, observed=True, sort=True):
        if not isinstance(key, tuple):
            key = (key,)
        T = _tensor(sub)
        R, U, _ = T.shape
        S = T.sum(axis=1)                     # (R, K)
        N = S[:, _IDX["n"]].mean()
        if N < min_n:
            continue
        m = _metrics_from_sums(S)             # each (R,)
        row = dict(zip(keys, key))
        for name in ALL_METRICS:
            row[name] = float(_nanmean(m[name], axis=0))
        r = m["Pearson_r"]
        row["Pearson_r_sd"] = float(np.nanstd(r, ddof=1)) if (R > 1 and np.isfinite(r).sum() > 1) else np.nan
        row["N"] = int(round(N))
        row["N_orgs"] = U
        if n_boot > 0 and U >= min_units_ci:
            ci = _bootstrap(T, n_boot, rng)
        else:
            ci = {k: (np.nan, np.nan) for k in CI_METRICS}
        for k in CI_METRICS:
            row[f"{k}_lo"], row[f"{k}_hi"] = ci[k]
        rows.append(row)
    cols = keys + ALL_METRICS[:1] + ["Pearson_r_sd", "Pearson_r_lo", "Pearson_r_hi",
            "CCC", "CCC_lo", "CCC_hi", "MAE", "MAE_lo", "MAE_hi", "RMSE", "R2",
            "Within_org_r", "Within_org_r_lo", "Within_org_r_hi", "Between_org_r", "N", "N_orgs"]
    return pd.DataFrame(rows)[cols] if rows else pd.DataFrame(columns=cols)


def score_per_group(long, n_boot=1000, seed=0, min_samples=3):
    st = unit_stats(long, ["Group2", "Method", "Base"])
    return score_from_stats(st, ["Group2", "Method", "Base"], n_boot=n_boot, seed=seed, min_n=min_samples)


def score_pooled(long, n_boot=200, seed=0):
    st = unit_stats(long, ["Method", "Base"])
    return score_from_stats(st, ["Method", "Base"], n_boot=n_boot, seed=seed, min_n=2)


def best_per_group_dict(scores, value_col="Pearson_r"):
    best = {}
    for (g, base), sub in scores.dropna(subset=[value_col]).groupby(["Group2", "Base"], observed=True):
        best[(g, base)] = sub.sort_values(value_col, ascending=False).iloc[0]["Method"]
    return best


def best_pooled_dict(scores, value_col="Pearson_r"):
    best = {}
    for base, sub in scores.groupby("Base", observed=True):
        s = sub.dropna(subset=[value_col]).sort_values(value_col, ascending=False)
        if len(s):
            best[base] = s.iloc[0]["Method"]
    return best


def add_baseline_skill(scores, baseline_scores, on):
    """Adds MAE_skill_vs_M0 = 1 - MAE/MAE_M0 and Pearson_gain_vs_M0 = r - r_M0."""
    b = baseline_scores[list(on) + ["MAE", "Pearson_r"]].rename(
        columns={"MAE": "MAE_M0", "Pearson_r": "Pearson_r_M0"})
    out = scores.merge(b, on=list(on), how="left")
    out["MAE_skill_vs_M0"] = 1 - out["MAE"] / out["MAE_M0"]
    out["Pearson_gain_vs_M0"] = out["Pearson_r"] - out["Pearson_r_M0"]
    return out.drop(columns=["MAE_M0", "Pearson_r_M0"])
