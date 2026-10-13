"""
identify_control_loops.py -- Step A (Note 33): how does the AiCU climate computer react?
Data only (logger + VIP setpoints), no model. Put this file in GreenhouseSim/test/ and run from GreenhouseSim:
    python3 test/identify_control_loops.py
Estimates, per loop: static characteristic (dead band / proportional band), reaction delay
(lagged cross-correlation of increments) and first-order time constant (ARX fit).
  1. Vents   : U8 = mean(VentLee, Ventwind)/100  vs  e_v = Tair - VentSP (zone mean), daytime
  2. Humidity venting: night vent opening vs humidity deficit relative to its setpoint
  3. Heating : PipeLow (when logged > 5 C) vs e_h = HeatSP - Tair
  4. Closed loop: Tair response to heating-setpoint steps (pre-dawn ramps)
  5. Screens : EnScr / BlackScr open-close events vs radiation, hour, Tout
Needs numpy, pandas, matplotlib.
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "..", "..", "invernaderos2.0", "data"))
FIG  = os.path.join(HERE, "figures", "AiCU_control_loops.png")
DAY0, DT = 43326.0, 300.0

gh = pd.read_csv(os.path.join(DATA, "AiCU", "Greenhouse_climate.csv"))
vp = pd.read_csv(os.path.join(DATA, "AiCU", "vip.csv"))
mt = pd.read_csv(os.path.join(DATA, "meteo.csv"))
n = min(len(gh), len(vp), len(mt)); gh, vp, mt = gh.iloc[:n], vp.iloc[:n], mt.iloc[:n]
d = pd.DataFrame({"t": gh.GHtime})
d["day"] = np.floor(d.t - DAY0); d["hour"] = (d.t - np.floor(d.t)) * 24
d["Tair"] = gh.Tair; d["RH"] = gh.RHair; d["HumDef"] = gh.HumDef
d["U8"] = (gh.VentLee + gh.Ventwind) / 200; d["Lee"] = gh.VentLee / 100; d["Wind"] = gh.Ventwind / 100
d["EnScr"] = gh.EnScr; d["Black"] = gh.BlackScr; d["PipeLow"] = gh.PipeLow; d["PipeGrow"] = gh.PipeGrow
d["VentSP"] = vp[[c for c in vp.columns if c.startswith("VentLeew_Vip") or c.startswith("VentWind_Vip")]].mean(axis=1)
d["HeatSP"] = vp[[f"HeatTemp_Vip_{z}" for z in (1, 2, 3, 5, 6)]].mean(axis=1)
d["HumDefSP"] = vp["HumDef_Vip"]
d["Ig"] = mt.Iglob; d["Tout"] = mt.Tout
d = d[(d.day >= 0) & (d.day <= 114)].reset_index(drop=True)
d["ev"] = d.Tair - d.VentSP; d["eh"] = d.HeatSP - d.Tair
day = d.Ig > 50; night = d.Ig < 1

def pw_fit(e, u):
    """u = clamp((e - e0)/PB, 0, 1): grid least squares"""
    best = (np.inf, None, None)
    for e0 in np.arange(-3, 3.01, 0.05):
        for PB in np.arange(0.2, 12.01, 0.1):
            r = np.sum((np.clip((e - e0) / PB, 0, 1) - u) ** 2)
            if r < best[0]: best = (r, e0, PB)
    return best[1], best[2], 1 - best[0] / np.sum((u - u.mean()) ** 2)

def arx(y, x, valid, dmax=6):
    """y_k = a*y_{k-1} + b*x_{k-d} + c on contiguous valid samples; one fit per delay d"""
    out = []
    for dl in range(0, dmax + 1):
        k = np.arange(max(1, dl), len(y))
        m = valid[k] & valid[k - 1] & valid[k - dl] & ~np.isnan(y[k]) & ~np.isnan(y[k-1]) & ~np.isnan(x[k-dl])
        if m.sum() < 200: continue
        A = np.column_stack([y[k-1][m], x[k-dl][m], np.ones(m.sum())]); b = y[k][m]
        coef, *_ = np.linalg.lstsq(A, b, rcond=None)
        r2 = 1 - np.sum((A @ coef - b) ** 2) / np.sum((b - b.mean()) ** 2)
        out.append((dl, coef, r2, m.sum()))
    return out

def xcorr_lag(x, y, valid, kmax=12):
    """corr(dx_t, dy_{t+k}) for k = 0..kmax on valid increments"""
    dx = np.diff(x); dy = np.diff(y); v = valid[1:] & valid[:-1]
    res = []
    for k in range(0, kmax + 1):
        a = dx[:len(dx)-k]; b = dy[k:]; m = v[:len(v)-k] & v[k:] & ~np.isnan(a) & ~np.isnan(b)
        res.append(np.corrcoef(a[m], b[m])[0, 1] if m.sum() > 100 else np.nan)
    return np.array(res)

print("=" * 96); print("  AiCU control-loop identification (logger data, days 0-114, 5-min samples)"); print("=" * 96)

# ── 1. vents: static characteristic ──────────────────────────────────────────
print("\n[1] VENTS vs Tair - VentSP (daytime, Iglob > 50)")
for lab, rng in [("days 0-41", (0, 41)), ("days 42-114", (42, 114)), ("all", (0, 114))]:
    m = day & d.day.between(*rng) & d.ev.notna() & d.U8.notna()
    e0, PB, r2 = pw_fit(d.ev[m].values, d.U8[m].values)
    e0l, PBl, r2l = pw_fit(d.ev[m].values, d.Lee[m].values)
    e0w, PBw, r2w = pw_fit(d.ev[m].values, d.Wind[m].values)
    print(f"  {lab:12s} U8: start {e0:+5.2f} K, P-band {PB:5.1f} K (R2 {r2:.2f}) | lee: {e0l:+5.2f}/{PBl:4.1f} K | wind: {e0w:+5.2f}/{PBw:4.1f} K  (n={m.sum()})")
m = day & d.ev.notna() & d.U8.notna()
bins = np.arange(-3, 5.01, 0.5); dg = pd.cut(d.ev[m], bins)
tab = d[m].groupby(dg, observed=True).U8.agg(["count", "median", lambda s: s.quantile(.25), lambda s: s.quantile(.75)])
tab.columns = ["n", "median", "q25", "q75"]
print("  binned (all days):"); print("    " + tab.round(3).to_string().replace("\n", "\n    "))

# ── 1b. vents: dynamics ──────────────────────────────────────────────────────
vday = (day & d.ev.notna() & d.U8.notna()).values
xc = xcorr_lag(d.ev.values, d.U8.values, vday)
print("\n  increment cross-correlation corr(dE_t, dU_t+k), k in 5-min steps:")
print("   " + "  ".join(f"k={k}:{c:+.2f}" for k, c in enumerate(xc)))
vunsat = vday & (d.U8 > 0.02).values & (d.U8 < 0.98).values
fits = arx(d.U8.values, d.ev.values, vunsat)
print("  ARX  U_k = a U_k-1 + b e_k-d + c   (unsaturated daytime samples)")
for dl, c, r2, nn in fits:
    tau = -DT / np.log(c[0]) / 60 if 0 < c[0] < 1 else np.nan
    print(f"   d={dl} ({5*dl:2d} min): a={c[0]:.3f} (tau {tau:5.1f} min)  b={c[1]:+.4f} /K  static gain b/(1-a)={c[1]/(1-c[0]):+.3f} /K  R2={r2:.3f}  n={nn}")

# ── 2. humidity venting at night ─────────────────────────────────────────────
print("\n[2] NIGHT vents vs humidity deficit (HumDef - HumDef_SP; negative = too humid) and screen")
m = night & d.HumDef.notna() & d.HumDefSP.notna() & d.U8.notna()
eH = (d.HumDef - d.HumDefSP)[m]
tabh = d[m].groupby(pd.cut(eH, np.arange(-3, 6.01, 0.5)), observed=True).U8.agg(["count", "median", "mean"])
print("    " + tabh.round(3).to_string().replace("\n", "\n    "))
for lab, sm in [("screen < 0.5", d.EnScr < 50), ("screen >= 0.5", d.EnScr >= 50)]:
    mm = m & sm; print(f"  night mean U8 with {lab}: {d.U8[mm].mean():.3f}  (n={mm.sum()})")

# ── 3. heating ───────────────────────────────────────────────────────────────
print("\n[3] HEATING: PipeLow (logged > 5 C) vs e_h = HeatSP - Tair")
pv = (d.PipeLow > 5).values
print(f"  PipeLow valid: {pv.mean()*100:.0f}% of rows; PipeGrow > 5 C: {(d.PipeGrow > 5).mean()*100:.1f}% of rows")
for lab, cond in [("e_h > 0 (heating demand)", d.eh > 0), ("e_h <= 0", d.eh <= 0)]:
    print(f"  {lab:26s}: PipeLow logged on {(pv & cond.values).sum() / max(cond.sum(),1) * 100:5.1f}% of those rows")
m = pv & d.eh.notna().values
dp = (d.PipeLow - d.Tair)[m]
tabp = pd.DataFrame({"eh": d.eh[m], "dp": dp}).groupby(pd.cut(d.eh[m], np.arange(-4, 4.01, 0.5)), observed=True).dp.agg(["count", "median"])
print("  pipe - air [K] by e_h bin:"); print("    " + tabp.round(2).to_string().replace("\n", "\n    "))
fitsp = arx(d.PipeLow.values, d.eh.values, pv)
print("  ARX  Pipe_k = a Pipe_k-1 + b e_h,k-d + c   (logged-pipe samples)")
for dl, c, r2, nn in fitsp:
    tau = -DT / np.log(c[0]) / 60 if 0 < c[0] < 1 else np.nan
    print(f"   d={dl} ({5*dl:2d} min): a={c[0]:.3f} (tau {tau:5.1f} min)  b={c[1]:+.3f} K/K  R2={r2:.3f}  n={nn}")

# ── 4. closed loop: Tair response to heating-setpoint steps ──────────────────
print("\n[4] CLOSED LOOP: Tair after heating-setpoint steps >= +1 K (night/pre-dawn, Iglob < 20)")
sp = d.HeatSP.values; T = d.Tair.values; resp = []; times = []
for k in range(1, len(sp) - 36):
    step = sp[k] - sp[k-1]
    if step >= 1.0 and d.Ig.values[k] < 20 and not np.isnan(T[k-1]):
        seg = T[k-1:k+36] - T[k-1]
        if np.isnan(seg).any(): continue
        target = sp[k] - T[k-1]
        if target <= 0.3: continue
        frac = seg / target; resp.append(frac)
        hit = np.where(frac >= 0.63)[0]; times.append(hit[0] * 5 if len(hit) else np.nan)
times = np.array(times)
if len(times):
    print(f"  {len(times)} steps; time to 63% of the requested rise: median {np.nanmedian(times):.0f} min,"
          f" IQR [{np.nanpercentile(times,25):.0f}, {np.nanpercentile(times,75):.0f}] min, never reached in 3 h: {np.isnan(times).sum()}")
else:
    print("  no qualifying setpoint steps found")

# ── 5. screens ───────────────────────────────────────────────────────────────
print("\n[5] SCREENS: transitions across 50 %")
for col in ("EnScr", "Black"):
    s = d[col].values
    for lab, cond in [("opening (>=50 -> <50)", lambda a, b: a >= 50 and b < 50), ("closing (<50 -> >=50)", lambda a, b: a < 50 and b >= 50)]:
        ix = [k for k in range(1, len(s)) if not np.isnan(s[k-1]) and not np.isnan(s[k]) and cond(s[k-1], s[k])]
        if not ix: print(f"  {col:6s} {lab}: none"); continue
        h = d.hour.values[ix]; ig = d.Ig.values[ix]; to = d.Tout.values[ix]
        print(f"  {col:6s} {lab}: n={len(ix):3d}  hour median {np.nanmedian(h):4.1f} [{np.nanpercentile(h,10):4.1f}-{np.nanpercentile(h,90):4.1f}]"
              f"  Iglob median {np.nanmedian(ig):5.0f} W/m2  Tout median {np.nanmedian(to):4.1f} C")

# ── figure ───────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(2, 2, figsize=(12, 8))
m = day & d.ev.notna() & d.U8.notna()
ax[0,0].scatter(d.ev[m], d.U8[m], s=1, alpha=.15, c="tab:blue")
e0, PB, _ = pw_fit(d.ev[m].values, d.U8[m].values); ee = np.linspace(-3, 5, 200)
ax[0,0].plot(ee, np.clip((ee - e0) / PB, 0, 1), "r", lw=2, label=f"P-band fit: start {e0:+.2f} K, PB {PB:.1f} K")
ax[0,0].set(xlabel="Tair - vent SP [K]", ylabel="vent opening U8", title="Vents (daytime)", xlim=(-3, 5)); ax[0,0].legend()
ax[0,1].plot(np.arange(len(xc)) * 5, xc, "o-"); ax[0,1].set(xlabel="lag [min]", ylabel="corr(dE, dU)", title="Vent reaction delay")
m = pv & d.eh.notna().values
ax[1,0].scatter(d.eh[m], (d.PipeLow - d.Tair)[m], s=1, alpha=.15, c="tab:red")
ax[1,0].set(xlabel="heating SP - Tair [K]", ylabel="PipeLow - Tair [K]", title="Heating pipe (logged)", xlim=(-4, 4))
if resp:
    R = np.array(resp); tt = np.arange(R.shape[1]) * 5
    ax[1,1].plot(tt, R.T, color="gray", alpha=.2); ax[1,1].plot(tt, np.nanmedian(R, 0), "k", lw=3, label="median")
    ax[1,1].axhline(.63, color="r", ls="--"); ax[1,1].set(xlabel="minutes after SP step", ylabel="fraction of requested rise",
        title="Tair response to heating-SP steps", ylim=(-0.5, 1.6)); ax[1,1].legend()
os.makedirs(os.path.dirname(FIG), exist_ok=True)
plt.tight_layout(); plt.savefig(FIG, dpi=110); print("\nFigure:", FIG)
