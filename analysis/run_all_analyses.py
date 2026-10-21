# =====================================================================================
# Vietnamese anterolateral knee anatomy in 3D planning — full analysis (single Colab cell)
#   Part 1  Vietnamese cadaveric study: morphometry, comparison with reference values,
#           sectioning effects, measurement reliability            -> Table 1, Fig. 1
#   Part 2  Population rotation model (20,000 virtual knees), leave-one-study-out
#           validation, procedure outcomes                          -> Fig. 2c, Additional file 1
#   Part 3  Geometry on the template knee (Dryad S193761, CC0): literature checks,
#           Vietnamese origin ring, insertion definitions, fixation angle, tunnels
#                                                                   -> Fig. 2a,b,d, Additional file 1
#   Part 4  Figures (black and white, hatch patterns, no titles inside images)
# Input folder (BASE) must contain:
#   Dataset_DCC_DCTN_1sheet.xlsx   Vietnamese specimen data (available on request)
#   vn_triplicates.csv             three repeated readings per measurement (on request)
#   lit_params_v5.json             published parameters used by the model
#   template_knee_S193761.npz      template bones, landmarks and robot passive path (CC0-derived)
# Outputs are written to BASE/outputs/. Fig. 3 (screenshots) is produced from the web app.
# =====================================================================================
import os, json, numpy as np, pandas as pd, matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
try:
    from google.colab import drive; drive.mount('/content/drive')
    BASE = '/content/drive/MyDrive/Sim_DCC_DCTN/'
except ImportError:
    BASE = os.environ.get('SIM_BASE', './')
OUT = os.path.join(BASE, 'outputs/'); os.makedirs(OUT, exist_ok=True)
import builtins, datetime
_LOG = open(OUT + 'run_log.txt', 'w', encoding='utf-8'); _print = builtins.print
def print(*a, **k):
    _print(*a, **k); _print(*a, **{**k, 'file': _LOG}); _LOG.flush()
print('Run started', datetime.datetime.now().isoformat(timespec='seconds'))
rng = np.random.default_rng(2026)
pd.set_option('display.width', 200)

# ---------------------------------------------------------------- Part 1: Vietnamese cadaveric study
vn = pd.read_excel(BASE + 'Dataset_DCC_DCTN_1sheet.xlsx')
# entries judged to be transcription errors (see Methods, Data cleaning)
vn.loc[vn.sample_id == 1, 'all_d_lfe_mm'] = np.nan      # 79 mm
vn.loc[vn.sample_id == 15, 'all_thick_mm'] = np.nan     # 8.0 mm
vn.loc[vn.sample_id == 19, 'all_length_mm'] = np.nan    # 6.6 mm
vn['ins_ratio'] = vn.all_d_gerdy_mm / (vn.all_d_gerdy_mm + vn.all_d_fibhead_mm)

def describe(x):
    x = pd.Series(x).dropna(); n = len(x); m, s = x.mean(), x.std(ddof=1)
    lo, hi = stats.t.interval(0.95, n - 1, m, s / np.sqrt(n)); return n, m, s, lo, hi

def fmt_p(p): return f"{p:.2g}" if p >= 1e-3 else f"{p:.1e}"

refs = [  # variable, column, reference value, source label
    ('Origin to lateral epicondyle (mm)', 'all_d_lfe_mm', np.hypot(4.7, 2.8), 'Kennedy 2015'),
    ('Insertion to Gerdy tubercle (mm)', 'all_d_gerdy_mm', 21.6, 'Claes 2013'),
    ('Insertion to Gerdy tubercle (mm)', 'all_d_gerdy_mm', 23.9, 'Cai 2017'),
    ('Insertion to fibular head (mm)', 'all_d_fibhead_mm', 23.2, 'Claes 2013'),
    ('Insertion to fibular head (mm)', 'all_d_fibhead_mm', 23.8, 'Cai 2017'),
    ('Insertion ratio', 'ins_ratio', 21.6 / (21.6 + 23.2), 'Claes 2013'),
    ('Insertion ratio', 'ins_ratio', 23.9 / (23.9 + 23.8), 'Cai 2017'),
    ('ALL length (mm)', 'all_length_mm', 39.8, 'Hohmann 2025'),
    ('ALL length (mm)', 'all_length_mm', 44.8, 'Cai 2017'),
    ('Width at joint line (mm)', 'all_w_jl_mm', 6.7, 'Claes 2013'),
    ('Width at joint line (mm)', 'all_w_jl_mm', 5.6, 'Hohmann 2025'),
    ('Thickness (mm)', 'all_thick_mm', 1.6, 'Totlis 2024'),
    ('Thickness (mm)', 'all_thick_mm', 1.3, 'Claes 2013')]
rows = []
for lab, col, ref, src in refs:
    n, m, s, lo, hi = describe(vn[col]); p = stats.ttest_1samp(vn[col].dropna(), ref).pvalue
    rows.append(dict(variable=lab, n=n, mean=m, sd=s, ci_lo=lo, ci_hi=hi, reference=ref, source=src, diff=m - ref, p=p))
for lab, col in [('Width at femur (mm)', 'all_w_fem_mm'), ('Width at tibia (mm)', 'all_w_tib_mm'),
                 ('Femoral footprint (mm2)', 'all_fem_area_mm2'), ('Tibial footprint (mm2)', 'all_tib_area_mm2'),
                 ('ITB width (mm)', 'itb_width_mm'), ('ITB thickness (mm)', 'itb_thick_mm'), ('ITB footprint (mm2)', 'itb_gerdy_area_mm2')]:
    if col in vn:
        n, m, s, lo, hi = describe(vn[col]); rows.append(dict(variable=lab, n=n, mean=m, sd=s, ci_lo=lo, ci_hi=hi))
T1 = pd.DataFrame(rows); T1.round(4).to_csv(OUT + 'Table1_vietnamese_vs_reference.csv', index=False)
prev = stats.binomtest(int(vn.all_present.sum()) if 'all_present' in vn else 22, len(vn), 0.458)
k_all = int(vn.all_present.sum()) if 'all_present' in vn else 22
print('PART 1  Vietnamese cadaveric study')
print(f"  ALL present {k_all}/{len(vn)} (95% CI {prev.proportion_ci(method='exact').low*100:.1f}–{prev.proportion_ci(method='exact').high*100:.1f}%), "
      f"vs pooled Asian 45.8%: p = {fmt_p(prev.pvalue)}")
for c, lab in [('all_lm_rel', 'ALL attached to lateral meniscus'), ('kaplan', 'Kaplan fibres'), ('itb_deep_layer', 'Deep ITB layer'), ('itb_all_rel', 'ITB–ALL fibre continuity (code 2)')]:
    if c in vn: v = vn[c]; k = int((v == (2 if c == 'itb_all_rel' else 1)).sum()); print(f"  {lab}: {k}/{len(vn)}")
print(T1.assign(p=T1.p.map(lambda x: '' if pd.isna(x) else fmt_p(x))).round(3).to_string(index=False))
print(f"  10 mm LET strip uses {1000 / vn.itb_width_mm.mean():.0f}% of the ITB width")

# reliability of the three repeated readings
if os.path.exists(BASE + 'vn_triplicates.csv'):
    tri = pd.read_csv(BASE + 'vn_triplicates.csv')
    def icc31(Y):
        n, k = Y.shape; gm = Y.mean(); msr = k * ((Y.mean(1) - gm) ** 2).sum() / (n - 1)
        mse = ((Y - Y.mean(1, keepdims=True) - Y.mean(0, keepdims=True) + gm) ** 2).sum() / ((n - 1) * (k - 1))
        return (msr - mse) / (msr + (k - 1) * mse), np.sqrt(mse)
    rel = []
    for (test, ang, st), g in tri.groupby(['test', 'angle', 'state']):
        Y = g[['r1', 'r2', 'r3']].to_numpy(float); icc, sem = icc31(Y)
        rel.append(dict(test=test, angle=ang, state=st, n=len(Y), icc=icc, sem=sem, mdc95=1.96 * np.sqrt(2) * sem,
                        identical=np.mean(Y.max(1) == Y.min(1)), whole_numbers=np.mean(Y == np.round(Y))))
    REL = pd.DataFrame(rel); REL.round(3).to_csv(OUT + 'reliability.csv', index=False)
    ir = REL[REL.test == 'ir']; sem_pool = np.sqrt((ir['sem'] ** 2).mean())
    print(f"  Internal rotation reliability: ICC {ir.icc.min():.2f}–{ir.icc.max():.2f}; pooled SEM {sem_pool:.2f}°; "
          f"MDC95 {1.96*np.sqrt(2)*sem_pool:.2f}°; identical triplicates {ir.identical.mean()*100:.0f}%; whole numbers {ir.whole_numbers.mean()*100:.0f}%")
else:
    REL = pd.read_csv(BASE + 'outputs/reliability.csv') if os.path.exists(BASE + 'outputs/reliability.csv') else pd.read_csv(BASE + 'reliability.csv')
    ir = REL[REL.test == 'ir']; sem_pool = np.sqrt((ir['sem'] ** 2).mean()); print(f'  (triplicates not found; reliability read from reliability.csv, pooled SEM {sem_pool:.2f}°)')
MDC = 1.96 * np.sqrt(2) * sem_pool

# sectioning effects (ACL intact)
sect = []
for a in (0, 30, 60, 90):
    x = vn[[f'ir{a}_intact', f'ir{a}_d_itb', f'ir{a}_d_itball']].dropna()
    for lab, v in [('ITB sectioned', x[f'ir{a}_d_itb']), ('ALL added', x[f'ir{a}_d_itball'] - x[f'ir{a}_d_itb']), ('ITB + ALL', x[f'ir{a}_d_itball'])]:
        n, m, s, lo, hi = describe(v)
        sect.append(dict(angle=a, step=lab, n=n, mean_deg=m, ci_lo=lo, ci_hi=hi, p=stats.ttest_1samp(v, 0).pvalue,
                         n_above_MDC=int((v > MDC).sum()), pct_of_intact=(100 * v / x[f'ir{a}_intact']).mean()))
SE = pd.DataFrame(sect); SE.round(4).to_csv(OUT + 'sectioning_effects.csv', index=False)
print(SE.assign(p=SE.p.map(fmt_p)).round(2).to_string(index=False))

# ---------------------------------------------------------------- Part 2: population rotation model
L = json.load(open(BASE + 'lit_params_v5.json', encoding='utf-8')); N = 20000; Tq = 5.0
AS, G, H, PR = L['assumptions'], L['Geeslin2018'], L['Huser2017_ACLintact'], L['procedures']
ang = np.array(G['angles']); ix = {a: i for i, a in enumerate(ang)}
col = lambda S, k: np.array([r[0] for r in G[S][k]]); se_ci = lambda lo, hi, z: (hi - lo) / (2 * z)
thI = (col('seq1', 'Intact') + col('seq2', 'Intact')) / 2; thA = (col('seq1', 'ACL') + col('seq2', 'ACL')) / 2
sdI = (np.array([r[1] for r in G['seq1']['Intact']]) + np.array([r[1] for r in G['seq2']['Intact']])) / 2
thr, over = AS['clinical_thresholds_deg'][1], AS['overconstraint_threshold_deg']; sdd = AS['sd_paired_diff_deg']['value']
Lm = L['Lee2023_meta_ALL_in_ACLdef']
dALL_g = ((col('seq1', 'ACL+ALL') - col('seq1', 'ACL')) + (col('seq2', 'ACL+KAP+ALL') - col('seq2', 'ACL+KAP'))) / 2
dALL = np.where(ang < 30, dALL_g, np.interp(ang, Lm['angles'], [x[0] for x in Lm['dIR_ALL_after_ACL']]))
seALL = np.where(ang < 30, sdd / np.sqrt(20), np.interp(ang, Lm['angles'], [se_ci(x[1], x[2], 1.96) for x in Lm['dIR_ALL_after_ACL']]))
dKAP = ((col('seq1', 'ACL+ALL+KAP') - col('seq1', 'ACL+ALL')) + (col('seq2', 'ACL+KAP') - col('seq2', 'ACL'))) / 2
seKAP = np.full(len(ang), sdd / np.sqrt(20))
dIT = np.interp(ang, H['angles'], [x[0] for x in H['dIR_vs_intact']['ALL+ITB']])
seIT = np.interp(ang, H['angles'], [se_ci(x[1], x[2], 1.645) for x in H['dIR_vs_intact']['ALL+ITB']])
CV = AS['cv_between_knee']
frac = lambda prev, d: thI / prev - thI / (prev + d)
draw = lambda d, se, cv: np.clip(rng.normal(d, se, (N, len(ang))) * rng.normal(1, cv, (N, len(ang))), 0, None)
th0 = np.clip(rng.normal(thI, sdI, (N, len(ang))), 3, None); K0 = Tq / th0
fACL = np.tile(frac(thI, thA - thI), (N, 1))
fAL = frac(thA, draw(dALL, seALL, CV['ALL'])) + frac(thA + dALL, draw(dKAP, seKAP, CV['KAP']))
fIT = frac(thI, draw(dIT, seIT, CV['ITB']))
e = rng.uniform(*AS['ACLR_restore_fraction'], (N, 1)); z = np.clip(rng.normal(1, AS['graft_cv']['value'], (N, 1)), 0.1, None)
dth = lambda f: Tq / np.clip(K0 * (1 - f), 1e-3, None) - th0
hgrid = np.exp(np.linspace(-1.5, 1.5, 61))
curve = {nm: np.array([dth((1 - e) * fACL + h * fa).mean(axis=0) for h in hgrid]) for nm, fa in [('ALLKAP', fAL), ('ALLITB', fIT)]}
base_ACLR = dth((1 - e) * fACL).mean(axis=0)
def h_for(obs, nm, a=None, pooled=False, minus_base=False):
    cv = curve[nm].mean(axis=1) if pooled else curve[nm][:, ix[a]] - (base_ACLR[ix[a]] if minus_base else 0)
    return np.interp(obs, cv, hgrid)
G7, D20, NY = PR['Geeslin2017'], PR['Delaloye2020'], L['Noyes2017_ACLR']   # 'Geeslin2017' = Geeslin et al. AJSM 2018;46(3)
acl_g17 = (np.array(G7['fix30']['ACLR']) + np.array(G7['fix70']['ACLR'])) / 2
studies = {'Geeslin 2018b': np.log([h_for(acl_g17[ix[a]], 'ALLKAP', a) for a in (60, 90)]).mean(),
           'Delaloye 2020': np.log(h_for(D20['ACLR_residual_pooled_deg'][0], 'ALLKAP', pooled=True)),
           'Noyes 2017': np.log([h_for(NY['dIR_ALL+ITB_after_ACLR'][j][0], 'ALLITB', a, minus_base=True) for j, a in enumerate((60, 90))]).mean()}
Hs = pd.Series(studies); mu_h, sd_h = Hs.mean(), Hs.std(ddof=1)
obs_pts = {'Geeslin 2018b': [('ALLKAP', 60, acl_g17[ix[60]], False), ('ALLKAP', 90, acl_g17[ix[90]], False)],
           'Delaloye 2020': [('ALLKAP', None, D20['ACLR_residual_pooled_deg'][0], False)],
           'Noyes 2017': [('ALLITB', 60, NY['dIR_ALL+ITB_after_ACLR'][0][0], True), ('ALLITB', 90, NY['dIR_ALL+ITB_after_ACLR'][1][0], True)]}
loo = []
for s in studies:
    oth = Hs.drop(s); hs = np.exp(rng.normal(oth.mean(), max(oth.std(ddof=1), AS['between_study_sd_floor']['value']), 4000))
    for nm, a, ob, mb in obs_pts[s]:
        cv = curve[nm].mean(axis=1) if a is None else curve[nm][:, ix[a]] - (base_ACLR[ix[a]] if mb else 0)
        pr = np.interp(hs, hgrid, cv); lo, hi = np.percentile(pr, [2.5, 97.5])
        loo.append(dict(study=s, angle='pooled' if a is None else a, observed=ob, predicted=np.median(pr), pi_lo=lo, pi_hi=hi, within=lo <= ob <= hi))
LOO = pd.DataFrame(loo); LOO.round(2).to_csv(OUT + 'AdditionalFile1_leave_one_out.csv', index=False)
def c_from(r, a, i): return (1 / a - 1 / r) / (1 / i - 1 / r)
cA, mL = [], []
for fx in ('fix30', 'fix70'):
    for a in (60, 75, 90):
        j = ix[a]; r = thI[j] + G7[fx]['ACLR'][j]
        ca = c_from(r, thI[j] + G7[fx]['ALLR20'][j], thI[j]); cl = c_from(r, thI[j] + G7[fx]['LET20'][j], thI[j])
        cA.append(('Geeslin', ca)); mL.append(('Geeslin', cl / ca))
for j, a in enumerate((60, 90)):
    o = PR['Noyes2017_ALLR']['diff_vs_ACLR_state'][j][0]; r5 = NY['dIR_ALL+ITB_after_ACLR'][j][0]; ti = thI[ix[a]]
    cA.append(('Noyes', c_from(ti + r5, ti + o, ti)))
cA.append(('Delaloye', D20['c_ALLR_assumed'][0])); mL.append(('Delaloye', D20['m_LET_assumed'][0]))
cdf = pd.DataFrame(cA, columns=['s', 'c']).groupby('s').c.apply(lambda v: np.exp(np.log(v).mean()))
mdf = pd.DataFrame(mL, columns=['s', 'm']).groupby('s').m.apply(lambda v: np.exp(np.log(v).mean()))
mu_c, sd_c = np.log(cdf).mean(), np.log(cdf).std(ddof=1); mu_m, sd_m = np.log(mdf).mean(), max(np.log(mdf).std(ddof=1), 0.05)
h_s = np.exp(rng.normal(mu_h, sd_h, (N, 1))); c_s = np.exp(rng.normal(mu_c, sd_c, (N, 1))) * z; m_s = np.exp(rng.normal(mu_m, sd_m, (N, 1)))
fres = h_s * fAL
scen = {'ACLR': dth((1 - e) * fACL + fres), 'ACLR + ALLR': dth((1 - e) * fACL + fres * (1 - c_s)), 'ACLR + LET': dth((1 - e) * fACL + fres * (1 - c_s * m_s))}
OUTC = pd.DataFrame([dict(procedure=k, angle=a, median=np.median(d[:, j]), p2_5=np.percentile(d[:, j], 2.5), p97_5=np.percentile(d[:, j], 97.5),
                          under=(d[:, j] > thr).mean(), restored=((d[:, j] <= thr) & (d[:, j] >= over)).mean(), over=(d[:, j] < over).mean())
                     for k, d in scen.items() for j, a in enumerate(ang)])
OUTC.round(4).to_csv(OUT + 'procedure_outcomes.csv', index=False)
pd.DataFrame([('between-laboratory factor h (geometric mean)', np.exp(mu_h)), ('SD of log h', sd_h), ('graft effect c, ALLR (geometric mean)', np.exp(mu_c)),
              ('SD of log c', sd_c), ('LET multiplier m', np.exp(mu_m)), ('SD of log m', sd_m), ('virtual knees per procedure', N), ('torque (N·m)', Tq),
              ('under-constrained threshold (°)', thr), ('over-constrained threshold (°)', over)] +
             [(f'h, {k}', np.exp(v)) for k, v in Hs.items()] + [(f'c, {k}', v) for k, v in cdf.items()] + [(f'm, {k}', v) for k, v in mdf.items()],
             columns=['parameter', 'value']).to_csv(OUT + 'model_parameters.csv', index=False)
print('\nPART 2  Population rotation model')
print(f"  between-laboratory factor h = {np.exp(mu_h):.2f} (SD log {sd_h:.2f}); graft effect c = {np.exp(mu_c):.2f} (SD log {sd_c:.2f}); LET multiplier m = {np.exp(mu_m):.2f}")
print(LOO.round(2).to_string(index=False)); print(f"  within 95% PI: {LOO.within.sum()}/{len(LOO)}")
print(OUTC[OUTC.angle.isin([30, 60, 90])].pivot_table(index='procedure', columns='angle', values=['under', 'restored', 'over']).round(3).to_string())

# ---------------------------------------------------------------- Part 3: geometry on the template knee
TK = np.load(BASE + 'template_knee_S193761.npz')
FV, TV = TK['femur'].astype(float), TK['tibia'].astype(float); fT, tT = cKDTree(FV), cKDTree(TV)
LM = {k[3:]: TK[k] for k in TK.files if k.startswith('lm_')}
LFE, MFE, GER, FH = LM['LFE'], LM['MFE'], LM['GERDY'], LM['FH']
PIV = (LFE + MFE) / 2; AXV = (LFE - MFE) / np.linalg.norm(LFE - MFE); GFl = np.linalg.norm(FH - GER); GFu = (FH - GER) / GFl
def kin(f):
    x = min(70, max(0, f)); return -np.interp(x, TK['pflex_ang'], TK['pflex_dIE']), np.interp(x, TK['pflex_ang'], TK['pflex_dAP']), np.interp(x, TK['pflex_ang'], TK['pflex_dML'])
def radial(p): r = np.array([p[0], p[1], 0.0]); n = np.linalg.norm(r); return r / n if n > 1e-9 else np.array([1., 0, 0])
def femP(prox, post): v = FV[fT.query([LFE[0], LFE[1] - post, LFE[2] + prox])[1]]; return v + radial(v)
def tibP(fromG, below): p = GER + GFu * fromG; v = TV[tT.query([p[0], p[1], -below])[1]]; return v + radial(v)
def tib2w(p, f, extra_ir=0.0):
    ir, ap, ml = kin(f); q = Rotation.from_rotvec([0, 0, np.radians(ir + extra_ir)]).apply(p) + np.array([ml, ap, 0])
    return Rotation.from_rotvec(AXV * np.radians(-f)).apply(q - PIV) + PIV
def viaLCL(f): return LM['LCL_F'] + 0.25 * (tib2w(LM['LCL_D'], f) - LM['LCL_F']) + np.array([-3., 0, 0])
def glen(F, T, f, ir=0.0, deep=False):
    t = tib2w(T, f, ir)
    if not deep: return np.linalg.norm(F - t)
    V = viaLCL(f); return np.linalg.norm(F - V) + np.linalg.norm(V - t)
pc = lambda a, b: (a / b - 1) * 100
flex = np.arange(0, 91, 5)
r0, sdr = describe(vn.all_d_lfe_mm)[1], describe(vn.all_d_lfe_mm)[2]
ratio = vn.all_d_gerdy_mm.mean() / (vn.all_d_gerdy_mm.mean() + vn.all_d_fibhead_mm.mean())
Tvn = tibP(round(ratio * GFl, 1), 7)
def rng_mm(F, T=None):
    T = Tvn if T is None else T; v = [glen(F, T, f) for f in flex]; return max(v) - min(v), (max(v) - min(v)) / min(v) * 100
# literature checks (no fitted parameter)
T0 = tibP(round(16 / 31.8 * GFl, 1), 7); Lr = lambda F, f, ir=0, deep=False, T=None: glen(F, T0 if T is None else T, f, ir, deep)
Fepi, Fken, Fpp, Fad = femP(0, 0), femP(4.7, 2.8), femP(8, 8), femP(-4, -4)
h_ = pc(Lr(Fepi, 90), Lr(Fepi, 0)); zi = [pc(Lr(Fken, 90, 25), Lr(Fken, 0)), pc(Lr(Fepi, 90, 25), Lr(Fepi, 0))]
sPP, sAD = pc(Lr(Fpp, 90), Lr(Fpp, 0)), pc(Lr(Fad, 90), Lr(Fad, 0)); d60 = pc(Lr(Fken, 60), Lr(Fken, 0))
Tg, Flet = tibP(0, 20), femP(10, 8)
rg = lambda deep: (lambda v: (max(v) - min(v)) / v[0] * 100)([Lr(Flet, a, 0, deep, Tg) for a in (0, 30, 60, 90)])
VER = pd.DataFrame([
    ['Helito 2014: epicondylar ALL 0–90°', '+16.7%', f'{h_:+.1f}%', abs(h_ - 16.7) < 5],
    ['Zens 2015: per degree of flexion', '+0.15%/°', f'{h_/90:+.2f}%/°', abs(h_ / 90 - 0.15) <= 0.05],
    ['Zens 2015: 90° + 25° internal rotation', '+33.8 ± 9.6%', f'{min(zi):.0f} to {max(zi):.0f}%', min(zi) <= 43.4 and max(zi) >= 24.2],
    ['Imbert 2016: proximal-posterior / epicondylar / anterior-distal', '−, +, +', f'{sPP:+.1f} / {h_:+.1f} / {sAD:+.1f}%', sPP < 0 < h_ and sAD > 0],
    ['Dodds 2014: proximal-posterior 0–60°', 'near-isometric', f'{d60:+.1f}%', abs(d60) <= 3],
    ['Kittl 2015: LET deep vs superficial to LCL', 'deep smaller', f'{rg(True):.1f}% vs {rg(False):.1f}%', rg(True) < rg(False)]],
    columns=['check', 'reference', 'model', 'agrees'])
VER.to_csv(OUT + 'AdditionalFile1_geometry_checks.csv', index=False)
# Vietnamese origin ring
ring = {}
for lab, r in [('mean', r0), ('minusSD', r0 - sdr), ('plusSD', r0 + sdr)]:
    ring[lab] = pd.DataFrame([(phi, *rng_mm(femP(round(r * np.cos(np.radians(phi)), 2), round(r * np.sin(np.radians(phi)), 2))))
                              for phi in range(0, 360, 15)], columns=['direction_deg', 'range_mm', 'range_pct'])
pd.concat({k: v for k, v in ring.items()}).to_csv(OUT + 'ring_sweep.csv')
bphi = int(ring['mean'].loc[ring['mean'].range_mm.idxmin(), 'direction_deg'])
Fvn = femP(round(r0 * np.cos(np.radians(bphi)), 2), round(r0 * np.sin(np.radians(bphi)), 2))
ins = {'Vietnamese ratio': Tvn, 'Vietnamese 16 mm': tibP(16, 7), 'Western 21.6 mm': tibP(21.6, 7)}
curves = {}
for nm, T in ins.items():
    for onm, F in [('Vietnamese origin', Fvn), ('Western origin', Fken)]:
        L0 = glen(F, T, 0); curves[f'{onm} | {nm}'] = [pc(glen(F, T, f), L0) for f in flex]
CUR = pd.DataFrame(curves, index=flex); CUR.round(3).to_csv(OUT + 'insertion_curves.csv')
fixc = {fa: [pc(glen(Fvn, Tvn, f), glen(Fvn, Tvn, fa)) for f in flex] for fa in (0, 30, 60, 90)}
irc = {f: pc(glen(Fvn, Tvn, f, 15), glen(Fvn, Tvn, f)) for f in (30, 60, 90)}
pd.DataFrame(fixc, index=flex).add_prefix('fixed_at_').rename_axis('flexion_deg').round(3).to_csv(OUT + 'fixation_angle_curves.csv')
pd.DataFrame({'flexion_deg': flex, 'length_increase_with_15deg_IR_pct': [pc(glen(Fvn, Tvn, f, 15), glen(Fvn, Tvn, f)) for f in flex]}).round(3).to_csv(OUT + 'internal_rotation_effect.csv', index=False)
# tunnel convergence
def segdist(p1, q1, p2, q2):
    d1, d2, r = q1 - p1, q2 - p2, p1 - p2; a, e_, f_ = d1 @ d1, d2 @ d2, d2 @ r; c, b = d1 @ r, d1 @ d2; den = a * e_ - b * b
    s = np.clip((b * f_ - c * e_) / den, 0, 1) if den > 1e-9 else 0; t = (b * s + f_) / e_
    if t < 0: t, s = 0, np.clip(-c / a, 0, 1)
    elif t > 1: t, s = 1, np.clip((b - c) / a, 0, 1)
    return np.linalg.norm(p1 + d1 * s - (p2 + d2 * t))
tun = []
for x in (0, 8, 15):
    A1 = np.array([LFE[0] + 2, LFE[1] - 0.8 * x, LFE[2] + x])
    for a in range(0, 50, 5):
        for b in range(0, 50, 5):
            ar, br = np.radians(a), np.radians(b); d = np.array([-np.cos(ar) * np.cos(br), np.sin(ar) * np.cos(br), np.sin(br)]); d /= np.linalg.norm(d)
            tun.append(dict(acl_exit_mm=x, aim_anterior=a, aim_proximal=b, cortex_mm=np.linalg.norm(Fvn - A1) - 7.0,
                            deep_mm=segdist(LM['ACL_F'], A1, Fvn + d * 6, Fvn + d * 25) - 7.0))
TUN = pd.DataFrame(tun); TUN.round(2).to_csv(OUT + 'tunnel_convergence.csv', index=False)
print('\nPART 3  Geometry on the template knee')
print(VER.to_string(index=False))
m_ = ring['mean']; print(f"  Vietnamese ring {r0:.1f} mm: best {bphi}° ({m_.range_mm.min():.1f} mm, {m_.range_pct.min():.1f}%), worst "
      f"{int(m_.loc[m_.range_mm.idxmax(), 'direction_deg'])}° ({m_.range_mm.max():.1f} mm, {m_.range_pct.max():.1f}%); "
      f"≤2 mm: {(m_.range_mm <= 2).sum()}/24; ≤4 mm: {(m_.range_mm <= 4).sum()}/24")
for k in ('minusSD', 'plusSD'): print(f"  ring {k}: best {ring[k].range_mm.min():.1f} mm")
print(f"  Western origin: {rng_mm(Fken)[0]:.1f} mm ({rng_mm(Fken)[1]:.1f}%)")
print(f"  Gerdy–fibular head {GFl:.1f} mm; ratio insertion {ratio*GFl:.1f} mm from Gerdy; 16 mm point lies {np.linalg.norm(ins['Vietnamese 16 mm']-Tvn):.1f} mm away")
for c in CUR.columns: print(f"  {c}: range {CUR[c].max()-CUR[c].min():.1f}%")
print('  +15° internal rotation lengthens graft: ' + ', '.join(f'{f}° {v:.1f}%' for f, v in irc.items()))
for fa, v in fixc.items(): print(f"  fixation at {fa}°: max {max(v):+.1f}%, min {min(v):+.1f}%")
t8 = TUN[TUN.acl_exit_mm == 8]; sel = lambda a, b: t8[(t8.aim_anterior == a) & (t8.aim_proximal == b)].iloc[0]
print(f"  tunnels (ACL exit 8 mm): cortex {sel(0,0).cortex_mm:.1f} mm; deep {sel(0,0).deep_mm:.1f} / {sel(30,30).deep_mm:.1f} / {sel(45,45).deep_mm:.1f} mm at 0/30/45°")
for x in (0, 15):
    tx = TUN[TUN.acl_exit_mm == x]; print(f"  ACL exit {x} mm: deep {tx.deep_mm.min():.1f}–{tx.deep_mm.max():.1f} mm, cortex {tx.cortex_mm.iloc[0]:.1f} mm")

# ---------------------------------------------------------------- Part 4: figures (The Knee colour palette)
NAVY, BLUE, STEEL, LIGHT, GOLD, DGOLD, GREY, RED = '#081C55', '#305390', '#6F8BAD', '#99B2C7', '#EEC83C', '#D3AD31', '#D7DEDC', '#C8553D'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.linewidth': 0.8, 'hatch.linewidth': 0.5, 'hatch.color': NAVY, 'axes.edgecolor': NAVY, 'axes.labelcolor': NAVY,
                     'xtick.color': NAVY, 'ytick.color': NAVY, 'text.color': NAVY, 'axes.spines.top': False, 'axes.spines.right': False})
def letter(ax, l): ax.text(-0.14, 1.04, l, transform=ax.transAxes, fontsize=11, fontweight='bold', va='bottom', color=NAVY)
ci = lambda c: (lambda n, m, s, lo, hi: (m, [[m - lo], [hi - m]]))(*describe(vn[c]))
mk = {'Vietnamese (n=22)': ('o', GOLD), 'Claes 2013': ('s', NAVY), 'Cai 2017': ('^', BLUE), 'Zhang 2018': ('D', STEEL),
      'Hohmann 2025': ('v', LIGHT), 'Totlis 2024': ('P', DGOLD)}
def dots(ax, data, sp):
    for i, (m_, dd) in enumerate(data.items()):
        for j, (s, (v, e_)) in enumerate(dd.items()):
            y = i + (j - (len(dd) - 1) / 2) * sp; mm, fc = mk[s]
            ax.errorbar(v, y, xerr=e_, fmt=mm, mfc=fc, mec=NAVY, mew=0.6, ecolor=NAVY, ms=5.5 if s.startswith('Viet') else 4.5, capsize=2, lw=0.8)
    ax.set_yticks(range(len(data))); ax.set_yticklabels(list(data)); ax.invert_yaxis(); ax.set_xlabel('mm'); ax.grid(axis='x', color=GREY, lw=0.5)
fig = plt.figure(figsize=(7.1, 6.6)); gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.35)
ax = fig.add_subplot(gs[0, 0]); letter(ax, 'a')
dots(ax, {'To Gerdy tubercle': {'Vietnamese (n=22)': ci('all_d_gerdy_mm'), 'Claes 2013': (21.6, None), 'Cai 2017': (23.9, 3.3)},
          'To fibular head': {'Vietnamese (n=22)': ci('all_d_fibhead_mm'), 'Claes 2013': (23.2, None), 'Cai 2017': (23.8, 4.0)},
          'ALL length': {'Vietnamese (n=22)': ci('all_length_mm'), 'Cai 2017': (44.8, 5.1), 'Zhang 2018': (35.8, 7.2), 'Hohmann 2025': (39.8, None)}}, 0.2)
ax = fig.add_subplot(gs[0, 1]); letter(ax, 'b')
dots(ax, {'Width at joint line': {'Vietnamese (n=22)': ci('all_w_jl_mm'), 'Claes 2013': (6.7, None), 'Zhang 2018': (6.2, 2.3), 'Hohmann 2025': (5.6, None)},
          'Thickness': {'Vietnamese (n=22)': ci('all_thick_mm'), 'Claes 2013': (1.3, None), 'Zhang 2018': (1.4, 0.9), 'Hohmann 2025': (1.5, None),
                        'Totlis 2024': (1.6, [[0.4], [0.4]])}}, 0.16)
fig.legend(handles=[Line2D([], [], marker=mk[s][0], mfc=mk[s][1], mec=NAVY, ls='none', ms=5, label=s) for s in mk],
           frameon=False, fontsize=6.2, ncol=6, loc='upper center', bbox_to_anchor=(0.5, 0.515))
ax = fig.add_subplot(gs[1, 0]); letter(ax, 'c'); jit = np.random.default_rng(1)
for i, a in enumerate((0, 30, 60, 90)):
    x = vn[[f'ir{a}_d_itb', f'ir{a}_d_itball']].dropna()
    for v, off, m_, fc in [(x[f'ir{a}_d_itb'].to_numpy(), -0.17, 'o', BLUE), ((x[f'ir{a}_d_itball'] - x[f'ir{a}_d_itb']).to_numpy(), 0.17, 's', GOLD)]:
        ax.plot(i + off + jit.uniform(-0.05, 0.05, len(v)), v, m_, mfc=fc, mec=NAVY, mew=0.4, ms=3, lw=0, alpha=0.9)
        n, mu, s_, lo, hi = describe(v); ax.errorbar(i + off + 0.11, mu, yerr=[[mu - lo], [hi - mu]], fmt='_', color=NAVY, ms=9, capsize=3, lw=1.2)
ax.set_xticks(range(4)); ax.set_xticklabels(['0°', '30°', '60°', '90°']); ax.set_xlabel('Knee flexion'); ax.set_ylabel('Increase in internal rotation (°)')
ax.axhline(MDC, color=RED, lw=0.8, ls=(0, (3, 2))); ax.text(-0.45, MDC - 1.1, f'MDC {MDC:.1f}°', fontsize=5.8, color=RED)
ax.legend(handles=[Line2D([], [], marker='o', mfc=BLUE, mec=NAVY, ls='none', ms=4.5, label='ITB sectioned'),
                   Line2D([], [], marker='s', mfc=GOLD, mec=NAVY, ls='none', ms=4.5, label='ALL added (incremental)')], frameon=False, fontsize=6.2, loc='upper left')
ax = fig.add_subplot(gs[1, 1]); letter(ax, 'd')
pct = lambda a, st: SE[(SE.angle == a) & (SE.step == st)].pct_of_intact.iloc[0]
bars = [([pct(60, 'ITB sectioned'), 14.0, 25.1, np.nan], BLUE, 'ITB, 60°'), ([pct(90, 'ITB sectioned'), 11.0, 31.1, np.nan], LIGHT, 'ITB, 90°'),
        ([pct(60, 'ITB + ALL'), 20.9, np.nan, np.nan], GOLD, 'ITB + ALL, 60°'), ([pct(90, 'ITB + ALL'), 19.5, np.nan, 63.0], DGOLD, 'ITB + ALL, 90°')]
for k_, (vals, col, lab) in enumerate(bars): ax.bar(np.arange(4) + (k_ - 1.5) * 0.2, vals, 0.2, color=col, edgecolor=NAVY, lw=0.5, hatch=['', '////', '', '....'][k_], label=lab)
ax.set_xticks(range(4)); ax.set_xticklabels(['Vietnamese\n(goniometer)', 'Huser 2017\n(robot)', 'Vap 2017\n(robot)', 'Sonnery-\nCottet 2016\n(navigation)'], fontsize=6.2)
ax.set_ylabel('Increase vs intact knee (%)'); ax.legend(frameon=False, fontsize=6.2, loc='upper left'); ax.grid(axis='y', color=GREY, lw=0.5); ax.set_axisbelow(True)
fig.savefig(OUT + 'Fig1.tif', dpi=600, bbox_inches='tight', pil_kwargs={'compression': 'tiff_lzw'}); fig.savefig(OUT + 'Fig1.png', dpi=300, bbox_inches='tight'); plt.show()

fig = plt.figure(figsize=(7.1, 6.6)); gs = fig.add_gridspec(2, 2, hspace=0.5, wspace=0.32)
ax = fig.add_subplot(gs[0, 0], projection='polar'); ax.text(-0.14, 1.04, 'a', transform=ax.transAxes, fontsize=11, fontweight='bold', color=NAVY)
for lab, st, m_, col, nm in [('mean', '-', 'o', NAVY, f'{r0:.1f} mm (mean)'), ('minusSD', '--', 's', STEEL, f'{r0-sdr:.1f} mm (−1 SD)'), ('plusSD', ':', '^', DGOLD, f'{r0+sdr:.1f} mm (+1 SD)')]:
    d = ring[lab]; th = np.radians(list(d.direction_deg) + [0]); rr = list(d.range_mm) + [d.range_mm.iloc[0]]
    ax.plot(th, rr, st, color=col, lw=1.2, marker=m_, ms=3, mfc=col, mec=col, label=nm)
ax.fill(np.linspace(0, 2 * np.pi, 100), [2] * 100, color=GOLD, alpha=0.35, lw=0)
ax.set_theta_zero_location('N'); ax.set_theta_direction(-1); ax.set_xticks(np.radians([0, 90, 180, 270])); ax.set_xticklabels(['Proximal', 'Posterior', 'Distal', 'Anterior'], fontsize=6.5)
ax.set_rlabel_position(112); ax.set_yticks([4, 8, 12, 16]); ax.set_yticklabels(['4', '8', '12', '16 mm']); ax.tick_params(axis='y', labelsize=5.8); ax.set_ylim(0, 18); ax.grid(color=GREY, lw=0.5)
ax.legend(frameon=False, fontsize=6, loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=3)
ax = fig.add_subplot(gs[0, 1]); letter(ax, 'b')
for (c, col, m_, ls_) in [('Vietnamese origin | Vietnamese ratio', NAVY, 'o', '-'), ('Vietnamese origin | Vietnamese 16 mm', GOLD, 's', '-'),
                          ('Vietnamese origin | Western 21.6 mm', BLUE, '^', '--'), ('Western origin | Vietnamese ratio', STEEL, 'D', ':')]:
    ax.plot(flex, CUR[c], ls_, color=col, marker=m_, ms=3, mfc=col, mec=NAVY if col == GOLD else col, mew=0.4, lw=1.3, label=c.replace(' | ', ' + '))
ax.axhline(0, color=NAVY, lw=0.4); ax.axvspan(70, 90, color=GREY, alpha=0.6, lw=0); ax.text(80, 17.3, 'kinematics\nheld at 70°', fontsize=5.6, ha='center', va='top')
ax.set_xlabel('Knee flexion (°)'); ax.set_ylabel('Length change vs extension (%)'); ax.set_xlim(0, 90); ax.set_ylim(-6, 18)
ax.legend(frameon=False, fontsize=5.6, loc='upper center', bbox_to_anchor=(0.5, -0.16), ncol=2)
ax = fig.add_subplot(gs[1, 0]); letter(ax, 'c'); pos, lab, i = [], [], 0
for p in scen:
    for a in (30, 60, 90):
        r = OUTC[(OUTC.procedure == p) & (OUTC.angle == a)].iloc[0]
        ax.bar(i, r.under * 100, 0.75, color=BLUE, edgecolor='white', lw=0.5, hatch='////', label='Under-constrained (> +2°)' if i == 0 else None)
        ax.bar(i, r.restored * 100, 0.75, bottom=r.under * 100, color=LIGHT, edgecolor='white', lw=0.5, label='Restored (±2°)' if i == 0 else None)
        ax.bar(i, r.over * 100, 0.75, bottom=(r.under + r.restored) * 100, color=GOLD, edgecolor=NAVY, lw=0.4, hatch='....', label='Over-constrained (< −2°)' if i == 0 else None)
        pos.append(i); lab.append(f'{a}°'); i += 1
    i += 0.6
ax.set_xticks(pos); ax.set_xticklabels(lab, fontsize=6.5); ax.set_ylabel('Virtual knees (%)'); ax.set_ylim(0, 100)
for k_, p in enumerate(scen): ax.text(np.mean(pos[k_ * 3:k_ * 3 + 3]), -12, p, ha='center', fontsize=7)
ax.legend(frameon=False, fontsize=5.8, loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3)
ax = fig.add_subplot(gs[1, 1]); letter(ax, 'd')
A_ = sorted(t8.aim_anterior.unique()); B_ = sorted(t8.aim_proximal.unique())
Z = np.array([[sel(a, b).deep_mm for a in A_] for b in B_])
from matplotlib.colors import LinearSegmentedColormap
cm = LinearSegmentedColormap.from_list('knee', [GOLD, '#F6F2E4', LIGHT, BLUE])
cf = ax.contourf(A_, B_, Z, levels=np.arange(4, 8.6, 0.5), cmap=cm)
cs = ax.contour(A_, B_, Z, levels=[4, 5, 6, 7, 8], colors=NAVY, linewidths=0.7); ax.clabel(cs, fmt='%d mm', fontsize=6)
ax.plot([0, 30, 45], [0, 30, 45], 'o', color=NAVY, ms=3.5)
for a in (0, 30, 45): ax.annotate(f'{sel(a, a).deep_mm:.1f} mm', (a, a), textcoords='offset points', xytext=((5, -9) if a < 45 else (-34, -10)), fontsize=6)
ax.set_xlabel('ALL socket aim, anterior (°)'); ax.set_ylabel('ALL socket aim, proximal (°)')
cb = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.03); cb.set_label('Clearance inside bone (mm)', fontsize=6.5); cb.ax.tick_params(labelsize=6)
fig.savefig(OUT + 'Fig2.tif', dpi=600, bbox_inches='tight', pil_kwargs={'compression': 'tiff_lzw'}); fig.savefig(OUT + 'Fig2.png', dpi=300, bbox_inches='tight'); plt.show()
print('\nAll outputs saved to', OUT); _LOG.close()
