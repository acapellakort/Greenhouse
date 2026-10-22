"""
check_tex.py -- the equations of docs/greenhouse_model_n79.tex, written in Python from the document, evaluated at the
states and inputs dumped by test/dump_terms.jl, and compared with the Julia terms (model_terms).
Usage (from GreenhouseTwin):  python3 tools/check_tex.py [test/terms_dump.csv]
Each row of the report: term, max relative difference over the dump; PASS if < 1e-9.
"""
import csv, json, math, sys

P = {k: v["val"] for k, v in json.load(open("configfiles/climate_n79.json")).items()}
rows = list(csv.DictReader(open(sys.argv[1] if len(sys.argv) > 1 else "test/terms_dump.csv")))
rows = [{k: float(v) for k, v in r.items()} for r in rows]
CO2o, Tmc = 834.7, 276.0
pi = math.pi
clamp = lambda x, a, b: min(max(x, a), b)
pos = lambda x: max(x, 0.0)

def Ps(T):                                   # Section "Psychrometrics"
    Tc = T - 273.15
    if T > 273.15:
        return 611.21 * math.exp((18.678 - Tc / 234.5) * (Tc / (257.14 + Tc)))
    return 611.15 * math.exp((23.036 - Tc / 333.7) * (Tc / (279.82 + Tc)))

def column(layers, rf):                      # Section "Solar radiation" (adding method)
    (t1, r1), (t2, r2), (t3, r3), (t4, r4) = layers
    R4 = r4 + t4 ** 2 * rf / (1 - r4 * rf); R3 = r3 + t3 ** 2 * R4 / (1 - r3 * R4)
    R2 = r2 + t2 ** 2 * R3 / (1 - r2 * R3); R1 = r1 + t1 ** 2 * R2 / (1 - r1 * R2)
    D2 = t1 / (1 - r1 * R2); D3 = t2 * D2 / (1 - r2 * R3); D4 = t3 * D3 / (1 - r3 * R4); D5 = t4 * D4 / (1 - r4 * rf)
    return dict(R1=R1, A1=(1 - t1 - r1) * (1 + R2 * D2), A2=(1 - t2 - r2) * (D2 + R3 * D3), A3=(1 - t3 - r3) * (D3 + R4 * D4),
                A4=(1 - t4 - r4) * (D4 + rf * D5), Af=(1 - rf) * D5, D4=D4, D5=D5)

def model(r):
    p = P; T1, T2, C1, V1, eint, T5, Td, Tm = (r[k] for k in ("T1", "T2", "C1", "V1", "eint", "T5", "Td", "Tm"))
    I = lambda k: r["in_" + k]
    Ig, Tsky, To, v, Vo, L = I("Iglob"), I("Tsky"), I("Tout"), I("Wind"), I("VPout"), I("LAI")
    U1, Ubl, U2, U3, U4, U5, U6, U7 = I("EnScr"), I("BlackScr"), I("FanPad"), I("MechCool"), I("Heater"), I("Shade"), I("SideVent"), I("ForcedVent")
    Ulee, Uwin, U9, U10, U11, U12, Tlog, Rc, FW = I("VentLee"), I("VentWind"), I("Fog"), I("CO2Valve"), I("HeatSP"), I("Lamps"), I("Tpipe"), I("CropResp"), I("CropFW")
    o = {}
    # solar column
    tr = (1 - p["eta1"]) * p["tau1"]; rr = 1 - tr - p["eta1"] - p["cov_abs"]
    es, bl = clamp(U1, 0, 1), clamp(Ubl, 0, 1)
    def band(tes, res, tbl, rbl, k, rinf, rf):
        tc = math.exp(-k * max(L, 0)); rcn = rinf * (1 - tc)
        return column([(tr, rr), (1 - es * (1 - tes), es * res), (1 - bl * (1 - tbl), bl * rbl), (tc, rcn)], rf)
    Pb = band(p["scr_es_tau_par"], p["scr_es_rho_par"], p["scr_bl_tau_par"], p["scr_bl_rho_par"], p["beta1"], p["rho1"], p["rho2"])
    Nb = band(p["scr_es_tau_nir"], p["scr_es_rho_nir"], p["scr_bl_tau_nir"], p["scr_bl_rho_nir"], p["beta3"], p["rho_can_nir"], 1 - p["alpha7"])
    IP, IN = p["eta2"] * Ig, p["eta3"] * Ig
    Qroof = IP * Pb["A1"] + IN * Nb["A1"]; fs = min(p["eta1"] / (1 - tr - rr), 1)
    Qst = Qroof * fs; Qscr = IP * (Pb["A2"] + Pb["A3"]) + IN * (Nb["A2"] + Nb["A3"])
    o["so_Qcan"] = IP * Pb["A4"] + IN * Nb["A4"]; o["so_Qflr"] = IP * Pb["Af"] + IN * Nb["Af"]
    o["so_Qair"] = Qst + p["scr_abs_air"] * Qscr; o["so_Qcov"] = (Qroof - Qst) + (1 - p["scr_abs_air"]) * Qscr
    o["so_Qout"] = IP * Pb["R1"] + IN * Nb["R1"]
    # plant light and FvCB
    PARp = IP * Pb["D4"] + IP * p["rho2"] * Pb["D5"] + p["lamp_par_W"] * clamp(U12, 0, 1); o["PARp"] = PARp
    T = T2; Rd = p["Rd_day"] * 2 ** ((T1 - 298.15) / 10)
    Iq = 4.6e-6 * PARp; Ci = 0.509e-6 * C1; g = p["g_s"] * p["g_m"] / (p["g_s"] + p["g_m"])
    Gs = 0.5 * p["O_a"] / (p["Sco25"] * math.exp((T - 298.15) / T * p["E_Soc"] / (298.15 * p["Rgas"])))
    a, J = p["alpha"], p["J_max"]; th = p["theta"]
    Jt = ((a * Iq + J) - math.sqrt((a * Iq + J) ** 2 - 4 * th * a * Iq * J)) / (2 * th)
    Vc = p["V_cmax25"] * p["Q10_Vcmax"] ** (0.1 * (T - 298.15)) / (1 + math.exp(0.128 * (T - 315.15)))
    KC = p["K_C25"] * p["Q10_KC"] ** (0.1 * (T - 298.15)); KO = p["K_O25"] * p["Q10_KO"] ** (0.1 * (T - 298.15))
    def root(b, c):                         # smaller root of x^2 + b x + c = 0 (document: x^2 - (...) x + ... = 0)
        d = b * b - 4 * c
        return 0.5 * (-b - math.sqrt(d)) if d >= 0 else -Rd
    Ac = root(-(Vc + g * (Ci + KC * (1 + p["O_a"] / KO)) - Rd), g * (Vc * pos(Ci - Gs) - (Ci + KC * (1 + p["O_a"] / KO)) * Rd))
    Aj = root(-(0.25 * Jt - Rd + g * (Ci + 2 * Gs)), 0.25 * g * Jt * pos(Ci - Gs) - g * (Ci + 2 * Gs) * Rd)
    Ap = 1.5 * Vc - Rd
    kp = p["k_lai_photo"]; o["A"] = 44010 * (1 - math.exp(-kp * L)) / kp * min(Aj, Ac, Ap)
    # heating
    e = U11 - T2; dmd = p["hc_Kp"] * e + p["hc_Ki"] * max(eint, 0); room = max(p["hc_Tmax"] - T2, 0)
    warm = Tlog > T2 + 1
    Tp = Tlog if warm else max(Tlog, T2 + clamp(dmd, 0, room)); o["Tp"] = Tp
    o["deint"] = 0.0 if (warm or (dmd >= room and e > 0) or (eint <= 0 and e < 0)) else e
    # stomata
    q2 = Ps(T1); q7 = 1 / (1 + math.exp((PARp - p["delta1"]) / p["gamma5"]))
    q8 = p["delta5"] + (p["delta4"] - p["delta5"]) * q7; q9 = p["delta7"] + (p["delta6"] - p["delta7"]) * q7
    q10 = (PARp + p["delta2"]) / (PARp + p["delta3"])
    q3 = p["gamma4"] * q10 * (1 + q8 * (C1 * T2 / 536.4 - 200) ** 2) * (1 + q9 * (q2 - V1) ** 2); o["q3"] = q3
    base = 2 * p["rho3"] * p["alpha5"] * L / (p["gamma"] * p["gamma2"])
    q1 = base / p["gamma3"] if q2 < V1 else base / (p["gamma3"] + q3); o["q1"] = q1
    # cover node
    kR = 4 * p["epsil3"] * p["sigma"] * (0.5 * (T2 + Tsky)) ** 3; h = p["h_cover_out"]
    Tc0 = (kR * (T2 + Tsky) + h * To + o["so_Qcov"]) / (2 * kR + h); fcov = 1 - U1 * (1 - p["vent_scr_perm"])
    Tc = Tc0
    for _ in range(4):
        Tc = Tc0 + p["gamma2"] * fcov * 6.4e-9 * 1.7 * abs(T2 - Tc) ** 0.33 * pos(V1 - Ps(Tc)) / (2 * kR + h)
    o["Tc"] = Tc
    # long-wave
    sg = p["sigma"]; g1 = 0.49 * (1 - math.exp(-p["beta3"] * L)); g2 = 1 - U1 * (1 - p["tau3"])
    pf = 1 - 0.49 * pi * p["gamma1"] * p["phi1"]; g3 = p["tau1"] * g2 * pf * math.exp(-p["beta3"] * L)
    F = p["epsil3"] * sg * (T2 ** 4 - Tc ** 4)
    o["r6"] = p["alpha3"] * p["epsil1"] * p["epsil2"] * g1 * sg * (Tp ** 4 - T1 ** 4)
    o["r7"] = g1 / 0.49 * p["epsil2"] * p["epsil3"] * g2 * sg * (T1 ** 4 - Tc ** 4)
    o["r11"] = p["epsil4"] * g3 * p["epsil3"] * sg * (T5 ** 4 - Tc ** 4)
    o["r12"] = p["eps_screen"] * U1 * p["tau2"] * F; o["r13"] = p["cover_area"] * p["epsil6"] * F
    o["r14"] = p["epsil4"] * p["epsil2"] * pf * sg * (T1 ** 4 - T5 ** 4)
    # air exchange
    n3 = p["nu2"] * (1 - p["eta11"] * U5); n1 = p["nu1"] * (1 - p["eta10"] * U5)
    B = p["omega1"] * p["nu6"] * (T2 - To) / (T2 + To); w = n3 * v ** 2; ex = p["vs_exp"]
    Ul, Uw = clamp(Ulee, 0, 1), clamp(Uwin, 0, 1)
    Phi = 0.5 * (Ul ** ex * math.sqrt(pos(B + p["vs_clee"] * w)) + p["vs_gwind"] * Uw ** ex * math.sqrt(pos(B + p["vs_cwind"] * w)))
    U8 = Phi / math.sqrt(pos(B + w)) if pos(B + w) > 1e-12 else 0.5 * (Ul ** ex + p["vs_gwind"] * Uw ** ex); o["U8"] = U8
    fscr = 1 - U1 * (1 - p["vent_scr_perm"])
    f7 = U8 * fscr * p["nu5"] * n1 / (2 * p["alpha6"]) * math.sqrt(pos(B + w)); o["f7"] = f7
    f4 = p["eta6"] * f7; f5 = n1 * p["nu3"] * U6 * v * math.sqrt(n3) / (2 * p["alpha6"]); f2 = p["eta6"] * U1 * f5
    f3 = U7 * p["phi8"] / p["alpha6"]; f6 = p["nu4"] * math.sqrt(max(v, 0.25) ** 2 + p["leak_cs"] * abs(T2 - To))
    o.update(f2=f2, f3=f3, f4=f4, f6=f6); ft = f2 + f3 + f4 + f6
    # lamps, canopy, air
    Ql = p["lamp_Pel"] * U12; Qlc = p["lamp_fcan"] * (1 - math.exp(-p["beta1"] * L)) * Ql; Qla = Ql - Qlc
    o["Qlc"], o["Qla"] = Qlc, Qla
    h1 = 2 * p["alpha4"] * L * (T1 - T2); l1 = p["gamma2"] * q1 * (q2 - V1); o["h1"], o["l1"] = h1, l1
    o["dT1"] = (o["so_Qcan"] + Qlc + o["r6"] - h1 - l1 - o["r7"] - o["r14"]) / (p["alpha1"] * L)
    a6 = p["alpha6"]
    h2 = U2 * p["phi7"] / a6 * (p["rho3"] * p["alpha5"] * To - p["gamma2"] * p["rho3"] * p["eta5"] * (p["phi5"] - p["phi6"]))
    h4 = p["n_pipes"] * 1.99 * pi * p["phi1"] * p["gamma1"] * abs(Tp - T2) ** 0.32 * (Tp - T2)
    h7 = p["rho3"] * p["alpha5"] * ft * (T2 - To); h11 = 2 * p["nu7"] * (T2 - T5) / (p["phi2"] + p["nu8"])
    l2 = p["gamma2"] * U9 * p["phi9"] / a6; Qaf = p["k_ground"] * (T2 - T5); Qam = p["h_mass"] * (T2 - Tm)
    o.update(h2=h2, h4=h4, h7=h7, h11=h11, l2=l2, Qaf=Qaf, Qam=Qam)
    cair = p["phi2"] * p["rho3"] * p["alpha5"]
    o["dT2"] = (h1 + h2 + h4 - h7 - h11 - o["r12"] - o["r13"] - l2 + o["so_Qair"] - Qaf - Qam + Qla) / (p["air_cap_mult"] * cair)
    # vapour
    Tscr = 0.5 * (To + T2); ps1, om2 = p["psi1"], p["omega2"]
    p1 = q1 * (q2 - V1); p2 = p["rho3"] * U2 * p["phi7"] / a6 * (p["eta5"] * (p["phi5"] - p["phi6"]) + p["phi6"])
    p3 = U9 * p["phi9"] / a6; p4 = p["eta12"] * U4 * p["lamb4"] / a6; p6 = U2 * p["phi7"] / a6 * ps1 / om2 * V1 / T2
    p5 = ps1 / om2 * (V1 / T2 - Vo / To) * ft; p7 = U3 * 6.4e-9 * pos(V1 - Ps(Tmc))
    p8 = 6.4e-9 * 1.7 * U1 * abs(T2 - Tscr) ** 0.33 * pos(V1 - Ps(Tscr))
    p9 = 6.4e-9 * 1200 * abs(T1 - T2) ** 0.33 * pos(V1 - Ps(T2))
    p10 = p["cover_area"] * fcov * 6.4e-9 * 1.7 * abs(T2 - Tc) ** 0.33 * pos(V1 - Ps(Tc))
    o.update(p1=p1, p2=p2, p3=p3, p4=p4, p5=p5, p6=p6, p7=p7, p8=p8, p9=p9, p10=p10)
    k3 = ps1 * p["phi2"] / (om2 * T2)
    o["dV1"] = (p1 + p2 + p3 + p4 - p5 - p6 - p7 - p8 - p9 - p10) / k3
    # CO2
    o1 = p["eta13"] * U4 * p["lamb4"] / a6; o2 = U10 * p["psi2"] / a6; o3 = U2 * p["phi7"] / a6 * (C1 - CO2o)
    o5 = (C1 - CO2o) * ft; oR = Rc * 2 ** ((T1 - 293.15) / 10); o.update(o1=o1, o2=o2, o3=o3, o5=o5, oR=oR)
    o["dC1"] = (o1 + o2 + o3 - o["A"] - o5 + oR) / p["phi2"]
    # floor, deep soil, mass
    Qfd = p["k_deep"] * (T5 - Td); Qfs = (1 - p["f_mass"]) * o["so_Qflr"]; Qms = p["f_mass"] * o["so_Qflr"]
    o.update(Qfd=Qfd, Qfs=Qfs, Qms=Qms)
    o["dT5"] = (Qaf - Qfd + h11 + o["r14"] - o["r11"] + Qfs) / p["C_soil"]
    o["dTd"] = (T5 - Td) / p["sd_tau"]
    o["dTm"] = (Qam + Qms) / (p["C_mass"] + p["crop_cp"] * FW)
    return o

worst = {}
for r in rows:
    for k, v in model(r).items():
        ref = r[k]; den = max(abs(ref), abs(v), 1e-300)
        d = 0.0 if ref == v else abs(ref - v) / den
        sc = max(abs(r2[k]) for r2 in rows)
        d = abs(ref - v) / max(sc, 1e-300)
        worst[k] = max(worst.get(k, 0.0), d)
bad = 0
for k in sorted(worst):
    ok = worst[k] < 1e-9; bad += not ok
    print(f"{k:10s} max |doc - code| / max|code| = {worst[k]:.2e}  {'PASS' if ok else 'FAIL'}")
print(f"\n{len(worst) - bad} PASS, {bad} FAIL  ({len(rows)} states)")
sys.exit(1 if bad else 0)
