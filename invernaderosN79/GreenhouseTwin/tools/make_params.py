"""
make_params.py -- builds the flag-free n79 parameter files of GreenhouseTwin from the V2.2 JSONs.

Provenance script (run once; kept for traceability). Value of every parameter =
  V2.2 JSON value  ->  overridden by the n79 fixed list (run_climate_calibration_cl.jl, PRIORS=n79, SOLARC=1)
                   ->  inferred parameters at their n79 prior centres.
Usage: python3 tools/make_params.py <invernaderos root> <output configfiles dir>
"""
import json, sys, os

ROOT, OUT = sys.argv[1], sys.argv[2]
cl = json.load(open(os.path.join(ROOT, "invernaderos2.0/code/configfiles/constants_climate.json")))
cp = json.load(open(os.path.join(ROOT, "invernaderos2.0/code/configfiles/constants_cropphoto.json")))
gr = json.load(open(os.path.join(ROOT, "invernaderosV2.2/GreenhouseSim/configfiles/cucumber_growth.json")))
v = lambda d, k: float(d[k]["val"])

# --- n79 fixed values (FIXN/FIXV of run_climate_calibration_cl.jl, PRIORS=n79) ----------------------------------
N79_FIXED = {"phi2": 5.4, "alpha5": 1005.0, "nu5": 3000.0, "gamma1": 1.0, "n_pipes": 1.0, "tau1": 0.65 / 0.9,
             "psi2": 35556.0, "cover_area": 1.6, "lamp_Pel": 117.0, "crop_cp": 3800.0, "sd_tau": 14 * 86400.0,
             "hc_Kp": 10.0, "hc_Ki": 10.0 / 1800, "vs_clee": 1.0}
# --- n79 priors [min, centre, max] (CQ_dict after the n79 block, SOLARC=1) ----------------------------------------
N79_PRIORS = {"gamma4": [40.0, 82.0, 200.0], "vent_scr_perm": [0.0, 0.4, 1.0], "k_ground": [1.0, 8.0, 40.0],
              "nu1": [0.3, 0.65, 1.0], "nu2": [0.005, 0.05, 0.3], "C_soil": [1.5e5, 6e5, 2e6], "k_deep": [0.5, 1.5, 3.0],
              "C_mass": [2e4, 5e4, 1.2e5], "h_mass": [1.0, 5.0, 40.0], "f_mass": [0.0, 0.3, 1.0],
              "lamp_fcan": [0.05, 0.3, 0.9], "vs_gwind": [0.5, 1.5, 4.0], "vs_cwind": [0.5, 1.0, 10.0],
              "vs_exp": [0.5, 0.9, 1.2], "nu4": [1e-5, 4.5e-5, 1e-4], "h_cover_out": [10.0, 25.0, 50.0],
              "air_cap_mult": [1.0, 3.0, 12.0], "leak_cs": [0.0, 0.3, 1.0]}
NOISE = {"sigma_T": [0.2, 1.5, 6.0], "sigma_RH": [0.5, 6.0, 25.0], "sigma_C": [10.0, 60.0, 400.0]}

# name: (group, units, description, source)   -- climate + photosynthesis
C = {
 # geometry and structure
 "alpha6": ("geometry", "m2", "floor area used to express absolute capacities per m2 floor", "V2.0"),
 "phi2": ("geometry", "m", "mean air height (volume per floor area)", "Note 67 (compartment drawing)"),
 "cover_area": ("geometry", "m2 m-2", "exterior cover area per floor area (roof + one side wall)", "Note 67"),
 "nu5": ("geometry", "m2", "roof-vent area (0.30 of the floor, insect nets)", "Notes 78, 92"),
 "nu6": ("geometry", "m", "vertical dimension of a roof-vent opening", "V2.0"),
 "nu3": ("geometry", "m2", "side-wall vent surface (side vents, input U6)", "V2.0"),
 "gamma1": ("geometry", "m m-2", "heating-pipe length per floor area (rail pipe, one circuit)", "Note 78"),
 "phi1": ("geometry", "m", "outer diameter of the heating pipe", "V2.0 / Note 78"),
 "n_pipes": ("geometry", "-", "number of heating circuits (crop pipe folded into the rail pipe by the inputs)", "Note 54/79"),
 "alpha3": ("geometry", "m2 m-2", "pipe surface seen by the canopy (long-wave)", "V2.0 (audit: geometric 0.16)"),
 "nu8": ("geometry", "m", "floor-to-ground distance in the air-floor conduction length", "V2.0"),
 # physical constants
 "alpha5": ("constants", "J kg-1 K-1", "specific heat of air", "Note 78"),
 "rho3": ("constants", "kg m-3", "air density", "V2.0"),
 "gamma": ("constants", "Pa K-1", "psychrometric constant", "V2.0"),
 "gamma2": ("constants", "J kg-1", "latent heat of evaporation", "V2.0"),
 "sigma": ("constants", "W m-2 K-4", "Stefan-Boltzmann constant", "V2.0"),
 "omega1": ("constants", "m s-2", "gravity", "V2.0"),
 "omega2": ("constants", "J kmol-1 K-1", "molar gas constant", "V2.0"),
 "psi1": ("constants", "kg kmol-1", "molar mass of water", "V2.0"),
 # canopy
 "alpha1": ("canopy", "J K-1 m-2 leaf", "canopy heat capacity per leaf area", "V2.0"),
 "alpha4": ("canopy", "W m-2 K-1", "leaf-air convective coefficient (per leaf side)", "Vanthoor 2011"),
 "gamma3": ("canopy", "s m-1", "canopy boundary-layer resistance to vapour", "V2.0"),
 "gamma5": ("canopy", "W m-2", "inverse slope of the day/night switch of the stomatal resistance", "V2.0"),
 "delta1": ("canopy", "W m-2", "PAR defining day/night for the stomatal resistance", "V2.0"),
 "delta2": ("canopy", "W m-2", "stomatal light response parameter", "V2.0"),
 "delta3": ("canopy", "W m-2", "stomatal light response parameter", "V2.0"),
 "delta4": ("canopy", "ppm-2", "stomatal CO2 response, day", "V2.0"),
 "delta5": ("canopy", "ppm-2", "stomatal CO2 response, night", "V2.0"),
 "delta6": ("canopy", "Pa-2", "stomatal vapour-pressure response, day", "V2.0"),
 "delta7": ("canopy", "Pa-2", "stomatal vapour-pressure response, night", "V2.0"),
 "beta1": ("canopy", "-", "canopy PAR extinction coefficient", "V2.0"),
 "beta3": ("canopy", "-", "canopy NIR extinction coefficient (also used for the long-wave terms)", "V2.0 (audit: K_FIR 0.94)"),
 "rho1": ("canopy", "-", "canopy PAR reflection at high LAI", "V2.0"),
 "rho_can_nir": ("canopy", "-", "canopy NIR reflection at high LAI", "Vanthoor 2011 / Note 103"),
 "epsil2": ("canopy", "-", "canopy long-wave emissivity", "V2.0"),
 "k_lai_photo": ("canopy", "-", "extinction in the effective leaf area (1-exp(-kL))/k of the climate photosynthesis", "V2.2 K_EXT"),
 # floor and soil
 "rho2": ("floor", "-", "floor PAR reflection", "V2.0"),
 "alpha7": ("floor", "-", "floor NIR absorption (reflection 1 - alpha7)", "V2.0"),
 "epsil4": ("floor", "-", "floor long-wave emissivity", "V2.0"),
 "nu7": ("floor", "W m-1 K-1", "soil conductivity in the air-floor conduction h11", "V2.0"),
 "sd_tau": ("floor", "s", "time constant of the deep-soil temperature", "Note 47"),
 # cover and roof
 "tau1": ("cover", "-", "cover PAR transmission; (1-eta1) tau1 = 0.65 roof transmission", "Note 78 (organisers 0.60-0.70)"),
 "eta1": ("cover", "-", "share of the global radiation absorbed by the structure (booked to the air)", "V2.0"),
 "eta2": ("cover", "-", "PAR share of the global radiation", "V2.0"),
 "eta3": ("cover", "-", "NIR share of the global radiation", "V2.0"),
 "cov_abs": ("cover", "-", "share of the global radiation absorbed by the glass (booked to the cover)", "Note 46"),
 "epsil3": ("cover", "-", "sky/cover long-wave emissivity", "V2.0"),
 "epsil6": ("cover", "-", "cover long-wave emissivity (air-cover term r13)", "V2.0"),
 "tau2": ("cover", "-", "cover long-wave transmission in the screen term r12", "V2.0"),
 # screens
 "tau3": ("screens", "-", "energy-screen long-wave transmission", "V2.0"),
 "eps_screen": ("screens", "-", "energy-screen long-wave emissivity (upper face, r12)", "V2.2"),
 "scr_es_tau_par": ("screens", "-", "energy screen PAR transmission at full closure (PROVISIONAL)", "organisers 25 % interception"),
 "scr_es_rho_par": ("screens", "-", "energy screen PAR reflection (PROVISIONAL, assumed)", "Note 101"),
 "scr_es_tau_nir": ("screens", "-", "energy screen NIR transmission (PROVISIONAL, = PAR)", "Note 101"),
 "scr_es_rho_nir": ("screens", "-", "energy screen NIR reflection (PROVISIONAL, = PAR)", "Note 101"),
 "scr_bl_tau_par": ("screens", "-", "blackout screen PAR transmission at full closure (PROVISIONAL)", "organisers 98 % interception"),
 "scr_bl_rho_par": ("screens", "-", "blackout screen PAR reflection (PROVISIONAL, assumed)", "Note 101"),
 "scr_bl_tau_nir": ("screens", "-", "blackout screen NIR transmission (PROVISIONAL, = PAR)", "Note 101"),
 "scr_bl_rho_nir": ("screens", "-", "blackout screen NIR reflection (PROVISIONAL, = PAR)", "Note 101"),
 "scr_abs_air": ("screens", "-", "share of screen-absorbed solar released to the crop air (rest to the cover)", "Note 101 (PROVISIONAL)"),
 # heating
 "epsil1": ("heating", "-", "pipe long-wave emissivity", "V2.0"),
 "hc_Kp": ("heating", "K K-1", "PI heating controller proportional gain", "Note 34"),
 "hc_Ki": ("heating", "s-1", "PI heating controller integral gain (tau_i 30 min)", "Note 34"),
 "hc_Tmax": ("heating", "K", "maximum pipe temperature (65 C)", "Note 33"),
 # ventilation and leakage
 "eta6": ("ventilation", "-", "ventilation reduction factor (Vanthoor 8.70, zeta 0.33)", "V2.0"),
 "eta10": ("ventilation", "-", "shading effect on the discharge coefficient", "V2.0"),
 "eta11": ("ventilation", "-", "shading effect on the wind-pressure coefficient", "V2.0"),
 "vs_clee": ("ventilation", "-", "leeward wind-pressure multiplier", "Note 58"),
 # actuators kept as inputs (zero input in the AiCU season)
 "eta5": ("actuators", "-", "fan-pad efficiency (U2)", "V2.0"),
 "phi5": ("actuators", "kg kg-1", "vapour content of fan-pad air (U2)", "V2.0"),
 "phi6": ("actuators", "kg kg-1", "vapour content of outside air for the fan-pad (U2)", "V2.0"),
 "phi7": ("actuators", "m3 s-1", "fan-pad air flow capacity (U2)", "V2.0"),
 "phi8": ("actuators", "m3 s-1", "forced-ventilation capacity (U7)", "V2.0"),
 "phi9": ("actuators", "kg s-1", "fog capacity (U9)", "V2.0"),
 "lamb4": ("actuators", "W", "direct air heater capacity (U4)", "V2.0"),
 "eta12": ("actuators", "kg J-1", "vapour released per J of direct-heater heat (U4)", "V2.0"),
 "eta13": ("actuators", "mg J-1", "CO2 released per J of direct-heater heat (U4)", "V2.0"),
 "psi2": ("actuators", "mg s-1", "CO2 valve capacity (12.8 g m-2 h-1 x alpha6)", "Note 78"),
 # internal mass
 "crop_cp": ("mass", "J kg-1 K-1", "specific heat of fresh plant tissue", "Note 91"),
 # lamps
 "lamp_Pel": ("lamps", "W m-2", "lamp electrical power at 100 % (187 umol / 1.6 umol J-1)", "Note 90"),
 "lamp_par_W": ("lamps", "W m-2", "lamp PAR at the canopy at 100 % (187 umol / 4.6)", "Note 71"),
 # photosynthesis (FvCB)
 "O_a": ("photosynthesis", "mol m-3", "O2 concentration", "V2.0"),
 "Sco25": ("photosynthesis", "-", "Rubisco CO2/O2 specificity at 25 C", "V2.0"),
 "E_Soc": ("photosynthesis", "J mol-1", "activation energy of Sco", "V2.0"),
 "Rgas": ("photosynthesis", "J mol-1 K-1", "gas constant", "V2.0"),
 "alpha": ("photosynthesis", "-", "quantum efficiency of electron transport", "V2.0"),
 "J_max": ("photosynthesis", "mol m-2 s-1", "maximum electron transport rate (climate model; the crop calibration infers its own)", "V2.0"),
 "theta": ("photosynthesis", "-", "curvature of the light response", "V2.0"),
 "V_cmax25": ("photosynthesis", "mol m-2 s-1", "maximum carboxylation rate at 25 C", "V2.0"),
 "Q10_Vcmax": ("photosynthesis", "-", "Q10 of Vcmax", "V2.0"),
 "K_C25": ("photosynthesis", "mol m-3", "Michaelis constant for CO2 at 25 C", "V2.0"),
 "Q10_KC": ("photosynthesis", "-", "Q10 of K_C", "V2.0"),
 "K_O25": ("photosynthesis", "mol m-3", "Michaelis constant for O2 at 25 C", "V2.0"),
 "Q10_KO": ("photosynthesis", "-", "Q10 of K_O", "V2.0"),
 "Rd_day": ("photosynthesis", "mol m-2 s-1", "leaf dark respiration at 25 C (x Q10 2 on canopy T)", "V2.0 / Note 86"),
 "g_s": ("photosynthesis", "mol m-2 s-1", "stomatal conductance in the constant effective conductance", "Bush 2023 (V2.0)"),
 "g_m": ("photosynthesis", "mol m-2 s-1", "mesophyll conductance in the constant effective conductance", "Bush 2023 (V2.0)"),
}
CI = {  # inferred (n79): group, units, description
 "gamma4": ("canopy", "s m-1", "minimum stomatal resistance"),
 "vent_scr_perm": ("ventilation", "-", "share of roof ventilation and cover contact reaching the crop air through a closed energy screen"),
 "k_ground": ("floor", "W m-2 K-1", "air-floor convective conductance"),
 "nu1": ("ventilation", "-", "roof-vent discharge coefficient"),
 "nu2": ("ventilation", "-", "roof-vent wind-pressure coefficient"),
 "C_soil": ("floor", "J m-2 K-1", "heat capacity of the active floor layer"),
 "k_deep": ("floor", "W m-2 K-1", "floor to deep-soil conductance"),
 "C_mass": ("mass", "J m-2 K-1", "heat capacity of the internal mass (structure, slabs, water)"),
 "h_mass": ("mass", "W m-2 K-1", "air-internal mass heat-transfer coefficient"),
 "f_mass": ("mass", "-", "share of the floor-level solar absorbed by the internal mass"),
 "lamp_fcan": ("lamps", "-", "share of lamp power absorbed by a closed canopy"),
 "vs_gwind": ("ventilation", "-", "windward discharge multiplier relative to leeward"),
 "vs_cwind": ("ventilation", "-", "windward wind-pressure multiplier"),
 "vs_exp": ("ventilation", "-", "curvature of opening -> flow"),
 "nu4": ("ventilation", "m s-1 per m s-1", "leakage coefficient"),
 "h_cover_out": ("cover", "W m-2 K-1", "cover-outside convective coefficient"),
 "air_cap_mult": ("air", "-", "effective air heat-capacity multiplier (fast masses lumped with the air)"),
 "leak_cs": ("ventilation", "m2 s-2 K-1", "stack (temperature-difference) term of the leakage"),
}
src_all = {**cl, **cp}
special = {"k_lai_photo": 0.7, "g_s": 0.3, "g_m": 0.25}
out = {}
for k, (grp, un, ds, so) in C.items():
    val = N79_FIXED[k] if k in N79_FIXED else special[k] if k in special else float(src_all[k]["val"])
    out[k] = {"val": val, "units": un, "group": grp, "status": "fixed", "desc": ds, "source": so}
for k, (grp, un, ds) in CI.items():
    out[k] = {"val": N79_PRIORS[k][1], "units": un, "group": grp, "status": "inferred", "prior": N79_PRIORS[k],
              "desc": ds, "source": "n79 prior (Note 78; f_mass Note 104)"}
for k, pr in NOISE.items():
    out[k] = {"val": pr[1], "units": {"sigma_T": "K", "sigma_RH": "%", "sigma_C": "mg m-3"}[k], "group": "noise",
              "status": "inferred", "prior": pr, "desc": "likelihood noise", "source": "n79 prior"}
json.dump(out, open(os.path.join(OUT, "climate_n79.json"), "w"), indent=1)

# --- crop (n105) ----------------------------------------------------------------------------------------------------
CROP_PRIORS = {"W_pot": [20.0, 60.7, 120.0], "E_fruit": [150.0, 275.0, 450.0], "app_beta1": [0.05, 0.0882, 0.13],
               "V20": [1.0, 3.2, 8.0], "k_m_leaf": [0.01, 0.03, 0.06], "k_m_stem": [0.005, 0.015, 0.03],
               "k_m_fruit": [0.003, 0.01, 0.02], "gen_leaf_target": [20.0, 25.0, 30.0], "gen_N0": [2.0, 4.0, 8.0],
               "set_r_low": [0.0, 0.3, 1.0], "set_r_high": [0.2, 1.0, 2.5], "DMC_fruit": [0.025, 0.035, 0.045],
               "J_max": [1.0e-4, 2.0e-4, 3.5e-4],
               "sigma_y": [0.02, 0.3, 3.0], "sigma_leaf": [0.5, 3.0, 20.0], "sigma_fd": [0.5, 3.0, 15.0], "sigma_n": [0.3, 3.0, 20.0]}
CD = {
 "W_pot": ("M94", "g DM fruit-1", "potential fruit dry mass"), "E_fruit": ("M94", "C d", "fruit duration (degree-days above 10 C)"),
 "app_beta1": ("M94", "nodes stem-1 d-1 C-1", "node (flower) appearance slope"), "V20": ("M94", "g DM stem-1 d-1", "vegetative demand at 20 C"),
 "k_m_leaf": ("respiration", "g g-1 d-1", "leaf maintenance at 25 C"), "k_m_stem": ("respiration", "g g-1 d-1", "stem maintenance at 25 C"),
 "k_m_fruit": ("respiration", "g g-1 d-1", "fruit maintenance at 25 C"), "gen_leaf_target": ("canopy", "leaves stem-1", "retained-leaf target (pruning)"),
 "gen_N0": ("canopy", "nodes", "nodes at transplant"), "set_r_low": ("fruit set", "-", "source/sink ratio with no fruit set"),
 "set_r_high": ("fruit set", "-", "source/sink ratio with full fruit set"), "DMC_fruit": ("quality", "-", "fruit dry-matter content"),
 "J_max": ("photosynthesis", "mol m-2 s-1", "maximum electron transport rate (crop path)"),
}
CF = {  # fixed crop values: name -> (value, group, units, desc, source)
 "stem_density": (3.6, "canopy", "stems m-2", "stem density (one stem = one modelled plant)", "Hemming et al. 2019 (Note 92)"),
 "gen_fL": (0.59, "canopy", "-", "leaf share of vegetative growth", "Aaron (Note 94)"),
 "frac_stem": (v(gr, "frac_stem") if isinstance(gr["frac_stem"], dict) else float(gr["frac_stem"]), "canopy", "-", "stem weight in the stem:root split", "V2.2"),
 "frac_root": (float(gr["frac_root"]), "canopy", "-", "root weight in the stem:root split", "V2.2"),
 "sla_mean": (float(gr["sla_mean"]), "canopy", "cm2 g-1", "seasonal SLA mean", "Heuvelink 1999 eq. 3"),
 "sla_amp": (float(gr["sla_amp"]), "canopy", "cm2 g-1", "seasonal SLA amplitude", "Heuvelink 1999 eq. 3"),
 "sla_phase": (float(gr["sla_phase"]), "canopy", "d", "seasonal SLA phase", "Heuvelink 1999 eq. 3"),
 "gen_doy0": (226.0, "canopy", "day of year", "transplant day (14 Aug)", "Note 94"),
 "LAI0": (0.25, "canopy", "m2 m-2", "LAI at transplant (assumed)", "Note 94"),
 "gen_flower_node": (5.0, "architecture", "node", "first node with a flower", "Aaron (Note 94)"),
 "set_start_day": (15.0, "architecture", "d", "first day with fruit set", "Note 66"),
 "top_day": (float(gr["top_day"]), "architecture", "d", "topping: no node appearance after this day", "V2.2"),
 "k_m_root": (float(gr["k_m_root"]), "respiration", "g g-1 d-1", "root maintenance at 25 C", "V2.2"),
 "Q10_resp": (float(gr["Q10_resp"]), "respiration", "-", "Q10 of maintenance", "V2.2"),
 "asrq": (float(gr["asrq"]), "respiration", "g CH2O g-1 DM", "construction cost", "Aaron (1.4)"),
 "k_ext": (float(gr["k_ext"]), "photosynthesis", "-", "extinction in the effective leaf area of the crop path", "V2.2 (climate uses 0.7)"),
 "harvest_FW": (323.0, "harvest", "g FW", "picking size (observed mean A+B fruit)", "Note 92"),
 "transmission": (0.65, "light", "-", "greenhouse PAR transmission in the measured-light formula", "Note 78"),
 "lamp_umol": (187.0, "light", "umol m-2 s-1", "lamp PAR at 100 %", "AGC ReadMe"),
 "VT": (0.25, "M94", "g DM stem-1 d-1 C-1", "temperature slope of vegetative demand", "Marcelis 1994"),
 "sink_s": (0.0170 * 275.0, "M94", "-", "fruit sink shape s = 0.0170 x 275", "Marcelis 1994 via Aaron"),
 "sink_p": (1.0 + 1.0 / 0.0111, "M94", "-", "fruit sink exponent p = 1 + 1/0.0111", "Marcelis 1994 via Aaron"),
 "sink_peak": (131.0 / 275.0, "M94", "-", "stage of peak sink = 131/275", "Marcelis 1994 via Aaron"),
 "app_base": (0.75 / 0.09, "M94", "C", "base temperature of node appearance (beta0 = -beta1 x 0.75/0.09)", "Marcelis 1994"),
 "dev_base": (10.0, "M94", "C", "base temperature of fruit development and of the sink amplitude", "Marcelis 1994"),
 "par_ref": (6.8, "M94", "MJ m-2 d-1", "PAR at which M94 fitted the appearance rate", "Marcelis 1994"),
 "dmc_stem": (0.06, "series", "-", "stem dry-matter content (standing fresh mass)", "Note 91 (assumed)"),
 "dmc_root": (0.07, "series", "-", "root dry-matter content (standing fresh mass)", "Note 91 (assumed)"),
}
outc = {}
for k, pr in CROP_PRIORS.items():
    grp, un, ds = CD[k] if k in CD else ("noise", "-", "likelihood noise")
    outc[k] = {"val": pr[1], "units": un, "group": grp, "status": "inferred", "prior": pr, "desc": ds, "source": "n94/n105 prior (Notes 74, 94, 98)"}
for k, (val, grp, un, ds, so) in CF.items():
    outc[k] = {"val": val, "units": un, "group": grp, "status": "fixed", "desc": ds, "source": so}
json.dump(outc, open(os.path.join(OUT, "crop_n105.json"), "w"), indent=1)
print(len(out), "climate +", len(outc), "crop parameters written")
