"""
audit.py -- static checks of GreenhouseTwin (no Julia needed)
  1. no V2.2 flag name or `hasproperty(p, ...)` switch is left in src/
  2. every parameter the climate code reads (p.<name> or destructured from p) is in climate_n79.json, and every
     climate_n79.json key is read by the code (noise terms by the calibration)
  3. same for the crop (gp.<name>) and crop_n105.json
Usage: python3 tools/audit.py   (from the GreenhouseTwin folder)
"""
import json, re, glob, os, sys

src = {f: open(f).read() for f in glob.glob("src/*.jl")}
code = "\n".join(src.values())
FLAGS = """use_cover_node heat_ctrl vent_screen stomata_fix vapour_fix cond_heat vent_ctrl pipe_lag black_scr vc_night_log
vc_min_night hc_log_pipe vc_min_log soil_dyn floor_full solar_split deep_k mass_node lamp_full vent_split lamp_scr_loss
s_rad dday phys_light photo_lai cover_cond co2_resp dew_bl leak_shared rd_q10 leak_stack crop_cap scr_solar solar_cons
solar_mass_share par_below gen_canopy gen_parlight gen_zclamp gen_clampT gen_beta gen_kappa gen_A gen_V gen_node
nu4_min U8_min_leak vent_exp scr_order roof_rho_par roof_rho_nir roof_tau_nir stomata_model""".split()
ok = True
def strip_comments(s):
    out = []
    for line in s.splitlines():
        q = line.split("#", 1)[0]
        out.append(q)
    return "\n".join(out)
codenc = "\n".join(strip_comments(s) for s in src.values())
found = [f for f in FLAGS if re.search(r"\b" + f + r"\b", codenc)]
print("1. flag names in src/ (outside comments):", found or "none")
hp = re.findall(r"hasproperty\(\s*p\b", codenc)
print("   hasproperty(p, ...) switches:", len(hp))
ok &= not found and not hp

def used(prefix, text):
    names = set(re.findall(r"\b" + prefix + r"\.([A-Za-z_][A-Za-z0-9_]*)", text))
    for blk in re.findall(r"\(;\s*([^)]*?)\)\s*=\s*" + prefix + r"\b", text, re.S):
        names |= set(x.strip() for x in blk.replace("\n", " ").split(",") if x.strip())
    return names

clim = json.load(open("configfiles/climate_n79.json"))
crop = json.load(open("configfiles/crop_n105.json"))
climate_files = ["src/" + f for f in ("solar.jl", "photosynthesis.jl", "climate.jl", "crop_measured.jl")]
pc = set()
for f in climate_files:
    pc |= used("p", strip_comments(src[f]))
pc -= {"t"}
missing = sorted(pc - set(clim))
noise = {"sigma_T", "sigma_RH", "sigma_C"}
unused = sorted(set(clim) - pc - noise)
print("2. climate parameters read by the code:", len(pc), "| JSON keys:", len(clim))
print("   read but not in JSON:", missing or "none")
print("   in JSON but never read:", unused or "none")
ok &= not missing and not unused

gc = set()
for f in ("src/crop.jl", "src/crop_measured.jl"):
    gc |= used("gp", strip_comments(src[f]))
cropnoise = {"sigma_y", "sigma_leaf", "sigma_fd", "sigma_n", "J_max", "W_pot", "E_fruit", "app_beta1", "V20",
             "transmission", "lamp_umol"}
missing = sorted(gc - set(crop))
unused = sorted(set(crop) - gc - cropnoise)
print("3. crop parameters read by the code:", len(gc), "| JSON keys:", len(crop))
print("   read but not in JSON:", missing or "none")
print("   in JSON but never read (beyond inferred-by-name and light-formula keys):", unused or "none")
ok &= not missing and not unused
print("\nAUDIT", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
