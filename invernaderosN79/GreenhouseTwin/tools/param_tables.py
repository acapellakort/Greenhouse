"""
param_tables.py -- LaTeX parameter tables of the n79 model, generated from the JSON files (no hand transcription).
Usage: python3 tools/param_tables.py > ../../docs/n79_param_tables.tex   (from the GreenhouseTwin folder)
"""
import json, re

def tex(s):
    s = str(s)
    for a, b in (("\\", r"\textbackslash{}"), ("_", r"\_"), ("%", r"\%"), ("&", r"\&"), ("#", r"\#"), ("~", r"\~{}")):
        s = s.replace(a, b)
    s = re.sub(r"\bm(\d)", r"m$^{\1}$", s)
    s = s.replace("-1", "$^{-1}$").replace("-2", "$^{-2}$").replace("-3", "$^{-3}$").replace("-4", "$^{-4}$")
    return s

def num(v):
    if v == 0:
        return "0"
    a = abs(v)
    if 1e-3 <= a < 1e5:
        s = f"{v:.6g}"
    else:
        m, e = f"{v:.4e}".split("e")
        m = m.rstrip("0").rstrip(".")
        s = f"${m}\\times10^{{{int(e)}}}$"
    return s

ORDER_C = ["geometry", "constants", "air", "canopy", "photosynthesis", "floor", "mass", "cover", "screens", "heating",
           "ventilation", "lamps", "actuators", "noise"]
ORDER_K = ["M94", "canopy", "architecture", "respiration", "fruit set", "quality", "harvest", "photosynthesis", "light",
           "series", "noise"]

def table(path, order, caption, label):
    d = json.load(open(path))
    out = [r"{\small", r"\begin{longtable}{>{\raggedright\arraybackslash}p{0.17\textwidth}>{\raggedright\arraybackslash}p{0.17\textwidth}"
           r">{\raggedright\arraybackslash}p{0.13\textwidth}>{\raggedright\arraybackslash}p{0.43\textwidth}}",
           r"\caption{" + caption + r"}\label{" + label + r"}\\", r"\toprule",
           r"name & value (prior) & units & meaning; source\\", r"\midrule\endfirsthead",
           r"\toprule name & value (prior) & units & meaning; source\\ \midrule\endhead"]
    groups = sorted(set(v["group"] for v in d.values()), key=lambda g: order.index(g) if g in order else 99)
    for g in groups:
        out.append(r"\multicolumn{4}{l}{\textit{" + tex(g) + r"}}\\")
        for k in sorted(k for k, v in d.items() if v["group"] == g):
            v = d[k]
            val = r"\textbf{inferred} [" + ", ".join(num(x) for x in v["prior"]) + "]" if v["status"] == "inferred" else num(v["val"])
            out.append(r"\texttt{" + tex(k) + "} & " + val + " & " + tex(v["units"]) + " & " + tex(v["desc"]) + "; "
                       + tex(v["source"]) + r"\\")
    out += [r"\bottomrule", r"\end{longtable}}"]
    return "\n".join(out)

print(table("configfiles/climate_n79.json", ORDER_C,
            "Climate and photosynthesis parameters (\\texttt{configfiles/climate\\_n79.json}). Inferred parameters: "
            "[min, centre, max] of the uniform prior; the centre is the default value.", "tab:pclim"))
print()
print(table("configfiles/crop_n105.json", ORDER_K,
            "Crop parameters (\\texttt{configfiles/crop\\_n105.json}).", "tab:pcrop"))
