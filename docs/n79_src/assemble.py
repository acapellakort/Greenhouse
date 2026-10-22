# assemble.py -- builds greenhouse_model_n79.tex from template_n79.tex, parts/*.tex and the generated parameter tables
import subprocess
s = open('template_n79.tex').read()
tables = subprocess.run(['python3', 'tools/param_tables.py'], cwd='../../invernaderosN79/GreenhouseTwin', capture_output=True, text=True, check=True).stdout
rep = {'CROPSECTION': open('parts/crop.tex').read(), 'PARAMTABLES': tables, 'INFERENCE': open('parts/inference.tex').read(),
       'GAPS': open('parts/gaps.tex').read(), 'VERIFY': open('parts/verify.tex').read()}
for k, v in rep.items():
    assert s.count(k) == 1, k
    s = s.replace(k, v)
open('../greenhouse_model_n79.tex', 'w').write(s)
print('wrote greenhouse_model_n79.tex')
