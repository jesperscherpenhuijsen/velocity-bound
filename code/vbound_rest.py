import json, sys
sys.path.insert(0,'/home/claude/ftl_check')
import vbound as V
configs = [
    ('C=0.50', 10.0, 20.0, 0.50, 4000),
    ('C=0.70', 10.0, 20.0, 0.70, 4000),
    ('R1/R2=0.25', 5.0, 20.0, 1/3, 4000),
    ('R1/R2=0.75', 15.0, 20.0, 1/3, 4000),
    ('sF=2000 (sharper shift)', 10.0, 20.0, 1/3, 2000),
    ('sF=8000 (smoother shift)', 10.0, 20.0, 1/3, 8000),
]
res = json.load(open('/home/claude/ftl_check/vbound_results_part1.json'))
for c in configs:
    res.append(V.analyse_config(*c))
    json.dump(res, open('/home/claude/ftl_check/vbound_results.json','w'), indent=1)
print('ALL DONE')
