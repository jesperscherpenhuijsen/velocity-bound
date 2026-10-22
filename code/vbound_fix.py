import os
RESULTS = os.environ.get('FTL_RESULTS', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results'))
os.makedirs(RESULTS, exist_ok=True)
import json, sys
import vbound as V
name = sys.argv[1]
cfg = {'R1/R2=0.25': ('R1/R2=0.25', 5.0, 20.0, 1/3, 4000, 0.0, True),
       'R1/R2=0.75': ('R1/R2=0.75', 15.0, 20.0, 1/3, 4000, 0.0, False),
       'sF=2000': ('sF=2000 (sharper shift)', 10.0, 20.0, 1/3, 2000, 0.0, False),
       'sF=8000': ('sF=8000 (smoother shift)', 10.0, 20.0, 1/3, 8000, 0.0, False)}[name]
res = json.load(open(RESULTS + '/vbound_results.json'))
new = V.analyse_config(*cfg)
for i, r in enumerate(res):
    if r['name'] == new['name']:
        if not cfg[6]:
            new['v_NEC_bisection'] = r['v_NEC_bisection']   # keep the earlier valid bisection
        res[i] = new
json.dump(res, open(RESULTS + '/vbound_results.json','w'), indent=1)
print(name, 'v1*=%.4f at r=%.1f, vpar=%.4f, bisect=%s' % (new['v_first_order'], new['v_first_order_location']['r'], new['v_first_order_e_parallel_J'], new['v_NEC_bisection']))
