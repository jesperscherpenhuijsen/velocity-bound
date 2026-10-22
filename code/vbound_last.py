import os
RESULTS = os.environ.get('FTL_RESULTS', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results'))
os.makedirs(RESULTS, exist_ok=True)
import json, sys
import vbound as V
res = json.load(open(RESULTS + '/vbound_results.json'))
res.append(V.analyse_config('sF=8000 (smoother shift)', 10.0, 20.0, 1/3, 8000))
json.dump(res, open(RESULTS + '/vbound_results.json','w'), indent=1)
