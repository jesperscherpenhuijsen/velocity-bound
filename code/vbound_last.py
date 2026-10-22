import json, sys
sys.path.insert(0,'/home/claude/ftl_check')
import vbound as V
res = json.load(open('/home/claude/ftl_check/vbound_results.json'))
res.append(V.analyse_config('sF=8000 (smoother shift)', 10.0, 20.0, 1/3, 8000))
json.dump(res, open('/home/claude/ftl_check/vbound_results.json','w'), indent=1)
