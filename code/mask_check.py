import json, sys
sys.path.insert(0,'/home/claude/ftl_check')
import vbound as V
cfg = {'pub': ('published', 10.0, 20.0, 1/3, 4000), 'c07': ('C=0.70', 10.0, 20.0, 0.70, 4000)}[sys.argv[1]]
V.MASK = float(sys.argv[2])
r = V.analyse_config(*cfg, vlo=0.0, vhi=0.08)
json.dump(r, open('/home/claude/ftl_check/mask_%s_%s.json' % (sys.argv[1], sys.argv[2]), 'w'))
