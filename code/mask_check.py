import os
RESULTS = os.environ.get('FTL_RESULTS', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results'))
os.makedirs(RESULTS, exist_ok=True)
import json, sys
import vbound as V
cfg = {'pub': ('published', 10.0, 20.0, 1/3, 4000), 'c07': ('C=0.70', 10.0, 20.0, 0.70, 4000)}[sys.argv[1]]
V.MASK = float(sys.argv[2])
r = V.analyse_config(*cfg, vlo=0.0, vhi=0.08)
json.dump(r, open(RESULTS + '/mask_%s_%s.json' % (sys.argv[1], sys.argv[2]), 'w'))
