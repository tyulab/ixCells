import numpy as np

TIERS_THRESHOLD = 1.5
Z_SCORE_THRESHOLD = 1.8
naive = 'Naive'
neg_control = '(Ionis676630|Ionis 676630).*_10'
RENORMALIZE = False
ROUND = 'Round 2.3'

# bins to check intervals on. nomenclature: "10" means 0.1 <= x < 0.2, "90" means 0.9 <= x < inf, etc...
BINS = np.linspace(0.0,1.0,11) # right = False
names = list(map(lambda x: '0' + str(x),np.arange(0,100,10)))
names[0] = '000'
names.append('100')
# the name assigned to each bin by corresponding indices
NAMES = names
# key starting from 1, names at NAMES
DICT = dict(enumerate(NAMES, 1))