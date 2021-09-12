import numpy as np

# Folder for Round with data (also used in plot titles)
ROUND = 'Round 3.1'
# True if normalized to negative
RENORMALIZE = True


# Initial Z score threshold (remove entries in output csv and send to output_dropped)
Z_SCORE_THRESHOLD = 1.8
# Z score range threshold used when creating tiers (send to tiers_dropped)
TIERS_THRESHOLD = 1.5

# Show plots while running
SHOW_PLOTS = False
# naive = 'Naive'
# neg_control = '(Ionis676630|Ionis 676630).*_10'

# Bins/names for Tiers (these can also be entered by hand)
# list of bins for tier intervals. nomenclature: "10" means 0.1 <= x < 0.2, "90" means 0.9 <= x < inf, etc...
BINS = np.linspace(0.0,1.0,11) # right = False
names = list(map(lambda x: '0' + str(x),np.arange(0,100,10)))
names[0] = '000'
names.append('100')
# the name assigned to each bin by corresponding indices
NAMES = names
# create key starting from 1, names at NAMES
DICT = dict(enumerate(NAMES, 1))

# Change file type for certain figures
SAVE_FIG_AS = 'pdf' # png, pdf
