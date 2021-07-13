#!/usr/bin/env python

# Exp zscore plots and dose response step- flag ASO if avg 10mmol exp value > avg 3mmol exp value

# Calculate mean exp for each biological replicate
# TODO: take out naive/controls
from run import create_avg_exp, create_tiers

if __name__ == "__main__":
    create_avg_exp()
    create_tiers()