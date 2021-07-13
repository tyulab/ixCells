#!/usr/bin/env python
# Loads csvs, gets zscore, average zscore by sample and plate, and puts to output csv

# set outliers above a certain z score threshold to na (preserve shape)
from run import create_output

if __name__ == "__main__":
    create_output()
