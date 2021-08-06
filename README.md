## iXCells scripts

store script(s) for figures and tables from iXCells data

### How to use
#### Input files required
- `/Round#` (contain all the data and output for the round)
  - `/data` (contain all the data to be searched/used as input)

Example folder structure within repo directory:
```shell
C:.
└───Round 2.3
    └───data
        ├───R2.3-1,2 summary 07262021
        ├───R2.3-3,4 summary 07262021
        ├───R2.3-5,6 summary 07272021
        ├───R2.3-7,8 summary 07272021
        └───R2.3-9 summary 07272021
```
#### Running the script
1. Check/install requirements: `pip install -r requirements.txt`
1. `config.py` open to change options for running figures
   - Must specify folder with round data (e.g. Round 1, Round 2...)
1. `python run.py` to run all figures. To create select figures comment out the irrelevant functions.
    - `run_renata_data.py` run figures on Renata's sheets and similar formatted files
#### Scripts
- `tasks.py`: store functions that load csvs and create individual tables/figures
- `/utils`
  - `plots.py`: contains functions to generate/save figures
  - `stats.py`: functions that perform mathematical operations
  - `functions.py`: misc helper functions for dataframe manipulation
- `normalized_naive` or `normalized_neg`: folders separated by normalization containing output data
  - `Control plots`: store plots and figures

### Todo:
- [x] Plots to pdf
  - [x] config setting for pdf vs png
  - [x] control plates histograms
  - [x] r squared
- [ ] Add col for biological replicates B1/B2
  - [ ] rewrite all functions that use separate bio replicates to use new column
    - - [ ] Avg Exp in avg_exp
    - - [ ] ctrl+f for // 3
    - - [ ] r^2 plots
  - [ ] fix functions to reformat renata's data to work on same tasks (avg exp, output excluded)
- [ ] Round prefixes to plot titles
- [ ] check changes on older rounds, error handling 
  
#### Future improvements
- [ ] Rework some variables into dictionaries
  - [ ] config options
  - [ ] column names
- [ ] Clarify z score column naming
- [ ] Script or function for data validiation?
  
Other possible todos in scripts

### Notes:
- Creating initial output requires correct column names (or some functions to make them correct)
  - Necessary columns:
    - Experiment Name
      - Type (MT, WT, Total) is assigned from searching values
    - SampleName
      - used with Experiment Name in various groupbys
    - ASO Microsynth ID (or ASO Short Description)
        - plate # is found by searching the last number from values
    - CrossingPoint
    - CrossingPoint.1
      - crossing point order doesn't matter since it is changed to abs value
    - dCt (Round 1)
- Check letter case for type/control/naive search errors