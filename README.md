# recycle_model
 model of phototroph/heterotroph interaction

This repository contains SW and notebooks used in the creation of:


**Four Potential Mechanisms Underlying Phytoplankton-Bacteria Interactions Assessed Using Experiments and Models**

Osnat Weissberg, Dikla Aharonovich, Zhen Wu, Michael J. Follows, Daniel Sher

An ecological model of the interactions between *Prochlorococcus* (cyanobacteria) and marine heterotrophs. 
Prochlorococcus is marine autotroph, ubiquitous in the world oligotrophic oceans. When grown in coculture with heterotrophic bacteria, different outcomes are observed, both positive and negative.
The goal of this repo is to model this interaction, explore possible mechanism underlying the observed outcomes, associated organism traits and possible biogeochemical implications.
The code contains an implementation of the model, code for running simulations using randomly selected parameter values, a classification engine to classify to the different experimental outcomes and jupyter notebooks analyzing the results.



# System Requirments
* Python 3.11
* Additional required python packages are listed in requirements.txt
The code has been tested on both windows and Linux.

# Installation

Clone the github repo.
Install additional required python packages: 
 pip install -r requirements.txt

Compilation is not required. Python is used for running the code in the repo.
Typical install time is less than a minute.

# Running the code
The published work is contained in the folder STORE_model.
```
 cd STORE_model
```

The basic script for run the models, randomly select parameter values, simulate the model, compare to experimental reference results and output simulation results, parameter values and comparison to experimental results.
```
model_equations_separate_NC_store_numba.py [-h] [--ref_csv REF_CSV]
                                                  [--ref_pro99_csv REF_PRO99_CSV]
                                                  [--json JSON [JSON ...]]
                                                  [--jsondpath JSONDPATH]
                                                  [--maxday MAXDAY]
                                                  [--outdpath OUTDPATH]
                                                  --run_id RUN_ID --model
                                                  {MIN,MIXOTROPH,OVERFLOW,ROS,EXOENZYME,ROS-MIXOTROPH-OVERFLOW-EXOENZYME,ROS-MIXOTROPH,EXOENZYME-ROS,EXOENZYME-MIXOTROPH,OVERFLOW-MIXOTROPH,OVERFLOW-EXOENZYME,OVERFLOW-ROS,EXOENZYME-ROS-MIXOTROPH,OVERFLOW-EXOENZYME-MIXOTROPH,OVERFLOW-ROS-MIXOTROPH,OVERFLOW-ROS-EXOENZYME}
                                                  [--which_organism {ponly,honly,all}]
                                                  [--pro99_mode]
                                                  [--t_eval T_EVAL [T_EVAL ...]]
                                                  [--param_sensitivity PARAM_SENSITIVITY]
                                                  [--organism_to_tune {PRO,HET}]
                                                  [--number_of_runs NUMBER_OF_RUNS]
                                                  [--sobol SOBOL]
                                                  [--monte_max_params MONTE_MAX_PARAMS]
                                                  [--monte_add_noise]
                                                  [--monte]
Run models - nutrients recycle with separate N/C and quotas.

options:
  -h, --help            show this help message and exit
  --ref_csv REF_CSV     reference CSV
  --ref_pro99_csv REF_PRO99_CSV
                        reference pro99 CSV
  --json JSON [JSON ...]
                        json with param vals
  --jsondpath JSONDPATH
                        directory with json files
  --maxday MAXDAY       max day of simulation
  --outdpath OUTDPATH   output dir
  --run_id RUN_ID       run id
  --model {MIN,MIXOTROPH,OVERFLOW,ROS,EXOENZYME,ROS-MIXOTROPH-OVERFLOW-EXOENZYME,ROS-MIXOTROPH,EXOENZYME-ROS,EXOENZYME-MIXOTROPH,OVERFLOW-MIXOTROPH,OVERFLOW-EXOENZYME,OVERFLOW-ROS,EXOENZYME-ROS-MIXOTROPH,OVERFLOW-EXOENZYME-MIXOTROPH,OVERFLOW-ROS-MIXOTROPH,OVERFLOW-ROS-EXOENZYME}
                        model to run
  --which_organism {ponly,honly,all}
                        which organism to run
  --pro99_mode          run on pro99 media
  --t_eval T_EVAL [T_EVAL ...]
  --param_sensitivity PARAM_SENSITIVITY
                        index of param to update (0 based)
  --organism_to_tune {PRO,HET}
                        which organism to tune
  --number_of_runs NUMBER_OF_RUNS
                        number of simulations to run
  --sobol SOBOL         run sobol (will run 2^<sobol> simulation
  --monte_max_params MONTE_MAX_PARAMS
                        max number of params to update per monte simulations
  --monte_add_noise     add noise to monte fixed params
  --monte               run monte carlo
```


Post process the results of model simulations. The simulator creates a multiple files per each simulation. This scripts concatenate these into one large file.
```
results_postprocess.py [-h] --dpath DPATH [DPATH ...]
                              [--outdpath OUTDPATH] --run_id RUN_ID
                              [--ref_csv REF_CSV]

post process results

options:
  -h, --help            show this help message and exit
  --dpath DPATH [DPATH ...]
                        paths to load from
  --outdpath OUTDPATH   output dir
  --run_id RUN_ID       run id
  --ref_csv REF_CSV     reference CSV
```


Analyze Prochlorococcus monoculture simulations. This script cleans up problematic simulations (i.e. those where the simulator ran into mathematical difficulties or where the biomass becomse negative (very small % of simulations), compares the simulations results to the experimental monocultures using 2 different media (each randomly generated virtual Prochlorococcus - vPro - is simulated on two different media, and both simulated are compared to the experimental growth. Select vPros that have a low root mean square error (RMSE) vs the lab on both media. 

```
run_ponly_clean_find_vpros.py [-h] [--outdpath OUTDPATH]
                                     [--vprooutdpath VPROOUTDPATH]
                                     [--extendvpro] [--prefixvpro PREFIXVPRO]
                                     [--indpath INDPATH] --run_id RUN_ID

cleanup PONLY runs and produce VPROS.

options:
  -h, --help            show this help message and exit
  --outdpath OUTDPATH   output dir
  --vprooutdpath VPROOUTDPATH
                        output dir
  --extendvpro          this is a run extending existing vpros
  --prefixvpro PREFIXVPRO
                        prefix to remove from idx when constracting VPRO name
  --indpath INDPATH     input dir
  --run_id RUN_ID       run id
```

Once the coculture simulations are done and postprocessed, this script cleans up problematic simulations (i.e. those where the simulator ran into mathematical difficulties or where the biomass becomse negative (very small % of simulations), classifies each simulation into one of the coculture outcomes, and identifies versatile vPros, vPros with both positive and negative outcomes.
```
run_het_clean_classify_veratiles.py [-h] [--outdpath OUTDPATH]
                                           [--indpath INDPATH] --run_id RUN_ID

cleanup coculture runs and produce VPROS.

options:
  -h, --help           show this help message and exit
  --outdpath OUTDPATH  output dir
  --indpath INDPATH    input dir
  --run_id RUN_ID      run id
```


Go over the classification results and select the simulations with the highest classification probability, per model and per outcome.

```
run_get_top_df.py [-h] [--outdpath OUTDPATH] --topids TOPIDS
                         [--infpath INFPATH]
create df with only the top runs

options:
  -h, --help           show this help message and exit
  --outdpath OUTDPATH  output dir
  --topids TOPIDS      csv file with a list of top runs
  --infpath INFPATH    input dir
```

Go over simulation results and integrate to estimate the total fluxes. Also select the final simulation outcomes.
```
run_het_biogeo_create_files.py [-h] [--outdpath OUTDPATH]
                                      [--indpath INDPATH] --run_id RUN_ID

Integrate to identify biogeo impacts.

options:
  -h, --help           show this help message and exit
  --outdpath OUTDPATH  output dir
  --indpath INDPATH    input dir
  --run_id RUN_ID      run id
```

Identify the portion of the simulations were growth is balanced.
``` 
run_het_stage_df.py [-h] [--outdpath OUTDPATH] [--infpath INFPATH]

Identify the portion (stage) of the simulations were growth is balanced.

options:
  -h, --help           show this help message and exit
  --outdpath OUTDPATH  output dir
  --infpath INFPATH    input dir
```

# Notebooks used to analyze the simulation results

* STORE_model/notebooks/CC10_phylogenetic_tree.ipynb : create phylogenetic tree of the heterotroph bacteria used in the experiment (figure 1). 
* STORE_model/notebooks/model_store_analyze_monte_ponly_curves.ipynb - analyze monoculuture simulations (Figure 1)
* STORE_model/notebooks/model_store_analyze_monte_curves.ipynb - analyze coculuture simulations (Figure 1)
* STORE_model/notebooks/model_store_analyze_monte_ponly_params.ipynb - analyze vPro parameters in monoculuture simulations (Figure S5, S6)
* STORE_model/notebooks/verastile_vpros.ipynb - distribution of outcomes and versatility (Figure 2, S12)
* STORE_model/notebooks/verastile_vpros_multi.ipynb - distribution of outcomes and versatility for mechanism combinations (Figure S15, S16)
* STORE_model/notebooks/model_store_analyze_monte_correlation_to_biomass.ipynb - correlation between param values and biomass (Figure 3)
* STORE_model/notebooks/model_store_analyze_monte_curves_top.ipynb - analyze top select coculuture simulations (Figure 4)
* STORE_model/notebooks/model_store_analyze_monte_biogeo.ipynb - analyze biogeo results (figure 5)
* STORE_model/notebooks/model_store_analyze_monte_biogeo_multi.ipynb - analyze biogeo results for mechanism combinations
* STORE_model/notebooks/create param table.ipynb : Create the table of parameters model (supplemantary table S2)
* STORE_model/notebooks/model_store_analyze_monte_compare_MSE_prediction.ipynb - compare two alternative classification methods (ML and RMSE) (Figure S10)



  

