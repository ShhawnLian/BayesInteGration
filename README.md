# BayesInteGration
An illustrative example for our paper "*Bayesian integrative detection of structural variations with FDR control*"

#### Step 1: Obtain the SV calls from individual tools 
  We require the variant call format (VCF) files from different detection tools and merged them into `Data/merged.csv` files.
  
  | Tool | Available Quality Scores | Record in VCF files |
  |----------|----------|----------|
  | [cuteSV](https://github.com/tjiangHIT/cuteSV) | None | None |
  | [pbsv](https://github.com/PacificBiosciences/pbsv) | None | None |
  | [Sniffles](https://github.com/fritzsedlazeck/Sniffles) | Mean mapping quality of supporting reads | QUAL |
  | [DeBreak](https://github.com/Maggi-Chen/DeBreak) | Mean mapping quality of supporting reads | INFO.MAPQ |
  | [SVIM](https://github.com/eldariont/svim) | Scores that accounts for supporting read count, span deviation, and position deviation | QUAL |
  | ... | ... | ... |

#### Step 2: Index the SVs

  We employ a tool-aware merging procedure by running the python file `index.py` to generate the indexed matrix $\boldsymbol{Y}$ and the score matrix $\boldsymbol{S}$.

#### Step 3: Bayesian integration modeling

  Users defined values for hyperparameters, boundaries, initial values, and so on. Then run the R file `BayesianModel.R`

#### Step 4: FDR control

  After the MCMC step get converged, output the results at any specific thresholds. 
  
  Default setting in `BayesianModel.R` will output three levels: Model-0.950, Model-0.990, and Model-0.999, with FDR at 0.05, 0.01, and 0.001, respectively.
