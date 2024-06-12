# BayesInteGration
Illustrative example for our paper "*Bayesian integrative detection of structural variations with FDR control*"

#### Step 1: Obtain the SV calls from individual tools 

  | Tool | Available Quality Scores | Column 3 |
  |----------|----------|----------|
  | [cuteSV](https://github.com/tjiangHIT/cuteSV) | None | Cell 3   |
  | [pbsv](https://github.com/PacificBiosciences/pbsv) | None | Cell 6   |
  | [Sniffles](https://github.com/fritzsedlazeck/Sniffles) | Mean mapping quality of supporting reads | Cell 9   |
  | [DeBreak](https://github.com/Maggi-Chen/DeBreak) | Mean mapping quality of supporting reads | Cell 9   |
  | [SVIM](https://github.com/eldariont/svim) | Scores that accounts for supporting read count, span deviation, and position
deviation | Cell 9   |
  | ... | ... | ... |


#### Step 2: Tool-aware merging

  We required the VCF files from different detection tools and merged them into 'Data/merged.csv' files.

  Then run python file `index.py` to generate the indexed matrix $\boldsymbol{Y}$ and the score matrix $\boldsymbol{S}$.

#### Step 2: Bayesian integration modeling

  Users defined values for hyperparameters, boundaries, initial values, and so on. Then run the R file `BayesianModel.R`

#### Step 3: FDR control

  After the MCMC step get converged, output the results at any specific thresholds. 
  
  Default setting in `BayesianModel.R` will output three levels: Model-0.950, Model-0.990, and Model-0.999, with FDR at 0.05, 0.01, and 0.001, respectively.
