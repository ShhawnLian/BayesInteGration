# BayesInteGration
Example for our paper "*Bayesian integrative detection of structural variations with FDR control*"

#### Step 1: Tool-aware merging

  We required the VCF files from different detection tools and merged them into 'Data/merged.csv' files.

  Then run python file `index.py` to generated the indexed matrix $\boldsymbol{Y}$ and the score matrix $\boldsymbol{S}$.

#### Step 2: Bayesian integration modeling

  Users defined values for priors, boundaries, initial values, and so on. Then run the Bayesian model via the R file `BayesianModel.R`

#### Step 3: FDR control

  After the MCMC step get converged, output the results at any specific thresholds. 
  
  Default setting in `BayesianModel.R` will output three levels: Model-0.950, Model-0.990, and Model-0.999, with FDR at 0.05, 0.01, and 0.001, respectively.
