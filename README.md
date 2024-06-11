# BayesInteGration
Example for our paper "*Bayesian integrative detection of structural variations with FDR control*"

#### Step 1: Tool-aware merging

  We required the VCF files from different detection tools and merged them into 'Data/merged.csv' files.

  Then run `python index.py` to generated the indexed matrix $\boldsymbol{Y}$ and the score matrix $\boldsymbol{S}$.

#### Step 2: Bayesian integration modeling

  User defined values for priors, boundaries, initial values, and so on, and run the Bayesian model via the R programming `BayesianModel.R`

#### Step 3: FDR control

  After the MCMC step get converged, output the results given any thresholds. Default setting in `BayesianModel.R`.
