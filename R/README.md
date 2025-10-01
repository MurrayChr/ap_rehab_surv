## R files to perform analyses

Files that all relate to the same model share a common numerical prefic from 01 to 08. Files starting with the prefix 00 mostly either contain functions, or perform various data-related tasks. Lastly files with a prefix in the 20's contain code to make various figures.   

### Files with 00 prefix

#### Functions for goodness-of-fit testing
- `00_functions_get_gof_tables.R` constructs component contingency tables for tests 3G.SR (transience), M.ITEC (trap-dependence) and WBWA (memory) from the capture histories (directly, not via an m-array).
- `00_function_get_expected_frequencies` calculates expected cell frequencies for a contingency table under null hypotheses of no difference between the rows (used for tests 3G.SR and WBWA) or a mixture model (used for test M.ITEC)
- `00_function_fit_mixtures` returns maximum likelihood estimates of multinomial cell probabilities and mixture probabilities for test M.ITEC of mixtures
- `00_function_get_ft` calculates Freeman-Tukey discrepancy measure between observed and expected frequencies
 
#### Function to create m-arrays
- `00_function_get_marray.R` converts suitably encoded capture history matrix into m-array representation (handles single and multistate cases)

#### Files to prepare mark-recapture datasets
- `00a_multisite_birdyears.R` finds all years in which a bird was encountered at more than one site ('multi-site bird-years'), and assigns a single site.
- `00b_create_cmr_data.R` creates a mark-recapture dataset encoded for a multi-age, multi-site model incorporating trap-dependence in the adult states.
- `00c_cmr_data_plots_and_summaries.R` various numerical and graphical summaries of the mark-recapture data

#### Other 
- `00_beta_spike_and_slab_engineering.R` simulates data from the "beta-spike-and-slab" prior used on residency probabilities in model 06 (see below)

### Files with prefices 01 to 08

Here is the correspondence between the model names in the manuscript and the file prefices used to analyse them:

| File prefix | Model | Trap-dependence | Transience | Hand-rearing | Immature |  Further description |
| :----: | :-----------: | :-------: | :-------: |:-------: |:-------: | :----------- |
| 01 | $\boldsymbol{\mathcal{M}}$ | no | no | no | no | multi-age, multi-site model |
| 02 | $\boldsymbol{\mathcal{M}}_{\text{td}}$ | yes | no | no | no |  |
| 03 | $\boldsymbol{\mathcal{M}}_{\text{td}+\text{tr}}$ | yes | yes | no | no |  |
| 04 | $\boldsymbol{\mathcal{M}}_{\text{td}+\text{tr}}^{\text{hr}}$ | yes | yes | yes | no |  |
| 05 | $\boldsymbol{\mathcal{M}}_{\text{td}+\text{tr}+\text{imm}}^{\text{hr}}$ | yes | yes | yes | yes | immature survival not assumed to be the same as adult survival|
| 06 | $\boldsymbol{\mathcal{M}}_{\text{td}+\text{tr-bss}+\text{imm}}^{\text{hr}}$ | yes | yes | yes | yes | spike-and-slab prior on residency probabilities |
| 07 | $\boldsymbol{\mathcal{M}}_{\text{Stony}}^{\text{hr}}$ | yes | - | yes | no | multi-age, single-site model fit to a subset of data consisting of birds marked as juveniles at Stony |
| 08 | - | yes | yes| yes | no | estimates survival separately for each combination of site, age class, year and hand-rearing level |

For the first six of these models there are a few files associated with each:
- Models 01 to 03 each have a file to fit the model (e.g. `01a_fit_*`), a file containing functions to simulate replicate data under the model (e.g. `01b_fn_sim_rep_data_*`), and a file to do posterior predictive checking of the model (e.g. `01c_ppc_*`).
- Models 04 and 05 each have a file to fit the model, a file containing functions to simulate data under the model, and a file to evaluate the model using simulation (e.g. `04c_sim_*`). 
- Model 04 has two model fitting files, one using the EPM likelihood and one the HHMM likelihood

For the last two models there are only files to fit the model.




