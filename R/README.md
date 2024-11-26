### R files to perform analyses

#### File naming conventions
Files that contribute to a common task share a numerical prefix, e.g. "00" files are all used to create mark-recapture datasets from the raw encounter data.

#### Functions
- **00_function_get_expected_frequencies** calculates expected cell frequencies for a contingency table under the null hypotheses of no difference between the rows (used for tests 3G.SR and WBWA)
- **00_function_get_ft** calculates Freeman-Tukey discrepancy measure between observed and expected frequencies
- **00_function_get_marray.R** converts suitably encoded capture history matrix into m-array representation (handles single and multistate cases)
- **00_functions_get_gof_tables.R** construct component contingency tables for tests 3G.SR (transience), M.ITEC (trap-dependence) and WBWA (memory) from the capture histories (directly, not via an m-array).

#### Files to prepare mark-recapture datasets
- **00a_multisite_birdyears.R** finds all years in which a bird was encountered at more than one site ('multi-site bird-years'), and assigns a single site.
- **00b_create_cmr_data.R** creates a mark-recapture dataset encoded for a multi-age, multi-site model incorporating trap-dependence in the adult states.

#### Models without a hand-rearing covariate

- **01a_fit_multiage_multisite_trap-dep.R** fits a multi-age, multi-site model incorporating trap-dependence (but not transience) in the 
adult states.

- **01b_fn_sim_rep_data_multiage_multisite_trap-dep.R** contains three functions: `sim_cmr_data` simulates data under model `01_multiage_multisite_trap-dep.stan` given parameter values and a dataset structure, `get_pars_from_posterior` extracts and formats a posterior sample from a fit of this model, and `sim_rep_cmr_data` calls the first two functions to simulate a replicate dataset under the posterior predictive distribution of a particular fit of this model.

-  **02a_fit_multiage_multisite_trap-dep_trans.R** fits a multi-age, multi-site model incorporating trap-dependence and transience in the 
adult states.
