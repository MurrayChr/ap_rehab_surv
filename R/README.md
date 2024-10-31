### R files to perform analyses

#### File naming conventions
Files that contribute to a common task share a numerical prefix, e.g. "00" files are all used to create mark-recapture datasets from the raw encounter data.

#### Files to prepare mark-recapture datasets
-**00_function_get_marray.R** converts suitably encoded capture history matrix into m-array representation (handles single and multistate cases)
- **00a_multisite_birdyears.R** finds all years in which a bird was encountered at more than one site ('multi-site bird-years'), and assigns a single site.
- **00b_create_cmr_data.R** creates a mark-recapture dataset encoded for a multi-age, multi-site model incorporating trap-dependence in the adult states.

#### Models without a hand-rearing covariate

- **01a_fit_multiage_multisite_trap-dep.R** fits a multi-age, multi-site model incorporating trap-dependence (but not transience) in the 
adult states.
