## R files to perform analyses

### File naming conventions
Files that contribute to a common task share a numerical prefix, e.g. "00" files are all used to create mark-recapture datasets from the raw encounter data.

- **00a_multisite_birdyears.R** finds all years in which a bird was encountered at more than one site ('multi-site bird-years'), and assigns a single site.
- **00b_create_cmr_data.R** creates a mark-recapture dataset encoded for a multi-age, multi-site model incorporating trap-dependence in the adult states.
