# Data

### Raw data on encounters and penguins

- `encounters_wca.RDS` contains raw encounter data for the region studied ("west of Cape Agulhas")
- `penguins.RDS` contains a variable indicating whether a penguin was hand-reared or not

These are the data used in the paper, complete until end of 2024. These tidied, formatted files are created in a separate repo 
`ap_database_wrangling`, where they are in the folder `data/2025_04_24`. Older versions of this data, complete up until some time in 2023, are in `encounter_wca_old.RDS` and `penguins_old.RDS`.


### Derived capture-mark-recapture datasets

`00b_cmr_data_multisite_multiage.RDS` and `00b_cmr_data_multisite_multiage_trapdep.RDS` are mark-recapture datasets derived from the raw encounter data in the script `R\00b_create_cmr_data.R`. They differ in the state encoding, the former excluding trap-dependent states which the latter includes.
