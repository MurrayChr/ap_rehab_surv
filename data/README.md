# Dataframes of encounters or penguin-level information

`encounter_wca_old.RDS` and `penguins_old.RDS` were the original data that the 
analysis worked with, which included encounters until some time in 2023.

`encounters_wca.RDS` and `penguins.RDS` are the updated versions complete until
end of 2024. These tidied, formatted files are created in the repo 
`ap_database_wrangling`, where they are in the folder `data/2025_04_24`.

`00b_cmr_data_multisite_multiage.RDS` and `00b_cmr_data_multisite_multiage_trapdep.RDS` are mark-recapture datasets derived from the raw encounter data in the script `R\00b_create_cmr_data.R`. They differ in the state encoding, the former excluding trap-dependent states which the latter includes.
