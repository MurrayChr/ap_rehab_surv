### Stan models

This folder contains Stan code for all the mark-recapture models, plus one model `00_multinomial_mixture.stan` used in the goodness-of-fit test of trap-dependence (M.ITEC). The file prefices correspond to those of the R scripts where the mark-recapture models are fit:

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

Model 04 has two Stan files, `04_td_tr_hr.stan` using the EPM likelihood and `04_td_tr_hr_hhmm.stan` using the HHMM likelihood.


