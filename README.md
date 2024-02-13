# ALPHA

(Data & Codes for:)<br>
**Phenotypic evolution of SARS-CoV-2: a statistical inference approach**

Wakinyan Benhamou, Sébastien Lion, Rémi Choquet and Sylvain Gandon
*CEFE, CNRS, Univ Montpellier, EPHE, IRD, Montpellier, France*<br>

Evolution, Volume 77, Issue 10, October 2023, Pages 2213–2223, [https://doi.org/10.1093/evolut/qpad133](https://doi.org/10.1093/evolut/qpad133)

## Data (sources and short descriptions)

This is previously publicly available data that were downloaded from the Internet. All data used in the scripts are in the 'Data' folder (and may also be found from the sources mentioned below). 

### covid-stringency-index.csv

&#11169;&emsp; Stringency Index (from 21/01/2020 to 23/03/2021), computed by the Oxford COVID-19 Government Response Tracker (OxCGRT) [[Hale *et al.*, 2021](https://doi.org/10.1038/s41562-021-01079-8)]. This file was originally downloaded from the website [*Our World in Data*](https://ourworldindata.org/) but does not seem to be available there anymore... However, values from old databases can still be found in a GitHub repositery of OxCGRT: '[covid-policy-tracker-legacy](https://github.com/OxCGRT/covid-policy-tracker-legacy)' (see 'Legacy OxCGRT dataset (retired July 2022)'). Although 184 countries are included in this file, we only used the Stringency Index in the UK (`Entity == "United Kingdom"`).

### daily-tests-and-daily-new-confirmed-covid-cases.csv

&#11169;&emsp; Daily number of tests with new cases tested positive in the UK (from 07/04/2020 to 28/04/2021), downloaded from *Our World in Data*:
https://ourworldindata.org/grapher/daily-tests-and-daily-new-confirmed-covid-cases?country=~GBR

### data_daily_deaths_UK_2021_May_05.csv

&#11169;&emsp; Daily new fatality cases in the UK (from 02/02/2020 to 23/04/2021), downloaded from *GOV.UK* (visited 05/05/2021):
https://coronavirus.data.gov.uk/details/deaths (see 'Daily deaths with COVID-19 on the death certificate by date of death').<br>
<sub>(licenced under [Open Government Licence v3.0](https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/))</sub>

### Variant_of_Concern_VOC_202012_01_Technical_Briefing_5_Data_England.ods

[DOWNLOAD](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/957631/Variant_of_Concern_VOC_202012_01_Technical_Briefing_5_Data_England.ods) (ods file)

&#11169;&emsp; Fequencies of S-Gene Target Failures (SGTF) in the 9 regions of England; underlying data from *Public Health England* [*Technical Briefing 5: Investigation of novel SARS-CoV-2 variant - Variant of Concern 202012/01*](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/959426/Variant_of_Concern_VOC_202012_01_Technical_Briefing_5.pdf).

<sub>All PHE (now  replaced by UKHSA) technical briefings dealing with the investigation of SARS-CoV-2 variants may be found on *GOV.UK*:</sub><br>
<sub>Technical briefings 1 to 23: https://www.gov.uk/government/publications/investigation-of-novel-sars-cov-2-variant-variant-of-concern-20201201</sub><br>
<sub>From technical briefing 24: https://www.gov.uk/government/publications/investigation-of-sars-cov-2-variants-technical-briefings</sub>

## R codes

All the scripts are in the root directory of this repositery. We developed a two-step approach based on an SEIR model and the analysis of a combination of epidemiological and evolutionary information:

(1) In the first step, we infer how the Stringency Index reduces the amount of viral transmission - assuming that it wouldn't affect the duration of infectiousness - just before the emergence of the Alpha Variant;

Scripts:
- **STEP1_functions.R** (script with the custom functions used in the first step)
- **STEP1_Script_UK_data.Rmd** (analyses using the UK data)
- **STEP1_Identifiability_profiles.Rmd** (build and plot identifiability profiles of the model)

(2) In the second step, based on a novel theoretical derivation of the selection gradient in an SEIR model, we infer the phenotype of the Alpha variant in terms of transmission and recovery from the analysis of the change in its logit-frequency.

Scripts:
- **STEP2_functions.R** (script with the custom functions used in the second step)
- **STEP2_Script_England_data.Rmd** (analyses - including the mixed-effects models - using the data from England)

---

In addition:

- **Selection_coefficient_vs_Stringency_Index.R** explores the correlation between the selection coefficient of the Alpha variant and the Stringency Index in the UK; the script generates Fig. S2
- **Fig_2step_analysis_with_real_data.R** generates Fig. 1
- **Fig_Spread_variant_B117_UK.R** generates Fig. S1
- **selection_vs_Reproduction_number.R** plots the relationship between the growth rate of the epidemic and the effective reproduction number using the framework of  [Blanquart *et al.*, 2022](https://doi.org/10.7554/eLife.75791); the script generates Fig. S13

---

## Outputs

The 'Outputs' folder hosts a part of what was generated by the scripts (csv and rds files). Most of them were generated after too long running times not to be stored - e.g. non-linear optimizations - and/or are reused elsewhere - e.g. a file generated in the first step and imported in the second step.

Non-linear optimizations

- Initial values: 
**Initial_values_v4_gamma01_kappa02_pS09.csv**

- Parameter estimates: 
**Estim_tab_v4_gamma01_kappa02_pS09_R025.csv**

- Best estimates from the previous file: 
**Best_estimates_phase1_v4_gamma01_kappa02_pS09_R025.csv**

Wild bootstraps

- Using Mammen's 2-points distribution:

     (estimates) **bootstrap_wild_v4_Mammen2_gamma01_kappa02_pS09_R025.csv**

     (simulations) **simul_obs_bootstrap_wild_v4_Mammen2_gamma01_kappa02_pS09_R025.csv**

- Using Rademacher distribution:

     (estimates) **bootstrap_wild_v4_Rademacher_gamma01_kappa02_pS09_R025.csv**

     (simulations) **simul_obs_bootstrap_wild_v4_Rademacher_gamma01_kappa02_pS09_R025.csv**

Identifiability profiles

- Parameter E(t0)/N: 
**Identifiability_profile_pE_pS09_gamma01_kappa02_R025_norm_noise_sd5000_500_5.csv**
- Parameter omega: 
**Identifiability_profile_omega_pS09_gamma01_kappa02_R025_norm_noise_sd5000_500_5.csv**
- Parameter p: 
**Identifiability_profile_p_pS09_gamma01_kappa02_R025_norm_noise_sd5000_500_5.csv**
- Parameter alpha: 
**Identifiability_profile_alpha_pS09_gamma01_kappa02_R025_norm_noise_sd5000_500_5.csv**
- Parameter k: 
**Identifiability_profile_k_pS09_gamma01_kappa02_R025_norm_noise_sd5000_500_5.csv**
- Parameter a: 
**Identifiability_profile_a_pS09_gamma01_kappa02_R025_norm_noise_sd5000_500_5.csv**
- Parameter eta: 
**Identifiability_profile_eta_pS09_gamma01_kappa02_R025_norm_noise_sd5000_500_5.csv**
- Parameter mu: 
**Identifiability_profile_mu_pS09_gamma01_kappa02_R025_norm_noise_sd5000_500_5.csv**

Effects of small variations in the fixed parameters (phase 1)

- Varying gamma: 
**Optim_v4_vary_gamma_nstarts500.csv**
- Varying kappa: 
**Optim_v4_vary_kappa_nstarts500.csv**
- Varying R0 (basic reproduction number): 
**Optim_v4_vary_R0_nstarts500.csv**
- Varying S(t0)/N: 
**Optim_v4_vary_SN0_nstarts500.csv**
- Overview: 
**Estimates_with_pertubed_parameters_phase1.csv**
- Figure (R file): 
**Fig_optim_var_step1.rds**
