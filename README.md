# ALPHA

**Phenotypic evolution of SARS-CoV-2: a statistical inference approach**

Wakinyan Benhamou, Sébastien Lion, Rémi Choquet and Sylvain Gandon

CEFE, CNRS, Univ Montpellier, EPHE, IRD, Montpellier, France

## Data

This is publicly available data that were downloaded from the Internet. All data used in the scripts are in the 'Data' folder. 

### covid-stringency-index.csv

&#11169;&emsp; Stringency Index (from 21/01/2020 to 23/03/2021), computed by the Oxford COVID-19 Government Response Tracker (OxCGRT) [[Hale *et al.*, 2021](https://doi.org/10.1038/s41562-021-01079-8)]. This file was originally downloaded from the website *Our World in Data* but does not seem to be available there anymore... However, values from old databases can still be found in a GitHub repositery of OxCGRT: '[covid-policy-tracker-legacy](https://github.com/OxCGRT/covid-policy-tracker-legacy)'.

### daily-tests-and-daily-new-confirmed-covid-cases.csv

&#11169;&emsp; Daily number of tests with new cases tested positive (from 07/04/2020 to 28/04/2021), downloaded from *Our World in Data*:

https://ourworldindata.org/grapher/daily-tests-and-daily-new-confirmed-covid-cases?country=~GBR

### data_daily_deaths_UK_2021_May_05.csv

&#11169;&emsp; Daily new fatality cases (from 02/02/2020 to 23/04/2021), downloaded from *GOV.UK* (visited 05/05/2021):

https://coronavirus.data.gov.uk/details/deaths

(see 'Daily deaths with COVID-19 on the death certificate by date of death')

### Variant_of_Concern_VOC_202012_01_Technical_Briefing_5_Data_England.ods

[DOWNLOAD](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/957631/Variant_of_Concern_VOC_202012_01_Technical_Briefing_5_Data_England.ods) (ods file)

&#11169;&emsp; Fequencies of SGTF (S-Gene Target Failures) in the 9 regions of England, from *Public Health England Technical Briefing 5*:

https://www.gov.uk/government/publications/investigation-of-novel-sars-cov-2-variant-variant-of-concern-20201201

## R codes

We developed a two-step approach based on a SEIR model and the analysis of a combination of epidemiological and evolutionary information:

(1) In the first step, we infer how the Stringency Index reduces the amount of viral transmission just before the emergence of the Alpha Variant;

Scripts:
- **STEP1_functions.R** (script with functions used in the first step)
- **STEP1_Script_UK_data.Rmd** (analyses using the UK data)
- **STEP1_Identifiability_profiles.Rmd** (build and plot identifiability profile of our model)

(2) In the second step, based on a novel theoretical derivation of the selection gradient in a SEIR model, we infer the phenotype of the Alpha variant from the analysis of the change in its logit-frequency.

Scripts:
- **STEP2_functions.R** (script with functions used in the second step)
- **STEP2_Script_England_data.Rmd** (analyses - including the mixed-effects model - using the data from England)

------------------------------------------------------------------------------------------------

In addition:

- **Selection_coefficient_vs_Stringency_Index.R** explores the correlation between the selection coefficient of the Alpha variant and the Stringency Index in the UK; the script generates Fig. S2;

- **Fig_2step_analysis_with_real_data.R** generates Fig. 1;

- **Fig_Spread_variant_B117_UK.R** generates Fig. S1.

------------------------------------------------------------------------------------------------

The 'Outputs' folder hosts a part of what was generated by the scripts and are reused elsewhere in codes - e.g. a file generated in the first step and imported in the second step.
