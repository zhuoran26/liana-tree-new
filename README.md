# liana-tree-comp
Model code and code for figures and sensitivity analysis for Willson et al. 2021

Please direct questions to the corresponding author, David Medvigy, at dmedvigy@nd.edu

Short descriptions of the contents of the file, along with all dependencies, are described below.

Main_figures.R - code to create Figures 1-4 in the main text of the publication

x full_met_analysis_data_12_July.csv - compiled data from meta-analysis of hydraulic traits and those accepted from TRY. Not available but can be reproduced by following the steps in Supplementary Methods: TRY analysis & Meta-analysis

x param.input.RData - parameters & model code. Created in NPP_models.R and Input_parameter_est.R

x bci_met_mxh.RData - soil water potential and vapor pressure deficit model forcings from Barro Colorado Island, Panama. Created in BCI_met_process.R

x horizontes_met_mxh.RData - soil water potential and vapor pressure deficit model forcings from Horizontes, Costa Rica. Created in Horizontes_met_process.R

x Interpolation_mxh_100.RData - soil water potential and vapor pressure deficit interpolated indices. Dimensions as follows: soil water potential = 12 months x 100 indices, vapor pressure deficit = 12 months x 24 hours x 100 indices. Created in Sensitivity_met_interpolation.R

x tree_req_future_established_200.RData - Required tree stem-specific whole-plant hydraulic conductivity under scenarios of future increasing vapor pressure deficit (10% to 100% increase from present annual vapor pressure deficit at Barro Colorado Island and Horizontes). Created in Supplementary_figures.R

x liana_req_future_established_200.RData - Required liana stem-specific whole-plant hyraulic conductivity under scenarios of future increasing vapor pressure deficit (10% to 100% increase from present annual vapor pressure deficit at Barro Colorado Island and Horizontes). Created in Supplementary_figures.R

Supplementary_figures.R - code to create Supplementary Figures 1-23 in the supplement of the publication

x Supplement_parameters.csv - spreadsheet created in Excel to make an organized table of parameters. Formatted in R to create Supplementary Figure 1.

x full_met_analysis_data_8_Feb.csv - compiled data from meta-analysis of hydraulic traits and those accepted from TRY. Not available but can be reproduced by following the steps in Supplementary Methods: TRY analysis & Meta-analysis

x param.input.RData - parameters & model code. Created in NPP_models.R and Input_parameter_est.R

x bci_met_mxh.RData - soil water potential and vapor pressure deficit model forcings from Barro Colorado Island, Panama. Created in BCI_met_process.R

x horizontes_met_mxh.RData - soil water potential and vapor pressure deficit model forcings from Horizontes, Costa Rica. Created in Horizontes_met_process.R

x output_9oct.csv - data from TRY database. Includes mean trait values for all traits to the species level and assignment of each species to a growth form: "tree," "liana," or "other" based on the most frequently used growth form attribute for the species. Created in TRY_analysis_1.R and TRY_analysis_2.R

x Interpolation_mxh_100.RData - soil water potential and vapor pressure deficit interpolated indices. Dimensions as follows: soil water potential = 12 months x 100 indices, vapor pressure deficit = 12 months x 24 hours x 100 indices. Created in Sensitivity_met_interpolation.R
Sensitivity_analysis.R - code to create Supplementary Figures 24-27 in the supplement of the publication. This is only the extended sensitivity analysis of each parameter. In some sections, the competition models are re-written so as to allow a parameter to vary from the default value without changing the model used for the rest of the analysis.

x param.input.RData - parameters & model code. Created in NPP_models.R and Input_parameter_est.R

x bci_met_mxh.RData - soil water potential and vapor pressure deficit model forcings from Barro Colorado Island, Panama. Created in BCI_met_process.R

x horizontes_met_mxh.RData - soil water potential and vapor pressure deficit model forcings from Horizontes, Costa Rica. Created in Horizontes_met_process.R
NPP_models.R - competition model code. Tree and liana models are separate and coupled by the inputs. These models should be saved with the products of Input_parameter_est.R in one RData file for use in Main_figures.R, Supplementary figures.R, and Sensitivity_analysis.R

Input_paramter_est.R - processes data from multiple sources for use as inputs to the competition model. The products of this script should be saved with the models from NPP_models.R in one RData file for use in Main_figures.R, Supplementary figures.R, and Sensitivity_analysis.R

x full_met_analysis_data_12_July.csv - compiled data from meta-analysis of hydraulic traits and those accepted from TRY. Not avaialble but can be reproduced by following the steps in Supplementary Methods: TRY analysis & Meta-analysis

x SSD_TRY.csv - subset of data from TRY (from output_9oct.csv) that includes only species averages of stem specific density by growth form

x Horizontes_liana_tree_biomass_Martin.csv - data from Chris Smith-Martin et al. 2020 (New Phytologist)

x LianaSurveyMay2018DataForAlyssa.csv - unpublished liana survey data from Horizontes, Costa Rica from Chris Smith-Martin
TRY_analysis_1.R and TRY_analysis_2.R - processing script for all traits downloaded from the TRY plant trait database. Some traits included in this script were excluded from the final analysis due to there being too few liana species or because the definition of the trait was called into question. The data are not available, but the analysis can be repeated by following the steps outlined in Supplementary Methods: TRY analysis. All data were read in as txt files with the original name assigned from the TRY database. Species-averaged trait values and corresponding growth forms were compiled in the file output_9oct.csv

BCI_met_process.R - creates vector of average monthly soil water potential and matrix of average monthly x hourly vapor pressure deficit for Barro Colorado Island, Panama. Data are saved in an RData file to use as climate forcings for the competition model in Main_figures.R, Supplementary figures.R, and Sensitivity_analysis.R

x bci_elect_48m_at/bci_lutz48m_at_elect.csv - air temperature data downloaded from the STRI website at the Lutz Tower site for 48 m canopy height. Used in calculating vapor pressure deficit  following the procedure outlined in Supplementary Methods: Climate data

x bci_elect_48m_rh/bci_lutz48m_rh_elect.csv - relative humidity data downloaded from the STRI website at the Lutz Tower site for 48 m canopy height. Used in calculating vapor pressure deficit following the procedure outlined in Supplementary Methods: Climate data

x SWP_BCI.csv - soil water potential data for multiple soil depths at Barro Colorado Island from Levy-Varon et al. 2019 (Nature Communications)
Horizontes_met_process.R - creates vector of average monthly soil water potential and matrix of average monthly x hourly vapor pressure deficit for Horizontes, Costa Rica. Data are saved in an RData file to use as climate forcings for the competition model in Main_figures.R, Supplementary figures.R, and Sensitivity_analysis.R

x vpd_output.h5 - vapor pressure deficit data from the ERA5 re-analysis data product

x SWP_CostaRica.csv - soil water potential data for multiple soil depths at Horizontes from Medvigy et al. 2019 (New Phytologist)
Sensitivity_met_interpolation.R - creates a matrix of soil water potential values ranging from the wettest site (Barro Colorado Island) to the driest site (Horizontes) with dimensions 12 months x 100 indices and an array of vapor pressure deficit values ranging from the wettest site (Barro Colorado Island) to the driest site (Horizonte) with dimensios 12 months x 24 hours x 100 indices. Data are saved in an RData file to use as climate forcings for figures related to climate sensitivity in Main_figures.R and Supplementary figures.R

x bci_met_mxh.Rdata - soil water potential and vapor pressure deficit model forcings from Barro Colorado Island, Panama. Created in BCI_met_process.R

x horizontes_met_mxh.RData - soil water potential and vapor pressure deficit model forcings from Horizontes, Costa Rica. Created in Horizontes_met_process.R
