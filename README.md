# liana-tree-comp
## Model code and code for figures and sensitivity analysis for Willson et al. 2022

### Please direct questions to the corresponding author, David Medvigy, at dmedvigy@nd.edu

### Short descriptions of the contents of the file, along with all dependencies, are described below.

**Main_figures.R** - code to create Figures 1-5 in the main text of the publication

* *full_met_analysis_data_8_Feb.csv* - compiled data from meta-analysis of hydraulic traits and those accepted from TRY. Download from figshare repository associated with the manuscript

* *param.input.RData* - parameters. Created in Compile_parameters.R

* *bci_met_mxh.RData* - soil water potential and vapor pressure deficit model forcings from Barro Colorado Island, Panama. Created in BCI_met_process.R

* *horizontes_met_mxh.RData* - soil water potential and vapor pressure deficit model forcings from Horizontes, Costa Rica. Created in Horizontes_met_process.R

* *Interpolation_mxh_100.RData* - soil water potential and vapor pressure deficit interpolated indices. Dimensions as follows: soil water potential = 12 months x 100 indices, vapor pressure deficit = 12 months x 24 hours x 100 indices. Created in Sensitivity_met_interpolation.R

* *tree.NPP.mxh* & *liana.NPP.mxh* - tree and liana models. Called from NPP_models.R

* *tree.NPP.mxh.ca* & *liana.NPP.mxh.ca* - tree and liana models with additional input paramter atmospheric carbon dioxide concentration. Used to increase atmospheric carbon dioxide concentration for future simulations. Called from NPP_models_wCO2.R

**Supplementary_figures.R** - code to create Supplementary Figures & Tables of the manuscript

* *full_met_analysis_data_8_Feb.csv* - compiled data from meta-analysis of hydraulic traits and those accepted from TRY. Download from the figshare repository associated with this manuscript

* *param.input.RData* - parameters. Created in Compile_parameters.R

* *bci_met_mxh.RData* - soil water potential and vapor pressure deficit model forcings from Barro Colorado Island, Panama. Created in BCI_met_process.R

* *horizontes_met_mxh.RData* - soil water potential and vapor pressure deficit model forcings from Horizontes, Costa Rica. Created in Horizontes_met_process.R

* *filtered_TRY_analysis_16-03-21.csv* - data from TRY database. Includes mean trait values for all traits to the species level and assignment of each species to a growth form: "tree," "liana," or "other" based on the most frequently used growth form attribute for the species. Created in TRY_analysis_1.R and TRY_analysis_2.R and available in the fighsare repository associated with this manuscript

* *Interpolation_mxh_100.RData* - soil water potential and vapor pressure deficit interpolated indices. Dimensions as follows: soil water potential = 12 months x 100 indices, vapor pressure deficit = 12 months x 24 hours x 100 indices. Created in Sensitivity_met_interpolation.R 

* *tree.NPP.mxh* & *liana.NPP.mxh* - tree and liana models. Sourced from NPP_models.R

* *tree.NPP.mxh.ca* & *liana.NPP.mxh.ca* - tree and liana models with additional input parameter atmospheric carbon dioxide concentration. Used to increase atmospheric carbon dioxide concentration for future simulations. Sourced from NPP_models_wCO2.R

**Sensitivity_analysis.R** - code to create Supplementary Figures & Tables from the manuscript. This is only the extended sensitivity analysis of each parameter. In some sections, the competition models are re-written so as to allow a parameter to vary from the default value without changing the model used for the rest of the analysis.

* *param.input.RData* - parameters. Created in Compile_parameters.R

* *bci_met_mxh.RData* - soil water potential and vapor pressure deficit model forcings from Barro Colorado Island, Panama. Created in BCI_met_process.R

* *horizontes_met_mxh.RData* - soil water potential and vapor pressure deficit model forcings from Horizontes, Costa Rica. Created in Horizontes_met_process.R

* *tree.NPP.mxh* & *liana.NPP.mxh* - tree and liana models. Sourced from NPP_models.R

**NPP_models.R** - competition model code. Tree and liana models are separate and coupled by the inputs. These models are used in Main_figures.R, Supplementary_figures.R, and Sensitivity_analysis.R

**NPP_models_wCO2.R** - competition model code with variable carbon dioxide concentration input. Otherwise, identical to NPP_models.R. These models are used in Main_figures.R and Supplementary_figures.R

**Compile_parameters.R** - Contains data informed parameters for the competition model. The products of this script should be saved in one RData file for use in Main_figures.R, Supplementary figures.R, and Sensitivity_analysis.R

* *full_met_analysis_data_8_Feb.csv* - compiled data from meta-analysis of hydraulic traits and those accepted from TRY. Available from the figshare repository associated with this manuscript

**TRY_analysis_1.R** & **TRY_analysis_2.R** - processing script for all traits downloaded from the TRY plant trait database. Some traits included in this script were excluded from the final analysis due to there being too few liana species or because the definition of the trait was called into question. The data are not available, but the analysis can be repeated by following the steps outlined in Supplementary Methods: TRY analysis. All data were read in as txt files with the original name assigned from the TRY database. Species-averaged trait values and corresponding growth forms were compiled in the file *output_9oct.csv* (TRY_analysis_1.R) and *filtered_TRY_analysis_16-03-21.csv* (TRY_analysis_2.R)

**BCI_met_process.R** - creates vector of average monthly soil water potential and matrix of average monthly x hourly vapor pressure deficit for Barro Colorado Island, Panama. Data are saved in an RData file to use as climate forcings for the competition model in Main_figures.R, Supplementary figures.R, and Sensitivity_analysis.R

* *bci_elect_48m_at/bci_lutz48m_at_elect.csv* - air temperature data downloaded from the STRI website at the Lutz Tower site for 48 m canopy height. Used in calculating vapor pressure deficit following the procedure outlined in Supplementary Methods: Climate data

* *bci_elect_48m_rh/bci_lutz48m_rh_elect.csv* - relative humidity data downloaded from the STRI website at the Lutz Tower site for 48 m canopy height. Used in calculating vapor pressure deficit following the procedure outlined in Supplementary Methods: Climate data

* *SWP_BCI.csv* - soil water potential data for multiple soil depths at Barro Colorado Island from Levy-Varon et al. 2019 (Nature Communications) 

**Horizontes_met_process.R** - creates vector of average monthly soil water potential and matrix of average monthly x hourly vapor pressure deficit for Horizontes, Costa Rica. Data are saved in an RData file to use as climate forcings for the competition model in Main_figures.R, Supplementary figures.R, and Sensitivity_analysis.R

* *vpd_output.h5* - vapor pressure deficit data from the ERA5 re-analysis data product

* *SWP_CostaRica.csv* - soil water potential data for multiple soil depths at Horizontes from Medvigy et al. 2019 (New Phytologist) 

**Sensitivity_met_interpolation.R** - creates a matrix of soil water potential values ranging from the wettest site (Barro Colorado Island) to the driest site (Horizontes) with dimensions 12 months x 100 indices and an array of vapor pressure deficit values ranging from the wettest site (Barro Colorado Island) to the driest site (Horizonte) with dimensios 12 months x 24 hours x 100 indices. Data are saved in an RData file to use as climate forcings for figures related to climate sensitivity in Main_figures.R and Supplementary figures.R

* *bci_met_mxh.Rdata* - soil water potential and vapor pressure deficit model forcings from Barro Colorado Island, Panama. Created in BCI_met_process.R

* *horizontes_met_mxh.RData* - soil water potential and vapor pressure deficit model forcings from Horizontes, Costa Rica. Created in Horizontes_met_process.R
