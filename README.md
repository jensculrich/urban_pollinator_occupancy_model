# urban_pollinator_occupancy_model

## data and code to accompany the manuscript: "Habitat restoration promotes long-term diversity in an urban pollinator metacommunity"

#### prepared by Jens Ulrich, May 2024

## the paper contains two key sets of analyses:
1) multi-species dyanmic occupancy model for estimating changes in pollinator occurrence through time.
2) logistic regression model for estimating effects of habitat enhancements on pollen limitation. 

## data are stored in ./data/:
- flower_resources_herb_quadrat.csv (the data that was collected on flower counts in the herb enhancements or analogous lawn space)
- flower_resources_woody.csv (the data that was collected on flower counts from woody plants)
- land_cover_by_site_reduced.csv (lat/long and landscape buffer impervious surface and tree canopy cover summaries for each of the 18 park sites)
- pollinator_data.csv (pollinator detections and the plants that they were recorded interacting with)
- metada.doc (contains detailed descriptions of column names from each data file)
- Pollen limitation data are stored in a separate folder. Navigate to this directory to access: ./pollen_limitation_experiment/data/clarkia_pollination_data_2022.csv.
- Appendix S2 (species detections, specialization, and phenology metrics) is available at: https://github.com/jensculrich/urban_pollinator_occupancy_model/blob/main/dynamic_occupancy_model/model_outputs/appendix_s2.xlsx. The metadata are included on the second page of the excel document.

## To conduct the occupancy model analysis:
- navigate to ./dynamic_occupancy_model/run_model/run_model.R
- the run_model.R file will access a prep_data function from prep_data.R to organize the data into the appropriate types including array of detection/non-detection for [species,site,year,survey] and the covariates (including herbaceous enhancement, woody plants, specialization).
- specifiy the appropriate stan model to estimate posterior distributions for unknown parameters and then call stan using rstan to run the model. The model used for the submitted version of the manuscript is labelled: ./dynamic_occupancy_model/models/final_model.stan
- model diagnostics and posterior distribution summaries can be pulled out at the end of the run_model.R file.
- navigate to ./dynamic_occupancy_model/PPCs/PPCs.R to run posterior predictive checks.
- navigate to ./dynamic_occupancy_model/make_figures/ to recreate figures 3 and 4 from the paper.

## To conduct the pollen limitation analysis:
- navigate to ./pollen_limitation_experiment/analysis/run_model.R
- the run_model.R file will access a prep_data function from prep_data.R to organize the data into the appropriate types including array of pollen limitation outcomes and the covariates (including herb enhancement).
- the prep_data function will access the pollen limitation data which is stored at ./pollen_limitation_experiment/data/clarkia_pollination_data_2022.csv.
- specifiy the appropriate stan model to estimate posterior distributions for unknown parameters and then call stan using rstan to run the model. The model used for the submitted version of the manuscript is labelled: ./pollen_limitation_experiment/models/logistic_model_binary_covs.stan
- model diagnostics and posterior distribution summaries can be pulled out at the end of the run_model.R file.
- posterior predictive checks can be pulled out at the end of the run_model.R file..
- navigate to ./pollen_limitation_experiment/figures/make_figures.R to recreate pollen limitation figure from the paper (Figure S23).

## Other stuff:
#### ./static_occupancy_model/
I originally built a static occupancy model (that does not tease apart changes in occurrence through time). I expanded on this by transforming the simpler structure into the dynamic structure used for the analysis presented in the manuscript.

#### ./other_analyses_and_comparisons/
- includes information on how I classified landscape context (QGIS Landscape Classification Protocol.doc)
- includes scripts to compare pollinator abundance and diversity using conventional GLMs
- includes scripts to summarize flower resources in the different parks.

#### ./figures/
- storage folder for figures used in the manuscript and supplements.
