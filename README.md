R scripts used for data analysis in the manuscript: 'Extreme compound climate hazards surge by 2040: governance is critical for the world’s most vulnerable'

by Craparo, A.C.W., Minoarivelo, O. H., Basel, A.M., Nguyen, K.T., Dao., H., Montes, C., Melly, B., Birner, J., Yonetani, M., Antinoja, E.S., Sanchez Torres, D.G., Brown, O.D., Pacillo., G., Harper, A., Läderach., P.

The codes were run through RStudio 2023.06.1 using R version 4.3.1 on a device with a Windows 11 64-bit OS, an Intel Core i7-10510U @1.80GHz 2.30 GHz processor and a 16Gb RAM. The processing time for all codes to run is approzimately 50 minutes once all R required packages and their dependencies are installed correctly.

While raw input data are too big to be stored as part of the repository, climate data that are results of averaging over model ensembles are available as input within 'Ouput/Intermediate' folder.  

In order to replicate the analytics of the manuscript, the codes should be executed in the following order:

(1) The script named 'build_individual_hazard_index.R' creates the individual hazards (heat, drought, flood) indices.

(2) The script named 'build_composite_hazard_index.R' creates the composite index.

(3) The script named 'FDP_exposure_compoud.R' performs the statistical analyses with regard to the number of forcibly displaced people (FDP) within each class of compound climate hazard, both with and without governance.

(4) The script named 'FDP_exposure_compound_WB_governance.R' performs the statistical analysis related to the number of FDP exposed to climate hazards when the governance index from world bank is used instead of the BTI governance index.

(5) The scrip named 'FDP_exposre_individual_hazards.R' performs the statistical analyses with regard to the number of FDP within each class of individual climate hazards (heat, drought and flood). 

(6) The scrip named 'SSP126_analysis.R' combines the steps (1) to (4) but using projected climated hazards under SSP1-2.6.

The script 'build_mean_ensemble_from_raw_data.R' computes the (weighted) mean across model ensembles and should be the first code to be run if raw data are accessible. 
