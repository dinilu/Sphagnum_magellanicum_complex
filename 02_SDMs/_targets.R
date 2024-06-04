# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  packages = c("tibble", 
               "xlsx", 
               "tidyverse",
               "magrittr",
               "sf",
               "terra", 
               "flexsdm",
               "targets",
               "tidyterra",
               "qs"
  ), # packages that your targets need to run
  format = "rds", # default storage format
  # Set other options as needed.
  memory = "transient",
  garbage_collection = TRUE
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# source("other_functions.R") # Source other scripts as needed. # nolint

# Replace the target list below with your own:
list(
  # // DEFINE SPECIES NAMES
  tar_target(species,
             c("diabolicum", "divinum", "magniae", "medium")),
  # // LOAD OCCURRENCE DATA
  tar_target(file,
             "../Datos/S_mag_RADseqANDbarcoding_Cleaned_for_Diego2.xlsx",
             format = "file"),
  tar_target(occ_data,
             load_occurrence_data(file, "Sheet1", species),
             pattern = map(species)),
  # // PREPARE ENVIRONMENTAL DATA AND LOAD THEM
  tar_target(vars,
             c("bio_4", "bio_5", "bio_12", "bio_15")),
  tar_target(env_file,
             # Data should be downloaded from this website: http://www.paleoclim.org/
             prep_env_data(env_folder = "../Data/PaleoClim/CHELSA_cur_V1_2B_r30s/30sec/",
                           vars = vars,
                           mask_file = NULL),
             format = "qs"),
  tar_target(env_data,
             terra::rast(env_file) %>% terra::wrap(),
             format = "qs"),
  # // DEFINE CALIBRATION AREA
  tar_target(cal_area, calc_calib_area(occ_data, env_data),
             pattern = map(occ_data),
             format = "qs"), # 500 km buffer around occurrence points
  # // FILTER DATA USING ENVIRONMENTAL FILTERING
  tar_target(occ_filt_data,
             filter_occurrence_data(occ_data, env_data),
             pattern = map(occ_data)),
  # // PARTITION DATA FOR TRAINING-TESTING
  tar_target(occ_filt_part_data,
             calc_part_occurrences(occ_filt_data, env_data),
             pattern = map(occ_filt_data),
             format = "qs"),
  tar_target(occ_filt_part,
             occ_filt_part_data$part,
             pattern = map(occ_filt_part_data)),
  # // CREATE PARTITION BLOCKS LAYER FOR BACKGROUND AND PSEUDO-ABSENCES
  tar_target(blocks_layer,
             get_blocks_layer(occ_filt_part_data, env_data),
             pattern = map(occ_filt_part_data),
             format = "qs"),
  # // CREATE BACKGROUND AND PSEUDOABSENCE DATA
  tar_target(bg,
             get_bg_data(occ_filt_part, blocks_layer, cal_area, 1:4, env_data),
             pattern = map(occ_filt_part, blocks_layer, cal_area)),
  tar_target(psa,
             get_psa_data(occ_filt_part, blocks_layer, cal_area, 1:4, env_data),
             pattern = map(occ_filt_part, blocks_layer, cal_area)),
  tar_target(species_pa,
             combine_pres_abse_data(occ_filt_part, psa, env_data),
             pattern = map(occ_filt_part, psa)),
  # // CALIBRATE MODELS
  tar_target(max_mods,
             adjust_max(species_pa, env_data, bg, vars),
             pattern = map(species_pa, bg)),
  tar_target(glm_mods,
             adjust_glm(species_pa, vars),
             pattern = map(species_pa)),
  tar_target(gau_mods,
             adjust_gau(species_pa, vars),
             pattern = map(species_pa)),
  tar_target(gam_mods,
             adjust_gam(species_pa, vars),
             pattern = map(species_pa)),
  tar_target(gbm_mods,
             adjust_gbm(species_pa, vars),
             pattern = map(species_pa)),
  tar_target(net_mods,
             adjust_net(species_pa, vars),
             pattern = map(species_pa)),
  tar_target(svm_mods,
             adjust_svm(species_pa, vars),
             pattern = map(species_pa)),
  # // ENSEMBLE MODELS
  tar_target(ens_mods,
             fit_ensemble(models = list(max_mods, 
                                        gbm_mods, 
                                        net_mods, 
                                        # svm_mods,
                                        glm_mods,
                                        # gam_mods,
                                        gau_mods),
                          ens_method = "meanw",
                          thr = c("max_sens_spec"),
                          thr_model = "max_sens_spec",
                          metric = "TSS"),
             pattern = map(max_mods, gbm_mods, net_mods, glm_mods, gau_mods)),
  # // EVALUATE MODELS (AND ENSEMBLE) PERFORMANCE
  tar_target(model_performance,
             sdm_summarize(list(max_mods, gbm_mods, net_mods, glm_mods, gau_mods, ens_mods)),
             pattern = map(max_mods, gbm_mods, net_mods, glm_mods, gau_mods, ens_mods)),
  tar_target(write_performance,
             write.xlsx(model_performance,
                        "output/model_performance.xlsx")),
  # // MAKE PREDICTIONS
  # UNDER CURRENT CONDITIONS
  tar_target(present_prediction,
             make_prediction(ens_mods, env_data, species),
             pattern = map(ens_mods, species),
             format = "qs"),
  tar_target(current_maps,
             plot_prediction(present_prediction, 
                             species,
                             "output/current_maps.pdf")),
  # UNDER PAST (LGM) CONDITIONS
  tar_target(past_env_file,
             # Data should be downloaded from this website: http://www.paleoclim.org/
             prep_past_env_data(env_folder = "../Data/PaleoClim/chelsa_LGM_v1_2B_r30s/30sec/",
                                vars = vars,
                                mask_file = NULL),
             format = "qs"),
  tar_target(lgm_env_data,
             terra::rast(past_env_file) %>% terra::wrap(),
             format = "qs"),
  tar_target(lgm_prediction,
             make_prediction(ens_mods, lgm_env_data),
             pattern = map(ens_mods),
             format = "qs"),
  tar_target(lgm_maps,
             plot_prediction(lgm_prediction, 
                             species,
                             "output/lgm_maps.pdf")),
  
  # // FUTURE PREDICTIONS

  # FUT V1.2

  # Data should be downloaded from this website: https://chelsa-climate.org/
  
  tar_target(fut_rcp26_env_data,
             prep_fut1_env_data(env_folder = "/media/NAS/Public/Data/Chelsa/v1.2/2061-80/bio/",
                                rcp = "rcp26",
                                vars = vars),
             format = "qs"),
  tar_target(fut_rcp45_env_data,
             prep_fut1_env_data(env_folder = "/media/NAS/Public/Data/Chelsa/v1.2/2061-80/bio/",
                                rcp = "rcp45",
                                vars = vars),
             format = "qs"),
  tar_target(fut_rcp60_env_data,
             prep_fut1_env_data(env_folder = "/media/NAS/Public/Data/Chelsa/v1.2/2061-80/bio/",
                                rcp = "rcp60",
                                vars = vars),
             format = "qs"),
  tar_target(fut_rcp85_env_data,
             prep_fut1_env_data(env_folder = "/media/NAS/Public/Data/Chelsa/v1.2/2061-80/bio/",
                                rcp = "rcp85",
                                vars = vars),
             format = "qs"),
  tar_target(fut_rcp26_prediction,
             make_fut_prediction(ens_mods,
                                 fut_rcp26_env_data %>%
                                   magrittr::extract2(1) %>%
                                   terra::unwrap()),
             pattern = map(ens_mods),
             format = "qs"),
  tar_target(fut_rcp45_prediction,
             make_fut_prediction(ens_mods,
                                 fut_rcp45_env_data %>%
                                   magrittr::extract2(1) %>%
                                   terra::unwrap()),
             pattern = map(ens_mods),
             format = "qs"),
  tar_target(fut_rcp60_prediction,
             make_fut_prediction(ens_mods,
                                 fut_rcp60_env_data %>%
                                   magrittr::extract2(1) %>%
                                   terra::unwrap()),
             pattern = map(ens_mods),
             format = "qs"),
  tar_target(fut_rcp85_prediction,
             make_fut_prediction(ens_mods,
                                 fut_rcp85_env_data %>%
                                   magrittr::extract2(1) %>%
                                   terra::unwrap()),
             pattern = map(ens_mods),
             format = "qs"),
  tar_target(fut_maps, 
             plot_fut_prediction(fut_rcp26_prediction, 
                                 fut_rcp45_prediction, 
                                 fut_rcp60_prediction,
                                 fut_rcp85_prediction,
                                 species,
                                 "output/fut_maps.pdf",
                                 tr = 0.4))
  
)
  

