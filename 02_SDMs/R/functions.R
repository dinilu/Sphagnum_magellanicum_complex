load_occurrence_data <- function(file, sheet, species) {
  data <- read.xlsx(file, sheet) %>%
    select(1:7) %>%
    filter(Species == species) %>% 
    rename(Loc = Loc.abbr) %>%
    rename(Sample = Unique_Coll_ID) %>%
    rename(Cover = Vegetation.cover) %>%
    rename(x = Lon) %>% 
    rename(y = Lat) %>% 
    mutate(Hydrology = str_to_title(Hydrology)) %>%
    mutate(Cover = str_to_title(Cover)) %>%
    mutate(Genus = "Sphagnum") %>%
    mutate(across(where(is_character), as_factor)) %>% 
    mutate(pr_ab = 1) %>% 
    mutate(id = row_number()) %>% 
    select(id, x, y, pr_ab, Species) 
    # group_by(.by = Species) %>% 
    # nest()
  return(data)
}
 

prep_env_data <- function(env_folder, vars, mask_file) {
  loc_sf <- ext(-169, -50, 10, 75)

  chelsa <- env_folder %>%
    list.files(full.names = TRUE) %>%
    terra::rast()
 
  # names(chelsa) <- chelsa %>%
  #   names() %>%
  #   str_remove("CHELSA_") %>%
  #   str_remove("_1981-2010") %>%
  #   str_remove("_V.2.1")
 
  chelsa <- chelsa[[vars]] %>% 
    # mask(mask) %>% 
    crop(loc_sf)

  chelsa <- chelsa %>% aggregate(fact = 2, fun = mean)

  writeRaster(chelsa, "output/env_data.tif", overwrite = TRUE)

  return("output/env_data.tif")
  }


# # https://sjevelazco.github.io/flexsdm/articles/v04_Red_fir_example.html
calc_calib_area <- function(occ, env){
  ca <- flexsdm::calib_area(occ, x = 'x', y = 'y',
                              method =  c('buffer', width = 500000),
                              crs = terra::crs(terra::rast(env))) %>% 
    terra::wrap()
  ca                            
}

plot_data <- function(env, cal, occ, occ_f, blayer, bg, psa) {
# env <- targets::tar_load(env_data)
# cal <- targets::tar_load(cal_area, 1)
# occ <- targets::tar_load(occ_data, 1)
# occ_f <- targets::tar_load(occ_filt_part, 1)
# blayer <- targets::tar_load(blocks_layer, 1)
# bg <- targets::tar_load(bg, 1)
# psa <- targets:: tar_load(psa, 1)
  
    # mask <- env %>% terra::rast() %>% terra::subset(1)
  mask <- env %>% terra::unwrap() %>% terra::subset(subset = 1)
  mask[!is.na(mask)] <- 1
  mask %>% plot(col = "gray80", legend = FALSE)
  # blayer %>% terra::unwrap() %>% plot(col=cl, legend=FALSE, axes=FALSE)
  cal %>% terra::vect() %>% plot(add = TRUE)
  occ %>% st_as_sf(coords = c("x", "y")) %>% points(cex = 1.2)
  
  cl <- c("#280B50", "#9E2962", "#F47C15", "#FCFFA4")

  bg %>% st_as_sf(coords = c("x", "y")) %>% points(col = cl[bg$.part], cex=0.8, pch = 1) # Background points
  psa %>% st_as_sf(coords = c("x", "y")) %>% points(bg = cl[psa$.part], cex=0.8, pch=21) # Pseudo-absences
  # 
  occ_f %>% st_as_sf(coords = c("x", "y")) %>% points(bg = "green", pch=21)
}


filter_occurrence_data <- function(occ, env){
  data <- occfilt_env(data = occ,
              x = "x",
              y = "y",
              id = "id",
              nbins = 100,
              env_layer = rast(env)) %>% 
    left_join(occ, by = c("id", "x", "y"))
  return(data)
}


calc_part_occurrences <- function(occ_f, env){
  set.seed(1)
  data <- part_sblock(
    data = occ_f,
    env_layer = rast(env),
    pr_ab = "pr_ab",
    x = "x",
    y = "y",
    n_part = 4,
    min_occ = 4
  )
  data$grid <- terra::wrap(data$grid)
  return(data)
}

# # Transform best block partition to a raster layer with same resolution and extent than
# # predictor variables
get_blocks_layer <- function(occ_part, env) {
  layer <- get_block(env_layer = rast(env), best_grid = terra::unwrap(occ_part$grid))
  terra::wrap(layer)
}


# # Number of presences per block
# lapply(species_data_filtered, FUN = function(x) {
#   x %>% 
#   dplyr::group_by(.part) %>%
#   dplyr::count()
# })
# # Additional information of the best block
# lapply(occ_part, FUN = function(x) x$best_part_info)

# # Spatial blocks where species occurs
# # Sample background points throughout study area with random method, allocating 10X the number of presences a background
get_bg_data <- function(sp, blayer, c, nblocks, env){
  set.seed(10)
  lapply(nblocks, 
         FUN = function(i){
           sample_background(data = sp,
                             x = "x",
                             y = "y",
                             n = sum(sp$.part == i) * 10,
                             method = "random",
                             rlayer = terra::unwrap(blayer),
                             maskval = nblocks,
                             calibarea = terra::vect(c)
           )
         }
  ) %>%
    bind_rows() %>% 
    sdm_extract(x = "x", y = "y", env_layer = terra::unwrap(blayer)) %>% 
    sdm_extract(x = "x", y = "y", env_layer = terra::rast(env), filter_na = TRUE)
}

# # Sample a number of pseudo-absences equal to the presence in each partition
get_psa_data <- function(sp, blayer, c, nblocks, env){
  set.seed(10)
  lapply(nblocks, 
         FUN = function(i) {
           sample_pseudoabs(data = sp,
                            x = "x",
                            y = "y",
                            n = sum(sp$.part == i),
                            method = "random",
                            rlayer = terra::unwrap(blayer),
                            maskval = nblocks,
                            calibarea = terra::vect(c)
           )
         }
  ) %>%
    bind_rows() %>% 
    sdm_extract(x = "x", y = "y", env_layer = terra::unwrap(blayer)) %>% 
    sdm_extract(x = "x", y = "y", env_layer = terra::rast(env), filter_na = TRUE)
} 


combine_pres_abse_data <- function(occ_filt_part, psa, env){
  occ_filt_part %>% sdm_extract(x = "x", y = "y", env_layer = terra::rast(env)) %>% 
    bind_rows(occ_filt_part, psa)
}


adjust_max <- function(species, env, bg, vars){
  tune_max(
    data = species,
    response = "pr_ab",
    predictors = vars,
    background = bg,
    partition = ".part",
    grid = expand.grid(
      regmult = seq(0.5, 1, 1.5),
      classes = c("l", "lq", "lqhpt")
    ),
    thr = c("max_sens_spec"),
    metric = c("TSS"),
    clamp = TRUE,
    pred_type = "cloglog"
  )
}


adjust_gau <- function(species, vars){
  fit_gau(data = species,
          response = "pr_ab",
          predictors = vars,
          partition = ".part",
          thr = c("max_sens_spec")
  )
} 


adjust_gam <- function(species, vars){
  fit_gam(data = species,
          response = "pr_ab",
          predictors = vars,
          partition = ".part",
          thr = c("max_sens_spec")
  )
} 


adjust_glm <- function(species, vars){
  fit_glm(
    data = species,
    response = "pr_ab",
    predictors = vars,
    partition = ".part",
    thr = c("max_sens_spec"),
    poly = 2
  )
}              


adjust_gbm <- function(species, vars){
  
  tune_grid <- expand.grid(
    n.trees = c(500, 1000, 1500),
    shrinkage = c(0.1, 0.5, 1),
    n.minobsinnode = c(1, 3, 5, 7, 9)
  )
  
  tune_gbm(
    data = species,
    response = "pr_ab",
    predictors = vars,
    partition = ".part",
    grid = tune_grid,
    thr = c("max_sens_spec"),
    metric = "TSS"
  )
}              


adjust_net <- function(species, vars){
  
  tune_grid <- expand.grid(
    size = c(2, 4, 6, 8, 10),
    decay = c(0.001, 0.05, 0.1, 1, 5, 10)
  )
  
  tune_net(
    data = species,
    response = "pr_ab",
    predictors = vars,
    partition = ".part",
    grid = tune_grid,
    thr = c("max_sens_spec"),
    metric = "TSS"
  )
}              


adjust_svm <- function(species, vars){

  tune_grid <-
    expand.grid(
      C = c(2, 4, 8, 16, 20),
      sigma = c(0.01, 0.1, 0.2, 0.3, 0.4)
    )
  
  tune_svm(
    data = species,
    response = "pr_ab",
    predictors = vars,
    partition = ".part",
    grid = tune_grid,
    thr = c("max_sens_spec"),
    metric = "TSS"
  )
}              




# lapply(ens_m, function(x) x$performance)

make_prediction <- function(model, env, species){
  pred <- sdm_predict(
    models = model,
    pred = terra::rast(env),
    thr = "max_sens_spec",
    con_thr = TRUE
  )
  # writeRaster(pred$meanw, paste0("output/", species, "_current.tif"), overwrite = TRUE)
  pred$meanw <- wrap(pred$meanw)
  return(pred)
} 


make_fut_prediction <- function(model, env){
  pred <- sdm_predict(
    models = model,
    pred = env,
    thr = "max_sens_spec",
    con_thr = TRUE
  )
  pred$meanw <- wrap(pred$meanw)
  return(pred)
} 


plot_prediction <- function(pred, species, filename, tr = 0.4){
  # pred <- present_prediction
  pred <- lapply(pred, FUN = terra::unwrap) %>% 
    lapply(FUN = terra::subset, subset = 1) %>% 
    rast() %>% 
    aggregate(fact = 10, fun = "mean")
  
  names(pred) <- paste0("S. ", species)
  
  pred <- pred %>% as.data.frame(xy = TRUE) %>% 
    pivot_longer(cols = all_of(paste0("S. ", species)),
                 names_to = "species",
                 values_to = "suitability") %>%
    rename(Longitude = x) %>% 
    rename(Latitude = y)

  # pred_wo_zero <- pred %>% filter(suitability != 0)  
  
  maps <- ggplot() +
    # geom_tile(data = pred, aes(x = x, y = y), fill = "grey") +
    geom_tile(data = pred, aes(x = Longitude, y = Latitude, col = suitability)) +
    # geom_tile(data = pred_wo_zero, aes(x = x, y = y, col = suitability)) +
    # scale_fill_viridis_c(limits = c(0.4, 1)) +
    scale_color_viridis_c(limits = c(tr, 1), na.value = "grey", name = "Suitability") +
    facet_wrap(~species, ncol = 2) +
    theme_bw() +
    coord_equal() + 
    theme(strip.text = element_text(face = "italic"))
  
  ggsave(filename, 
         plot = maps,
         width = 8,
         height = 4.5)
  
  return("Done")
}


prep_past_env_data <- function(env_folder, vars, mask_file) {
  
  loc_sf <- ext(-169, -50, 10, 75)

  # mask <- terra::rast(mask_file) > 0 
  # NAflag(mask) <- 0
  
  chelsa <- env_folder %>%
    list.files(full.names = TRUE) %>%
    terra::rast()
  
  # names(chelsa) <- chelsa %>%
  #   names() %>%
  #   str_remove("CHELSA_TraCE21k_") %>%
  #   str_remove("_-200") %>%
  #   str_remove("_V1.0")
  
  chelsa <- chelsa[[vars]] %>% 
    # mask(mask) %>% 
    crop(loc_sf)

  # names(chelsa) <- c("bio_4", "bio_5", "bio_12", "bio_15")
  
  chelsa <- chelsa %>% aggregate(fact = 2, fun = mean)

  writeRaster(chelsa, "output/lgm_env_data.tif", overwrite = TRUE)
  
  return("output/lgm_env_data.tif")
}


prep_fut_env_data <- function(env_folder, gcm, ssp, vars, mask_file) {

  loc_sf <- ext(-169, -50, 10, 75)
  
  env_folder <- paste0(env_folder, "/", gcm, "/", ssp, "/bio")
  
  chelsa <- env_folder %>%
    list.files(full.names = TRUE) %>%
    terra::rast()

  names(chelsa) <- chelsa %>%
    names() %>%
    str_remove("CHELSA_") %>%
    str_remove("_2071-2100_(.*)_V.2.1")

  chelsa <- chelsa[[paste0("bio", 1:19)]]
  names(chelsa) <- paste0("bio_", 1:19)
  
  mask <- mask_file %>% rast()
  
  chelsa <- chelsa[[vars]] %>% 
    mask(mask) %>% 
    crop(loc_sf)
  
  chelsa <- chelsa %>% aggregate(fact = 2, fun = mean)
  
  chelsa <- list(chelsa %>% terra::wrap())

  names(chelsa) <- paste0(gcm, "-", ssp)
  
  chelsa
  # writeRaster(chelsa, "output/lgm_env_data.tif", overwrite = TRUE)
  
  # return("output/lgm_env_data.tif")
}



prep_fut1_env_data <- function(env_folder, rcp, vars) {
  # env_folder <- "/media/NAS/Public/Data/Chelsa/v1.2/2061-80/bio/"
  # rcp <- "rcp26"
  # targets::tar_load(vars)
  
  loc_sf <- ext(-169, -50, 10, 75)
  
  ivars <- vars %>% str_remove("bio_") %>% as.numeric()
  
  env_files <- lapply(ivars, FUN = function(i) list.files(env_folder, pattern = paste0(rcp, "_(.*)_", i, "_"), full.names = TRUE))
  
  .load_rast <- function(x, loc) {
    x %>% terra::rast() %>% 
      terra::crop(loc) %>%
      terra::aggregate(fact = 2, fun = mean) %>% 
      terra::mean()
  }
  
  env_rast <- lapply(env_files, FUN = .load_rast, loc = loc_sf) %>% 
    Reduce(c, .)
  
  names(env_rast) <- vars
  
  env_rast <- list(env_rast %>% terra::wrap())
  
  names(env_rast) <- rcp

  env_rast  
}


plot_fut_prediction <- function(pred1, pred2, pred3, pred4, species, filename, tr = 0.4){
  
  .rast_to_long_df <- function(pred, sp, rcp){
    pred <- lapply(pred, FUN = terra::unwrap) %>% 
      lapply(FUN = terra::subset, subset = 1) %>% 
      rast() %>% 
      aggregate(fact = 10, fun = "mean")
    
    names(pred) <- paste0("S. ", sp)
    
    pred <- pred %>% as.data.frame(xy = TRUE) %>% 
      pivot_longer(cols = all_of(paste0("S. ", sp)),
                   names_to = "species",
                   values_to = "suitability") %>% 
      mutate(rcp = rcp)  %>%
      rename(Longitude = x) %>% 
      rename(Latitude = y)
    pred
  }
  
  pred1 <- .rast_to_long_df(pred1, species, "rcp26")
  pred2 <- .rast_to_long_df(pred2, species, "rcp45")
  pred3 <- .rast_to_long_df(pred3, species, "rcp60")
  pred4 <- .rast_to_long_df(pred4, species, "rcp85")
  
  pred <- rbind(pred1, pred2, pred3, pred4)
  
  maps <- ggplot() +
    # geom_tile(data = pred, aes(x = x, y = y), fill = "grey") +
    geom_tile(data = pred, aes(x = Longitude, y = Latitude, col = suitability)) +
    # geom_tile(data = pred_wo_zero, aes(x = x, y = y, col = suitability)) +
    # scale_fill_viridis_c(limits = c(0.4, 1)) +
    scale_color_viridis_c(limits = c(tr, 1), na.value = "grey", name = "Suitability") +
    facet_grid(rcp~species) +
    theme_bw() +
    coord_equal() + 
    theme(strip.text.x = element_text(face = "italic"))
  
  
  ggsave(filename, 
         plot = maps,
         width = 8,
         height = 5)
  
  return("Done")
}
