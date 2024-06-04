library(xlsx)
library(tidyverse)
library(sf)
library(terra)
library(flexsdm)
library(targets)
# library(ape)
# library(ggfortify)

data <- read.xlsx("../Datos/S_mag_RADseqANDbarcoding_Cleaned_for_Diego2.xlsx", "Sheet1") %>%
  select(1:7) %>%
  rename(Loc = Loc.abbr) %>%
  rename(Sample = Unique_Coll_ID) %>%
  rename(Cover = Vegetation.cover) %>%
  rename(x = Lon) %>% 
  rename(y = Lat) %>% 
  mutate(Hydrology = str_to_title(Hydrology)) %>%
  mutate(Cover = str_to_title(Cover)) %>%
  mutate(Genus = "Sphagnum") %>%
  mutate(across(where(is_character), as_factor)) %>% 
  mutate(pr_ab = 1)

species_data <- data %>% select(x, y, pr_ab, Species) %>% 
  group_by(.by = Species) %>% 
  nest()

chelsa <- "/media/NAS/Public/Data/Chelsa/v2.1/1981-2010/bio/" %>%
  list.files(full.names = TRUE) %>%
  terra::rast()

names(chelsa) <- chelsa %>%
  names() %>%
  str_remove("CHELSA_") %>%
  str_remove("_1981-2010") %>%
  str_remove("_V.2.1")

somevar <- chelsa[[c("bio4", "bio12", "gdd5")]]

mask <- rast("/media/NAS/Public/Data/Chelsa/v1.2/bio/CHELSA_bio10_01.tif") 
  
somevar <- mask(somevar, mask)

loc_sf <- data %>%
  # select(genus, Lon, Lat) %>%
  st_as_sf(coords = c("x", "y"), crs = 4326)

loc_buffer <- loc_sf %>%
  st_union() %>%
  st_convex_hull() %>%
  st_buffer(dist = 100000)

somevar <- somevar %>%
  crop(loc_buffer)

rm(chelsa)

# https://sjevelazco.github.io/flexsdm/articles/v04_Red_fir_example.html

ca <- species_data[[2]] %>%
  lapply(X = ., 
         FUN = calib_area, 
         x = 'x',
         y = 'y',
         method =  c('buffer', width = 500000),
         crs = crs(somevar)
  ) # create a calibration area with 500 km buffer around occurrence points

# visualize the species occurrences
layer1 <- somevar[[1]]
layer1[!is.na(layer1)] <- 1

plot(layer1, col="gray80", legend=FALSE, axes=FALSE)
plot(crop(ca[[3]], layer1), add=TRUE)
points(species_data[[2]][[3]][,c("x", "y")], col = "#00000480")

species_data <- lapply(species_data[[2]], FUN = function(x){
  x$id <- 1:nrow(x)
  return(x)
})

species_data_filtered <- lapply(species_data, FUN = function(x){
  data <- occfilt_env(data = x,
              x = "x",
              y = "y",
              id = "id",
              nbins = 100,
              env_layer = somevar) 
  data <-  left_join(data, x, by = c("id", "x", "y"))
  return(data)
})


plot(layer1, col="gray80", legend=FALSE, axes=FALSE)
plot(crop(ca[[4]], layer1), add=TRUE)
points(species_data[[4]][,c("x", "y")], col = "#00000480")
points(species_data_filtered[[4]][,c("x", "y")], col = "red")


set.seed(10)
occ_part <- lapply(species_data_filtered,
                   FUN = function(x){
                     part_sblock(
                       data = x,
                       env_layer = somevar,
                       pr_ab = "pr_ab",
                       x = "x",
                       y = "y",
                       n_part = 4,
                       min_occ = 5
                     )
                   }
)

species_data_filtered <- lapply(occ_part, FUN = function(x) x$part)

# Transform best block partition to a raster layer with same resolution and extent than
# predictor variables
block_layer <- lapply(occ_part, FUN = function(x) get_block(env_layer = somevar, best_grid = x$grid))

cl <- c("#64146D", "#9E2962", "#F47C15", "#FCFFA4")
plot(block_layer[[4]], col=cl, legend=FALSE, axes=FALSE)
points(species_data_filtered[[4]][,c("x", "y")])




# Number of presences per block
lapply(species_data_filtered, FUN = function(x) {
  x %>% 
  dplyr::group_by(.part) %>%
  dplyr::count()
})
# Additional information of the best block
lapply(occ_part, FUN = function(x) x$best_part_info)


# Spatial blocks where species occurs
# Sample background points throughout study area with random method, allocating 10X the number of presences a background
set.seed(10)
bg <- mapply(FUN = function(sp, blayer, c, nblocks){
  lapply(nblocks, function(x, y, z, w) {
    sample_background(
      data = y,
      x = "x",
      y = "y",
      n = sum(y$.part == x) * 10,
      method = "random",
      rlayer = z,
      maskval = x,
      calibarea = w
    )
  }, sp, blayer, c) %>%
    bind_rows()
  }, species_data_filtered, block_layer, ca, MoreArgs = list(nblocks = 1:4), SIMPLIFY = FALSE)

bg <- mapply(FUN = function(x, y){sdm_extract(data = x, x = "x", y = "y", env_layer = y)}, bg, block_layer, SIMPLIFY = FALSE)

# Sample a number of pseudo-absences equal to the presence in each partition
set.seed(10)
psa <- mapply(FUN = function(sp, blayer, c, nblocks){
  lapply(nblocks, function(x, y, z, w) {
    sample_pseudoabs(
      data = y,
      x = "x",
      y = "y",
      n = sum(y$.part == x),
      method = "random",
      rlayer = z,
      maskval = x,
      calibarea = w
    )
  }, sp, blayer, c) %>%
    bind_rows()
  }, species_data_filtered, block_layer, ca, MoreArgs = list(nblocks = 1:4), SIMPLIFY = FALSE)

psa <- mapply(function(x, y) sdm_extract(data = x, x = "x", y = "y", env_layer = y), psa, block_layer, SIMPLIFY = FALSE)

cl <- c("#280B50", "#9E2962", "#F47C15", "#FCFFA4")
plot(block_layer[[4]], col="gray80", legend=FALSE, axes=FALSE)
points(bg[[4]][,c("x", "y")], col=cl[bg[[4]]$.part], cex=0.8) # Background points
points(psa[[4]][,c("x", "y")], bg=cl[psa[[4]]$.part], cex=0.8, pch=21) # Pseudo-absences




# Bind a presences and pseudo-absences
species_pa <- mapply(bind_rows, species_data_filtered, psa, SIMPLIFY = FALSE)
species_pa # Presence-Pseudo-absence database
bg # Background points

lapply(species_pa, summary)
lapply(species_pa, dim)


species_pa <- lapply(species_pa, function(x) sdm_extract(
    data = x,
    x = "x",
    y = "y",
    env_layer = somevar,
    filter_na = TRUE
  ))
bg <- lapply(bg, function(x) sdm_extract(
    data = x,
    x = "x",
    y = "y",
    env_layer = somevar,
    filter_na = TRUE
  ))




t_max <- mapply(function(x, y) tune_max(
  data = x,
  response = "pr_ab",
  predictors = names(somevar),
  background = y,
  partition = ".part",
  grid = expand.grid(
    regmult = seq(0.1, 0.5, 1, 3),
    classes = c("l", "lq", "lqhpt")
  ),
  thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
  metric = c("TSS"),
  clamp = TRUE,
  pred_type = "cloglog"
), species_pa, bg, SIMPLIFY = FALSE)


f_gau <- lapply(species_pa, function(x) fit_gau(
  data = x,
  response = "pr_ab",
  predictors = names(somevar),
  partition = ".part",
  thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
))

f_glm <- lapply(species_pa, function(x) fit_glm(
  data = x,
  response = "pr_ab",
  predictors = names(somevar),
  partition = ".part",
  thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
  poly = 2
))


ens_m <- mapply(function(t_max, f_gau, f_glm) fit_ensemble(
  models = list(t_max, f_gau, f_glm),
  ens_method = "meanw",
  thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
  thr_model = "max_sens_spec",
  metric = "TSS"
), t_max, f_gau, f_glm, SIMPLIFY = FALSE)

lapply(ens_m, function(x) x$performance)


model_perf <- mapply(function(t_max, f_gau, f_glm, ens_m) sdm_summarize(list(t_max, f_gau, f_glm, ens_m)), t_max, f_gau, f_glm, ens_m, SIMPLIFY = FALSE)


pr_1 <- lapply(ens_m, function(x) sdm_predict(
  models = x,
  pred = somevar,
  thr = "max_sens_spec",
  con_thr = TRUE,
  predict_area = NULL
))

unconstrained <- lapply(pr_1, function(x) x$meanw[[1]])
lapply(unconstrained, function(x) {names(x) <- "unconstrained"; return(x)})

cl <- c("#FDE725", "#B3DC2B", "#6DCC57", "#36B677", "#1F9D87", "#25818E", "#30678D", "#3D4988", "#462777", "#440154")
plot(unconstrained[[4]], col=cl, legend=TRUE, axes=FALSE)
lapply(model_perf)


thr_val <- lapply(ens_m, function(x){
  x$performance %>%
  dplyr::filter(threshold == "max_sens_spec") %>%
  pull(thr_value)
  }
)

m_pres <- mapply(function(x, y, z){
  msdm_posteriori(
  records = x,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  cont_suit = y$meanw[[1]],
  method = c("obr"),
  thr = c("sensitivity", sens = z),
  buffer = NULL
  )
}, species_data, pr_1, thr_val, SIMPLIFY = FALSE)

constrained <- lapply(m_pres, function(x) x$meanw[[1]])
constrained <- lapply(constrained, function(x) {
  names(x) <- "constrained"
  x})

cl <- c("#FDE725", "#B3DC2B", "#6DCC57", "#36B677", "#1F9D87", "#25818E", "#30678D", "#3D4988", "#462777", "#440154")
plot(constrained[[4]], col=cl, legend=FALSE, axes=FALSE)
