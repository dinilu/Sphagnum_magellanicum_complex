
# // LOAD LIBRARIES

library(xlsx)
library(tidyverse)
library(factoextra)
library(sf)
library(terra)
library(ape)
library(ggfortify)
library(Hmsc)
library(ggpubr)


# // LOAD OCCURRENCE DATA

data <- read.xlsx("../Data/S_mag_RADseqANDbarcoding_Cleaned_for_Diego2.xlsx", "Sheet1") %>%
  select(1:7) %>%
  rename(Loc = Loc.abbr) %>%
  rename(Sample = Unique_Coll_ID) %>%
  rename(Cover = Vegetation.cover) %>%
  mutate(Hydrology = str_to_title(Hydrology)) %>%
  mutate(Cover = str_to_title(Cover)) %>%
  mutate(genus = "Sphagnum") %>%
  mutate(across(where(is_character), as_factor)) %>%
  mutate(Cover = fct_relevel(Cover, c("Partially Shady", "Shady", "Sunny"))) %>%
  mutate(Species = fct_relevel(Species, c("diabolicum", "divinum", "magniae", "medium")))


# // LOAD PHYLOGENETIC DATA

tree <- read.tree(file = "../Data/Mag.newick")
species <- unique(sapply(strsplit(tree$tip.label,"_"),function(x) x[2]))
ii <- sapply(species, function(x,y) {grep(x,y)[1]}, y = tree$tip.label)
tree <- drop.tip(tree, setdiff(tree$tip.label, tree$tip.label[ii]))
tree$tip.label <- species
plot(tree)
rm(species, ii)


# // SELECT SITES WITH MORE THAN 5 SAMPLES

sel_loc <- data %>%
  group_by(Loc) %>%
  summarize(n = n()) %>%
  filter(n >= 5)

data <- data %>%
  filter(Loc %in% sel_loc$Loc)


# // TRANSFORM OCCURRENCE DATA FRAME INTO A COMMUNITY MATRIX

cm <- data %>%
  select(Species, Sample) %>%
  mutate(Occ = 1) %>%
  pivot_wider(names_from = Species, values_from = Occ,
              values_fill = 0, values_fn = sum) %>%
  column_to_rownames(var = "Sample") %>%
  mutate(across(colnames(.), ~ ifelse(.x > 0, 1, 0))) %>%
  select(colnames(.) %>% sort()) %>%
  as.matrix()


# // EXTRACT GEOGRAPHIC INFORMATION AT OCCURRENCE LOCATIONS

loc_sf <- data %>%
  # select(genus, Lon, Lat) %>%
  st_as_sf(coords = c("Lon", "Lat"), crs = 4326)


# // CREATE A BUFFER AREA AROUND LOCATIONS

loc_buffer <- loc_sf %>%
  st_union() %>%
  st_convex_hull() %>%
  st_buffer(dist = 100000)


# // DOWNLOAD ENVIRONMENTAL DATA FROM CHELSA WEBSITE (RUN ONCE)

# setwd("../Data")
# system2("wget", args = c("--no-host-directories", "--force-directories", "--input-file=envidatS3paths.txt"))
# chelsa.files <- list.files("envicloud/chelsa/chelsa_V1/GLOBAL/climatologies/bio/", full.names = TRUE)
# dir.create("Chelsa")
# unlink("Data/envicloud/", recursive = TRUE)
# file.rename(from = chelsa.files, to = paste0("Chelsa/", basename(chelsa.files)))
# setwd("../01_HMSC")


# // LOAD CLIMATE DATA AND TAILOR THEM

# chelsa <- "/media/NAS/Public/Data/Chelsa/v1.2/bio/" %>%
chelsa <- "../Data/Chelsa/v1.2/bio/" %>%
  list.files(full.names = TRUE) %>%
  terra::rast()

names(chelsa) <- chelsa %>%
  names() %>%
  str_remove("CHELSA_") %>%
  str_remove("10_[0]{0,1}")

chelsa <- chelsa %>%
  crop(loc_buffer)

plot(chelsa[["bio12"]])
plot(loc_buffer, add = TRUE)
plot(loc_sf, col = "black", add = TRUE)


# // EXTRACT CLIMATE VALUES AT LOCATIONS AND REMOVE WHEN THERE ARE NAs

env <- chelsa %>%
  extract(loc_sf) %>%
  dplyr::select(where(~!any(is.na(.))))


# // CALCULATE PCA ON CLIMATE DATA

env_pca <- env %>%
  select(-ID) %>%
  prcomp(scale.=TRUE)


# // PLOT PCAs

# Define transparency for all variables but a few selected (BIO4, BIO12, and CMI_MIN)
# pca_alpha <- rep(0.2, 64)
# pca_alpha[c(5, 15, 27)] <- 1
pca_alpha <- rep(0.2, 19)
pca_alpha[c(4, 5, 12, 15)] <- 1
# pca_col <- rep("grey", 64)
# pca_col[c(5, 15, 27)] <- "black"
pca_col <- rep("grey", 19)
pca_col[c(4, 5, 12, 15)] <- "black"

pca_a <- autoplot(env_pca,
                  loadings = TRUE,
                  loadings.label = TRUE,
                  loadings.colour = pca_col,
                  loadings.label.colour  = "black",
                  loadings.label.alpha = pca_alpha)

pca_b <- autoplot(env_pca,
                  loadings = TRUE,
                  loadings.label = TRUE,
                  loadings.colour = pca_col,
                  loadings.label.colour = "black",
                  y = 3,
                  loadings.label.alpha = pca_alpha)

pdf("output/Supplementary/pca.pdf", width = 7, height = 3)
  ggarrange(pca_a, pca_b, ncol = 2)
dev.off()



# // ADD FIRST 4 PCA AXIS TO ENVIRONMENTAL DATASET, AS WELL AS HABITAT VARIABLES

evplot <- function(ev)
{
  # Broken stick model (MacArthur 1957)
  n <- length(ev)
  bsm <- data.frame(j=seq(1:n), p=0)
  bsm$p[1] <- 1/n
  for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n + 1 - i))
  bsm$p <- 100*bsm$p/n
  # Plot eigenvalues and % of variation for each axis
  op <- par(mfrow=c(2,1))
  barplot(ev, main="Eigenvalues", col="bisque", las=2)
  abline(h=mean(ev), col="red")
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE,
          main="% variation", col=c("bisque",2), las=2)
  legend("topright", c("% eigenvalue", "Broken stick model"),
         pch=15, col=c("bisque",2), bty="n")
  par(op)
}

evplot(env_pca$sdev^2)

env <- env_pca$x %>%
  as.data.frame() %>%
  select(1:4) %>%
  bind_cols(env, .) %>%
  mutate(Hydrology = data$Hydrology, Cover = data$Cover) %>%
  mutate(across(where(is_character), as_factor)) %>%
  mutate(Sample = data$Sample) %>%
  distinct(Sample, .keep_all = TRUE)


# // PREPARE ALL DATA FOR HMSC WITH COMPLETE CASES

# Complete cases index
cc_index <- env %>%
  complete.cases()

# Community matrix
cm_cc <- cm %>%
  .[cc_index,]

# Environmental dataset
env_cc <- env %>%
  filter(cc_index)

# Original dataset for random effect levels
data_cc <- data %>%
  filter(Sample %in% env_cc$Sample) %>%
  mutate(Sample = as.factor(Sample)) %>%
  droplevels()

# Random effect at the sample level
rl_sample <- data_cc$Sample %>%
  HmscRandomLevel(units = .)

# Random effect at the location level
rl_loc <- data_cc$Loc %>%
  HmscRandomLevel(units = .)

# HMSC formula
# form <- ~bio4+I(bio4^2)+gdd5+I(gdd5^2)+bio12+I(bio12^2)+Hydrology+Cover # Tested this
# form <- ~PC1+I(PC1^2)+PC2+I(PC2^2)+PC3+I(PC3^2)+Hydrology+Cover # and this...
form <- ~PC1+PC2+PC3+PC4+Hydrology+Cover # Finally kept this.


# // DEFINE HMSC MODEL

m <- Hmsc(Y = cm_cc, XData = env_cc, XFormula = form, studyDesign = select(data_cc, Sample, Loc), ranLevels = list("Sample" = rl_sample, "Loc" = rl_loc), phyloTree = tree, distr = "probit")


# // RUN CALIBRATE HMSC POSTERIORS

m <- sampleMcmc(m, samples = 1000, transient = 10000, thin = 100, nChains = 2, nParallel = 2)


# // CONVERT TO CODA OBJECT

mpost <- convertToCodaObject(m)


# // EVALUATE MODEL CONVERGENCE

# Effective sample size
es <- effectiveSize(mpost$Beta)

# Gelman diagnostics
gd <- gelman.diag(mpost$Beta, multivariate = FALSE)$psrf

# Plot ES and GD histograms
pdf(file = "output/Supplementary/Beta_convergence.pdf", width = 6, height = 4)
par(mfrow=c(1,2))
hist(es, main = "", xlab = "Effective sample size")
hist(gd, main = "", xlab = "Gelman diagnostic")
dev.off()

# Plot chain convergence
pdf(file = "output/Supplementary/Chain_convergence.pdf", height = 9)
plot(mpost$Beta)
dev.off()


# // EVALUATE MODEL PERFORMANCE

# On training dataset
preds_train <- computePredictedValues(m)
evaluateModelFit(hM = m, predY = preds_train) %>%
  as.data.frame() %>%
  mutate(species = c("diabolicum", "divinum", "magniae", "medium")) %>%
  write.xlsx("output/training_evaluation.xlsx")

# On four random testing partitions
partition <- createPartition(m, nfolds = 4)
preds_test <- computePredictedValues(m, partition = partition, nParallel = 2)
evaluateModelFit(hM = m, predY = preds_test) %>%
  as.data.frame() %>%
  mutate(species = c("diabolicum", "divinum", "magniae", "medium")) %>%
  write.xlsx("output/test_evaluation.xlsx")


# // STUDY ENVIRONMENTAL VARIABLES EFFECT ON SPECIES OCCURRENCE

postBeta <- getPostEstimate(m, parName = "Beta")
pdf(file = "output/Betas.pdf", width = 6, height = 6)
plotBeta(m,
         post = postBeta,
         param = "Mean",
         supportLevel = 0.95,
         plotTree = TRUE,
         spNamesNumbers = c(1, 0),
         covNamesNumbers = c(1, 0),
         split = 0.25,
         mar = c(10,1,1,0.5),
         marTree = c(11,1.5,2,0.5))
dev.off()

round(summary(mpost$Rho, quantiles = c(0.025, 0.5, 0.975))[[1]],2)


# // ANALYZE RESIDUAL CO-OCCURRENCE MATRIX

hist(effectiveSize(mpost$Omega[[1]]), main = "ess(omega)")
hist(gelman.diag(mpost$Omega[[1]], multivariate = FALSE)$psrf, main = "psrf(omega)")

plot_omega <- function(m, n_random, supportLevel = 0.95) {
  OmegaCor <- computeAssociations(m)

  tosample <- ((OmegaCor[[n_random]]$support > supportLevel)
            + (OmegaCor[[n_random]]$support < (1 - supportLevel)) > 0) * OmegaCor[[n_random]]$mean

  corrplot::corrplot(tosample, method = "color",
           col = colorRampPalette(c("blue", "white", "red"))(200),
           title = paste("Random Effect:", m$rLNames[n_random]),
           mar = c(0,0,1,0),
           tl.col = "black")
}

plot_omega(m, n_random = 1)
pdf(file = "output/Omegas.pdf", width = 4, height = 4)
plot_omega(m, n_random = 2)
dev.off()


# // ANALYZE VARIANCE PARTITIONING

par(mfrow=c(1,1))
VP <- computeVariancePartitioning(m)
plotVariancePartitioning(m, VP = VP)
VP <- computeVariancePartitioning(m, group = c(1, 2, 3, 4, 5, 6, 6, 7, 7), groupnames = c("intercept", "PC1", "PC2", "PC3", "PC4", "Hydrology", "Cover"))
plotVariancePartitioning(m, VP = VP)

VP <- computeVariancePartitioning(m, group = c(1, 1, 1, 1, 1, 2, 2, 2, 2), groupnames = c("Climate", "Habitat"))

pdf(file = "output/VariancePartitioning.pdf", width = 7, height = 5)
plotVariancePartitioning(m, VP = VP, args.legend = list(x = 2.3, y = 0.97, bg = "white", cex = 0.8))
dev.off()


# // PLOT RESPONSE CURVES

plot_response_curve <- function(var, m) {

  source("custom_functions.R")
  Gradient <- constructGradient(m,
                                focalVariable = var,
                                non.focalVariables = 1)

  predY <- predict(m,
                   XData = Gradient$XDataNew,
                   studyDesign = Gradient$studyDesignNew,
                   ranLevels = Gradient$rLNew,
                   expected = TRUE)

  if(var == "Hydrology"){
    a <- plotGradient(m, Gradient, pred=predY, measure="Y", index=1, showData = TRUE, jigger = 0.2) + scale_x_discrete(limits = c("Hollow", "Low Hummock", "High Hummock"))
    b <- plotGradient(m, Gradient, pred=predY, measure="Y", index=2, showData = TRUE, jigger = 0.2) + scale_x_discrete(limits = c("Hollow", "Low Hummock", "High Hummock"))
    c <- plotGradient(m, Gradient, pred=predY, measure="Y", index=3, showData = TRUE, jigger = 0.2) + scale_x_discrete(limits = c("Hollow", "Low Hummock", "High Hummock"))
    d <- plotGradient(m, Gradient, pred=predY, measure="Y", index=4, showData = TRUE, jigger = 0.2) + scale_x_discrete(limits = c("Hollow", "Low Hummock", "High Hummock"))
    ggarrange(a, b, c, d, ncol = 2, nrow = 2)
  } else if(var == "Cover") {
    a <- plotGradient(m, Gradient, pred=predY, measure="Y", index=1, showData = TRUE, jigger = 0.2) + scale_x_discrete(limits = c("Shady", "Partially Shady", "Sunny"))
    b <- plotGradient(m, Gradient, pred=predY, measure="Y", index=2, showData = TRUE, jigger = 0.2) + scale_x_discrete(limits = c("Shady", "Partially Shady", "Sunny"))
    c <- plotGradient(m, Gradient, pred=predY, measure="Y", index=3, showData = TRUE, jigger = 0.2) + scale_x_discrete(limits = c("Shady", "Partially Shady", "Sunny"))
    d <- plotGradient(m, Gradient, pred=predY, measure="Y", index=4, showData = TRUE, jigger = 0.2) + scale_x_discrete(limits = c("Shady", "Partially Shady", "Sunny"))
    ggarrange(a, b, c, d, ncol = 2, nrow = 2)
  }else {
    a <- plotGradient(m, Gradient, pred = predY, measure = "Y", index = 1, showData = TRUE)
    b <- plotGradient(m, Gradient, pred = predY, measure = "Y", index = 2, showData = TRUE)
    c <- plotGradient(m, Gradient, pred = predY, measure = "Y", index = 3, showData = TRUE)
    d <- plotGradient(m, Gradient, pred = predY, measure = "Y", index = 4, showData = TRUE)
    ggarrange(a, b, c, d, ncol = 2, nrow = 2)
  }
}

pdf(file = "output/PC1_rc.pdf", width = 8, height = 5.5)
plot_response_curve("PC1", m)
dev.off()

pdf(file = "output/PC2_rc.pdf", width = 8, height = 5.5)
plot_response_curve("PC2", m)
dev.off()

pdf(file = "output/PC3_rc.pdf", width = 8, height = 5.5)
plot_response_curve("PC3", m)
dev.off()

pdf(file = "output/PC4_rc.pdf", width = 8, height = 5.5)
plot_response_curve("PC4", m)
dev.off()

pdf(file = "output/Hydrology_rc.pdf", width = 8, height = 5.5)
plot_response_curve("Hydrology", m)
dev.off()

pdf(file = "output/Cover_rc.pdf", width = 8, height = 5.5)
plot_response_curve("Cover", m)
dev.off()


