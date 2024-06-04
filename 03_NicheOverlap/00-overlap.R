# install.packages("ENMTools")
library(ENMTools)
library(tidyverse)
library(terra)
library(ggplot2)
library(ggtext)
library(xlsx)
source("R/functions.R")

env_var <- rast("../02_SDMs/output/env_data.tif")

env_var <- check.env(env_var)

myspecies <- c("diabolicum", "divinum", "magniae", "medium")

file <- "../Datos/S_mag_RADseqANDbarcoding_Cleaned_for_Diego2.xlsx"

format_species_data <- function(species, species_name, env) {
  # species <- diabolicum
  # env <- env_var
  points <- vect(species, geom = c("x", "y"))
  crs(points) <- crs(env)
  
  range <- background.buffer(points, 
                             500000, 
                             return.type = "raster", 
                             mask = env)
  
  bg.points <- background.buffer(points,
                                 500000, 
                                 return.type = "points",
                                 n = 1000, 
                                 mask = env)
  
  species.enmtools <- enmtools.species(species.name = species_name, 
                                       presence.points = points,
                                       range = range,
                                       background.points = bg.points)
  
  species.enmtools <- check.species(species.enmtools)
}

diabolicum <- load_occurrence_data(file, "Sheet1", "diabolicum")
divinum <- load_occurrence_data(file, "Sheet1", "divinum")
magniae <- load_occurrence_data(file, "Sheet1", "magniae")
medium <- load_occurrence_data(file, "Sheet1", "medium")

species_enmtools <- list(diabolicum, divinum, magniae, medium)

species_enmtools <- mapply(FUN = format_species_data, species_enmtools, myspecies, MoreArgs = list(env_var), SIMPLIFY = FALSE)

idtest <- list()
for(i in 1:(length(species_enmtools)-1)){
  tmp <- list()
  for(j in 2:length(species_enmtools)){
    if(i < j){
      cat(i, ":", j)
      tmp[[j - 1]] <- enmtools.ecospat.id(species_enmtools[[i]], species_enmtools[[j]], env = env_var, nreps = 999, bg.source = "range")
    }
  }
  names(tmp) <- myspecies[-1]
  idtest[[i]] <- tmp
}
names(idtest) <- myspecies[1:3]

bgtest <- list()
for(i in 1:(length(species_enmtools)-1)){
  tmp <- list()
  for(j in 2:length(species_enmtools)){
    if(i < j){
      cat(i, ":", j)
      tmp[[j - 1]] <- enmtools.ecospat.bg(species_enmtools[[i]], species_enmtools[[j]], env = env_var, nreps = 999, bg.source = "range")
    }
  }
  names(tmp) <- myspecies[-1]
  bgtest[[i]] <- tmp
}
names(bgtest) <- myspecies[1:3]


get_test_data <- function(obj) {
  if(!is.null(obj)){
    # obj$sp1.bg.plot$data |> rename(sp1.bg = Density) |> 
    #   left_join(obj$sp2.bg.plot$data |> rename(sp2.bg = Density)) |> 
    #   left_join(as.data.frame(obj$sp1.niche$z.cor, xy = TRUE) |> rename(sp1.niche = lyr.1, X = x, Y = y)) |> 
    #   left_join(as.data.frame(obj$sp2.niche$z.cor, xy = TRUE) |> rename(sp2.niche = lyr.1, X = x, Y = y))
    obj1 <- obj$sp1.bg.plot$data |> mutate(species = "sp1") |> 
      rbind(obj$sp2.bg.plot$data |> mutate(species = "sp2"))
    obj2 <- as.data.frame(obj$sp1.niche$z.cor, xy = TRUE) |> mutate(species = "sp1") |> rename(Density.niche = lyr.1, X = x, Y = y) |> 
      rbind(as.data.frame(obj$sp2.niche$z.cor, xy = TRUE) |> mutate(species = "sp2") |> rename(Density.niche = lyr.1, X = x, Y = y))
    obj <- obj1 |> left_join(obj2)
  }
}

get_overlap_data <- function(obj) {
  if(!is.null(obj)){
    d <- obj$test.results$obs$D
    i <- obj$test.results$obs$I
    d.p.value <- obj$p.values[["D"]]
    i.p.value <- obj$p.values[["I"]]
    df <- data.frame(d = d, i = i, d.p.value = d.p.value, i.p.value = i.p.value)
    df
  }
}

plot_test_data <- lapply(bgtest, function(x){lapply(x, get_test_data)}) |> 
  lapply(function(x){discard(x, is.null)})

plot_id_overlap_data <- lapply(idtest, function(x){lapply(x, get_overlap_data)}) |> 
  lapply(function(x){discard(x, is.null)})
plot_bg_overlap_data <- lapply(bgtest, function(x){lapply(x, get_overlap_data)}) |> 
  lapply(function(x){discard(x, is.null)})

plot_test_data <- plot_test_data |> reshape2::melt(id.vars = c("X", "Y", "species", "Density", "Density.niche"))
colnames(plot_test_data) <- c("X", "Y", "species", "density", "density.niche", "sp2", "sp1")

plot_id_overlap_data <- plot_id_overlap_data |> reshape2::melt(id.vars = c("i", "d", "i.p.value", "d.p.value"))
colnames(plot_id_overlap_data) <- c("i", "d", "i.id.p.value", "d.id.p.value", "sp2", "sp1")

plot_bg_overlap_data <- plot_bg_overlap_data |> reshape2::melt(id.vars = c("i", "d", "i.p.value", "d.p.value"))
colnames(plot_bg_overlap_data) <- c("i", "d", "i.bg.p.value", "d.bg.p.value", "sp2", "sp1")

plot_overlap_data <- plot_id_overlap_data |> left_join(plot_bg_overlap_data)

ggplot(plot_test_data) +
  geom_contour(aes(x = X, y = Y, z = density, colour = species),
               breaks = .001,
               linetype = 2,
               size = 0.5) +
  geom_contour_filled(aes(x = X, y = Y, z = density.niche, fill = species), 
                      breaks = c(.001, 1),
                      alpha = 0.5,) +
  geom_richtext(data = plot_overlap_data, 
                aes(label = paste0("D = ", format(round(d, digits = 2), nsmall = 2), 
                                   "; P<sub>id</sub> = ", format(round(d.id.p.value, 3), nsmall = 3), 
                                   "; P<sub>bg</sub> = ", format(round(d.bg.p.value, 3), nsmall = 3))), 
                x = -1,
                y = 3, 
                size = 3,
                hjust = 0,
                fill = NA, 
                label.color = NA) +
  geom_richtext(data = plot_overlap_data, 
                aes(label = paste0("I = ", format(round(i, 2), nsmall = 2), 
                                   "; P<sub>id</sub> = ", format(round(i.id.p.value, 3), nsmall = 3), 
                                   "; P<sub>bg</sub> = ", format(round(i.bg.p.value, 3), nsmall = 3))), 
                x = -1,
                y = 3.75, 
                size = 3,
                hjust = 0,
                fill = NA, 
                label.color = NA) +
  facet_grid(sp1 ~ sp2) +
  theme_bw() +
  theme(legend.position = "right",
        strip.text.x = element_text(colour = '#DC3220'),
        strip.text.y = element_text(colour = '#005AB5')) +
  scale_fill_manual(values=c("#005AB5", "#DC3220")) +
  scale_colour_manual(values=c("#005AB5", "#DC3220")) +
  coord_cartesian(ylim = c(-2.5, 4.5)) +
  labs(y = "PC2 (environmental space)", x = "PC1 (environmental space)")

ggsave("output/overlap.pdf", width = 9, height = 6.5)
