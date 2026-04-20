##### Fit floral foraging distance models
### April 15, 2026
### J Melanson


### Load in packages
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(ggplot2)
library(rstan)
library(sf)
source("src/ff_helper.R")

bombus_path = "/Users/jenna1/Documents/bombus_project/"

# Load in data
mixtus_sibs2022 = read.csv("cleandata/siblingships/mixtus_sibships_2022.csv")
mixtus_sibs2023 = read.csv("cleandata/siblingships/mixtus_sibships_2023.csv")
impatiens_sibs2022 = read.csv("cleandata/siblingships/impatiens_sibships_2022.csv")
impatiens_sibs2023 = read.csv("cleandata/siblingships/impatiens_sibships_2023.csv")
effort2022 = read.csv("cleandata/fielddata/2022sampledata.csv")
effort2023 = read.csv("cleandata/fielddata/2023sampledata.csv")
veg2022 = read.csv("cleandata/fielddata/2022vegetationdata.csv")
veg2023 = read.csv("cleandata/fielddata/2023vegetationdata.csv")
landscapemets = read.csv("analysis/landscapemetrics_nocrops.csv")


# Prep data for stan
data_mix = prep_stan_floralforaging(sibships1 = mixtus_sibs2022,
                         sibships2 = mixtus_sibs2023,
                         effort1 = effort2022,
                         effort2 = effort2023,
                         veg1 = veg2022,
                         veg2 = veg2023,
                         landscapemetrics = landscapemets)

data_imp = prep_stan_floralforaging(sibships1 = impatiens_sibs2022,
                                sibships2 = impatiens_sibs2023,
                                effort1 = effort2022,
                                effort2 = effort2023,
                                veg1 = veg2022,
                                veg2 = veg2023,
                                landscapemetrics = landscapemets)
stan_data_mix = data_mix[[1]]
stan_data_imp = data_imp[[1]]
CKT_mix = data_mix[[2]]
CKT_imp = data_imp[[2]]
               
# Fit models in stan
stanfile = "models/floral_foraging.stan"
stan_data_mix$deltaprior = 1
stan_data_imp$deltaprior = 1

fitmix_nocrops = stan(file = stanfile,
           data = stan_data_mix, seed = 5838299,
           chains = 4, cores = 4,
           iter = 4000, warmup = 1000,
           verbose = TRUE)
saveRDS(fitmix_nocrops, "analysis/floralfit_mixtus_nocrops.rds")


fitimp_nocrops = stan(file = stanfile,
                      data = stan_data_imp, seed = 5838299,
                      chains = 4, cores = 4,
                      iter = 4000, warmup = 1000,
                      verbose = TRUE)
saveRDS(fit, "analysis/floralfit_impatiens_nocrops.rds")

mixfit = readRDS("analysis/floralfit_mixtus.rds")
View(summary(mixfit)$summary)
