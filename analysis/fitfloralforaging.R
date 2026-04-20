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

#bombus_path = "/Users/jenna1/Documents/UBC/bombus_project/"
#bombus_path = "/Users/jenna1/Documents/bombus_project/"
bombus_path = "~/projects/def-ckremen/melanson"

# Get task ID from slurm manager
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

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
floralprot = read.csv("cleandata/plantdata/proteinbyfamily.csv")


# All four vars
if (task_id == 1){
  # Stan file
  stanfile = "models/floral_foraging.stan"
  
  # Get data
  data_mix = prep_stan_floralforaging(sibships1 = mixtus_sibs2022,
                                      sibships2 = mixtus_sibs2023,
                                      effort1 = effort2022,
                                      effort2 = effort2023,
                                      veg1 = veg2022,
                                      veg2 = veg2023,
                                      landscapemetrics = landscapemets,
                                      pollenprotein = floralprot)
  stan_data_mix = data_mix[[1]]
  CKT_mix = data_mix[[2]]
  stan_data_mix$deltaprior = 1
  
  
  # Fitmodel
  fitmix_all = stan(file = stanfile,
                        data = stan_data_mix, seed = 5838299,
                        chains = 4, cores = 4,
                        iter = 4000, warmup = 1000,
                        verbose = TRUE)
  saveRDS(fitmix_all, "analysis/stanfits/floralfit_mixtus_all.rds")
  
} else if(task_id == 2){
  stanfile = "models/floral_foraging.stan"
  
  # Get impatiens data
  data_imp = prep_stan_floralforaging(sibships1 = impatiens_sibs2022,
                                      sibships2 = impatiens_sibs2023,
                                      effort1 = effort2022,
                                      effort2 = effort2023,
                                      veg1 = veg2022,
                                      veg2 = veg2023,
                                      landscapemetrics = landscapemets,
                                      pollenprotein = floralprot)
  
  stan_data_imp = data_imp[[1]]
  CKT_imp = data_imp[[2]]
  stan_data_imp$deltaprior = 1
  
  # Fit model
  fitimp_all = stan(file = stanfile,
                        data = stan_data_imp, seed = 5838299,
                        chains = 4, cores = 4,
                        iter = 4000, warmup = 1000,
                        verbose = TRUE)
  saveRDS(fitimp_all, "analysis/stanfits/floralfit_impatiens_all.rds")
  
} else if (task_id == 3){
  
  ## No composition
  
  # Get stan file
  stanfile = "models/floral_foraging_nocomp.stan"
  
  #Get mixtus data
  data_mix = prep_stan_floralforaging(sibships1 = mixtus_sibs2022,
                                      sibships2 = mixtus_sibs2023,
                                      effort1 = effort2022,
                                      effort2 = effort2023,
                                      veg1 = veg2022,
                                      veg2 = veg2023,
                                      landscapemetrics = landscapemets,
                                      pollenprotein = floralprot)
  stan_data_mix = data_mix[[1]]
  CKT_mix = data_mix[[2]]
  stan_data_mix$deltaprior = 1
  stan_data_mix$comp = NULL
  
  # Fit model
  fitmix_nocomp = stan(file = stanfile,
                        data = stan_data_mix, seed = 5838299,
                        chains = 4, cores = 4,
                        iter = 4000, warmup = 1000,
                        verbose = TRUE)
  saveRDS(fitmix_nocomp, "analysis/stanfits/floralfit_mixtus_nocomp.rds")
  
} else if(task_id == 4){
  stanfile = "models/floral_foraging_nocomp.stan"
  
  
  # Get impatiens data
  data_imp = prep_stan_floralforaging(sibships1 = impatiens_sibs2022,
                                      sibships2 = impatiens_sibs2023,
                                      effort1 = effort2022,
                                      effort2 = effort2023,
                                      veg1 = veg2022,
                                      veg2 = veg2023,
                                      landscapemetrics = landscapemets,
                                      pollenprotein = floralprot)
  
  stan_data_imp = data_imp[[1]]
  CKT_imp = data_imp[[2]]
  stan_data_imp$deltaprior = 1
  stan_data_imp$comp = NULL
  
  fitimp_nocomp = stan(file = stanfile,
                        data = stan_data_imp, seed = 5838299,
                        chains = 4, cores = 4,
                        iter = 4000, warmup = 1000,
                        init = init, verbose = TRUE)
  saveRDS(fitimp_nocomp, "analysis/stanfits/floralfit_impatiens_nocomp.rds")
  
} else if (task_id == 5){
  
  ## No composition or protein
  stanfile = "models/floral_foraging_nocomp_noprot.stan"
  data_mix = prep_stan_floralforaging(sibships1 = mixtus_sibs2022,
                                      sibships2 = mixtus_sibs2023,
                                      effort1 = effort2022,
                                      effort2 = effort2023,
                                      veg1 = veg2022,
                                      veg2 = veg2023,
                                      landscapemetrics = landscapemets,
                                      pollenprotein = floralprot)
  stan_data_mix = data_mix[[1]]
  CKT_mix = data_mix[[2]]
  stan_data_mix$deltaprior = 1
  stan_data_mix$comp = NULL
  stan_data_mix$mp = NULL
  
  fitmix_nocomporpoll = stan(file = stanfile,
                        data = stan_data_mix, seed = 5838299,
                        chains = 4, cores = 4,
                        iter = 4000, warmup = 1000,
                        init = init, verbose = TRUE)
  saveRDS(fitmix_nocomporpoll, "analysis/floralfit_mixtus_nocomporpoll.rds")
  
} else if(task_id == 6){
  stanfile = "models/floral_foraging_nocomp_noprot.stan"
  # Fit impatiens model
  data_imp = prep_stan_floralforaging(sibships1 = impatiens_sibs2022,
                                      sibships2 = impatiens_sibs2023,
                                      effort1 = effort2022,
                                      effort2 = effort2023,
                                      veg1 = veg2022,
                                      veg2 = veg2023,
                                      landscapemetrics = landscapemets,
                                      pollenprotein = floralprot)
  
  stan_data_imp = data_imp[[1]]
  CKT_imp = data_imp[[2]]
  stan_data_imp$deltaprior = 1
  stan_data_imp$comp = NULL
  stan_data_imp$mp = NULL
  
  fitimp_nocomporpoll = stan(file = stanfile,
                        data = stan_data_imp, seed = 5838299,
                        chains = 4, cores = 4,
                        iter = 4000, warmup = 1000,
                        init = init, verbose = TRUE)
  saveRDS(fitimp_nocomporpoll, "analysis/floralfit_impatiens_nocomporpoll.rds")
  
} else if (task_id == 7){
  
  # Use crop inclusive floral habitat! 
  landscapemets = read.csv("analysis/landscapemetrics.csv")
  
  # Stan file
  stanfile = "models/floral_foraging.stan"
  
  # Get data
  data_mix = prep_stan_floralforaging(sibships1 = mixtus_sibs2022,
                                      sibships2 = mixtus_sibs2023,
                                      effort1 = effort2022,
                                      effort2 = effort2023,
                                      veg1 = veg2022,
                                      veg2 = veg2023,
                                      landscapemetrics = landscapemets,
                                      pollenprotein = floralprot)
  stan_data_mix = data_mix[[1]]
  CKT_mix = data_mix[[2]]
  stan_data_mix$deltaprior = 1
  
  
  # Fitmodel
  fitmix_pheno = stan(file = stanfile,
                    data = stan_data_mix, seed = 5838299,
                    chains = 4, cores = 4,
                    iter = 4000, warmup = 1000,
                    verbose = TRUE)
  saveRDS(fitmix_pheno, "analysis/stanfits/floralfit_mixtus_pheno.rds")
  
} else if(task_id == 8){
  # Use crop inclusive floral habitat! 
  landscapemets = read.csv("analysis/landscapemetrics.csv")
  
  stanfile = "models/floral_foraging.stan"
  
  # Get impatiens data
  data_imp = prep_stan_floralforaging(sibships1 = impatiens_sibs2022,
                                      sibships2 = impatiens_sibs2023,
                                      effort1 = effort2022,
                                      effort2 = effort2023,
                                      veg1 = veg2022,
                                      veg2 = veg2023,
                                      landscapemetrics = landscapemets,
                                      pollenprotein = floralprot)
  
  stan_data_imp = data_imp[[1]]
  CKT_imp = data_imp[[2]]
  stan_data_imp$deltaprior = 1
  
  # Fit model
  fitimp_pheno = stan(file = stanfile,
                    data = stan_data_imp, seed = 5838299,
                    chains = 4, cores = 4,
                    iter = 4000, warmup = 1000,
                    verbose = TRUE)
  saveRDS(fitimp_pheno, "analysis/stanfits/floralfit_impatiens_pheno.rds")
  
}





