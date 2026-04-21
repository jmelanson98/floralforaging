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

#bombus_path = "/Users/jenna1/Documents/UBC/bombus_project/"
#bombus_path = "/Users/jenna1/Documents/bombus_project/"
bombus_path = "~/projects/def-ckremen/melanson"
setwd("~/projects/def-ckremen/melanson/floralforaging")
source("src/ff_helper.R")

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
                        verbose = TRUE)
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
                        verbose = TRUE)
  saveRDS(fitmix_nocomporpoll, "analysis/stanfits/floralfit_mixtus_nocomporpoll.rds")
  
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
                        verbose = TRUE)
  saveRDS(fitimp_nocomporpoll, "analysis/stanfits/floralfit_impatiens_nocomporpoll.rds")
  
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




# Make some posterior predictions
fitmix = readRDS("analysis/floralfit_mixtus.rds")
fitimp = readRDS("analysis/floralfit_impatiens.rds")

landscapemets = read.csv("analysis/landscapemetrics.csv")
floralmean = mean(CKT_mix$floral_abundance)
floralvals = seq(min(CKT_mix$floral_abundance), max(CKT_mix$floral_abundance), by = 0.01)
ijivals = seq(min(CKT_mix$iji_centered), max(CKT_mix$iji_centered), by = 0.02)
compvals = seq(min(CKT_mix$floweringpercent_centered), max(CKT_mix$floweringpercent_centered), by = 0.1)

mixdraws = as.data.frame(fitmix)
mixrhos = mixdraws[c("rhomax", "rho0", "rho1", "rho2", "rho3")]

impdraws = as.data.frame(fitimp)
imprhos = impdraws[c("rhomax", "rho0", "rho1", "rho2", "rho3")]

mixrhos$species = "B. mixtus"
imprhos$species = "B. impatiens"

rhos = rbind(mixrhos, imprhos)
rhos_long = rhos %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(cols =c(rho1, rho2, rho3),
               names_to= "variable",
               values_to = "value")

# Set color scheme for plots
faded_pale = "#D2E4D4"
faded_light = "#B6D3B8"
faded_medium = "#609F65"
faded_strong = "#4E8353"
faded_green = "#355938"
faded_dark = "#1B2D1C"


light_gold = "#F7EAC0"
lm_gold = "#F2DC97"
medium_gold = "#ECCF6F"
gold = "#E2B41D"
dark_gold = "#B99318"
darker_gold = "#907313"

# histogram of parameters
ggplot(rhos_long, aes(x = value, colour = species, fill = species)) +
  scale_colour_manual(values = c(dark_gold, faded_green)) +
  scale_fill_manual(values = c(lm_gold, faded_light)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  facet_grid(species~variable) +
  theme_bw()

# how rho changes with covariates (floral abundance)
mix_df = cbind(data.frame(floral_abundance =floralvals,
           as.data.frame(matrix(NA, nrow = length(floralvals), ncol = 100))))
imp_df = cbind(data.frame(floral_abundance =floralvals,
                          as.data.frame(matrix(NA, nrow = length(floralvals), ncol = 100))))
           
           
for (draw in 1:100){
  row = sample(1:nrow(mixrhos),1)
  mix_df[,draw+1] = mixrhos$rhomax[row] * inv.logit(mixrhos$rho0[row] + mixrhos$rho1[row]*floralvals)
  imp_df[,draw+1] = imprhos$rhomax[row] * inv.logit(imprhos$rho0[row] + imprhos$rho1[row]*floralvals)
}


# get summary stats
y_mix = as.matrix(mix_df[ , -1])
mix_df$mean  <- rowMeans(y_mix, na.rm = TRUE)
mix_df$lower <- apply(y_mix, 1, quantile, probs = 0.025, na.rm = TRUE)
mix_df$upper <- apply(y_mix, 1, quantile, probs = 0.975, na.rm = TRUE)

y_imp = as.matrix(imp_df[ , -1])
imp_df$mean  <- rowMeans(y_imp, na.rm = TRUE)
imp_df$lower <- apply(y_imp, 1, quantile, probs = 0.025, na.rm = TRUE)
imp_df$upper <- apply(y_imp, 1, quantile, probs = 0.975, na.rm = TRUE)

ggplot() +
  geom_line(data = mix_df, aes(x = floral_abundance, y = mean), color = faded_green, size = 1) +
  geom_ribbon(data = mix_df, aes(x = floral_abundance, ymin = lower, ymax = upper), fill = faded_medium, alpha = 0.4) +
  geom_line(data = imp_df, aes(x = floral_abundance, y = mean), color = dark_gold, size = 1) +
  geom_ribbon(data = imp_df, aes(x = floral_abundance, ymin = lower, ymax = upper), fill = medium_gold, alpha = 0.4) +
  ylab(expression(rho)) +
  xlab("Floral abundance (log-scaled)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

# how rho changes with covariates (IJI)
mix_df = cbind(data.frame(iji =ijivals,
                          as.data.frame(matrix(NA, nrow = length(ijivals), ncol = 100))))
imp_df = cbind(data.frame(iji =ijivals,
                          as.data.frame(matrix(NA, nrow = length(ijivals), ncol = 100))))


for (draw in 1:100){
  row = sample(1:nrow(mixrhos),1)
  mix_df[,draw+1] = mixrhos$rhomax[row] * inv.logit(mixrhos$rho0[row] + mixrhos$rho1[row]*floralmean + mixrhos$rho3[row]*ijivals)
  imp_df[,draw+1] = imprhos$rhomax[row] * inv.logit(imprhos$rho0[row] + imprhos$rho1[row]*floralmean + imprhos$rho3[row]*ijivals)
}


# get summary stats
y_mix = as.matrix(mix_df[ , -1])
mix_df$mean  <- rowMeans(y_mix, na.rm = TRUE)
mix_df$lower <- apply(y_mix, 1, quantile, probs = 0.025, na.rm = TRUE)
mix_df$upper <- apply(y_mix, 1, quantile, probs = 0.975, na.rm = TRUE)

y_imp = as.matrix(imp_df[ , -1])
imp_df$mean  <- rowMeans(y_imp, na.rm = TRUE)
imp_df$lower <- apply(y_imp, 1, quantile, probs = 0.025, na.rm = TRUE)
imp_df$upper <- apply(y_imp, 1, quantile, probs = 0.975, na.rm = TRUE)

ggplot() +
  geom_line(data = mix_df, aes(x = iji, y = mean), color = faded_green, size = 1) +
  geom_ribbon(data = mix_df, aes(x = iji, ymin = lower, ymax = upper), fill = faded_medium, alpha = 0.4) +
  geom_line(data = imp_df, aes(x = iji, y = mean), color = dark_gold, size = 1, linetype = "dashed") +
  geom_ribbon(data = imp_df, aes(x = iji, ymin = lower, ymax = upper), fill = "lightgrey", alpha = 0.4) +
  ylab(expression(rho)) +
  xlab("Interspersion and juxtaposition index") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))


# how distance decay changes with covariates (floral abundance)
dist = seq(0,2, by = 0.01)
mix_df1 = cbind(data.frame(distance = dist,
                          as.data.frame(matrix(NA, nrow = length(dist), ncol = 100))))
mix_df2 = cbind(data.frame(distance = dist,
                           as.data.frame(matrix(NA, nrow = length(dist), ncol = 100))))
mix_df3 = cbind(data.frame(distance = dist,
                           as.data.frame(matrix(NA, nrow = length(dist), ncol = 100))))
imp_df1 = cbind(data.frame(distance = dist,
                          as.data.frame(matrix(NA, nrow = length(dist), ncol = 100))))
imp_df2 = cbind(data.frame(distance = dist,
                           as.data.frame(matrix(NA, nrow = length(dist), ncol = 100))))
imp_df3 = cbind(data.frame(distance = dist,
                           as.data.frame(matrix(NA, nrow = length(dist), ncol = 100))))


for (draw in 1:100){
  row = sample(1:nrow(mixrhos),1)
  mix_df1[,draw+1] = exp(-dist/(mixrhos$rhomax[row] * inv.logit(mixrhos$rho0[row] + mixrhos$rho1[row]*0)))
  mix_df2[,draw+1] = exp(-dist/(mixrhos$rhomax[row] * inv.logit(mixrhos$rho0[row] + mixrhos$rho1[row]*2)))
  mix_df3[,draw+1] = exp(-dist/(mixrhos$rhomax[row] * inv.logit(mixrhos$rho0[row] + mixrhos$rho1[row]*4)))
  imp_df1[,draw+1] = exp(-dist/(imprhos$rhomax[row] * inv.logit(imprhos$rho0[row] + imprhos$rho1[row]*0)))
  imp_df2[,draw+1] = exp(-dist/(imprhos$rhomax[row] * inv.logit(imprhos$rho0[row] + imprhos$rho1[row]*2)))
  imp_df3[,draw+1] = exp(-dist/(imprhos$rhomax[row] * inv.logit(imprhos$rho0[row] + imprhos$rho1[row]*4)))
}


# get summary stats
y_mix1 = as.matrix(mix_df1[ , -1])
mix_df1$mean  <- rowMeans(y_mix1, na.rm = TRUE)
mix_df1$lower <- apply(y_mix1, 1, quantile, probs = 0.025, na.rm = TRUE)
mix_df1$upper <- apply(y_mix1, 1, quantile, probs = 0.975, na.rm = TRUE)

y_mix2 = as.matrix(mix_df2[ , -1])
mix_df2$mean  <- rowMeans(y_mix2, na.rm = TRUE)
mix_df2$lower <- apply(y_mix2, 1, quantile, probs = 0.025, na.rm = TRUE)
mix_df2$upper <- apply(y_mix2, 1, quantile, probs = 0.975, na.rm = TRUE)

y_mix3 = as.matrix(mix_df3[ , -1])
mix_df3$mean  <- rowMeans(y_mix3, na.rm = TRUE)
mix_df3$lower <- apply(y_mix3, 1, quantile, probs = 0.025, na.rm = TRUE)
mix_df3$upper <- apply(y_mix3, 1, quantile, probs = 0.975, na.rm = TRUE)

ggplot() +
  geom_line(data = mix_df1, aes(x = distance, y = mean), color = faded_strong, size = 1) +
  geom_ribbon(data = mix_df1, aes(x = distance, ymin = lower, ymax = upper), fill = faded_pale, alpha = 0.4) +
  geom_line(data = mix_df2, aes(x = distance, y = mean), color = faded_green, size = 1) +
  geom_ribbon(data = mix_df2, aes(x = distance, ymin = lower, ymax = upper), fill = faded_medium, alpha = 0.4) +
  geom_line(data = mix_df3, aes(x = distance, y = mean), color = faded_dark, size = 1) +
  geom_ribbon(data = mix_df3, aes(x = distance, ymin = lower, ymax = upper), fill = faded_green, alpha = 0.4) +
  ylab(expression(rho)) +
  xlab("Distance (km)") +
  xlim(c(0,1.5)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))


# get summary stats
y_imp1 = as.matrix(imp_df1[ , -1])
imp_df1$mean  <- rowMeans(y_imp1, na.rm = TRUE)
imp_df1$lower <- apply(y_imp1, 1, quantile, probs = 0.025, na.rm = TRUE)
imp_df1$upper <- apply(y_imp1, 1, quantile, probs = 0.975, na.rm = TRUE)

y_imp2 = as.matrix(imp_df2[ , -1])
imp_df2$mean  <- rowMeans(y_imp2, na.rm = TRUE)
imp_df2$lower <- apply(y_imp2, 1, quantile, probs = 0.025, na.rm = TRUE)
imp_df2$upper <- apply(y_imp2, 1, quantile, probs = 0.975, na.rm = TRUE)

y_imp3 = as.matrix(imp_df3[ , -1])
imp_df3$mean  <- rowMeans(y_imp3, na.rm = TRUE)
imp_df3$lower <- apply(y_imp3, 1, quantile, probs = 0.025, na.rm = TRUE)
imp_df3$upper <- apply(y_imp3, 1, quantile, probs = 0.975, na.rm = TRUE)

ggplot() +
  geom_line(data = imp_df1, aes(x = distance, y = mean), color = gold, size = 1) +
  geom_ribbon(data = imp_df1, aes(x = distance, ymin = lower, ymax = upper), fill = light_gold, alpha = 0.4) +
  geom_line(data = imp_df2, aes(x = distance, y = mean), color = dark_gold, size = 1) +
  geom_ribbon(data = imp_df2, aes(x = distance, ymin = lower, ymax = upper), fill = medium_gold, alpha = 0.4) +
  geom_line(data = imp_df3, aes(x = distance, y = mean), color = darker_gold, size = 1) +
  geom_ribbon(data = imp_df3, aes(x = distance, ymin = lower, ymax = upper), fill = dark_gold, alpha = 0.4) +
  ylab(expression(rho)) +
  xlab("Distance (km)") +
  xlim(c(0,1.5)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))




# how distance decay changes with covariates (IJI)
dist = seq(0,2, by = 0.01)
mix_df1 = cbind(data.frame(distance = dist,
                           as.data.frame(matrix(NA, nrow = length(dist), ncol = 100))))
mix_df2 = cbind(data.frame(distance = dist,
                           as.data.frame(matrix(NA, nrow = length(dist), ncol = 100))))
mix_df3 = cbind(data.frame(distance = dist,
                           as.data.frame(matrix(NA, nrow = length(dist), ncol = 100))))
imp_df1 = cbind(data.frame(distance = dist,
                           as.data.frame(matrix(NA, nrow = length(dist), ncol = 100))))
imp_df2 = cbind(data.frame(distance = dist,
                           as.data.frame(matrix(NA, nrow = length(dist), ncol = 100))))
imp_df3 = cbind(data.frame(distance = dist,
                           as.data.frame(matrix(NA, nrow = length(dist), ncol = 100))))


for (draw in 1:100){
  row = sample(1:nrow(mixrhos),1)
  mix_df1[,draw+1] = exp(-dist/(mixrhos$rhomax[row] * inv.logit(mixrhos$rho0[row] + mixrhos$rho1[row]*floralmean + mixrhos$rho3[row]*-5)))
  mix_df2[,draw+1] = exp(-dist/(mixrhos$rhomax[row] * inv.logit(mixrhos$rho0[row] + mixrhos$rho1[row]*floralmean+ mixrhos$rho3[row]*0)))
  mix_df3[,draw+1] = exp(-dist/(mixrhos$rhomax[row] * inv.logit(mixrhos$rho0[row] + mixrhos$rho1[row]*floralmean + mixrhos$rho3[row]*5)))
  imp_df1[,draw+1] = exp(-dist/(imprhos$rhomax[row] * inv.logit(imprhos$rho0[row] + imprhos$rho1[row]*floralmean + imprhos$rho3[row]*-5)))
  imp_df2[,draw+1] = exp(-dist/(imprhos$rhomax[row] * inv.logit(imprhos$rho0[row] + imprhos$rho1[row]*floralmean + imprhos$rho3[row]*0)))
  imp_df3[,draw+1] = exp(-dist/(imprhos$rhomax[row] * inv.logit(imprhos$rho0[row] + imprhos$rho1[row]*floralmean + imprhos$rho3[row]*5)))
}


# get summary stats
y_mix1 = as.matrix(mix_df1[ , -1])
mix_df1$mean  <- rowMeans(y_mix1, na.rm = TRUE)
mix_df1$lower <- apply(y_mix1, 1, quantile, probs = 0.025, na.rm = TRUE)
mix_df1$upper <- apply(y_mix1, 1, quantile, probs = 0.975, na.rm = TRUE)

y_mix2 = as.matrix(mix_df2[ , -1])
mix_df2$mean  <- rowMeans(y_mix2, na.rm = TRUE)
mix_df2$lower <- apply(y_mix2, 1, quantile, probs = 0.025, na.rm = TRUE)
mix_df2$upper <- apply(y_mix2, 1, quantile, probs = 0.975, na.rm = TRUE)

y_mix3 = as.matrix(mix_df3[ , -1])
mix_df3$mean  <- rowMeans(y_mix3, na.rm = TRUE)
mix_df3$lower <- apply(y_mix3, 1, quantile, probs = 0.025, na.rm = TRUE)
mix_df3$upper <- apply(y_mix3, 1, quantile, probs = 0.975, na.rm = TRUE)

ggplot() +
  geom_line(data = mix_df1, aes(x = distance, y = mean), color = faded_strong, size = 1) +
  geom_ribbon(data = mix_df1, aes(x = distance, ymin = lower, ymax = upper), fill = faded_pale, alpha = 0.4) +
  geom_line(data = mix_df2, aes(x = distance, y = mean), color = faded_green, size = 1) +
  geom_ribbon(data = mix_df2, aes(x = distance, ymin = lower, ymax = upper), fill = faded_medium, alpha = 0.4) +
  geom_line(data = mix_df3, aes(x = distance, y = mean), color = faded_dark, size = 1) +
  geom_ribbon(data = mix_df3, aes(x = distance, ymin = lower, ymax = upper), fill = faded_green, alpha = 0.4) +
  ylab(expression(rho)) +
  xlab("Distance (km)") +
  xlim(c(0,1.5)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))


# get summary stats
y_imp1 = as.matrix(imp_df1[ , -1])
imp_df1$mean  <- rowMeans(y_imp1, na.rm = TRUE)
imp_df1$lower <- apply(y_imp1, 1, quantile, probs = 0.025, na.rm = TRUE)
imp_df1$upper <- apply(y_imp1, 1, quantile, probs = 0.975, na.rm = TRUE)

y_imp2 = as.matrix(imp_df2[ , -1])
imp_df2$mean  <- rowMeans(y_imp2, na.rm = TRUE)
imp_df2$lower <- apply(y_imp2, 1, quantile, probs = 0.025, na.rm = TRUE)
imp_df2$upper <- apply(y_imp2, 1, quantile, probs = 0.975, na.rm = TRUE)

y_imp3 = as.matrix(imp_df3[ , -1])
imp_df3$mean  <- rowMeans(y_imp3, na.rm = TRUE)
imp_df3$lower <- apply(y_imp3, 1, quantile, probs = 0.025, na.rm = TRUE)
imp_df3$upper <- apply(y_imp3, 1, quantile, probs = 0.975, na.rm = TRUE)

ggplot() +
  geom_line(data = imp_df1, aes(x = distance, y = mean), color = gold, size = 1) +
  geom_ribbon(data = imp_df1, aes(x = distance, ymin = lower, ymax = upper), fill = light_gold, alpha = 0.4) +
  geom_line(data = imp_df2, aes(x = distance, y = mean), color = dark_gold, size = 1) +
  geom_ribbon(data = imp_df2, aes(x = distance, ymin = lower, ymax = upper), fill = medium_gold, alpha = 0.4) +
  geom_line(data = imp_df3, aes(x = distance, y = mean), color = darker_gold, size = 1) +
  geom_ribbon(data = imp_df3, aes(x = distance, ymin = lower, ymax = upper), fill = dark_gold, alpha = 0.4) +
  ylab(expression(rho)) +
  xlab("Distance (km)") +
  xlim(c(0,1.5)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
