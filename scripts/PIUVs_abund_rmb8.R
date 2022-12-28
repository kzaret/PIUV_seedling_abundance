#Hierarchical linear mixed models of Pilgerodendron uviferum seedling abundance at the control (Raul Marin Balmaceda) site; data from 2017 and 2018.

#===========================================================================
# SETUP
#===========================================================================

library(here)
library(dplyr)
library(tidyverse)
library(rstanarm)
library(rstan)
library(brms)
library(bayesplot)
library(feasts)
library(lubridate)
library(hms)
library(shinystan)
library(tidybayes)
library(bayesrules)
library(broom.mixed) # use tidy function to tailor summary output
library(loo) # cross-validation
#library(posterior) #masks ff from package:stats -- mad, sd, var
library(RColorBrewer) #plots
library(modelr) #plots

rstan_options(auto_write = TRUE)

bayesplot_theme_set(ggplot2::theme_bw())
color_scheme_set("brightblue")

#===============================================================================
# READ IN & TIDY DATA
#===============================================================================

# read in subplot dataset
sp <- read.csv(here("data/Subplot_PIUVcounts_ExplanVar_2017to2018.csv"), header = TRUE)

# tidy date & time
sp$date <- as.POSIXlt(sp$date, format = '%m/%d/%Y') #convert date field to date
sp$par_time <- as_hms(sp$par_time)


# tidy par
sp <- sp %>% mutate (log_par0.3 = log10(par0.3),
                     log_par1.3 = log10(par1.3))


# tidy dwt
sp.dwt <- sp %>% 
  pivot_longer(c("dwt1", "dwt2", "dwt3"), names_to = "dwt_msmnt", values_to = "dwt_cm")

sp.dwt.s <- sp.dwt %>%
  group_by(subplot) %>%
  summarise (
    avgdwt = mean(dwt_cm, na.rm=T),
    mindwt = min(dwt_cm, na.rm=T),
    maxdwt = max(dwt_cm, na.rm=T),
    rangedwt = (maxdwt - mindwt),
  ) %>%
  mutate(avgdwt = round(avgdwt, 0))

sp <- left_join(sp, sp.dwt.s, by = "subplot")

# convert continuous avg dwt values to ordinal categories (low, <= 0; medium, 1 to 25; high, >25)
sp <- sp%>%
  mutate(DWT = cut(avgdwt, breaks=c(-Inf, 0, 25, Inf), labels=c("low", "medium", "high")))
sp$DWT <- factor(sp$DWT, c("low", "medium", "high"), ordered = TRUE)


#filter & set factor levels by site
sp.rmb <- sp %>% filter(site=="RMB")
sp.rmb$patch <- factor(sp.rmb$patch, c("Forest", "Bog Forest", "Bog"))
sp.rmb$plot = as.factor(sp.rmb$plot)


# scale variables by site for use w/ brms models
sp.rmb <- sp.rmb %>%
  mutate (log_PAR0.3.rs = as.vector(scale(log_par0.3)),
          log_PAR1.3.rs = as.vector(scale(log_par1.3)),
          Pool.rs = as.vector(scale(pool)),
          Plot_PIUVcnt.rs = as.vector(scale(Plot_PIUVcnt)),
          Trunk.rs = as.vector(scale(trunk)),
          Sphagnum.rs = as.vector(scale(Sphagnum)),
          Cushion_spp.rs = as.vector(scale(cush_spp)),
          Dicran_spp.rs = as.vector(scale(Dicran)),
          OtherMoss_spp.rs = as.vector(scale(om)),
          Graminoids.rs = as.vector(scale(graminoids)),
          Forbs.rs = as.vector(scale(forbs)),
          Ferns.rs = as.vector(scale(ferns)),
          LowShrubs.rs = as.vector(scale(low_shrubs)),
          OtherShrubs.rs = as.vector(scale(other_shrubs)),
          Shrubs.rs = as.vector(scale(all_shrubs)),
          US_Trees.rs = as.vector(scale(us_trees)),
          US_PIUV.rs = as.vector(scale(us_PIUV)),
          OtherTrees.rs = as.vector(scale(OS_great130cm)))
          

#verify & summarize tidied data by site
glimpse(sp.rmb)
summary(sp.rmb)
View(sp.rmb)

#write.csv(sp.rmb, "C:/Users/Marsh Hawk Zaret/Documents/Big_D/Data_Analysis/RQ1_Neo/Data/Ordination/rmb_explanvars_Mar2022.csv")

#################################################################################### DETERMINING OPTIMAL HIERARCHICAL STRUCTURE & LIKELIHOOD
################################################################################

#===============================================================================
# Intercept-only models
#===============================================================================

#------------
# Fit models
#------------

# note:  'k' refers to set of priors [see sensitiviy script]

#nbnk -- nb intercept-only without random effects
rmsa.glmer.nbnk <- brm(bf(PIUV_less5cmDBH_total ~ 1), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior ("gamma(1, 0.1)", class = "shape"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nbnk, file="F:/_R_objects/rmsa.glmer.nbnk.rda")

#nb0k -- nb intercept-only with plot as random effects
rmsa.glmer.nb0k <- brm(bf(PIUV_less5cmDBH_total ~ 1 + (1 | plot)), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior ("gamma(1, 0.1)", class = "shape"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb0k, file="F:/_R_objects/rmsa.glmer.nb0k.rda")

#nb0k.2 -- nb intercept-only with plot as random effects & shape modeled by patch
rmsa.glmer.nb0k.2 <- brm(bf(PIUV_less5cmDBH_total ~ 1 + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0,5)", class = "b", dpar = "shape"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb0k.2, file="F:/_R_objects/rmsa.glmer.nb0k.2.rda")

#zinb0k -- zinb intercept-only w/ plot as random effects
rmsa.glmer.zinb0k <- brm(bf(PIUV_less5cmDBH_total ~ 1 + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior ("normal(0, 5)", class = "b", dpar="zi") + set_prior ("gamma(1, 0.1)", class = "shape"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb0k, file="F:/_R_objects/rmsa.glmer.zinb0k.rda")


load(file="F:/_R_objects/rmsa.glmer.zinb0k.rda")
load(file="F:/_R_objects/rmsa.glmer.nb0k.rda")
load(file="F:/_R_objects/rmsa.glmer.nb0k.2.rda")
load(file="F:/_R_objects/rmsa.glmer.nbnk.rda")


#--------------
# Check models
#--------------

mod <- rmsa.glmer.nb1k.2
mod_name <- "rmsa.glmer.nb1k.2"

prior_summary(mod)
summary(mod)
variables(mod)


launch_shinystan(mod)

bayesplot_theme_set(ggplot2::theme_bw())
color_scheme_set("brightblue")

# density & trace plots
dev.new()
plot(mod)

dev.new()
mcmc_trace(mod)

# posterior vs. prior
variables(mod)
param <- "b_shape_Intercept"
prior_summary(mod)    # inspect prior for param
pdf <- function(x) dnorm(x, 0, 5) # density function for prior on param
pdf <- function(x) dstudent_t(x, 3, 0, 2.5)
pdf <- function(x) dgamma(x, 1, 0.1)
pdf <- function(x) dlogis(x, 0, 1)


hist(as.matrix(mod, variable = param), 25, prob = TRUE,
     col = "gray", border = "white", 
     las = 1, cex.axis = 1.2, cex.lab = 1.5, cex.main = 1.5,
     xlim = c(-10, 10),
     #ylim = c(0, .3),
     xlab = param, ylab = "Probability", main = paste("prior vs. posterior", mod_name))
curve(pdf(x), col = "blue", lwd = 2, add = TRUE)

exp(5)
log(10)
plogis(2.5) #inverse logit


# posterior predictive checks
pp_check(mod, ndraws = 250) + 
  xlab("Seedling counts") + ggtitle(mod_name)

y <- sp.rmb$PIUV_less5cmDBH_total # vector of outcome values
yrep <- posterior_predict(mod, draws = 500) # matrix of draws from the ppd

prop_zero <- function(x) mean(x == 0)

ppc_stat(y, yrep, stat = "prop_zero", binwidth = 0.005) +
  ggtitle(mod_name)
ppc_stat_grouped(y, yrep, sp.rmb$patch, stat = "prop_zero", binwidth = 0.005) +
  ggtitle(mod_name)
ppc_stat(y, yrep, stat = "max", binwidth = 5) + ggtitle(mod_name)
ppc_stat_grouped(y, yrep, sp.rmb$patch, stat = "max", binwidth = 5) + 
  ggtitle(mod_name)
ppc_stat(y, yrep, stat = "sd", binwidth = 5) +
  ggtitle(mod_name)
ppc_stat_grouped(y, yrep, sp.rmb$patch, stat = "sd", binwidth = 5) +
  ggtitle(mod_name)


# interval estimates
get_variables(mod)
dev.new()
mcmc_intervals(mod, prob = 0.90, prob_outer = 0.95,pars=c("b_Intercept", "b_zi_Intercept", "b_log_PAR1.3.rs", "b_patchBogForest", "b_patchBog", "b_log_PAR1.3.rs:patchBogForest", "b_log_PAR1.3.rs:patchBog", "b_zi_patchBogForest", "b_zi_patchBog", "sd_plot__Intercept", "shape")) + ggtitle(mod_name)

dev.new()
mcmc_intervals(mod, prob = 0.90, prob_outer = 0.95,pars=c("b_Intercept", "b_log_PAR1.3.rs", "b_patchBogForest", "b_patchBog", "b_log_PAR1.3.rs:patchBogForest", "b_log_PAR1.3.rs:patchBog", "sd_plot__Intercept", "b_shape_Intercept", "b_shape_patchBogForest", "b_shape_patchBog", "Intercept_shape")) + ggtitle(mod_name)

# pairwise plots
dev.new()
mcmc_pairs(as.matrix(mod),pars = c("b_Intercept", "b_zi_Intercept", "b_log_PAR1.3.rs", "b_patchBogForest", "b_patchBog", "b_log_PAR1.3.rs:patchBogForest", "b_log_PAR1.3.rs:patchBog", "b_zi_patchBogForest", "b_zi_patchBog", "sd_plot__Intercept", "shape")) + ggtitle(mod_name)

posterior <- as.matrix(mod)

dev.new()
mcmc_pairs(posterior, pars = c("b_Intercept", "b_zi_Intercept", "b_log_PAR1.3.rs", "b_patchBogForest", "b_patchBog", "b_log_PAR1.3.rs:patchBogForest", "b_log_PAR1.3.rs:patchBog", "b_zi_patchBogForest", "b_zi_patchBog", "sd_plot__Intercept", "shape")) + ggtitle(mod_name)


# conditional effects
plot(conditional_effects(mod), ask = FALSE) + ggtitle(mod_name)



#--------------------------
# Validate & compare models
#--------------------------

loo_rmsa.nbnk <- loo(rmsa.glmer.nbnk, save_psis = TRUE, moment_match = TRUE)
save(loo_rmsa.nbnk, file="F:/_R_objects/loo_rmsa.nbnk.rda")

loo_rmsa.nb0k <- loo(rmsa.glmer.nb0k, save_psis = TRUE, moment_match = TRUE)
save(loo_rmsa.nb0k, file="F:/_R_objects/loo_rmsa.nb0k.rda")

loo_rmsa.nb0k.2 <- loo(rmsa.glmer.nb0k.2, save_psis = TRUE, moment_match = TRUE)
save(loo_rmsa.nb0k.2, file="F:/_R_objects/loo_rmsa.nb0k.2.rda")

loo_rmsa.zinb0k <- loo(rmsa.glmer.zinb0k, save_psis = TRUE, moment_match = TRUE)
save(loo_rmsa.zinb0k, file="F:/_R_objects/loo_rmsa.zinb0k.rda")


load(file="F:/_R_objects/loo_rmsa.nbnk.rda")
load(file="F:/_R_objects/loo_rmsa.nb0k.rda")
load(file="F:/_R_objects/loo_rmsa.nb0k.2.rda")
load(file="F:/_R_objects/loo_rmsa.zinb0k.rda")


loo_compare(loo_rmsa.nbnk, loo_rmsa.nb0k, loo_rmsa.nb0k.2, loo_rmsa.zinb0k)


#===============================================================================
# Model 1 - Patch
#===============================================================================

#------------
# Fit models
#------------

#zinb1k
rmsa.glmer.zinb1k <- brm(bf(PIUV_less5cmDBH_total ~ patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b") + set_prior ("normal(0, 5)", class = "b", dpar="zi") + set_prior ("gamma(1, 0.1)", class = "shape"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb1k, file="F:/_R_objects/rmsa.glmer.zinb1k.rda")

#nb1k
rmsa.glmer.nb1k <- brm(bf(PIUV_less5cmDBH_total ~ patch + (1 | plot)), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b") + set_prior ("gamma(1, 0.1)", class = "shape"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb1k, file="F:/_R_objects/rmsa.glmer.nb1k.rda")

#nb1k.2 [w/ shape modeled by patch]
rmsa.glmer.nb1k.2 <- brm(bf(PIUV_less5cmDBH_total ~ patch + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb1k.2, file="F:/_R_objects/rmsa.glmer.nb1k.2.rda")


load(file="F:/_R_objects/rmsa.glmer.zinb1k.rda")
load(file="F:/_R_objects/rmsa.glmer.nb1k.rda")
load(file="F:/_R_objects/rmsa.glmer.nb1k.2.rda")


mod <- rmsa.glmer.nb1k.2
mod_name <- "rmsa.glmer.nb1k.2"

prior_summary(mod)
summary(mod)
variables(mod)


#---------------------------
# Validate & compare models
#---------------------------

loo_rmsa.zinb1k <- loo(rmsa.glmer.zinb1k, save_psis = TRUE, moment_match = TRUE)
save(loo_rmsa.zinb1k, file="F:/_R_objects/loo_rmsa.zinb1k.rda")

loo_rmsa.nb1k <- loo(rmsa.glmer.nb1k, save_psis = TRUE, moment_match = TRUE)
save(loo_rmsa.nb1k, file="F:/_R_objects/loo_rmsa.nb1k.rda")

loo_rmsa.nb1k.2 <- loo(rmsa.glmer.nb1k.2, save_psis = TRUE, moment_match = TRUE)
save(loo_rmsa.nb1k.2, file="F:/_R_objects/loo_rmsa.nb1k.2.rda")


load(file="F:/_R_objects/loo_rmsa.zinb1k.rda")
load(file="F:/_R_objects/loo_rmsa.nb1k.rda")
load(file="F:/_R_objects/loo_rmsa.nb1k.2.rda")

loo_compare(loo_rmsa.zinb0k, loo_rmsa.zinb1k, loo_rmsa.nb0k, loo_rmsa.nb0k.2, loo_rmsa.nb1k, loo_rmsa.nb1k.2)


#===============================================================================
# Model 2 - PAR & Patch w/ interaction; compare likelihood specification
#===============================================================================

#------------------------------------------------
# Fit and compare models w/ ZINB & NB likelihoods
#------------------------------------------------


#zinb2k
rmsa.glmer.zinb2k <- brm(bf(PIUV_less5cmDBH_total ~ log_PAR1.3.rs + patch + log_PAR1.3.rs:patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b") + set_prior ("normal(0, 5)", class = "b", dpar="zi") + set_prior ("gamma(1, 0.1)", class = "shape"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb2k, file="F:/_R_objects/rmsa.glmer.zinb2k.rda")


#nb2k
rmsa.glmer.nb2k <- brm(bf(PIUV_less5cmDBH_total ~ log_PAR1.3.rs + patch + log_PAR1.3.rs:patch + (1 | plot)), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b") + set_prior ("gamma(1, 0.1)", class = "shape"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb2k, file="F:/_R_objects/rmsa.glmer.nb2k.rda")


#nb2k.2 [w/ shape modeled by patch]
rmsa.glmer.nb2k.2 <- brm(bf(PIUV_less5cmDBH_total ~ log_PAR1.3.rs + patch + log_PAR1.3.rs:patch + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb2k.2, file="F:/_R_objects/rmsa.glmer.nb2k.2.rda")


load(file="F:/_R_objects/rmsa.glmer.zinb0k.rda")
load(file="F:/_R_objects/rmsa.glmer.zinb2k.rda")
load(file="F:/_R_objects/rmsa.glmer.nb0k.rda")
load(file="F:/_R_objects/rmsa.glmer.nb0k.2.rda")
load(file="F:/_R_objects/rmsa.glmer.nb2k.2.rda")


#--------------
# Check models
#--------------

mod <- rmsa.glmer.nb2k.2
mod_name <- "rmsa.glmer.nb2k.2"

prior_summary(mod)
summary(mod)
variables(mod)


launch_shinystan(mod)

bayesplot_theme_set(ggplot2::theme_bw())
color_scheme_set("brightblue")

# density & trace plots
dev.new()
plot(mod)

dev.new()
mcmc_trace(mod)

# posterior vs. prior
variables(mod)
param <- "b_shape_Intercept"
prior_summary(mod)    # inspect prior for param
pdf <- function(x) dnorm(x, 0, 5) # density function for prior on param
pdf <- function(x) dstudent_t(x, 3, 0, 2.5)
pdf <- function(x) dgamma(x, 1, 0.1)
pdf <- function(x) dlogis(x, 0, 1)


hist(as.matrix(mod, variable = param), 25, prob = TRUE,
     col = "gray", border = "white", 
     las = 1, cex.axis = 1.2, cex.lab = 1.5, cex.main = 1.5,
     xlim = c(-10, 10),
     #ylim = c(0, .3),
     xlab = param, ylab = "Probability", main = paste("prior vs. posterior", mod_name))
curve(pdf(x), col = "blue", lwd = 2, add = TRUE)

exp(5)
log(10)
plogis(2.5) #inverse logit


# posterior predictive checks
pp_check(mod, ndraws = 250) + 
  xlab("Seedling counts") + ggtitle(mod_name)

y <- sp.rmb$PIUV_less5cmDBH_total # vector of outcome values
yrep <- posterior_predict(mod, draws = 500) # matrix of draws from the ppd

prop_zero <- function(x) mean(x == 0)

ppc_stat(y, yrep, stat = "prop_zero", binwidth = 0.005) +
  ggtitle(mod_name)
ppc_stat_grouped(y, yrep, sp.rmb$patch, stat = "prop_zero", binwidth = 0.005) +
  ggtitle(mod_name)
ppc_stat(y, yrep, stat = "max", binwidth = 5) + ggtitle(mod_name)
ppc_stat_grouped(y, yrep, sp.rmb$patch, stat = "max", binwidth = 20) + 
  ggtitle(mod_name)
ppc_stat(y, yrep, stat = "sd", binwidth = 10) +
  ggtitle(mod_name)
ppc_stat_grouped(y, yrep, sp.rmb$patch, stat = "sd", binwidth = 10) +
  ggtitle(mod_name)


# interval estimates
get_variables(mod)
dev.new()
mcmc_intervals(mod, prob = 0.90, prob_outer = 0.95,pars=c("b_Intercept", "b_zi_Intercept", "b_log_PAR1.3.rs", "b_patchBogForest", "b_patchBog", "b_log_PAR1.3.rs:patchBogForest", "b_log_PAR1.3.rs:patchBog", "b_zi_patchBogForest", "b_zi_patchBog", "sd_plot__Intercept", "shape")) + ggtitle(mod_name)

dev.new()
mcmc_intervals(mod, prob = 0.90, prob_outer = 0.95,pars=c("b_Intercept", "b_log_PAR1.3.rs", "b_patchBogForest", "b_patchBog", "b_log_PAR1.3.rs:patchBogForest", "b_log_PAR1.3.rs:patchBog", "sd_plot__Intercept", "b_shape_Intercept", "b_shape_patchBogForest", "b_shape_patchBog", "Intercept_shape")) + ggtitle(mod_name)

# pairwise plots
dev.new()
mcmc_pairs(as.matrix(mod),pars = c("b_Intercept", "b_zi_Intercept", "b_log_PAR1.3.rs", "b_patchBogForest", "b_patchBog", "b_log_PAR1.3.rs:patchBogForest", "b_log_PAR1.3.rs:patchBog", "b_zi_patchBogForest", "b_zi_patchBog", "sd_plot__Intercept", "shape")) + ggtitle(mod_name)

posterior <- as.matrix(mod)

dev.new()
mcmc_pairs(posterior, pars = c("b_Intercept", "b_zi_Intercept", "b_log_PAR1.3.rs", "b_patchBogForest", "b_patchBog", "b_log_PAR1.3.rs:patchBogForest", "b_log_PAR1.3.rs:patchBog", "b_zi_patchBogForest", "b_zi_patchBog", "sd_plot__Intercept", "shape")) + ggtitle(mod_name)


# conditional effects
plot(conditional_effects(mod), ask = FALSE) + ggtitle(mod_name)


#---------------------------
# Validate & compare models
#---------------------------

loo_rmsa.zinb2k <- loo(rmsa.glmer.zinb2k, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.zinb2k, file="F:/_R_objects/loo_rmsa.zinb2k.rda")
# one observation with pareto k > 0.7 [obs 29]

loo_rmsa.nb2k.2 <- loo(rmsa.glmer.nb2k.2, save_psis = TRUE, moment_match = TRUE)
save(loo_rmsa.nb2k.2, file="F:/_R_objects/loo_rmsa.nb2k.2.rda")

load(file="F:/_R_objects/loo_rmsa.zinb0k.rda")
load(file="F:/_R_objects/loo_rmsa.zinb2k.rda")
load(file="F:/_R_objects/loo_rmsa.nb0k.rda")
load(file="F:/_R_objects/loo_rmsa.nb0k.2.rda")
load(file="F:/_R_objects/loo_rmsa.nb2k.2.rda")

loo_compare(loo_rmsa.zinb0k, loo_rmsa.zinb2k,loo_rmsa.nb0k, loo_rmsa.nb0k.2, loo_rmsa.nb2k.2)


#################################################################################### FIT & EVALUATE MODELS 3 - x USING OPTIMAL LIKELIHOOD SPECIFICATION
################################################################################

#================================================
# Model 3 - DWT w/ and w/o patch and interaction
#================================================

#nb3k [w/ shape modeled by patch; change in object name format (no longer .2)]
rmsa.glmer.nb3k <- brm(bf(PIUV_less5cmDBH_total ~ DWT + patch + DWT:patch + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb3k, file="F:/_R_objects/rmsa.glmer.nb3k.rda")

#nb3k.b
rmsa.glmer.nb3k.b <- brm(bf(PIUV_less5cmDBH_total ~ DWT + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb3k.b, file="F:/_R_objects/rmsa.glmer.nb3k.b.rda")

#zinb3k
rmsa.glmer.zinb3k <- brm(bf(PIUV_less5cmDBH_total ~ DWT + patch + DWT:patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b") + set_prior ("normal(0, 5)", class = "b", dpar="zi") + set_prior ("gamma(1, 0.1)", class = "shape"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb3k, file="F:/_R_objects/rmsa.glmer.zinb3k.rda")

#===============================================
# Model 4 - pool & patch w/o and w/ interaction
#===============================================

#nb4k.a
rmsa.glmer.nb4k.a <- brm(bf(PIUV_less5cmDBH_total ~ Pool.rs + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb4k.a, file="F:/_R_objects/rmsa.glmer.nb4k.a.rda")

#nb4k.b
rmsa.glmer.nb4k.b <- brm(bf(PIUV_less5cmDBH_total ~ Pool.rs + patch + Pool.rs:patch + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb4k.b, file="F:/_R_objects/rmsa.glmer.nb4k.a.rda")

# from loo_compare, .b is almost significantly better at predicting previously unseen data than .a, but this is likely due to the inclusion of patch

#===============================================
# Model 5 -- PIUV cnt (seed trees)
#===============================================

#nb5k
rmsa.glmer.nb5k <- brm(bf(PIUV_less5cmDBH_total ~ Plot_PIUVcnt.rs + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb5k, file="F:/_R_objects/rmsa.glmer.nb5k.rda")

# median estimate differs from zero with 95% probability

#====================================================================
# Model 6 -- percent cover Sphagnum w/o and w/ interaction w/ patch
#====================================================================

#nb6k.a
rmsa.glmer.nb6k.a <- brm(bf(PIUV_less5cmDBH_total ~ Sphagnum.rs + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb6k.a, file="F:/_R_objects/rmsa.glmer.nb6k.a.rda")

#nb6k.b
rmsa.glmer.nb6k.b <- brm(bf(PIUV_less5cmDBH_total ~ Sphagnum.rs + patch + Sphagnum.rs:patch + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb6k.b, file="F:/_R_objects/rmsa.glmer.nb6k.b.rda")

#Sphagnum almost has parameter estimates outside the 95 CI in .a; .a has sig lower predictive density than .b

#nb6k.c -- patch included, but no interaction
rmsa.glmer.nb6k.c <- brm(bf(PIUV_less5cmDBH_total ~ Sphagnum.rs + patch + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb6k.c, file="F:/_R_objects/rmsa.glmer.nb6k.c.rda")


#===============================================
# Model 7 -- percent cover Cushion spp.
#===============================================

#nb7k
rmsa.glmer.nb7k <- brm(bf(PIUV_less5cmDBH_total ~ Cushion_spp.rs + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb7k, file="F:/_R_objects/rmsa.glmer.nb7k.rda")

#===============================================
# Model 8 -- percent cover Dicranaloma spp.
#===============================================

#nb8k
rmsa.glmer.nb8k <- brm(bf(PIUV_less5cmDBH_total ~ Dicran_spp.rs + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb8k, file="F:/_R_objects/rmsa.glmer.nb8k.rda")

# Dicran might have an estimate that differs from zero w/ 90% probability


#==========================================================
# Model 9 -- patch, percent cover Other Moss & interaction
#==========================================================

#nb9k
rmsa.glmer.nb9k <- brm(bf(PIUV_less5cmDBH_total ~ patch + OtherMoss_spp.rs + OtherMoss_spp.rs:patch + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb9k, file="F:/_R_objects/rmsa.glmer.nb9k.rda")

#=====================================================================
# Model 10 -- percent cover Other Moss spp., DWT and their interaction
#=====================================================================

#nb10k
rmsa.glmer.nb10k <- brm(bf(PIUV_less5cmDBH_total ~ DWT + OtherMoss_spp.rs + OtherMoss_spp.rs:DWT + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb10k, file="F:/_R_objects/rmsa.glmer.nb10k.rda")

#=========================================================================
# Model 11 -- percent cover understory PIUV w/o & w/ interaction w/ patch
#=========================================================================

#nb11k
rmsa.glmer.nb11k <- brm(bf(PIUV_less5cmDBH_total ~ US_PIUV.rs + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb11k, file="F:/_R_objects/rmsa.glmer.nb11k.rda")

#US_PIUV.rs may have a positive effect that differs from zero with 90% probability

#nb11k.b
rmsa.glmer.nb11k.b <- brm(bf(PIUV_less5cmDBH_total ~ US_PIUV.rs + patch + US_PIUV.rs:patch + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb11k.b, file="F:/_R_objects/rmsa.glmer.nb11k.b.rda")

#no apparent effect of interaction

#============================================================================
# Model 12 -- percent cover understory non-PIUV trees w/o and w/ patch 
#============================================================================

#nb12k
rmsa.glmer.nb12k <- brm(bf(PIUV_less5cmDBH_total ~ US_Trees.rs + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb12k, file="F:/_R_objects/rmsa.glmer.nb12k.rda")

#US_Trees.rs may have a negative effect that differs from zero with 90% probability

#nb12k.b -- interaction w/ patch
rmsa.glmer.nb12k.b <- brm(bf(PIUV_less5cmDBH_total ~ US_Trees.rs + patch + US_Trees.rs:patch + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb12k.b, file="F:/_R_objects/rmsa.glmer.nb12k.b.rda")

#nb12k.c -- no interaction w/ patch
rmsa.glmer.nb12k.c <- brm(bf(PIUV_less5cmDBH_total ~ US_Trees.rs + patch + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb12k.c, file="F:/_R_objects/rmsa.glmer.nb12k.c.rda")

#============================================================================
# Model 13 -- percent cover Other Shrubs w/o and w/ interaction w/ patch
#============================================================================

#nb13k.a
rmsa.glmer.nb13k.a <- brm(bf(PIUV_less5cmDBH_total ~ OtherShrubs.rs + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb13k.a, file="F:/_R_objects/rmsa.glmer.nb13k.a.rda")

#no effect of other shrubs on seedling counts

#nb13k.b
rmsa.glmer.nb13k.b <- brm(bf(PIUV_less5cmDBH_total ~ OtherShrubs.rs + patch + OtherShrubs.rs:patch + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb13k.b, file="F:/_R_objects/rmsa.glmer.nb13k.b.rda")

#negative effect of shrubs on seedlings counts in the Forest patch, but the patchBog estimate includes zero and the CI is very wide leading to unreasonable counts.

#nb13k.c -- model shape as a function of patch, OtherShrubs & interaction
rmsa.glmer.nb13k.c <- brm(bf(PIUV_less5cmDBH_total ~ OtherShrubs.rs + patch + OtherShrubs.rs:patch + (1 | plot), shape ~ patch + OtherShrubs.rs + OtherShrubs.rs:patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb13k.c, file="F:/_R_objects/rmsa.glmer.nb13k.c.rda")
#this took a long time to run & there were many errors

#nb13k.d -- zinb version
rmsa.glmer.zinb13k.b <- brm(bf(PIUV_less5cmDBH_total ~ OtherShrubs.rs + patch + OtherShrubs.rs:patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b") + set_prior ("normal(0, 5)", class = "b", dpar="zi") + set_prior ("gamma(1, 0.1)", class = "shape"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb13k.b, file="F:/_R_objects/rmsa.glmer.zinb13k.b.rda")

#============================================================================
# Model 14 -- percent cover Low Shrubs w/o and w/ interaction w/ patch
#============================================================================

#nb14k.a
rmsa.glmer.nb14k.a <- brm(bf(PIUV_less5cmDBH_total ~ LowShrubs.rs + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb14k.a, file="F:/_R_objects/rmsa.glmer.nb14k.a.rda")

#nb14k.b
rmsa.glmer.nb14k.b <- brm(bf(PIUV_less5cmDBH_total ~ LowShrubs.rs + patch + LowShrubs.rs:patch + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb14k.b, file="F:/_R_objects/rmsa.glmer.nb14k.b.rda")

#============================================================================
# Model 15 -- percent cover Ferns w/o and w/ interaction w/ patch
#============================================================================

#nb15k.a
rmsa.glmer.nb15k.a <- brm(bf(PIUV_less5cmDBH_total ~ Ferns.rs + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb15k.a, file="F:/_R_objects/rmsa.glmer.nb15k.a.rda")

#nb15k.b
rmsa.glmer.nb15k.b <- brm(bf(PIUV_less5cmDBH_total ~ Ferns.rs + patch + Ferns.rs:patch + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb15k.b, file="F:/_R_objects/rmsa.glmer.nb15k.b.rda")
# some pareto k diagnostic values too high


#============================================================================
# Model 16 -- percent cover Forbs w/o and w/ interaction w/ patch
#============================================================================

#nb16k.a
rmsa.glmer.nb16k.a <- brm(bf(PIUV_less5cmDBH_total ~ Forbs.rs + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb16k.a, file="F:/_R_objects/rmsa.glmer.nb16k.a.rda")

#nb16k.b
rmsa.glmer.nb16k.b <- brm(bf(PIUV_less5cmDBH_total ~ Forbs.rs + patch + Forbs.rs:patch + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb16k.b, file="F:/_R_objects/rmsa.glmer.nb16k.b.rda")


#============================================================================
# Model 17 -- Number stems other tree spp. (height > 130 cm)
#============================================================================

#nb17k.a
rmsa.glmer.nb17k.a <- brm(bf(PIUV_less5cmDBH_total ~ OtherTrees.rs + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb17k.a, file="F:/_R_objects/rmsa.glmer.nb17k.a.rda")
#1 divergent transition

#nb17k.b
rmsa.glmer.nb17k.b <- brm(bf(PIUV_less5cmDBH_total ~ OtherTrees.rs + patch + OtherTrees.rs:patch + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb17k.b, file="F:/_R_objects/rmsa.glmer.nb17k.b.rda")


#### Global Models
#============================================================================
# Model 18
#============================================================================

#nb18k.a
rmsa.glmer.nb18k.a <- brm(bf(PIUV_less5cmDBH_total ~ patch + DWT + log_PAR1.3.rs + Plot_PIUVcnt.rs + Sphagnum.rs + Dicran_spp.rs + DWT:patch + log_PAR1.3.rs:patch + Sphagnum.rs:patch + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb18k.a, file="F:/_R_objects/rmsa.glmer.nb18k.a.rda")

#nb18k.b
rmsa.glmer.nb18k.b <- brm(bf(PIUV_less5cmDBH_total ~ patch + DWT + log_PAR1.3.rs + Plot_PIUVcnt.rs + Sphagnum.rs + Dicran_spp.rs + DWT:patch + log_PAR1.3.rs:patch + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb18k.b, file="F:/_R_objects/rmsa.glmer.nb18k.b.rda")

#nb18k.c
rmsa.glmer.nb18k.c <- brm(bf(PIUV_less5cmDBH_total ~ patch + DWT + log_PAR1.3.rs + Plot_PIUVcnt.rs + Dicran_spp.rs + DWT:patch + log_PAR1.3.rs:patch + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb18k.c, file="F:/_R_objects/rmsa.glmer.nb18k.c.rda")

#nb18k.d
rmsa.glmer.nb18k.d <- brm(bf(PIUV_less5cmDBH_total ~ DWT + log_PAR1.3.rs + Plot_PIUVcnt.rs + patch + DWT:patch + log_PAR1.3.rs:patch + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb18k.d, file="F:/_R_objects/rmsa.glmer.nb18k.d.rda")

#nb18k.e
rmsa.glmer.nb18k.e <- brm(bf(PIUV_less5cmDBH_total ~ patch + DWT + Plot_PIUVcnt.rs + DWT:patch + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb18k.e, file="F:/_R_objects/rmsa.glmer.nb18k.e.rda")

#nb18k.f
rmsa.glmer.nb18k.f <- brm(bf(PIUV_less5cmDBH_total ~ patch + log_PAR1.3.rs + Plot_PIUVcnt.rs + log_PAR1.3.rs:patch + (1 | plot), shape ~ patch), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(2, 3)", class = "Intercept") + set_prior("normal(0, 5)", class = "b"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb18k.f, file="F:/_R_objects/rmsa.glmer.nb18k.f.rda")


# ppc ok
load(file="F:/_R_objects/rmsa.glmer.nb1k.2.rda")
load(file="F:/_R_objects/rmsa.glmer.nb2k.2.rda")
load(file="F:/_R_objects/rmsa.glmer.nb6k.b.rda")
load(file="F:/_R_objects/rmsa.glmer.nb9k.rda")

# ppc wonky
load(file="F:/_R_objects/rmsa.glmer.nb0k.2.rda")
load(file="F:/_R_objects/rmsa.glmer.nb3k.rda")
load(file="F:/_R_objects/rmsa.glmer.nb4k.a.rda")
load(file="F:/_R_objects/rmsa.glmer.nb5k.rda")
load(file="F:/_R_objects/rmsa.glmer.nb6k.a.rda")
load(file="F:/_R_objects/rmsa.glmer.nb7k.rda")
load(file="F:/_R_objects/rmsa.glmer.nb8k.rda")
load(file="F:/_R_objects/rmsa.glmer.nb12k.b.rda")
load(file="F:/_R_objects/rmsa.glmer.nb13k.b.rda")
load(file="F:/_R_objects/rmsa.glmer.nb14k.b.rda")
load(file="F:/_R_objects/rmsa.glmer.nb15k.b.rda")
load(file="F:/_R_objects/rmsa.glmer.nb16k.b.rda")
load(file="F:/_R_objects/rmsa.glmer.nb18k.a.rda")
load(file="F:/_R_objects/rmsa.glmer.nb18k.b.rda")
load(file="F:/_R_objects/rmsa.glmer.nb18k.c.rda")
load(file="F:/_R_objects/rmsa.glmer.nb18k.d.rda")


#scrap
load(file="F:/_R_objects/rmsa.glmer.nb3k.b.rda")
load(file="F:/_R_objects/rmsa.glmer.nb4k.b.rda")
load(file="F:/_R_objects/rmsa.glmer.nb6k.c.rda")
load(file="F:/_R_objects/rmsa.glmer.nb10k.rda")
load(file="F:/_R_objects/rmsa.glmer.nb11k.rda")
load(file="F:/_R_objects/rmsa.glmer.nb11k.b.rda")
load(file="F:/_R_objects/rmsa.glmer.nb12k.rda")
load(file="F:/_R_objects/rmsa.glmer.nb12k.c.rda")
load(file="F:/_R_objects/rmsa.glmer.nb13k.a.rda")
load(file="F:/_R_objects/rmsa.glmer.zinb13k.b.rda")
load(file="F:/_R_objects/rmsa.glmer.nb14k.a.rda")
load(file="F:/_R_objects/rmsa.glmer.nb15k.a.rda")
load(file="F:/_R_objects/rmsa.glmer.nb16k.a.rda")
load(file="F:/_R_objects/rmsa.glmer.nb17k.a.rda")
load(file="F:/_R_objects/rmsa.glmer.nb17k.b.rda")
load(file="F:/_R_objects/rmsa.glmer.nb18k.e.rda")
load(file="F:/_R_objects/rmsa.glmer.nb18k.f.rda")




#-------------------------------------------------------------------------------
# Check models
#-------------------------------------------------------------------------------

mod <- rmsa.glmer.nb9k
mod_name <- "rmsa.glmer.nb9k"

prior_summary(mod)
summary(mod)
variables(mod)


launch_shinystan(mod)

bayesplot_theme_set(ggplot2::theme_bw())
color_scheme_set("brightblue")

# density & trace plots
dev.new()
plot(mod)

dev.new()
mcmc_trace(mod)

# posterior vs. prior
variables(mod)
param <- "b_shape_patchBogForest"
prior_summary(mod)    # inspect prior for param
pdf <- function(x) dnorm(x, 0, 5) # density function for prior on param
pdf <- function(x) dstudent_t(x, 3, 0, 2.5)
pdf <- function(x) dgamma(x, 1, 0.1)
pdf <- function(x) dlogis(x, 0, 1)


hist(as.matrix(mod, variable = param), 25, prob = TRUE,
     col = "gray", border = "white", 
     las = 1, cex.axis = 1.2, cex.lab = 1.5, cex.main = 1.5,
     xlim = c(-10, 10),
     #ylim = c(0, .3),
     xlab = param, ylab = "Probability", main = paste("prior vs. posterior", mod_name))
curve(pdf(x), col = "blue", lwd = 2, add = TRUE)

exp(5)
log(10)
plogis(2.5) #inverse logit


# posterior predictive checks
y <- sp.rmb$PIUV_less5cmDBH_total # vector of outcome values
yrep <- posterior_predict(mod, draws = 500) # matrix of draws from the ppd

prop_zero <- function(x) mean(x == 0)

pp_check(mod, ndraws = 250) + 
  xlab("Seedling counts") + ggtitle(mod_name)

ppc_stat_grouped(y, yrep, sp.rmb$patch, stat = "prop_zero", binwidth = 0.005) +
  ggtitle(mod_name)
ppc_stat_grouped(y, yrep, sp.rmb$patch, stat = "max", binwidth = 20) + 
  ggtitle(mod_name) #+ coord_cartesian(xlim = c(0, 1000))
ppc_stat_grouped(y, yrep, sp.rmb$patch, stat = "sd", binwidth = 10) +
  ggtitle(mod_name) + coord_cartesian(xlim = c(0, 100))

ppc_stat(y, yrep, stat = "prop_zero", binwidth = 0.005) +
  ggtitle
ppc_stat(y, yrep, stat = "max", binwidth = 5) + ggtitle
ppc_stat(y, yrep, stat = "sd", binwidth = 5) +
  ggtitle(mod_name)


# interval estimates
get_variables(mod)
dev.new()
mcmc_intervals(mod, prob = 0.90, prob_outer = 0.95,pars=c("b_Intercept", "b_patchBogForest", "b_patchBog", "sd_plot__Intercept", "b_shape_Intercept", "b_shape_patchBogForest", "b_shape_patchBog")) + ggtitle(mod_name)

get_variables(mod)
dev.new()
mcmc_intervals(mod, prob = 0.90, prob_outer = 0.95,pars=c("b_Intercept", "b_OtherShrubs.rs", "b_patchBogForest", "b_patchBog", "b_OtherShrubs.rs:patchBogForest", "b_OtherShrubs.rs:patchBog", "sd_plot__Intercept", "b_shape_Intercept", "b_shape_patchBogForest", "b_shape_patchBog")) + ggtitle(mod_name)

get_variables(mod)
dev.new()
mcmc_intervals(mod, prob = 0.90, prob_outer = 0.95,pars=c("b_Intercept", "b_OtherShrubs.rs", "b_patchBogForest", "b_patchBog", "b_OtherShrubs.rs:patchBogForest", "b_OtherShrubs.rs:patchBog", "sd_plot__Intercept", "b_shape_Intercept", "b_shape_patchBogForest", "b_shape_patchBog")) + ggtitle(mod_name)

# pairwise plots
dev.new()
mcmc_pairs(as.matrix(mod),pars = c("b_Intercept", "b_zi_Intercept", "b_log_PAR1.3.rs", "b_patchBogForest", "b_patchBog", "b_log_PAR1.3.rs:patchBogForest", "b_log_PAR1.3.rs:patchBog", "b_zi_patchBogForest", "b_zi_patchBog", "sd_plot__Intercept", "shape")) + ggtitle(mod_name)


# conditional effects
plot(conditional_effects(mod), ask = FALSE) + ggtitle(mod_name)


#######
#point summaries and intervals for posterior draws of models
get_variables(rmsa.glmer.nb16k.a)
modsum <- rmsa.glmer.nb16k.a %>%
  gather_draws(b_Intercept,
               b_shape_Intercept,
               b_Forbs.rs,
               b_shape_patchBogForest, b_shape_patchBog,
               sd_plot__Intercept) %>%
  median_qi(na.rm = TRUE, .width = c(.95, .9, .8))

write.table(modsum, "clipboard", sep="\t", row.names=TRUE)
View(modsum) 

get_variables(rmsa.glmer.nb18k.f)
modsum <- rmsa.glmer.nb18k.f %>%
  gather_draws(b_Intercept,
               b_shape_Intercept,
               b_patchBogForest, b_patchBog,
               b_log_PAR1.3.rs, b_Plot_PIUVcnt.rs,
               `b_patchBogForest:log_PAR1.3.rs`, `b_patchBog:log_PAR1.3.rs`,
               b_shape_patchBogForest, b_shape_patchBog,
               sd_plot__Intercept) %>%
  median_qi(na.rm = TRUE, .width = c(.95, .9, .8))

write.table(modsum, "clipboard", sep="\t", row.names=TRUE)
View(modsum)


summary(rmsa.glmer.nb14k.b)

#---------------------------
# Validate & compare models
#---------------------------

loo_rmsa.nb3k <- loo(rmsa.glmer.nb3k, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb3k, file="F:/_R_objects/loo_rmsa.nb3k.rda")

loo_rmsa.nb3k.b <- loo(rmsa.glmer.nb3k.b, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb3k.b, file="F:/_R_objects/loo_rmsa.nb3k.b.rda")

loo_rmsa.nb4k.a <- loo(rmsa.glmer.nb4k.a, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb4k.a, file="F:/_R_objects/loo_rmsa.nb4k.a.rda")

loo_rmsa.nb4k.b <- loo(rmsa.glmer.nb4k.b, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb4k.b, file="F:/_R_objects/loo_rmsa.nb4k.b.rda")

loo_rmsa.nb5k <- loo(rmsa.glmer.nb5k, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb5k, file="F:/_R_objects/loo_rmsa.nb5k.rda")

loo_rmsa.nb6k.a <- loo(rmsa.glmer.nb6k.a, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb6k.a, file="F:/_R_objects/loo_rmsa.nb6k.a.rda")

loo_rmsa.nb6k.b <- loo(rmsa.glmer.nb6k.b, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb6k.b, file="F:/_R_objects/loo_rmsa.nb6k.b.rda")

loo_rmsa.nb6k.c <- loo(rmsa.glmer.nb6k.c, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb6k.c, file="F:/_R_objects/loo_rmsa.nb6k.c.rda")

loo_rmsa.nb7k <- loo(rmsa.glmer.nb7k, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb7k, file="F:/_R_objects/loo_rmsa.nb7k.rda")

loo_rmsa.nb8k <- loo(rmsa.glmer.nb8k, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb8k, file="F:/_R_objects/loo_rmsa.nb8k.rda")

loo_rmsa.nb9k <- loo(rmsa.glmer.nb9k, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb9k, file="F:/_R_objects/loo_rmsa.nb9k.rda")

loo_rmsa.nb10k <- loo(rmsa.glmer.nb10k, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb10k, file="F:/_R_objects/loo_rmsa.nb10k.rda")

loo_rmsa.nb11k <- loo(rmsa.glmer.nb11k, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb11k, file="F:/_R_objects/loo_rmsa.nb11k.rda")

loo_rmsa.nb11k.b <- loo(rmsa.glmer.nb11k.b, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb11k.b, file="F:/_R_objects/loo_rmsa.nb11k.b.rda")

loo_rmsa.nb12k <- loo(rmsa.glmer.nb12k, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb12k, file="F:/_R_objects/loo_rmsa.nb12k.rda")

loo_rmsa.nb12k.b <- loo(rmsa.glmer.nb12k.b, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb12k.b, file="F:/_R_objects/loo_rmsa.nb12k.b.rda")

loo_rmsa.nb12k.c <- loo(rmsa.glmer.nb12k.c, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb12k.c, file="F:/_R_objects/loo_rmsa.nb12k.c.rda")

loo_rmsa.nb13k.a <- loo(rmsa.glmer.nb13k.a, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb13k.a, file="F:/_R_objects/loo_rmsa.nb13k.a.rda")

loo_rmsa.nb13k.b <- loo(rmsa.glmer.nb13k.b, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb13k.b, file="F:/_R_objects/loo_rmsa.nb13k.b.rda")

loo_rmsa.nb14k.a <- loo(rmsa.glmer.nb14k.a, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb14k.a, file="F:/_R_objects/loo_rmsa.nb14k.a.rda")

loo_rmsa.nb14k.b <- loo(rmsa.glmer.nb14k.b, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb14k.b, file="F:/_R_objects/loo_rmsa.nb14k.b.rda")

loo_rmsa.nb15k.a <- loo(rmsa.glmer.nb15k.a, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb15k.a, file="F:/_R_objects/loo_rmsa.nb15k.a.rda")

loo_rmsa.nb15k.b <- loo(rmsa.glmer.nb15k.b, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb15k.b, file="F:/_R_objects/loo_rmsa.nb15k.b.rda")
# one problematic obs; some  pareto k diagnostic values too high

loo_rmsa.nb16k.a <- loo(rmsa.glmer.nb16k.a, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb16k.a, file="F:/_R_objects/loo_rmsa.nb16k.a.rda")

loo_rmsa.nb16k.b <- loo(rmsa.glmer.nb16k.b, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb16k.b, file="F:/_R_objects/loo_rmsa.nb16k.b.rda")


loo_rmsa.nb18k.a <- loo(rmsa.glmer.nb18k.a, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb18k.a, file="F:/_R_objects/loo_rmsa.nb18k.a.rda")

loo_rmsa.nb18k.b <- loo(rmsa.glmer.nb18k.b, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb18k.b, file="F:/_R_objects/loo_rmsa.nb18k.b.rda")

loo_rmsa.nb18k.c <- loo(rmsa.glmer.nb18k.c, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb18k.c, file="F:/_R_objects/loo_rmsa.nb18k.c.rda")

loo_rmsa.nb18k.d <- loo(rmsa.glmer.nb18k.d, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb18k.d, file="F:/_R_objects/loo_rmsa.nb18k.d.rda")

loo_rmsa.nb18k.e <- loo(rmsa.glmer.nb18k.e, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb18k.e, file="F:/_R_objects/loo_rmsa.nb18k.e.rda")

loo_rmsa.nb18k.f <- loo(rmsa.glmer.nb18k.f, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.nb18k.f, file="F:/_R_objects/loo_rmsa.nb18k.f.rda")

loo_rmsa.nb0k.2
load(file="F:/_R_objects/loo_rmsa.nb0k.2.rda")
load(file="F:/_R_objects/loo_rmsa.nb1k.2.rda")
load(file="F:/_R_objects/loo_rmsa.nb2k.2.rda")
load(file="F:/_R_objects/loo_rmsa.nb3k.rda")
load(file="F:/_R_objects/loo_rmsa.nb3k.b.rda")
load(file="F:/_R_objects/loo_rmsa.nb4k.a.rda")
load(file="F:/_R_objects/loo_rmsa.nb4k.b.rda")
load(file="F:/_R_objects/loo_rmsa.nb5k.rda")
load(file="F:/_R_objects/loo_rmsa.nb6k.a.rda")
load(file="F:/_R_objects/loo_rmsa.nb6k.b.rda")
load(file="F:/_R_objects/loo_rmsa.nb7k.rda")
load(file="F:/_R_objects/loo_rmsa.nb8k.rda")
load(file="F:/_R_objects/loo_rmsa.nb9k.rda")
load(file="F:/_R_objects/loo_rmsa.nb10k.rda")
load(file="F:/_R_objects/loo_rmsa.nb11k.rda")
load(file="F:/_R_objects/loo_rmsa.nb11k.b.rda")
load(file="F:/_R_objects/loo_rmsa.nb12k.rda")
load(file="F:/_R_objects/loo_rmsa.nb12k.b.rda")
load(file="F:/_R_objects/loo_rmsa.nb12k.c.rda")
load(file="F:/_R_objects/loo_rmsa.nb13k.a.rda")
load(file="F:/_R_objects/loo_rmsa.nb13k.b.rda")
load(file="F:/_R_objects/loo_rmsa.nb14k.a.rda")
load(file="F:/_R_objects/loo_rmsa.nb14k.b.rda")
load(file="F:/_R_objects/loo_rmsa.nb15k.a.rda")
load(file="F:/_R_objects/loo_rmsa.nb15k.b.rda")
load(file="F:/_R_objects/loo_rmsa.nb16k.a.rda")
load(file="F:/_R_objects/loo_rmsa.nb16k.b.rda")
load(file="F:/_R_objects/loo_rmsa.nb18k.a.rda")
load(file="F:/_R_objects/loo_rmsa.nb18k.b.rda")
load(file="F:/_R_objects/loo_rmsa.nb18k.c.rda")
load(file="F:/_R_objects/loo_rmsa.nb18k.d.rda")
load(file="F:/_R_objects/loo_rmsa.nb18k.e.rda")
load(file="F:/_R_objects/loo_rmsa.nb18k.f.rda")


# compare all models w/ hypotheses but w/o problematic results
loo_compare(loo_rmsa.nb0k.2, loo_rmsa.nb1k.2, loo_rmsa.nb2k.2, loo_rmsa.nb3k, loo_rmsa.nb4k.a,  loo_rmsa.nb5k, loo_rmsa.nb6k.a, loo_rmsa.nb6k.b, loo_rmsa.nb7k, loo_rmsa.nb8k, loo_rmsa.nb9k, loo_rmsa.nb12k.b, loo_rmsa.nb18k.a, loo_rmsa.nb18k.b, loo_rmsa.nb18k.c, loo_rmsa.nb18k.d)

# compare select models
loo_compare(loo_rmsa.nb0k.2, loo_rmsa.nb1k.2, loo_rmsa.nb2k.2, loo_rmsa.nb3k, loo_rmsa.nb5k, loo_rmsa.nb6k.b, loo_rmsa.nb8k, loo_rmsa.nb18k.a, loo_rmsa.nb18k.b, loo_rmsa.nb18k.c, loo_rmsa.nb18k.d)


# compare select models
loo_compare(loo_rmsa.nb0k.2, loo_rmsa.nb1k.2, loo_rmsa.nb18k.a, loo_rmsa.nb18k.b, loo_rmsa.nb18k.c, loo_rmsa.nb18k.d)



#################################################################################### CONSTRUCT & PLOT A FINAL STACKED MODEL
################################################################################

#===============================================================================
# Determine model averaging weights using Bayesian stacking
#===============================================================================

# When we have more than two models, it can be easier to understand the uncertainties in the comparison by looking at the model averaging weights based on Bayesian stacking (Yao et al., 2018) which optimizes the model weights so that the combined predictive distribution maximizes the estimate leave-one-out cross-validation performance. https://avehtari.github.io/modelselection/rats_kcv.html#4_Leave-one-out_cross-validation

# all models
lpd_point <- cbind(
  loo_rmsa.nb0k.2$pointwise[,"elpd_loo"],
  loo_rmsa.nb1k.2$pointwise[,"elpd_loo"],
  loo_rmsa.nb2k.2$pointwise[,"elpd_loo"],
  loo_rmsa.nb3k$pointwise[,"elpd_loo"],
  loo_rmsa.nb4k.a$pointwise[,"elpd_loo"],
  loo_rmsa.nb5k$pointwise[,"elpd_loo"],
  loo_rmsa.nb6k.a$pointwise[,"elpd_loo"],
  loo_rmsa.nb6k.b$pointwise[,"elpd_loo"],
  loo_rmsa.nb7k$pointwise[,"elpd_loo"],
  loo_rmsa.nb8k$pointwise[,"elpd_loo"],
  loo_rmsa.nb9k$pointwise[,"elpd_loo"],
  loo_rmsa.nb12k.b$pointwise[,"elpd_loo"],
  loo_rmsa.nb18k.a$pointwise[,"elpd_loo"],
  loo_rmsa.nb18k.b$pointwise[,"elpd_loo"],
  loo_rmsa.nb18k.c$pointwise[,"elpd_loo"],
  loo_rmsa.nb18k.d$pointwise[,"elpd_loo"])
stacking_weights(lpd_point)

# select models
lpd_point <- cbind(
  loo_rmsa.nb0k.2$pointwise[,"elpd_loo"],
  loo_rmsa.nb1k.2$pointwise[,"elpd_loo"],
  loo_rmsa.nb2k.2$pointwise[,"elpd_loo"],
  loo_rmsa.nb3k$pointwise[,"elpd_loo"],
  loo_rmsa.nb5k$pointwise[,"elpd_loo"],
  loo_rmsa.nb6k.a$pointwise[,"elpd_loo"],
  loo_rmsa.nb6k.b$pointwise[,"elpd_loo"],
  loo_rmsa.nb8k$pointwise[,"elpd_loo"],
  loo_rmsa.nb18k.a$pointwise[,"elpd_loo"],
  loo_rmsa.nb18k.b$pointwise[,"elpd_loo"],
  loo_rmsa.nb18k.c$pointwise[,"elpd_loo"],
  loo_rmsa.nb18k.d$pointwise[,"elpd_loo"])
stacking_weights(lpd_point)


?stacking_weights
#===============================================================================
# ADD POSTERIOR DRAWS FROM MODELS (BASED ON STACKING WEIGHTS) TO SINGLE TIBBLE
#===============================================================================

load(file="F:/_R_objects/rmsa.glmer.nb1k.2.rda")
load(file="F:/_R_objects/rmsa.glmer.nb3k.rda")
load(file="F:/_R_objects/rmsa.glmer.nb5k.rda")
load(file="F:/_R_objects/rmsa.glmer.nb6k.b.rda")
load(file="F:/_R_objects/rmsa.glmer.nb8k.rda")
load(file="F:/_R_objects/rmsa.glmer.nb12k.b.rda")
load(file="F:/_R_objects/rmsa.glmer.nb18k.d.rda")



# extract draws from posterior distributions of each model; set number of draws to stacking weight * 10000
mod <- rmsa.glmer.nb18k.d
get_variables(mod) #note:  will need backticks (``) around names that use colons

post.mod1 <- rmsa.glmer.nb1k.2 %>%
  spread_draws(b_Intercept, b_shape_Intercept,
               b_patchBogForest, b_patchBog,
               b_shape_patchBogForest, b_shape_patchBog,
               sd_plot__Intercept,
               ndraws=4120)

post.mod3 <- rmsa.glmer.nb3k %>%
  spread_draws(b_Intercept, b_shape_Intercept,
               b_patchBogForest, b_patchBog,
               b_DWT.L, b_DWT.Q,
               `b_DWT.L:patchBogForest`, `b_DWT.L:patchBog`,
               `b_DWT.Q:patchBogForest`, `b_DWT.Q:patchBog`,
               b_shape_patchBogForest, b_shape_patchBog,
               sd_plot__Intercept,
               ndraws=210)

post.mod5 <- rmsa.glmer.nb5k %>%
  spread_draws(b_Intercept, b_shape_Intercept,
               b_Plot_PIUVcnt.rs,
               b_shape_patchBogForest, b_shape_patchBog,
               sd_plot__Intercept,
               ndraws=320)

post.mod6b <- rmsa.glmer.nb6k.b %>%
  spread_draws(b_Intercept, b_shape_Intercept,
               b_patchBogForest, b_patchBog,
               b_Sphagnum.rs,
               `b_Sphagnum.rs:patchBogForest`, `b_Sphagnum.rs:patchBog`,
               b_shape_patchBogForest, b_shape_patchBog,
               sd_plot__Intercept,
               ndraws=10)

post.mod8 <- rmsa.glmer.nb8k %>%
  spread_draws(b_Intercept, b_shape_Intercept,
               b_Dicran_spp.rs,
               b_shape_patchBogForest, b_shape_patchBog,
               sd_plot__Intercept,
               ndraws=1100)

post.mod12b <- rmsa.glmer.nb12k.b %>%
  spread_draws(b_Intercept, b_shape_Intercept,
               b_patchBogForest, b_patchBog,
               b_US_Trees.rs,
               `b_US_Trees.rs:patchBogForest`, `b_US_Trees.rs:patchBog`,
               b_shape_patchBogForest, b_shape_patchBog,
               sd_plot__Intercept,
               ndraws=830)

post.mod18d <- rmsa.glmer.nb18k.d %>%
  spread_draws(b_Intercept, b_shape_Intercept,
               b_patchBogForest, b_patchBog,
               b_DWT.L, b_DWT.Q,
               b_log_PAR1.3.rs, b_Plot_PIUVcnt.rs,
               `b_DWT.L:patchBogForest`, `b_DWT.L:patchBog`,
               `b_DWT.Q:patchBogForest`, `b_DWT.Q:patchBog`,
               `b_log_PAR1.3.rs:patchBogForest`, `b_log_PAR1.3.rs:patchBog`,
               b_shape_patchBogForest, b_shape_patchBog,
               sd_plot__Intercept,
               ndraws=3410)


glimpse(post.mod18d)


# bind posterior draws from all models into a single tibble & tidy

post.rmb.stack <- bind_rows(post.mod1, post.mod3, post.mod5, post.mod6b, post.mod8, post.mod12b, post.mod18d)

post.rmb.stack <- post.rmb.stack %>%
  mutate(draw_orig = .draw) %>% # save original draw numbers just in case
  mutate(.draw = 1:10000) # each .draw needs to be unique

glimpse(post.rmb.stack)

save(post.rmb.stack, file="F:/_R_objects/post.rmb.stack.rda")
load(file="F:/_R_objects/post.rmb.stack.rda")


#===============================================================================
# SUMMARIZE & VISUALIZE FINAL STACKED MODEL
#===============================================================================

get_variables(post.rmb.stack)


#### point summaries and intervals for posterior draws of stacked models
post.rmb.stack.sum <- post.rmb.stack %>%
  gather_draws(b_Intercept, b_shape_Intercept,
               b_patchBogForest, b_patchBog,
               b_DWT.L, b_DWT.Q,
               b_log_PAR1.3.rs, b_Plot_PIUVcnt.rs,
               b_Sphagnum.rs, b_Dicran_spp.rs,
               b_US_Trees.rs,
               `b_DWT.L:patchBogForest`, `b_DWT.L:patchBog`,
               `b_DWT.Q:patchBogForest`, `b_DWT.Q:patchBog`,
               `b_log_PAR1.3.rs:patchBogForest`, `b_log_PAR1.3.rs:patchBog`,
               `b_Sphagnum.rs:patchBogForest`, `b_Sphagnum.rs:patchBog`,
               `b_US_Trees.rs:patchBogForest`, `b_US_Trees.rs:patchBog`,
               b_shape_patchBogForest, b_shape_patchBog,
               sd_plot__Intercept,
               ) %>%
  mutate(.variable = str_remove(.variable, "b_")) %>% 
  mutate(.variable = str_remove(.variable, ".rs")) %>%
  mutate(.variable = str_remove(.variable, "patch")) %>%
  median_qi(na.rm = TRUE, .width = c(.95, .9, .8))


write.table(post.rmb.stack.sum, "clipboard", sep="\t", row.names=TRUE)
View(post.rmb.stack.sum) 

?median_qi

#----------
# plots
#----------

#point intervals of parameter estimates from posterior distribution

post.rmb.stack %>%
  gather_draws(b_Intercept, b_shape_Intercept,
               b_patchBogForest, b_patchBog,
               b_DWT.L, b_DWT.Q,
               b_log_PAR1.3.rs, b_Plot_PIUVcnt.rs,
               b_Sphagnum.rs, b_Dicran_spp.rs,
               b_US_Trees.rs,
               `b_DWT.L:patchBogForest`, `b_DWT.L:patchBog`,
               `b_DWT.Q:patchBogForest`, `b_DWT.Q:patchBog`,
               `b_log_PAR1.3.rs:patchBogForest`, `b_log_PAR1.3.rs:patchBog`,
               `b_Sphagnum.rs:patchBogForest`, `b_Sphagnum.rs:patchBog`,
               `b_US_Trees.rs:patchBogForest`, `b_US_Trees.rs:patchBog`,
               b_shape_patchBogForest, b_shape_patchBog,
               sd_plot__Intercept) %>%
  median_qi(.width = c(.95, .8), na.rm = TRUE) %>%
  mutate(.variable = str_remove(.variable, "b_")) %>% 
  mutate(.variable = str_remove(.variable, ".rs")) %>%
  mutate(.variable = str_remove(.variable, "patch")) %>%
  ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper)) +
  geom_pointinterval(color="purple") +
  scale_y_discrete(limits = c("sd_plot__Intercept", "shape_Bog",
                              "shape_BogForest", "shape_Intercept",
                              "US_Trees:Bog", "US_Trees:BogForest",
                              "US_Trees", "Dicran_spp",
                              "Sphagnum:Bog", "Sphagnum:BogForest",
                              "Sphagnum", "Plot_PIUVcnt",
                              "DWT.Q:Bog", "DWT.Q:BogForest", "DWT.Q",
                              "DWT.L:Bog", "DWT.L:BogForest", "DWT.L",
                              "log_PAR1.3:Bog", "log_PAR1.3:BogForest",
                              "log_PAR1.3",
                              "Bog", "BogForest",
                              "Intercept"),
  labels = c("BogForest" = "Bog Forest",
             "log_PAR1.3" = "log(PAR 1.3)",
             "log_PAR1.3:BogForest" = "log(PAR 1.3) x Bog Forest",
             "log_PAR1.3:Bog" = "log(PAR 1.3) x Bog",
             "DWT.L" = "DWT.med",
             "DWT.L:BogForest" = "DWT.med x BogForest",
             "DWT.L:Bog" = "DWT.med x Bog",
             "DWT.Q" = "DWT.high",
             "DWT.Q:BogForest" = "DWT.high x Bog Forest",
             "DWT.Q:Bog" = "DWT.high x Bog",
             "Sphagnum:BogForest" = "Sphagnum x Bog Forest",
             "Sphagnum:Bog" = "Sphagnum x Bog",
             "US_Trees:BogForest" = "US_Trees x Bog Forest",
             "US_Trees:Bog" = "US_Trees x Bog",
             "shape_BogForest" = "shape_Bog Forest",
                              "Plot_PIUVcnt" = "mature_PIUV",
                              "Dicran_spp" = "Dicranaloma_spp")) +
  #coord_cartesian(xlim = c(-2.5, 4.0)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold")) +
  geom_vline(xintercept=0,linetype="dashed") +
  xlab("Stacked Model Estimate") + ylab("Parameter") +
  ggtitle("Control Site")
