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
          Plot_PIUVcnt.rs = as.vector(scale(Plot_PIUVcnt)),
          Trunk.rs = as.vector(scale(trunk)),
          Sphagnum.rs = as.vector(scale(Sphagnum)),
          Cushion_spp.rs = as.vector(scale(cush_spp)),
          Dicran_spp.rs = as.vector(scale(Dicran)),
          OtherMoss_spp.rs = as.vector(scale(om)),
          US_Emshr.rs = as.vector(scale(us_emshr)),
          US_Shrubs.rs = as.vector(scale(us_shrubs)),
          US_Trees.rs = as.vector(scale(us_trees)))
          

#verify & summarize tidied data by site
glimpse(sp.rmb)
summary(sp.rmb)
View(sp.rmb)

#write.csv(sp.rmb, "C:/Users/Marsh Hawk Zaret/Documents/Big_D/Data_Analysis/RQ1_Neo/Data/Ordination/rmb_explanvars_Mar2022.csv")

#################################################################################### DETERMINING OPTIMAL LIKELIHOOD & HIERARCHICAL STRUCTURE
###################################################################################

#===============================================================================
# Intercept-only models
#===============================================================================

#------------
# Fit models
#------------
# Negative binomial likelihood, rstanarm, no random effects
rmsa.glmer.nbn <- stan_glm.nb(PIUV_less5cmDBH_total ~ 1, data = sp.rmb, prior_intercept = normal(0, 1, autoscale = TRUE), na.action = na.omit, chains = 4,iter = 5000, warmup = 1000)
save(rmsa.glmer.nbn, file="F:/_R_objects/rmsa.glmer.nbn.rda")

# Negative binomial likelihood, rstanarm, random effects
rmsa.glmer.nb0 <- stan_glmer.nb(PIUV_less5cmDBH_total ~ 1 + (1 | plot), data = sp.rmb, prior_intercept = normal(0, 1, autoscale = TRUE), na.action = na.omit, chains = 4,iter = 5000, warmup = 1000)
save(rmsa.glmer.nb0, file="F:/_R_objects/rmsa.glmer.nb0.rda")

# Negative binomial likelihood, brms, random effects
rmsa.glmer.nb0.brm <- brm(PIUV_less5cmDBH_total ~ 1 + (1 | plot), family = negbinomial("log"), data = sp.rmb, prior = c(prior(normal(0, 1), class = "Intercept"), prior(normal(0, 1), class = "sd")), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb0.brm, file="F:/_R_objects/rmsa.glmer.nb0.brm.rda")

# Zero-inflated negative binomial, brms, random effects
rmsa.glmer.zinb0 <- brm(PIUV_less5cmDBH_total ~ 1 + (1 | plot), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = c(prior(normal(0, 1), class = "Intercept"), prior(normal(0, 1), class = "sd")),iter=7000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb0, file="F:/_R_objects/rmsa.glmer.zinb0.rda")


prior_summary(rmsa.glmer.zinb0)
summary(rmsa.glmer.zinb0)


load(file="F:/_R_objects/rmsa.glmer.nbn.rda")
load(file="F:/_R_objects/rmsa.glmer.nb0.rda")
load(file="F:/_R_objects/rmsa.glmer.nb0.brm.rda")
load(file="F:/_R_objects/rmsa.glmer.zinb0.rda")


#--------------
# Check models
#--------------

launch_shinystan(rmsa.glmer.zinb0)

dev.new()
plot(rmsa.glmer.zinb0)

dev.new()
mcmc_trace(rmsa.glmer.zinb0)

pp_check(rmsa.glmer.zinb0, nreps = 100) + 
  xlab("Seedling counts") + ggtitle("rmsa.glmer.zinb0")

pp_check(rmsa.glmer.zinb0, ndraws = 100) + 
  xlab("Seedling counts") + ggtitle("rmsa.glmer.zinb0")

y <- sp.rmb$PIUV_less5cmDBH_total # vector of outcome values
yrep <- posterior_predict(rmsa.glmer.nbn, draws = 500) # matrix of draws from the ppd

prop_zero <- function(x) mean(x == 0)

ppc_stat(y, yrep, stat = "prop_zero", binwidth = 0.005) +
  ggtitle("rmsa.glmer.nbn")
ppc_stat_grouped(y, yrep, sp.rmb$patch, stat = "prop_zero", binwidth = 0.005) +
  ggtitle("rmsa.glmer.zinb0")
ppc_stat(y, yrep, stat = "max", binwidth = 5) + ggtitle("rmsa.glmer.zinb0")
ppc_stat_grouped(y, yrep, sp.rmb$patch, stat = "max", binwidth = 5) + 
  ggtitle("rmsa.glmer.zinb0")
ppc_stat(y, yrep, stat = "sd", binwidth = 5) + ggtitle("rmsa.glmer.zinb0")
ppc_stat_grouped(y, yrep, sp.rmb$patch, stat = "sd", binwidth = 5) +
  ggtitle("rmsa.glmer.zinb0")


#--------------------------
# Validate & compare models
#--------------------------

loo_rmsa.nbn <- loo(rmsa.glmer.nbn, save_psis = TRUE, k_threshold = 0.7)
save(loo_rmsa.nbn, file="F:/_R_objects/loo_rmsa.nbn.rda")

loo_rmsa.nb0 <- loo(rmsa.glmer.nb0, save_psis = TRUE, k_threshold = 0.7)
#6 problematic obs
save(loo_rmsa.nb0, file="F:/_R_objects/loo_rmsa.nb0.rda")

loo_rmsa.nb0.brm <- loo(rmsa.glmer.nb0.brm, save_psis = TRUE, moment_match = TRUE)
save(loo_rmsa.nb0.brm, file="F:/_R_objects/loo_rmsa.nb0.brm.rda")

loo_rmsa.zinb0 <- loo(rmsa.glmer.zinb0, save_psis = TRUE, moment_match = TRUE)
save(loo_rmsa.zinb0, file="F:/_R_objects/loo_rmsa.zinb0.rda")

load(file="F:/_R_objects/loo_rmsa.nbn.rda")
load(file="F:/_R_objects/loo_rmsa.nb0.rda")
load(file="F:/_R_objects/loo_rmsa.nb0.brm.rda")
load(file="F:/_R_objects/loo_rmsa.zinb0.rda")

plot(loo_rmsa.nb0, label_points = TRUE)


psis  <- loo_rmsa.nb0$psis_object
lw     <- weights(psis)
y     <- sp.rmb$PIUV_less5cmDBH_total
yrep  <- posterior_predict(rmsa.glmer.nb0)

ppc_loo_intervals(y, yrep, psis, prob_outer= 0.95)
ppc_loo_ribbon(y, yrep, lw, psis, prob_outer = 0.95)


loo_compare(loo_rmsa.nbn, loo_rmsa.nb0, loo_rmsa.nb0.brm, loo_rmsa.zinb0)
#nbn has significantly less predictive density than other models, which are not statistically different


#===============================================================================
# Model 1 - Patch
#===============================================================================

#------------
# Fit models
#------------

?set_prior

# Negative binomial likelihood, rstanarm, random effects
rmsa.glmer.nb1 <- stan_glmer.nb(PIUV_less5cmDBH_total ~ patch + (1 | plot), data = sp.rmb, na.action = na.omit, chains = 4,iter = 5000, warmup = 1000, prior_intercept = normal(0, 1, autoscale = TRUE), prior = normal(0, 1, autoscale = TRUE))
save(rmsa.glmer.nb1, file="F:/_R_objects/rmsa.glmer.nb1.rda")


# Negative binomial likelihood, brms, random effects, predictors scaled before fit
rmsa.glmer.nb1.brm <- brm(PIUV_less5cmDBH_total ~ patch + (1 | plot), family = negbinomial("log"), data = sp.rmb, prior = set_prior("normal(0, 1)", class = "Intercept") + set_prior("normal(0, 0.5)", class = "b") + set_prior ("normal(0, 1)", class = "sd"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.nb1.brm, file="F:/_R_objects/rmsa.glmer.nb1.brm.rda")


# Zero-inflated negative binomial, brms, random effects, predictors scaled before fit
rmsa.glmer.zinb1 <- brm(bf(PIUV_less5cmDBH_total ~ patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(0, 1)", class = "Intercept") + set_prior("normal(0, 0.5)", class = "b") + set_prior ("normal(0, 1)", class = "sd"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb1, file="F:/_R_objects/rmsa.glmer.zinb1.rda")


load(file="F:/_R_objects/rmsa.glmer.nb1.rda")
load(file="F:/_R_objects/rmsa.glmer.nb1.brm.rda")
load(file="F:/_R_objects/rmsa.glmer.zinb1.rda")


prior_summary(rmsa.glmer.zinb1)
summary(rmsa.glmer.zinb1)
variables(rmsa.glmer.zinb1)

#--------------
# Check models
#--------------

launch_shinystan(rmsa.glmer.zinb1)

dev.new()
plot(rmsa.glmer.zinb1)

dev.new()
mcmc_trace(rmsa.glmer.zinb1)

pp_check(rmsa.glmer.zinb1, nreps = 100) + 
  xlab("Seedling counts") + ggtitle("rmsa.glmer.zinb1")

pp_check(rmsa.glmer.zinb1, ndraws = 100) + 
  xlab("Seedling counts") + ggtitle("rmsa.glmer.zinb1")

y <- sp.rmb$PIUV_less5cmDBH_total # vector of outcome values
yrep <- posterior_predict(rmsa.glmer.zinb1, draws = 500) # matrix of draws from the ppd

prop_zero <- function(x) mean(x == 0)

ppc_stat(y, yrep, stat = "prop_zero", binwidth = 0.005) +
  ggtitle("rmsa.glmer.zinb1")
ppc_stat_grouped(y, yrep, sp.rmb$patch, stat = "prop_zero", binwidth = 0.005) +
  ggtitle("rmsa.glmer.zinb1")
ppc_stat(y, yrep, stat = "max", binwidth = 5) + ggtitle("rmsa.glmer.zinb1")
ppc_stat_grouped(y, yrep, sp.rmb$patch, stat = "max", binwidth = 5) + 
  ggtitle("rmsa.glmer.zinb1")
ppc_stat(y, yrep, stat = "sd", binwidth = 5) + ggtitle("rmsa.glmer.zinb1")
ppc_stat_grouped(y, yrep, sp.rmb$patch, stat = "sd", binwidth = 5) +
  ggtitle("rmsa.glmer.zinb1")


get_variables(rmsa.glmer.zinb1)

# density plots of posterior interval estimates from MCMC draws
# + ggplot2::coord_cartesian(xlim = c(-1.5, 1.6))
dev.new()
mcmc_areas(rmsa.glmer.zinb1, prob = 0.95, pars=c("b_Intercept", "b_patchBogForest", "b_patchBog")) + ggtitle("rmsa.glmer.zinb1")

mcmc_areas(rmsa.glmer.zinb1, prob = 0.95, pars=c("b_zi_Intercept", "b_zi_patchBogForest", "b_zi_patchBog", "sd_plot__Intercept", "shape","Intercept", "Intercept_zi")) + ggtitle("rmsa.glmer.zinb1")

# interval estimates
dev.new()
mcmc_intervals(rmsa.glmer.zinb1, prob = 0.90, prob_outer = 0.95, pars=c("b_Intercept", "b_patchBogForest", "b_patchBog")) + ggtitle("rmsa.glmer.zinb1")

mcmc_intervals(rmsa.glmer.zinb1, prob = 0.90, prob_outer = 0.95, pars=c("b_zi_Intercept", "b_zi_patchBogForest", "b_zi_patchBog", "sd_plot__Intercept", "shape","Intercept", "Intercept_zi")) + ggtitle("rmsa.glmer.zinb1")

# pairwise plots
mcmc_pairs(as.matrix(rmsa.glmer.zinb1),pars = c("(Intercept)","sqrt_roach1","treatment","senior"))

# conditional effects
plot(conditional_effects(rmsa.glmer.zinb1), ask = FALSE)

#--------------------------
# Validate & compare models
#--------------------------

loo_rmsa.nb1 <- loo(rmsa.glmer.nb1, save_psis = TRUE, k_threshold = 0.7)
save(loo_rmsa.nb1, file="F:/_R_objects/loo_rmsa.nb1.rda")
#3 problematic obs

loo_rmsa.nb1.brm <- loo(rmsa.glmer.nb1.brm, save_psis = TRUE, moment_match = TRUE)
save(loo_rmsa.nb1.brm, file="F:/_R_objects/loo_rmsa.nb1.brm.rda")
#some slightly high Pareto k values

loo_rmsa.zinb1 <- loo(rmsa.glmer.zinb1, save_psis = TRUE, moment_match = TRUE)
save(loo_rmsa.zinb1, file="F:/_R_objects/loo_rmsa.zinb1.rda")
#some slightly high Pareto k values

load(file="F:/_R_objects/loo_rmsa.nb1.rda")
load(file="F:/_R_objects/loo_rmsa.nb1.brm.rda")
load(file="F:/_R_objects/loo_rmsa.zinb1.rda")

plot(loo_rmsa.nb0, label_points = TRUE)


psis  <- loo_rmsa.nb0$psis_object
lw     <- weights(psis)
y     <- sp.rmb$PIUV_less5cmDBH_total
yrep  <- posterior_predict(rmsa.glmer.nb0)

ppc_loo_intervals(y, yrep, psis, prob_outer= 0.95)
ppc_loo_ribbon(y, yrep, lw, psis, prob_outer = 0.95)


loo_compare(loo_rmsa.nb1, loo_rmsa.nb1.brm, loo_rmsa.zinb1) #none of these have significantly different epld values

loo_compare(loo_rmsa.zinb0, loo_rmsa.zinb1) #slight statistically significant difference between the intercept-only and full model

loo_compare(loo_rmsa.nbn, loo_rmsa.nb0, loo_rmsa.nb0.brm, loo_rmsa.zinb0, loo_rmsa.nb1, loo_rmsa.nb1.brm, loo_rmsa.zinb1) #the nbn model is statistically different (w/ lower ELPD) than the other models; the zinb0 is almost statistically different from the zinb1


#===============================================================================
# Model 2 - PAR & Patch w/ interaction; explore prior specification
#===============================================================================

#------------
# Fit models
#------------

#zinb2a
rmsa.glmer.zinb2a <- brm(bf(PIUV_less5cmDBH_total ~ log_PAR1.3.rs + patch + log_PAR1.3.rs:patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(0, 1)", class = "Intercept") + set_prior("normal(0, 0.5)", class = "b") + set_prior ("normal(0, 1)", class = "sd"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb2a, file="F:/_R_objects/rmsa.glmer.zinb2a.rda")

#zinb2b
rmsa.glmer.zinb2b <- brm(bf(PIUV_less5cmDBH_total ~ log_PAR1.3.rs + patch + log_PAR1.3.rs:patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(0, 1)", class = "Intercept") + set_prior("normal(0, 0.5)", class = "b") + set_prior ("normal(0, 1)", class = "sd") + set_prior ("gamma(0.5, 0.5)", class = "shape"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb2b, file="F:/_R_objects/rmsa.glmer.zinb2b.rda")

#zinb2c
rmsa.glmer.zinb2c <- brm(bf(PIUV_less5cmDBH_total ~ log_PAR1.3.rs + patch + log_PAR1.3.rs:patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(0, 1)", class = "Intercept") + set_prior("normal(0, 0.5)", class = "b") + set_prior ("normal(0, 1)", class = "sd") + set_prior ("logistic(0, 0.5)", class = "Intercept", dpar="zi"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb2c, file="F:/_R_objects/rmsa.glmer.zinb2c.rda")

#zinb2d
rmsa.glmer.zinb2d <- brm(bf(PIUV_less5cmDBH_total ~ log_PAR1.3.rs + patch + log_PAR1.3.rs:patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(0, 1)", class = "Intercept") + set_prior("normal(0, 0.5)", class = "b") + set_prior ("normal(0, 1)", class = "sd") + set_prior ("logistic(0, 0.5)", class = "Intercept", dpar="zi") + set_prior ("gamma(0.5, 0.5)", class = "shape"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb2d, file="F:/_R_objects/rmsa.glmer.zinb2d.rda")

#zinb2e
rmsa.glmer.zinb2e <- brm(bf(PIUV_less5cmDBH_total ~ log_PAR1.3.rs + patch + log_PAR1.3.rs:patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(0, 0.5)", class = "Intercept") + set_prior("normal(0, 0.5)", class = "b") + set_prior ("normal(0, 0.5)", class = "sd") + set_prior ("logistic(0, 0.5)", class = "Intercept", dpar="zi") + set_prior ("gamma(0.5, 0.5)", class = "shape"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb2e, file="F:/_R_objects/rmsa.glmer.zinb2e.rda")

#zinb2f
rmsa.glmer.zinb2f <- brm(bf(PIUV_less5cmDBH_total ~ log_PAR1.3.rs + patch + log_PAR1.3.rs:patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(0, 1)", class = "Intercept") + set_prior("normal(0, 0.5)", class = "b") + set_prior ("normal(0, 1)", class = "sd") + set_prior("logistic(0,0.1)", class = "Intercept", dpar="zi"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb2f, file="F:/_R_objects/rmsa.glmer.zinb2f.rda")

#zinb2g
rmsa.glmer.zinb2g <- brm(bf(PIUV_less5cmDBH_total ~ log_PAR1.3.rs + patch + log_PAR1.3.rs:patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(0, 0.5)", class = "Intercept") + set_prior("normal(0, 0.5)", class = "b") + set_prior ("normal(0, 0.5)", class = "sd") + set_prior ("logistic(0, 0.1)", class = "Intercept", dpar="zi") + set_prior ("gamma(0.05, 0.05)", class = "shape"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb2g, file="F:/_R_objects/rmsa.glmer.zinb2g.rda")

#zinb2h
rmsa.glmer.zinb2h <- brm(bf(PIUV_less5cmDBH_total ~ log_PAR1.3.rs + patch + log_PAR1.3.rs:patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(0, 1)", class = "Intercept") + set_prior("normal(0, 0.5)", class = "b") + set_prior ("normal(0, 0.5)", class = "sd") + set_prior ("logistic(0, 0.5)", class = "Intercept", dpar="zi") + set_prior ("gamma(0.5, 0.5)", class = "shape"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb2h, file="F:/_R_objects/rmsa.glmer.zinb2h.rda")


load(file="F:/_R_objects/rmsa.glmer.zinb2a.rda")
load(file="F:/_R_objects/rmsa.glmer.zinb2b.rda")
load(file="F:/_R_objects/rmsa.glmer.zinb2c.rda")
load(file="F:/_R_objects/rmsa.glmer.zinb2d.rda")
load(file="F:/_R_objects/rmsa.glmer.zinb2e.rda")
load(file="F:/_R_objects/rmsa.glmer.zinb2f.rda")
load(file="F:/_R_objects/rmsa.glmer.zinb2g.rda")
load(file="F:/_R_objects/rmsa.glmer.zinb2h.rda")


prior_summary(rmsa.glmer.zinb2h)
summary(rmsa.glmer.zinb2e)
variables(rmsa.glmer.zinb2h)

?set_prior

#--------------
# Check models
#--------------

mod <- rmsa.glmer.zinb2h

launch_shinystan(mod)

bayesplot_theme_set(ggplot2::theme_bw())
color_scheme_set("brightblue")

# posterior vs. prior
variables(rmsa.glmer.zinb2i)
param <- "b_Intercept"
prior_summary(rmsa.glmer.zinb2i)    # inspect prior for param
pdf <- function(x) dnorm(x, 0, 1) # density function for prior on param
pdf <- function(x) dgamma(x, 1, 0.1)
pdf <- function(x) dlogis(x, 0, 1)

hist(as.matrix(rmsa.glmer.zinb2i, variable = param), 25, prob = TRUE,
     col = "gray", border = "white", 
     las = 1, cex.axis = 1.2, cex.lab = 1.5, cex.main = 1.5,
     xlim = c(0, 5),
     xlab = param, ylab = "Probability", main = "prior vs. posterior, mod 2i")
curve(pdf(x), col = "blue", lwd = 2, add = TRUE)

# density & trace plots
dev.new()
plot(mod)

dev.new()
mcmc_trace(mod)

?pp_check
pp_check(mod, ndraws = 250) + 
  xlab("Seedling counts") + ggtitle("mod")

y <- sp.rmb$PIUV_less5cmDBH_total # vector of outcome values
yrep <- posterior_predict(mod, draws = 500) # matrix of draws from the ppd

prop_zero <- function(x) mean(x == 0)

ppc_stat(y, yrep, stat = "prop_zero", binwidth = 0.005) +
  ggtitle("mod")
ppc_stat_grouped(y, yrep, sp.rmb$patch, stat = "prop_zero", binwidth = 0.005) +
  ggtitle("mod")
ppc_stat(y, yrep, stat = "max", binwidth = 5) + ggtitle("mod")
ppc_stat_grouped(y, yrep, sp.rmb$patch, stat = "max", binwidth = 5) + 
  ggtitle("mod")
ppc_stat(y, yrep, stat = "sd", binwidth = 5) + ggtitle("mod")
ppc_stat_grouped(y, yrep, sp.rmb$patch, stat = "sd", binwidth = 5) +
  ggtitle("mod")


# interval estimates
get_variables(mod)
dev.new()
mcmc_intervals(mod, prob = 0.90, prob_outer = 0.95,pars=c("b_Intercept", "b_zi_Intercept", "b_log_PAR1.3.rs", "b_patchBogForest", "b_patchBog", "b_log_PAR1.3.rs:patchBogForest", "b_log_PAR1.3.rs:patchBog", "b_zi_patchBogForest", "b_zi_patchBog")) + ggtitle("mod")

mcmc_intervals(mod, prob = 0.90, prob_outer = 0.95, pars=c("sd_plot__Intercept", "shape","Intercept", "Intercept_zi")) + ggtitle("mod")

# pairwise plots
dev.new()
mcmc_pairs(as.matrix(mod),pars = c("b_Intercept","b_log_PAR1.3.rs", "b_patchBogForest", "b_patchBog")) + ggtitle("mod")

# conditional effects
plot(conditional_effects(mod), ask = FALSE) + ggtitle("mod")





######################## OLD / TO EDIT ###################################################
#################################################################################### FIT & EVALUATE MODELS 2 - 12 USING OPTIMAL LIKELIHOOD SPECIFICATION
###################################################################################

#===============================================================================
# Fit models
#===============================================================================

#zinb2
rmsa.glmer.zinb2 <- brm(bf(PIUV_less5cmDBH_total ~ minwtd.rs + log_par1.3.rs + subNMDS_a1.rs + subNMDS_a2.rs + subNMDS_a3.rs + Plot_PIUVcnt.rs + patch + minwtd.rs:patch + log_par1.3.rs:patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(0, 1)", class = "Intercept") + set_prior("normal(0, 0.5)", class = "b") + set_prior ("normal(0, 1)", class = "sd"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb2, file="F:/_R_objects/rmsa.glmer.zinb2.rda")


#zinb3
rmsa.glmer.zinb3 <- brm(bf(PIUV_less5cmDBH_total ~ minwtd.rs + log_par1.3.rs + Sphagnum.rs + ns_moss.rs + Plot_PIUVcnt.rs + patch + minwtd.rs:patch + log_par1.3.rs:patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(0, 1)", class = "Intercept") + set_prior("normal(0, 0.5)", class = "b") + set_prior ("normal(0, 1)", class = "sd"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb3, file="F:/_R_objects/rmsa.glmer.zinb3.rda")


#zinb4
rmsa.glmer.zinb4 <- brm(bf(PIUV_less5cmDBH_total ~ minwtd.rs + log_par1.3.rs + Sphagnum.rs + Plot_PIUVcnt.rs + patch + minwtd.rs:patch + log_par1.3.rs:patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(0, 1)", class = "Intercept") + set_prior("normal(0, 0.5)", class = "b") + set_prior ("normal(0, 1)", class = "sd"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb4, file="F:/_R_objects/rmsa.glmer.zinb4.rda")


#zinb5
rmsa.glmer.zinb5 <- brm(bf(PIUV_less5cmDBH_total ~ minwtd.rs + log_par1.3.rs + patch + minwtd.rs:patch + log_par1.3.rs:patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(0, 1)", class = "Intercept") + set_prior("normal(0, 0.5)", class = "b") + set_prior ("normal(0, 1)", class = "sd"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb5, file="F:/_R_objects/rmsa.glmer.zinb5.rda")


#zinb6
rmsa.glmer.zinb6 <- brm(bf(PIUV_less5cmDBH_total ~ minwtd.rs + patch + minwtd.rs:patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(0, 1)", class = "Intercept") + set_prior("normal(0, 0.5)", class = "b") + set_prior ("normal(0, 1)", class = "sd"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb6, file="F:/_R_objects/rmsa.glmer.zinb6.rda")


#zinb7
rmsa.glmer.zinb7 <- brm(bf(PIUV_less5cmDBH_total ~ log_par1.3.rs + patch + log_par1.3.rs:patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(0, 1)", class = "Intercept") + set_prior("normal(0, 0.5)", class = "b") + set_prior ("normal(0, 1)", class = "sd"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb7, file="F:/_R_objects/rmsa.glmer.zinb7.rda")


#zinb8
rmsa.glmer.zinb8 <- brm(bf(PIUV_less5cmDBH_total ~ NMDS_axis1.rs + NMDS_axis2.rs + NMDS_axis3.rs + patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(0, 1)", class = "Intercept") + set_prior("normal(0, 0.5)", class = "b") + set_prior ("normal(0, 1)", class = "sd"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb8, file="F:/_R_objects/rmsa.glmer.zinb8.rda")


#zinb9
rmsa.glmer.zinb9 <- brm(bf(PIUV_less5cmDBH_total ~ subNMDS_a1.rs + subNMDS_a2.rs + subNMDS_a3.rs + patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(0, 1)", class = "Intercept") + set_prior("normal(0, 0.5)", class = "b") + set_prior ("normal(0, 1)", class = "sd"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb9, file="F:/_R_objects/rmsa.glmer.zinb9.rda")


#zinb10
rmsa.glmer.zinb10 <- brm(bf(PIUV_less5cmDBH_total ~ Sphagnum.rs + patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(0, 1)", class = "Intercept") + set_prior("normal(0, 0.5)", class = "b") + set_prior ("normal(0, 1)", class = "sd"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb10, file="F:/_R_objects/rmsa.glmer.zinb10.rda")


#zinb11
rmsa.glmer.zinb11 <- brm(bf(PIUV_less5cmDBH_total ~ Plot_PIUVcnt.rs + patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(0, 1)", class = "Intercept") + set_prior("normal(0, 0.5)", class = "b") + set_prior ("normal(0, 1)", class = "sd"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb11, file="F:/_R_objects/rmsa.glmer.zinb11.rda")


#zinb12
rmsa.glmer.zinb12 <- brm(bf(PIUV_less5cmDBH_total ~ patch + (1 | plot), zi ~ patch), family = zero_inflated_negbinomial("log"), data = sp.rmb, prior = set_prior("normal(0, 1)", class = "Intercept") + set_prior("normal(0, 0.5)", class = "b") + set_prior ("normal(0, 1)", class = "sd"), iter=5000, save_pars = save_pars(all=TRUE))
save(rmsa.glmer.zinb12, file="F:/_R_objects/rmsa.glmer.zinb12.rda")


load(file="F:/_R_objects/rmsa.glmer.zinb1.rda")
load(file="F:/_R_objects/rmsa.glmer.zinb2.rda")
load(file="F:/_R_objects/rmsa.glmer.zinb3.rda")
load(file="F:/_R_objects/rmsa.glmer.zinb4.rda")
load(file="F:/_R_objects/rmsa.glmer.zinb5.rda")
load(file="F:/_R_objects/rmsa.glmer.zinb6.rda")
load(file="F:/_R_objects/rmsa.glmer.zinb7.rda")
load(file="F:/_R_objects/rmsa.glmer.zinb8.rda")
load(file="F:/_R_objects/rmsa.glmer.zinb9.rda")
load(file="F:/_R_objects/rmsa.glmer.zinb10.rda")
load(file="F:/_R_objects/rmsa.glmer.zinb11.rda")
load(file="F:/_R_objects/rmsa.glmer.zinb12.rda")


#===============================================================================
# Check models
#===============================================================================

prior_summary(rmsa.glmer.zinb2)
summary(rmsa.glmer.zinb2)

launch_shinystan(rmsa.glmer.zinb2)

dev.new()
plot(rmsa.glmer.zinb2)

dev.new()
mcmc_trace(rmsa.glmer.zinb2)

pp_check(rmsa.glmer.zinb2, nreps = 100) + 
  xlab("Seedling counts") + ggtitle("rmsa.glmer.zinb2")

pp_check(rmsa.glmer.zinb2, ndraws = 100) + 
  xlab("Seedling counts") + ggtitle("rmsa.glmer.zinb2")

y <- sp.rmb$PIUV_less5cmDBH_total # vector of outcome values
yrep <- posterior_predict(rmsa.glmer.zinb2, draws = 500) # matrix of draws from the ppd

prop_zero <- function(x) mean(x == 0)

ppc_stat(y, yrep, stat = "prop_zero", binwidth = 0.005) +
  ggtitle("rmsa.glmer.zinb2")
ppc_stat_grouped(y, yrep, sp.rmb$patch, stat = "prop_zero", binwidth = 0.005) +
  ggtitle("rmsa.glmer.zinb2")
ppc_stat(y, yrep, stat = "max", binwidth = 5) + ggtitle("rmsa.glmer.zinb2")
ppc_stat_grouped(y, yrep, sp.rmb$patch, stat = "max", binwidth = 5) + 
  ggtitle("rmsa.glmer.zinb2")
ppc_stat(y, yrep, stat = "sd", binwidth = 5) + ggtitle("rmsa.glmer.zinb2")
ppc_stat_grouped(y, yrep, sp.rmb$patch, stat = "sd", binwidth = 5) +
  ggtitle("rmsa.glmer.zinb2")


get_variables(rmsa.glmer.zinb2)
dev.new()

# density plots of posterior interval estimates from MCMC draws
mcmc_areas(rmsa.glmer.zinb2, prob = 0.95, pars=c("b_Intercept", "b_minwtd.rs", "b_log_par1.3.rs", "b_Sphagnum.rs", "b_ns_moss.rs", "b_Plot_PIUVcnt.rs", "b_patchBogForest", "b_patchBog", "b_minwtd.rs:patchBogForest", "b_minwtd.rs:patchBog", "b_log_par1.3.rs:patchBogForest", "b_log_par1.3.rs:patchBog")) + ggtitle("rmsa.glmer.zinb2")

mcmc_areas(rmsa.glmer.zinb2, prob = 0.95, pars=c("b_zi_Intercept", "b_zi_patchBogForest", "b_zi_patchBog", "sd_plot__Intercept", "shape","Intercept", "Intercept_zi")) + ggtitle("rmsa.glmer.zinb2")


#===============================================================================
# Validate & compare models
#===============================================================================

loo_rmsa.zinb2 <- loo(rmsa.glmer.zinb2, save_psis = TRUE, moment_match = TRUE)
save(loo_rmsa.zinb2, file="F:/_R_objects/loo_rmsa.zinb2.rda")

loo_rmsa.zinb3 <- loo(rmsa.glmer.zinb3, save_psis = TRUE, moment_match = TRUE)
save(loo_rmsa.zinb3, file="F:/_R_objects/loo_rmsa.zinb3.rda")

loo_rmsa.zinb4 <- loo(rmsa.glmer.zinb4, save_psis = TRUE, moment_match = TRUE)
save(loo_rmsa.zinb4, file="F:/_R_objects/loo_rmsa.zinb4.rda")

loo_rmsa.zinb5 <- loo(rmsa.glmer.zinb5, save_psis = TRUE, moment_match = TRUE)
save(loo_rmsa.zinb5, file="F:/_R_objects/loo_rmsa.zinb5.rda")

loo_rmsa.zinb6 <- loo(rmsa.glmer.zinb6, save_psis = TRUE, moment_match = TRUE)
save(loo_rmsa.zinb6, file="F:/_R_objects/loo_rmsa.zinb6.rda")

loo_rmsa.zinb7 <- loo(rmsa.glmer.zinb7, save_psis = TRUE, moment_match = TRUE)
save(loo_rmsa.zinb7, file="F:/_R_objects/loo_rmsa.zinb7.rda")

loo_rmsa.zinb8 <- loo(rmsa.glmer.zinb8, save_psis = TRUE, moment_match = TRUE, reloo = TRUE)
save(loo_rmsa.zinb8, file="F:/_R_objects/loo_rmsa.zinb8.rda")

loo_rmsa.zinb9 <- loo(rmsa.glmer.zinb9, save_psis = TRUE, moment_match = TRUE)
save(loo_rmsa.zinb9, file="F:/_R_objects/loo_rmsa.zinb9.rda")

loo_rmsa.zinb10 <- loo(rmsa.glmer.zinb10, save_psis = TRUE, moment_match = TRUE)
save(loo_rmsa.zinb10, file="F:/_R_objects/loo_rmsa.zinb10.rda")

loo_rmsa.zinb11 <- loo(rmsa.glmer.zinb11, save_psis = TRUE, moment_match = TRUE)
save(loo_rmsa.zinb11, file="F:/_R_objects/loo_rmsa.zinb11.rda")

loo_rmsa.zinb12 <- loo(rmsa.glmer.zinb12, save_psis = TRUE, moment_match = TRUE)
save(loo_rmsa.zinb12, file="F:/_R_objects/loo_rmsa.zinb12.rda")


load(file="F:/_R_objects/loo_rmsa.zinb0.rda")
load(file="F:/_R_objects/loo_rmsa.zinb1.rda")
load(file="F:/_R_objects/loo_rmsa.zinb2.rda")
load(file="F:/_R_objects/loo_rmsa.zinb3.rda")
load(file="F:/_R_objects/loo_rmsa.zinb4.rda")
load(file="F:/_R_objects/loo_rmsa.zinb5.rda")
load(file="F:/_R_objects/loo_rmsa.zinb6.rda")
load(file="F:/_R_objects/loo_rmsa.zinb7.rda")
load(file="F:/_R_objects/loo_rmsa.zinb8.rda")
load(file="F:/_R_objects/loo_rmsa.zinb9.rda")
load(file="F:/_R_objects/loo_rmsa.zinb10.rda")
load(file="F:/_R_objects/loo_rmsa.zinb11.rda")
load(file="F:/_R_objects/loo_rmsa.zinb12.rda")


plot(loo_rmsa.nb0, label_points = TRUE)


psis  <- loo_rmsa.zinb2$psis_object
lw     <- weights(psis)
y     <- sp.rmb$PIUV_less5cmDBH_total
yrep  <- posterior_predict(rmsa.glmer.zinb2)

ppc_loo_intervals(y, yrep, psis, prob_outer= 0.95)
ppc_loo_ribbon(y, yrep, lw, psis, prob_outer = 0.95)


loo_compare(loo_rmsa.zinb0, loo_rmsa.zinb1, loo_rmsa.zinb2, loo_rmsa.zinb3, loo_rmsa.zinb4, loo_rmsa.zinb5, loo_rmsa.zinb6, loo_rmsa.zinb7, loo_rmsa.zinb8, loo_rmsa.zinb9, loo_rmsa.zinb10, loo_rmsa.zinb11, loo_rmsa.zinb12) #zinb5 has the greatest predictive density, but only statistically different from zinb0, perhaps statistically different from zinb9 (which had the least predictive density)


#################################################################################### CONSTRUCT & PLOT A FINAL STACKED MODEL
###################################################################################

#===============================================================================
# Determine model averaging weights using Bayesian stacking
#===============================================================================

# When we have more than two models, it can be easier to understand the uncertainties in the comparison by looking at the model averaging weights based on Bayesian stacking (Yao et al., 2018) which optimizes the model weights so that the combined predictive distribution maximizes the estimate leave-one-out cross-validation performance. https://avehtari.github.io/modelselection/rats_kcv.html#4_Leave-one-out_cross-validation

# full models that vary by substrate predictors
lpd_point <- cbind(
  loo_rmsa.zinb1$pointwise[,"elpd_loo"],
  loo_rmsa.zinb2$pointwise[,"elpd_loo"],
  loo_rmsa.zinb3$pointwise[,"elpd_loo"],
  loo_rmsa.zinb4$pointwise[,"elpd_loo"])
stacking_weights(lpd_point)

# models with single predictors (except zinb5), interaction terms if applicable, & patch (zinb12 just patch)
lpd_point <- cbind(
  loo_rmsa.zinb5$pointwise[,"elpd_loo"],
  loo_rmsa.zinb6$pointwise[,"elpd_loo"],
  loo_rmsa.zinb7$pointwise[,"elpd_loo"],
  loo_rmsa.zinb8$pointwise[,"elpd_loo"],
  loo_rmsa.zinb9$pointwise[,"elpd_loo"],
  loo_rmsa.zinb10$pointwise[,"elpd_loo"],
  loo_rmsa.zinb11$pointwise[,"elpd_loo"],
  loo_rmsa.zinb12$pointwise[,"elpd_loo"])
stacking_weights(lpd_point)

# all models
lpd_point <- cbind(
  loo_rmsa.zinb1$pointwise[,"elpd_loo"],
  loo_rmsa.zinb2$pointwise[,"elpd_loo"],
  loo_rmsa.zinb3$pointwise[,"elpd_loo"],
  loo_rmsa.zinb4$pointwise[,"elpd_loo"],
  loo_rmsa.zinb5$pointwise[,"elpd_loo"],
  loo_rmsa.zinb6$pointwise[,"elpd_loo"],
  loo_rmsa.zinb7$pointwise[,"elpd_loo"],
  loo_rmsa.zinb8$pointwise[,"elpd_loo"],
  loo_rmsa.zinb9$pointwise[,"elpd_loo"],
  loo_rmsa.zinb10$pointwise[,"elpd_loo"],
  loo_rmsa.zinb11$pointwise[,"elpd_loo"],
  loo_rmsa.zinb12$pointwise[,"elpd_loo"])
stacking_weights(lpd_point)


#===============================================================================
# ADD POSTERIOR DRAWS FROM MODELS (BASED ON STACKING WEIGHTS) TO SINGLE TIBBLE
#===============================================================================

variables <- as.list(get_variables(rmsa.glmer.zinb1))
variables #note:  need backticks (``) around names that use colons

# extract draws from posterior distributions of each model; set number of draws to stacking weight * 10000
post.rmsa.glmer.zinb1 <- rmsa.glmer.zinb1 %>%
  spread_draws(b_Intercept, b_zi_Intercept, b_minwtd.rs, b_log_par1.3.rs,
               b_NMDS_axis1.rs, b_NMDS_axis2.rs, b_NMDS_axis3.rs,
               b_Plot_PIUVcnt.rs,
               b_patchBogForest, b_patchBog,
               `b_minwtd.rs:patchBogForest`, `b_minwtd.rs:patchBog`, 
               `b_log_par1.3.rs:patchBogForest`, `b_log_par1.3.rs:patchBog`,
               b_zi_patchBogForest, b_zi_patchBog, sd_plot__Intercept, shape,
               Intercept, Intercept_zi, ndraws=3070)

post.rmsa.glmer.zinb5 <- rmsa.glmer.zinb5 %>%
  spread_draws(b_Intercept, b_zi_Intercept, b_minwtd.rs, b_log_par1.3.rs,
               b_patchBogForest, b_patchBog,
               `b_minwtd.rs:patchBogForest`, `b_minwtd.rs:patchBog`, 
               `b_log_par1.3.rs:patchBogForest`, `b_log_par1.3.rs:patchBog`,
               b_zi_patchBogForest, b_zi_patchBog, sd_plot__Intercept, shape,
               Intercept, Intercept_zi, ndraws=3130)

post.rmsa.glmer.zinb7 <- rmsa.glmer.zinb7 %>%
  spread_draws(b_Intercept, b_zi_Intercept, b_log_par1.3.rs,
               b_patchBogForest, b_patchBog,
               `b_log_par1.3.rs:patchBogForest`, `b_log_par1.3.rs:patchBog`,
               b_zi_patchBogForest, b_zi_patchBog, sd_plot__Intercept, shape,
               Intercept, Intercept_zi, ndraws=1490)

post.rmsa.glmer.zinb12 <- rmsa.glmer.zinb12 %>%
  spread_draws(b_Intercept, b_zi_Intercept,
               b_patchBogForest, b_patchBog,
               b_zi_patchBogForest, b_zi_patchBog, sd_plot__Intercept, shape,
               Intercept, Intercept_zi, ndraws=2310)


glimpse(post.rmsa.glmer.zinb1)


# bind posterior draws from all models into a single tibble & tidy

post.rmsa.glmer.zinb.stack <- bind_rows(post.rmsa.glmer.zinb1, post.rmsa.glmer.zinb5, post.rmsa.glmer.zinb7, post.rmsa.glmer.zinb12)

post.rmsa.glmer.zinb.stack <- post.rmsa.glmer.zinb.stack %>%
  mutate(draw_orig = .draw) %>% # save original draw numbers just in case
  mutate(.draw = 1:10000) # each .draw needs to be unique

glimpse(post.rmsa.glmer.zinb.stack)

save(post.rmsa.glmer.zinb.stack, file="F:/_R_objects/post.rmsa.glmer.zinb.stack.rda")
load(file="F:/_R_objects/post.rmsa.glmer.zinb.stack.rda")


#===============================================================================
# SUMMARIZE & VISUALIZE FINAL STACKED MODEL
#===============================================================================

variables2 <- get_variables(post.rmsa.glmer.zinb.stack)
variables2

#### point summaries and intervals for posterior draws of stacked models
post.rmsa.glmer.zinb.stack.sum <- post.rmsa.glmer.zinb.stack %>%
  gather_draws(b_Intercept, b_zi_Intercept, b_minwtd.rs, b_log_par1.3.rs,
               b_NMDS_axis1.rs, b_NMDS_axis2.rs, b_NMDS_axis3.rs,
               b_Plot_PIUVcnt.rs,
               b_patchBogForest, b_patchBog,
               `b_minwtd.rs:patchBogForest`, `b_minwtd.rs:patchBog`, 
               `b_log_par1.3.rs:patchBogForest`,`b_log_par1.3.rs:patchBog`,
               b_zi_patchBogForest, b_zi_patchBog, sd_plot__Intercept,
               shape, Intercept, Intercept_zi) %>%
  median_qi(na.rm = TRUE, .width = c(.95, .9, .8, .5))

View(post.rmsa.glmer.zinb.stack.sum) 
write.table(post.rmsa.glmer.zinb.stack.sum, "clipboard", sep="\t", row.names=TRUE)

?median_qi
#----------
# plots
#----------

#point intervals of parameter estimates from posterior distribution

#non-zi parameters
post.rmsa.glmer.zinb.stack %>%
  gather_draws(b_Intercept, b_minwtd.rs, b_log_par1.3.rs,
               b_NMDS_axis1.rs, b_NMDS_axis2.rs, b_NMDS_axis3.rs,
               b_Plot_PIUVcnt.rs,
               b_patchBogForest, b_patchBog,
               `b_minwtd.rs:patchBogForest`, `b_minwtd.rs:patchBog`, 
               `b_log_par1.3.rs:patchBogForest`,`b_log_par1.3.rs:patchBog`, 
               sd_plot__Intercept) %>%
  median_qi(.width = c(.95, .9), na.rm = TRUE) %>%
  mutate(.variable = str_remove(.variable, "b_")) %>% 
  mutate(.variable = str_remove(.variable, ".rs")) %>% 
  ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper)) +
  geom_pointinterval(color="purple") +
  #geom_slab() +
  scale_y_discrete(limits = c("sd_plot__Intercept", "log_par1.3:patchBog", "log_par1.3:patchBogForest",  "minwtd:patchBog",  "minwtd:patchBogForest", "patchBog", "patchBogForest", "Plot_PIUVcnt", "NMDS_axis3", "NMDS_axis2", "NMDS_axis1", "log_par1.3", "minwtd", "Intercept")) +
  coord_cartesian(xlim = c(-2.5, 4.0)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold")) +
  geom_vline(xintercept=0,linetype="dashed") +
  xlab("Stacked Model Estimate") + ylab("Parameter") +
  ggtitle("Control Site")

#zi (and shape) parameter
post.rmsa.glmer.zinb.stack %>%
  gather_draws(b_zi_Intercept, b_zi_patchBogForest, b_zi_patchBog,
               shape, Intercept_zi, Intercept) %>%
  median_qi(.width = c(.95, .9), na.rm = TRUE) %>%
  mutate(.variable = str_remove(.variable, "b_")) %>% 
  ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper)) +
  geom_pointinterval(color="purple") +
  scale_y_discrete(limits = c("shape", "Intercept", "Intercept_zi", "zi_patchBog", "zi_patchBogForest", "zi_Intercept")) +
  coord_cartesian(xlim = c(-15, 15)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold")) +
  geom_vline(xintercept=0,linetype="dashed") +
  xlab("Stacked Model Estimate") + ylab("Parameter") +
  ggtitle("Control Site")
