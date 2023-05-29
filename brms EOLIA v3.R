###############################
#EOLIA re-analysis using brms#
##############################
# the EOLIA trial was a multcentre randomized trial of artificial lung support (ECMO) vs conventional treatment in patients with severe pnumonia
# the observed effect size was a 9% decrease in absolute mortality in the ECMO group
# the trial had a small sample size and was under powered
# the p-value on the primary mortality outcome was 0.09


# create data frame based on published results
# 0 = alive, 1 = dead
# 57/125 dead in control
# 44/125 dead in ECMO group

# create vectors
treat <- rep(c("Control","ECMO"),times=c(125, 125))
dead <- rep(c(0,1,0,1),times=c(68,57,81,44))

# convert to factors
treat <- as.factor(treat)
dead <- as.factor(dead)

# check levels
levels(treat)
levels(dead)


# create data frame
EOLIA <- data.frame(treat, dead)

# check values are correct
table(EOLIA)


# run frequentest logistic regression model
EOLIA_logistic <- glm(dead ~ treat, 
                      family = 'binomial', data = EOLIA)
summary(EOLIA_logistic)


# obtain point estimate (95% confidence interval) for OR the treat = ECMO relative to control
exp(coef(EOLIA_logistic)["treatECMO"])
exp(confint(EOLIA_logistic,"treatECMO"))

# do the Bayesian analysis in brms

#load packages
library(tidyverse)
library(brms)
library(ProbBayes)
library(marginaleffects)
library(collapse)
library(ggplot2)
library(ggdist)
library(qwraps2)
library(bayestestR)

# brms model default flat priors
EOLIA_flat <- brm(dead ~ treat, 
                  bernoulli(link = "logit"), 
                  data = EOLIA,
                  chains = 2, # nb of chains
                  iter = 5000, # nb of iterations, including burnin
                  warmup = 1000, # burnin
                  thin = 1)
summary(EOLIA_flat)

# obtain odds ratios an associated 95% credible interval
# first extract the fixed effects from the model
fixedEffects_flat <- fixef(EOLIA_flat)
fixedEffects_flat

# exponentiate the fixed effect coefficients to obtain odds ratios
OR_EOLIA_flat <- exp(fixef(EOLIA_flat))
OR_EOLIA_flat

# comparison of frequentest and brms model with flat priors:
# frequentest - OR (95% confidence interval): 0.65 (0.39, 1.08) 
# flat priors - OR (95% credible interval) 0.65 (0.38, 1.08)

# plot the OR and 95% credible interval
mcmc_plot(EOLIA_flat, 
          type = "areas",
          prob = 0.95,
          transformations = "exp") +
  geom_vline(xintercept = 1, color = "grey")

#####################################
## Choosing more informative priors ##
######################################

# extract priors for EOLIA_brms_flat
prior_summary(EOLIA_flat)

# prior for the control (expected event rate of 0.5 with 90% density less than 0.6)
# use beta.select from the learn package
beta.select(list(x = 0.5, p = 0.5),
            list(x = 0.6, p = 0.9))
# gives a beta distribution with parameters 20.37 and 20.37

# convert to normal distribution
set.seed(2)
p_sim_neutral <- rbeta(1000, 20.37, 20.37)
theta_sim_neutral <- logit(p_sim_neutral)
c(mean(theta_sim_neutral), sd(theta_sim_neutral))
# gives a normal distribution with mean  0.01159056 and sd 0.33132171 - use for neutral prior for control and treatment groups

#  set the strong prior for the ECMO group - expected event rate of 0.4 with 90% of density less than 0.5
beta.select(list(x = 0.4, p = 0.5),
            list(x = 0.5, p = 0.9))
# gives gives a beta distribution with shape parameters 16.59 and 24.72

# convert to normal distribution
set.seed(3)
p_sim_strong <- rbeta(1000, 16.59, 24.72)
theta_sim_stong <- logit(p_sim_strong)
c(mean(theta_sim_stong), sd(theta_sim_stong))
# gives a normal distribution with mean -0.4168554 and sd of  0.3177808 - use as strong prior for the ECMO group


# set the very strong prior for the ECMO group - expected event rate of 0.3 with 90% density less than 0.4

beta.select(list(x = 0.3, p = 0.5),
            list(x = 0.4, p = 0.9))
# gives gives a beta distribution with shape parameters 11.68 and 26.81

# convert to normal distribution
set.seed(4)
p_sim_vstrong <- rbeta(1000, 11.68, 26.81)
theta_sim_vstong <- logit(p_sim_vstrong)
c(mean(theta_sim_vstong), sd(theta_sim_vstong))
# gives a normal distribution with mean -0.8612911 and sd of  0.3417656 - use as strong prior for the ECMO group




# define neutral prior using identical neutral priors for the control and ECMO groups
prior_neutral <- c(
  set_prior("normal(0, 0.33132171)",class = "b", coef = ""),
  set_prior("normal(0, 0.33132171)", class = "b", coef = "treatECMO"),
  set_prior("student_t(3, 0, 2.5)",class = "Intercept", coef = "") # use default for this prior
)


# define strong prior using neutral prior for control and strong prior for ECMO group
prior_strong <-c(
  set_prior("(0, 0.33132171)",class = "b", coef = ""),
  set_prior("normal(-0.4168554, 0.3177808)", class = "b", coef = "treatECMO"),
  set_prior("student_t(3, 0, 2.5)",class = "Intercept", coef = "") # don't change this prior
)

# define very strong prior using neutral prior for control and very strong prior for ECMO group
prior_vstrong <-c(
  set_prior("(0, 0.33132171)",class = "b", coef = ""),
  set_prior("normal(-0.8612911, 0.3417656)", class = "b", coef = "treatECMO"),
  set_prior("student_t(3, 0, 2.5)",class = "Intercept", coef = "") # don't change this prior
)

# define 3 more models with neutral, strong and very strong priors

# model with a neutral prior
EOLIA_neutral <- brm(dead ~ treat, 
                     bernoulli(link = "logit"), 
                     data = EOLIA,
                     prior = prior_neutral,
                     chains = 2, # nb of chains
                     iter = 5000, # nb of iterations, including burnin
                     warmup = 1000, # burnin
                     thin = 1)

plot(EOLIA_neutral)
# model with a strong prior
EOLIA_strong <- brm(dead ~ treat, 
                    bernoulli(link = "logit"), 
                    data = EOLIA,
                    prior = prior_strong,
                    chains = 2, # nb of chains
                    iter = 5000, # nb of iterations, including burnin
                    warmup = 1000, # burnin
                    thin = 1)
plot(EOLIA_strong)

# model with a very strong prior
EOLIA_vstrong <- brm(dead ~ treat, 
                    bernoulli(link = "logit"), 
                    data = EOLIA,
                    prior = prior_vstrong,
                    chains = 2, # nb of chains
                    iter = 5000, # nb of iterations, including burnin
                    warmup = 1000, # burnin
                    thin = 1)
plot(EOLIA_vstrong)

###################################################
## Both models give warning message after compiling##
######################################################
## The global prior 'normal(0, 0.3)' of class 'b' will not be used in the model,as all related coefficients have individual priors already.
# so extract the priors from the two models and check
prior_summary(EOLIA_neutral)
prior_summary(EOLIA_strong)
prior_summary(EOLIA_vstrong)
# the output suggests the models are using the priors as specified


# get the absolute treatment effect size from the posteriors

# extract the absolute treatment effect using the average slopes function (marginal effects package)
avg_slopes_logistic <- avg_slopes(EOLIA_logistic)
avg_slopes_flat <- avg_slopes(EOLIA_flat)
avg_slopes_neutral <- avg_slopes(EOLIA_neutral)
avg_slopes_strong <- avg_slopes(EOLIA_strong)
avg_slopes_vstrong <- avg_slopes(EOLIA_vstrong)

avg_slopes_logistic # point estimate of treatment effect and 95% confidence interval
avg_slopes_flat # point estimate of treatment effect and 95% credible interval
avg_slopes_neutral # point estimate of treatment effect and 95% credible interval
avg_slopes_strong # point estimate of treatment effect and 95% credible interval
avg_slopes_vstrong # point estimate of treatment effect and 95% credible interval

# get a vector of random draw from the treatment effect to calculate probabilities for brms models
# use the posterior draws function in the marginal effects package to take random samples from posteriors
post_flat <- posterior_draws(avg_slopes_flat)
post_neutral <- posterior_draws(avg_slopes_neutral)
post_strong <- posterior_draws(avg_slopes_strong)
post_vstrong <- posterior_draws(avg_slopes_vstrong)


# extract the column headed 'draw'
post_flat_draw <- post_flat$draw
post_neutral_draw <- post_neutral$draw
post_strong_draw <- post_strong$draw
post_vstrong_draw <- post_vstrong$draw

# check the credible intervals are the same as given by the avg_slopes function
# use ci function in learnBayes

ci(post_flat_draw)
ci(post_neutral_draw)
ci(post_strong_draw)
ci(post_vstrong_draw)
# gives pretty similar results to the earlier method

# now calculate probabilities
# probability the treatment effect is beneficial (ie less than 0)
prob_flat_eff_ben <- ecdf(post_flat_draw)(0)
prob_neutral_eff_ben <- ecdf(post_neutral_draw)(0)
prob_strong_eff_ben <- ecdf(post_strong_draw)(0)
prob_vstrong_eff_ben <- ecdf(post_vstrong_draw)(0)

prob_flat_eff_ben
prob_neutral_eff_ben
prob_strong_eff_ben
prob_vstrong_eff_ben

# probability the treatment effect is larger than a theoretical mcid of -5%
prob_flat_eff_5 <- ecdf(post_flat_draw)(-0.05)
prob_neutral_eff_5 <- ecdf(post_neutral_draw)(-0.05)
prob_strong_eff_5 <- ecdf(post_strong_draw)(-0.05)
prob_vstrong_eff_5 <- ecdf(post_vstrong_draw)(-0.05)

prob_flat_eff_5
prob_neutral_eff_5
prob_strong_eff_5
prob_vstrong_eff_5


# probability the treatment effect is larger than a theoretical mcid of -10%

prob_flat_eff_10 <- ecdf(post_flat_draw)(-0.1)
prob_neutral_eff_10 <- ecdf(post_neutral_draw)(-0.1)
prob_strong_eff_10 <- ecdf(post_strong_draw)(-0.1)
prob_vstrong_eff_10 <- ecdf(post_vstrong_draw)(-0.1)

prob_flat_eff_10
prob_neutral_eff_10
prob_strong_eff_10
prob_vstrong_eff_10
 





# plots

posterior_dist <- data.frame(post_flat_draw, post_neutral_draw, post_strong_draw, post_vstrong_draw)
head(posterior_dist)



ggplot(posterior_dist, aes(post_flat_draw)) +
  geom_density(data = posterior_dist, aes(post_neutral_draw, color = 'brown4'), size = 0.3,fill = 'brown1', alpha = 0.2) +
  geom_density(data = posterior_dist, aes(post_strong_draw, col = 'blue4'), size = 0.3, fill = 'cadetblue', alpha = 0.2) +
  geom_density(data = posterior_dist, aes(post_vstrong_draw, col = 'darkgreen'), size = 0.3, fill = 'lightgreen', alpha = 0.2) +
  theme_bw() +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  labs(title = "Posterior distributions for EOLIA models",
       x= "Absolute Effect Size",
       y = "Density") +
  scale_color_identity(name = "Priors",
                       breaks = c("brown4", "blue4", "darkgreen"),
                       labels = c("Neutral", "Strong", "Very strong"),
                       guide = "legend") +
  geom_vline(xintercept = mean(post_neutral_draw), color = 'red', linetype = 'dotted', size = 0.5)+
  geom_vline(xintercept = mean(post_strong_draw), color = 'blue', linetype = 'dotted', size = 0.5)+
  geom_vline(xintercept = mean(post_vstrong_draw), color = 'darkgreen', linetype = 'dotted', size = 0.5) +
  geom_vline(xintercept = 0, color = 'black', size = 0.3)

  
# composite image
library(patchwork)

Post_brms_uniform <- ggplot(posterior_dist, aes(post_flat_draw)) +
  geom_density(data = posterior_dist, aes(post_flat_draw), color = 'darkgreen',
               size = 0.3,fill = 'chartreuse', alpha = 0.5) +
  geom_vline(xintercept = mean(post_flat_draw), color = 'darkgreen', linetype = 'dotted', size = 0.5)+
  theme_bw() +
  xlim(-0.35, 0.15) +
  ylim(0, 8.5) +
  geom_vline(xintercept = 0, color = 'black', size = 0.3) +
  labs(title = "Uniform",
    x= "Effect Size",
       y = "Density") +
  theme(axis.title=element_text(size=8), plot.title=element_text(size=12)) #change font size of plot title)

Post_brms_uniform
  


Post_brms_neut <- ggplot(posterior_dist, aes(post_neutral_draw)) +
  geom_density(data = posterior_dist, aes(post_neutral_draw), color = 'darkgreen',
               size = 0.3,fill = 'chartreuse', alpha = 0.5) +
  geom_vline(xintercept = mean(post_neutral_draw), color = 'darkgreen', linetype = 'dotted', size = 0.5)+
  theme_bw() +
  xlim(-0.35, 0.15) +
  ylim(0, 8.5) +
  geom_vline(xintercept = 0, color = 'black', size = 0.3) +
  labs(title = "Neutral",
       x= "Effect Size",
       y = "Density") +
  theme(axis.title=element_text(size=8), plot.title=element_text(size=12)) #change font size of plot title)

    
Post_brms_neut
  
Post_brms_strong <- ggplot(posterior_dist, aes(post_strong_draw)) +
  geom_density(data = posterior_dist, aes(post_strong_draw), color = 'darkgreen',
               size = 0.3,fill = 'chartreuse', alpha = 0.5) +
  geom_vline(xintercept = mean(post_strong_draw), color = 'darkgreen', linetype = 'dotted', size = 0.5)+
  theme_bw() +
  xlim(-0.35, 0.15) +
  ylim(0, 8.5) +
  geom_vline(xintercept = 0, color = 'black', size = 0.3) +
  labs(title = "Enthusiastic",
       x= "Effect Size",
       y = "Density") +
  theme(axis.title=element_text(size=8), plot.title=element_text(size=12)) #change font size of plot title)

Post_brms_strong

Post_brms_vstrong <- ggplot(posterior_dist, aes(post_vstrong_draw)) +
  geom_density(data = posterior_dist, aes(post_vstrong_draw), color = 'darkgreen',
               size = 0.3,fill = 'chartreuse', alpha = 0.5) +
  geom_vline(xintercept = mean(post_vstrong_draw), color = 'darkgreen', linetype = 'dotted', size = 0.5)+
  theme_bw() +
  xlim(-0.35, 0.15) +
  ylim(0, 8.5) +
  geom_vline(xintercept = 0, color = 'black', size = 0.3) +
  labs(title = "Very enthusiastic",
       x= "Effect Size",
       y = "Density") +
  theme(axis.title=element_text(size=8), plot.title=element_text(size=12)) #change font size of plot title)

Post_brms_vstrong

Post_brms_uniform + Post_brms_neut + Post_brms_strong + Post_brms_vstrong + plot_layout(ncol = 2, byrow = TRUE)
