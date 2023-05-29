rm(list = ls())

# Comparing proportions using the Beta-binomial model
library(tidyverse)
library(bayestestR) # ci function
library(qwraps2) # logit and invlogit functions
library(ProbBayes) # beta.select function
library(LearnBayes)
library(patchwork)



##################
## trial data ###
#################

# 57/125 dead in control
# 44/125 dead in ECMO group


# treatment group
# y2 = number experiencing an event
# n2 = number in control group
y1 = 44
n1 = 125

# control group
# y1 = number experiencing an event
# n1 = number in control group
y2 = 57
n2 = 125


####################
# Flat prior###
###################

# vector of random draws for the prior for the control group
prior_con_draws_bb <- rbeta(10000, 1, 1)
mean(prior_con_draws_bb)

# vector of random draws for the prior for the treatment group
prior_treat_draws_bb <- rbeta(10000, 1, 1)

# vector of random draws for the posterior for the control group
post_con_draws_bb <- rbeta(10000, 1 + y2, n2-y2 + 1)


# vector of random draws for the posterior for the treatment group
post_treat_draws_bb <- rbeta(10000, 1 + y1, n1-y1 + 1)



# vector for the prior treatment effect
prior_effect_dist_bb <- prior_treat_draws_bb - prior_con_draws_bb

# vector for posterior treatment effect
post_effect_dist_bb <- post_treat_draws_bb - post_con_draws_bb

# 95% credible interval for the treatment effect
ci(post_effect_dist_bb)

# probabilities for the treatment effect
ecdf(post_effect_dist_bb)(0)
ecdf(post_effect_dist_bb)(-0.05)
ecdf(post_effect_dist_bb)(-0.1)

###################
## neutral ##
###################
# beta prior for the control group
# anticipated median event rate in the control group
x1_con_neut = 0.5
# anticipated 90th centile for the event rate in the control group
x2_con_neut = 0.6
# anticipated median event rate in the treatment group
x1_treat_neut = 0.5
# anticipated 90th centile for the event rate in the treatment group
x2_treat_neut = 0.6


quCon1_neut=list(p=.5,x=x1_con_neut)
quCon2_neut=list(p=.9,x=x2_con_neut)
prior_con_neut <- beta.select(quCon1_neut, quCon2_neut)
prior_con_neut[1]
prior_con_neut[2]

quTreat1_neut=list(p=.5,x=x1_treat_neut)
quTreat2_neut=list(p=.9,x=x2_treat_neut)
prior_treat_neut <- beta.select(quTreat1_neut, quTreat2_neut)
prior_treat_neut[1]
prior_treat_neut[2]

# vector of random draws for the prior for the control group
prior_con_draws_neut <- rbeta(10000, prior_con_neut[1], prior_con_neut[2])
mean(prior_con_draws_neut)

# vector of random draws for the prior for the treatment group
prior_treat_draws_neut <- rbeta(10000, prior_treat_neut[1], prior_treat_neut[2])

# vector of random draws for the posterior for the control group
post_con_draws_neut <- rbeta(10000, prior_con_neut[1] + y2, n2-y2 + prior_con_neut[2])


# vector of random draws for the posterior for the treatment group
post_treat_draws_neut <- rbeta(10000, prior_treat_neut[1] + y1, n1-y1 + prior_treat_neut[2])



# vector for the prior treatment effect
prior_effect_dist_neut <- prior_treat_draws_neut - prior_con_draws_neut


# vector for posterior treatment effect
post_effect_dist_neut <- post_treat_draws_neut - post_con_draws_neut


# 95% credible interval for the treatment effect
ci(post_effect_dist_neut)

# probabilities for the treatment effect
ecdf(post_effect_dist_neut)(0)
ecdf(post_effect_dist_neut)(-0.05)
ecdf(post_effect_dist_neut)(-0.1)


###################
## strong ########
###################

# beta prior for the control group
# anticipated median event rate in the control group
x1_con_strong = 0.5
# anticipated 90th centile for the event rate in the control group
x2_con_strong = 0.6
# anticipated median event rate in the treatment group
x1_treat_strong = 0.4
# anticipated 90th centile for the event rate in the treatment group
x2_treat_strong = 0.5


quCon1_strong=list(p=.5,x=x1_con_strong)
quCon2_strong=list(p=.9,x=x2_con_strong)
prior_con_strong <- beta.select(quCon1_strong, quCon2_strong)
prior_con_strong[1]
prior_con_strong[2]

quTreat1_strong=list(p=.5,x=x1_treat_strong)
quTreat2_strong=list(p=.9,x=x2_treat_strong)
prior_treat_strong <- beta.select(quTreat1_strong, quTreat2_strong)
prior_treat_strong[1]
prior_treat_strong[2]

# vector of random draws for the prior for the control group
prior_con_draws_strong <- rbeta(10000, prior_con_strong[1], prior_con_strong[2])
mean(prior_con_draws_strong)

# vector of random draws for the prior for the treatment group
prior_treat_draws_strong <- rbeta(10000, prior_treat_strong[1], prior_treat_strong[2])

# vector of random draws for the posterior for the control group
post_con_draws_strong <- rbeta(10000, prior_con_strong[1] + y2, n2-y2 + prior_con_strong[2])


# vector of random draws for the posterior for the treatment group
post_treat_draws_strong <- rbeta(10000, prior_treat_strong[1] + y1, n1-y1 + prior_treat_strong[2])



# vector for the prior treatment effect
prior_effect_dist_strong <- prior_treat_draws_strong - prior_con_draws_strong


# vector for posterior treatment effect
post_effect_dist_strong <- post_treat_draws_strong - post_con_draws_strong


# 95% credible interval for the treatment effect
ci(post_effect_dist_strong)

# probabilities for the treatment effect
ecdf(post_effect_dist_strong)(0)
ecdf(post_effect_dist_strong)(-0.05)
ecdf(post_effect_dist_strong)(-0.1)


###################
## very strong ########
###################

# beta prior for the control group
# anticipated median event rate in the control group
x1_con_vstrong = 0.5
# anticipated 90th centile for the event rate in the control group
x2_con_vstrong = 0.6
# anticipated median event rate in the treatment group
x1_treat_vstrong = 0.3
# anticipated 90th centile for the event rate in the treatment group
x2_treat_vstrong = 0.4


quCon1_vstrong=list(p=.5,x=x1_con_vstrong)
quCon2_vstrong=list(p=.9,x=x2_con_vstrong)
prior_con_vstrong <- beta.select(quCon1_vstrong, quCon2_vstrong)
prior_con_vstrong[1]
prior_con_vstrong[2]

quTreat1_vstrong=list(p=.5,x=x1_treat_vstrong)
quTreat2_vstrong=list(p=.9,x=x2_treat_vstrong)
prior_treat_vstrong <- beta.select(quTreat1_vstrong, quTreat2_vstrong)
prior_treat_vstrong[1]
prior_treat_vstrong[2]

# vector of random draws for the prior for the control group
prior_con_draws_vstrong <- rbeta(10000, prior_con_vstrong[1], prior_con_vstrong[2])
mean(prior_con_draws_vstrong)

# vector of random draws for the prior for the treatment group
prior_treat_draws_vstrong <- rbeta(10000, prior_treat_vstrong[1], prior_treat_vstrong[2])

# vector of random draws for the posterior for the control group
post_con_draws_vstrong <- rbeta(10000, prior_con_vstrong[1] + y2, n2-y2 + prior_con_vstrong[2])


# vector of random draws for the posterior for the treatment group
post_treat_draws_vstrong <- rbeta(10000, prior_treat_vstrong[1] + y1, n1-y1 + prior_treat_vstrong[2])



# vector for the prior treatment effect
prior_effect_dist_vstrong <- prior_treat_draws_vstrong - prior_con_draws_vstrong


# vector for posterior treatment effect
post_effect_dist_vstrong <- post_treat_draws_vstrong - post_con_draws_vstrong


# 95% credible interval for the treatment effect
ci(post_effect_dist_vstrong)

# probabilities for the treatment effect
ecdf(post_effect_dist_vstrong)(0)
ecdf(post_effect_dist_vstrong)(-0.05)
ecdf(post_effect_dist_vstrong)(-0.1)


# posterior vectors
post_effect_dist_bb
post_effect_dist_neut
post_effect_dist_strong
post_effect_dist_vstrong

# plots

post_distributions_betabinomial <- data.frame(post_effect_dist_bb, post_effect_dist_neut, 
                                              post_effect_dist_strong, post_effect_dist_vstrong)
head(post_distributions_betabinomial)

# composite image

Post_bb_uniform <- ggplot(post_distributions_betabinomial, aes(post_effect_dist_bb)) +
  geom_density(data = post_distributions_betabinomial, aes(post_effect_dist_bb), color = 'darkgreen',
               size = 0.3,fill = 'chartreuse', alpha = 0.5) +
  geom_vline(xintercept = mean(post_effect_dist_bb), color = 'darkgreen', linetype = 'dotted', size = 0.5)+
  theme_bw() +
  xlim(-0.35, 0.15) +
  ylim(0, 8.5) +
  geom_vline(xintercept = 0, color = 'black', size = 0.3) +
  labs(title = "Uniform",
       x= "Effect Size",
       y = "Density") +
  theme(axis.title=element_text(size=8), plot.title=element_text(size=12)) #change font size of plot title)

Post_bb_uniform


Post_bb_neut <- ggplot(post_distributions_betabinomial, aes(post_effect_dist_neut)) +
  geom_density(data = post_distributions_betabinomial, aes(post_effect_dist_neut), color = 'darkgreen',
               size = 0.3,fill = 'chartreuse', alpha = 0.5) +
  geom_vline(xintercept = mean(post_effect_dist_neut), color = 'darkgreen', linetype = 'dotted', size = 0.5)+
  theme_bw() +
  xlim(-0.35, 0.15) +
  ylim(0, 8.5) +
  geom_vline(xintercept = 0, color = 'black', size = 0.3) +
  labs(title = "Neutral",
       x= "Effect Size",
       y = "Density") +
  theme(axis.title=element_text(size=8), plot.title=element_text(size=12)) #change font size of plot title)


Post_bb_neut

Post_bb_strong <- ggplot(post_distributions_betabinomial, aes(post_effect_dist_strong)) +
  geom_density(data = post_distributions_betabinomial, aes(post_effect_dist_strong), color = 'darkgreen',
               size = 0.3,fill = 'chartreuse', alpha = 0.5) +
  geom_vline(xintercept = mean(post_effect_dist_strong), color = 'darkgreen', linetype = 'dotted', size = 0.5)+
  theme_bw() +
  xlim(-0.35, 0.15) +
  ylim(0, 8.5) +
  geom_vline(xintercept = 0, color = 'black', size = 0.3) +
  labs(title = "Enthusiastic",
       x= "Effect Size",
       y = "Density") +
  theme(axis.title=element_text(size=8), plot.title=element_text(size=12)) #change font size of plot title)

Post_bb_strong

Post_bb_vstrong <- ggplot(post_distributions_betabinomial, aes(post_effect_dist_vstrong)) +
  geom_density(data = post_distributions_betabinomial, aes(post_effect_dist_vstrong), color = 'darkgreen',
               size = 0.3,fill = 'chartreuse', alpha = 0.5) +
  geom_vline(xintercept = mean(post_effect_dist_vstrong), color = 'darkgreen', linetype = 'dotted', size = 0.5)+
  theme_bw() +
  xlim(-0.35, 0.15) +
  ylim(0, 8.5) +
  geom_vline(xintercept = 0, color = 'black', size = 0.3) +
  labs(title = "Very enthusiastic",
       x= "Effect Size",
       y = "Density") +
  theme(axis.title=element_text(size=8), plot.title=element_text(size=12)) #change font size of plot title)

Post_bb_vstrong

Post_bb_uniform + Post_bb_neut + Post_bb_strong + Post_bb_vstrong + plot_layout(ncol = 2, byrow = TRUE)

