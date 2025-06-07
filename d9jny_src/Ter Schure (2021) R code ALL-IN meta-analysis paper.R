### Author: Judith ter Schure
### E-mail: schure@cwi.nl
### Institute: CWI - Machine Learning group
###            CWI is the national research institute 
###            for mathematics and computer science in the Netherlands
###            Science Park 123, 1098 XG Amsterdam, NETHERLANDS
### Date: 22 September 2021
### Licence: CC0

### This code is part of the publication:
### Ter Schure & Grünwald (2021) ALL-IN meta-analysis: breathing life into living systematic reviews

### Details on the OS, system and software that produced the results
### in the paper:

# Platform: x86_64-w64-mingw32/x64 (64-bit)
# OS: Microsoft Windows
# System: Windows 10 x64 (build 19042)
# R version: 3.6.1 (2019-07-05)

### The following packages were used:
# devtools version 2.2.2
# usethis version 1.6.1
# safestats version 0.8.6
# ggplot2 version 3.3.5
# ggside version 0.1.1
# ggpointdensity version 0.1.0
# stringr version 1.4.3

# stats version 3.6.1
# graphics version 3.6.1
# grDevices version 3.6.1
# utils version 3.6.1
# methods version 3.6.1


######### Content #########
######### 
######### 
######
###### 1. Settings of the FDA game
###### 2. Pfizer/BioNtech results in the FDA game
###### 3. CureVac results in terms of confidence interval and the FDA game
###### 4. Plotting the FDA game betting scores under the null hypothesis (Figure 1)
###### 5. Plotting the expected sample size for various strategies in the FDA game (Figure 2)
###### 6. Plotting the anytime-valid confidence sequence for a random ordering of the CureVac data (Figure 3)
###### 7. Plotting the implied target of the CureVac design in the FDA game (Figure 4 & 5)
###### 8. Plotting a p-value of 0.05 (Figure 6; Appendix)
######
#########
#########
###########################

library(devtools)
devtools::install_github("AlexanderLyNL/Safestats", ref = "logrank")
library(safestats)
library(ggplot2)
library(ggside)
library(ggpointdensity)
library(stringr)

# Colours for plots
f1000blue <- "#006a89"
f1000orange <- "#cc6638"
f1000grey <- "#374e59"



######### 
######
###### 1. Settings of the FDA game
###### 
#########

calcPfromVE <- function(VE) {
  (100 - VE)/(100 + (100 - VE))  # calculate the probability that the next event
}                                # is in the vaccine group based on the Vaccine Efficacy (VE)
calcVEfromP <- function(p) {  # calculate VE from probability of next event in vaccine group
  RR = p/(1 - p)     # relative risk or (constant) hazard ratio
  VE = (1 - RR)*100  # Vaccine Efficacy (VE)
}

VE0 <- 30 # <= 30% VE under the null hypothesis
pVacc0 <- calcPfromVE(VE0)
pPlac0 <- 1 - pVacc0

pVacc0  # smallest probability of the next event in the vaccine group under the null hypothesis
pPlac0  # largest probability of the next event in the placebo group under the null hypothesis

multiplyerVacc <- 1/pVacc0
multiplyerPlac <- 1/pPlac0

multiplyerVacc  # analogue for the doubling factor in roulette: multiplyer in the FDA game for vaccine event
multiplyerPlac  # analogue for the doubling factor in roulette: multiplyer in the FDA game for placebo event

VE1 <- 50 # >=50% VE under the alternative hypothesis
pVacc1 <- calcPfromVE(VE1)  # probabiliity of the next event in the vaccine group under the alt hypothesis
pPlac1 <- 1 - pVacc1   # probabiliity of the next event in the placebo group under the alt hypothesis

pVacc1  # largest  probability of the next event in the vaccine group under the alternative hypothesis
pPlac1  # smallest probability of the next event in the placebo group under the alternative hypothesis


# you multiply your initial capital by this likelihood ratio, if you bet the FDA game according to p1 (= pVacc1)
calcLRfromP <- function(group, p1 = pVacc1, p0 = pVacc0) {
  if (group == "vaccine") p1/p0
  else (1 - p1)/(1 - p0)
}
calcLRfromVE <- function(group, VE1, VE0) {
  calcLRfromP(group, p1 = calcPfromVE(VE1), p0 = calcPfromVE(VE0))
}

calcLR <- function(group){
  calcLRfromVE(group, 50, 30)
}


######### 
######
###### 2. Pfizer/BioNtech results in the FDA game
######
#########

# https://www.nejm.org/doi/full/10.1056/NEJMoa2034577
PfiBioNobs1 <- 8    # number of observed events of Covid-19 infection in the vaccine group
PfiBioNobs0 <- 162  # number of observed events of Covid-19 infection in the placebo control group

# Betting score Pfizer/BioNTech trial
PfizerBioNtech <- calcLR("vaccine")^PfiBioNobs1 *
                  calcLR("placebo")^PfiBioNobs0
PfizerBioNtech/1000000  # in million ???



######### 
######
###### 3. Curevac results in terms of confidence interval and the FDA game
######
#########

# https://www.curevac.com/en/2021/06/30/curevac-final-data-from-phase-2b-3-trial-of-first-generation-covid-19-vaccine-candidate-cvncov-demonstrates-protection-in-age-group-of-18-to-60/
CureVobs1 <- 83  # number of observed events of Covid-19 infection in the vaccine group
CureVobs0 <- 145  # number of observed events of Covid-19 infection in the placebo control group

# Quote below from p. 124 of protocol:
# https://www.curevac.com/wp-content/uploads/2021/06/HERALD_CV-NCOV-004-Protocol.pdf

#"
# VE = 1 - RR = 1 - (ARV/ARP) = 1 - {p / r (1-p)} where 
# ARV = attack rate in vaccinated group = nv/Nv = number of subjects reporting at least one 
# COVID-19 episode in the CVnCoV group / total follow-up time of evaluable subjects in the 
# CVnCoV group (number of person-month). 
# ARP = attack rate in placebo group = np/Np = number of subjects reporting at least one 
# COVID-19 episode in the placebo group / total follow-up time of evaluable subjects in the 
# placebo group (number of person-month). 
# RR = relative risk = ARV/ARP 
# p = proportion of COVID-19 cases (according to primary case definition) coming from the 
# CVnCoV group among all cases = nv/(nv+np). 
# r = ratio of total follow-up time of evaluable subjects in the CVnCoV group over total follow-up time 
# of evaluable subjects in the placebo group = Nv/Np
#"

# We assume that the risk set is balanced troughout the trial, so r = 1 and Nv = Np, 
# so RR = (nv/Nv)/(np/Np) = nv/np. nv = CureVobs1 and np = CureVobs0. So:

CureVRR <- CureVobs1/CureVobs0
CureVVE <- (1 - CureVRR)*100
CureVVE  # estimated VE based on the CureVac data
# The press release reports a VE of 48%, so uses a different r (ratio of follow-up time in the two groups). 
# In such large trials with 50:50 randomization r can be assumed to stay close to 1, 
# so we set it to 1 to make all calculations simpler.

CureValphaFinal <- 0.02281  # from the CureVac protocol linked above (Table 8)
CureVZalphaFinal <- abs(qnorm(CureValphaFinal))

CureVp <- calcPfromVE(CureVVE)  # estimated probability that the next event is in the vaccine group
CureVse.p <- sqrt(CureVp*(1 - CureVp)/(CureVobs1 + CureVobs0))
CureVinterval.p <- c(CureVp - CureVZalphaFinal*CureVse.p, CureVp + CureVZalphaFinal*CureVse.p)
CureVintervalVE <- sort(calcVEfromP(CureVinterval.p))
CureVintervalVE

# Betting score CureVac trial
CureVac <- calcLR("vaccine")^CureVobs1 *
           calcLR("placebo")^CureVobs0
CureVac


######### 
######
###### 4. Plotting the FDA game betting scores under the null hypothesis (Figure 1)
######
#########

set.seed(2021)

#numSim <- 1000
numSim <- 100  # for faster plotting
numBets <- PfiBioNobs1 + PfiBioNobs0  # sample size from the Pfizer/BioNTech trial
alpha = 0.025


sampleEvents <- function(numBets) {
  sample(c("vaccine", "control"), 
         size = numBets, 
         replace = TRUE, 
         prob = c(pVacc0, pPlac0))
}

betSamplePayout <- function(numBets) {
  sapply(sampleEvents(numBets), calcLR)
}

# betting scores by round
betScByRound <- data.frame(simulation = rep(paste("sim", 1:numSim), each = numBets),
                           bettingScores = unlist(lapply(1:numSim, function(x) cumprod(betSamplePayout(numBets)))),
                           bettingRound = rep(1:numBets, times = numSim))

# final betting scores after numBets rounds
finbetScs <- data.frame(finalbettingScore = betScByRound$bettingScores[betScByRound$bettingRound == numBets],
                        bettingRound = rep(numBets, times = numSim))

ggplot() + 
  geom_line(aes(y = bettingScores, 
                x = bettingRound, 
                group = simulation), 
            colour = f1000blue, lwd = 0.05,
            data = betScByRound, show.legend = FALSE) +
  stat_pointdensity(aes(y = bettingScores, 
                        x = bettingRound), 
                    adjust = min(betScByRound$bettingScores),
                    size = 0.1, data = betScByRound) +
  scale_color_gradientn(colors = c(f1000grey, f1000blue, "black", "white")) +
  geom_jitter(aes(y = finalbettingScore, 
                  x = bettingRound), 
              colour = f1000orange, size = 0.5,
              width = 0.001, data = finbetScs) +
  geom_ysidehistogram(aes(y = finalbettingScore), fill = f1000orange, colour = NA,
                      binwidth = (max(finbetScs$finalbettingScore) - min(finbetScs$finalbettingScore))/max(finbetScs$finalbettingScore)/2,
                      data = finbetScs) +
  geom_hline(colour = f1000grey, yintercept = 1/alpha, size = 1, lty = 2) +
  scale_y_continuous(trans = 'log2', breaks = 2^(seq(from = -20, to = 20, by = 2)), 
                     labels = c(paste("1/", 2^(seq(from = 20, to = 2, by = -2)), sep = ""), 1, 
                                2^seq(from = 2, to = 20, by = 2)),
                     sec.axis = sec_axis(trans = function(x) x, breaks = 2^(seq(from = -20, to = 20, by = 2)), 
                                         labels = c(paste("1/", 2^(seq(from = 20, to = 2, by = -2)), sep = ""), 1, 
                                                    2^seq(from = 2, to = 20, by = 2)))) +
  scale_ysidex_continuous(labels = NULL) +
  labs(colour = "density", x = "betting round / number of events (n)", y = bquote("??? /" ~ LR^(n)),
       subtitle = bquote("Simulated betting scores ??? /" ~ LR^(n) ~ " over betting rounds (events) n under the null hypothesis of 30% VE")) +
  guides(colour = guide_colorbar(title = "density", label = FALSE, barheight = unit(10, "mm"))) +
  theme(text = element_text(size = 9), legend.position = c(0.02, 0.05), legend.justification = c(0.02, 0.05))

ggsave(filename = "bettingScoreSequence.pdf",
       height   = 70,
       width    = 210 - 1.4*25.4, 
       units    = "mm")

# This plot is not shown in the paper
ggplot() +
  geom_histogram(aes(x = finalbettingScore), 
                 fill = f1000orange, colour = NA,
                 binwidth = (max(finbetScs$finalbettingScore) - min(finbetScs$finalbettingScore))/max(finbetScs$finalbettingScore)/2, 
                 data = finbetScs) +
  scale_x_continuous(trans = 'log2', breaks = 2^(seq(from = -20, to = 20, by = 2)), 
                     labels = c(paste("1/", 2^(seq(from = 20, to = 2, by = -2)), sep = ""), 1, 
                                2^seq(from = 2, to = 20, by = 2))) +
  scale_y_continuous(labels = NULL) +
  geom_vline(xintercept = 1/alpha, size = 1, lty = 2, colour = f1000grey) +
  labs(y = NULL, x = bquote("??? /" ~ LR^(170)), subtitle = bquote("Distribution of simulated final betting scores ??? /" ~ LR^(170) ~ "after 170 events under the null hypothesis of 30% VE")) +
  theme(text = element_text(size = 9))

ggsave(filename = "bettingScoreDistribution.pdf",
       height   = 50,
       width    = 210 - 1.5*25.4, 
       units    = "mm")


# proportion simulated type-I errors at 170 events, and before 170 events
mean(finbetScs$finalbettingScore > 1/alpha)
mean(sapply(unique(betScByRound$simulation), function(sim) 
  any(betScByRound$bettingScores[betScByRound$simulation == sim] > 1/alpha)))



######### 
######
###### 5. Plotting the expected sample size for various strategies in the FDA game (Figure 2)
######
#########

trueVEs <- list("40" = 40:100,
                "50" = 42.8:100,
                "60" = 48:100)
VE1s <- c(40, 50, 60)

calcElogLR <- function(trueVE, VE1) {
  calcPfromVE(trueVE)*log(calcLR("vaccine")) +
  (1 - calcPfromVE(trueVE))*log(calcLR("placebo"))  
}

calcN <- function(trueVE, VE1, alpha) {
  log(1/alpha)/calcElogLR(trueVE, VE1)  # Using Wald's identity
}

data <- data.frame("trueVE" = unlist(trueVEs),
                   "N" = unlist(lapply(VE1s, function(VE1) 
                     sapply(trueVEs[[as.character(VE1)]], function(trueVE) calcN(trueVE, VE1, alpha = 0.025)))),
                   "VE1" = as.factor(c(rep(paste(VE1s[1], "%", sep = ""), times = length(trueVEs[["40"]])),
                                       rep(paste(VE1s[2], "%", sep = ""), times = length(trueVEs[["50"]])),
                                       rep(paste(VE1s[3], "%", sep = ""), times = length(trueVEs[["60"]])))))

ggplot(data) +
  geom_line(aes(x = trueVE, y = N, col = VE1, linetype = VE1)) +
  scale_x_continuous(breaks = c(40, 50, 60, 70, 80, 90, 100), 
                     labels = paste(c(40, 50, 60, 70, 80, 90, 100), "%", sep = "")) +
  scale_y_continuous(breaks = c(250, 500, 750, 1000, 1250)) +
  scale_colour_manual(values = c(f1000blue, f1000orange, f1000grey)) +
  labs(y = bquote(N[alpha]), x = "true VE", colour = bquote(VE[1]), linetype = bquote(VE[1]),
       subtitle = bquote("Expected" ~ N ~ "needed to achieve" ~ 1/alpha ~ "= 40")) +
  theme(text = element_text(size = 9))


ggsave(filename = "expectedSampleSize.pdf",
       height   = 50,
       width    = 85, 
       units    = "mm")

######### 
######
###### 6. Plotting the anytime-valid confidence sequence for a random ordering of the CureVac data (Figure 3)
######
#########

set.seed(21)

# sample a random sequence of the CureVac observations
CureVobs1s <- sample(c(rep(1, times = CureVobs1), rep(0, times = CureVobs0)), replace = FALSE)

# calculate logrank Z-statistics assuming balanced risk set (in O-E, E is always 0.5)
zs <- sapply(1:(CureVobs1 + CureVobs0), function(n) (sum(CureVobs1s[1:n]) - n/2)/sqrt(n/4))

# specify the design of the safe/e-value logrank test used to optimize the confidence sequence
designObj <- designSafeLogrank(hrMin = 0.5, exact = FALSE)  # based on the Gaussian approximation (exact = FALSE)

confInt <- function(z, n, hrMin = 0.5) {
  confIntHR <- safeLogrankTestStat(z = z, nEvents = n, designObj = designObj, ciValue = 0.90)$confSeq
  return(confIntHR)  # Peto estimator based on the logrank statistic
}

confSeq <- t(sapply(1:(CureVobs1 + CureVobs0), function(n) confInt(z = zs[n], n = n, hrMin = 0.5)))
confSeq <- as.data.frame(cbind(confSeq, 1:(CureVobs1 + CureVobs0)))
colnames(confSeq) <- c("lower", "upper", "n")
confSeq$lowerRunning <- cummax(confSeq$lower)
confSeq$upperRunning <- cummin(confSeq$upper)

          
# final confidence interval
tail(confSeq)
sort((1 - confSeq[CureVobs1 + CureVobs0, c("lowerRunning", "upperRunning")])*100)

myColours <- c("AVconfInt" = f1000blue, 
               "runInt" = f1000orange)
# remove very wide initial intervals for plotting purposes:               
confSeq$lower[confSeq$lower < 0.1] <- 0.1
confSeq$lowerRunning[confSeq$lowerRunning < 0.1] <- NA
confSeq$upper[confSeq$upper > 20] <- 20
confSeq$upperRunning[confSeq$upperRunning > 20] <- NA

ggplot(confSeq) + 
  geom_segment(aes(x = n, xend = n, y = lower, yend = upper, colour = "AVconfInt"), lwd = 0.1) +
  geom_line(aes(x = n, y = lowerRunning, colour = "runInt")) + 
  geom_line(aes(x = n, y = upperRunning, colour = "runInt")) +
  scale_y_continuous(trans = 'log2', #limits = c(1/10^20, 21),
                     breaks = c(0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 1, 
                                round(1/c(0.7, 0.5, 0.35, 0.2, 0.1, 0.05), digits = 1)),
                     labels = c(0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 1, 
                                round(1/c(0.7, 0.5, 0.35, 0.2, 0.1, 0.05), digits = 1)),
                     sec.axis = sec_axis(trans = function(x) x,
                                         breaks = c(0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 1),
                                         labels = c("95%", "90%", "80%", "65%", "50%", "30%", "0%"),
                                         name = "estimated VE")) +
  scale_colour_manual(values = myColours, labels = c(str_wrap("90% anytime-valid confidence interval", 20),
                                                     str_wrap("running intersection", 20))) +
                        #c(str_wrap("90% anytime-valid confidence interval", 15), str_wrap("running intersection", 60))) +
  labs(y = bquote("estimated hazard ratio"), colour = "",
       x = "betting round / number of events (n)") +
  theme(text = element_text(size = 9), legend.position = c(0.7, 0.9),
              legend.key.size = unit(5, 'mm'))
        
ggsave(filename = "confSeq.pdf",
               height   = 75,
               width    = 80, 
               units    = "mm")


######### 
######
###### 7. Plotting the implied target of the CureVac design in the FDA game (Figure 4 & 5)
######
#########

set.seed(2020)

#numSim <- 1000
numSim <- 100  # for faster plotting
numBets <- 160

sampleEvents <- function(numBets) {
  sample(c("vaccine", "control"), 
         size = numBets, 
         replace = TRUE, 
         prob = c(calcPfromVE(60), 1 - calcPfromVE(60)))
}

betSamplePayout <- function(numBets) {
  sapply(sampleEvents(numBets), function(group) calcLR(group = group))
}

betScByRound <- data.frame(simulation = rep(paste("sim", 1:numSim), each = numBets),
                             bettingScores = unlist(lapply(1:numSim, function(x)
                               cumprod(betSamplePayout(numBets)))),
                             bettingRound = rep(1:numBets, times = numSim))

finbetScs <- data.frame(finalbettingScore = betScByRound$bettingScores[betScByRound$bettingRound == numBets],
                           bettingRound = rep(numBets, times = numSim))

impliedTarget <- exp(160*((40/140)*log((50/150)/(70/170)) + (100/140)*log((100/150)/(100/170))))

# which calculates:
exp(160 * (     calcPfromVE(60)  * log(calcLRfromVE("vaccine", 50, 30)) +
            (1 - calcPfromVE(60)) * log(calcLRfromVE("placebo", 50, 30))
          )
    )

ggplot() + 
  geom_line(aes(y = bettingScores, 
                x = bettingRound, 
                group = simulation), 
            colour = f1000orange, lwd = 0.05,
            data = betScByRound) +
  stat_pointdensity(aes(y = bettingScores, x = bettingRound), adjust = min(betScByRound$bettingScores),
                    size = 0.1, data = betScByRound) +
  scale_color_gradientn(colors = c(f1000grey, f1000orange, "black", "white")) +
  geom_jitter(aes(y = finalbettingScore, 
                  x = bettingRound),
              colour = f1000blue, size = 0.5,
              width = 0.001, data = finbetScs) +
  geom_ysidehistogram(aes(y = finalbettingScore), 
                      fill = f1000blue, colour = NA,
                      binwidth = (max(finbetScs$finalbettingScore) - min(finbetScs$finalbettingScore))/max(finbetScs$finalbettingScore)/2,
                      data = finbetScs) +
  geom_hline(colour = f1000grey,
             yintercept = 1/alpha, size = 1, lty = 2) +
  geom_segment(aes(x = min_x, 
                   y = min_y, 
                   xend = max_x,
                   yend = max_y), data = data.frame(min_x = 0, 
                                                    min_y = 1, 
                                                    max_x = 160,
                                                    max_y = impliedTarget),
               size = 1, lty = "solid", colour = f1000orange) +
  scale_y_continuous(trans = 'log2', breaks = 2^(seq(from = -20, to = 20, by = 2)), 
                     labels = c(paste("1/", 2^(seq(from = 20, to = 2, by = -2)), sep = ""), 1, 
                                2^seq(from = 2, to = 20, by = 2)),
                     sec.axis = sec_axis(trans = function(x) x, breaks = 2^(seq(from = -20, to = 20, by = 2)), 
                                         labels = c(paste("1/", 2^(seq(from = 20, to = 2, by = -2)), sep = ""), 1, 
                                                    2^seq(from = 2, to = 20, by = 2)))) +
  scale_ysidex_continuous(labels = NULL) +
  labs(colour = "", x = "betting round / number of events (n)", y = bquote("??? /" ~ LR^(n)),
       subtitle = bquote("Simulated betting scores ??? /" ~ LR^(n) ~ " over betting rounds (events) n under the alternative hypothesis of 60% VE")) +
  guides(colour = guide_colorbar(title = "density", label = FALSE, barheight = unit(10, "mm"))) +
  theme(text = element_text(size = 9), legend.position = c(0.05, 0.95), legend.justification = c(0.05, 0.95))


ggsave(filename = "impliedTargetSequence.pdf",
       height   = 70,
       width    = 210 - 1.4*25.4, 
       units    = "mm")

ggplot() +
  geom_histogram(aes(x = finalbettingScore), 
                 fill = f1000blue, colour = NA,
                 binwidth = (max(finbetScs$finalbettingScore) - min(finbetScs$finalbettingScore))/max(finbetScs$finalbettingScore)/2,
                 data = finbetScs) +
  scale_x_continuous(trans = 'log2', breaks = 2^(seq(from = -20, to = 20, by = 2)), 
                     labels = c(paste("1/", 2^(seq(from = 20, to = 2, by = -2)), sep = ""), 1, 
                                2^seq(from = 2, to = 20, by = 2))) +
  scale_y_continuous(labels = NULL) +
  geom_vline(xintercept = 1/alpha, size = 1, lty = 2, colour = f1000grey) +
  geom_vline(xintercept = exp(((40/140)*log((50/150)/(70/170)) + (100/140)*log((100/150)/(100/170)))*160), 
             size = 1, lty = "solid", colour = f1000orange) +
  labs(y = NULL, x = bquote("??? /" ~ LR^(160)), subtitle = bquote("Distribution of simulated final betting scores ??? /" ~ LR^(160) ~ "after 160 events under the alternative hypothesis of 60% VE")) +
  theme(text = element_text(size = 9))

ggsave(filename = "impliedTargetDistribution.pdf",
       height   = 50,
       width    = 210 - 1.4*25.4, 
       units    = "mm")

exp(mean(log(finbetScs$finalbettingScore)))  # simulated implied target
mean(finbetScs$finalbettingScore > 1/alpha)  # power at 160 events and before 160 events:
mean(sapply(unique(betScByRound$simulation), function(sim) any(betScByRound$bettingScores[betScByRound$simulation == sim] > 1/alpha)))



######### 
######
###### 8. Plotting a p-value of 0.05 (Figure 6; Appendix)
######
#########


alpha <- 0.025
myxlim <- c(-4, 4)

ZsL <- seq(from = -4, to = -1.96, by = 0.01)
ZsR <- seq(from = 1.96, to = 4, by = 0.01)
pValueLeft <- data.frame(x = c(ZsL, rev(ZsL)),
                         y = c(dnorm(ZsL), rep(0, length(ZsL))))
pValueRight <- data.frame(x = c(ZsR, rev(ZsR)),
                          y = c(dnorm(ZsR), rep(0, length(ZsR))))
ggplot() +
  geom_function(fun = dnorm, colour = f1000orange, lwd = 1) +
  xlim(-4, 4) +
  geom_polygon(aes(x, y, fill = "pValue"), data = pValueLeft) +
  geom_polygon(aes(x, y, fill = "pValue"), data = pValueRight) +
  scale_fill_manual(values = c("pValue" = f1000blue), label = "p-value") +
  labs(x = "Z-statistic", y = bquote(phi[0]), fill = "") +
  theme(text = element_text(size = 9), legend.position = c(0.86, 0.95),
        legend.key.size = unit(5, 'mm'))

ggsave(filename = "p-value.pdf",
       height   = 45,
       width    = 80, 
       units    = "mm")