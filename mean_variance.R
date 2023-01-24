library(rstan)
library(hBayesDM)
library(bayestestR)
library(bayesplot)
library(dplyr)
library(lmerTest)
library(lme4)
library(jtools)
library(stargazer)

setwd("/Users/kseniapanidi/Desktop/KPanidi_remote/PPC/")
data <- read.csv2("ALLdata_v3.csv",header = T, sep = ",", dec=",")
N = unique(data[,'subject_ID'])
exclude = c("S4")
N = N[-which(N%in%exclude)]
myfun <- function(subj, i){
  s = sum(data[which(data$subject_ID==subj & data$session==i),'rt'])
}

df <- data.frame(subj=rep(unique(data$subject_ID),each=3), session=rep(c(0,1,2), length(unique(data$subject_ID))))
df1 <- data.frame()
df1 <- mapply(FUN=myfun, subj=df$subj, i=df$session)
df1 <- df1/1000
df1 <- df1/60
df <- cbind(df,df1)
mean(df[which(df$subj%in%N),3])
sd(df[which(df$subj%in%N),3])
df <- mutate(df, total_rt = sum(data[which(subject_ID==df$subj & session==df$session),'rt']))
############## add personality scales ##########

scales <- read.csv('ankety_PPC.csv')
ss_key <- c(1,1,2,2,2,2,2,2,1,2,1,2,2,2,1,2)
myf <- function(subj) { 
  s = length(which(scales[which(scales$Subject_ID==subj),which( colnames(scales)=="S1" ):which( colnames(scales)=="S16" )]==ss_key))
  return(s)}
scales$SS = mapply(myf, scales$Subject_ID)
data$SS = mapply(myf, data$subject_ID)
scales$Drive = scales$B3+scales$B9+scales$B12+scales$B21
scales$Fun = scales$B5+scales$B10+scales$B15+scales$B20
scales$Reward = scales$B4+scales$B7+scales$B14+scales$B18+scales$B23
scales$BIS = scales$B2+scales$B8+scales$B13+scales$B16+scales$B19+scales$B22+scales$B24

myf1 <- function(subj) {
  s1 = scales[which(scales$Subject_ID==subj),'Drive']
  s2 = scales[which(scales$Subject_ID==subj),'Fun']
  s3 = scales[which(scales$Subject_ID==subj),'Reward']
  s4 = scales[which(scales$Subject_ID==subj),'BIS']
  s5 = scales[which(scales$Subject_ID==subj),'Discomf1']
  s6 = scales[which(scales$Subject_ID==subj),'Discomf2']
  s7 = scales[which(scales$Subject_ID==subj),'Discomf3']
  if (scales[which(scales$Subject_ID==subj),'Gender']=='f')
    {s8=1}
  else {s8=0}
  
  return(cbind(s1,s2,s3,s4,s5,s6,s7,s8))
}

suppl= mapply(myf1, data$subject_ID)
rownames(suppl)=c('Drive', 'Fun', 'Reward','BIS', 'Discomf1', 'Discomf2', 'Discomf3', 'Gender')
data <- cbind(data, t(suppl))
session_order <- read.csv2("PPC_session_order.csv", sep=',')
names(session_order)=c("subject_ID", "session1", "session2", "session3")
data <- mutate(data, session_temp = ifelse(session==1, "A", ifelse(session==2,"B", "C")))
data$order = NA
for (i in 1: nrow(data)){
  s = data[i,'session_temp']
  data[i,'order'] = which(session_order[which(session_order[,1]==data[i,'subject_ID']),]==s)-1
}
data <- mutate(data, trial = rep(1:85, 3*length(unique(data[,'subject_ID']))))

data <- mutate(data, discomf = ifelse( order==1, Discomf1, ifelse(order==2, Discomf2, Discomf3)))

data <- mutate(data, EV_A = outcome1*prob1 + outcome2*prob2)
data <- mutate(data, EV_B = outcome3*prob1 + outcome4*prob2)
data <- mutate(data, sd_A = sqrt(((outcome1-EV_A)^2)*prob1 + ((outcome2-EV_A)^2)*prob2))
data <- mutate(data, sd_B = sqrt(((outcome3-EV_B)^2)*prob1 + ((outcome4-EV_B)^2)*prob2))

data <- mutate(data, diff_m = ifelse( (sd_A>sd_B), EV_A-EV_B, EV_B-EV_A))
data <- mutate(data, diff_sd = ifelse( (sd_A>sd_B), sd_A-sd_B, sd_B-sd_A))

data <- mutate(data, choice_sd = ifelse(sd_A>sd_B, 1, 2) ) 
data <- mutate(data, choice_is_sd = ifelse(choice==choice_sd, 1, 0))

data <- mutate(data, risky = ifelse((sd_A>sd_B & choice==1) | (sd_A<sd_B & choice==2), 1, 0))

attach(data)


####### Mean-variance model #########

library(lmerTest)
library(lme4)

m0  <- lme4::glmer(choice_is_sd  ~ as.factor(session)+I(diff_m/diff_sd)+diff_m+diff_sd+(1|subject_ID), family = 'binomial', data=data[which(prob2!=0 & prob1!=0 & subject_ID%in%N),])
m01 <- lme4::glmer(choice_is_sd  ~ as.factor(session)+I(diff_m/diff_sd)+ diff_m+diff_sd+ trial+order+discomf+(1|subject_ID), family = 'binomial', data=data[which(prob2!=0 & prob1!=0 & subject_ID%in%N),])
m02 <- lme4::glmer(choice_is_sd  ~ as.factor(session)*I(diff_m/diff_sd)+ diff_m+diff_sd+ trial+order+discomf+(1|subject_ID), family = 'binomial', data=data[which(prob2!=0 & prob1!=0 & subject_ID%in%N),])
m03 <- lme4::glmer(choice_is_sd  ~ as.factor(session)*I(diff_m/diff_sd)+ diff_m+diff_sd + trial+order + as.factor(session)*diff_m+as.factor(session)*diff_sd+discomf+(1|subject_ID), family = 'binomial', data=data[which(prob2!=0 & prob1!=0 & subject_ID%in%N),])

stargazer(m0, m01, m02,m03, type="html",
          out="PPC_reg_tables1.doc", intercept.bottom = F, intercept.top = T, report = "vc*"  , digits=4, 
          star.char = c("*", "**", "***"), 
          star.cutoffs = c(0.05, 0.01, 0.001),
          notes = c("* p<0.05", "** p<0.01", "*** p<0.001"), notes.append = F)

ss <- getME(m0,c("theta","fixef"))
m0 <- update(m0,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))




####==================== Bayes tables and plots ======================


#fit <- readRDS("fit_sigma1_v3_16_28_3_6.RDS")
fit <- readRDS("fit_narrow1.RDS")
parvals <- rstan::extract(fit, permuted=TRUE)
varnames = c('mu_delta1_rho', 'mu_delta1_prw', 'mu_delta1_tau', 'mu_delta2_rho', 'mu_delta2_prw', 'mu_delta2_tau')
#prior_intervals = data.frame(low = c(-0.5, -0.5, -5, -0.5, -0.5, -5 ), high=c(0.5, 0.5, 5, 0.5, 0.5, 5))
prior_intervals = data.frame(low = c(-1, -1, -5, -1, -1, -5 ), high=c(1, 1, 5, 1, 1, 5))

rownames(prior_intervals)=varnames

extract_ci <- function(varname){
  hdi = ci(parvals[[varname]], method="HDI")
  hdi = paste0("[",round(hdi[['CI_low']], digits=2) , '; ', round(hdi[['CI_high']], digits=2), "]")
  return(hdi)
}

extract_BF <- function(varname){
  dx <- density(parvals[[varname]], n=1200000)
  xnew <- 0
  new_dens_at_zero <- approx(dx$x,dx$y,xout=xnew) 
  old_dens_at_zero <- 1/(prior_intervals[varname, 'high']-prior_intervals[varname, 'low'])
  BF = old_dens_at_zero/new_dens_at_zero$y
  return(BF)
}


means = sapply(parvals[varnames], FUN=mean)
HDI = sapply(varnames, FUN=extract_ci)
BF = sapply(varnames, FUN=extract_BF)
df<- data.frame(varnames, round(means,3), HDI, round(BF,2))
rownames(df)=NULL
names(df)=c('Parameter', 'Mean', '95% HDI', 'BF')

library(tidyverse) #needed for str_replace() function
df1 <- rbind(c('TMS effects for right PPC', '','',''), df[1:3,], c('TMS effects for left PPC', '','',''), df[4:6,])

output <- stargazer(df1, type='html',summary=FALSE, style='jpam', out='effects_narrow.doc', digits=3)

output1 <- str_replace(output, "mu_delta1_rho", 'Δ risk aversion ( <i>&mu;</i><sub>Δr<sup><i>right</i></sup></sub>)')
output1 <- str_replace(output1, "mu_delta1_prw", 'Δ prob.weighting ( <i>&mu;</i><sub> Δ &gamma;<sup><i>right</i></sup></sub>)')
output1 <- str_replace(output1, "mu_delta1_tau", 'Δ consistency ( <i>&mu;</i><sub> Δ &tau;<sup><i>right</i></sup></sub>)')
output1 <- str_replace(output1, "mu_delta2_rho", 'Δ risk aversion ( <i>&mu;</i><sub>Δr<sup><i>left</i></sup></sub>)')
output1 <- str_replace(output1, "mu_delta2_prw", 'Δ prob.weighting ( <i>&mu;</i><sub> Δ &gamma;<sup><i>left</i></sup></sub>)')
output1 <- str_replace(output1, "mu_delta2_tau", 'Δ consistency ( <i>&mu;</i><sub> Δ &tau;<sup><i>left</i></sup></sub>)')

write_lines(output1, "effects_narrow.doc")


plotHDI(parvals[['mu_delta1_rho']], credMass=0.95, binSize = 50)
plotHDI(parvals[['mu_delta2_rho']], credMass=0.95, binSize = 50)
plotHDI(parvals[['mu_delta1_prw']], credMass=0.95, binSize = 50)
plotHDI(parvals[['mu_delta2_prw']], credMass=0.95, binSize = 50)
plotHDI(parvals[['mu_delta1_tau']], credMass=0.95, binSize = 50)
plotHDI(parvals[['mu_delta2_tau']], credMass=0.95, binSize = 50)

dx <- density(parvals[['mu_delta2_tau']], n=1200000)
xnew <- 0
approx(dx$x,dx$y,xout=xnew) 


plotHDI(parvals[['mu_delta1_tau']], credMass=0.89, binSize = 50)
plotHDI(parvals[['mu_delta2_tau']], credMass=0.89, binSize = 50)

#plotHDI(parvals[['mu_delta2_prw']][which(parvals[['mu_delta2_rho']]>0)], credMass=0.89, binSize=50)



posterior = data.frame(mu_rho=parvals[['mu_rho']], 
                       mu_prw=parvals[['mu_prw']],
                       mu_tau = parvals[['mu_tau']])
names(posterior)=c("Baseline RA", "Baseline Gamma", "Baseline Tau")


posteriors = data.frame(mu_rho=parvals[['mu_rho']], 
                        mu_prw=parvals[['mu_prw']], 
                        mu_tau=parvals[['mu_tau']], 
                        mu_delta1_rho=parvals[['mu_delta1_rho']],
                        mu_delta1_prw=parvals[['mu_delta1_prw']],
                        mu_delta1_tau=parvals[['mu_delta1_tau']],
                        mu_delta2_rho=parvals[['mu_delta2_rho']],
                        mu_delta2_prw=parvals[['mu_delta2_prw']],
                        mu_delta2_tau=parvals[['mu_delta2_tau']])
post = as.list(posteriors)
#labels= c("Risk aversion", "Probability weighting", "Consistency")
#  \u0394 is a Unicode for a greek letter Delta
library(stringr) #to wrap long text in labels
#labels= c("Risk aversion","Prob. weighting", "Consistency", "\u0394 risk aversion (right vs. sham)", "\u0394 prob. weighting (right vs. sham)", "\u0394 consistency (right vs.sham)","\u0394 risk aversion (left vs. sham)", "\u0394 prob. weighting (left vs. sham)", "\u0394 consistency (left vs.sham)" )
labels= c("Risk aversion","Prob. weighting", "Consistency", "\u0394 risk aversion", "\u0394 prob. weighting", "\u0394 consistency","\u0394 risk aversion", "\u0394 prob. weighting", "\u0394 consistency" )
plot_labels = data.frame( vars = names(post)[1:9], labels=str_wrap(labels, width=20))

####Draw MCMC plots with mcmc_dens()
#draw_HDI <- function(x){
  #par_plot <- mcmc_hist(as.data.frame(post[[x]]), binwidth=0.05, freq=TRUE)
#  par_plot <- mcmc_dens(as.data.frame(post[[x]]), alpha=0.5)
#  ci_int_89 <- bayestestR::ci(post[[x]], ci=0.89, method='HDI')
#  ci_int_95 <- bayestestR::ci(post[[x]], ci=0.95, method='HDI')
#  xlab0 <- plot_labels[x,2]
  #par_plot <-par_plot+geom_segment(aes(x=ci_int[,2],xend=ci_int[,3],y=0,yend=0), colour="red", size=2)+yaxis_text(on=TRUE)+ylab("frequency") + xlab(xlab0)
  #par_plot <-par_plot+geom_segment(aes(x=ci_int_89[,2],xend=ci_int_89[,3],y=0,yend=0), colour="black", size=4)+geom_segment(aes(x=ci_int_95[,2],xend=ci_int_95[,3],y=0,yend=0), colour="black", size=1)+yaxis_text(on=TRUE)+ylab("density") + xlab(xlab0)
#  par_plot <-par_plot+geom_segment(aes(x=ci_int_89[,2],xend=ci_int_89[,3],y=0,yend=0),size=4)+geom_segment(aes(x=ci_int_95[,2],xend=ci_int_95[,3],y=0,yend=0), colour="black", size=1)+yaxis_text(on=TRUE)+ylab("density") + xlab(xlab0)
#  return(par_plot)
#}

####Draw MCMC plots manually
draw_HDI <- function(i){
  data <- data.frame(mcmc=posteriors[,i])
  cutoff1 <- bayestestR::ci(data$mcmc, ci=0.95, method='HDI')[,2]
  cutoff2 <- bayestestR::ci(data$mcmc, ci=0.95, method='HDI')[,3]
  xlab0 <- plot_labels[i,2]
  
  d <- density(data$mcmc) 
  d = data.frame(mcmc=d[['x']], dens=d[['y']])
  #hist.y=mutate(hist.y, area = x >= cutoff1 & x<=cutoff2)
  data1 <- subset(d, d$mcmc>cutoff1 & d$mcmc<cutoff2)
  
  p<-ggplot(data) + geom_density(aes(x=mcmc), colour='#404B6F',fill='#8AA2F0', alpha=0.5)+
    geom_area(data=data1, aes(x=mcmc,y=dens), colour='#404B6F',fill='#788DCF',alpha=0.8)+
    geom_vline(xintercept =mean(data$mcmc), color='#404B6F', linetype='dashed', size=0.7)+
    yaxis_text(on=TRUE)+ylab("density") + xlab(xlab0)+
    theme_classic()
  return(p)
}



#color_scheme_set("brightblue")
#my_color_scheme <- c("sienna", "sienna1",
#                     "sienna2", "sienna3",
#                     "sienna4", "sienna1")
#color_scheme_set(my_color_scheme)
#color_scheme_set('brewer-Blues') #sets a colorscheme from RColorBrewer
pars = as.list(c(1:9))
posterior_plots <- lapply(pars, FUN=draw_HDI)
posterior_plots[[1]]


library(gridExtra)
library(ggpubr)
library(cowplot)
theme_set(theme_cowplot(font_size=12, font_family = "Helvetica") + 
            theme(text = element_text(colour = "black")))
plots1 <- ggarrange(posterior_plots[[1]],posterior_plots[[2]], posterior_plots[[3]], ncol=3, nrow=1, labels=c("A.1", "A.2", "A.3"), font.label = list(size=9), vjust=1, hjust=0.05)
#plots1 <- annotate_figure(plots1, fig.lab.pos="top.left", fig.lab="Effect of the right PPC stimulation", fig.lab.size = 10)
plots2 <- ggarrange(posterior_plots[[4]],posterior_plots[[5]], posterior_plots[[6]], ncol=3, nrow=1, labels=c("B.1", "B.2", "B.3"), font.label = list(size=9), vjust=1.4, hjust=0.5)
plots3 <- ggarrange(posterior_plots[[7]],posterior_plots[[8]], posterior_plots[[9]], ncol=3, nrow=1, labels=c("C.1", "C.2", "C.3"), font.label = list(size=9), vjust=1.2, hjust=0.5)
plots <- plot_grid(plotlist=list(plots2, plots3), scale=c(0.8,0.8), vjust=1.05, hjust=-0.3, nrow=2, ncol=1, labels=c("TMS effects: right PPC", "TMS effects: left PPC"), label_size=10 )
plots
ggsave("baseline_plots_PPC.tiff", plot = plots1, dpi=300, compression = 'lzw', height=5, width=18,units="cm", scale=0.8)
ggsave("delta_plots_PPC.tiff", plot = plots, dpi=300, compression = 'lzw', height=11, width=18,units="cm")

#ggsave("baseline_plots_PPC.tiff",dpi=300, compression = 'lzw', height=5, width=15,units="cm")
#ggsave("delta_plots_PPC.tiff",dpi=300, compression = 'lzw', height=10, width=15,units="cm")

######========== Posterior predictive check ==========

#parvals <- rstan::extract(fit, permuted=TRUE)
attach(data)
subjList <- unique(data[, "subject_ID"]) 
#exclude = c("S4", "S9")
#N = unique(data[,'subject_ID'])
exclude = c("S4")
N = subjList[-which(subjList%in%exclude)]

pred0=parvals[['y_pred']]
raw_data <- data[which(subject_ID%in%N),]


library(progress) #to add progress bar in the cycle
LB<- c()
HB <- c()
v1 <- sample.int(16000, 4000, replace=FALSE)
n_iter <- length(v1) # Number of iterations of the loop
pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = n_iter,
                       complete = "=",   # Completion bar character
                       incomplete = "-", # Incomplete bar character
                       current = ">",    # Current bar character
                       clear = FALSE,    # If TRUE, clears the bar when finish
                       width = 100)      # Width of the progress bar
mean_ppp <- c()
for (v in v1){
  # curSubj <- subjList[i]
  #tmp <- raw_data [ which (raw_data$subject!=14 & raw_data$subject!=16 ), ]
  pb$tick()
  tmp <- raw_data
  choice0 <- tmp[1:nrow(tmp), "choice"]
  choice<- (choice0-2)*(-1)
  ppp_ij <- c()
  #subjList <- unique(raw_data[, "subject_ID"]) 
  #subjList <- subjList[c(-14,-16)]
  #numSubjs=length(subjList)
  #maxTrials=48
  #Tsubj=as.vector(rep(maxTrials, numSubjs))
  #choice0 <- array(0, c(numSubjs, maxTrials))
  #pred = pred_right[v,,]
  pred=pred0[v,,]
  pred <- as.vector(t(pred))
  for (j in 1:nrow(tmp)){
    if (choice[j]==pred[j]){ppp_ij <- cbind(ppp_ij, 1)  
    } else{ppp_ij <- cbind(ppp_ij, 0)}
  }
  
  mean_ppp <- cbind(mean_ppp, mean(ppp_ij))
}

LB = rbind(LB, HDIofMCMC(mean_ppp)[1])
HB = rbind(HB, HDIofMCMC(mean_ppp)[2])

cbind(LB, HB)
plotHDI(t(mean_ppp))
posterior=data.frame(ppp=t(mean_ppp))
names(posterior)[1]="Mean PPP"


