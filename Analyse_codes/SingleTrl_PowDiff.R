##### run linear mixed model to investigate the differnece between powers at the tagging frequency bands and the noise level

# load in library
install.packages("readxl")
library(readxl)
library(lme4)
library(lmerTest)
library(ggplot2)

# read in dataset
fovea_para <- read_excel("SingleTrl_PowDiff.xlsx",sheet = "fovea_para")
fovea_noise <- read_excel("SingleTrl_PowDiff.xlsx",sheet = "fovea_noise")
para_noise <- read_excel("SingleTrl_PowDiff.xlsx",sheet = "para_noise")
minipow_noise <- read_excel("SingleTrl_PowDiff.xlsx",sheet = "minipow_noise")

## single_trial power differences between fovea & para
model <- lmer(PowDiff ~ 1 + (1 | Subjects), data = fovea_para)
summary(model)
model <- lmer(PowDiff ~ 1 + (1 | Subjects), data = fovea_para)
summary(model)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
#Random effects:
#  Groups   Name        Variance Std.Dev.
#Subjects (Intercept) 0.21373  0.4623  
#Residual             0.04412  0.2100  
#Number of obs: 7930, groups:  Subjects, 29
#Fixed effects:
#  Estimate Std. Error       df t value Pr(>|t|)
#(Intercept) -0.06464    0.08588 27.99589  -0.753    0.458
g = ggplot(fovea_para, aes(x = PowDiff)) +
      geom_histogram(bins = 30, fill ='gray50') +
      facet_wrap(~ Subjects) +
      xlab("Normalized single-trial power differences between fovea and parafovea (fovea-parafovea)") +
      ylab("Histogram (# of counts)")
g+theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
ggsave(file="fovea_para.svg", plot=g, width=10, height=8)

##
model <- lmer(PowDiff ~ 1 + (1 | Subjects), data = fovea_noise)
summary(model)
#Random effects:
#  Groups   Name        Variance Std.Dev.
#Subjects (Intercept) 0.24527  0.4952  
#Residual             0.02186  0.1478  
#Number of obs: 7930, groups:  Subjects, 29
#Fixed effects:
#  Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)  0.71680    0.09198 27.99566   7.793 1.73e-08 ***

model <- lmer(PowDiff ~ 1 + (1 | Subjects), data = para_noise)
summary(model)
#Random effects:
#  Groups   Name        Variance  Std.Dev.
#Subjects (Intercept) 0.0002914 0.01707 
#Residual             0.0082090 0.09060 
#Number of obs: 7930, groups:  Subjects, 29
#Fixed effects:
#  Estimate Std. Error        df t value Pr(>|t|)    
#(Intercept)  0.896275   0.003333 28.218591   268.9   <2e-16 ***

model <- lmer(PowDiff ~ 1 + (1 | Subjects), data = minipow_noise)
summary(model)
#Random effects:
#  Groups   Name        Variance  Std.Dev. 
#Subjects (Intercept) 4.479e-50 2.116e-25
#Residual             1.732e-49 4.162e-25
#Number of obs: 7930, groups:  Subjects, 29
#Fixed effects:
#  Estimate Std. Error        df t value Pr(>|t|)    
#(Intercept) 3.126e-25  3.958e-26 2.797e+01   7.896 1.35e-08 ***













