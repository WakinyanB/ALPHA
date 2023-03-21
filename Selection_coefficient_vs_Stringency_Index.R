# Phenotypic evolution of SARS-CoV-2: a statistical inference approach

# -> Explores the correlation between the selection coefficient of the Alpha variant and the Stringency Index in the UK.
# -> Generates Fig. S2. The selection coefficient of the Alpha variant in England is negatively correlated
#                       with the Stringency Index (fall - winter 2020-2021).

rm(list=ls())

library(tidyverse)
library(plyr)
library(ggplot2)
library(readODS)
library(lubridate)
library(gridExtra)
library(cowplot)
library(deSolve)
library(scales)
library(knitr)
library(boot)
library(lme4)

Sys.setlocale("LC_TIME", "English")

# Data import, processing and graphical visualisation (UK)

## Data import

### Demographic and genetic data

# at national scale
data <- read_ods("Data/Variant_of_Concern_VOC_202012_01_Technical_Briefing_5_Data_England.ods",
                 sheet = 5, col_names = TRUE)
data$n_Taqpath <- read_ods("Data/Variant_of_Concern_VOC_202012_01_Technical_Briefing_5_Data_England.ods",
                           sheet = 3, col_names = TRUE)[-1,3] # Nombre de tests TaqPath
# at regional scale
data_region <- read_ods("Data/Variant_of_Concern_VOC_202012_01_Technical_Briefing_5_Data_England.ods",
                        sheet = 6, col_names = TRUE)

### Stringency Index

# (reflect the amount of NPIs implemented to mitigate the spread of the epidemic)
Stringency_Index_global <- read.csv("Data/covid-stringency-index.csv")
Stringency_Index_global <- Stringency_Index_global[Stringency_Index_global$Entity == "United Kingdom",-2]

## Data process

# Regions : class factor
data_region$Region <- as.factor(data_region$Region)

# Date : class dates (package "lubridate")
data$week <- dmy(data$week)
data_region$week <- dmy(data_region$week)
Stringency_Index_global$Day <- ymd(Stringency_Index_global$Day)

# We keep only values of Stringency Index for the studied dates
Stringency_Index <- Stringency_Index_global[Stringency_Index_global$Day >= (min(data$week)-7) &
                                              Stringency_Index_global$Day <= (max(data$week)+6),]
# Weekly average
StrInd_mean <- data.frame("Week" = c(min(data$week)-7, data$week),
                          "stringency_index_mean" = seq(1,dim(Stringency_Index)[1]-6,7) %>%
                            sapply(function(x){Stringency_Index[x:(x+6),3] %>% mean}))

Fig_stringency_index <- ggplot(data = Stringency_Index, aes(x = Day, stringency_index)) +
  geom_step(cex = 0.9, col = "darkred") +
  ylab("Stringency Index (UK)\n") +
  ylim(c(0,100)) +
  scale_x_date(date_labels = "%b %Y")+
  theme_bw() +
  theme(axis.text.x = element_text(size=11, hjust=0),
        axis.text.y = element_text(size=11),
        axis.title.x = element_blank())

Fig_stringency_index_mean <- ggplot(data = StrInd_mean, aes(x = Week, stringency_index_mean)) +
  geom_step(cex = 1, col ="blue") +
  ylab("Stringency Index (UK, weekly average) \n") +
  ylim(c(0,100)) +
  scale_x_date(date_labels = "%b %Y")+
  theme_bw() +
  theme(axis.text.x = element_text(size=11, hjust=0),
        axis.text.y = element_text(size=11),
        axis.title.x = element_blank())

grid.arrange(Fig_stringency_index, Fig_stringency_index_mean,
             ncol = 1, nrow= 2, bottom = paste("Time series of Stringency Index (top) or weekly average (bottom) in the UK
(from", min(Stringency_Index$Day), "to", max(Stringency_Index$Day), ")"))

# How does the Stringency Index impact the selection coefficient of the variant B.1.1.7 ?

threshold = 0.1

# National scale

dim <- dim(data)

data_fm <- data[, c(1,5)]
data_fm$`percent_Confirmed SGTF` <- data_fm$`percent_Confirmed SGTF`/100
colnames(data_fm)[2] <- "fm"
data_fm$logit_fm <- data_fm$fm %>% logit

data_fm <- data_fm[data_fm$fm > threshold,]
s_weeks <- data_fm$week[-c(1, dim(data_fm)[1])]
s = data.frame("Week" = s_weeks,
               "s" = (data_fm$logit_fm[-c(1,2)]-data_fm$logit_fm[-(dim(data_fm)[1]-c(1,0))])/2,
               "StrInd" = StrInd_mean[match(s_weeks, StrInd_mean$Week),2])

# Regional scale

Regions <- unique(data_region$Region) %>% as.character
N_regions <- length(Regions)
dim_region <- dim(data_region)

data_fm_region <- data_region[,c(1,2,6)]
data_fm_region$`percent_Confirmed SGTF` <- data_fm_region$`percent_Confirmed SGTF`/100
colnames(data_fm_region)[3] <- "fm"
data_fm_region$logit_fm <- data_fm_region$fm %>% logit

Fig_logit_fm_region <- ggplot(data = data_fm_region, aes(x = week, y = logit_fm, col = Region)) +
  geom_line(cex = 0.5) +
  geom_point(pch = 19, cex = 1.5) +
  ylab("logit-frequency of the\n Alpha variant (England)\n") +
  scale_x_date(date_labels = "%b %Y", limits = c(min(StrInd_mean$Week), max(StrInd_mean$Week)))+
  theme_bw() +
  theme(axis.text.x = element_text(size=11, hjust = 0),
        axis.text.y = element_text(size=11),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = logit(threshold), lty = 'dashed', lwd = 0.4)

data_fm_region <- data_fm_region[data_fm_region$fm > threshold,]
s_region <- data_fm_region %>% ddply(~Region, function(x){N_weeks <- dim(x)[1]
data.frame("Week" = x[-c(1, N_weeks),2],
           "s" = (x[-c(1,2),4]-
                    x[-(N_weeks + c(-1,0)),4])/2)})
s_region$StrInd <- StrInd_mean[match(s_region$Week, StrInd_mean$Week),2]

## Graphical vizualisations

plot_s_vs_StrInd <- function(data, ylim = c(0,5)){
  
  # data = data frame with colums :
  #               - "s" = selection coefficient
  #               - "StrInd" = corresponding Stringency Index
  # (FACULTATIVE) - "Region" = names of the regions (whether we want to distinguish them)
  
  colnames(data)[grep("StrInd", colnames(data))] <- "StrInd"
  
  if(!is.null(data$Region)){
    Fig <- ggplot(data = data, aes(x = StrInd, y = s, col = Region)) +
      geom_line(cex = 0.5)
  }else{
    Fig <- ggplot(data = data, aes(x = StrInd, y = s))
  }
  Fig <- Fig +
    geom_point(cex = 1.7, pch = 19) +
    ylim(ylim) +
    xlab(paste("Stringency Index (UK, weekly average)")) +
    ylab("Selection coefficient of the variant s(t)") +
    theme_bw() +
    theme(axis.text.x = element_text(size=11),
          axis.text.y = element_text(size=11))
  
  return(Fig)
}

ymax <- max(s$s, s_region$s)*1.01

s_fct_StrInd <- plot_s_vs_StrInd(data = s[,-4], ylim = c(0,ymax)) # national scale
s_fct_StrInd_region <- plot_s_vs_StrInd(data = s_region[,-5], ylim = c(0,ymax)) # regional scale

# graphical visualization

plot_grid(s_fct_StrInd + geom_smooth(method = "lm"), s_fct_StrInd_region,
          ncol = 2, nrow = 1, rel_widths = c(0.4,0.6), labels=c("National scale", "Regional scale"),
          label_size=12, label_x=0.05, label_y=0.3)

# Selection_coefficient_vs_Stringency_Index
# SVG: 950 x 600
# PDF: 9.5 x 6
plot_grid(Fig_stringency_index + ylim(c(50,100)) + theme(axis.title.x = element_blank()),
          s_fct_StrInd + ylab("Selection coefficient s(t)\nof the Alpha variant (per week)\n"),
          Fig_logit_fm_region + theme(axis.text.x = element_text(size = 10, angle = -30)),
          s_fct_StrInd_region + theme(legend.position = 'none') +
            ylab("Selection coefficient s(t)\nof the Alpha variant (per week)\n"),
          ncol = 2, rel_widths = c(0.6,0.4), labels = paste0(LETTERS[1:4],")"),
          label_x = -0.012, label_y = 1.012)

## Correlations between s(t) and the Stringency Index

cor_national <- cor.test(s$s, s$StrInd)

cor_regions <- ddply(s_region, ~Region, function(x){
  r <- cor.test(x$s, x$StrInd)
  return(c(r$estimate, r$conf.int, r$p.value))
})

colnames(cor_regions)[3:5] <- c("CI95.lower", "CI95.upper", "p.value")

p <- data.frame("signif" = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                "code" = c("***", "**", "*", ".", rep(" ", 2)))

kable(cbind(c("UK (national)", " ", data_region$Region %>% unique %>% as.character),
            c(cor_national$estimate %>% unname %>% round(digits=3), " ",
              round(cor_regions$cor, digits=3)),
            c(paste0("[", round(cor_national$conf.int[1], digits=3), "; ",
                     round(cor_national$conf.int[2], digits=3),"]"), " ",
              paste0("[", round(cor_regions$CI95.lower, digits=3), "; ",
                     round(cor_regions$CI95.upper, digits=3), "]")),
            c(signif(cor_national$p.value, digits=3), " ",
              signif(cor_regions$p.value, digits=3)),
            c(p$code[which(cor_national$p.value-p$signif <= 0)[1]-1], " ",
              sapply(cor_regions$p.value, function(x){p$code[which(x-p$signif<=0)[1]-1]}))),
      col.names = c("Geographical location", "Correlation", "95% CI", "p-value", ""))

## Linear model (regional level)

reg1 <- lm(s ~ Region, data = s_region)
reg2 <- lm(s ~ Region + StrInd, data = s_region)
reg3 <- lm(s ~ Region*StrInd, data = s_region)

anova(reg1,reg2)
anova(reg2,reg3)

par(mfrow=c(2,2))
plot(reg3)
par(mfrow=c(1,1))

reg3 %>% residuals %>% shapiro.test
car::durbinWatsonTest(reg3)

summary(reg3)