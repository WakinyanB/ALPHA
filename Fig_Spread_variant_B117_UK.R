# Phenotypic evolution of SARS-CoV-2: a statistical inference approach

# -> Generates Fig. S1. Epidemiological and genetic data from the COVID-19 outbreak in England
#                   between September 2020 and January 2021 for the resident strain of SARS-CoV-2
#                   and for the Alpha variant.

rm(list=ls())

library(readODS)
library(lubridate)
library(tidyverse)
library(ggpubr)
library(plyr)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(boot)
library(scales)
library(RColorBrewer)

Sys.setlocale("LC_TIME", "English")

# Spread of the variant B.1.1.7 in the UK

## National scale
data <- read_ods("Data/Variant_of_Concern_VOC_202012_01_Technical_Briefing_5_Data_England.ods",
                 sheet = 5, col_names = TRUE)
data$n_Taqpath <- read_ods("Data/Variant_of_Concern_VOC_202012_01_Technical_Briefing_5_Data_England.ods",
                           sheet = 3, col_names = TRUE)[-1,3] # Number of TaqPath tests
data$week <- dmy(data$week)

## Regional scale
data_region <- read_ods("Data/Variant_of_Concern_VOC_202012_01_Technical_Briefing_5_Data_England.ods",
                        sheet = 6, col_names = TRUE)
data_region$Region <- as.factor(data_region$Region)
data_region$week <- dmy(data_region$week)

# Stringency Index (control measures to curb the epidemic of COVID-19)
Stringency_Index_global <- read.csv("Data/covid-stringency-index.csv")
Stringency_Index_global <- Stringency_Index_global[Stringency_Index_global$Entity == "United Kingdom",-2]
Stringency_Index_global$Day <- ymd(Stringency_Index_global$Day)
Stringency_Index <- Stringency_Index_global[Stringency_Index_global$Day >= min(data$week)-7 & 
                                              Stringency_Index_global$Day <= (max(data$week)+6),] # Stringency Index for the studied time period
SI_mean <- data.frame("Week" = c(min(data$week)-7, data$week),
                      "strigency_index_mean" = seq(1,dim(Stringency_Index)[1]-6,7) %>%
                        sapply(function(x){Stringency_Index[x:(x+6),3] %>% mean})) # Weekly averaged

dim <- dim(data)
dim_region <- dim(data_region)
Regions <- unique(data_region$Region) %>% as.character
N_regions <- length(Regions)
threshold.1 = 0.05
threshold.2 = 0.1

plot_tab <- data.frame("Week" = data$week %>% rep(3),
                       "Numbers" = c(data$`n_Confirmed S-gene`,
                                     data$`n_Confirmed SGTF`,data$n_Taqpath),
                       "Percentage" = c(data$`percent_Confirmed S-gene`/100,
                                        data$`percent_Confirmed SGTF`/100, rep(NA, dim[1])),
                       "Legend" = (1:3) %>% rep(each=dim[1]) %>%
                         factor(labels=c("Resident strain", "Alpha variant", "Total tests")))

Fig_tests <- ggplot(data = plot_tab, aes(x = Week, y = Numbers, color = Legend)) +
  geom_line(cex = 0.6) +
  geom_point(cex = 1.4, pch = 19) +
  labs(y = "Confirmed numbers", color = element_blank()) +
  scale_x_date(date_labels = "%b %Y")+
  scale_color_manual(values = c("#fa810f", "#c21919", "black")) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=11, hjust = 0),
        axis.text.y = element_text(size=11, angle=45),
        legend.position = "top")

Fig_tests2 <- ggplot(data = plot_tab[-which(plot_tab$Legend == "Total tests"),],
                     aes(x = Week, y = Numbers, color = Legend)) +
  geom_line(cex = 0.6) +
  geom_point(cex = 1.4, pch = 19) +
  labs(y = "Confirmed numbers", color = element_blank()) +
  scale_x_date(date_labels = "%b %Y")+
  scale_color_manual(values = c("#fa810f", "#c21919")) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=11, hjust = 0),
        axis.text.y = element_text(size=11, angle=45),
        legend.position = "top")

data_qm <- data[, c(1,5)]
data_qm$`percent_Confirmed SGTF` <- data_qm$`percent_Confirmed SGTF`/100
colnames(data_qm)[2] <- "qm"
data_qm$logit_qm <- data_qm$qm %>% logit

Fig_percent <- ggplot(data = subset(plot_tab, Legend != "Total tests"),
                      aes(x = Week, y = Percentage, fill = Legend)) +
  geom_col(cex = 0.5, width = 5.5, col = "black") +
  labs(fill = element_blank()) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),
                     breaks = seq(0,1,0.1), expand = c(0,0)) +
  scale_x_date(date_labels = "%b %Y", expand = c(0,0))+
  scale_fill_manual(values = c("#fa810f", "#c21919", "black")) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        legend.position = "top",
        axis.line.x.bottom=element_line(size=1),
        axis.line.y.left=element_line(size=1)) +
  geom_hline(yintercept = threshold.1, lwd = 0.6)

plot_tab_region <- data.frame("Region" = data_region$Region %>% rep(2),
                              "Week" = data_region$week %>% rep(2),
                              "Numbers" = c(data_region$`n_Confirmed S-gene`,
                                            data_region$`n_Confirmed SGTF`),
                              "Percentage" = c(data_region$`percent_Confirmed S-gene`/100,
                                               data_region$`percent_Confirmed SGTF`/100),
                              "Legend" = 1:2 %>% rep(each=dim_region[1]) %>%
                                factor(labels=c("Resident strain", "Alpha variant")))

Fig_tests_region <- ggplot(data = plot_tab_region,
                           aes(x = Week, y = Numbers, color = Legend)) +
  facet_wrap(~ Region, ncol=3) +
  geom_line(cex = 0.4) +
  geom_point(cex = 0.6, pch = 19) +
  theme_bw() +
  labs(y = "Cases tested positive by TaqPath laboratories", color = element_blank()) +
  scale_color_manual(values = c("#fa810f", "#c21919", "black")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=9.5),
        axis.text.y = element_text(size=9.5, angle=45),
        legend.position = "top")

data_qm_region <- data_region[,c(1,2,6)]
data_qm_region$`percent_Confirmed SGTF` <- data_qm_region$`percent_Confirmed SGTF`/100
colnames(data_qm_region)[3] <- "qm"
data_qm_region$logit_qm <- data_qm_region$qm %>% logit

Fig_percent_region <- ggplot(data = plot_tab_region, aes(x = Week, y = Percentage, fill = Legend)) +
  facet_wrap(~ Region, ncol=3) +
  geom_col(col = 'black', lwd = 0.01) +
  theme_classic() +
  labs(fill = element_blank()) +
  # scale_x_date(expand=c(0,0))+
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = c("#fa810f", "#c21919", "black")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=9.5),
        axis.text.y = element_text(size=9.5),
        legend.position = "top",
        strip.background = element_rect(fill="grey85")) +
  geom_hline(yintercept = threshold.2, col ='white', lwd = 0.7)

#####################################################################################
# Spread_variant_B117_UK
# SVG (Width = 1100, Height = 700)
# pdf 11 x 8

plot_grid(Fig_tests, Fig_percent, Fig_tests_region + theme(legend.position = 'none'),
          Fig_percent_region + theme(legend.position = 'none'),
          labels = c("A)", "B)", "C)", "D)"), rel_heights = c(0.47,0.53))

#####################################################################################