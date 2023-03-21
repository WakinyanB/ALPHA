# Phenotypic evolution of SARS-CoV-2: a statistical inference approach

# -> Generates Fig. 1. The two consecutive phases of the spread of the Alpha variant

rm(list=ls())

library(tidyverse)
library(plyr)
library(ggplot2)
library(readODS)
library(lubridate)
library(gridExtra)
library(cowplot)
library(scales)
library(knitr)
library(boot)
library(pBrackets)

Sys.setlocale("LC_TIME", "English")

# at national scale
data <- read_ods("Data/Variant_of_Concern_VOC_202012_01_Technical_Briefing_5_Data_England.ods",
                 sheet = 5, col_names = TRUE)
data$n_Taqpath <- read_ods("Data/Variant_of_Concern_VOC_202012_01_Technical_Briefing_5_Data_England.ods",
                           sheet = 3, col_names = TRUE)[-1,3] # Nombre de tests TaqPath
# at regional scale
data_region <- read_ods("Data/Variant_of_Concern_VOC_202012_01_Technical_Briefing_5_Data_England.ods",
                        sheet = 6, col_names = TRUE)

# (reflètent les mesures de contrôles misent en place)
Stringency_Index_global <- read.csv("Data/covid-stringency-index.csv")
Stringency_Index_global <- Stringency_Index_global[Stringency_Index_global$Entity == "United Kingdom",-2]

# Regions : class factor
data_region$Region <- as.factor(data_region$Region)

# Date : class dates (package "lubridate")
data$week <- dmy(data$week) + 3
data_region$week <- dmy(data_region$week) + 3
Stringency_Index_global$Day <- ymd(Stringency_Index_global$Day)

# We keep only values of Stringency Index for the studied dates
Stringency_Index <- Stringency_Index_global[Stringency_Index_global$Day >= (min(data$week)-7) &
                                              Stringency_Index_global$Day <= (max(data$week)+6),]

dim <- dim(data)

threshold <- 5 # %
# threshold of frequency of the variant among positive cases below which genetic diversity is neglected
t_threshold <- data$week[which(data$`percent_Confirmed SGTF` > threshold)[1]]

plot_tab_S_gene <- Stringency_Index[,-1] %>% cbind("cases" = NA)
plot_tab_S_gene$cases[match(data$week, Stringency_Index$Day)] <- data$`n_Confirmed S-gene`
plot_tab_S_gene$Strain <- factor("Resident strain")

plot_tab_SGTF <- Stringency_Index[,-1] %>% cbind("cases" = NA)
plot_tab_SGTF$cases[match(data$week, Stringency_Index$Day)] <- data$`n_Confirmed SGTF`
plot_tab_SGTF$Strain <- factor("Alpha variant")

plot_tab <- rbind(plot_tab_S_gene, plot_tab_SGTF)

x_limits <- which(!is.na(plot_tab$cases))
x_limits <- plot_tab$Day[c(x_limits[1], x_limits[length(x_limits)])]

y_limits <- c(0,130000)

Fig <- ggplot(data = plot_tab, aes(x = Day, y = cases)) +
  geom_rect(aes(xmin = Day, xmax = Day+1, ymin = 0, ymax = Inf, fill = stringency_index)) +
  geom_hline(yintercept = plot_tab_SGTF[plot_tab_SGTF$Day == t_threshold, 3],
             linetype = 'longdash') +
  geom_vline(xintercept = t_threshold, linetype = 'longdash') +
  geom_line(data = na.omit(plot_tab), aes(group = Strain, col = Strain), cex = 1.3) +
  labs(x = "\nEngland (2020-2021)", fill = "Stringency Index",
       y = "Weekly number of COVID-19 cases\n (by TaqPath laboratories)",
       col = "Viral Strain") +
  scale_color_manual(values = c("#fa810f", "#c21919")) +
  scale_fill_gradient(low = 'grey95', high = 'grey40') +
  theme_classic() +
  scale_x_date(expand = c(0,0), limits = x_limits) +
  scale_y_continuous(breaks = seq(y_limits[1],y_limits[2],25000),
                     limits = y_limits, expand = c(0,0)) +
  theme(axis.text.x = element_text(size=11), axis.text.y = element_text(size=11, angle=45),
        plot.margin = unit(c(1,1,1,1), "lines"), # add space outside the plot with plot.margin
        legend.position = "right") +
  coord_cartesian(clip = "off")

# 8 x 4.5 (pdf)
plot_grid(plot_grid(ggdraw() + draw_label("Phase 1", fontface = 'bold',
                                          x = 0.58, hjust = 0.5, size = 12),
                    ggdraw() + draw_label("Phase 2", fontface = 'bold',
                                          x = 0.25, hjust = 0.5, size = 12), ncol = 2),
          Fig, ncol = 1, rel_heights = c(0.125, 0.875))
grid.brackets(x1=105, y1=66, x2=337, y2=66, lwd = 2)
grid.brackets(x1=350, y1=65, x2=609, y2=66, lwd = 2)