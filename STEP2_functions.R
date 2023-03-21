#########################################################################
#######           Phenotypic evolution of SARS-CoV-2:             #######
#######             a statistical inference approach              #######
#########################################################################
#######            STEP 2 - R Script with functions               #######
#########################################################################

# Packages
# (may first require installations: install.packages())

# Packages (may first require installations: install.packages())
library(readODS) # to read ODS files
library(lubridate) # to convert into dates
## Data manipulation
library(tidyverse)
library(plyr)
## Plots
library(ggplot2)
library(scales)
library(gridExtra)
library(cowplot)
library(ggpubr)
## Colors and palettes
library(RColorBrewer)
library(lattice)
library(viridis)
##
library(deSolve) # to (numerically) integrate a system of ODEs
library(boot) # for functions 'logit' and 'inv.logit'
library(knitr) # to draw nice tables
library(nlme) # to estimate with correlated residuals


# Estimations of phenotypic differences in transmission (Delta beta) and in recovory (Delta gamma)
# between the Alpha variant and the WT strain (Mixed-Effects Model (MEM) or ANCOVA)

Estim_Phase2_SEIR <- function(data, k_a_df, kappa = 0.2, beta_w = 0.25, gamma_w = 0.1, SN = 0.75,
                              Stringency_Index, subdivisions = 200, method = "lmer"){
  # data = data frame with a column "Region", a column "t" (time)
  #                    and a column "logit_fm" (logit of the frequency of the variant)
  # k_a_df = data frame with previous estimates for parameters (respectively) k and a (cf. phase 1)
  # kappa = transition rate from the exposed to the infectious state
  # Parameter values of the WT: beta_w = transmission rate ; gamma_w = recovery rate
  # SN = Fraction of the population that is susceptible (S(t)/N), approx as a constant through time
  # Stringency_Index = values of Stringency Index ; note that element at position i must be 
  #                    the Stringency Index for time points i-1 <= t < i (e.g. 0 <= t < 1 --> index 1)
  # subdivisions = number of subdivisions used by 'integrate'
  # method = 'lmer' (mixed-effects model with variable 'Region' treated as a random effect)
  #       or 'lm'   (ANCOVA) 
  #                                   
  Linear_regressions <- list()
  
  # f1 and f2, functions to integrate:
  f1 <- function(t, kappa, beta, gamma, SN, k, a, stringency_index){ # (in front of Delta beta)
    c <- k*(stringency_index[trunc(t)+1]/100)^a
    return((1-c)/sqrt((kappa-gamma)^2 + 4*kappa*(1-c)*beta*SN))
  }
  f2 <- function(t, kappa, beta, gamma, SN, k, a, stringency_index){ # (in front of Delta gamma)
    c <- k*(stringency_index[trunc(t)+1]/100)^a
    return(1/sqrt((kappa-gamma)^2 + 4*kappa*(1-c)*beta*SN))
  }
  
  for(i in 1:nrow(k_a_df)){
    
    k <- k_a_df$k[i]
    a <- k_a_df$a[i]
    
    data_lm <- data %>% ddply(~Region, function(tab){
      t0 <- tab$t[1]
      tab$t %>% sapply(FUN = function(t,t0,kappa,beta,gamma,SN,k,a,stringency_index,subdivisions){
        
        coef_Delta_beta <- kappa*SN*
          integrate(f1, lower=t0, upper=t, subdivisions = subdivisions,
                    kappa=kappa, beta=beta, gamma=gamma, SN=SN, k=k, a=a,
                    stringency_index = stringency_index)$value
        
        coef_Delta_gamma <- -0.5*((kappa - gamma)*
                                    integrate(f2, lower=t0, upper=t, subdivisions = subdivisions,
                                              kappa=kappa, beta=beta, gamma=gamma, SN=SN, k=k, a=a,
                                              stringency_index = stringency_index)$value + t-t0)
        
        return(c("t" = t, "Delta_beta" = coef_Delta_beta, "Delta_gamma" = coef_Delta_gamma))
        
      }, t0 = t0, kappa = kappa, beta = beta_w, gamma = gamma_w, SN = SN, k = k, a = a,
      stringency_index = Stringency_Index, subdivisions = subdivisions) %>% data.frame %>% t
    }) %>% cbind("logit_fm" = data$logit_fm)
    
    if(method == "lm"){
      
      lm.model <- lm(logit_fm ~ Region + Delta_beta + Delta_gamma, data = data_lm)
      
      lm.fit <- data.frame("Region" = data_lm$Region, "Date" = data_lm$t+min(data_region$week),
                           "obs_logit_fm" = data_lm$logit_fm, "fitted_values" = lm.model$fitted.values)
      
      lm.estimates <- cbind("estimate" = lm.model$coefficients, confint(lm.model))
      
      Linear_regressions[[i]] <- list("model" = lm.model, "data_model" = data_lm,
                                      "fit" = lm.fit, "estimates" = lm.estimates)
      
    }else if(method == "lmer"){
      
      mem.model <- lme4::lmer(logit_fm ~ Delta_beta + Delta_gamma + (1|Region), data = data_lm)
      
      mem.fit <- data.frame("Region" = data_lm$Region, "Date" = data_lm$t+min(data_region$week),
                            "obs_logit_fm" = data_lm$logit_fm, "fitted_values" = fitted(mem.model))
      
      mem.estimates <- cbind("estimate" = c(unique(coef(mem.model)$Region$Delta_beta),
                                            unique(coef(mem.model)$Region$Delta_gamma)),
                             confint(mem.model)[4:5,])
      
      Linear_regressions[[i]] <- list("model"=mem.model, "data_model" = data_lm,
                                      "fit"=mem.fit, "estimates"=mem.estimates)
    }else{
      stop("Argument 'method' must be either 'lm' (for ANCOVA) or 'lmer' (for mixed-effects model)")
    }
  }
  return(Linear_regressions)
}


# Density plot for a given parameters

plot_density <- function(name, data, col = "red", lwd = 1, alpha = 0.05){
  # name: name of the parameter (-> name of the x-axis)
  # data: vector containing all the values of the parameter
  # col: color of the curb (red by default)
  # alpha: order for defining quantiles (default: 0.05 / 5%)
  
  if(alpha < 0 | alpha >= 1){stop("ERROR: 'alpha' must belong to the interval [0 ; 1 [")}
  data <- data[!is.na(data)] %>% as.data.frame
  
  quantile <- data %>% apply(2, function(x){quantile(x, probs = c(alpha/2, 1-alpha/2))})
  
  return(ggplot(data, aes(x = data[,1])) +
           geom_density(cex = lwd, col = col, fill = col, alpha = 0.1) +
           xlim(c(quantile[1], quantile[2])) +
           xlab(name) +
           theme_bw() +
           ylab(paste0("Density (", 100*(1-alpha), "%) \n")) +
           theme(axis.text.x = element_text(size=11), axis.text.y = element_text(size=11),
                 axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 12)))
}

# Plot estimations and 95% CI of Delta gamma against Delta beta

plot_coord <- function(coord, col = 'RSS', density = FALSE, lwd = 0.4, point_size = 1, digits = 2){
  
  # coord = data.frame with coordinates for estimates and confidence interval of parameters
  #         Delta beta (columns 'Delta_beta', 'Delta_beta_2.5', 'Delta_beta_97.5')
  #     and Delta gamma (columns 'Delta_gamma', 'Delta_gamma_2.5', 'Delta_gamma_97.5')
  #         optional: column 'RSS' (Residual Sum of Squares)
  
  # density = logical indicating whether density plots for Delta beta and Delta gamma should be added
  # col = color: either a color code or (default) 'RSS'
  #                                          (need a column 'RSS' (Residual Sum of Squares) in coord)
  # lwd = line width
  # point_size = point size in ggplot (argument 'cex')
  # accuracy = accuracy for numbers on x- and y-axis
  
  if(col == 'RSS'){
    Fig <- ggplot(data = coord, aes(col = RSS)) +
      geom_segment(aes(x = Delta_beta_2.5, xend = Delta_beta_97.5,
                       y = Delta_gamma, yend = Delta_gamma), size = lwd) +
      geom_segment(aes(x = Delta_beta, xend = Delta_beta,
                       y = Delta_gamma_2.5, yend = Delta_gamma_97.5), size = lwd) +
      geom_point(aes(x = Delta_beta, y = Delta_gamma), cex = point_size) +
      scale_color_gradientn(colors = brewer.pal(9, "RdPu")[-(1:2)] %>% rev) +
      labs(x = expression(paste("Transmission effect  ", Delta, beta)),
           y = expression(paste("Recovery effect  ", Delta, gamma)),
           col = "Residual Sum\n of Squares")
  }else{
    Fig <- ggplot(data = coord)+
      geom_segment(aes(x = Delta_beta_2.5, xend = Delta_beta_97.5,
                       y = Delta_gamma, yend = Delta_gamma), col = col, size = lwd) +
      geom_segment(aes(x = Delta_beta, xend = Delta_beta,
                       y = Delta_gamma_2.5, yend = Delta_gamma_97.5), col = col, size = lwd) +
      geom_point(aes(x = Delta_beta, y = Delta_gamma), col = col, cex = point_size) +
      labs(x = expression(paste("Transmission effect  ", Delta, beta)),
           y = expression(paste("Recovery effect  ", Delta, gamma)))
  }
  Fig <- Fig +
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    scale_x_continuous(position = 'top', labels = function(x){add_plus_sign(x, digits = digits)}) +
    scale_y_continuous(labels = function(x){add_plus_sign(x, digits = digits)}) +
    theme_bw() +
    theme(axis.text.x = element_text(size=11), axis.text.y = element_text(size=11),
          axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))
  
  if(density == TRUE){
    return(
      plot_grid(Fig + theme(legend.position = 'none'),
                plot_density(name=element_blank(), data=coord$Delta_gamma,
                             col = 'black', alpha=0, lwd=0.4) +
                  labs(y = "\n") +
                  xlim(ggplot_build(Fig)$layout$panel_scales_y[[1]]$range$range) +
                  theme_void()+
                  theme(axis.title.x = element_text(size = 13),
                        axis.text.x = element_blank(), axis.text.y = element_blank(),
                        axis.line.x = element_line(color = 'transparent'),
                        axis.line.y = element_line(color = 'transparent')) +
                  coord_flip() + scale_y_continuous(position = 'right'),
                get_legend(Fig),
                plot_density(name=element_blank(), data=coord$Delta_beta,
                             col = 'black', alpha=0, lwd=0.4) +
                  labs(y = "\n") +
                  xlim(ggplot_build(Fig)$layout$panel_scales_x[[1]]$range$range) +
                  scale_y_reverse() +
                  theme_classic() +
                  theme(axis.title.y = element_text(size = 20),
                        axis.text.y = element_blank(), axis.text.x = element_blank(),
                        axis.line.x = element_line(color = 'transparent'),
                        axis.line.y = element_line(color = 'transparent'),
                        axis.ticks = element_blank()),
                NULL, NULL,
                ncol = 3, rel_widths = c(1, 0.25, 0.3), rel_heights = c(1, 0.3))
    )
  }else{
    return(Fig)
  }
}

# Same but for the part with small variations in fixed parameters

plot_coord.var <- function(coord, ref, legend = NULL, colors = NULL, lwd = 0.4, point_size = 1.3, digits = 2){
  
  # coord = data.frame with coordinates for estimates and confidence interval of parameters
  #         Delta beta (columns 'Delta_beta', 'Delta_beta_2.5', 'Delta_beta_97.5')
  #     and Delta gamma (columns 'Delta_gamma', 'Delta_gamma_2.5', 'Delta_gamma_97.5')
  #         optional: column 'RSS' (Residual Sum of Squares)
  # ref = reference values for the parameter that is varied
  # legend = label for the legend (colors)
  # lwd = line width
  # point_size = point size in ggplot (argument 'cex')
  # accuracy = accuracy for numbers on x- and y-axis
  
  i <- which(!colnames(coord) %in% c("Delta_beta", "Delta_gamma", "Delta_beta_2.5",
                                     "Delta_gamma_2.5", "Delta_beta_97.5", "Delta_gamma_97.5"))
  if(is.null(legend)){legend <- colnames(coord)[i]}
  colnames(coord)[i] <- "col"
  
  var <- paste0("(", add_plus_sign(100*(coord$col - ref)/ref, digits = 0), "%)")
  var[var== "( 0%)"] <- ""
  coord$col <- coord$col %>% paste(var) %>% as.factor
  
  Fig <- ggplot(data = coord)+
    geom_segment(aes(x = Delta_beta_2.5, xend = Delta_beta_97.5,
                     y = Delta_gamma, yend = Delta_gamma, col = col), size = lwd) +
    geom_segment(aes(x = Delta_beta, xend = Delta_beta,
                     y = Delta_gamma_2.5, yend = Delta_gamma_97.5, col = col), size = lwd) +
    geom_point(aes(x = Delta_beta, y = Delta_gamma, col = col), cex = point_size) +
    labs(x = expression(paste("Transmission effect  ", Delta, beta)),
         y = expression(paste("Recovery effect  ", Delta, gamma)),
         col = legend) +
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    scale_x_continuous(position = 'top', labels = function(x){add_plus_sign(x, digits = digits)}) +
    scale_y_continuous(labels = function(x){add_plus_sign(x, digits = digits)}) +
    theme_bw() +
    theme(axis.text.x = element_text(size=11), axis.text.y = element_text(size=11),
          axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))
  
  if(is.null(colors)){
    return(Fig)
  }else{
    return(Fig + scale_color_manual(values = colors))
  }
}

## For graphical purpose: Add a '+' before positive numbers ('-' before negative numbers, nothing before 0)

add_plus_sign <- function(x, digits = 2){
  if(length(x) > 1){
    return(sapply(x, add_plus_sign, digits = digits))
  }
  if(isTRUE(all.equal(x,0))){
    return(sprintf(paste0("% .",digits,"f"),x))
  }else{
    return(sprintf(paste0("%+.",digits,"f"),x))
  }
}
