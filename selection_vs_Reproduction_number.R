rm(list=ls())

library(tidyverse)
library(cowplot)

selection_gradient <- function(R, mu, sigma, d1=0, d2=0, d3=0){
  cv2 <- (sigma/mu)^2
  return((d2-2*d3)*(R^cv2-1)*mu/sigma^2 + R^cv2*(d1+2*(d3-d2)*log(R))/mu)
}

d <- 0.3 # phenotypic difference
mu <- 10 # mean generation time

Rw <- seq(0, 2, 0.01) # Effective reproduction number of the WT strain

selection <- data.frame("Rw" = Rw,
                        "s" = c(selection_gradient(R=Rw, mu=mu, sigma=mu, d1=d),
                                selection_gradient(R=Rw, mu=mu, sigma=mu, d2=d),
                                selection_gradient(R=Rw, mu=mu, sigma=mu, d3=d)))

selection$effect <- rep(1:3, each=length(Rw))

Fig <- ggplot(selection, aes(x=Rw, y=s, col=as.factor(effect))) +
  geom_hline(yintercept = 0, lwd = 0.5) +
  geom_vline(xintercept = 1, lwd = 0.75, color = 'grey') +
  geom_line(lwd=1) +
  labs(x = expression(paste(R[w](t), ", effective reproduction number of the WT strain")),
       y = "s(t), selection gradient\n",
       col = "Variant with:") +
  scale_color_manual(values = c("#428CCD", "#FF8000", "#36A213"),
                     labels = c(expression(paste(delta[1]>0, " (higher reproduction number)")),
                                expression(paste(delta[2]>0, " (longer mean generation time)")),
                                expression(paste(delta[3]>0, " (higher variance in generation time)")))) +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(plot.margin=unit(c(2.5,0.5,0.5,0.5),"cm"),
        axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12, vjust=-2),
        axis.title.y = element_text(size = 12),
        legend.text.align = 0, legend.position = c(0.31, 1.065))

# (selection gradient for d2 variant is max. when Rw = -exp(-0.5*(mu/sigma)^2))

arrow <- ggplot(data.frame("x"=c(0.5,1.5), "y" = c(0,0)), aes(x=x, y=y)) +
  geom_line(arrow = arrow(length=unit(0.30,"cm"), ends="first", type = "closed")) +
  xlim(c(0.3,1.7)) + ylim(c(-0.5,0.1)) +
  annotate(geom='text', label = "c(t), effectiveness of control measures", x = 1, y = -0.25, size = 4) +
  theme_void()

# pdf : landscape (6.34 x 5.60)
plot_grid(Fig, arrow, ncol = 1, align = 'v', rel_heights = c(0.9, 0.1))

# Approximation with Rw not too far from 1

selection_gradient_approx <- function(R, mu, sigma, d1=0, d2=0, d3=0){
  cv2 <- (sigma/mu)^2
  return((d2-2*d3)*(R^cv2-1)*mu/sigma^2 + R^cv2*(d1+2*(d3-d2)*(R-1))/mu) # log(R) -> R-1
}

selection$s_approx <- c(selection_gradient_approx(R=Rw, mu=mu, sigma=mu, d1=d),
                        selection_gradient_approx(R=Rw, mu=mu, sigma=mu, d2=d),
                        selection_gradient_approx(R=Rw, mu=mu, sigma=mu, d3=d))

plot_grid(Fig + geom_line(data=selection, aes(y=s_approx)) + ylim(layer_scales(Fig)$y$range$range), arrow,
          ncol = 1, align = 'v', rel_heights = c(0.9, 0.1))