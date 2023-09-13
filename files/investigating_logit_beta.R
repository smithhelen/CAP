## investigating the logit line for different values of beta

library(boot)
library(tidyverse)
library(patchwork)

beta <- seq(from=1, to=20, by=1)
x <- seq(from=-1, to=1, by=0.01)
y <- map(beta, function(z) inv.logit(z*x))
p <- map2(y, beta, function(z,w) ggplot(data.frame(x=x,y=z),aes(x=x,y=z))+geom_line()+ylim(0,1)+ggtitle(paste0("beta=",w)))
names(p) <- paste0("p",factor(1:20))
wrap_plots(p)

png("beta.png", width=640, height=480)
wrap_plots(p)
dev.off()
