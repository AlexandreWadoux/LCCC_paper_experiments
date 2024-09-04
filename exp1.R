
# First limitation: The $\rho_c$ does not inform on the individual contribution 
# of correlation and bias.

########################################
########################################
library(ggplot2)
library(DescTools)
library(sp)
library(ggsci)
library(latex2exp)
library(ggpubr)

# make a function to calculate the CCC
eval <- function(x, y){
  
  # concordance correlation coefficient
  n <- length(x)
  sdx <- sd(x, na.rm = T)
  sdy <- sd(y, na.rm = T)
  r <- stats::cor(x, y, method = 'pearson', use = 'pairwise.complete.obs')
  # scale shift
  v <- sdx / sdy
  sx2 <- var(x, na.rm = T) * (n - 1) / n
  sy2 <- var(y, na.rm = T) * (n - 1) / n
  # location shift relative to scale
  u <- (mean(x, na.rm = T) - mean(y, na.rm = T)) / ((sx2 * sy2)^0.25)
  Cb <- ((v + 1 / v + u^2)/2)^-1
  rCb <- r * Cb
  rhoC <- round(rCb, digits = 2)
  
  Cb <- round(Cb, digits = 2)
  r <- round(r, digits = 2)
  
  # return the results
  evalRes <- data.frame(rhoC = rhoC, Cb = Cb, r=r)
  
  return(evalRes)
}

# Reference field 
grid <- expand.grid(x1 = seq(0, 100, length.out = 100),
                    x2 = seq(0, 100, length.out = 100))

#compute spatial trend; x1 and x2 ares used as covariate, the linear model has an intercept
grid$mu <- 5 +  0.1*grid$x1 + 0.05*grid$x2

#define covariance function for simulation of residuals
covfun <- function(sill, range, Dist) {
  sill * exp(-Dist / range)
}

#compute matrix with distances between simulation nodes
distx<-outer(grid$x1,grid$x1,FUN="-")
disty<-outer(grid$x2,grid$x2,FUN="-")
Dist<-sqrt(distx^2+disty^2)
#Dist <- as.matrix(dist(grid[c('x1', 'x2')]^2))

#compute matrix with mean covariances
sill <- 5
range <- 10
C <- covfun(sill, range, Dist = Dist)

#simulate values for residuals by Cholesky decomposition
set.seed(31415)
Upper <- chol(C)
G <- rnorm(n=nrow(grid),0,1) #simulate random numbers from standard normal distribution
grid$residuals <- crossprod(Upper,G)

# add the trend and the Cholesky decomposed values
grid$y <- grid$mu+grid$residuals


sim1 <- -0.82460249 +2.15179855*grid$y-0.06611936*grid$y^2
sim2 <- 0.738113 + grid$y*1.305521 
sim3 <- 0.006630962*(grid$y-2.247695894)^3

# summarize; 
eval(x=grid$y, y=sim1)
eval(x=grid$y, y=sim2)
eval(x=grid$y, y=sim3)

reference <- data.frame(x1=grid$x1, x2=grid$x2, Value=grid$y, name = "Reference")
sim1 <- data.frame(x1=grid$x1, x2=grid$x2, Value=sim1, name = "sim1")
sim2 <- data.frame(x1=grid$x1, x2=grid$x2, Value=sim2, name = "sim2")
sim3 <- data.frame(x1=grid$x1, x2=grid$x2, Value=sim3, name = "sim3")


grid.plot <- rbind(reference, sim1, sim2, sim3)
grid.plot$name <- factor(grid.plot$name ,
                         levels=c("Reference", "sim1","sim2","sim3"),
                         labels=c("Reference",
                                  "Prediction 1",
                                  "Prediction 2",
                                  "Prediction 3"))


ggplot(grid.plot) + geom_tile(aes(x1, x2, fill = Value)) +
  scale_fill_distiller(
    palette = "Spectral",
    limits = c(-5, 35),
    na.value = "#BD0026")+ 
  facet_grid(~name)+
  theme_minimal() + 
  coord_fixed()+
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    strip.text.x = element_text(size = 10)
  )



datsim1 <- data.frame(y=grid$y, sim = -0.82460249 +2.15179855*grid$y-0.06611936*grid$y^2)
p1 <- ggplot(datsim1, aes(y, sim))+
  geom_point(shape = 16, size = 1.5, show.legend = FALSE, alpha = .5) +
  theme_bw() +  
  theme(text=element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)))+
  scale_x_continuous(expand = c(0,0), limits = c(0,35)) +  
  scale_y_continuous(expand = c(0,0), limits = c(0,35))+
  geom_abline()+
  xlab("Reference") + ylab("Prediction")
p1

datsim2 <- data.frame(y=grid$y, sim = 0.738113 + grid$y*1.305521 )
p2 <- ggplot(datsim2, aes(y, sim))+
  geom_point(shape = 16, size = 1.5, show.legend = FALSE, alpha = .5) +
  theme_bw() +  
  theme(text=element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)))+
  scale_x_continuous(expand = c(0,0), limits = c(0,35)) +  
  scale_y_continuous(expand = c(0,0), limits = c(0,35))+
  geom_abline()+
  xlab("Reference") + ylab("")     
p2


datsim4 <- data.frame(y=grid$y, sim = 0.006630962*(grid$y-2.247695894)^3)
# plot scatterplot
p4 <- ggplot(datsim4, aes(y, sim))+
  geom_point(shape = 16, size = 1.5, show.legend = FALSE, alpha = .5) +
  theme_bw() +  
  theme(text=element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)))+
  scale_x_continuous(expand = c(0,0), limits = c(0,35)) +  
  scale_y_continuous(expand = c(0,0), limits = c(0,35))+
  geom_abline()+
  xlab("Reference") + ylab("")    
p4

ggarrange(p1, p2, p4, nrow=1)
