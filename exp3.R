
# Third limitation: The $\rho_c$ is prone to the same problems as other linear correlation statistics.

library(viridis)

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


# first scatterplot
x = seq(1, 100, 0.1)
y1 <- 0.3*x
y2 <- 0.7*x
y3 <- 1.2*x
y4 <- 1.39*x

datacal1 <- rbind(data.frame(x=x, y=y1, val='sim1'),
                  data.frame(x=x, y=y2, val='sim2'),
                  data.frame(x=x, y=y3, val='sim3'),
                  data.frame(x=x, y=y4, val='sim4'))

p1 <- ggplot(datacal1, aes(x, y, color=val))+
  geom_line(linewidth = 1.5, show.legend = T, alpha = .8) +
  theme_bw() +  
  theme(text=element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        legend.position = c(.01, .99),
        legend.title = element_blank(),
        legend.justification = c("left", "top"),
        legend.box.just = "top",
        legend.margin = margin(6, 6, 6, 6),
        plot.margin = unit(c(0.3, 0.1, 0, 0), "cm"))+
  scale_x_continuous(expand = c(0,0), limits = c(0,100)) +  
  scale_y_continuous(expand = c(0,0), limits = c(0,100))+
  geom_abline()+
  xlab("Measured") + ylab("Predicted")+
  scale_colour_manual(labels  = c(TeX(r"(\it{\hat{z} = 0.3z})"), 
                                  TeX(r"(\it{\hat{z} = 0.7z})"), 
                                  TeX(r"(\it{\hat{z} = 1.2z})"), 
                                  TeX(r"(\it{\hat{z} = 1.4z})")),
                      values = viridis(4)) +     
  annotate("text", x = 75, y = 12, label = TeX(r"(\textit{r} = 1; $\rho_c$ = 0.23)"), color='black') + 
  annotate("text", x = 60, y = 30, label = TeX(r"(\textit{r} = 1; $\rho_c$ = 0.79)"), color='black') + 
  annotate("text", x = 50, y = 45, label = TeX(r"(\textit{r} = 1; $\rho_c$ = 0.94)"), color='black')+ 
  annotate("text", x = 38, y = 72, label = TeX(r"(\textit{r} = 1; $\rho_c$ = 0.82)"), color='black')
p1


# second scatterplot
x = seq(1, 100, 0.1)
y1 <- mean(x) + 0.1*(x- mean(x)) 
y2 <- mean(x) + 0.3*(x- mean(x)) 
y3 <- mean(x) + 0.8*(x- mean(x)) 

datacal2 <- rbind(data.frame(x=x, y=y1, val='sim1'),
                  data.frame(x=x, y=y2, val='sim2'),
                  data.frame(x=x, y=y3, val='sim3'))

p2 <- ggplot(datacal2, aes(x, y, color=val))+
  geom_line(linewidth = 1.5, show.legend = T, alpha = .8) +
  theme_bw() +  
  theme(text=element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        legend.position = c(.01, .99),
        legend.title = element_blank(),
        legend.justification = c("left", "top"),
        legend.box.just = "top",
        legend.margin = margin(6, 6, 6, 6),
        plot.margin = unit(c(0.3, 0.1, 0, 0), "cm"))+
  scale_x_continuous(expand = c(0,0), limits = c(0,100)) +  
  scale_y_continuous(expand = c(0,0), limits = c(0,100))+
  geom_abline()+
  xlab("Measured") + ylab("")+
  scale_colour_manual(labels  = c(TeX(r"(\it{\hat{z} = \bar{z} + 0.1(z- \bar{z})})"), 
                                  TeX(r"(\it{\hat{z} = \bar{z} + 0.3(z- \bar{z})})"), 
                                  TeX(r"(\it{\hat{z} = \bar{z} + 0.8(z- \bar{z})})")),
                      values = viridis(3)) +     
  annotate("text", x = 12, y = 50, label = TeX(r"(\textit{r} = 1; $\rho_c$ = 0.2)"), color='black') + 
  annotate("text", x = 14, y = 35, label = TeX(r"(\textit{r} = 1; $\rho_c$ = 0.55)"), color='black') + 
  annotate("text", x = 26, y = 20, label = TeX(r"(\textit{r} = 1; $\rho_c$ = 0.98)"), color='black')
p2


# third scatterplot
x = seq(1, 100, 0.1)
y1 <-  x - 15
y2 <-  x - 5  
y3 <-  x + 25

datacal3 <- rbind(data.frame(x=x, y=y1, val='sim1'),
                  data.frame(x=x, y=y2, val='sim2'),
                  data.frame(x=x, y=y3, val='sim3'))

p3 <- ggplot(datacal3, aes(x, y, color=val))+
  geom_line(linewidth = 1.5, show.legend = T, alpha = .8) +
  theme_bw() +  
  theme(text=element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        legend.position = c(.01, .99),
        legend.title = element_blank(),
        legend.justification = c("left", "top"),
        legend.box.just = "bottom",
        legend.margin = margin(6, 6, 6, 6),
        legend.text.align = 0,
        plot.margin = unit(c(0.3, 0.5, 0, 0), "cm"))+
  scale_x_continuous(expand = c(0,0), limits = c(0,100)) +  
  scale_y_continuous(expand = c(0,0), limits = c(0,100))+
  geom_abline()+
  xlab("Measured") + ylab("")+
  scale_colour_manual(labels  = c(TeX(r"(\it{\hat{z} = z - 15})"), 
                                  TeX(r"(\it{\hat{z} = z - 5})"), 
                                  TeX(r"(\it{\hat{z} = z + 25})")),
                      values = viridis(3)) +     
  annotate("text", x = 55, y = 28, label = TeX(r"(\textit{r} = 1; $\rho_c$ = 0.88)"), color='black') + 
  annotate("text", x = 43, y = 52, label = TeX(r"(\textit{r} = 1; $\rho_c$ = 0.98)"), color='black') + 
  annotate("text", x = 32, y = 72, label = TeX(r"(\textit{r} = 1; $\rho_c$ = 0.72)"), color='black')
p3
