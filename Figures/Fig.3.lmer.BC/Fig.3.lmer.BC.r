library(data.table)
library(ggplot2)
library(patchwork)
library(sf)
library(lme4)
library(performance)

setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")
cor.df<-readRDS("../Figures/Fig.3.t.test.BC/lmer.rda")
cor.df<-cor.df[var!="intercept" & type=="2m_temperature"]
cor.df[e.days==90, e.days:=91]
e.days_levels <- unique(cor.df$e.days)

cor.df$e.days <- factor(cor.df$e.days)
cor.df$var <- factor(cor.df$var, levels = c("accu_temp", "sd"))
rect_data <- cor.df[, .SD[which.max(R2_conditional)], 
                    by = .(phenology, growthform)]

rect_data <- rect_data[, .(phenology, growthform, e.days)]

rect_data[, `:=`(
  x_pos = as.numeric(factor(e.days, levels = e.days_levels)),
  ymin = 0.5,
  ymax = 2.5
)]

rect_data[, `:=`(
  xmin = x_pos - 0.5,
  xmax = x_pos + 0.5
)]


p <- ggplot(
  data = cor.df,
  aes(x = e.days, y = var, fill = Estimate)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_rect(
    data = rect_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    color = "grey", 
    linetype = "dashed",
    linewidth = 1,
    fill = NA
  ) +
  geom_text(
    aes(label = round(R2_conditional, 2)), 
    color = "black",
    size = 3
  )+

  facet_grid(growthform ~ phenology) +
  scale_fill_gradient2(
    low = "deepskyblue4", 
    mid = "white", 
    high = "firebrick", 
    midpoint = 0,
    name = "Slope"
  ) +
  coord_equal()+
  theme(
    legend.position = "right" 
  )
p

cor.df<-readRDS("../Figures/Fig.3.t.test.BC/lmer.rda")

ggplot(cor.df[var=="sd"])+
  geom_point(aes(x=e.days, y=R2_conditional, color=type))+
  facet_grid(growthform ~ phenology)
