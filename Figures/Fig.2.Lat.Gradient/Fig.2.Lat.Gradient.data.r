library(ggalluvial)
library(data.table)
library(ggplot2)
library(patchwork)
library(sf)
setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")
dist<-10
paired.records.filtered<-readRDS(sprintf("../Data/paired.records.clean.%dkm.rda", dist))

paired.differ<-unique(paired.records.filtered[, 
                                              c("e.invasive", "e.status", "e.scientific_name",
                                                "e.decimalLatitude", "e.decimalLongitude",
                                                "e.phenology", "e.bot_country", 
                                                "n.invasive", "n.status", "n.scientific_name",
                                                "n.phenology", "n.bot_country", 
                                                "growthform", "e.pheno_arc_fixed",
                                                "n.pheno_arc_fixed")])
paired.differ$pheno.diff<-paired.differ$e.pheno_arc_fixed - paired.differ$n.pheno_arc_fixed
hist(paired.differ$pheno.diff)
window.size<-2
step<-1
lats<-seq(floor(min(abs(paired.differ$e.decimalLatitude)))+window.size, 
          ceiling(max(abs(paired.differ$e.decimalLatitude)))-window.size,
          by=step)
lat<-lats[1]
all<-list()
for (lat in lats){
  item<-paired.differ[between(abs(e.decimalLatitude), lat-window.size, lat+window.size)]
  flowing<-item[e.phenology %in% c(1, 12)]
  flowing$phenology<-"Flowering"
  fruiting<-item[e.phenology %in% c(2, 12)]
  fruiting$phenology<-"Fruiting"
  item<-rbindlist(list(flowing, fruiting))
  item.N.all<-item[, .(pheno.diff=mean(pheno.diff),
                       pheno.diff.sd=sd(pheno.diff),
                       pheno.diff.ci.low = if (.N >= 2) t.test(pheno.diff)$conf.int[1] else NA_real_,
                       pheno.diff.ci.high = if (.N >= 2) t.test(pheno.diff)$conf.int[2] else NA_real_,
                       p = if (.N >= 2) t.test(pheno.diff)$p.value else NA_real_,
                       N=.N,
                       lat=lat,
                       growthform="All"),
                   by=list(phenology)]
  item.N.growthform<-item[, .(pheno.diff=mean(pheno.diff),
                       pheno.diff.sd=sd(pheno.diff),
                       pheno.diff.ci.low = if (.N >= 2) t.test(pheno.diff)$conf.int[1] else NA_real_,
                       pheno.diff.ci.high = if (.N >= 2) t.test(pheno.diff)$conf.int[2] else NA_real_,
                       p = if (.N >= 2) t.test(pheno.diff)$p.value else NA_real_,
                       N=.N,
                       lat=lat),
                   by=list(phenology, growthform)]
  item.result<-rbindlist(list(item.N.all, item.N.growthform), use.names=TRUE)
  all[[length(all)+1]]<-item.result
}
all.df<-rbindlist(all)
saveRDS(all.df, "../Figures/Fig.2.Lat.Gradient/Fig.2.Lat.Gradient.rda")
p<-ggplot(all.df[N>20])+
  geom_ribbon(aes(x=lat, ymin=pheno.diff.ci.low, ymax=pheno.diff.ci.high, fill=phenology), alpha=0.2)+
  geom_line(aes(x=lat, y=pheno.diff, color=phenology))+
  facet_wrap(~growthform, nrow=3)+
  theme_bw()
p

ggsave(p, filename="../Figures/Fig.2.Lat.Gradient/Fig.2.Lat.Gradient.png", width=6, height=5)
check<-paired.differ[between(abs(e.decimalLatitude), 22-window.size, 22+window.size)]
ggplot(check)+geom_histogram(aes(x=pheno.diff, fill=growthform), bins=50)+
  facet_grid(e.bot_country~e.invasive)

for (g in c("All", unique(paired.differ$growthform))){
  if (g!="All"){
    item<-paired.differ[growthform==g]
  }else{
    item<-paired.differ
  }
  flowing<-item[e.phenology %in% c(1, 12)]
  flowing$phenology<-"Flowering"
  fruiting<-item[e.phenology %in% c(2, 12)]
  fruiting$phenology<-"Fruiting"
  dt<-rbindlist(list(flowing, fruiting))
  dt$lat<-abs(dt$e.decimalLatitude)
  models <- dt[, {
    fit <- lm(pheno.diff ~ lat)
    r2 <- summary(fit)$r.squared
    pval <- summary(fit)$coefficients[2,4]
    list(model = list(fit), r2 = r2, pval = pval)
  }, by = phenology]
  
  label_pos <- dt[, .SD[which.max(lat)-5], by = phenology]
  
  label_pos <- merge(label_pos, models[, .(phenology, r2, pval)], by = "phenology")
  label_pos[, label := paste0("R² = ", round(r2, 2), 
                              "\nP = ", format.pval(pval, digits = 2, eps = 1e-3))]
  
  p<-ggplot(dt[sample(nrow(dt), 1e3)], aes(x = lat, y = pheno.diff, color = phenology)) +
    geom_point(alpha = 0.7, size=0.5) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_text(data = label_pos, 
              aes(label = label), 
              hjust = 0, vjust = 1, 
              size = 3.5, show.legend = FALSE) +
    labs(title = g,
         x = "Lat", y = "pheno.diff")
  ggsave(p, 
         filename=sprintf("../Figures/Fig.2.Lat.Gradient/Fig.2.Lat.Gradient.LM.%s.png", g), 
         width=6, height=4)
  
}

