library(tidyverse)
library(sf)
library(castor)
library(lubridate)
library(U.PhyloMaker)
library(do)
library(metafor)
library(ggbeeswarm)
library(stringr)
library(knitr)
library(broom)
library(kableExtra)

setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")
iNadata <- read_csv("../Data/Similar_phenology_invasive_code&data_20250218/iNadata.csv") 
iNadata[which(iNadata$id==548648),]
glimpse(iNadata)

#Select species with penological information
iNadata <- iNadata %>% filter(phenology != 0)

#Correct the species names (aligning with Zanne’s tree from Zanne et al., Nature: 506, 89–92 (2014) (https://doi.org/10.1038/nature12872))
iNadata$scientific_name <- gsub("Chamaenerion angustifolium", "Epilobium angustifolium", iNadata$scientific_name)
iNadata$scientific_name <- gsub("Chamaenerion latifolium", "Epilobium latifolium", iNadata$scientific_name)
iNadata$scientific_name <- gsub("Polygonella articulata", "Polygonum articulatum", iNadata$scientific_name)
iNadata$scientific_name <- gsub("Fallopia cilinodis", "Polygonum cilinode", iNadata$scientific_name)
iNadata$scientific_name <- gsub("Polygonella americana", "Polygonum serotinum", iNadata$scientific_name)
iNadata$scientific_name <- gsub("Millettia pinnata", "Pongamia pinnata", iNadata$scientific_name)

#2.1.2 Remove images before 2008

year_dat <- iNadata %>% filter(phenology != 0)  %>%
  group_by(LEVEL3_COD, invasive, scientific_name, status) %>%
  filter(n() >= 20) %>% mutate(Year = as.numeric(format(as.Date(observed_on) , "%Y"))) %>%
  mutate(Year = ifelse(Year < 2008, 2008, Year)) 
glimpse(year_dat)

#Since iNaturalist was founded in 2008, we only included images taken after 2008.

iNadata$observed_on <- ymd(iNadata$observed_on)
iNadata_before_2008 <- iNadata %>% filter(observed_on < as.Date("2008-01-01"))
iNadata_ori <- dim(iNadata)[1]
iNadata <- iNadata %>% filter(observed_on > as.Date("2008-01-01"))

#Calculate number of observations before 2008

dim(iNadata_before_2008)[1]/iNadata_ori

## [1] 0.00507492


#Only around 0.5% images were taken before 2008.

year_distribution <- ggplot(data = year_dat, aes(x = Year)) + 
  geom_histogram(aes(y = ..count../sum(..count..)), binwidth = 1, col = I("black"), fill = "lightgray", alpha = 0.5) +
  labs(x = "Observation year", y = "Proportion") +
  theme_classic() + 
  scale_y_continuous(limits = c(0, 0.25),
                     labels = scales::percent_format(accuracy = 1)) + 
  scale_x_continuous(breaks = seq(2008, 2023, 1),
                     labels = c("Before 2008", seq(2009, 2023, 1))) + 
  theme(axis.text = element_text(colour = "black", size = 18), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(margin = margin(r = 15)),  
        legend.background = element_rect(fill = "transparent"),
        title = element_text(colour = "black", size = 18),
        axis.ticks = element_line(color = "black"),
        axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5),
        plot.margin = margin(10, 20, 10, 10)) +
  geom_vline(xintercept = quantile(year_dat$Year, 0.10), linetype = 2)
year_distribution


#Covert the date to an angular measure

iNadata <- iNadata %>% 
  mutate(pheno_arc = as.numeric(format(as.Date(observed_on, format = "%Y-%m-%d"), format = "%j")) * 360 / 365)

#Defining the function to calculate the sd and mean of circular vector

vector_mean <- function(x) {
  c <- sum(cos(x * pi / 180)) / length(x)
  s <- sum(sin(x * pi / 180)) / length(x)
  angle <- atan2(s, c) * 180 / pi
  if (angle < 0) { angle <- angle + 360 }
  m_d <- angle
}

vector_sd <- function(x) {
  c <- sum(cos(x * pi / 180)) / length(x)
  s <- sum(sin(x * pi / 180)) / length(x)
  angle <- atan2(s, c) * 180 / pi
  s.d <- sqrt(-2 * log(sqrt(s^2 + c^2))) * 180 / pi
  d <- s.d
}


#2.1.4 Flowering data

#Extract flowering data

flower <- iNadata %>% 
  filter(phenology %in% c(1,12)) %>%
  group_by(LEVEL3_COD, invasive, scientific_name, status) %>%
  summarise(n = n(), mean = vector_mean(pheno_arc), sd = vector_sd(pheno_arc)) %>%
  filter(n >= 20)

#Flowering of native species

native_flower <- flower %>% filter(status == "native") %>% as.data.frame()
glimpse(native_flower)

#Rename columns for clarification

names(native_flower)[3] <- "n.species"
names(native_flower)[5:7] <- c("n.n", "n.mean", "n.sd")
glimpse(native_flower)

#Flowering of exotic

exotic_flower <- flower %>% filter(status == "exotic") %>% as.data.frame()
glimpse(exotic_flower)

#Rename columns for clarity

names(exotic_flower)[3] <- "e.species"
names(exotic_flower)[5:7] <- c("e.n", "e.mean", "e.sd")
glimpse(exotic_flower)

#Pair the invasive species with their relative natives

n <- 0
flower_dis <- tibble()
for (i in 1:length(exotic_flower$LEVEL3_COD)) {
  for (j in 1:length(native_flower$LEVEL3_COD)) {
    if (exotic_flower$LEVEL3_COD[i] == native_flower$LEVEL3_COD[j] && 
        exotic_flower$invasive[i] == native_flower$invasive[j]) {
      n <- n + 1
      flower_dis[n, 1:6] <- exotic_flower[i, c(1:3, 5:7)]
      flower_dis[n, 7:10] <- native_flower[j, c(3, 5:7)]
    }
  }
}

#Defining the sequence of phenology (short arc)

for (i in 1:length(flower_dis$LEVEL3_COD)) {
  if (abs(flower_dis$e.mean[i] - flower_dis$n.mean[i]) <= 
      360 - abs(flower_dis$e.mean[i] - flower_dis$n.mean[i])) {
    flower_dis$e.vector[i] <- flower_dis$e.mean[i]
    flower_dis$n.vector[i] <- flower_dis$n.mean[i]
    
  } else {
    if (flower_dis$e.mean[i] > flower_dis$n.mean[i]) { 
      flower_dis$e.vector[i] <- flower_dis$e.mean[i]
      flower_dis$n.vector[i] <- flower_dis$n.mean[i] + 360
    } else {
      flower_dis$e.vector[i] <- flower_dis$e.mean[i] + 360
      flower_dis$n.vector[i] <- flower_dis$n.mean[i]
    }
  }
}

#Display the flowering data

glimpse(flower_dis)

#2.1.5 Fruiting data

#Extract fruiting phenological data

fruit <- iNadata %>%
  filter(phenology %in% c(2,12)) %>%
  group_by(LEVEL3_COD, invasive, scientific_name, status) %>%
  summarise(n = n(), mean = vector_mean(pheno_arc), sd = vector_sd(pheno_arc)) %>%
  filter(n >= 20)

#Fruiting of native species

native_fruit <- fruit %>% filter(status == "native") %>% as.data.frame()
glimpse(native_fruit)

#Rename columns for clarification

names(native_fruit)[3] <- "n.species"
names(native_fruit)[5:7] <- c("n.n", "n.mean", "n.sd")
glimpse(native_fruit)


#Fruiting of exotic species

exotic_fruit <- fruit %>% filter(status == "exotic") %>% as.data.frame()
glimpse(exotic_fruit)

#Rename columns for clarity

names(exotic_fruit)[3] <- "e.species"
names(exotic_fruit)[5:7] <- c("e.n", "e.mean", "e.sd")
glimpse(exotic_fruit)

#Pair the invasive species with their relative natives

n <- 0
fruit_dis <- tibble()
for (i in 1:length(exotic_fruit$LEVEL3_COD)) {
  for (j in 1:length(native_fruit$LEVEL3_COD)) {
    if (exotic_fruit$LEVEL3_COD[i] == native_fruit$LEVEL3_COD[j] && 
        exotic_fruit$invasive[i] == native_fruit$invasive[j]) {
      n <- n + 1
      fruit_dis[n, 1:6] <- exotic_fruit[i, c(1:3, 5:7)]
      fruit_dis[n, 7:10] <- native_fruit[j, c(3, 5:7)]
    }
  }
}

#Defining the sequence of phenology (short arc)

for (i in 1:length(fruit_dis$LEVEL3_COD)) {
  if (abs(fruit_dis$e.mean[i] - fruit_dis$n.mean[i]) <= 
      360 - abs(fruit_dis$e.mean[i] - fruit_dis$n.mean[i])) {
    fruit_dis$e.vector[i] <- fruit_dis$e.mean[i]
    fruit_dis$n.vector[i] <- fruit_dis$n.mean[i]
  } else {
    if (fruit_dis$e.mean[i] > fruit_dis$n.mean[i]) { 
      fruit_dis$e.vector[i] <- fruit_dis$e.mean[i]
      fruit_dis$n.vector[i] <- fruit_dis$n.mean[i] + 360
    } else {
      fruit_dis$e.vector[i] <- fruit_dis$e.mean[i] + 360
      fruit_dis$n.vector[i] <- fruit_dis$n.mean[i]
    }
  }
}
glimpse(fruit_dis)

#2.1.6 Get continent information

#Load WCVP checklist

checklist <- read.table("../Data/Similar_phenology_invasive_code&data_20250218/wcvp_names.csv", sep = "|", header = TRUE, quote = "", 
                        fill = TRUE, encoding = "UTF-8") %>% 
  dplyr::select("plant_name_id", "taxon_status", "family", 
                "genus", "species", "lifeform_description", 
                "taxon_name", "hybrid_formula") %>% 
  arrange(plant_name_id)

distribution <- read.table("../Data/Similar_phenology_invasive_code&data_20250218/wcvp_distribution.csv", sep = "|", header = TRUE, quote = "", 
                           fill = TRUE, encoding = "UTF-8") %>% 
  arrange(plant_name_id)

#Merge checklist and distribution data

plantcountry <- left_join(distribution, checklist, by = "plant_name_id") %>% 
  arrange(area) %>% 
  filter(area != "" & hybrid_formula == "" & species != "") %>% 
  mutate(scientific_name = paste0(genus, " ", species)) %>% 
  dplyr::select("scientific_name", "taxon_name", everything())

plant <- plantcountry %>% dplyr::select("scientific_name", "genus", "family") %>%  
  mutate(specie.relative = "", genus.relative = "") %>% 
  unique()

#Merge the continental information

continent <- plantcountry %>% 
  group_by(continent, area_code_l3) %>% 
  summarise(n = n())

#Merge the continent information

flower_dis <- left_join(flower_dis, continent, by = c("LEVEL3_COD" = "area_code_l3"))
fruit_dis <- left_join(fruit_dis, continent, by = c("LEVEL3_COD" = "area_code_l3")) 


#2.2 Extract coordinates
#2.2.1 Load the data

#The original data of botanical country boundary (https://www.tdwg.org/standards/wgsrpd/)

bot_country <- read_sf("../Data/Similar_phenology_invasive_code&data_20250218/geographic boundary data/bot_country_ori.shp")

#2.2.2 Flowering

for (i in 1:length(flower_dis$LEVEL3_COD)) {
  cor <- bot_country[which(bot_country$LEVEL3_COD == flower_dis$LEVEL3_COD[i]), ]
  flower_dis$longitude[i] <- st_coordinates(st_centroid(cor$geometry))[1]
  flower_dis$latitude[i] <- st_coordinates(st_centroid(cor$geometry))[2]
}

#2.2.3 Fruiting

for (i in 1:length(fruit_dis$LEVEL3_COD)) {
  cor <- bot_country[which(bot_country$LEVEL3_COD == fruit_dis$LEVEL3_COD[i]), ]
  fruit_dis$longitude[i] <- st_coordinates(st_centroid(cor$geometry))[1]
  fruit_dis$latitude[i] <- st_coordinates(st_centroid(cor$geometry))[2]
}

#2.3 Get the growth form
#2.3.1 Flowering

herbaceous <- c("Reynoutria japonica", 
                "Hedychium gardnerianum",
                "Pueraria montana",
                "Mikania micrantha", 
                "Lythrum salicaria",
                "Chromolaena odorata",
                "Sphagneticola trilobata")

flower_dis$growthform <- if_else(flower_dis$invasive %in% herbaceous, "Herbaceous", "Woody")

#2.3.2 Fruiting

herbaceous <- c("Reynoutria japonica", 
                "Hedychium gardnerianum",
                "Pueraria montana",
                "Mikania micrantha", 
                "Lythrum salicaria",
                "Chromolaena odorata",
                "Sphagneticola trilobata")

fruit_dis$growthform <- if_else(fruit_dis$invasive %in% herbaceous, "Herbaceous", "Woody")

#2.4 Final data tabel
#2.4.1 Flowering

#Standardise the data

flower_dis$n.vector <- as.numeric(flower_dis$n.vector)
flower_dis$e.vector <- as.numeric(flower_dis$e.vector)
flower_dis$id <- 1:length(flower_dis$LEVEL3_COD)
flower_dis$latitude.abs <- abs(flower_dis$latitude)

#Calculate effect size

flower_dis <- escalc(
  measure = "SMD", 
  m1i = n.vector, 
  m2i = e.vector,
  sd1i = n.sd, 
  sd2i = e.sd, 
  n1i = n.n, 
  n2i = e.n, 
  data = flower_dis
)

#2.4.2 Fruiting

#Standardise the data

fruit_dis$n.vector <- as.numeric(fruit_dis$n.vector)
fruit_dis$e.vector <- as.numeric(fruit_dis$e.vector)
fruit_dis$id <- 1:length(fruit_dis$LEVEL3_COD)
fruit_dis$latitude.abs <- abs(fruit_dis$latitude)

#Calculate effect size

fruit_dis <- escalc(
  measure = "SMD", 
  m1i = n.vector, 
  m2i = e.vector,
  sd1i = n.sd, 
  sd2i = e.sd, 
  n1i = n.n, 
  n2i = e.n, 
  data = fruit_dis
)

#2.4.3 Combine and save the data

#Combine flowering and fruiting data into a single table

meta <- rbind(
  flower_dis %>% mutate(phenology = "Flowering"), 
  fruit_dis %>% mutate(phenology = "Fruiting")
) %>% 
  mutate(id = 1:length(LEVEL3_COD)) %>% 
  dplyr::select(id, -e.mean, -n.mean, everything()) %>% 
  dplyr::select(-e.mean, -n.mean)

#Save the data

meta$continent <- str_to_title(meta$continent)
write_csv(meta, "Analysis_data.csv")

#3 Data analyses
#3.1 Overall phenological difference
#3.1.1 Number of corresponding native species for each invasive species (Flowering)
invasive_native <- meta %>% group_by(invasive,n.species, phenology) %>% summarise(n = n()) %>% 
  group_by(invasive, phenology) %>% summarise(n = n())

invasive_native_flow <- ggplot(data = invasive_native %>% 
                                 filter(phenology == "Flowering"), aes(n)) + 
  geom_histogram(binwidth = 1, col = I("black"), fill = "#FBB4AE", alpha = 0.5) +
  labs(x = "Number of invasive-native pairs", y = "Frequency") +
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4.5)) +
  scale_x_continuous(limits = c(0, 80))+
  theme(axis.text = element_text(colour = "black", size = 18), 
        legend.background = element_rect("transparent"),
        title = element_text(colour = "black", size = 18),
        axis.ticks = element_line(color = "black"),
        axis.text.y = element_text(angle = 90)) +
  geom_vline(xintercept = median(invasive_native[invasive_native$phenology == "Flowering", ]$n), linetype = 2)

#Number of relative natives for each invasive plants for flowering

invasive_native_flow 

#3.1.2 Number of corresponding native species for each invasive species (Fruiting)
invasive_native_frui <- ggplot(data = invasive_native %>% filter(phenology == "Fruiting"), aes(n)) + 
  geom_histogram(binwidth = 1, col = I("black"), fill = "#FFD92F", alpha = 0.5) +
  labs(x = "Number of invasive-native pairs", y = "Frequency") +
  theme_classic() + 
  scale_y_continuous(limits = c(0, 4.5), breaks = c(0, 1, 2, 3, 4), expand = c(0, 0)) +
  theme(axis.text = element_text(colour = "black", size = 18), 
        legend.background = element_rect("transparent"),
        title = element_text(colour = "black", size = 18),
        axis.text.y = element_text(angle = 90)) +
  geom_vline(xintercept = median(invasive_native[invasive_native$phenology == "Fruiting", ]$n), linetype = 2)
#Number of relative natives for each invasive plants for fruiting

invasive_native_frui

#3.1.3 The meta regression model
#Load the data

phenology_level <- c("Flowering", "Fruiting")
meta <- read.csv("Analysis_data.csv") %>% 
  mutate(phenology = factor(phenology, levels = rev(phenology_level)))
#Run the model

## Perform multilevel meta-analysis on phenology (Flowering and Fruiting)
ina_phenology <- rma.mv(
  yi, vi, 
  mods = ~phenology - 1,  ## Model phenology as a fixed effect (no intercept)
  random = list(           ## Random effects structure
    ~ 1 | LEVEL3_COD / id, ## Nested random effect for LEVEL3_COD and id
    ~ 1 | n.species,       ## Random effect for native species
    ~ 1 | e.species        ## Random effect for invasive species
  ),
  data = meta              ## Use the meta data for fitting the model
)
summary(ina_phenology)


#Format the results

## Create a summary of the model results for Flowering and Fruiting phenology
ina_result <- data.frame(
  phenology = c("Flowering", "Fruiting"),           ## Phenology types
  estimate = ina_phenology$b,                       ## Effect estimates (beta coefficients)
  lower_95_CI = ina_phenology$ci.lb,                ## Lower bound of the 95% CI
  upper_95_CI = ina_phenology$ci.ub,                ## Upper bound of the 95% CI
  test_statistic = ina_phenology$pval,              ## P-value for significance testing
  row.names = NULL
) %>% 
  mutate(phenology = factor(phenology, levels = rev(phenology_level))) ## Set factor levels for the plot
#Phenological difference between invasive species and native relatives

kable(ina_result, digits = 3, align = "c") %>%
  kable_styling(bootstrap_options = "bordered", font_size = 15, latex_options = "scale_down") %>%
  row_spec(0, background = "gray", color = "white") %>%
  column_spec(1, background = "gray", color = "white")


#3.1.4 Visualise overall phenological difference
#Plot

bot_plot <- ggplot(meta, aes(x = phenology, y = yi)) + 
  # Adding a horizontal line at y = 0
  geom_hline(yintercept = 0, linetype = 2) +
  
  # Adding a quasirandom scatter plot with fill by phenology
  geom_quasirandom(aes(fill = phenology), shape = 21, alpha = 0.5) + 
  
  # Adding point ranges from the results (estimate and confidence intervals)
  geom_pointrange(data = ina_result, 
                  mapping = aes(x = phenology, 
                                y = estimate, 
                                ymin = lower_95_CI,
                                ymax = upper_95_CI, 
                                fill = phenology), 
                  shape = 21, 
                  color = "black", 
                  size = 0.75) +
  
  # Set the Y-axis limits and label
  scale_y_continuous(name = "Phenological difference", limits = c(-10, 10)) + 
  
  # Set the fill color for different phenological stages
  scale_fill_manual(values = c("Flowering" = "#FBB4AE", 
                               "Fruiting" = "#FFD92F")) +  
  
  # Using classic theme for clean plot design
  theme_classic() + 
  
  # Customize axis text and titles
  theme(
    axis.text = element_text(size = 18, color = "black"),
    axis.title.x.bottom = element_text(size = 18),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 18, angle = 90, hjust = 0.5)
  ) +
  
  # Flip the coordinates to make it horizontal
  coord_flip() + 
  
  # Remove legends for fill and color
  guides(fill = FALSE, color = FALSE)
#Flowering and fruiting differences between invasive plants and their native relatives

bot_plot

#3.2 Different invasive species
#3.2.1 Flowering
#Function to extract the abbreviation of invasive species

extract_letters <- function(species_name) {
  parts <- str_split(species_name, " ")[[1]]
  genus <- str_sub(parts[1], 1, 2)
  species <- str_sub(parts[2], 1, 2)
  return(str_c(genus, species, sep = ""))
}
#Run the model

invasive_flow <- rma.mv(yi,vi, mods = ~invasive - 1,
                        random = list(~ 1 | LEVEL3_COD / id,
                                      ~ 1 | n.species),
                        data = meta,
                        subset = phenology == "Flowering")
summary(invasive_flow)

#Format the results

invasive_flow_result <- left_join(
  data.frame(invasive = Replace(row.names(as.data.frame(invasive_flow$b)), "invasive",""),
             estimate = invasive_flow$b,
             lower_95_CI = invasive_flow$ci.lb,
             upper_95_CI = invasive_flow$ci.ub,
             test_statistic = invasive_flow$pval,
             row.names = NULL), meta %>% 
    dplyr::select(invasive, growthform) %>% 
    unique, by = "invasive") %>% 
  mutate(growthforms = factor(growthform, levels = c("Woody", "Herbaceous"))) %>% 
  dplyr::select(-growthform) %>% 
  arrange(growthforms, estimate) %>% 
  mutate(spec_code = as.character(sapply(invasive, extract_letters))) %>% 
  dplyr::select(invasive, spec_code, growthforms, estimate, 
                lower_95_CI, upper_95_CI, test_statistic)

invasive_flow_result$spec_code = factor(invasive_flow_result$spec_code, 
                                        levels = invasive_flow_result$spec_code)
kable(invasive_flow_result, digits = 3, align = "c") %>%
  kable_styling(bootstrap_options = "bordered", font_size = 15, latex_options = "scale_down") %>%
  row_spec(0, background = "gray", color = "white") %>%
  column_spec(1, background = "gray", color = "white")


#Visualise flowering difference across invasive species

spec_flow <- ggplot() + 
  geom_hline(yintercept = 0, linetype = 2) +
  geom_quasirandom(meta %>% 
                     mutate(spec_code = factor(sapply(invasive, extract_letters), 
                                               levels = invasive_flow_result$spec_code)) %>% 
                     arrange(yi) %>% filter(phenology == "Flowering") %>% 
                     mutate(growthforms = growthform),
                   mapping = aes(x = spec_code, y = yi, 
                                 fill = growthforms, shape = growthforms), alpha = 0.5) + 
  geom_pointrange(data = invasive_flow_result, 
                  mapping = aes(x = spec_code, y = estimate, ymin = lower_95_CI, 
                                ymax = upper_95_CI, fill = growthforms, shape = growthforms), 
                  color = "black", size = 0.75) +
  scale_y_continuous(name = "Phenological difference", limits = c(-11, 11)) + 
  scale_fill_manual(values = c("Herbaceous" = "#66A61E", "Woody" = "#A6761D")) + 
  scale_shape_manual(values = c("Herbaceous" = 23, "Woody" = 21)) + 
  theme_classic() + 
  theme(plot.title = element_text(size = 18, hjust = 0.68),
        axis.text = element_text(size = 18, color = "black"),
        axis.title.x.bottom = element_text(size = 18),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank()) +
  coord_flip() + 
  guides(color = F, fill = F, shape = F)
#Flowering difference across invasive species

spec_flow

#3.2.2 Fruiting
#Run the model

invasive_frui <- rma.mv(yi,vi, mods = ~invasive - 1,
                        random = list(~ 1 | LEVEL3_COD / id,
                                      ~ 1 | n.species),
                        data = meta,
                        subset = phenology == "Fruiting")
#Format the results

invasive_frui_result <- left_join(data.frame(
  invasive = Replace(row.names(as.data.frame(invasive_frui$b)), "invasive",""),
  estimate = invasive_frui$b,
  lower_95_CI = invasive_frui$ci.lb,
  upper_95_CI = invasive_frui$ci.ub,
  test_statistic = invasive_frui$pval,
  row.names = NULL), meta %>% 
    dplyr::select(invasive, growthform) %>% 
    unique, by = "invasive") %>% 
  mutate(growthforms = factor(growthform, levels = c("Woody", "Herbaceous"))) %>% 
  dplyr::select(-growthform) %>% 
  arrange(growthforms, estimate) %>% 
  mutate(spec_code = as.character(sapply(invasive, extract_letters)))

invasive_frui_result$spec_code = factor(invasive_frui_result$spec_code, 
                                        levels = invasive_frui_result$spec_code)

kable(invasive_frui_result, digits = 3, align = "c") %>%
  kable_styling(bootstrap_options = "bordered", font_size = 15, latex_options = "scale_down") %>%
  row_spec(0, background = "gray", color = "white") %>%
  column_spec(1, background = "gray", color = "white")

#Visualise fruiting difference across invasive species

spec_frui <- ggplot() + 
  geom_hline(yintercept = 0, linetype = 2) +
  geom_quasirandom(meta %>% 
                     mutate(spec_code = factor(sapply(invasive, extract_letters), 
                                               levels = invasive_frui_result$spec_code)) %>% 
                     arrange(yi) %>% filter(phenology == "Fruiting") %>% 
                     mutate(growthforms = growthform),
                   mapping = aes(x = spec_code, y = yi, fill = growthforms, shape = growthforms),
                   alpha = 0.5) + 
  geom_pointrange(data = invasive_frui_result, 
                  mapping = aes(x = spec_code, y = estimate, ymin = lower_95_CI, 
                                ymax = upper_95_CI, fill = growthforms, shape = growthforms), 
                  color = "black", size = 0.75) +
  scale_y_continuous(name = "Phenological difference", limits = c(-11, 11)) + 
  scale_fill_manual(values = c("Herbaceous" = "#66A61E", "Woody" = "#A6761D")) + 
  scale_shape_manual(values = c("Herbaceous" = 23, "Woody" = 21)) +
  theme_classic() + 
  theme(plot.title = element_text(size = 18, hjust = 0.68),
        axis.text = element_text(size = 18, color = "black"),
        axis.title.x.bottom = element_text(size = 18),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank()) +
  coord_flip() + guides(color = F, fill = F, shape = F)
#Fruiting difference across invasive species

spec_frui

#3.3 Botanical country patterns
#3.3.1 Pair numerber of botanical countries (Flowering)
#Plot

pairs <- meta %>% group_by(LEVEL3_COD, phenology) %>% summarise(n = n())

pairs_flow <- ggplot(data = pairs %>% filter(phenology == "Flowering"), aes(n)) + 
  geom_histogram(binwidth = 1, col = I("black"), fill = "#FBB4AE", alpha = 0.5) +
  labs(x = "Number of invasive-native pairs", y = "Frequency") +
  theme_classic() + 
  scale_y_continuous(limits = c(0, 35), expand = c(0, 0)) +
  theme(axis.text = element_text(colour = "black", size = 18), 
        legend.background = element_rect("transparent"),
        title = element_text(colour = "black", size = 18),
        axis.ticks = element_line(color = "black"),
        axis.text.y = element_text(angle = 90)) +
  geom_vline(xintercept = median(pairs[pairs$phenology == "Flowering", ]$n), linetype = 2)
#Number of studied species pairs for flowring difference in each botanical country

pairs_flow

#3.3.2 Pair numerber of botanical countries (Fruiting)
#Plot

pairs_frui <- ggplot(data = pairs %>% filter(phenology == "Fruiting"), aes(n)) + 
  geom_histogram(binwidth = 1, col = I("black"), fill = "#FFD92F", alpha = 0.5) +
  labs(x = "Number of invasive-native pairs", y = "Frequency") +
  theme_classic() + 
  scale_y_continuous(limits = c(0, 35), expand = c(0, 0)) +
  theme(axis.text = element_text(colour = "black", size = 18), 
        legend.background = element_rect("transparent"),
        title = element_text(colour = "black", size = 18),
        axis.text.y = element_text(angle = 90)) +
  geom_vline(xintercept = median(pairs[pairs$phenology == "Fruiting", ]$n), linetype = 2)
#Number of studied species pairs for fruiting difference in each botanical country

pairs_frui

#3.3.3 Frequency of phenology shifting direction
#Aggregrate data into botanical country level

bot_level <- meta %>% 
  mutate(phenology = factor(phenology, levels = rev(phenology_level))) %>% 
  mutate(low_CI = yi - qnorm(0.975)*sqrt(vi),
         up_CI = yi + qnorm(0.975)*sqrt(vi)) %>% 
  arrange(yi) %>% 
  mutate(fill = ifelse(low_CI < 0 & up_CI > 0, "neutral", ifelse(low_CI >= 0, "positive", "negative"))) %>% 
  group_by(phenology) %>% 
  mutate(rank = 1:(length(fill)))
#Count the frequency

Phenology_bar_data <- as.data.frame(bot_level %>% 
                                      group_by(phenology, fill) %>% summarise(n = n()))

bar_data <- Phenology_bar_data %>% 
  group_by(phenology) %>% 
  mutate(phenology = factor(phenology, levels = c("Flowering", "Fruiting")),
         proportion = n/sum(n))
#Frequency of phenology shifting direction in different botanical country

kable(bar_data, digits = 3, align = "c") %>%
  kable_styling(bootstrap_options = "bordered", font_size = 15) %>%
  row_spec(0, background = "gray", color = "white")
#Plot

bar <- ggplot() +
  geom_bar(bar_data, 
           mapping = aes(x = phenology, y = proportion, fill = fill), 
           stat = "identity", width = 0.5, color = "black", position = "stack") +  
  scale_fill_manual(values = c("neutral" = "#FFFFBF",
                               "positive" = "#91CF60",
                               "negative" = "#FC8D59"), guide = F) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_x_discrete(expand = c(0.7, 0)) +
  labs(x = NULL, y = "Frequency of invasive-native pairs") + 
  theme(axis.ticks.x = element_blank(),
        axis.text = element_text(size = 18, color = "black"),
        title = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(color = "black", fill = "white", linewidth = 1)) +
  annotate("text", x = c(2, 2, 2, 1, 1, 1), 
           y = c(bar_data[bar_data$phenology == "Fruiting",]$proportion[3]/2,
                 bar_data[bar_data$phenology == "Fruiting",]$proportion[3] + 
                   bar_data[bar_data$phenology == "Fruiting",]$proportion[2]/2,
                 1 - bar_data[bar_data$phenology == "Fruiting",]$proportion[1]/2,
                 bar_data[bar_data$phenology == "Flowering",]$proportion[3]/2,
                 bar_data[bar_data$phenology == "Flowering",]$proportion[3] + 
                   bar_data[bar_data$phenology == "Flowering",]$proportion[2]/2,
                 1 - bar_data[bar_data$phenology == "Flowering",]$proportion[1]/2), 
           label = c(bar_data$n[3], bar_data$n[2], bar_data$n[1], bar_data$n[6], bar_data$n[5], bar_data$n[4]),
           size = 18/.pt)
#Frequency of phenology shifting direction in different botanical country

bar

#3.3.4 Overall
#Run the model

pattern <- rma.mv(yi, vi, mods = ~ LEVEL3_COD - 1,
                  random = list(~ 1 | n.species, 
                                ~ 1 | e.species),
                  data = meta)
summary(pattern)

#Load the map data

bot_country <- st_read("../Data/Similar_phenology_invasive_code&data_20250218/geographic boundary data/bot_country_robin.shp") ## Transformed boundary of botanical country (Robinson projection)
bbox <- read_sf("../Data/Similar_phenology_invasive_code&data_20250218/geographic boundary data/worldmap_boundary_box.shp", layer = "worldmap_boundary_box") ## World map boundary
bbox_robin <- st_transform(bbox, crs=st_crs("+proj=robin"))

#Formatting the data

patern_result <- data.frame(diff = pattern$b) %>% 
  rownames_to_column("LEVEL3_COD")

patern_result$LEVEL3_COD <- substr(patern_result$LEVEL3_COD, 11, 13)
bot <- left_join(patern_result, bot_country, by = "LEVEL3_COD")
glimpse(bot)

map <- ggplot() + 
  ## Add the bounding box for the map area (background color: light cyan)
  geom_sf(data=bbox_robin, fill = "#E0FFFF") +
  ## Add the bot_country data with light gray fill to show the countries on the map
  geom_sf(data = bot_country, mapping = aes(geometry = `geometry`), fill = "#F7F7F7") + 
  ## Add the bot data with color scale showing the phenological difference (green to red)
  geom_sf(data = bot, mapping = aes(geometry = `geometry`, fill = diff)) +
  ## Apply a color gradient for the phenological difference
  scale_fill_gradient2(midpoint = 0,
                       low = "#FC8D59",  # Red for negative difference
                       mid = "#FFFFBF",  # Yellow for neutral (0 difference)
                       high = "#91CF60",  # Green for positive difference
                       limits = c(-4, 2.22499103),  # Set the scale range
                       space = "Lab",
                       name = "Phenological difference") +
  ## Customize the appearance of the map (no background, no axis labels)
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.position = "top",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
#Overall phenological differences across the botanical countries

map


#3.3.5 Flowering
#Run the model

flow_pattern <- rma.mv(yi, vi, mods = ~ LEVEL3_COD - 1,
                       random = list(~ 1 | n.species, 
                                     ~ 1 | e.species),
                       data = meta, subset = phenology == "Flowering")
summary(flow_pattern)


#Format the data

## Create a data frame with the results for flowering phenology
patern_flower <- data.frame(diff = flow_pattern$b) %>% 
  rownames_to_column("LEVEL3_COD")

## Extract the 3-character LEVEL3_COD from the full code
patern_flower$LEVEL3_COD <- substr(patern_flower$LEVEL3_COD, 11, 13)

## Merge the pattern results with the botanical country data
bot_flower <- left_join(patern_flower, bot_country, by = "LEVEL3_COD")
glimpse(bot_flower)

#Plot the flowering pattern

## Create a map for flowering phenology differences
map_flow <- ggplot() + 
  ## Add background bounding box layer
  geom_sf(data=bbox_robin, fill = "#E0FFFF") +
  
  ## Add botanical country boundary layer
  geom_sf(data = bot_country, mapping = aes(geometry = `geometry`), fill = "#F7F7F7") + 
  
  ## Add the flowering phenology pattern layer
  geom_sf(data = bot_flower, mapping = aes(geometry = `geometry`, fill = diff)) +
  
  ## Customize the color scale for phenological differences
  scale_fill_gradient2(midpoint = 0,
                       low = "#FC8D59",
                       mid = "#FFFFBF",
                       high = "#91CF60",
                       limits = c(-4.5, 2),
                       space = "Lab",
                       name = "Phenological difference") +
  
  ## Customize the plot theme for better presentation
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.position = "top",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
#Flowering differences across the botanical countries

map_flow

#3.3.6 Fruiting
#Run the model

frui_pattern <- rma.mv(yi, vi, mods = ~ LEVEL3_COD - 1,
                       random = list(~ 1 | n.species, 
                                     ~ 1 | e.species),
                       data = meta, subset = phenology == "Fruiting")
summary(frui_pattern)

#Format the results

## Create a data frame with the results for fruiting phenology
patern_fruit <- data.frame(diff = frui_pattern$b) %>% 
  rownames_to_column("LEVEL3_COD")

## Extract the 3-character LEVEL3_COD from the full code
patern_fruit$LEVEL3_COD <- substr(patern_fruit$LEVEL3_COD, 11, 13)

## Merge the pattern results with the botanical country data
bot_fruit <- left_join(patern_fruit, bot_country, by = "LEVEL3_COD")
glimpse(bot_fruit)

#Plot of fruiting difference pattern

## Create a map for fruiting phenology differences
map_frui <- ggplot() + 
  ## Add background bounding box layer
  geom_sf(data=bbox_robin, fill = "#E0FFFF") +
  
  ## Add botanical country boundary layer
  geom_sf(data = bot_country, mapping = aes(geometry = `geometry`), fill = "#F7F7F7") + 
  
  ## Add the fruiting phenology pattern layer
  geom_sf(data = bot_fruit, mapping = aes(geometry = `geometry`, fill = diff)) +
  
  ## Customize the color scale for phenological differences
  scale_fill_gradient2(
    midpoint = 0,
    low = "#FC8D59",
    mid = "#FFFFBF",
    high = "#91CF60",
    limits = c(-6, 6),
    space = "Lab",
    name = "Phenological difference"
  ) +
  
  ## Customize the plot theme for better presentation
  theme(
    panel.background = element_blank(),
    axis.title = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.position = "top",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )
#Fruiting difference pattern across botanical countries

map_frui

#3.4 Effects of botanical country area
#3.4.1 Calculate botanical country area
area_dat <- left_join(meta, 
                      bot %>% mutate(log_area_km2 = log2(as.numeric(st_area(st_sf(bot)))/1000000)),
                      by = "LEVEL3_COD")
#3.4.2 Flowering
#Run the model

area_flow <- rma.mv(yi,vi, mods = ~log_area_km2,
                    random = list(~ 1 | LEVEL3_COD / id,
                                  ~ 1 | n.species, 
                                  ~ 1 | e.species),
                    data = area_dat, subset = phenology == "Flowering")
summary(area_flow)

#Predictive values

area_flow_pred <- as.data.frame(predict(area_flow, 
                                        newmods = seq(min(area_dat[area_dat$phenology == "Flowering",
                                        ]$log_area_km2), 
                                        max(area_dat[area_dat$phenology == "Flowering",
                                        ]$log_area_km2),
                                        0.1), 
                                        addx = TRUE)) %>% 
  dplyr::select(X.log_area_km2, pred, ci.lb, ci.ub)
glimpse(area_flow_pred)

#Plot

area_flow_plot <- ggplot() + 
  geom_point(area_dat %>% select(phenology == "Flowering"), 
             mapping = aes(x = log_area_km2, y = yi, size = 1/vi),
             alpha = 0.5, shape = 21, color = "black", fill = "#FBB4AE")  +
  geom_ribbon(data = area_flow_pred, 
              aes(x = X.log_area_km2, y = pred, ymin = ci.lb, ymax = ci.ub), 
              fill = "gray", alpha = 0.5) +
  geom_line(data = area_flow_pred, 
            mapping = aes(x = X.log_area_km2, y = pred), 
            color = "black", alpha = 0.5, linetype = "dashed") + 
  scale_size(range = c(0.1, 5)) +
  scale_y_continuous(limits = c(-12, 7)) +
  scale_x_continuous(limits = c(6, 22), breaks = c(6, 11, 16, 22)) +
  theme_bw() + 
  ggtitle("Flowering") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        text = element_text(size = 18),
        axis.text= element_text(color = "black", size = 18),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black"),
        axis.text.y = element_text(angle = 90)) +
  labs(x = "Botanical country area (km^2)\n[log scale]", y = "Phenological difference") +
  guides(size = F)
#Effects of bountanical country area on flowring difference

area_flow_plot

#3.4.3 Fruiting
#Run the model

area_frui <- rma.mv(yi,vi, mods = ~log_area_km2,
                    random = list(~ 1 | LEVEL3_COD / id,
                                  ~ 1 | n.species, 
                                  ~ 1 | e.species),
                    data = area_dat, subset = phenology == "Fruiting")
summary(area_frui)

#Predictive values

area_frui_pred <- as.data.frame(predict(area_frui, 
                                        newmods = seq(min(area_dat[area_dat$phenology == "Fruiting",
                                        ]$log_area_km2), 
                                        max(area_dat[area_dat$phenology == "Fruiting",
                                        ]$log_area_km2),
                                        0.1), 
                                        addx = TRUE)) %>% 
  dplyr::select(X.log_area_km2, pred, ci.lb, ci.ub)
glimpse(area_frui_pred)

#Plot

area_frui_plot <- ggplot() + 
  geom_point(area_dat %>% select(phenology == "Fruiting"), 
             mapping = aes(x = log_area_km2, y = yi, size = 1/vi),
             alpha = 0.5, shape = 21, color = "black", fill = "#FFD92F")  +
  geom_ribbon(data = area_frui_pred, 
              aes(x = X.log_area_km2, y = pred, ymin = ci.lb, ymax = ci.ub), 
              fill = "gray", alpha = 0.5) +
  geom_line(data = area_frui_pred, 
            mapping = aes(x = X.log_area_km2, y = pred), 
            color = "black", alpha = 0.5, linetype = "dashed") + 
  scale_size(range = c(0.1, 5)) +
  theme_bw() + 
  scale_y_continuous(limits = c(-10, 6)) +
  scale_x_continuous(limits = c(8, 23), breaks = c(8, 13, 18, 23)) +
  ggtitle("Fruiting") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        text = element_text(size = 18),
        axis.text= element_text(color = "black", size = 18),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black"),
        axis.text.y = element_text(angle = 90)) +
  labs(x = "Botanical country area (km^2)\n[log scale]", y = "Phenological difference") +
  guides(size = F)
#Effects of bountanical country area on fruiting difference

area_frui_plot

#3.4.4 Result tabel
#Combine the results

area_result <- data.frame(Phenology = c("Flowering", "Fruiting"),
                          Slope = c(area_flow$b[2], area_frui$b[2]),
                          Lower_95_CI = c(area_flow$ci.lb[2], area_frui$ci.lb[2]),
                          Upper_95_CI = c(area_flow$ci.ub[2], area_frui$ci.ub[2]),
                          P_value = c(area_flow$pval[2], area_frui$pval[2])) 
#Phenological difference across botanical country area

kable(area_result,
      digits = 3,align = "c") %>%
  kable_styling(bootstrap_options = "bordered", font_size = 15,latex_options = "scale_down") %>%
  row_spec(0, background = "gray", color = "white") %>% 
  column_spec(1, background = "gray", color = "white")

#3.5 Latitudinal patterns
#3.5.1 Flowering
#Run the model

lat_form_flow <- rma.mv(yi, vi, mods = ~latitude.abs + growthform + 
                          latitude.abs: growthform -1 -latitude.abs,
                        random = list(~ 1 | LEVEL3_COD / id,
                                      ~ 1 | n.species, 
                                      ~ 1 | e.species),
                        data = meta, subset = phenology == "Flowering")
summary(lat_form_flow)

#Predictive values

form_flow_pred <- as.data.frame(predict(lat_form_flow, 
                                        newdata = meta[meta$phenology == "Flowering",], 
                                        addx = TRUE))

glimpse(form_flow_pred)

#Plot

lat_flow_pattern <- ggplot() + 
  geom_point(meta[meta$phenology == "Flowering",], 
             mapping = aes(x = latitude.abs, y = yi, size = 1/vi, 
                           fill = growthform, shape = growthform, color = growthform), 
             alpha = 0.5)  +
  geom_ribbon(data = form_flow_pred[form_flow_pred$X.growthformWoody == 1,], 
              aes(x = X.latitude.abs.growthformWoody, y = pred, ymin = ci.lb, ymax = ci.ub), 
              fill = "gray", alpha = 0.5) +
  geom_line(data = form_flow_pred[form_flow_pred$X.growthformWoody == 1,], 
            mapping = aes(x = X.latitude.abs.growthformWoody, y = pred), 
            size = 1, color = "#A6761D", linetype = "dashed") +
  geom_ribbon(data = form_flow_pred[form_flow_pred$X.growthformHerbaceous == 1, ], 
              aes(x = X.latitude.abs.growthformHerbaceous, y = pred,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "gray", alpha = 0.5) +
  geom_line(data = form_flow_pred[form_flow_pred$X.growthformHerbaceous == 1, ], 
            mapping = aes(x = X.latitude.abs.growthformHerbaceous, y = pred), 
            size = 1, color = "#66A61E") +
  scale_size(range = c(0.1, 5)) +
  scale_x_continuous(limits = c(0, 60), breaks = seq(0,60,20),
                     name = "Absolute latitude") +
  scale_y_continuous(limits = c(min(meta[meta$phenology == "Flowering",]$yi),
                                max(meta[meta$phenology == "Flowering",]$yi)),
                     name = "Phenological difference") +
  scale_fill_manual(values = c("Herbaceous" = "#66A61E", "Woody" = "#A6761D"), 
                    guide = guide_legend(override.aes = list(size = 3), reverse = F)) +
  scale_color_manual(values = c("Herbaceous" = "#66A61E", "Woody" = "#A6761D"), 
                     guide = guide_legend(override.aes = list(size = 3), reverse = F)) +
  scale_shape_manual(values = c("Herbaceous" = 23, "Woody" = 21), 
                     guide = guide_legend(override.aes = list(size = 3), reverse = F)) +
  theme_bw() +
  ggtitle("Flowering") + 
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 18, color = "black"),
        panel.grid = element_blank(),
        axis.title = element_text(size = 18, color = "black"),
        panel.border = element_rect(color = "black"),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.125),
        legend.background = element_blank()) +
  guides(size = F)
#Latitudinal patterns of flowering differences in herbaceous and woody plants

lat_flow_pattern

#3.5.2 Fruiting
#Run the model

lat_form_frui <- rma.mv(yi, vi, mods = ~latitude.abs + growthform + 
                          latitude.abs : growthform -1 -latitude.abs,
                        random = list(~ 1 | LEVEL3_COD / id,
                                      ~ 1 | n.species, 
                                      ~ 1 | e.species),
                        data = meta, subset = phenology == "Fruiting")
summary(lat_form_frui)

#Predictive values
form_frui_pred <- as.data.frame(predict(lat_form_frui, 
                                        newdata = meta[meta$phenology == "Fruiting",], 
                                        addx = TRUE))
glimpse(form_frui_pred)

#Plot

lat_frui_pattern <- ggplot() + 
  geom_point(meta[meta$phenology == "Fruiting",], 
             mapping = aes(x = latitude.abs, y = yi, size = 1/vi, 
                           fill = growthform, shape = growthform, color = growthform), 
             alpha = 0.5)  +
  geom_ribbon(data = form_frui_pred[form_frui_pred$X.growthformWoody == 1,], 
              aes(x = X.latitude.abs.growthformWoody, y = pred, ymin = ci.lb, ymax = ci.ub), 
              fill = "gray", alpha = 0.5) +
  geom_line(data = form_frui_pred[form_frui_pred$X.growthformWoody == 1,], 
            mapping = aes(x = X.latitude.abs.growthformWoody, y = pred), 
            size = 1, color = "#A6761D", linetype = "dashed") +
  geom_ribbon(data = form_frui_pred[form_frui_pred$X.growthformHerbaceous == 1, ], 
              aes(x = X.latitude.abs.growthformHerbaceous, y = pred, ymin = ci.lb, ymax = ci.ub),
              fill = "gray", alpha = 0.5) +
  geom_line(data = form_frui_pred[form_frui_pred$X.growthformHerbaceous == 1, ], 
            mapping = aes(x = X.latitude.abs.growthformHerbaceous, y = pred),
            size = 1, color = "#66A61E") +
  scale_size(range = c(0.1, 5)) + 
  scale_y_continuous(limits=c(min(meta[meta$phenology == "Flowering",]$yi),
                              max(meta[meta$phenology == "Flowering",]$yi)),
                     name = "Phenological difference") +
  scale_x_continuous(limits=c(0, 60), breaks = seq(0, 60, 20),
                     name = "Absolute latitude") +
  scale_fill_manual(values = c("Herbaceous" = "#66A61E", "Woody" = "#A6761D"), 
                    guide = guide_legend(override.aes = list(size = 3), reverse = F)) + 
  scale_color_manual(values = c("Herbaceous" = "#66A61E", "Woody" = "#A6761D"), 
                     guide = guide_legend(override.aes = list(size = 3), reverse = F)) +
  scale_shape_manual(values = c("Herbaceous" = 23, "Woody" = 21), 
                     guide = guide_legend(override.aes = list(size = 3), reverse = F)) +
  theme_bw() + 
  ggtitle("Fruiting") + 
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        panel.border = element_rect(color = "black"),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.125),
        legend.background = element_blank()) +
  guides(size = F)
#Latitudinal patterns of flowering differences in herbaceous and woody plants

lat_frui_pattern

#3.5.3 Result table
#Combine the results

lat_form <- data.frame(phenology = rep(c("Flowering", "Fruiting"), c(2, 2)), 
                       growthform = rep(c("Herbaceous", "Woody"), 2),
                       slope = c(lat_form_flow$b[3], lat_form_flow$b[4],
                                 lat_form_frui$b[3], lat_form_frui$b[4]),
                       lower_95_CI = c(lat_form_flow$ci.lb[3], lat_form_flow$ci.lb[4],
                                       lat_form_frui$ci.lb[3], lat_form_frui$ci.lb[4]),
                       upper_95_CI = c(lat_form_flow$ci.ub[3], lat_form_flow$ci.ub[4],
                                       lat_form_frui$ci.ub[3], lat_form_frui$ci.ub[4]),
                       test_statistic = c(lat_form_flow$pval[3], lat_form_flow$pval[4],
                                          lat_form_frui$pval[3], lat_form_frui$pval[4]))  # Create final dataframe with results for both phenologies
#Latitudinal patterns of phenological differences in herbaceous and woody plants

kable(lat_form,
      digits = 3, align = "c") %>%
  kable_styling(bootstrap_options = "bordered", 
                font_size = 15, latex_options = "scale_down") %>%
  row_spec(0, background = "gray", color = "white") %>% 
  column_spec(1:2, background = "gray", color = "white") %>%
  row_spec(c(1, 3), bold = TRUE)

#3.6 Continental patterns
#3.6.1 Flowering
#Run the model

cont_flow <- rma.mv(yi,vi, mods = ~continent-1,
                    random = list(~ 1 | n.species, 
                                  ~ 1 | e.species),
                    data = meta, subset = phenology == "Flowering")
summary(cont_flow)

#Format the results

cont_flow_result <- data.frame(
  continent = Replace(row.names(as.data.frame(cont_flow$b)), "continent",""),
  estimate = cont_flow$b,
  lower_95_CI = cont_flow$ci.lb,
  upper_95_CI = cont_flow$ci.ub,
  test_statistic = cont_flow$pval,
  row.names = NULL) %>% 
  mutate(continent = factor(continent, levels = c("Europe", "Africa","Asia-Temperate", 
                                                  "Asia-Tropical","Australasia",
                                                  "Pacific","Northern America",
                                                  "Southern America"
  )))
kable(cont_flow_result, digits = 3, align = "c") %>%
  kable_styling(bootstrap_options = "bordered", font_size = 15, latex_options = "scale_down") %>%
  row_spec(0, background = "gray", color = "white") %>%
  column_spec(1, background = "gray", color = "white")

#Plot

cont_flow_plot <- ggplot() + 
  geom_hline(yintercept = 0, linetype = 2) +
  geom_quasirandom(meta %>% 
                     mutate(continent = factor(continent, 
                                               levels = c("Europe", "Africa","Asia-Temperate", 
                                                          "Asia-Tropical","Australasia",
                                                          "Pacific","Northern America",
                                                          "Southern America"
                                               ))) %>% 
                     filter(phenology == "Flowering"),
                   mapping = aes(x = continent, y = yi), fill = "#F7F7F7",
                   alpha = 0.5) + 
  geom_pointrange(data = cont_flow_result, 
                  mapping = aes(x = continent, y = estimate, ymin = lower_95_CI, 
                                ymax = upper_95_CI), 
                  shape = 21, color = "black",size = 0.75, fill = "white") +
  scale_y_continuous(name = "Phenological difference", limits = c(-10, 10)) + 
  theme_classic() + 
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 18, color = "black"),
        axis.title.x.bottom = element_text(size = 18),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank()) +
  coord_flip() + guides(color = F, fill = F, shape = F)

#Flowering difference patterns across continents

cont_flow_plot

#3.6.2 Fruiting
#Run the model

cont_frui <- rma.mv(yi,vi, mods = ~continent-1,
                    random = list(~ 1 | n.species, 
                                  ~ 1 | e.species),
                    data = meta, subset = phenology == "Fruiting")
summary(cont_frui)

#Format the results

cont_frui_result <- data.frame(
  continent = Replace(row.names(as.data.frame(cont_frui$b)), "continent",""),
  estimate = cont_frui$b,
  lower_95_CI = cont_frui$ci.lb,
  upper_95_CI = cont_frui$ci.ub,
  test_statistic = cont_frui$pval,
  row.names = NULL) %>% 
  mutate(continent = factor(continent, levels = c("Europe", "Africa","Asia-Temperate", 
                                                  "Asia-Tropical","Australasia",
                                                  "Pacific","Northern America",
                                                  "Southern America"
  )))
kable(cont_frui_result, digits = 3, align = "c") %>%
  kable_styling(bootstrap_options = "bordered", font_size = 15, latex_options = "scale_down") %>%
  row_spec(0, background = "gray", color = "white") %>%
  column_spec(1, background = "gray", color = "white")


#Plot

cont_frui_plot <- ggplot() + 
  geom_hline(yintercept = 0, linetype = 2) +
  geom_quasirandom(meta %>% 
                     mutate(continent = factor(continent, levels = 
                                                 c("Europe", "Africa","Asia-Temperate", 
                                                   "Asia-Tropical","Australasia",
                                                   "Pacific","Northern America",
                                                   "Southern America"
                                                 ))) %>% 
                     filter(phenology == "Fruiting"),
                   mapping = aes(x = continent, y = yi),fill = "#F7F7F7",
                   alpha = 0.5) + 
  geom_pointrange(data = cont_frui_result, 
                  mapping = aes(x = continent, y = estimate, ymin = lower_95_CI, 
                                ymax = upper_95_CI), 
                  shape = 21, color = "black",size = 0.75, fill = "white") +
  scale_y_continuous(name = "Phenological difference", limits = c(-10, 10)) + 
  theme_classic() + 
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 18, color = "black"),
        axis.title.x.bottom = element_text(size = 18),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank()) +
  coord_flip() + guides(color = F, fill = F, shape = F)
#Fruiting difference patterns across continents

cont_frui_plot

#3.7 Effects of growth forms
#3.7.1 FLowering
#Run the model

form_flow <- rma.mv(yi, vi, mods = ~growthform - 1,
                    random = list(~ 1 | LEVEL3_COD / id,
                                  ~ 1 | n.species, 
                                  ~ 1 | e.species),
                    data = meta, subset = phenology == "Flowering")
summary(form_flow)
#ANOVA test

anova(form_flow, L = c(1, -1))
#3.7.2 Fruiting
#Run the model

form_frui <- rma.mv(yi, vi, mods = ~growthform - 1,
                    random = list(~ 1 | LEVEL3_COD / id,
                                  ~ 1 | n.species, 
                                  ~ 1 | e.species),
                    data = meta, subset = phenology == "Fruiting")
summary(form_frui)
#ANOVA test

anova(form_frui, L = c(1, -1))
#3.7.3 Visulise Phenological difference across different growth forms
#Format the results

form_result <- data.frame(
  phenology = rep(c("Flowering", "Fruiting"), c(2, 2)),
  growthform = rep(c("Herbaceous", "Woody"), 2),
  estimate = c(form_flow$b, form_frui$b),
  lower_95_CI = c(form_flow$ci.lb, form_frui$ci.lb),
  upper_95_CI = c(form_flow$ci.ub, form_frui$ci.ub),
  test_statistic = c(form_flow$pval, form_frui$pval),
  row.names = NULL
) %>% 
  mutate(phenology = factor(phenology, levels = rev(phenology_level))) %>% 
  mutate(growthform = factor(growthform, levels = c("Woody", "Herbaceous")))
#Phenological difference across different growth forms

kable(form_result, 
      digits = 3, align = "c") %>%
  kable_styling(bootstrap_options = "bordered", 
                font_size = 15, latex_options = "scale_down") %>%
  row_spec(0, background = "gray", color = "white") %>% 
  column_spec(1:2, background = "gray", color = "white") %>%
  row_spec(c(1, 3), bold = TRUE)

#Plot

form <- ggplot(form_result, mapping = aes(x = estimate, y = phenology)) + 
  geom_vline(xintercept = 0, 
             color = 'black', 
             linetype = 'longdash') +  # Add a vertical line at x = 0 as a reference
  geom_pointrange(aes(xmin = lower_95_CI, xmax = upper_95_CI, 
                      fill = growthform, shape = growthform), 
                  position = position_dodge(width = 0.8), color = "black", size = 1) +  # Add point-range plot showing estimates and 95% CI
  scale_x_continuous(limits = c(-3, 3), 
                     breaks = seq(-4, 4, 2),
                     name = "Phenological difference", 
                     expand = c(0, 0)) +  # Set x-axis limits, breaks, and label
  scale_fill_manual(values = c("Herbaceous" = "#66A61E", "Woody" = "#A6761D"), 
                    guide = guide_legend(reverse = TRUE)) +  # Set custom colors for the fill (growthform)
  scale_shape_manual(values = c("Herbaceous" = 23, "Woody" = 21), 
                     guide = guide_legend(reverse = TRUE)) +  # Set custom shapes for the points based on growthform
  theme_classic() +  # Apply classic theme for a clean look
  theme(axis.text.x = element_text(size = 18, color = "black"), 
        axis.title.x = element_text(size = 18),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 18, color = "black", angle = 90, hjust = 0.5),
        axis.text = element_text(color = "black"),
        legend.justification = "top",
        legend.title = element_blank(),
        legend.position = c(0.15, 1.06),
        legend.background = element_rect(fill = "transparent")) +  # Customize text and legend styles
  annotate("segment", x = -3, xend = 3, y = 0.4, yend = 0.4, linewidth = 1)  # Add a horizontal line to separate the groups visually
#Flowering and fruiting differences in herbaceous and woody plants.

form

#3.8 Effects of phylogenetic distance
#3.8.1 Establish phylogenetic tree
#Species checklist analyzed in the study

iNa_species <- unique(data.frame(scientific_name = unique(flower$scientific_name)) %>% 
                        rbind(data.frame(scientific_name = unique(fruit$scientific_name))))
iNa_species <- left_join(iNa_species, plant, by = "scientific_name")
#Get the genu list

iNa_genuslist <- unique(iNa_species[2:3])
#Load the mega tree from Zanne et al., Nature: 506, 89–92 (2014)

zanne <- read.tree("../Data/Similar_phenology_invasive_code&data_20250218/Zanne_mega_tree.tre")
#Construct the phylogenetic tree

tree <- phylo.maker(iNa_species, zanne, iNa_genuslist, 
                    nodes.type = 1, scenario = 3)
tree <- compute.brlen(tree$phylo)
#3.8.2 Flowering
#Getting the phylogenetic distance

flower_dis <- meta[meta$phenology == "Flowering",] %>% dplyr::select(-n)

for (i in 1:length(flower_dis$LEVEL3_COD)) {
  flower_dis$dis[i] <- get_pairwise_distances(
    tree, 
    Replace(flower_dis$e.species[i], " ", "_"),
    Replace(flower_dis$n.species[i], " ", "_")
  )
}
#Run the model

phylo_flow <- rma.mv(yi, vi, mods = ~dis,
                     random = list(~ 1 | LEVEL3_COD / id,
                                   ~ 1 | n.species, 
                                   ~ 1 | e.species),
                     data = flower_dis)
summary(phylo_flow)
#Predictive values

phylo_flow_pred <- as.data.frame(predict(
  phylo_flow, newmods = seq(min(flower_dis$dis), 
                            max(flower_dis$dis), 0.001), 
  addx = TRUE)) %>% 
  dplyr::select(X.dis, pred, ci.lb, ci.ub)

glimpse(phylo_flow_pred)

#Plot

phylo_flow_plot <- ggplot() + 
  geom_point(flower_dis, mapping = aes(x = dis, y = yi, size = 1/vi),
             alpha = 0.5, shape = 21, color = "black", fill = "#FBB4AE")  +
  geom_ribbon(data = phylo_flow_pred, 
              aes(x = X.dis, y = pred, ymin = ci.lb, ymax = ci.ub), 
              fill = "gray", alpha = 0.5) +
  geom_line(data = phylo_flow_pred, 
            mapping = aes(x = X.dis, y = pred), 
            color = "black", alpha = 0.5, linetype = "dashed") + 
  scale_y_continuous( limits = c(min(flower_dis$yi), 
                                 max(flower_dis$yi)),
                      name = "Phenological difference") + 
  scale_x_continuous( name = "Phylogenetic distance") + 
  scale_size(range = c(0.1, 5)) +
  theme_bw() + 
  ggtitle("Flowering") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        text = element_text(size = 18),
        axis.text= element_text(color = "black", size = 18),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black")) +
  guides(size = F)
#Flowering differences across phylogenetic distance

phylo_flow_plot

#3.8.3 Fruiting
#Getting the phylogenetic distance

fruit_dis <- meta[meta$phenology == "Fruiting",] %>% 
  dplyr::select(-n)

for (i in 1:length(fruit_dis$LEVEL3_COD)) {
  fruit_dis$dis[i] <- get_pairwise_distances(
    tree, 
    Replace(fruit_dis$e.species[i], " ", "_"),
    Replace(fruit_dis$n.species[i], " ", "_")
  )
}
#Run the model

phylo_frui <- rma.mv(yi, vi, mods = ~dis,
                     random = list(~ 1 | LEVEL3_COD / id,
                                   ~ 1 | n.species, 
                                   ~ 1 | e.species),
                     data = fruit_dis)
summary(phylo_frui)

#Predictive values

phylo_frui_pred <- as.data.frame(predict(
  phylo_frui, newmods = seq(min(fruit_dis$dis), 
                            max(fruit_dis$dis), 0.001), 
  addx = TRUE)) %>% 
  dplyr::select(X.dis, pred, ci.lb, ci.ub)
glimpse(phylo_frui_pred)

#Plot

phylo_frui_plot <- ggplot() + 
  geom_point(fruit_dis, mapping = aes(x = dis, y = yi, size = 1/vi),
             alpha = 0.5, shape = 21, color = "black", fill = "#FFD92F")  +
  geom_ribbon(data = phylo_frui_pred, 
              aes(x = X.dis, y = pred, ymin = ci.lb, ymax = ci.ub), 
              fill = "gray", alpha = 0.5) +
  geom_line(data = phylo_frui_pred, 
            mapping = aes(x = X.dis, y = pred), 
            color = "black", alpha = 0.5, linetype = "dashed") + 
  scale_y_continuous(limits = c(min(fruit_dis$yi), 
                                max(fruit_dis$yi)),
                     name = "Phenological difference") + 
  scale_x_continuous( name = "Phylogenetic distance") + 
  scale_size(range = c(0.1, 5)) +
  theme_bw() + 
  ggtitle("Fruiting") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        text = element_text(size = 18),
        axis.text= element_text(color = "black", size = 18),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black")) +
  guides(size = F)
#Fruiting differences across phylogenetic distance

phylo_frui_plot

#3.8.4 Result table
#Combine the results

phylo_result <- data.frame(
  phenology = c("Flowering", "Fruiting"),
  slope = c(phylo_flow$b[2], phylo_frui$b[2]),
  lower_95_CI = c(phylo_flow$ci.lb[2], phylo_frui$ci.lb[2]),
  upper_95_CI = c(phylo_flow$ci.ub[2], phylo_frui$ci.ub[2]),
  test_statistic = c(phylo_flow$pval[2], phylo_frui$pval[2])
)

#Phenology difference across phylogenetic distance

kable(phylo_result,
      digits = 3, align = "c") %>%
  kable_styling(bootstrap_options = "bordered", font_size = 15, latex_options = "scale_down") %>%
  row_spec(0, background = "gray", color = "white") %>% 
  column_spec(1, background = "gray", color = "white")
