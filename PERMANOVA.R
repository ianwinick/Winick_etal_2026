library(tidyverse)
library(vegan)
library(FD)
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(ggordiplots)

# COMBINE AND CLEAN COMMUNITY MATRICES

# Read in cover for each year, add year column
cov.20<-read_csv("CoverSep2020.csv")
cov.20$Year<-2020

cov.21<-read_csv("CoverSep2021.csv")
cov.21$Year<-2021

cov.22<-read_csv("CoverSep2022.csv")
cov.22$Year<-2022

cov.23<-read_csv("CoverSep2023.csv")
cov.23$Year<-2023

cov.24<-read_csv("CoverSep2024.csv")
cov.24$Year<-2024

#combine for total cover across years and treatments
library(plyr) 
cov.CLIFF<-rbind.fill(cov.20, cov.21, cov.22, cov.23, cov.24)
detach("package:plyr", unload = TRUE) ##### included this because i (ian) have had issues with plyr functions conflicting with functions in other packages 
cov.CLIFF$Year<-as.factor(cov.CLIFF$Year)
cov.CLIFF$plot<-as.factor(cov.CLIFF$plot)
cov.CLIFF$severity<-as.factor(cov.CLIFF$severity)
cov.CLIFF[is.na(cov.CLIFF)] <- 0 #Make NAs 0

year <- cov.CLIFF$Year
plot <- cov.CLIFF$plot
severity <- cov.CLIFF$severity

#read in dominant (i.e., not rare) species
traits <- read_csv("TraitTable.csv") %>% 
  mutate(spp = case_when(spp == "MUVI" ~ "MUST", .default = spp)) %>% 
  arrange(spp) %>%
  column_to_rownames(var="spp") %>%
  select(sla, height, seedmass, resprouting) %>%
  mutate(resprouting=resprouting+1)

spp_pool <- rownames(traits)
cov.CLIFF <- cov.CLIFF[,(which(names(cov.CLIFF) %in% spp_pool))]
cov.CLIFF <- cov.CLIFF[,order(names(cov.CLIFF))]
cov.CLIFF$Year <- year ; cov.CLIFF$plot <- plot 
cov.CLIFF$severity <- severity

# For loop that relativizes cover across severity classes within years
year <- c(2020, 2021, 2022, 2023, 2024)
cov.stand <- NULL
for(i in year){
  cov.CLIFF_wisoncsin <- cov.CLIFF %>%
    filter(Year==i)
  wisconsin <- cov.CLIFF_wisoncsin %>%
    select(-c(plot, Year, severity)) %>%
    wisconsin()
  cov.stand<-rbind(cov.stand, wisconsin)
}

# REPEATED MEASURES PERMANOVA ##################################################

#This chunk of code controls how permutations are carried out
CTRL.t <- how(within = Within(type = "free"),
              plots = Plots(type = "none"),
              blocks=cov.CLIFF$plot,
              nperm = 999,
              observed = TRUE)


adonis.out<-adonis2(cov.stand~plot+severity*Year,
                    data=cov.CLIFF,
                    method="bray",
                    permutations=CTRL.t,
                    by="margin")

adonis.out #interaction is significant, create new factor for pairwise adonis

# Pairwise repeated-measures PERMANOVA
cov.CLIFF$sev.year<-paste(cov.CLIFF$severity,cov.CLIFF$Year)

pairwise.out<-pairwise.adonis2(cov.stand~sev.year,
                               data=cov.CLIFF,
                               method="bray",
                               by="margin")
pairwise.out #multiple significant pairwise comparisons

# REPEATED MEASURES TEST FOR BETA DISPERSION ###################################

permutest(betadisper(vegdist(cov.stand, method="bray"), severity, type="centroid"),pairwise=T, permutations=CTRL.t)

# ENVFIT ANALYSIS AND NMDS PLOTS WITH TRAIT VECTORS ############################

fd.out<-dbFD(traits, cov.stand)#run dbFD to get CWM
cwm.traits<-fd.out$CWM#isolate CWM, inspect
cwm.traits$resprouting<-cwm.traits$resprouting-1#subtract 1 from resprouting to make it binary

nmds <- metaMDS(cov.stand, distance="bray", k=2, trymax=100)
nmds$stress

envfit_results <- envfit(nmds, cwm.traits)

env <- data.frame(envfit(nmds_5yr, cwm.traits)$vectors$arrows, 
                 traits = c("sla", "height", "seedmass", "resprouting"))


env_0 <- data.frame(NMDS1 = c(0, 0, 0, 0), NMDS2 = c(0, 0, 0, 0), 
                      traits = c("sla", "height", "seedmass", "resprouting"))
env_lines <- rbind(env_5yr, env_5yr_0)

data_nmds <- as_tibble(scores(nmds_5yr, display="sites"))

plot_data <- select(cov.CLIFF, Year, severity)
tax_scores <- cbind(data_nmds, plot_data)

tax_means <- tax_scores %>% 
  group_by(Year, severity) %>% 
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))

severity_colors_points <- c("U" = "#5bbcd695", "L" = "#f9840295", "H" = "#fb040495")
severity_colors_centroids <- c("U" = "#5bbcd6", "L" = "#f98402", "H" = "#fb0404")

ggplot(tax_scores, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(color = severity), data = tax_scores, alpha = 0.5)+
  scale_color_manual(values = severity_colors_points)+
  scale_fill_manual(values = severity_colors_points)+
  geom_point(aes(fill = severity), size = 3, alpha = 1,
             data = tax_means, pch = 21)+
  stat_ellipse(aes(color = severity, linetype = severity), data = tax_scores,
               linewidth = 1)+
  geom_line(data = env_lines, aes(group = traits))+
  geom_label(data = env_5yr, aes(label = traits))+ # use this to figure out which env line is which
  facet_wrap(~Year, ncol = 1)+
  theme_minimal(base_size = 20)+
  labs(fill = "Severity", color = "Severity", linetype = "Severity")
 ggsave("outputs/5yrNMDS.png", width = 8, height = 20, units = "in", dpi = 600)
