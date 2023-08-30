################################################################
# dose effect FvFm
# Create a growth curve at the selectesd day for DMSO vs AZ data
################################################################

# Prevent from scientific notation use
options(scipen=999)

#### Part to modify #######################
#===============================================================================
# Variables 
WORKDIR <- "~/partage/bio-informatique/figures_papier_romain"

Data <- "data_FvFm_dose_effect_AZD.txt"


# Number of day
NB_DAY <- 5

# Day of interest
DAY <- 4

# Reference line
RefLine <- "WT"


# Target order for order on the graph
target_order <- c("WT", "tor-15")
az_order <- c("AZD.0.003.µM", "AZD.0.01.µM", "AZD.0.03.µM", "AZD.0.1.µM",  "AZD.0.3.µM", "AZD.1.µM", "AZD.3.µM", "AZD.10.µM")

# Variables de personnalisation du graphique pour la courbe de croissance

TITLE <- ""
Y_AXIS <- ""


# Variables de personnalisation du graphique pour l'analyse statistique

# Définir la couleur des points
COULEUR <- "#dbdad7"



if (!require(ggplot2)) { install.packages("ggplot2") }
library(ggplot2)
# Orientationd des étiquettes de l'axe X.
# Pour des étiquettes horizontales : angle = 0, vjust = 0, hjust=0.
# Pour des étiquettes tournées de 90° : angle = 90, vjust = 0.5, hjust=1
orientation_xlabels <- theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5))

# Position and appearance of group markers
# The default values ​​produce a graph where the parts of the graph are separated with the factor
# of grouping in a frame inside the graph.
# You can modify it like this:
# joining parts: panel.spacing = unit(0, "lines")
# No frame around grouping factor: strip.background = element_rect(colour=NA, fill=NA
# grouping factor outside the graph :strip.placement = "outside"
strip_pos <- theme(panel.spacing = unit(0.3, "lines"), 
                   strip.background = element_rect(colour="black", fill=NA),
                   strip.placement = "inside")

################################################################################################################################

# Fonction utilisée dans le script

#=======================================================================================================================================
# Vérifie la normalité des données. Sort de la boucle si au moins un des groupes de données ne suit 
# pas une loi normale.
#Retourne TRUE si les données suivent une loi normale
#=======================================================================================================================================
check_normality <- function(shapiro_df) {
  # On suppose que les données sont normales
  flag_normal <- TRUE
  
  for (i in 1 : nrow(shapiro_df)) {
    if(shapiro_df[i, 4] > 0.05) {
      # print(paste0("les données ",shapiro_df$grouping_factor[i],"-", 
      # shapiro_df$plant_line[i], " suivent une loi normale"), quote = FALSE)
      
    } else {
      # print(paste0("les données ",shapiro_df$grouping_factor[i],"-", 
      #          shapiro_df$plant_line[i], " ne suivent pas une loi normale"), quote = FALSE)
      
      # En fait les données ne sont pas normales, pas besoin d'aller plus loin
      flag_normal <- FALSE
      break
    }
  }
  return(flag_normal)
}

#=======================================================================================================================================
# Prend les résultats du test Anova
# Retourne TRUE si au moins une moyenne est significativement différente
#=======================================================================================================================================
check_anova <- function(anova_results) {
  flag_anova <- FALSE
  
  for (i in 1 : nrow(anova_results)) {
    if (anova_results$p[i] < 0.05) {
      print (paste0("Le test Anova compare les moyennes, la pvalue pour le jour ", anova_results$traitement[i] , " est < 0.05 ce qui indique qu’au moins 1 des moyennes est différentes des autres, on réalise un test post hoc de Tukey"))
      flag_anova <- TRUE
    } else {
      print (paste0("Le test Anova compare les moyennes, la pvalue pour le jour ", anova_results$traitement[i] , " est > 0.05 ce qui indique qu’il n'y a pas de différence entre les moyennes"))
    }
  }    
  
  return(flag_anova)
}

#=======================================================================================================================================
# Fait le test de Kruskal-Wallis
# Retourne TRUE si au moins une médianes est significativement différente
#=======================================================================================================================================
check_kruskal <- function(kruskal_pval) {
  flag_kruskal <- FALSE
  
  for (i in 1 : nrow(kruskal_pval)) {
    if (kruskal_pval$p[i] < 0.05) {
      print (paste0("Le test de Kruskall Wallis compare les médianes, la pvalue pour le groupe ", kruskal_pval$traitement[i] , " est < 0.05 ce qui indique qu’au moins 1 des médianes est différentes des autres, on réalise un test post hoc de Dunn"))
      flag_kruskal <- TRUE
    } else {
      print (paste0("Le test de Kruskall Wallis compare les médianes, la pvalue pour le groupe ", kruskal_pval$traitement[i] , " est > 0.05 ce qui indique qu’il n'y a pas de différence entre les médianes"))
    }
  }    
  
  return(flag_kruskal)
}

#=======================================================================================================================================
# Fait le test de Dunn
# retourne les pvalue dans un dataframe
#=======================================================================================================================================
test_dunn <- function(df_data) {
  pval <- as.data.frame(df_data%>%
                          mutate(lines = fct_relevel(lines, target_order)) %>% 
                          group_by(traitement) %>% dunn_test(DMSO_percent ~ lines, p.adjust.method = "BH"))
  #print(df_data %>% group_by(Day) %>% dunn_test(Length ~ line, p.adjust.method = "BH"))
  return(pval)
}

#=======================================================================================================================================
# Réalise le graphique lorsque les données suivent une loi normale. L'intervalle de confiance est 
# construit autour de la moyenne
#=======================================================================================================================================
plot_normal <- function(df, my_colours, my_summary) {
  p <- df %>%
    mutate(lines = fct_relevel(lines, 
                               target_order)) %>%
    ggplot( aes(x=lines, y=DMSO_percent)) +
    geom_quasirandom(dodge.width=0.8,alpha = 0.6, colour=COULEUR) +
    geom_pointrange(data=my_summary%>% mutate(lines = fct_relevel(lines,  target_order)), 
                    aes(ymin=DMSO_percent - ci, ymax=DMSO_percent + ci),
                    position=position_dodge(width=0.8)) + 
    scale_colour_manual(values="black") +
    scale_x_discrete(labels = c("WT", expression(italic("tor-15")))) +
    scale_y_continuous(name = "FvFm",
                       expand = expansion(mult = c(.1, .1))) +    # échelle augmentée de 10% au dessus du point le plus haut et en dessous du point le plus bas
    facet_wrap(~ traitement, strip.position = "bottom", scales = "free") +
    theme_classic() + 
    strip_pos +
    orientation_xlabels +
    theme(legend.title=element_blank())
  print(p)
}


#=======================================================================================================================================
# Réalise le graphique lorsque les données ne suivent pas une loi normale. L'intervalle de confiance est 
# construit autour de la médiane
#=======================================================================================================================================
plot_not_normal <- function(df, my_colours, conf_int) {
  # La colonne "median" doit porter le nom "mesured_value" pour le ggplot
  names(conf_int)[4] <- "DMSO_percent"
  
  p <- df %>%
    mutate(lines = fct_relevel(lines, 
                               target_order)) %>%
    ggplot( aes(x=lines, y=DMSO_percent)) +
    geom_quasirandom(dodge.width=0.8, alpha = 0.6, colour=COULEUR) +
    geom_pointrange(data=conf_int%>% mutate(lines = fct_relevel(lines,  target_order)), 
                    aes(ymin=Percentile.lower, ymax=Percentile.upper), 
                    position=position_dodge(width=0.8)) +
    scale_colour_manual(values="black") +
    scale_x_discrete(labels = c("WT", expression(italic("tor-15")))) +
    scale_y_continuous(name = "FvFm",
                       expand = expansion(mult = c(.1, .1))) + # échelle augmentée de 10% au dessus du point le plus haut et en dessous du point le plus bas
    facet_wrap(~ traitement, strip.position = "bottom", scales = "free") +
    theme_classic() +  
    strip_pos +
    orientation_xlabels +
    theme(legend.title=element_blank())
  print(p)
}
################################################################################################################################

# Installation des pacakges

if (!require(tidyr)) { install.packages("tidyr") }
if (!require(dplyr)) { install.packages("dplyr") }
if (!require(tidyverse)) { install.packages("tidyverse") }
if (!require(devtools)) { install.packages("devtools") }
if (!require(rstatix)) { install.packages("rstatix", repos = "https://cloud.r-project.org") }
if (!require(ggbeeswarm)) { install.packages("ggbeeswarm") }
if (!require(RColorBrewer)) { install.packages("RColorBrewer") }
if (!require(rcompanion)) { install.packages("rcompanion") }
if (!require(scales)) {install.packages("scales")}
if (!require(ggExtra)) {install.packages("ggExtra")}
if (!require(zoo)) {install.packages("zoo")}

################################################################################################################################

################################################################################################################################

# Set working directory
setwd(WORKDIR)

# Load data

data_df <- read.table(Data, header = TRUE, sep = "\t")


# Split the dataframe in dataframe of 2 columns and store in a list
list <- lapply(seq(1, ncol(data_df), by=2), function(i)  
  data_df[i: pmin((i+1), ncol(data_df))])

# Name the element of the list by first column name of each dataframe
NAMES <- c()

for (i in 1:length(list)) {
  NAMES <- c(NAMES, colnames(list[[i]])[1])
}
names(list) <- NAMES

# Change columns name
for (i in 1:length(list)) {
  colnames(list[[i]]) <- c("WT", "tor-15")
}     

#Supprimer la première ligne

for (i in 1:length(list)) {
  list[[i]] <- list[[i]][-1, ]
} 

# Ajouter l'info du traitement dans une colonne

for (i in 1:length(list)) {
  list[[i]]$traitement <- rep(names(list[i]), ncol(list[[i]]))
} 

# join the element of the list
data <- do.call(rbind.data.frame, list)

# add missing data

# 5 mesures per root
data$day <- rep(c(0:(NB_DAY-1)), nrow(data)/NB_DAY)

# Reoder columns
data <- data[, c(3, 4, 1, 2)]

# Change class of the datas
data$WT <- as.numeric(data$WT)
data$traitement <- as.factor(data$traitement)
data$day <- as.factor(data$day)

# keep only value for the day of interest DAY
data <- data[data$day==DAY, ]

# make data tidy
data <- pivot_longer(data, where(is.numeric), names_to = "lines", values_to = ("FvFm"))
data <- data[order(data$lines),]


# The curve is perfporm à day 4
# Compute mean for DMSO at day 4 for both lines
summary_DMSO <- data[data$traitement == "DMSO", ] %>% dplyr::group_by(lines) %>%
  summarise(mean = mean(FvFm, na.rm = TRUE), sd= sd(FvFm, na.rm = TRUE))

# Create a new dataframe with the result of Az/DMSO*100

AZ_df <- data.frame(traitement = c(data$traitement[data$traitement != "DMSO"]), 
                    lines = c(data$lines[data$traitement != "DMSO"]))
# Compute % AZ/DMSO for each mesure
AZ_df$DMSO_percent <- (data$FvFm[data$traitement != "DMSO"]/
                         summary_DMSO$mean[match(AZ_df$lines, summary_DMSO$lines)])*100

# Compute mean and SD for AZ_df$DMSO_percent 

summary_AZ <- AZ_df %>% dplyr::group_by(lines, traitement) %>%
  summarise(mean = mean(DMSO_percent, na.rm = TRUE), sd= sd(DMSO_percent, na.rm = TRUE))

# Add standard error
# Compute the real number of roots that have been mesured
REAL_ROOT_NUMBER <- ((aggregate(FvFm ~ lines + traitement, data=data[data$traitement != "DMSO", ], function(x) {sum(!is.na(x))}, na.action = NULL))) 
colnames(REAL_ROOT_NUMBER) <- c("lines", "traitement", "number")

REAL_ROOT_NUMBER <- REAL_ROOT_NUMBER[order(REAL_ROOT_NUMBER$lines),]

summary_AZ$SE <- summary_AZ$sd/(sqrt(REAL_ROOT_NUMBER$number))


#___________________________________________________________________________________________________________________________________
# Statistical analysis
#___________________________________________________________________________________________________________________________________

library(ggbeeswarm) # pour la fonction geom_quasirandom
library(RColorBrewer) # pour la définition des couleurs
library(rstatix)  # pour les tests statistiques
library(rcompanion) #pour le calcul de l'intervalle de confiance pour les données non paramétriques
library(forcats)
library(scales)

# s'assurer que les package rmisc et plyr ne sont pas chargés`
#detach(package:Rmisc)
#detach(package:plyr)


# Define color

my_colours <- c("#4d4c4a", "#dbdad7")

# Créer un objet target_order pour forcer l'ordre dans les graph 
# Inhib$lines <- as.factor(Inhib$lines)
# L <- levels(Inhib$lines)
# L <- L[!(L %in% RefLine)]
# target_order <- c(RefLine, L)

# create a data frame with no NA
AZ_df_No_NA <- subset(AZ_df, subset =! is.na(AZ_df$DMSO_percent))

# Determining data normality status
shapiro_df <- AZ_df_No_NA %>%
  dplyr::group_by(lines, traitement) %>%
  summarise(statistic = shapiro.test(DMSO_percent)$statistic, 
            p.value = shapiro.test(DMSO_percent)$p.value)

flag_normal <- check_normality(shapiro_df)

# Data treatement according to normality status
if(flag_normal == TRUE) {
  print("Les données suivent une loi normale")
  
  # Summary
  if (!require(Rmisc)) {install.packages("Rmisc")}
  library(plyr) # dépendence de rmisc
  library(Rmisc) # pour la commande summarySE
  my_summary <- summarySE(AZ_df_No_NA, measurevar="DMSO_percent", 
                          groupvars=c("lines", "traitement"), na.rm = TRUE)
  
  detach(package:Rmisc)
  detach(package:plyr)
  
  # Plot
  AZ_df_No_NA$traitement <- factor(AZ_df_No_NA$traitement, levels= az_order)
  plot_normal(AZ_df_No_NA, my_colours, my_summary)
  
  # Stats
  anova_results <- AZ_df_No_NA %>% subset(., subset =! traitement== "AZD.0.003.µM") %>% 
    subset(., subset =! traitement== "AZD.0.01.µM") %>% 
    group_by(traitement) %>% anova_test(DMSO_percent ~ lines)
  flag_anova <- check_anova(anova_results)
  if (flag_anova == TRUE) {
    tukey_results <- as.data.frame(AZ_df_No_NA %>% subset(., subset =! traitement== "AZD.0.003.µM") %>% 
                                     subset(., subset =! traitement== "AZD.0.01.µM") %>% 
                                     mutate(lines = fct_relevel(lines, target_order)) %>%
                                     group_by(traitement) %>% tukey_hsd(DMSO_percent ~lines))
  }
  
  # Sauver les fichiers
  write.table(my_summary, file = "Summary.txt", 
              quote = FALSE, row.names = FALSE, sep = '\t')
  write.table(anova_results, file = "Anova.txt", 
              quote = FALSE, row.names = FALSE, sep = '\t')
  if (flag_anova == TRUE) {
    write.table(tukey_results[, c(1,3,4, 9, 10)], file = "Tukey.txt", 
                quote = FALSE, row.names = FALSE, sep = '\t')
  } 
  #ggsave("Plot.svg", width=4, height=5)
  
  
} else {
  print("Les données ne suivent pas une loi normale")
  
  # Summary
  conf_int <- groupwiseMedian(data = AZ_df_No_NA,
                              var = "DMSO_percent",
                              group = c("lines", "traitement"),
                              conf       = 0.95,
                              R          = 5000,
                              percentile = TRUE,
                              bca        = FALSE,
                              digits     = 3)
  
  
  # Plot
  AZ_df_No_NA$traitement <- factor(AZ_df_No_NA$traitement, levels= az_order)
  plot_not_normal(AZ_df_No_NA, my_colours, conf_int)
  
  # Stats
  kruskal_pval <- (AZ_df_No_NA %>% subset(., subset =! traitement== "AZD.0.003.µM") %>% 
                     subset(., subset =! traitement== "AZD.0.01.µM") %>% 
                     group_by(traitement) %>% kruskal_test(DMSO_percent ~ lines)) %>% 
    select(traitement, p)
  
  flag_kruskal <- check_kruskal(kruskal_pval)
  if (flag_kruskal == TRUE) { pval_dunn <- test_dunn(AZ_df) }
  
  
  # Sauver les fichiers
  write.table(conf_int, file = "Summary.txt", 
              quote = FALSE, row.names = FALSE, sep = '\t')
  write.table(kruskal_pval, file = "Kruskal.txt", 
              quote = FALSE, row.names = FALSE, sep = '\t')
  
  if (flag_kruskal == TRUE) {
    write.table(pval_dunn[, c(1, 3, 4, 9, 10)], file = "Dunn.txt", 
                quote = FALSE, row.names = FALSE, sep = '\t')
  } 
  #ggsave("Plot.svg", width=4, height=5)
}

# Curve
# order the dataframe numericaly
summary_AZ <- summary_AZ[order(factor(summary_AZ$traitement, levels=unique(az_order))),]
summary_AZ <- summary_AZ[order(summary_AZ$lines),]

# Replace traitement factors by numeric value
concentration <- as.numeric(c(0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10))
summary_AZ$traitement <- rep(concentration, 2)

x <- data.frame("tor-15", 0.001, 100, NA, NA)
y <- data.frame("WT", 0.001, 100, NA, NA)

summary_AZ <- rbind(setNames(x, names(summary_AZ)), summary_AZ)  
summary_AZ <- rbind(setNames(y, names(summary_AZ)), summary_AZ) 
summary_AZ <- summary_AZ[order(summary_AZ$lines),]

summary_AZ$mean[is.nan(summary_AZ$mean)]<-NA

library(zoo)
summary_AZ$mean <- na.approx(summary_AZ$mean)

library(ggExtra)
# Plot
p2 <- summary_AZ %>%
  mutate(lines = fct_relevel(lines, 
                             target_order)) %>%
  ggplot( aes(x = traitement, y = mean,  group = lines)) + 
  geom_line(aes(group = lines, color = lines, size = lines), show_guide=FALSE) +
  geom_point(aes(shape = lines, color = lines), size = 3) + 
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE, color = lines), width=.2,
                position=position_dodge(0.01)) +
  scale_x_continuous(trans= 'log10', breaks=c(0.001, 0.010, 0.100, 1.000, 10.000 ),labels=c(0, 0.01, 0.1, 1, 10)) +
  theme_light() +
  scale_colour_manual(values=my_colours, labels = c(expression("WT", italic("tor-15"))))+
  scale_shape_manual(values=c(17,16), labels = c(expression("WT", italic("tor-15"))))+
  scale_size_manual(values=c(1.5, 1.5))+
  ylim(0,105)+
  removeGridX() +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14),
        legend.key.width = unit(2, 'cm')) +
  labs(title=TITLE, x="", y = Y_AXIS)

print(p2)
p2 + annotation_logticks(sides = "b")


# Violin plot
library(Hmisc)

AZ_df$traitement <- factor(AZ_df$traitement, levels= az_order)

p4 <- AZ_df %>% mutate(lines = fct_relevel(lines, target_order)) %>% 
  ggplot(aes(x=lines, y=DMSO_percent, fill = lines, color = lines)) +
  scale_fill_manual(values=my_colours)+
  scale_color_manual(values = c(rep("black", 2)))+
  geom_violin(trim=FALSE)+
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red")+
  geom_jitter(shape=21, position=position_jitter(0.2))+
  theme_light()+
  removeGridX() +
  theme(legend.position="none",
        axis.text.x = element_text(size=14, angle=45,  hjust = 1),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(color="black"))+
  facet_wrap(~ traitement, strip.position = "bottom")+
  strip_pos +
  scale_x_discrete(labels= c("WT" = "WT", "tor-15" = expression(italic("tor-15"))))

p4
