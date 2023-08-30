#########################################################
#Script for Fv/Fm data analysis with statistical analysis
#########################################################

#### Part to modify #######################
#===============================================================================
# Variables 
##Reference line : the line againt which the statistical tests will be performed
RefLine <- "WT"  #a string between ""

# Order in wich the lines will appear on the plot
target_order <- c("WT", "tor-15", "TORp::TORG2268E")

# Working directory
WORKDIR <- "~/OneDrive-Aix-MarseilleUniversité/bio-informatique/figures_papier_romain"

# File to analyse
DATA <- "AZD1_FvFm.txt"
DMSO <- "DMSO_FvFm.txt"

# Nombre de Jours
NB_DAY <- 4

# Variables de personnalisation du graphique pour la courbe de croissance

TITLE <- ""
Y_AXIS <- ""

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
#==============================================================================
#Installation des packages
if (!require(tidyr)) { install.packages("tidyr") }
if (!require(dplyr)) { install.packages("dplyr") }
if (!require(tidyverse)) { install.packages("tidyverse") }
if (!require(devtools)) { install.packages("devtools") }
if (!require(rstatix)) { install.packages("rstatix", repos = "https://cloud.r-project.org") }
if (!require(ggbeeswarm)) { install.packages("ggbeeswarm") }
if (!require(RColorBrewer)) { install.packages("RColorBrewer") }
if (!require(rcompanion)) { install.packages("rcompanion") }
if (!require(ggExtra)) { install.packages("ggExtra") }
if (!require(svglite)) { install.packages("svglite") }

################################################################################################################################

# Fonction utilisée dans le script

#=======================================================================================================================================
# Check data normality. Get out of th loop if at least one of the group of data does not follow a normal law.
# return TRUE if data follow a normal law
#=======================================================================================================================================
check_normality <- function(shapiro_df) {
  # The normality is assumed to be true
  flag_normal <- TRUE
  
  for (i in 1 : nrow(shapiro_df)) {
    if(shapiro_df[i, 4] > 0.05) {
      # print(paste0("les données ",shapiro_df$grouping_factor[i],"-", 
      # shapiro_df$plant_line[i], " suivent une loi normale"), quote = FALSE)
      
    } else {
      # print(paste0("les données ",shapiro_df$grouping_factor[i],"-", 
      #          shapiro_df$plant_line[i], " ne suivent pas une loi normale"), quote = FALSE)
      
      # If data don't follow a normal law, don't go further
      flag_normal <- FALSE
      break
    }
  }
  return(flag_normal)
}

#=======================================================================================================================================
# Take the results of Anova test
# Return TRUE if at least one mean is significantly different 
#=======================================================================================================================================
check_anova <- function(anova_results) {
  flag_anova <- FALSE
  
  for (i in 1 : nrow(anova_results)) {
    if (anova_results$p[i] < 0.05) {
      print (paste0("Anova test compares the mean, pvalue for the group ", anova_results$Day[i] , " is < 0.05 this means that at least one of the medians is significantly different from the others, a post hoc Tukey test is performed"))
      flag_anova <- TRUE
    } else {
      print (paste0("Anova test compares the mean, pvalue for the group", anova_results$Day[i] , " is > 0.05 this means that there is no significant difference between means"))
    }
  }    
  
  return(flag_anova)
}

#=======================================================================================================================================
# Perform Kruskal-Wallis test
# Return TRUE if one of the median is significantly different
#=======================================================================================================================================
check_kruskal <- function(kruskal_pval) {
  flag_kruskal <- FALSE
  
  for (i in 1 : nrow(kruskal_pval)) {
    if (kruskal_pval$p[i] < 0.05) {
      print (paste0("Kruskall Wallis test compares medians, pvalue for the group ", kruskal_pval$Day[i] , " is < 0.05 this means that at least one of the medians is significantly different from the others, a post hoc Dunn test is performed"))
      flag_kruskal <- TRUE
    } else {
      print (paste0("Kruskall Wallis test compares medians, pvalue for the group ", kruskal_pval$Day[i] , " is > 0.05 this means that there is no significant difference between medians"))
    }
  }    
  
  return(flag_kruskal)
}

#=======================================================================================================================================
# Perform Dunn test
# Return the pvalue in a dataframe
#=======================================================================================================================================
test_dunn <- function() {
  pval <- as.data.frame(df_data_long %>%
                          mutate(lines = fct_relevel(lines, target_order)) %>% 
                          group_by(Day) %>% dunn_test(FvFm ~ lines, p.adjust.method = "BH"))
  #print(df_data_long %>% group_by(Day) %>% dunn_test(value ~ line, p.adjust.method = "BH"))
  return(pval)
}

#==============================================================================
#Define Working directory
setwd(WORKDIR)

#Load the data
df_data <- read.table(DATA, header = TRUE, sep = "\t")
colnames(df_data) <-c("WT", "tor-15", "TORp::TORG2268E")

#Convert to long format
library(tidyr)
library(tidyverse)
#https://tidyr.tidyverse.org/reference/pivot_longer.html
df_data_long <- df_data %>% pivot_longer(cols = everything(), names_to = "lines", values_to = "FvFm")
df_data_long <- df_data_long[order(df_data_long$lines),]

# Séparer les disques (4 mesure par racine)
df_data_long$Day <- rep(c(0:NB_DAY), nrow(df_data_long)/(NB_DAY+1))

# reoder column
df_data_long <- df_data_long[, c(3, 1, 2)]





#___________________________________________________________________________________________________________________________________
# Statistical analysis
#___________________________________________________________________________________________________________________________________

library(ggbeeswarm) # for funcnction geom_quasirandom
library(RColorBrewer) # to define colors
library(rstatix)  # for statistical test
library(rcompanion) # to compute confidence interval for non parametric data
library(forcats)


# Convert line's name into factor to force the order in the plot 
df_data_long$lines <- as.factor(df_data_long$lines)


# Determining data normality status
shapiro_df <- df_data_long %>%
  dplyr::group_by(Day, lines) %>%
  summarise(statistic = shapiro.test(FvFm)$statistic, 
            p.value = shapiro.test(FvFm)$p.value)

flag_normal <- check_normality(shapiro_df)

# Data treatement according to normality status
if(flag_normal == TRUE) {
  print("Datas follow a normal law")
  # Stats
  anova_results <- df_data_long %>% group_by(Day) %>% anova_test(FvFm ~ lines)
  flag_anova <- check_anova(anova_results)
  if (flag_anova == TRUE) {
    tukey_results <- as.data.frame(df_data_long %>% mutate(lines = fct_relevel(lines, target_order)) %>%
                                     group_by(Day) %>%  tukey_hsd(FvFm ~ lines))
  }
    
  # Save the files
  write.table(anova_results, file = "Anova.txt", 
              quote = FALSE, row.names = FALSE, sep = '\t')
  if (flag_anova == TRUE) {
    write.table(tukey_results[, c(1,3,4, 9, 10)], file = "Tukey.txt", 
                quote = FALSE, row.names = FALSE, sep = '\t')
  } else {
    print("Datas are not significantly different, The Tukey test was not performed.")
  } 
} else {
  print("Datas don't follow a normal law")
  
  # Stats
  kruskal_pval <- (df_data_long %>% group_by(Day)%>%kruskal_test(FvFm ~ lines)) %>% select(Day, p)
  
  flag_kruskal <- check_kruskal(kruskal_pval)
  if (flag_kruskal == TRUE) { pval_dunn <- test_dunn() }
  
  # Save the files
  write.table(kruskal_pval, file = "Kruskal.txt", 
              quote = FALSE, row.names = FALSE, sep = '\t')
  
  if (flag_kruskal == TRUE) {
    write.table(pval_dunn[, c(1, 3, 4, 8, 9, 10)], file = "Dunn.txt", 
                quote = FALSE, row.names = FALSE, sep = '\t')
  }
} 
  
#___________________________________________________________________________________________________________________________________
# Courbe de suivi dans le temps
#___________________________________________________________________________________________________________________________________
  
#calcul de la moyenne et Standard déviation
  
summary <- df_data_long %>% dplyr::group_by(Day, lines) %>%
  summarise(mean = mean(FvFm, na.rm = TRUE), sd= sd(FvFm, na.rm = TRUE))
  
# ajouter l"erreur standard
# Compter le nombre réel de racine mesurée
REAL_disc_NUMBER <- ((aggregate(FvFm ~ lines, data=df_data_long, function(x) {sum(!is.na(x))}, na.action = NULL)) %>% 
                         pull(FvFm, lines))/(NB_DAY+1)
  
  
summary$SE <- summary$sd/sqrt(REAL_disc_NUMBER)
  
  
# Courbe de suivi dans le temps
  
library(ggExtra)
# Plot

# Define color
my_colours = c("#4d4c4a", "#dbdad7", "#b8b6b0")

p2 <- summary %>%
  mutate(lines = fct_relevel(lines, 
                               target_order)) %>%
  ggplot( aes(x = Day, y = mean, group = 1, color = lines)) + 
  geom_line(aes(group = lines, color = lines, size = lines), show_guide=FALSE) +
  geom_point(aes(shape = lines, color = lines), size = 3) + 
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width=.2,
                  position=position_dodge(0.05)) +
  theme_light() +
  scale_colour_manual(values=my_colours, 
                      labels = c(expression("WT", italic("tor-15"), italic("TORp::TOR"^G2268E))))+
  scale_shape_manual(values=c(17, 16, 18), 
                      labels = c(expression("WT", italic("tor-15"), italic("TORp::TOR"^G2268E))))+
  scale_size_manual(values=c(1.5, 1.5, 1.5))+
  removeGridX() +
  scale_y_continuous(limits = c(0, 0.9), breaks = seq(0, 0.9, 0.1)) +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14),
        legend.key.width = unit(2, 'cm')) +
  theme(legend.position="bottom") +     
  labs(title=TITLE, x="Day", y = Y_AXIS)
  
print(p2)

library(svglite)  
ggsave("fvfm_plot.svg", width=8, height=3)
    
  
#___________________________________________________________________________________________________________________________________
# Courbe de suivi dans le temps en pourcentage par rapport au DMSO
#___________________________________________________________________________________________________________________________________

# Préparation des données DMSO
df_dmso <- read.table(DMSO, header = TRUE, sep = "\t")
colnames(df_dmso) <-c("WT", "tor-15", "TORp::TORG2268E")

#Convert to long format

df_dmso_long <- df_dmso %>% pivot_longer(cols = everything(), names_to = "lines", values_to = "FvFm")
df_dmso_long <- df_dmso_long[order(df_dmso_long$lines),]

# Séparer les disques (4 mesure par racine)
df_dmso_long$Day <- rep(c(0:NB_DAY), nrow(df_dmso_long)/(NB_DAY+1))

# reoder column
df_dmso_long <- df_dmso_long[, c(3, 1, 2)]

# Compute mean and SD for DMSO on day DAY

summary_dmso <- df_dmso_long %>% dplyr::group_by(lines, Day) %>%
  summarise(mean = mean(FvFm, na.rm = TRUE), sd= sd(FvFm, na.rm = TRUE))

# Compute growth on inhibitor in % of growth on DMSO
# split Inih dataframe in a list of dataframe
L_dmso <- split(summary_dmso, seq(nrow(summary_dmso)))


L <- split(df_data_long, f = list(df_data_long$lines, df_data_long$Day))
L <-  L[order(names(L))]

# Divide all element of the list by the corresponding value in summary_DMSO
for (l in 1:length(L)) {
  L[[l]]$pourcent <- L[[l]]$FvFm/L_dmso[[l]][1, 3, drop = TRUE]*100
}

# Join all element of the list in one dataframe
Data <- bind_rows(L, .id = "column_label")
Data <- Data[,-1]

#___________________________________________________________________________________________________________________________________
# Statistical analysis
#___________________________________________________________________________________________________________________________________


# Convert line's name into factor to force the order in the plot 
Data$lines <- as.factor(Data$lines)


# Determining data normality status
shapiro_df_pourcent <- Data %>%
  dplyr::group_by(Day, lines) %>%
  summarise(statistic = shapiro.test(pourcent)$statistic, 
            p.value = shapiro.test(pourcent)$p.value)

flag_normal_pourcent <- check_normality(shapiro_df_pourcent)

# Data treatement according to normality status
if(flag_normal_pourcent == TRUE) {
  print("Datas follow a normal law")
  # Stats
  anova_results_pourcent <- Data %>% group_by(Day) %>% anova_test(pourcent ~ lines)
  flag_anova_pourcent <- check_anova(anova_results_pourcent)
  if (flag_anova_pourcent == TRUE) {
    tukey_results_pourcent <- as.data.frame(Data %>% mutate(lines = fct_relevel(lines, target_order)) %>%
                                     group_by(Day) %>%  tukey_hsd(pourcent ~ lines))
  }
  
  # Save the files
  write.table(anova_results_pourcent, file = "Anova_pourcent.txt", 
              quote = FALSE, row.names = FALSE, sep = '\t')
  if (flag_anova_pourcent == TRUE) {
    write.table(tukey_results_pourcent[, c(1,3,4, 9, 10)], file = "Tukey_pourcent.txt", 
                quote = FALSE, row.names = FALSE, sep = '\t')
  } else {
    print("Datas are not significantly different, The Tukey test was not performed.")
  } 
} else {
  print("Datas don't follow a normal law")
  
  # Stats
  kruskal_pval_pourcent <- (Data %>% group_by(Day)%>%kruskal_test(pourcent ~ lines)) %>% select(Day, p)
  
  flag_kruskal_pourcent <- check_kruskal(kruskal_pval_pourcent)
  if (flag_kruskal == TRUE) { pval_dunn_pourcent <- test_dunn() }
  
  # Save the files
  write.table(kruskal_pval_pourcent, file = "Kruskal_pourcent.txt", 
              quote = FALSE, row.names = FALSE, sep = '\t')
  
  if (flag_kruskal_pourcent == TRUE) {
    write.table(pval_dunn_pourcent[, c(1, 3, 4, 8, 9, 10)], file = "Dunn_pourcent.txt", 
                quote = FALSE, row.names = FALSE, sep = '\t')
  }
} 

#calcul de la moyenne et Standard déviation pour le pourcentage

summary_pourcent <- Data %>% dplyr::group_by(Day, lines) %>%
  summarise(mean = mean(pourcent, na.rm = TRUE), sd= sd(pourcent, na.rm = TRUE))

# ajouter l"erreur standard

summary_pourcent$SE <- summary_pourcent$sd/sqrt(REAL_disc_NUMBER)



# plot data
p3 <- summary_pourcent %>%
  mutate(lines = fct_relevel(lines, 
                             target_order)) %>%
  ggplot( aes(x = Day, y = mean, group = 1, color = lines)) + 
  geom_line(aes(group = lines, color = lines, size = lines), show_guide=FALSE) +
  geom_point(aes(shape = lines, color = lines), size = 3) + 
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width=.2,
                position=position_dodge(0.05)) +
  theme_light() +
  
  scale_colour_manual(values=my_colours, 
                      labels = c(expression("WT", italic("tor-15"), italic("TORp::TOR"^G2268E))))+
  scale_shape_manual(values=c(17, 16, 18), 
                     labels = c(expression("WT", italic("tor-15"), italic("TORp::TOR"^G2268E))))+
  scale_size_manual(values=c(1.5, 1.5, 1.5))+
  ylim(0,105) +
  removeGridX() +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14),
        legend.key.width = unit(2, 'cm')) +
  labs(title=TITLE, x="Day", y = Y_AXIS)

print(p3)

# Distribution violin plot

library(Hmisc)
p4 <- Data %>% mutate(lines = fct_relevel(lines, target_order)) %>% 
  ggplot(aes(x=lines, y=pourcent, fill = lines, color = lines)) +
  scale_fill_manual(values=my_colours)+
  scale_color_manual(values = c(rep("black", 3)))+
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
  facet_wrap(~ Day, strip.position = "bottom")+
  strip_pos +
  scale_x_discrete(labels= c("WT" = "WT", "tor-15" = expression(italic("tor-15")), 
                             "TORp::TORG2268E" = expression(italic("TORp::TOR"^G2268E))))

p4

