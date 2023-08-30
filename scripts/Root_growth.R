##############################################
# Root growth analyses
##############################################
#### Part to modify #######################
#===============================================================================
# Variables 

# Number of days
NB_DAY <- 4

# Number of roots
NB_ROOT <- 24

# reference line
RefLine <- "nom de la lignée de référence" # exemple : "WT"

# Name of inhibitor
inhibiteur <- "nom de l'inhibiteur"  # exemple : "WYE 0.6µM"

# Target order for order on the graph
target_order <- c("WT", "tor-15", "TORp::TORG2268E")


# Variables to customize the growth curve

TITLE <- ""
Y_AXIS <- ""


# Variables to customize the confidence interval plot

# Color of points
COULEUR <-  "#dbdad7"
  

if (!require(ggplot2)) { install.packages("ggplot2") }
library(ggplot2)
# Orientation of labels for X axis.
# For horizontal labels : angle = 0, vjust = 0, hjust=0.
# For 90° labels : angle = 90, vjust = 0.5, hjust=1
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
#===============================================================================
# Set working directory
setwd("~/onedrive_amu/bio-informatique/figures_papier_romain")

# Load data

DMSO <- read.table("Data_curve_DMSO.txt", header = TRUE)
colnames(DMSO) <- c("WT", "tor-15", "TORp::TORG2268E")

Inhib <- read.table("Data_curve_inhibiteur.txt", header = TRUE, sep = '\t')
colnames(Inhib) <- c("WT", "tor-15", "TORp::TORG2268E")

################################################################################################################################

# Packages installation 

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

################################################################################################################################

# Functions

#=======================================================================================================================================
# Check normality of the data. Break the loop if at least one of the group of data does not follow 
# a normal law.
# Return TRUE if the data follow a normal law
#=======================================================================================================================================
check_normality <- function(shapiro_df) {
  # On suppose que les données sont normales
  flag_normal <- TRUE
  
  for (i in 1 : nrow(shapiro_df)) {
    if(shapiro_df[i, 4] > 0.05) {
      
    } else {
     
      # En fait les données ne sont pas normales, pas besoin d'aller plus loin
      flag_normal <- FALSE
      break
    }
  }
  return(flag_normal)
}

#=======================================================================================================================================
# Take Anova test results
# Return TRUE if at list one of the mean is different
#=======================================================================================================================================
check_anova <- function(anova_results) {
  flag_anova <- FALSE
  
  for (i in 1 : nrow(anova_results)) {
    if (anova_results$p[i] < 0.05) {
      print (paste0("Le test Anova compare les moyennes, la pvalue pour le jour ", anova_results$Day[i] , " est < 0.05 ce qui indique qu’au moins 1 des moyennes est différentes des autres, on réalise un test post hoc de Tukey"))
      flag_anova <- TRUE
    } else {
      print (paste0("Le test Anova compare les moyennes, la pvalue pour le jour ", anova_results$Day[i] , " est > 0.05 ce qui indique qu’il n'y a pas de différence entre les moyennes"))
    }
  }    
  
  return(flag_anova)
}

#=======================================================================================================================================
# Take Kruskal-Wallis results
# Return TRUE if at least one of the median is different
#=======================================================================================================================================
check_kruskal <- function(kruskal_pval) {
  flag_kruskal <- FALSE
  
  for (i in 1 : nrow(kruskal_pval)) {
    if (kruskal_pval$p[i] < 0.05) {
      print (paste0("Le test de Kruskall Wallis compare les médianes, la pvalue pour le groupe ", kruskal_pval$Day[i] , " est < 0.05 ce qui indique qu’au moins 1 des médianes est différentes des autres, on réalise un test post hoc de Dunn"))
      flag_kruskal <- TRUE
    } else {
      print (paste0("Le test de Kruskall Wallis compare les médianes, la pvalue pour le groupe ", kruskal_pval$Day[i] , " est > 0.05 ce qui indique qu’il n'y a pas de différence entre les médianes"))
    }
  }    
  
  return(flag_kruskal)
}

#=======================================================================================================================================
# Make Dunn test
# Return pvalues in a dataframe
#=======================================================================================================================================
test_dunn <- function(df_data) {
  pval <- as.data.frame(df_data%>%
                          mutate(lines = fct_relevel(lines, target_order)) %>% 
                          group_by(Day) %>% dunn_test(length ~ lines, p.adjust.method = "BH"))
  print(df_data %>% group_by(Day) %>% dunn_test(length ~ lines, p.adjust.method = "BH"))
  return(pval)
}

#=======================================================================================================================================
# Make the plot when data follow a normal law 
# Confidence interval is build around the mean
#=======================================================================================================================================
plot_normal <- function(df, my_colours, my_summary) {
  p <- df %>%
    mutate(lines = fct_relevel(lines, 
                              target_order)) %>%
    ggplot( aes(x=lines, y=length)) +
    geom_quasirandom(dodge.width=0.8,alpha = 0.6, colour=COULEUR) +
    geom_pointrange(data=my_summary%>% mutate(lines = fct_relevel(lines,  target_order)), 
                    aes(ymin=pourcent - ci, ymax=pourcent + ci),
                    position=position_dodge(width=0.8)) + 
    scale_colour_manual(values="black") +
    scale_x_discrete(labels = c(expression("WT", italic("tor-15"), italic("TORp::TOR"^G2268E)))) +
    scale_y_continuous(name = "length",
                       expand = expansion(mult = c(.1, .1))) +    # Extended scale 10% above the higher point and below the lower point 
    facet_wrap(~ Day, strip.position = "bottom", scales = "free") +
    theme_classic() + 
    strip_pos +
    orientation_xlabels +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))
  print(p)
}


#=======================================================================================================================================
# Make the plot when data does not follow a normal law 
# Confidence interval is build around the median
#=======================================================================================================================================
plot_not_normal <- function(df, my_colours, conf_int) {
  # Column "median" have to be named "mesured_value" for the ggplot
  names(conf_int)[4] <- "pourcent"
  
  p <- df %>%
    mutate(lines = fct_relevel(lines, 
                              target_order)) %>%
    ggplot( aes(x=lines, y=pourcent)) +
    geom_quasirandom(dodge.width=0.8, alpha = 0.6, colour=COULEUR) +
    geom_pointrange(data=conf_int%>% mutate(lines = fct_relevel(lines,  target_order)), 
                    aes(ymin=Percentile.lower, ymax=Percentile.upper), 
                    position=position_dodge(width=0.8)) +
    scale_colour_manual(values="black") +
    scale_x_discrete(labels = c(expression("WT", italic("tor-15"), italic("TORp::TOR"^G2268E)))) +
    scale_y_continuous(name = "length",
                       expand = expansion(mult = c(.1, .1))) + # Extended scale 10% above the higher point and below the lower point 
    facet_wrap(~ Day, strip.position = "bottom", scales = "free") +
    theme_classic() +  
    strip_pos +
    orientation_xlabels +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))
  print(p)
}
################################################################################################################################
#### Main #######################

# Make data tidy

DMSO <- pivot_longer(DMSO, everything(), names_to = "lines", values_to = ("length"))
DMSO <- DMSO[order(DMSO$lines),]

Inhib <- pivot_longer(Inhib, everything(), names_to = "lines", values_to = ("length"))
Inhib <- Inhib[order(Inhib$lines),]

library(dplyr) # for functions :  %>%, group_by, summarise, select

# Add missing data

# Make rrots group (4 mesure per root)
DMSO$Day <- rep(c(1:NB_DAY), nrow(DMSO)/NB_DAY)
Inhib$Day <- rep(c(1:NB_DAY), nrow(Inhib)/NB_DAY)

# Add treatment data
DMSO$traitement <- rep("DMSO", nrow(DMSO))
Inhib$traitement <- rep("Inhib", nrow(DMSO))

# reoder column
DMSO <- DMSO[, c(4, 3, 1, 2)]
Inhib <- Inhib[, c(4, 3, 1, 2)]

# Compute mean and SD for DMSO on day DAY

summary_DMSO <- DMSO %>% dplyr::group_by(lines, Day) %>%
  summarise(mean = mean(length, na.rm = TRUE), sd= sd(length, na.rm = TRUE))

# Compute growth on inhibitor in % of growth on DMSO
# split Inih dataframe in a list of dataframe
L_DMSO <- split(summary_DMSO, seq(nrow(summary_DMSO)))


L <- split(Inhib, f = list(Inhib$lines, Inhib$Day))
L <-  L[order(names(L))]

# Divide all element of the list by the corresponding value in summary_DMSO
for (l in 1:length(L)) {
  L[[l]]$pourcent <- L[[l]]$length/L_DMSO[[l]][1, 3, drop = TRUE]*100
}

# Join all element of the list in one dataframe
Data <- bind_rows(L, .id = "column_label")
Data <- Data[,-1]

#___________________________________________________________________________________________________________________________________
# Statistical analyse
#___________________________________________________________________________________________________________________________________

library(ggbeeswarm) # for geom_quasirandom()
library(RColorBrewer) # for colors
library(rstatix)  # for statistical test
library(rcompanion) #to compute confidence interval for non parametric data
library(forcats)

# be sure that packages rmisc and plyr are not loaded
#detach(package:Rmisc)
#detach(package:plyr)


# Define color
my_colours = c("#4d4c4a", "#dbdad7", "#b8b6b0")


# Determining data normality status
shapiro_df <- Data %>%
  dplyr::group_by(Day, lines) %>%
  summarise(statistic = shapiro.test(pourcent)$statistic, 
            p.value = shapiro.test(pourcent)$p.value)

flag_normal <- check_normality(shapiro_df)

# Data treatement according to normality status
if(flag_normal == TRUE) {
  print("Les données suivent une loi normale")
  
  # Summary
  if (!require(Rmisc)) {install.packages("Rmisc")}
  library(plyr) # dépendence de rmisc
  library(Rmisc) # pour la commande summarySE
  my_summary <- summarySE(Data, measurevar="pourcent", groupvars=c("Day", "lines"))
  
  detach(package:Rmisc)
  detach(package:plyr)
  
  # Plot
  
  plot_normal(Data, my_colours, my_summary)
  
  # Stats
  anova_results <- Data %>% group_by(Day) %>%  anova_test(pourcent ~ lines)
  flag_anova <- check_anova(anova_results)
  if (flag_anova == TRUE) {
    tukey_results <- as.data.frame(Data %>% mutate(lines = fct_relevel(lines, target_order)) %>%
                                     group_by(Day) %>%  tukey_hsd(pourcent ~ lines))
  }
  
  # Save files
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
  conf_int <- groupwiseMedian(data = Data,
                              var = "pourcent",
                              group = c("Day", "lines"),
                              conf       = 0.95,
                              R          = 5000,
                              percentile = TRUE,
                              bca        = FALSE,
                              digits     = 3)
  
  
  # Plot
  plot_not_normal(Data, my_colours, conf_int)
  
  # Stats
  kruskal_pval <- (Data %>% group_by(Day)%>%kruskal_test(pourcent ~ lines)) %>% select(Day, p)
  
  flag_kruskal <- check_kruskal(kruskal_pval)
  if (flag_kruskal == TRUE) { pval_dunn <- test_dunn(Data) }
  
  # Save files
  write.table(conf_int, file = "Summary.txt", 
              quote = FALSE, row.names = FALSE, sep = '\t')
  write.table(kruskal_pval, file = "Kruskal.txt", 
              quote = FALSE, row.names = FALSE, sep = '\t')
  
  if (flag_kruskal == TRUE) {
    write.table(pval_dunn[, c(1, 3, 4, 8, 9, 10)], file = "Dunn.txt", 
                quote = FALSE, row.names = FALSE, sep = '\t')
  } 
  #ggsave("Plot.svg", width=4, height=5)
}

#___________________________________________________________________________________________________________________________________
# Growth curve
#___________________________________________________________________________________________________________________________________

#Compute mean and Standard déviation

summary <- Data %>% dplyr::group_by(Day, lines) %>%
  summarise(mean = mean(pourcent, na.rm = TRUE), sd= sd(pourcent, na.rm = TRUE))

# Add Standard error
# Compute the real number of roots
REAL_ROOT_NUMBER <- ((aggregate(pourcent ~ lines, data=Data, function(x) {sum(!is.na(x))}, na.action = NULL)) %>% 
                       pull(pourcent, lines))/NB_DAY


summary$SE <- summary$sd/sqrt(REAL_ROOT_NUMBER)


data0 <- data.frame(c(0,0,0), c("WT", "tor-15", "TORp::TORG2268E"), c(100,100,100), c(0,0,0), c(0,0,0))
colnames(data0) <- colnames(summary)

summary <- rbind(data0, summary)

# Growth curve

library(ggExtra)
# Plot
p2 <- summary %>%
  mutate(lines = fct_relevel(lines, 
                            target_order)) %>%
  ggplot( aes(x = Day, y = mean, group = 1, color = lines)) + 
  geom_line(aes(group = lines, color = lines, size = lines), show_guide=FALSE) +
  geom_point(aes(shape = lines, color = lines), size = 3) + 
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width=.2,
                position=position_dodge(0.05)) +
  theme_light() +
  ylim(0,130) +
  scale_colour_manual(values=my_colours, 
                      labels = c(expression("WT", italic("tor-15"), italic("TORp::TOR"^G2268E))))+
  scale_shape_manual(values=c(17, 16, 18), 
                     labels = c(expression("WT", italic("tor-15"), italic("TORp::TOR"^G2268E))))+
  scale_size_manual(values=c(1.5, 1.5, 1.5))+
  removeGridX() +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14),
        legend.key.width = unit(2, 'cm')) +
  labs(title=TITLE, x="Day", y = Y_AXIS)

print(p2)

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

# Growth at Day 4
summary4<- summary[summary$Day== 4,]


# Preparation to add stars to the barplot. 
# the pairwise comparison must be done with reference line :RefLine
if(flag_normal == TRUE) {
  if(flag_anova == TRUE) {
    df_star <- dplyr::filter(tukey_results[tukey_results$Day== 4,], grepl(paste0("\\b", RefLine, "\\b"), group1))
  }else {
    print("les données ne sont pas significativement différentes, le test de Tukey n'a pas été réalisé.")
  }
  
}else {
  if(flag_kruskal == TRUE) {
    df_star <- dplyr::filter(pval_dunn[pval_dunn$Day== 4,], grepl(paste0("\\b", RefLine, "\\b"), group1))
  }else {
    print("les données ne sont pas significativement différentes, le test de Dunn n'a pas été réalisé.")
  }
}


x <- c("4", "length", RefLine, RefLine, rep(x="", 6))
df_star <- rbind(x, df_star)

# join 2 dataframes
summary4 <- cbind(summary4[match(target_order, summary4$lines),]  , df_star$p.adj.signif)
colnames(summary4)[6] <- "Labels"

p2 <- summary4 %>% mutate(lines = fct_relevel(lines, target_order)) %>% 
  ggplot(aes(x=lines, y=mean, fill = lines, color = lines)) +
  scale_fill_manual(values=my_colours)+
  scale_color_manual(values = c(rep("black", length(c(1 : nrow(summary4))))))+
  geom_bar(stat="identity", width = 0.3)+
  geom_errorbar(aes(x = as.factor(lines), ymin=mean-SE, ymax=mean+SE), width=.2)+
  geom_text(aes(x = as.factor(lines), y= mean+SE + 1 , 
                label = Labels),
            size = 5, inherit.aes = TRUE)+
  theme_light()+
  removeGridX() +
  theme(legend.position="none",
        axis.text.x = element_text(size=14, , angle=45,  hjust = 1),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

p2 + scale_x_discrete(labels= c("WT" = "WT", "tor-15" = expression(italic("tor-15")), 
                                "TORp::TORG2268E" = expression(italic("TORp::TOR"^G2268E))))  
                               
