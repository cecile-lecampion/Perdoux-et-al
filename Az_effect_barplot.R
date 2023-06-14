###########################################################
# Create a barplot at the selectesd day for DMSO vs AZ data
###########################################################
#### Part to modify #######################
#===============================================================================
# Variables 
WORKDIR <- "~/partage/bio-informatique/figures_papier_romain"

Data_DMSO <- "data_barplot_DMSO.txt"
Data_Inhib <- "data_barplot_AZ.txt"

# Number of day
NB_DAY <- 4

# Number of roots
NB_ROOT <- 24

# Day of interest
DAY <- 4

# Reference line
RefLine <- "WT"

# lines real names
lines_names <- c("WT", "tor-15", "TORp::TORG2268E.1", "TORp::TORG2268E.2", "TORp::TORG2268E.3" )

# Target order for order on the graph
target_order <- c("WT", "tor-15", "TORp::TORG2268E.1", "TORp::TORG2268E.2", "TORp::TORG2268E.3")

# Variables to customize the plot

TITLE <- ""
Y_AXIS <- ""


# Variables to customize the confidence interval plot

# Color of points
COULEUR <- "#dbdad7"
  


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
if (!require(ggtext)) { install.packages("ggtext") }
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
    if(shapiro_df[i, 3] > 0.05) {
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
# Take Kruskal-Wallis results
# Return TRUE if at least one of the median is different
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
# Take Anova test results
# Return TRUE if at list one of the mean is different
#=======================================================================================================================================
check_anova <- function(anova_results) {
  flag_anova <- FALSE
  
  for (i in 1 : nrow(anova_results)) {
    if (anova_results$p[i] < 0.05) {
      print (paste0("Le test Anova compare les moyennes, la pvalue pour le jour ", anova_results$lines[i] , " est < 0.05 ce qui indique qu’au moins 1 des moyennes est différentes des autres, on réalise un test post hoc de Tukey"))
      flag_anova <- TRUE
    } else {
      print (paste0("Le test Anova compare les moyennes, la pvalue pour le jour ", anova_results$lines[i] , " est > 0.05 ce qui indique qu’il n'y a pas de différence entre les moyennes"))
    }
  }    
  
  return(flag_anova)
}

#=======================================================================================================================================
# Make Dunn test
# Return pvalues in a dataframe
#=======================================================================================================================================
test_dunn <- function(df_data) {
  pval <- as.data.frame(df_data%>%
                          mutate(lines = fct_relevel(lines, target_order)) %>% 
                          dunn_test(DMSO_percent ~ lines, p.adjust.method = "BH"))
  #print(df_data %>% group_by(lines) %>% dunn_test(DMSO_percent, p.adjust.method = "BH"))
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
    ggplot( aes(x=lines, y=DMSO_percent)) +
    geom_quasirandom(dodge.width=0.8,alpha = 0.6, colour=COULEUR) +
    geom_pointrange(data=my_summary%>% mutate(lines = fct_relevel(lines,  target_order)), 
                    aes(ymin=DMSO_percent - ci, ymax=DMSO_percent + ci),
                    position=position_dodge(width=0.8)) + 
    scale_colour_manual(values="black") +
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "DMSO_percent",
                       expand = expansion(mult = c(.1, .1))) +    # Extended scale 10% above the higher point and below the lower point 
    facet_wrap(~ lines, strip.position = "bottom", scales = "free") +
    theme_classic() + 
    strip_pos +
    orientation_xlabels +
    theme(legend.position = "none", axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  print(p)
}


#=======================================================================================================================================
# Make the plot when data does not follow a normal law 
# Confidence interval is build around the median
#=======================================================================================================================================
plot_not_normal <- function(df, my_colours, conf_int) {
  # Column "median" have to be named "mesured_value" for the ggplot
  names(conf_int)[3] <- "DMSO_percent"
  
  p <- df %>%
    mutate(lines = fct_relevel(lines, 
                              target_order)) %>%
    ggplot( aes(x=lines, y=DMSO_percent)) +
    geom_quasirandom(dodge.width=0.8, alpha = 0.6, colour=COULEUR) +
    geom_pointrange(data=conf_int%>% mutate(lines = fct_relevel(lines,  target_order)), 
                    aes(ymin=Percentile.lower, ymax=Percentile.upper), 
                    position=position_dodge(width=0.8)) +
    scale_colour_manual(values="black") +
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "",
                       expand = expansion(mult = c(.1, .1))) + # Extended scale 10% above the higher point and below the lower point 
    facet_wrap(~ lines, strip.position = "bottom", scales = "free") +
    theme_classic() +  
    strip_pos +
    orientation_xlabels +
    theme(legend.position = "none", axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  print(p)
}
################################################################################################################################

# Set working directory
setwd(WORKDIR)

# Load data

DMSO <- read.table(Data_DMSO, header = TRUE, sep = "\t")
colnames(DMSO) <- lines_names

Inhib <- read.table(Data_Inhib, header = TRUE, sep = "\t")
colnames(Inhib) <- lines_names

# pivot long

DMSO <- pivot_longer(DMSO, everything(), names_to = "lines", values_to = ("length"))
DMSO <- DMSO[order(DMSO$lines),]


Inhib <- pivot_longer(Inhib, everything(), names_to = "lines", values_to = ("length"))
Inhib <- Inhib[order(Inhib$lines),]


# Add the day of mesure for each row.
DMSO$day <- rep(c(1:NB_DAY), nrow(DMSO)/NB_DAY)
DMSO <- DMSO[, c(1,3,2)]

# keep only value for the day of interest DAY
DMSO <- DMSO[DMSO$day==DAY, ]

Inhib$day <- rep(c(1:NB_DAY), nrow(Inhib)/NB_DAY)
Inhib <- Inhib[, c(1,3,2)]

# keep only value for the day of interest DAY
Inhib <- Inhib[Inhib$day==DAY, ]

# Compute mean and SD for DMSO on day DAY

summary_DMSO <- DMSO %>% dplyr::group_by(lines) %>%
  summarise(mean = mean(length, na.rm = TRUE), sd= sd(length, na.rm = TRUE))

# Compute % AZ/DMSO for each mesure
Inhib$DMSO_percent <- (Inhib$length/summary_DMSO$mean[match(Inhib$lines, summary_DMSO$lines)])*100

# Compute mean and SD for Inhib$DMSO_percent on day DAY

summary_Inhib <- Inhib %>% dplyr::group_by(lines) %>%
  summarise(mean = mean(DMSO_percent, na.rm = TRUE), sd= sd(DMSO_percent, na.rm = TRUE))

# Add standard error
# Compute the real number of roots that have been mesured
REAL_ROOT_NUMBER <- ((aggregate(length ~ lines, data=Inhib, function(x) {sum(!is.na(x))}, na.action = NULL))) 
colnames(REAL_ROOT_NUMBER) <- c("lines", "number")

summary_Inhib$SE <- summary_Inhib$sd/sqrt(REAL_ROOT_NUMBER$number)

#___________________________________________________________________________________________________________________________________
# Statistical analysis
#___________________________________________________________________________________________________________________________________

library(ggbeeswarm) # for geom_quasirandom()
library(RColorBrewer) # for colors
library(rstatix)  # for statistical test
library(rcompanion) #to compute confidence interval for non parametric data
library(forcats)
library(ggtext)

# be sure that packages rmisc and plyr are not loaded
#detach(package:Rmisc)
#detach(package:plyr)



# Define color

my_colours <- c("#4d4c4a", "#dbdad7", rep("#b8b6b0", length(c(3 : nrow(summary_Inhib)))))


# Determining data normality status
shapiro_df <- Inhib %>%
  dplyr::group_by(lines) %>%
  summarise(statistic = shapiro.test(DMSO_percent)$statistic, 
            p.value = shapiro.test(DMSO_percent)$p.value)

flag_normal <- check_normality(shapiro_df)

# Data treatement according to normality status
if(flag_normal == TRUE) {
  print("Les données suivent une loi normale")
  
  # Summary
  if (!require(Rmisc)) {install.packages("Rmisc")}
  library(plyr) # dépendence de rmisc
  library(Rmisc) # For summarySE
  my_summary <- summarySE(Inhib, measurevar="DMSO_percent", groupvars=c("lines"), na.rm = TRUE)
  
  detach(package:Rmisc)
  detach(package:plyr)
  
  # Plot
  new_labels <- c(c("WT" = "WT", "tor-15" = expression(italic("tor-15")), 
                    "TORp::TORG2268E.1" = expression(italic("TORp::TORG2268E.1")), 
                    "TORp::TORG2268E.2" = expression(italic("TORp::TORG2268E.2")), 
                    "TORp::TORG2268E.3" = expression(italic("TORp::TORG2268E.3"))))
  plot_normal(Inhib, my_colours, my_summary)
  
  # Stats
  anova_results <- Inhib %>% anova_test(DMSO_percent ~ lines)
  flag_anova <- check_anova(anova_results)
  if (flag_anova == TRUE) {
    tukey_results <- as.data.frame(Inhib %>% mutate(lines = fct_relevel(lines, target_order)) %>%
                                     tukey_hsd(DMSO_percent ~lines))
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
  conf_int <- groupwiseMedian(data = Inhib,
                              var = "DMSO_percent",
                              group = c("lines"),
                              conf       = 0.95,
                              R          = 5000,
                              percentile = TRUE,
                              bca        = FALSE,
                              digits     = 3)
  
  
  # Plot
    
  plot_not_normal(Inhib, my_colours, conf_int)
  
  # Stats
  kruskal_pval <- (Inhib %>% kruskal_test(DMSO_percent ~ lines)) %>% select(.y., p)
  
  if (kruskal_pval$p < 0.05) {
    print ("Le test de Kruskall Wallis compare les médianes, la pvalue est < 0.05 ce qui indique qu’au moins 1 des médianes est différentes des autres, on réalise un test post hoc de Dunn")
    flag_kruskal <- TRUE
  } else {
    print ("Le test de Kruskall Wallis compare les médianes, la pvalue est > 0.05 ce qui indique qu’il n'y a pas de différence entre les médianes")
  }
  
  if (flag_kruskal == TRUE) { pval_dunn <- test_dunn(Inhib) }
  
  # Save files
  write.table(conf_int, file = "Summary.txt", 
              quote = FALSE, row.names = FALSE, sep = '\t')
  write.table(kruskal_pval, file = "Kruskal.txt", 
              quote = FALSE, row.names = FALSE, sep = '\t')
  
  if (flag_kruskal == TRUE) {
    write.table(pval_dunn[, c(2, 3, 8, 9)], file = "Dunn.txt", 
                quote = FALSE, row.names = FALSE, sep = '\t')
  } 
  #ggsave("Plot.svg", width=4, height=5)
}


# Barplot 

# Preparation to add stars to the barplot. 
# the pairwise comparison must be done with reference line :RefLine
if(flag_normal == TRUE) {
  if(flag_anova == TRUE) {
    df_star <- dplyr::filter(tukey_results, grepl(paste0("\\b", RefLine, "\\b"), group1))
  }else {
    print("les données ne sont pas significativement différentes, le test de Tukey n'a pas été réalisé.")
  }
  
}else {
  if(flag_kruskal == TRUE) {
    df_star <- dplyr::filter(pval_dunn, grepl(paste0("\\b", RefLine, "\\b"), group1))
  }else {
    print("les données ne sont pas significativement différentes, le test de Dunn n'a pas été réalisé.")
  }
}


x <- c("line", RefLine, RefLine, rep(x="", 6))
df_star <- rbind(x, df_star)

# Join the 2 data frame
summary_Inhib <- cbind(summary_Inhib[match(target_order, summary_Inhib$lines),]  , df_star$p.adj.signif)
colnames(summary_Inhib)[5] <- "Labels"

library(ggExtra)
   

p2 <- summary_Inhib %>% mutate(lines = fct_relevel(lines, target_order)) %>% 
ggplot(aes(x=lines, y=mean, fill = lines, color = lines)) +
  scale_fill_manual(values=my_colours)+
  scale_color_manual(values = c(rep("black", length(c(1 : nrow(summary_Inhib))))))+
  geom_bar(stat="identity", width = 0.5)+
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
                                "TORp::TORG2268E.1" = expression(italic("TORp::TOR"^G2268E.1)),  
                                "TORp::TORG2268E.2" = expression(italic("TORp::TOR"^G2268E.2)), 
                                "TORp::TORG2268E.3" = expression(italic("TORp::TOR"^G2268E.3))))
 
