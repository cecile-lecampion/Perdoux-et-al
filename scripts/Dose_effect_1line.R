#############################################################
# Create a growth curve for dose effect : 1 line, 1 inhibitor
#############################################################


# Prevent from scientific notation use
options(scipen=999)


#### Part to modify #######################
#===============================================================================
# Variables 
WORKDIR <- "~/partage/bio-informatique/figures_papier_romain"

Data <- "Data_curve_INK.txt"

# Target order for order on the graph

az_order <- c("INK.0.01.µM", "INK.0.1.µM", "INK.1.µM", "INK.10.µM")

# Variables to customize the plot

TITLE <- ""
Y_AXIS <- ""



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
    if(shapiro_df[i, ncol(shapiro_df)] > 0.05) {
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
# Take Anova test results
# Return TRUE if at list one of the mean is different
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
# Make Dunn test
# Return pvalues in a dataframe
#=======================================================================================================================================
test_dunn <- function(df_data) {
  pval <- as.data.frame(df_data%>%
                          dunn_test(DMSO_percent ~ Concentration, p.adjust.method = "BH"))
  #print(df_data %>% group_by(Day) %>% dunn_test(Length ~ line, p.adjust.method = "BH"))
  return(pval)
}


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
if (!require(zoo)) {install.packages("zoo")}

################################################################################################################################

################################################################################################################################

# Set working directory
setwd(WORKDIR)

# Load data

data_df <- read.table(Data, header = TRUE, sep = "\t")

# make data tidy
data <- pivot_longer(data_df, where(is.numeric), names_to = "Concentration", values_to = ("length"))
data <- data[order(data$Concentration),]

# make 2 data frame
DMSO_df <- data[data$Concentration=="DMSO",]
INK_df <- data[data$Concentration!="DMSO",]

# Compute mean for DMSO 
summary_DMSO <- DMSO_df %>%
  summarise(mean = mean(length, na.rm = TRUE), sd= sd(length, na.rm = TRUE))

# Create a new column in INK_df with the result of INK/DMSO*100
# Compute % INK/DMSO for each mesure
INK_df$DMSO_percent <- (INK_df$length/
                         summary_DMSO$mean*100)

# Compute mean and SD for INK_df$DMSO_percent 

summary_INK <- INK_df %>% dplyr::group_by(Concentration) %>%
  summarise(mean = mean(DMSO_percent, na.rm = TRUE), sd= sd(DMSO_percent, na.rm = TRUE))

# Add standard error
# Compute the real number of roots that have been mesured
REAL_ROOT_NUMBER <- ((aggregate(length ~ Concentration, data=INK_df, function(x) {sum(!is.na(x))}, na.action = NULL))) 
colnames(REAL_ROOT_NUMBER) <- c("Concentration", "number")

summary_INK$SE <- summary_INK$sd/(sqrt(REAL_ROOT_NUMBER$number))


#___________________________________________________________________________________________________________________________________
# Statistical analysis
#___________________________________________________________________________________________________________________________________

library(ggbeeswarm) # for geom_quasirandom()
library(RColorBrewer) # for colors
library(rstatix)  # for statistical test
library(rcompanion) #to compute confidence interval for non parametric data
library(forcats)
library(scales)

# be sure that packages rmisc and plyr are not loaded
#detach(package:Rmisc)
#detach(package:plyr)


# Define color

my_colours <- c("#4d4c4a", "#dbdad7")

# create a data frame with no NA
INK_df_No_NA <- subset(INK_df, subset =! is.na(INK_df$DMSO_percent))

# Determining data normality status
shapiro_df <- INK_df_No_NA %>%
  dplyr::group_by(Concentration) %>%
  summarise(statistic = shapiro.test(DMSO_percent)$statistic, 
            p.value = shapiro.test(DMSO_percent)$p.value)

flag_normal <- check_normality(shapiro_df)

# Data treatement according to normality status
if(flag_normal == TRUE) {
  print("Les données suivent une loi normale")
  
  # Stats
  anova_results <- INK_df_No_NA %>% anova_test(DMSO_percent ~ Concentration)
  flag_anova <- check_anova(anova_results)
  if (flag_anova == TRUE) {
    tukey_results <- as.data.frame(INK_df_No_NA %>% tukey_hsd(DMSO_percent ~Concentration))
  }
  
  # Save the files
  write.table(anova_results, file = "Anova.txt", 
              quote = FALSE, row.names = FALSE, sep = '\t')
  if (flag_anova == TRUE) {
    write.table(tukey_results[, c(2,3,8,9)], file = "Tukey.txt", 
                quote = FALSE, row.names = FALSE, sep = '\t')
  } 
  
  
} else {
  print("Les données ne suivent pas une loi normale")
  
  # Stats
  kruskal_pval <- (INK_df_No_NA %>% kruskal_test(DMSO_percent ~ Concentration)) %>% 
    select(p)
  
  flag_kruskal <- check_kruskal(kruskal_pval)
  if (flag_kruskal == TRUE) { pval_dunn <- test_dunn(INK_df) }
  
  # Save the files
  write.table(kruskal_pval, file = "Kruskal.txt", 
              quote = FALSE, row.names = FALSE, sep = '\t')
  
  if (flag_kruskal == TRUE) {
    write.table(pval_dunn[, c(1, 3, 4, 9, 10)], file = "Dunn.txt", 
                quote = FALSE, row.names = FALSE, sep = '\t')
  } 
 
}


# Curve

# Replace Concentration factors by numeric value
concentration <- as.numeric(c(0.01, 0.1, 1, 10))
summary_INK$Concentration <- concentration

# add a line to be able to have 0 on the graph
x <- data.frame(0.001, 100, NA, NA)
summary_INK <- rbind(setNames(x, names(summary_INK)), summary_INK)

library(zoo)
summary_AZ$mean <- na.approx(summary_AZ$mean)

library(ggExtra)
# Plot
p2 <- summary_INK %>%
  ggplot( aes(x = Concentration, y = mean)) + 
  geom_line(size = 2, show_guide=FALSE) +
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width=.2,
                position=position_dodge(0.01)) +
  scale_x_continuous(trans= 'log10', breaks=c(0.0010, 0.010, 0.100, 1.000, 10.000 ),labels=c(0, 0.01, 0.1, 1, 10)) +
  theme_light() +
  scale_size_manual(values=c(1.5, 1.5))+
  #ylim(0,105)+
  removeGridX() +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14),
        legend.key.width = unit(2, 'cm')) +
  labs(title=TITLE, x="", y = Y_AXIS)

print(p2)
p2 + annotation_logticks(sides = "b")

