##%######################################################%##
#                                                          #
#### Principal component analysis for rust test 1 Using ####
####              LEA package Firza Riany               ####
####      (firzariany2@gmail.com) 1 September 2020      ####
#                                                          #
##%######################################################%##

# Libraries
library(LEA)  # For ancestry coefficients
library(RColorBrewer)  # To get plot colors
library(readxl)  # To load the population names
library(tidyr)  # To tidy the data set
library(ggplot2)  # To plot the admixture coefficients
library(forcats)
library(ggthemes)
library(patchwork)

# Setting working directory
setwd("C:/Users/firza/
      OneDrive - University of Eastern Finland/Internship (thesis)/
      Scripts and results/PCA-population structure")

# Loading and converting the SNP dataset ----
# Ped file should be converted because LEA can only read lfmm dataset
genot_data <- ped2lfmm(input.file = "SNP_rust_filter.ped", 
                       output.file = "SNP_rust_filter.lfmm", 
                       force = TRUE)

genot <- read.lfmm(genot_data)

# Determining the number of clusters ----
# We tried with K ranging from 1:15 to test the cross-entropy
# Our wild poplar is gathered from 12 populations, clearly a group, so we want
# To keep the original numbers of 12 populations within the range
#pop_structure <- snmf(input.file = genot_data, K = 1:15, 
                      #ploidy = 2, entropy = T, 
                      #alpha = 100, project = "new")
#summary(pop_structure)

pop_structure2 <- snmf(input.file = genot_data, K = 1:15,
                       ploidy = 2, entropy = T,
                       alpha = 100, repetitions = 5,
                       project = "new")
summary(pop_structure2)

plot(pop_structure2, col = "blue4", cex = 1.4, pch = 19)
# Knee point at 4, so 4 ancestral pops

# Creating a matrix with the values of admixture coefficients ----
qmatrix <- Q(object = pop_structure2, K = 6, run = 1)  
head(qmatrix)

# Assigning the names of ancestral populations K
colnames(qmatrix) <- c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6")

# Assigning individual names + population of origin 
df_qmatrix <- data.frame(qmatrix)

SNP <- read.table("SNP_rust_filter.ped")  # To get the individual names
ind_names <- as.vector(SNP$V1)
df_qmatrix$Ind <- ind_names

# To get the population names
pop <- read_xlsx("B4EST_Selection_tests_multisouches_Task1.3.xlsx",
                 sheet = 1)
pop_filter <- pop[which(pop$Nom_genotype %in% ind_names), ]
df_qmatrix$Pop <- pop_filter$Population
df_qmatrix <- df_qmatrix[order(df_qmatrix$Pop), ]

# Ordering the columns
df_qmatrix <- df_qmatrix[, c("Ind", "Pop", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6")]  # To be exported

# Tidying the dataset, from wide to long
df_qmatrix_long <- gather(df_qmatrix, key = "Coeff", 
                          value = "Q", -c("Ind", "Pop"))


# Plotting the admixture coefficients ----
Q6_plot <- 
  ggplot(data = df_qmatrix_long, aes(x = factor(Ind), y = Q, fill = factor(Coeff))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(Pop), switch = "x", scales = "free", space = "free") +
  theme_minimal() + 
  labs(x = "Individuals", title = "K = 6", y = "Admixture coefficients") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = expansion(add = 4)) +
  theme(
    panel.spacing.x = unit(x = 0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_gdocs(guide = FALSE)

# Saving the dataset for GWAS ----
write.table(x = df_qmatrix[, -2], 
            file = "covariate_popstruct.txt", 
            sep = "\t", 
            dec = ".",
            row.names = FALSE,
            quote = FALSE)