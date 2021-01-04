##%######################################################%##
#                                                          #
####  Preliminary analysis for Rust test 1 (conducted   ####
####             in 15/07/2020) Firza Riany             ####
####       (firzariany2@gmail.com) 12 August 2020       ####
#                                                          #
##%######################################################%##

# Objectives ----
# 1. Verify the pathotype of the strains (for each strain)
# 2. Understand the global distribution of each phenotype
# 3. Understand the variability in distribution of the phenotypes for each strain
# 4. Observe the block effects
# 5. Provide the visualization of the data set

# Libraries ----
library(plyr)
library(naniar)
library(ggplot2)
library(breedR)  # to adjust trait variation based on block and residual effects

# Setting working directory
setwd("~/Doc4Firza/Scripts/Rust_test1 (for thesis)")

# Loading the raw data set ----
discrim <- read.table("Tab_discrim_vir.txt",
                      header = TRUE,
                      row.names = NULL)

rust_test1 <- read.delim2("B4est_test1_LectureJ13_J14_commented.txt", na = "NA")
str(rust_test1)

# Selection and coercion of column ----
test1_selec <- rust_test1[ ,-c(12, 13, 15)]

# Subsetting the raw dataset by strain
test1_09AX27 <- test1_selec[which(test1_selec$Strain == "09AX27"),]
test1_93JE3 <- test1_selec[which(test1_selec$Strain == "93JE3"),]
test1_P72 <- test1_selec[which(test1_selec$Strain == "P72"),]

# Verification of Mlp virulences and pathotype ----
# For strain 09AX27, virulence = 4
check_vir_09AX27 <- test1_09AX27[which(test1_09AX27$Genotype %in% 
                                         discrim$Discrim), ]

# For strain 93JE3, virulence = 2, 4
check_vir_93JE3 <- test1_93JE3[which(test1_93JE3$Genotype %in%
                                       discrim$Discrim), ]

# For strain P72, virulence = 3
check_vir_P72 <- test1_P72[which(test1_P72$Genotype %in%
                                   discrim$Discrim), ]

# Compatiblity and incompatibility of the strains ----
# Create a new column to differentiate normal genotypes (subject for inspection)
# From control (Genotype == Robusta) and discriminants (Genotype = Discrim)

## For strain 09AX27
test1_09AX27[ ,13] <- "NA"
for (i in c(1:nrow(test1_09AX27))){
  if(test1_09AX27$Genotype[i] %in% discrim$Discrim)
    test1_09AX27[i,13] <- "discrim"
  if(test1_09AX27$Genotype[i] == "Robusta")
    test1_09AX27[i,13] <- "control"
}

## Graphical inspection
test1_09AX27_nodisc <- test1_09AX27[-which(test1_09AX27[ ,13] %in%
                                             c("control", "discrim")), ]
par(mfrow = c(2,2))
a <- barplot(c(length(which(test1_09AX27_nodisc$Laten_24.02.15_gp == -5)),
               length(which(test1_09AX27_nodisc$Laten_24.02.15_gp != -5))),
             names.arg = c("incompatible", "compatible"),
             main = "Compatibility of strain 09AX27")
text(a, labels = c(length(which(test1_09AX27_nodisc$Laten_24.02.15_gp == -5)),
                   length(which(test1_09AX27_nodisc$Laten_24.02.15_gp != -5))),
     adj = 0,
     pos = 3,
     offset = 3,
     cex = 1.5)

# For strain 93JE3
test1_93JE3[ ,13] <- "NA"
for (i in c(1:nrow(test1_93JE3))){
  if(test1_93JE3$Genotype[i] %in% discrim$Discrim)
    test1_93JE3[i,13] <- "discrim"
  if(test1_93JE3$Genotype[i] == "Robusta")
    test1_93JE3[i,13] <- "control"
}

# Graphical inspection
test1_93JE3_nodisc <- test1_93JE3[-which(test1_93JE3[ ,13] %in%
                                             c("control", "discrim")), ]

b <- barplot(c(length(which(test1_93JE3_nodisc$Laten_24.02.15_gp == -5)),
               length(which(test1_93JE3_nodisc$Laten_24.02.15_gp != -5))),
             names.arg = c("incompatible", "compatible"),
             main = "Compatibility of strain 93JE3")
text(b, labels = c(length(which(test1_93JE3_nodisc$Laten_24.02.15_gp == -5)),
                   length(which(test1_93JE3_nodisc$Laten_24.02.15_gp != -5))),
     adj = 0,
     pos = 3,
     offset = 3,
     cex = 1.5)

# For strain P72
test1_P72[ ,13] <- "NA"
for (i in c(1:nrow(test1_P72))){
  if(test1_P72$Genotype[i] %in% discrim$Discrim)
    test1_P72[i,13] <- "discrim"
  if(test1_P72$Genotype[i] == "Robusta")
    test1_P72[i,13] <- "control"
}

# Graphical inspection
test1_P72_nodisc <- test1_P72[-which(test1_P72[ ,13] %in%
                                       c("control", "discrim")), ]
c <- barplot(c(length(which(test1_P72_nodisc$Laten_24.02.15_gp == -5)),
               length(which(test1_P72_nodisc$Laten_24.02.15_gp != -5))),
             names.arg = c("incompatible", "compatible"),
             main = "Compatibility of strain P72")
text(c, labels = c(length(which(test1_P72_nodisc$Laten_24.02.15_gp == -5)),
                   length(which(test1_P72_nodisc$Laten_24.02.15_gp != -5))),
     adj = 0,
     pos = 3,
     offset = 3,
     cex = 1.5)

# Counting the genotypes with LP = -5 (incompatible genotypes) ----
incomp_09AX27 <- plyr::count(test1_09AX27_nodisc$Genotype[
  which(test1_09AX27_nodisc$Laten_24.02.15_gp == -5)])

incomp_93JE3 <- plyr::count(test1_93JE3_nodisc$Genotype[
  which(test1_93JE3_nodisc$Laten_24.02.15_gp == -5)])

incomp_P72 <- plyr::count(test1_P72_nodisc$Genotype[
  which(test1_P72_nodisc$Laten_24.02.15_gp == -5)])

# Global distribution of the raw genotypes ----
# Data preparation:
## Removing -5 in all variables
test1_selec_NA <- test1_selec
test1_selec_NA <-  test1_selec_NA %>% 
  replace_with_na(replace = list(Laten_24.02.15_gp = -5,
                                 NbSS_24.02.15_gp = -5,
                                 NbSO_24.02.15_gp = -5,
                                 NbPtN_24.02.15_gp = -5,
                                 Taill_05.08.20_gp = -5))

## Total uredinia number = adding NbSS and NbSoS
test1_selec_NA$total_NbSS <- test1_selec_NA$NbSS_24.02.15_gp + 
  test1_selec_NA$NbSO_24.02.15_gp

## Total black points = substracting NbPtN by point noires avant le test
test1_selec_NA$PtN_real <- test1_selec_NA$NbPtN_24.02.15_gp - 
  test1_selec_NA$points.noirs.avant.le.test

test1_selec_NA <- test1_selec_NA %>%
  naniar::replace_with_na(replace = list(PtN_real = c(-16:-1)))

## Removing the discriminants and the control
test1_selec_NA[ ,15] <- "NA"
for (i in c(1:nrow(test1_selec_NA))){
  if(test1_selec_NA$Genotype[i] %in% discrim$Discrim)
    test1_selec_NA[i,15] <- "discrim"
  if(test1_selec_NA$Genotype[i] == "Robusta")
    test1_selec_NA[i,15] <- "control"
}

test1_NA_nodisc <- test1_selec_NA[-which(test1_selec_NA[ ,15] %in%
                                           c("control", "discrim")), ]

## Histogram of each trait
par(mfrow = c(2,2))

test1_loop <- data.frame("LP" = test1_NA_nodisc$Laten_24.02.15_gp,
                         "UN" = test1_NA_nodisc$total_NbSS,
                         "US" = test1_NA_nodisc$Taill_05.08.20_gp,
                         "Black_Points" = test1_NA_nodisc$PtN_real)
loop.vector <- 1:4

for (i in loop.vector){
  x <- test1_loop[,i]
  title <- names(test1_loop[i])
  hist(x,
       main = title,
       xlab = i)
}

# Trait distribution of compatible genotypes per strain (data preparation)----
## Removing the discriminant and control
s09AX27_clean <- test1_NA_nodisc[which(test1_NA_nodisc$Strain == "09AX27"),]
s93JE3_clean <- test1_NA_nodisc[which(test1_NA_nodisc$Strain == "93JE3"),]
sP72_clean <- test1_NA_nodisc[which(test1_NA_nodisc$Strain == "P72"),]

## Graphical inspection for each strain
## Strain 09AX27
s09AX27_loop <- data.frame("LP" = s09AX27_clean$Laten_24.02.15_gp,
                           "UN" = s09AX27_clean$total_NbSS,
                           "US" = s09AX27_clean$Taill_05.08.20_gp,
                           "Black_Points" = s09AX27_clean$PtN_real)
loop.vector <- 1:4

for (i in loop.vector){
  x <- s09AX27_loop[,i]
  title <- names(s09AX27_loop[i])
  hist(x,
       main = paste("09AX27", title),
       xlab = "Scores")
}

# Strain 93JE3
s93JE3_loop <- data.frame("LP" = s93JE3_clean$Laten_24.02.15_gp,
                          "UN" = s93JE3_clean$total_NbSS,
                          "US" = s93JE3_clean$Taill_05.08.20_gp,
                          "Black_Points" = s93JE3_clean$PtN_real)

for (i in loop.vector){
  x <- s93JE3_loop[,i]
  title <- names(s93JE3_loop[i])
  hist(x,
       main = paste("93JE3", title),
       xlab = "Scores")
}

# Strain P72
sP72_loop <- data.frame("LP" = sP72_clean$Laten_24.02.15_gp,
                        "UN" = sP72_clean$total_NbSS,
                        "US" = sP72_clean$Taill_05.08.20_gp,
                        "Black_Points" = sP72_clean$PtN_real)

for (i in loop.vector){
  x <- sP72_loop[,i]
  title <- names(sP72_loop[i])
  hist(x,
       main = paste("P72", title),
       xlab = "Scores")
}

# Mean Genotype and Mean Block ----
## Strain 09AX27
str(s09AX27_clean)
s09AX27_clean$Block <- as.factor(s09AX27_clean$Block)
s09AX27_clean$Genotype <- as.factor(s09AX27_clean$Genotype)
summary(s09AX27_clean)

par(mfrow = c(2,2))
plot.design(s09AX27_clean, y = s09AX27_clean$Laten_24.02.15_gp, ylab = "Mean LP")
plot.design(s09AX27_clean, y = s09AX27_clean$total_NbSS, ylab = "Mean UN")
plot.design(s09AX27_clean, y = s09AX27_clean$Taill_05.08.20_gp, ylab = "Mean US")
plot.design(s09AX27_clean, y = s09AX27_clean$PtN_real, ylab = "Mean Black Points")

## Strain 93JE3
str(s93JE3_clean)
s93JE3_clean$Block <- as.factor(s93JE3_clean$Block)
s93JE3_clean$Genotype <- as.factor(s93JE3_clean$Genotype)

plot.design(s93JE3_clean, y = s93JE3_clean$Laten_24.02.15_gp, ylab = "Mean LP")
plot.design(s93JE3_clean, y = s93JE3_clean$total_NbSS, ylab = "Mean UN")
plot.design(s93JE3_clean, y = s93JE3_clean$Taill_05.08.20_gp, ylab = "Mean US")
plot.design(s93JE3_clean, y = s93JE3_clean$PtN_real, ylab = "Mean Black Points")

## Strain P72
sP72_clean$Block <- as.factor(sP72_clean$Block)
sP72_clean$Genotype <- as.factor(sP72_clean$Genotype)

plot.design(sP72_clean, y = sP72_clean$Laten_24.02.15_gp, ylab = "Mean LP")
plot.design(sP72_clean, y = sP72_clean$total_NbSS, ylab = "Mean UN")
plot.design(sP72_clean, y = sP72_clean$Taill_05.08.20_gp, ylab = "Mean US")
plot.design(sP72_clean, y = sP72_clean$PtN_real, ylab = "Mean Black Points")

# Mean Genotype and sd for Laten Period ----
# Function for genotypic mean
mean.laten <- function(Laten_24.02.15_gp, Genotype){
  data.frame(aggregate(Laten_24.02.15_gp, by = Genotype,
                       mean, na.rm = TRUE))
}

# Function for genotypic sd
sd.laten <- function(Laten_24.02.15_gp, Genotype){
  data.frame(aggregate(Laten_24.02.15_gp, by = Genotype,
                       sd, na.rm = TRUE))
}

# Creating list for loop
strain_list <- list(s09AX27_clean, s93JE3_clean, sP72_clean)

mean_laten <- list()
sd_laten <- list()
for (i in 1:length(strain_list)) {
  mean_laten[[i]] <- mean.laten(strain_list[[i]]$Laten_24.02.15_gp,
                                list(strain_list[[i]]$Genotype))
  dat_mean_laten <- as.data.frame(mean_laten)
  sd_laten[[i]] <- sd.laten(strain_list[[i]]$Laten_24.02.15_gp,
                            list(strain_list[[i]]$Genotype))
  dat_sd_laten <- as.data.frame(sd_laten)
}

names(dat_mean_laten) <- c("s09AX27", "Mean_09AX27",
                           "s93JE3", "Mean_93JE3",
                           "sP72", "Mean_P72")
dat_mean_laten[dat_mean_laten == "NaN"] <- NA

names(dat_sd_laten) <- c("s09AX27", "Sd_09AX27",
                         "s93JE3", "Sd_93JE3",
                         "sP72", "Sd_P72")

# Extract the genotypes for each strain
s09AX27_means <- data.frame(
  Genotype = as.character(levels(as.factor(s09AX27_clean$Genotype)))
)

s93JE3_means <- data.frame(
  Genotype = as.character(levels(as.factor(s93JE3_clean$Genotype)))
)

sP72_means <- data.frame(
  Genotype = as.character(levels(as.factor(sP72_clean$Genotype)))
)

# Join them (per strain) with sd and mean
s09AX27_means$LP_mean <- dat_mean_laten$Mean_09AX27
s09AX27_means$LP_sd <- dat_sd_laten$Sd_09AX27

s93JE3_means$LP_mean <- dat_mean_laten$Mean_93JE3
s93JE3_means$LP_sd <- dat_sd_laten$Sd_93JE3

sP72_means$LP_mean <- dat_mean_laten$Mean_P72
sP72_means$LP_sd <- dat_sd_laten$Sd_P72

## Graphical inspection of meanGenotype per strain
## 09AX27
ggplot(s09AX27_means) +
  geom_bar(aes(x = Genotype, y = LP_mean), 
           stat = "identity", 
           fill = "cyan3",
           na.rm = TRUE) +
  geom_errorbar(aes(x = Genotype, ymin = LP_mean - LP_sd,
                    ymax = LP_mean + LP_sd), color = "orange", size = 1) +
  labs(x = "\nGenotypes (09AX27)", y = "Mean Latent Periods (days)\n", 
       title = "Genotypic means across blocks") +
  ylim(0,15) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank())

## 93JE3
ggplot(s93JE3_means) +
  geom_bar(aes(x = Genotype, y = LP_mean), 
           stat = "identity", 
           fill = "cyan3",
           na.rm = TRUE) +
  geom_errorbar(aes(x = Genotype, ymin = LP_mean - LP_sd,
                    ymax = LP_mean + LP_sd), color = "orange", size = 1) +
  labs(x = "\nGenotypes (93JE3)", y = "Mean Latent Periods (days)\n", 
       title = "Genotypic means across blocks") +
  ylim(0,15) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank())

## P72
ggplot(sP72_means) +
  geom_bar(aes(x = Genotype, y = LP_mean), 
           stat = "identity", 
           fill = "cyan3",
           na.rm = TRUE) +
  geom_errorbar(aes(x = Genotype, ymin = LP_mean - LP_sd,
                    ymax = LP_mean + LP_sd), color = "orange", size = 1) +
  labs(x = "\nGenotypes (P72)", y = "Mean Latent Periods (days)\n", 
       title = "Genotypic means across blocks") +
  ylim(0,15) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank())

# Mean Genotype and sd for Uredinia Number ----
# Function for genotypic mean
mean.UN <- function(total_NbSS, Genotype){
  data.frame(aggregate(total_NbSS, by = Genotype,
                       mean, na.rm = TRUE))
}

# Function for genotypic sd
sd.UN <- function(total_NbSS, Genotype){
  data.frame(aggregate(total_NbSS, by = Genotype,
                       sd, na.rm = TRUE))
}

# List for loop = strain_list

mean_UN <- list()
sd_UN <- list()
for (i in 1:length(strain_list)) {
  mean_UN[[i]] <- mean.UN(strain_list[[i]]$total_NbSS,
                          list(strain_list[[i]]$Genotype))
  dat_mean_UN <- as.data.frame(mean_UN)
  sd_UN[[i]] <- sd.UN(strain_list[[i]]$total_NbSS,
                      list(strain_list[[i]]$Genotype))
  dat_sd_UN <- as.data.frame(sd_UN)
}

names(dat_mean_UN) <- c("s09AX27", "UN_Mean_09AX27",
                        "s93JE3", "UN_Mean_93JE3",
                        "sP72", "UN_Mean_P72")
names(dat_sd_UN) <- c("s09AX27", "Sd_09AX27",
                      "s93JE3", "Sd_93JE3",
                      "sP72", "Sd_P72")
# Join them (per strain) with sd and mean
s09AX27_means$UN_mean <- dat_mean_UN$UN_Mean_09AX27
s09AX27_means$UN_sd <- dat_sd_UN$Sd_09AX27

s93JE3_means$UN_mean <- dat_mean_UN$UN_Mean_93JE3
s93JE3_means$UN_sd <- dat_sd_UN$Sd_93JE3

sP72_means$UN_mean <- dat_mean_UN$UN_Mean_P72
sP72_means$UN_sd <- dat_sd_UN$Sd_P72

## Graphical inspection of meanGenotype per strain
## 09AX27
ggplot(s09AX27_means) +
  geom_bar(aes(x = Genotype, y = UN_mean), 
           stat = "identity", 
           fill = "firebrick3",
           na.rm = TRUE) +
  geom_errorbar(aes(x = Genotype, ymin = UN_mean - UN_sd,
                    ymax = UN_mean + UN_sd), color = "orange", size = 1) +
  labs(x = "\nGenotypes (09AX27)", y = "Mean Uredinia Number\n", 
       title = "UN genotypic across blocks") +
  ylim(0,85) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank())

## 93JE3
ggplot(s93JE3_means) +
  geom_bar(aes(x = Genotype, y = UN_mean), 
           stat = "identity", 
           fill = "firebrick3",
           na.rm = TRUE) +
  geom_errorbar(aes(x = Genotype, ymin = UN_mean - UN_sd,
                    ymax = UN_mean + UN_sd), color = "orange", size = 1) +
  labs(x = "\nGenotypes (93JE3)", y = "Mean Uredinia Number\n", 
       title = "UN genotypic means across blocks") +
  ylim(0,85) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank())

## P72
ggplot(sP72_means) +
  geom_bar(aes(x = Genotype, y = UN_mean), 
           stat = "identity", 
           fill = "firebrick3",
           na.rm = TRUE) +
  geom_errorbar(aes(x = Genotype, ymin = UN_mean - UN_sd,
                    ymax = UN_mean + UN_sd), color = "orange", size = 1) +
  labs(x = "\nGenotypes (P72)", y = "Mean Uredinia number\n", 
       title = "UN genotypic means across blocks") +
  ylim(0,85) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank())

# Mean Genotype and sd for Uredinia Size ----
# Function for genotypic mean
mean.US <- function(Taill_05.08.20_gp, Genotype){
  data.frame(aggregate(Taill_05.08.20_gp, by = Genotype,
                       mean, na.rm = TRUE))
}

# Function for genotypic sd
sd.US <- function(Taill_05.08.20_gp, Genotype){
  data.frame(aggregate(Taill_05.08.20_gp, by = Genotype,
                       sd, na.rm = TRUE))
}

# List for loop = strain_list

mean_US <- list()
sd_US <- list()
for (i in 1:length(strain_list)) {
  mean_US[[i]] <- mean.US(strain_list[[i]]$Taill_05.08.20_gp,
                          list(strain_list[[i]]$Genotype))
  dat_mean_US <- as.data.frame(mean_US)
  sd_US[[i]] <- sd.US(strain_list[[i]]$Taill_05.08.20_gp,
                      list(strain_list[[i]]$Genotype))
  dat_sd_US <- as.data.frame(sd_US)
}

names(dat_mean_US) <- c("s09AX27", "US_Mean_09AX27",
                        "s93JE3", "US_Mean_93JE3",
                        "sP72", "US_Mean_P72")
names(dat_sd_US) <- c("s09AX27", "Sd_09AX27",
                      "s93JE3", "Sd_93JE3",
                      "sP72", "Sd_P72")

# Join them (per strain) with sd and mean
s09AX27_means$US_mean <- dat_mean_US$US_Mean_09AX27
s09AX27_means$US_sd <- dat_sd_US$Sd_09AX27

s93JE3_means$US_mean <- dat_mean_US$US_Mean_93JE3
s93JE3_means$US_sd <- dat_sd_US$Sd_93JE3

sP72_means$US_mean <- dat_mean_US$US_Mean_P72
sP72_means$US_sd <- dat_sd_US$Sd_P72

## Graphical inspection of meanGenotype per strain
## 09AX27
ggplot(s09AX27_means) +
  geom_bar(aes(x = Genotype, y = US_mean), 
           stat = "identity", 
           fill = "darkseagreen4",
           na.rm = TRUE) +
  geom_errorbar(aes(x = Genotype, ymin = US_mean - US_sd,
                    ymax = US_mean + US_sd), color = "orange", size = 1) +
  labs(x = "\nGenotypes (09AX27)", y = "Mean Uredinia Size\n", 
       title = "US genotypic means across blocks") +
  ylim(0,5) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank())

## 93JE3
ggplot(s93JE3_means) +
  geom_bar(aes(x = Genotype, y = US_mean), 
           stat = "identity", 
           fill = "darkseagreen4",
           na.rm = TRUE) +
  geom_errorbar(aes(x = Genotype, ymin = US_mean - US_sd,
                    ymax = US_mean + US_sd), color = "orange", size = 1) +
  labs(x = "\nGenotypes (93JE3)", y = "Mean Uredinia Size\n", 
       title = "US genotypic means across blocks") +
  ylim(0,5) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank())

## P72
ggplot(sP72_means) +
  geom_bar(aes(x = Genotype, y = US_mean), 
           stat = "identity", 
           fill = "darkseagreen4",
           na.rm = TRUE) +
  geom_errorbar(aes(x = Genotype, ymin = US_mean - US_sd,
                    ymax = US_mean + US_sd), color = "orange", size = 1) +
  labs(x = "\nGenotypes (P72)", y = "Mean Uredinia Size\n", 
       title = "US genotypic means across blocks") +
  ylim(0,5) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank())

# Mean Genotype and sd for Black Points ----
# Function for genotypic mean
mean.BP <- function(PtN_real, Genotype){
  data.frame(aggregate(PtN_real, by = Genotype,
                       mean, na.rm = TRUE))
}

# Function for genotypic sd
sd.BP <- function(PtN_real, Genotype){
  data.frame(aggregate(PtN_real, by = Genotype,
                       sd, na.rm = TRUE))
}

# List for loop = strain_list

mean_BP <- list()
sd_BP <- list()
for (i in 1:length(strain_list)) {
  mean_BP[[i]] <- mean.BP(strain_list[[i]]$PtN_real,
                          list(strain_list[[i]]$Genotype))
  dat_mean_BP <- as.data.frame(mean_BP)
  sd_BP[[i]] <- sd.BP(strain_list[[i]]$PtN_real,
                      list(strain_list[[i]]$Genotype))
  dat_sd_BP <- as.data.frame(sd_BP)
}

names(dat_mean_BP) <- c("s09AX27", "BP_Mean_09AX27",
                        "s93JE3", "BP_Mean_93JE3",
                        "sP72", "BP_Mean_P72")
dat_mean_BP[dat_mean_BP == "NaN"] <- NA

names(dat_sd_BP) <- c("s09AX27", "Sd_09AX27",
                      "s93JE3", "Sd_93JE3",
                      "sP72", "Sd_P72")

# Join them (per strain) with sd and mean
s09AX27_means$BP_mean <- dat_mean_BP$BP_Mean_09AX27
s09AX27_means$BP_sd <- dat_sd_BP$Sd_09AX27

s93JE3_means$BP_mean <- dat_mean_BP$BP_Mean_93JE3
s93JE3_means$BP_sd <- dat_sd_BP$Sd_93JE3

sP72_means$BP_mean <- dat_mean_BP$BP_Mean_P72
sP72_means$BP_sd <- dat_sd_BP$Sd_P72

## Graphical inspection of meanGenotype per strain
## 09AX27
ggplot(s09AX27_means) +
  geom_bar(aes(x = Genotype, y = BP_mean), 
           stat = "identity", 
           fill = "navy",
           na.rm = TRUE) +
  geom_errorbar(aes(x = Genotype, ymin = BP_mean - BP_sd,
                    ymax = BP_mean + BP_sd), color = "orange", size = 1) +
  labs(x = "\nGenotypes (09AX27)", y = "Mean Black Points\n", 
       title = "BP genotypic means across blocks") +
  ylim(0,15) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank())

## 93JE3
ggplot(s93JE3_means) +
  geom_bar(aes(x = Genotype, y = BP_mean), 
           stat = "identity", 
           fill = "navy",
           na.rm = TRUE) +
  geom_errorbar(aes(x = Genotype, ymin = BP_mean - BP_sd,
                    ymax = BP_mean + BP_sd), color = "orange", size = 1) +
  labs(x = "\nGenotypes (93JE3)", y = "Mean Black Points\n", 
       title = "BP genotypic means across blocks") +
  ylim(0,15) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank())

## P72
ggplot(sP72_means) +
  geom_bar(aes(x = Genotype, y = BP_mean), 
           stat = "identity", 
           fill = "navy",
           na.rm = TRUE) +
  geom_errorbar(aes(x = Genotype, ymin = BP_mean - BP_sd,
                    ymax = BP_mean + BP_sd), color = "orange", size = 1) +
  labs(x = "\nGenotypes (P72)", y = "Mean Black Points\n", 
       title = "BP genotypic means across blocks") +
  ylim(0,15) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank())

# Plotting lab experimental design (strain 09AX27) ----
loop.vector <- 1:4

for (i in loop.vector){
  var <- s09AX27_loop[,i]
  title <- names(s09AX27_loop[i])
  plot <- qplot(X, Y, fill = var, geom = "tile", main = paste("09AX27", title), 
                data = s09AX27_clean, xlab = "X", ylab = "Y") + 
    scale_fill_gradient(low = "green", high = "red") 
  print(plot)
}

# Plotting lab experimental design (strain 93JE3) ----
for (i in loop.vector){
  var <- s93JE3_loop[,i]
  title <- names(s93JE3_loop[i])
  plot2 <- qplot(X, Y, fill = var, geom = "tile", main = paste("93JE3", title), 
                data = s93JE3_clean, xlab = "X", ylab = "Y") + 
    scale_fill_gradient(low = "green", high = "red") 
  print(plot2)
}

# Plotting lab experimental design (strain P72) ----
for (i in loop.vector){
  var <- sP72_loop[,i]
  title <- names(sP72_loop[i])
  plot3 <- qplot(X, Y, fill = var, geom = "tile", main = paste("P72", title), 
                 data = sP72_clean, xlab = "X", ylab = "Y") + 
    scale_fill_gradient(low = "green", high = "red") 
  print(plot3)
}

# Observing block effects with mixed effects model ----
# list: strain_list, 1 = 09AX27, 2 = 93JE3, 3 = P72
LP_mix <- list()
UN_mix <- list()
US_mix <- list()
BP_mix <- list()
for (i in 1:length(strain_list)) {
  LP_mix[[i]] <- remlf90(fixed = Laten_24.02.15_gp ~ 1,
                         random = ~ Genotype + Block,
                         data = strain_list[[i]],
                         method = "ai")
  UN_mix[[i]] <- remlf90(fixed = total_NbSS ~ 1,
                         random = ~ Genotype + Block,
                         data = strain_list[[i]],
                         method = "ai")
  US_mix[[i]] <- remlf90(fixed = Taill_05.08.20_gp ~ 1,
                         random = ~ Genotype + Block,
                         data = strain_list[[i]],
                         method = "ai")
  BP_mix[[i]] <- remlf90(fixed = PtN_real ~ 1,
                         random = ~ Genotype + Block,
                         data = strain_list[[i]],
                         method = "ai")
}
names(LP_mix) <- c("s09AX27", "s93JE3", "sP72")
names(UN_mix) <- c("s09AX27", "s93JE3", "sP72")
names(US_mix) <- c("s09AX27", "s93JE3", "sP72")
names(BP_mix) <- c("s09AX27", "s93JE3", "sP72")

s09AX27_adj_LP <- s09AX27_clean[, -c(6:15)]
s93JE3_adj_LP <- s93JE3_clean[, -c(6:15)]
sP72_adj_LP <- sP72_clean[, -c(6:15)]
adj_LP <- list(s09AX27_adj_LP, s93JE3_adj_LP, sP72_adj_LP)

# Adjusting LP ----
## Functions for LP
fun.coord <- function(X, Y){
  paste(X, Y, sep = "-")
}

intercept <- function(i){
    LP_mix[[i]]$fixed$Intercept[[1]][1, "value"]
}

BLUP.genot.LP <- function(x, i){
  merge(x, data.frame("Genotype" = rownames(LP_mix[[i]]$ranef$Genotype[[1]]),
                      "BLUP_Genotype" = LP_mix[[i]]$ranef$Genotype[[1]]$value),
        by = "Genotype")
}

BLUP.block.LP <- function(x, i){
  merge(x, data.frame("Block" = rownames(LP_mix[[i]]$ranef$Block[[1]]),
                      "BLUP¨.block" = LP_mix[[i]]$ranef$Block[[1]]$value),
        by = "Block")
}

BLUP.resid.LP <- function(x, i){
  merge(x, data.frame("coord" = fun.coord(x$X, x$Y),
                      "BLUP_resid" = residuals(LP_mix[[i]])),
        by = "coord")
}

adj.LP <- function(x){
  x$intercept_LP + x$BLUP_Genotype + x$BLUP_resid
}

# Adjusted LP strain 09AX27 ----
s09AX27_adj_LP$coord <- fun.coord(s09AX27_adj_LP$X, s09AX27_adj_LP$Y)
s09AX27_adj_LP$intercept_LP <- inter.LP("s09AX27")
s09AX27_adj_LP <- BLUP.genot.LP(s09AX27_adj_LP, "s09AX27")
s09AX27_adj_LP <- BLUP.block.LP(s09AX27_adj_LP, "s09AX27")
s09AX27_adj_LP <- BLUP.resid.LP(s09AX27_adj_LP, "s09AX27")
s09AX27_adj_LP$adj_LP <- adj.LP(s09AX27_adj_LP)

# Adjusted LP strain 93JE3 ----
s93JE3_adj_LP$coord <- fun.coord(s93JE3_adj_LP$X, s93JE3_adj_LP$Y)
s93JE3_adj_LP$intercept_LP <- inter.LP("s93JE3")
s93JE3_adj_LP <- BLUP.genot.LP(s93JE3_adj_LP, "s93JE3")
s93JE3_adj_LP <- BLUP.block.LP(s93JE3_adj_LP, "s93JE3")
s93JE3_adj_LP <- BLUP.resid.LP(s93JE3_adj_LP, "s93JE3")
s93JE3_adj_LP$adj_LP <- adj.LP(s93JE3_adj_LP)

# Adjusted LP strain P72 ----
sP72_adj_LP$coord <- fun.coord(sP72_adj_LP$X, sP72_adj_LP$Y)
sP72_adj_LP$intercept_LP <- inter.LP("sP72")
sP72_adj_LP <- BLUP.genot.LP(sP72_adj_LP, "sP72")
sP72_adj_LP <- BLUP.block.LP(sP72_adj_LP, "sP72")
sP72_adj_LP <- BLUP.resid.LP(sP72_adj_LP, "sP72")
sP72_adj_LP$adj_LP <- adj.LP(sP72_adj_LP)
