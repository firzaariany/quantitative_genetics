summary(LP_mix[[1]])
LP_mix[[i]] <- remlf90(fixed = Laten_24.02.15_gp ~ 1,
                       random = ~ Genotype + Block,
                       data = strain_list[[i]],
                       method = "ai")

LP_09_nr <- remlf90(fixed = Laten_24.02.15_gp ~ Genotype,
                    random =  ~ Block,
                    data = strain_list[[1]],
                    method = "ai")
summary(LP_09_nr)

#calculate the chi-squared stat for the log-likelihood ratio test
2*(LP_mix[[1]]$fit$`-2logL` - LP_09_nr$fit$`-2logL`)

#calculate the associated significance
1 - pchisq(2*(LP_mix[[1]]$fit$`-2logL` - LP_09_nr$fit$`-2logL`), 153)

# Genotype is significant?
LP_09_nr2 <- remlf90(fixed = Laten_24.02.15_gp ~ Block,
                     random =  ~ Genotype,
                     data = strain_list[[1]],
                     method = "ai")

#calculate the associated significance
1 - pchisq(2*(LP_mix[[1]]$fit$`-2logL` - LP_09_nr2$fit$`-2logL`), 4)

# UN 09AX27 ----
# With random = UN_mix[[1]]
# Null model, genotype fixed
UN_09_nr <- remlf90(fixed = total_NbSS ~ Genotype,
                    random =  ~ Block,
                    data = strain_list[[1]],
                    method = "ai")

#calculate the associated significance
1 - pchisq(2*(UN_mix[[1]]$fit$`-2logL` - UN_09_nr$fit$`-2logL`), 153)

# Null model, bloc = fixed
UN_09_nr2 <- remlf90(fixed = total_NbSS ~ Block,
                    random =  ~ Genotype,
                    data = strain_list[[1]],
                    method = "ai")

#calculate the associated significance
1 - pchisq(2*(UN_mix[[1]]$fit$`-2logL` - UN_09_nr2$fit$`-2logL`), 4)

# US 09AX27 ----
# With random = US_mix[[1]]
# Null model, fixed = genotype
US_09_nr <- remlf90(fixed = Taill_05.08.20_gp ~ Genotype,
                    random =  ~ Block,
                    data = strain_list[[1]],
                    method = "ai")

#calculate the associated significance
1 - pchisq(2*(US_mix[[1]]$fit$`-2logL` - US_09_nr$fit$`-2logL`), 153)

# Null model, fixed = block
US_09_nr2 <- remlf90(fixed = Taill_05.08.20_gp ~ Block,
                    random =  ~ Genotype,
                    data = strain_list[[1]],
                    method = "ai")

#calculate the associated significance
1 - pchisq(2*(US_mix[[1]]$fit$`-2logL` - US_09_nr2$fit$`-2logL`), 4)

# LP 93JE3 ----
# With random = LP_mix[[2]]
# Null model, fixed = genotype
LP_93_nr <- remlf90(fixed = Laten_24.02.15_gp ~ Genotype,
                    random =  ~ Block,
                    data = strain_list[[2]],
                    method = "ai")

#calculate the associated significance
1 - pchisq(2*(LP_mix[[2]]$fit$`-2logL` - LP_93_nr$fit$`-2logL`), 153)

# Null model, fixed = block
LP_93_nr2 <- remlf90(fixed = Laten_24.02.15_gp ~ Block,
                    random =  ~ Genotype,
                    data = strain_list[[2]],
                    method = "ai")

#calculate the associated significance
1 - pchisq(2*(LP_mix[[2]]$fit$`-2logL` - LP_93_nr2$fit$`-2logL`), 4)

# UN 93JE3 ----
# With random = UN_mix[[2]]
# Null model, fixed = genotype
UN_93_nr <- remlf90(fixed = total_NbSS ~ Genotype,
                    random =  ~ Block,
                    data = strain_list[[2]],
                    method = "ai")

# Calculate the associated significance
1 - pchisq(2*(UN_mix[[2]]$fit$`-2logL` - UN_93_nr$fit$`-2logL`), 153)

# Null model, fixed = block
UN_93_nr2 <- remlf90(fixed = total_NbSS ~ Block,
                    random =  ~ Genotype,
                    data = strain_list[[2]],
                    method = "ai")

#calculate the associated significance
1 - pchisq(2*(UN_mix[[2]]$fit$`-2logL` - UN_93_nr2$fit$`-2logL`), 4)

# US 93JE3 ----
# With random = US_mix[[2]]
# Null model, fixed = genotype
US_93_nr <- remlf90(fixed = Taill_05.08.20_gp ~ Genotype,
                    random =  ~ Block,
                    data = strain_list[[2]],
                    method = "ai")

# Calculate the associated significance
1 - pchisq(2*(US_mix[[2]]$fit$`-2logL` - US_93_nr$fit$`-2logL`), 153)

# Null model, fixed = block
US_93_nr2 <- remlf90(fixed = Taill_05.08.20_gp ~ Block,
                    random =  ~ Genotype,
                    data = strain_list[[2]],
                    method = "ai")

# Calculate the associated significance
1 - pchisq(2*(US_mix[[2]]$fit$`-2logL` - US_93_nr2$fit$`-2logL`), 4)

# LP P72 ----
# With random = LP_mix[[3]]
# Null model, fixed = genotype
LP_72_nr <- remlf90(fixed = Laten_24.02.15_gp ~ Genotype,
                    random =  ~ Block,
                    data = strain_list[[3]],
                    method = "ai")

# Calculate the associated significance
1 - pchisq(2*(LP_mix[[3]]$fit$`-2logL` - LP_72_nr$fit$`-2logL`), 153)

# Null model, fixed = block
LP_72_nr2 <- remlf90(fixed = Laten_24.02.15_gp ~ Block,
                    random =  ~ Genotype,
                    data = strain_list[[3]],
                    method = "ai")

# Calculate the associated significance
1 - pchisq(2*(LP_mix[[3]]$fit$`-2logL` - LP_72_nr2$fit$`-2logL`), 4)

# UN P72 ----
# With random = UN_mix[[3]]
# Null model, fixed = genotype
UN_72_nr <- remlf90(fixed = total_NbSS ~ Genotype,
                    random =  ~ Block,
                    data = strain_list[[3]],
                    method = "ai")

# Calculate the associated significance
1 - pchisq(2*(UN_mix[[3]]$fit$`-2logL` - UN_72_nr$fit$`-2logL`), 153)

# Null model, fixed = block
UN_72_nr2 <- remlf90(fixed = total_NbSS ~ Block,
                    random =  ~ Genotype,
                    data = strain_list[[3]],
                    method = "ai")

# Calculate the associated significance
1 - pchisq(2*(UN_mix[[3]]$fit$`-2logL` - UN_72_nr2$fit$`-2logL`), 4)

# US P72 ----
# With random = US_mix[[3]]
# Null model, fixed = genotype
US_72_nr <- remlf90(fixed = Taill_05.08.20_gp ~ Genotype,
                    random =  ~ Block,
                    data = strain_list[[3]],
                    method = "ai")

# Calculate the associated significance
1 - pchisq(2*(US_mix[[3]]$fit$`-2logL` - US_72_nr$fit$`-2logL`), 153)

# Null model, fixed = block
US_72_nr2 <- remlf90(fixed = Taill_05.08.20_gp ~ Block,
                    random =  ~ Genotype,
                    data = strain_list[[3]],
                    method = "ai")

# Calculate the associated significance
1 - pchisq(2*(US_mix[[3]]$fit$`-2logL` - US_72_nr2$fit$`-2logL`), 4)
