library(tidyverse)
library(ppmSuite)

# Controls
Nobs <- 200
Minformative <- 5
Muninformative <- 5

# simulate data

# Seed and empty dataframe
df <- tibble("row" = 1:Nobs)

# Create Minformative new variables that are random observations with apply function

for (i in 1:Minformative) {
  df <- df %>%
    mutate(
      !!paste0("Inf", i) := rnorm(Nobs),
      !!paste0("Uninf", i) := rnorm(Nobs)
    )
}

df <- df %>%
  select(-row) %>%
  # sum across just Inf columns
  mutate(pRisk = rowSums(select(., starts_with("Inf")))) %>%
  # Rescale pRisk to be between 0 and 1
  mutate(
    pRisk = (pRisk - min(pRisk)) / (max(pRisk) - min(pRisk)),
    pRisk = pRisk^2
  ) %>%
  # mutate(pRisk = ifelse(pRisk < 0, 0.25, 0.75)) %>%
  # Create binary varaible Risk group based on pRisk using rbinom
  mutate(riskGroups = rbinom(Nobs, 1, pRisk)) %>%
  # Create outcome variable Y that depends on all columns including riskGroups
  mutate(
    # add an interaction term between Inf varaibles and  Risk group that affects Y
    Y = (riskGroups + 1)^3 * rowSums(select(., starts_with("Inf"))) + rowSums(select(., starts_with("Uninf"))) + riskGroups * 10 + rnorm(Nobs),
  )

df %>%
  ggplot(aes(x = pRisk, y = Y, col = riskGroups)) +
  geom_point()
M <- 1e-20
test <- ppmSuite::gaussian_ppmx(y = df$Y, X = df[c(paste0("Inf", 1:Minformative), paste0("Uninf", 1:Muninformative))], draws = 1000, burn = 50, thin = 5, M = M, meanModel = 1, similarity_function = 1, verbose = TRUE)


RHS <- gsub(pattern = ",", replacement = " + ", toString(paste0("Uninf", 1:Muninformative)))

form <- as.formula(paste("Y ~", RHS))

test2df <- df %>%
  # residualize out the Uninf varaibles
  mutate(
    Y = lm(form, data = .)$residuals
  )
test2 <- ppmSuite::gaussian_ppmx(y = test2df$Y, X = test2df[paste0("Inf", 1:Minformative)], draws = 1000, burn = 50, thin = 5, M = M, meanModel = 1, similarity_function = 1, verbose = TRUE)

llike <- rowSums(log(test$like))
mle <- which(max(llike) == llike)
df$predicted <- test$fitted.values[mle, ]
df$label <- test$Si[mle, ]
table(df$label)


llike <- rowSums(log(test2$like))
mle <- which(max(llike) == llike)
df$predicted <- test2$fitted.values[mle, ]
df$label2 <- test2$Si[mle, ]
table(df$label2)


# R2 <- case_when(
#   proj == "False" ~ summary(lm(Z ~ predicted, data = df))$r.squared,
#   proj == "Confound" ~ summary(lm(Z_c ~ predicted, data = df))$r.squared,
#   proj == "C+Xc" ~ summary(lm(Z_xc ~ predicted, data = df))$r.squared
# )
# print(paste("R2 : ", R2))

Metrics::f1(df$riskGroups, df$label)
Metrics::f1(df$riskGroups, df$label2)
