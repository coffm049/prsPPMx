library(ppmSuite)
library(tidyverse)


df <- data.frame("x" = rnorm(500)) %>%
  mutate(
    group = rep(1:2, each = 250),
    Y = x * 3 + group * 2 + rnorm(500) / 10
  ) %>%
  # standardize x and y
  mutate(
    x = (x - mean(x)) / sd(x),
    Y = (Y - mean(Y)) / sd(Y)
  )


M <- 1
simParms <- c(0.0, 1.0, 0.1, 1.0, 2.0, 0.1, 1)
simParms <- c(0.0, 1.0, 1, 1.0, 2.0, 0.1, 1)
test <- ppmSuite::gaussian_ppmx(y = df$Y, X = df["x"], draws = 1000, burn = 50, thin = 5, M = M, meanModel = 2, similarity_function = 1, verbose = TRUE, simParms = simParms)

llike <- rowSums(log(test$like))
mle <- which(max(llike) == llike)
df$predicted <- test$fitted.values[mle, ]
df$label <- test$Si[mle, ]
table(df$label)


Metrics::f1(df$group, df$label)


df %>%
  ggplot(aes(x = x, y = Y, col = factor(label))) +
  geom_point()


M <- 1
simParms <- c(0.0, 1.0, 0.1, 1.0, 2.0, 0.1, 1)
test <- ppmSuite::gaussian_ppmx(
  y = df$Y, X = df["x"], draws = 11000, burn = 1000, thin = 10,
  M = M, meanModel = 2, similarity_function = 1,
  simParms = simParms, verbose = TRUE
)

plot(df$x, apply(test$fitted.values, 2, mean), col = test$Si[1, ])

simParms <- c(0.0, 1.0, 1, 1.0, 2.0, 0.1, 1)
test <- ppmSuite::gaussian_ppmx(
  y = df$Y, X = df["x"], draws = 11000, burn = 1000, thin = 10,
  M = M, meanModel = 2, similarity_function = 1,
  simParms = simParms, verbose = TRUE
)

plot(df$x, apply(test$fitted.values, 2, mean), col = test$Si[1, ])
