library(ppmSuite)
data(bear)
# plot length, sex, and weight of bears
ck <- c(4, 3, 2)
pairs(bear[, ck])
# response is weight
Y <- bear$weight
# Continuous Covariate is length of chest
# Categorical covariate is sex
X <- bear[, c("length", "sex")]
X$sex <- as.factor(X$sex)
# Randomly partition data into 44 training and 10 testing
set.seed(1)
trainObs <- sample(1:length(Y), 44, replace = FALSE)
Ytrain <- Y[trainObs]
Ytest <- Y[-trainObs]
Xtrain <- X[trainObs, , drop = FALSE]
Xtest <- X[-trainObs, , drop = FALSE]
simParms <- c(0.0, 1.0, 0.1, 1.0, 2.0, 0.1)
modelPriors <- c(0, 100^2, 0.5 * sd(Y), 100)
M <- 1.0
niter <- 100000
nburn <- 50000
nthin <- 50
nout <- (niter - nburn) / nthin
mh <- c(1, 10)
# Run MCMC algorithm for Gaussian PPMx model
out1 <- gaussian_ppmx(
  y = Ytrain, X = Xtrain, Xpred = Xtest,
  M = M, PPM = FALSE,
  meanModel = 1,
  similarity_function = 1,
  consim = 1,
  calibrate = 0,
  simParms = simParms,
  modelPriors = modelPriors,
  draws = niter, burn = nburn, thin = nthin,
  mh = mh
)
