options(tidyverse.quiet = TRUE, readr.show_col_types = FALSE)
library(tidyverse)


loadSim <- function(filepath = "temp/simulation") {

    df <-  paste0(filepath, ".Y", seq(0, 4), ".best") |>
      map(read_table) |>
      purrr::reduce(full_join, by = c("FID", "IID")) |>
      select(!contains("Regression"))
    colnames(df) <- c("FID", "IID", "PRS0", "PRS1", "PRS2", "PRS3", "PRS4")
    pcs <- read_table(paste0(filepath, ".eigenvec"))
    covars <- read_table(paste0(filepath,".covar"))
    outcome <- read_table(paste0(filepath, ".pheno"))

    df <- full_join(df, pcs, by = c("FID" = "#FID", "IID")) |>
      full_join(covars, by = c("FID", "IID")) |>
      full_join(outcome, by = c("FID", "IID")) |>
      mutate(PRS0 = (PRS0 - mean(PRS0)) / sd(PRS0), 
             PRS1 = (PRS1 - mean(PRS1)) / sd(PRS1),
             PRS2 = (PRS2 - mean(PRS2)) / sd(PRS2),
             PRS3 = (PRS3 - mean(PRS3)) / sd(PRS3),
             PRS4 = (PRS4 - mean(PRS4)) / sd(PRS4),
             Z = (Z- mean(Z)) / sd(Z), 
             PC1 = (PC1 - mean(PC1) / sd(PC1)), 
             Z = Z + riskGroups * 4,
             # Create a covariate that has probability 0.25 if riskGroup == 0 
             # and 0.75 if riskGroup == 1
             Confound = ifelse(
              riskGroups == 0,
              sample(c(0, 1), size= 1, prob = c(0.75, 0.25)),
              sample(c(0, 1), size= 1, prob = c(0.25, 0.75))
             ), 
             Z = Z + Confound * 3,
             PRS0_c = lm(PRS0 ~ Confound)$residuals,
             PRS1_c = lm(PRS1 ~ Confound)$residuals,
             PRS2_c = lm(PRS2 ~ Confound)$residuals,
             PRS3_c = lm(PRS3 ~ Confound)$residuals,
             PRS4_c = lm(PRS4 ~ Confound)$residuals,
             Z_c = lm(Z ~ Confound)$residuals,
             PRS0_xc = lm(PRS0 ~ Xc + Confound)$residuals,
             PRS1_xc = lm(PRS1 ~ Xc + Confound)$residuals,
             PRS2_xc = lm(PRS2 ~ Xc + Confound)$residuals,
             PRS3_xc = lm(PRS3 ~ Xc + Confound)$residuals,
             PRS4_xc = lm(PRS4 ~ Xc + Confound)$residuals,
             Z_xc = lm(Z ~ Xc + Confound)$residuals,
             )
  return(df)
}


ppmxsummary <- function(mod, df, proj="False", ...) {
  
  llike <- rowSums(log(mod$like))
  mle <- which(max(llike) == llike)
  
  df$predicted <- mod$fitted.values[mle,] 
  df$label <- mod$Si[mle,]
  
  #g <- df %>%
  #  ggplot(aes(!!!ensyms(...))) +
  #  geom_point() +
  #  geom_point(aes(y = predicted), color = "red") 
  #print(g)
 
  # print(table(df[c("label", "riskGroups", "subj_ancestries")]))

  #print(paste("f1 : ", Metrics::f1(df$riskGroups, df$label)))
  R2 <- case_when(
    proj == "False" ~ summary(lm(Z ~ predicted, data = df))$r.squared,
    proj == "Confound" ~ summary(lm(Z_c ~ predicted, data = df))$r.squared,
    proj == "C+Xc" ~ summary(lm(Z_xc ~ predicted, data = df))$r.squared
    )
  #print(paste("R2 : ", R2))

  return(data.frame("f1" = Metrics::f1(df$riskGroups, df$label), "R2" = R2, "proj" = proj)) 
}


runSims <- function(meanModel = 1, M=1e-20, similarity_function = 1, draws = 2000, burn = 500, thin = 10, out="temp/simulation"){
    df <- loadSim(filepath = out)
    m1 = ppmSuite::gaussian_ppmx(y = df$Z, X = df[c("PRS0", "PRS1", "PRS2", "PRS3", "PRS4", "Confound", "Xc")], draws = draws, burn = burn, thin = thin, M = M, meanModel= meanModel,similarity_function=similarity_function )
    m2 = ppmSuite::gaussian_ppmx(y = df$Z_c, X = df[c("PRS0_c", "PRS1_c", "PRS2_c","PRS3_c", "PRS4_c", "Xc")], meanModel = meanModel, M = M, similarity_function = similarity_function, draws = draws, burn = burn, thin = thin)
    m3 = ppmSuite::gaussian_ppmx(y = df$Z_xc, X = df[c("PRS0_xc", "PRS1_xc", "PRS2_xc","PRS3_xc", "PRS4_xc")], meanModel = meanModel, M = M, similarity_function = similarity_function, draws = draws, burn = burn, thin = thin)

    m1 <- ppmxsummary(m1, df,proj = "False", x ="PC1", y= "Z")
    m2 <- ppmxsummary(m2, df, x ="PRS1_c", y= "Z_c", proj = "Confound")
    m3 <- ppmxsummary(m3, df, x ="PRS1_xc", y= "Z_xc", proj ="C+Xc")

    results <- rbind(m1, m2, m3)
    write.table(results, file = paste0(out, ".results"), sep = "\t", row.names = F, append = T)
}

# Make an argparser to parse the filepath flag

argparser <- argparse::ArgumentParser(prog = "ppmx",
                                      description = "Run a simulation of the PPMX model")
argparser$add_argument("--filepath", help = "The path to the simulation files", type = "character")
args <- argparser$parse_args()
runSims()



