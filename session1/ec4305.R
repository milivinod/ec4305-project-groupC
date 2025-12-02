library(ISLR)
library(gsl)
library(copula)
library(VineCopula)
library(fitdistrplus)
library(dplyr)
data(Credit)

## data cleaning

Credit <- Credit %>%
  mutate(
    Gender  = as.integer(Gender  == "Female"),
    Student = as.integer(Student == "Yes"),
    Married = as.integer(Married == "Yes"),
    Asian = as.integer(Ethnicity == "Asian"),
    AfricanAmerican = as.integer(Ethnicity == "African American"),
    Caucasian = as.integer(Ethnicity == "Caucasian")
  ) %>%
  dplyr::select(-c(ID, Ethnicity))

## to remove scientific notations

options(scipen = 999) 

## summary statistics

summary_statistics <- data.frame(
  Mean = sapply(Credit, mean),
  Median = sapply(Credit, median),
  Variance = sapply(Credit, var)
)

summary_statistics <- signif(summary_statistics, 3)
covariance_matrix <- signif(cov(Credit), 3)
correlation_matrix <- signif(cor(Credit), 3)

## copula estimations

## only key financial variables selected
numeric_variables <- Credit %>% dplyr::select(Income, Limit, Rating, Balance)

# fitting normal margins
parametric_marginals <- lapply(numeric_variables, function(x) {
  fitdistrplus::fitdist(x, "norm")
})

# transform to uniform(0,1) distribution
uniform_data <- as.data.frame(mapply(function(x, fit) {
  pnorm(x, mean=fit$estimate["mean"], sd=fit$estimate["sd"])
}, numeric_variables, parametric_marginals))

# robust copula fitting function 
fit_all_copulas <- function(u, v) {
  tau_value <- cor(u, v, method = "kendall")  
  
  fit_models <- list()
  
  # Gaussian copula
  cop_gaussian <- normalCopula(param = sin(pi * tau_value / 2), dim = 2)
  fit_models$Gaussian <- suppressWarnings(
    tryCatch(fitCopula(cop_gaussian, cbind(u,v), method = "ml"), error=function(e) NULL)
  )
  
  # Student t copula
  cop_studentt <- tCopula(param = sin(pi * tau_value / 2), dim = 2, df = 4, df.fixed = FALSE)
  fit_models$StudentT <- suppressWarnings(
    tryCatch(fitCopula(cop_studentt, cbind(u,v), method = "ml"), error=function(e) NULL)
  )
  
  # Clayton (tau > 0)
  if (tau_value > 0) {
    cop_clayton <- claytonCopula(param = iTau(claytonCopula(), tau_value))
    fit_models$Clayton <- suppressWarnings(
      tryCatch(fitCopula(cop_clayton, cbind(u,v), method = "ml"), error=function(e) NULL)
    )
  }
  
  # Gumbel (tau >= 0)
  if (tau_value >= 0) {
    cop_gumbel <- gumbelCopula(param = iTau(gumbelCopula(), tau_value))
    fit_models$Gumbel <- suppressWarnings(
      tryCatch(fitCopula(cop_gumbel, cbind(u,v), method = "ml"), error=function(e) NULL)
    )
  }
  
  # Compile results
  results <- lapply(names(fit_models), function(name) {
    fitmodel <- fit_models[[name]]
    if (is.null(fitmodel)) {
      return(data.frame(
        Copula = name,
        Loglikelihood = NA,
        AIC = NA,
        BIC = NA,
        Params = NA,
        KendallTau = round(tau_value, 4),
        stringsAsFactors = FALSE
      ))
    } else {
      return(data.frame(
        Copula = name,
        Loglikelihood = as.numeric(logLik(fitmodel)),
        AIC = AIC(fitmodel),
        BIC = BIC(fitmodel),
        Params = paste(round(fitmodel@estimate, 4), collapse = ","),
        KendallTau = round(tau_value, 4),
        stringsAsFactors = FALSE
      ))
    }
  })
  
  do.call(rbind, results)
}

# pairwise loop
variable_names <- names(uniform_data)
pairwise_results <- list()
summary_table <- data.frame()

for (i in 1:(length(variable_names) - 1)) {
  for (j in (i+1):length(variable_names)) {
    var1 <- variable_names[i]
    var2 <- variable_names[j]
    
    fit_df <- fit_all_copulas(uniform_data[[var1]], uniform_data[[var2]])
    
    # add pair info
    fit_df$Pair <- paste(var1, var2, sep = "_")
    
    pairwise_results[[paste(var1, var2, sep = "_")]] <- fit_df
    summary_table <- rbind(summary_table, fit_df)
  }
}

summary_table


