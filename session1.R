#load libraries
library(ISLR)
library(nortest)
library(tseries) 
library(ggplot2)
library(gridExtra)
library(tidyr)
library(tidyverse)
library(vcd)
library(gsl)
library(copula)
library(VineCopula)
library(fitdistrplus)
library(dplyr)
data(Credit)

pkgs <- c("ISLR", "nortest", "tseries", "ggplot2", "gridExtra", "tidyr",
  "tidyverse", "vcd", "gsl", "copula", "VineCopula",
  "fitdistrplus", "dplyr"
)

# Install any that are missing
to_install <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if(length(to_install) > 0) install.packages(to_install)

# Load all packages
lapply(pkgs, library, character.only = TRUE)


#convert categorical variables to binary variables (0-1)
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

#define variable groups
continuous_vars <- c("Income", "Limit", "Rating", "Balance")
count_vars      <- c("Cards")
binary_vars     <- c("Gender", "Student", "Married", "Asian", "AfricanAmerican", "Caucasian")


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


##DISTRIBUTION TESTS##
results <- list()

#for continuous variables
for (var in continuous_vars) {
  x <- Credit[[var]]
  
  results[[var]] <- list(
    Shapiro_Wilk = shapiro.test(x),
    Anderson_Darling = ad.test(x),
    Jarque_Bera  = jarque.bera.test(x),
    Lilliefors       = lillie.test(x)
  )
}

#for discrete variable (Cards)
for (var in count_vars) {
  x <- Credit[[var]]
  
  results[[var]] <- list(
    Poisson_GOFA = goodfit(x, type = "poisson")
  )
}


#for binary 0-1 variables 
for (var in binary_vars) {
  x <- Credit[[var]]
  
  results[[var]] <- list(
    Binomial_Test = binom.test(sum(x), length(x)),
    Chi_Square_Test = chisq.test(table(x))
  )
}

#print the distribution test results for continuous variables
dist_tests <- data.frame(
  Variable = character(),
  Shapiro_Wilk_p = numeric(),
  Anderson_Darling_p = numeric(),
  Lilliefors_p = numeric(),
  Jarque_Bera_p = numeric(),
  stringsAsFactors = FALSE
)

for (var in continuous_vars) {
  x <- results[[var]]
  
  dist_tests <- rbind(dist_tests, data.frame(
    Variable = var,
    Shapiro_Wilk_p     = x$Shapiro_Wilk$p.value,
    Anderson_Darling_p = x$Anderson_Darling$p.value,
    Lilliefors_p       = x$Lilliefors$p.value,
    Jarque_Bera_p      = x$Jarque_Bera$p.value
  ))
}

# Show the results neatly
print(dist_tests) #all of them reject null hypothesis -> non-normal

##MGF Estimates##
mgf_hat <- function(x, t) {
  mean(exp(t * x))
}

#t_values_1 <- c(-1, -0.5, 0, 0.5, 1) using this got numerical overflow (Inf)

t_vals <- seq(-0.001, 0.001, length.out = 200)

#mgf estimates for each continuous variable
mgf_income  <- sapply(t_vals, function(t) mgf_hat(Credit$Income, t))
mgf_limit <- sapply(t_vals, function(t) mgf_hat(Credit$Limit, t))
mgf_rating <- sapply(t_vals, function(t) mgf_hat(Credit$Rating, t))
mgf_balance <- sapply(t_vals, function(t) mgf_hat(Credit$Balance, t))


# Plot mgf estimates
plot(t_vals, mgf_income, type = "l", 
     main = "Estimated MGF: Income", xlab = "t", ylab = "M_X(t)")
plot(t_vals, mgf_limit, type = "l", 
     main = "Estimated MGF: Limit", xlab = "t", ylab = "M_X(t)")
plot(t_vals, mgf_rating, type = "l", 
     main = "Estimated MGF: Rating", xlab = "t", ylab = "M_X(t)")
plot(t_vals, mgf_balance, type = "l", 
     main = "Estimated MGF: Balance", xlab = "t", ylab = "M_X(t)")


##PLOT##
png(filename = "Univariate Empirical Distributions.png", 
    width = 1600, height = 1000, res = 150)

par(mfrow=c(2,2), mar=c(4,4,4,1))
for (var in continuous_vars) {
  x <- Credit[[var]]
  hist(x, probability=TRUE, main=paste("Histogram of", var), xlab=var, col="lightgray", border="white")
  lines(density(x), col="blue", lwd=2)
}
dev.off()

#plots confirm that they are not normal


##Scatter plot of the integral transformed data##

png(filename = "Integral_Transformed_Scatter.png", 
    width = 1600, height = 1000, res = 150)
pairs_list <- list(
  c("Income", "Limit"), 
  c("Income", "Rating"), 
  c("Income", "Balance"), 
  c("Limit", "Rating"), 
  c("Limit", "Balance"), 
  c("Rating", "Balance"))

# Set layout: 2 rows x 3 columns
par(mfrow = c(2,3), mar = c(4,4,2,1))

# Loop through each pair and plot
for(pair in pairs_list){
  X <- Credit[[pair[1]]]
  Y <- Credit[[pair[2]]]
  
  U1 <- ecdf(X)(X)
  U2 <- ecdf(Y)(Y)
  
  plot(U1, U2,
       main = paste(pair[1], "vs", pair[2]),
       xlab = paste("U =", pair[1]),
       ylab = paste("U =", pair[2]),
       pch = 20, col = rgb(0,0,1,0.4))
}
dev.off()
