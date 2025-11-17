#load libraries
library(ISLR)
library(nortest)
library(tseries) 
library(ggplot2)
library(gridExtra)
library(tidyr)
library(tidyverse)
library(vcd)
library(GGally)
library(ggplot2)

#Data
data("Credit")

#convert categorical variables to binary variables (0-1)
Credit <- Credit %>%
  select(-ID) %>%
  mutate(
    Gender = as.integer(Gender == "Female"), 
    Student = as.integer(Student == "Yes"),
    Married = as.integer(Married == "Yes"), 
    Asian = as.integer(Ethnicity == "Asian"), 
    AfricanAmerican = as.integer(Ethnicity == "African American"), 
    Caucasian = as.integer(Ethnicity == "Caucasian")) %>%
    select(-Ethnicity) 

#define variable groups
continuous_vars <- c("Income", "Limit", "Rating", "Age", "Education", "Balance")
count_vars      <- c("Cards")
binary_vars     <- c("Gender", "Student", "Married", "Asian", "AfricanAmerican", "Caucasian")


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
mgf_estimate <- function(x, t_values) {
  sapply(t_values, function(t) mean(exp(t * x)))
}

t_values_1 <- c(-1, -0.5, 0, 0.5, 1) #using this got numerical overflow (Inf)
t_values_2 <- c(-0.01, -0.005, 0, 0.005, 0.01)
mgf_results <- list()

for (var in continuous_vars) {
  x <- Credit[[var]]
  mgf_results[[var]] <- mgf_estimate(x, t_values_2)
}

mgf_results

##PLOT##
par(mfrow=c(2,3), mar=c(4,4,2,1))
for (var in continuous_vars) {
  x <- Credit[[var]]
  hist(x, probability=TRUE, main=paste("Histogram of", var), xlab=var, col="lightgray", border="white")
  lines(density(x), col="blue", lwd=2)
}

#plots confirm that they are not normal


##Scatter plot of the integral transformed data##
pit_data <- as.data.frame(lapply(Credit[continuous_vars], function(x) {
  jitter(rank(x) / (length(x) + 1), factor = 0.01)
}))

col_points <- rgb(0, 0, 1, 0.5) 
pairs(pit_data,
      main = "Scatterplot Matrix of PIT-Transformed Data",
      pch = 19,          # solid dots
      cex = 0.5,         # smaller point size
      col = col_points)
