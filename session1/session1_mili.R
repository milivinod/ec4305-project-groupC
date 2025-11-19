#load libraries
library(ISLR)
library(nortest)
library(tseries) 
library(ggplot2)
library(gridExtra)
library(tidyr)
library(tidyverse)
library(vcd)

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
continuous_vars <- c("Income", "Limit", "Rating", "Balance")
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

###RETRYING WITH LOG TRANSFORMATION on right skewed variables###
log_vars <- c("Income", "Limit", "Rating", "Education", "Balance")

# Create a log-transformed dataset
Credit_log <- Credit
for (var in log_vars) {
  Credit_log[[var]] <- log(Credit[[var]])  
}