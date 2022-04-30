# ----- LIBRARIES -------
install.packages('comprehenr')
install.packages("forecast")
install.packages('latex2exp')
install.packages('tseries')
install.packages("viridis") 
library(comprehenr)
library(matlib)
library(ggplot2)
library(forecast)
library(latex2exp)
library(tseries)
library(viridis)
library(RColorBrewer)

# ---------------- Data prepping -------------------
# Loading the data
data <- data.frame(read.csv(file.choose()))[,1] # interest
data

# Mean and length
m <- mean(data)
N <- length(data)

# Mean corrected data
mean_corr_data <- data - m
mean_corr_data

# Dividing data in training and testing sets
train_data <- mean_corr_data[1:201]
test_data <- mean_corr_data[202:302]

#Plotting data
data_plot <- ggplot(data.frame(mean_corr_data), aes(x=1:N, y=mean_corr_data)) + 
  geom_line(cex = 0.8, col="#5c5c5c") + geom_hline(yintercept = 0) + 
  ylab("Interest") + xlab("Index") + ggtitle("Data Plot")
data_plot

# ----------------- FUNCTIONS ------------------------
# Sample ACVF
samp_acvf <- function(data, maxh){
  vals <- c()
  m <- mean(data); n <- length(data)
  for (h in 0:maxh){
    s <- to_vec(for (i in 1:(n-h)) (data[i]-m)*(data[i+h]-m))
    vals <- c(vals, (1/n)*sum(s))
  }
  vals
}

# Sample ACF
samp_acf <- function(data, maxh){
  vals <- samp_acvf(data, maxh)
  vals / vals[1]   # Divide by phi(0) to obtain ACF values
}

# Aux function to calculate the ACVF matrix (used for the sample PACF function)
acvf_matrix <- function(data, h){
  v <- samp_acvf(data, h)
  toeplitz(v[-length(v)])
}

# Sample PACF
samp_pacf <- function(data, maxh){
  vals <- c(1)
  acvf_mat <- acvf_matrix(data, maxh)
  acvf_vec <- samp_acvf(data, maxh)
  for (h in 1:maxh){
    m1 <- acvf_mat[(1:h), (1:h)]
    m2 <- acvf_vec[2:(h+1)]
    if (h != 1){
      phi <- inv(m1) %*% m2
    } else {
      phi <- (1/m1) * m2
    }
    vals <- c(vals, phi[length(phi)])
  }
  vals
}

# Calculating the ACVF of AR(1) using only phi and sigma
# where, phi and sigma are s.t. X_t =   X_{t-1}*phi + Z_t, Z_t ~ IID(0,sigma^2)
parametric_acvf <- function(phi, sigma, max_lag){
  vals <- c(sigma^2/(1-phi^2)) # acvf(0)
  for (h in 1:max_lag){
    vals <- c(vals, phi*vals[length(vals)])
  }
  vals
}

# Calculating the ACVF matrix using previous function
parametric_acvf_mat <- function(phi, sigma, max_lag){
  v <- parametric_acvf(phi, sigma, max_lag)
  toeplitz(v[-length(v)])
}

# ----------------- QUESTION 2 -----------------------
# Calculating the ACF and PACF values
max_lag <- 20
acf_vals <- samp_acf(mean_corr_data, max_lag)
acf_vals
pacf_vals <- samp_pacf(mean_corr_data, max_lag)
pacf_vals


# With the automatic functions (same values for: [X] ACF; [X] PACF) CHECK!
auto_acf_vals <- Acf(mean_corr_data, lag.max = max_lag)
auto_acf_vals
auto_pacf_vals <- Pacf(mean_corr_data, lag.max = max_lag)
auto_pacf_vals

# Plotting the ACF and PACF 
bounds <- c(-1.96/sqrt(N), 1.96/sqrt(N))

acf_plot <- ggplot(data.frame(acf_vals), aes(x=0:max_lag, acf_vals)) +
  geom_hline(yintercept = bounds, linetype='longdash', color='#ad0000', cex=1.2) +
  geom_segment(aes(x = 0:max_lag, xend = 0:max_lag, y=0, yend=acf_vals), cex=0.9, color = '#74767d') +
  geom_point(pch=19,cex=2,col="#003f5c") + geom_hline(yintercept = 0) +
  ylab("ACF") + xlab("Lag") + ggtitle("Sample ACF of Mean Corrected Data")
acf_plot

pacf_plot <- ggplot(data.frame(pacf_vals), aes(x=0:max_lag, pacf_vals)) + 
  geom_hline(yintercept = bounds, linetype='longdash', color='#ad0000', cex=1.2) +
  geom_segment(aes(x = 0:max_lag, xend = 0:max_lag, y=0, yend=pacf_vals), cex=0.9, color = '#74767d') +
  geom_point(pch=19,cex=2,col="#003f5c") + geom_hline(yintercept = 0) +
  ylab("PACF") + xlab("Lag") + ggtitle("Sample PACF of Mean Corrected Data")
pacf_plot


# ------------------ QUESTION 3 ---------------------
# ----- 3a) -----
# Calculating the pairs
pairs <- matrix(to_vec(for (i in 1:(N-1)) c(mean_corr_data[i], mean_corr_data[i+1])), ncol = 2, byrow = TRUE)
pairs 

# Plotting the pairs
pairs_plot <- ggplot(data.frame(pairs), aes(x=pairs[,1], y=pairs[,2], color = c(1:length(pairs[,1])))) +
  scale_colour_gradient(low="#0b03ff", high="#ffc403") + labs(color="Time t") +
  coord_fixed(ratio = 1) + geom_point(pch=19, cex=1.25) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  ylab("X_{t+1}") + xlab("X_t") + ggtitle("Plot of consecutive points")
pairs_plot

# The pair of observations seem indeed to fall under a straight line, 
# but a simple linear fit can't accurately predict an increasing or decreasing 
# tendency of the time series (otherwise the fit suggests a nearly constant TS)

# ----- 3b) -----
# Fitting the linear regression model (with no intercept)
lin_reg <- lm(pairs[,2] ~ pairs[,1] + 0, data = data.frame(pairs))
summary(lin_reg)

# phi and sigma coefficients
phi <- lin_reg$coefficients[[1]]
sigma <- sigma(lin_reg)

# Plotting the linear fit in the same graph as the pairs
pairs_lr_plot <- pairs_plot + 
  geom_abline(intercept = 0, slope = phi, color="#ff4603", size=1.15)
pairs_lr_plot

# Predictive data points
predicted <- predict(lin_reg, newdata=data.frame(pairs[,1]))
predicted

# Plotting random paths along with the mean corrected data
# Creating the data frame
n_paths <- 3;
plotpaths <- data.frame(matrix(data = rep(mean_corr_data[1], times=n_paths+1), nrow=N*(n_paths+1), ncol=3, byrow=FALSE))
colnames(plotpaths) <- c("Index", "Data", "Interest") 
plotpaths[,1] <- rep(1:N, n_paths+1)
plotpaths[,2] <- rep(c(paste0("Gen. Path ",1:n_paths), "Original"), each = N)
for (p in 2:N){
  prev_inds <- plotpaths['Index'] == (p-1); next_inds <- plotpaths['Index'] == p;
  noise <- rnorm(n_paths, mean = 0, sd = sigma)
  prev_row <- plotpaths[prev_inds, 3][-(n_paths+1)]
  add_row <- c((prev_row * phi) + noise, mean_corr_data[p])
  plotpaths[next_inds, 3] <- add_row
}
plotpaths

# Plotting 
paths_plot <- ggplot(plotpaths, aes(x=Index, y=Interest, group=Data, color=Data)) +
  geom_line(cex = 0.8) + geom_hline(yintercept = 0) + ggtitle("Mean corrected data with simulated paths") +
  scale_colour_manual(values = c(heat.colors(n_paths), "black"))
paths_plot


# ------------ 3c) ----------------
# Calculating the residuals (Z_t = S_t - phi*S_{t-1}, t=2:N)
residuals <- pairs[,2] - predicted
sigma_sq <- sum(residuals^2)/length(residuals)
sigma_sq # We get sigma^2 ~ sigma_sq, as expected

# Plotting the residuals
residuals_plot <- ggplot(data.frame(residuals), aes(x=1:length(residuals), y=residuals)) + 
  geom_point(pch=19, cex=1.2, col="#003f5c") + geom_hline(yintercept = 0) + 
  ylab("Residual") + xlab("Index") + ggtitle("Plot of residuals")
residuals_plot

# ACF of residuals
max_lag <- floor(N/4); n_res <- length(residuals)
acf_residuals <- samp_acf(residuals, max_lag)
acf_residuals

auto_acf_res <- Acf(residuals, lag.max = max_lag) # Matching!

# Plotting the residuals ACF
acf_res_plot <- ggplot(data.frame(acf_residuals), aes(x=0:max_lag, acf_residuals)) + 
  geom_hline(yintercept = bounds, linetype='longdash', color='#ad0000', cex=1.2) +
  geom_segment(aes(x = 0:max_lag, xend = 0:max_lag, y=0, yend=acf_residuals), cex=0.9, color = '#74767d') +
  geom_point(pch=19,cex=2,col="#003f5c") + geom_hline(yintercept = 0) +
  ylab("PACF") + xlab("Lag") + ggtitle("Sample PACF of Residuals")
acf_res_plot
# Less than 95% of the observations fall within the bounds - immediate indicator
# that it is unlikely to be iid noise.

# Ljung-Box test
pacf_max <- samp_pacf(residuals, max_lag)
add <- 0; test_stat <- c(); pvals <- c(); auto_pvals <- c()
for (h in 1:max_lag){
  # Manually done
  add <- add + ((n_res^2 + 2*n_res) * (pacf_max[h+1]^2 / (n_res - h)))
  test_stat <- c(test_stat, add)
  pvals <- c(pvals, 1 - pchisq(q = add, df = h))
  # Automatically done
  ljungbox_test <- Box.test(residuals, lag = h, type = 'Ljung-Box')
  auto_pvals <- c(auto_pvals, ljungbox_test$p.value)
}
pvals_df <- data.frame(cbind(c(1:max_lag), test_stat, pvals, auto_pvals))
colnames(pvals_df) <- c("Lag", "Test Stat", "Manual p-value", "Automatic p-value")
pvals_df
# All the calculated p-vals are equal or extremely close to zero, meaning there is
# very strong evidence that supports the ALTERNATIVE hypothesis, that is,
# E_i does not behave like iid noise


# -------------- QUESTION 4 -------------------
# Calculating sample ACVF from train data
max_lag <- 20
train_acvf <- samp_acvf(train_data, max_lag)
train_acvf_mat <- acvf_matrix(train_data, max_lag)

# Computing coefficients for best linear predictor
coeffs <- solve(train_acvf_mat, train_acvf[2:(max_lag+1)])

# Computing predictions
past_data <- train_data[(length(train_data)-max_lag+1):length(train_data)]
predicted_data <- c()
for (i in 1:length(test_data)){
  pred <- sum(coeffs * rev(past_data))
  predicted_data <- c(predicted_data, pred)
  past_data <- c(past_data[-1], test_data[i])
}
predicted_data

# Plotting predictions with test data
# Creating the data frame
pred_test_data <- data.frame(cbind(rep(c(1:length(test_data))),
                                   rep(c("Predicted", "Test"), each=length(test_data)),
                                   c(predicted_data, test_data))) 
colnames(pred_test_data) <- c("Index", "Data", "Value")
pred_test_data[,1] <- as.numeric(pred_test_data[,1]) # Indexes were being read as chars
pred_test_data[,3] <- as.numeric(pred_test_data[,3]) # Values were being read as chars
pred_test_data
# Plot
pt_plot <- ggplot(pred_test_data, aes(x=Index, y=Value, group=Data, color=Data)) +
  geom_hline(yintercept=0) + geom_line(cex = 0.8) + 
  scale_colour_manual(values = c("#13a600","#0285c2")) +
  xlab("Index") + ylab("Interest") + ggtitle("Predicted and Test data")
pt_plot

# ---- Error calculation ----
# Checking if the predictions are reasonable, i.e looking at the 20 historical data points used.
test_plot1 <- ggplot(data.frame(mean_corr_data[182:201]), aes(x=182:201,y=mean_corr_data[182:201])) + 
  geom_point(pch=19,cex=2,col="#003f5c") + geom_hline(yintercept = 0) +
  ylab("Interest") + xlab("Index") + ggtitle("The last 20 datapoints of the training set")
test_plot1

# Compute the error
pred_error_list <- (predicted_data - test_data)^2
pred_error <- (1/101)*sum(pred_error_list)
pred_error

# Compute the naive "prediction" error
naive_error <- (1/101)*sum(test_data^2)
naive_error

# Computing the square errors for each observation
SE <- data.frame(cbind(rep(c(1:length(test_data))),
                       rep(c("Square error prediction", "Square error from mean"), each=length(test_data)),
                       c(pred_error_list, naive_error_list))) 
colnames(SE) <- c("Index", "Data", "Value")
SE[,1] <- as.numeric(SE[,1]) # Indexes were being read as chars
SE[,3] <- as.numeric(SE[,3]) # Values were being read as chars
SE
# Plot
SE_plot <- ggplot(SE, aes(x=Index, y=Value, group=Data, color=Data)) +
  geom_hline(yintercept=0) + geom_line(cex = 0.8) + 
  scale_colour_manual(values = c("#13a600","#0285c2")) +
  xlab("Month") + ylab("Square Error") + ggtitle("Square Errors")
SE_plot

# Testing for errors in the first 12 months
pred_error_12 <- 1/12*sum(pred_error_list[1:12])
pred_error_12
naive_error_12 <- 1/12*sum(naive_error_list[1:12])
naive_error_12

cum_MSE <- data.frame(cbind(rep(c(1:length(test_data))),
                            rep(c("Square error prediction", "Square error from mean"), each=length(test_data)),
                            c(cumsum(pred_error_list)/101, cumsum(naive_error_list)/101))) 
colnames(cum_MSE) <- c("Index", "Data", "Value")
cum_MSE[,1] <- as.numeric(cum_MSE[,1]) # Indexes were being read as chars
cum_MSE[,3] <- as.numeric(cum_MSE[,3]) # Values were being read as chars
cum_MSE
# Plot
cum_MSE_plot <- ggplot(cum_MSE, aes(x=Index, y=Value, group=Data, color=Data)) +
  geom_hline(yintercept=0) + geom_line(cex = 0.8) + 
  scale_colour_manual(values = c("#13a600","#0285c2")) +
  xlab("Month") + ylab("Cumulative MSE") + ggtitle("Cumulative Mean Square Errors")
cum_MSE_plot

# Is the forecast better?

# On average, the prediction is considerably worse than just using the mean.
# However it is clear that the early predictions, you can graphically determine
# that the nine first predictions beat the average.
# Since the predictive data comes from a higher period than average and
# then actually dives into the lowest period of the whole period it is
# not strange that the predictions are worse than the average.

# -------------------- QUESTION 5 -----------------------
# Computing coefficients for best linear predictor
max_lag <- 20
par_acvf <- parametric_acvf(phi, sigma, max_lag)
par_acvf_mat <- parametric_acvf_mat(phi, sigma, max_lag)
coeffs_q5 <- solve(par_acvf_mat, par_acvf[2:(max_lag+1)])
coeffs_q5

# We get that the best linear predictor of X_{t+1} is simply = phi*X_t, as expected!
# Computing the predicted data
past_data_q5 <- train_data[(length(train_data)-max_lag+1):length(train_data)]
predicted_data_q5 <- c()
for (i in 1:length(test_data)){
  pred <- sum(coeffs_q5 * rev(past_data_q5))
  predicted_data_q5 <- c(predicted_data_q5, pred)
  past_data_q5 <- c(past_data_q5[-1], test_data[i])
}
predicted_data_q5

# Plotting predictions with test data
# Creating the data frame
pred_test_data_q5 <- data.frame(cbind(rep(c(1:length(test_data))),
                                      rep(c("Predicted", "Test"), each=length(test_data)),
                                      c(predicted_data_q5, test_data))) 
colnames(pred_test_data_q5) <- c("Index", "Data", "Value")
pred_test_data_q5[,1] <- as.numeric(pred_test_data_q5[,1]) # Indexes were being read as chars
pred_test_data_q5[,3] <- as.numeric(pred_test_data_q5[,3]) # Values were being read as chars
pred_test_data_q5
# Plot
pt_plot_q5 <- ggplot(pred_test_data_q5, aes(x=Index, y=Value, group=Data, color=Data)) +
  geom_hline(yintercept=0) + geom_line(cex = 0.8) + 
  scale_colour_manual(values = c("#13a600","#0285c2")) +
  xlab("Index") + ylab("Interest") + ggtitle("Predicted and Test data")
pt_plot_q5

# Plotting the data with both predictions from Q4 and Q5
pt_data_q45 <- data.frame(cbind(rep(c(1:length(test_data))),
                                rep(c("Non-parametric", "Parametric", "Test"), each=length(test_data)),
                                c(predicted_data, predicted_data_q5, test_data))) 
colnames(pt_data_q45) <- c("Index", "Data", "Value")
pt_data_q45[,1] <- as.numeric(pt_data_q45[,1]) # Indexes were being read as chars
pt_data_q45[,3] <- as.numeric(pt_data_q45[,3]) # Values were being read as chars
pt_data_q45
# Plot
pt_plot_q45 <- ggplot(pt_data_q45, aes(x=Index, y=Value, group=Data, color=Data)) +
  geom_hline(yintercept=0) + geom_line(cex = 0.8) + 
  scale_colour_manual(values = c("#fc2e17", "#13a600","#0285c2")) +
  xlab("Index") + ylab("Interest") + ggtitle("Predicted and Test data")
pt_plot_q45

# ----- Calculating errors -----
# Compute the best linear predictor error
pred_error_list_q5 <- (predicted_data_q5 - test_data)^2
pred_error_q5 <- (1/101)*sum(pred_error_list_q5)
pred_error_q5

# Compute the naive "prediction" error (same as last question)
naive_error_list <- test_data^2
naive_error <- (1/101)*sum(test_data^2)
naive_error

# Computing the MSE and cumulative MSE for both parametric and non-parametric 
# predictions along with the "naive" prediction of just using the mean
SE_q5 <- data.frame(cbind(rep(c(1:length(test_data))),
                       rep(c("Non-Parametric", "Parametric", "Mean"), each=length(test_data)),
                       c(pred_error_list, pred_error_list_q5, naive_error_list))) 
colnames(SE_q5) <- c("Index", "Data", "Value")
SE_q5[,1] <- as.numeric(SE_q5[,1]) # Indexes were being read as chars
SE_q5[,3] <- as.numeric(SE_q5[,3]) # Values were being read as chars
SE_q5
# Plot
SE_plot_q5 <- ggplot(SE_q5, aes(x=Index, y=Value, group=Data, color=Data)) +
  geom_hline(yintercept=0) + geom_line(cex = 0.8) + 
  scale_colour_manual(values = c("#fc2e17", "#13a600","#0285c2")) +
  xlab("Month") + ylab("Square Error") + ggtitle("Square Errors")
SE_plot_q5

cum_MSE_q5 <- data.frame(cbind(rep(c(1:length(test_data))),
                            rep(c("Non-Parametric", "Parametric", "Mean"), each=length(test_data)),
                            c(cumsum(pred_error_list)/101, cumsum(pred_error_list_q5)/101, 
                              cumsum(naive_error_list)/101))) 
colnames(cum_MSE_q5) <- c("Index", "Data", "Value")
cum_MSE_q5[,1] <- as.numeric(cum_MSE_q5[,1]) # Indexes were being read as chars
cum_MSE_q5[,3] <- as.numeric(cum_MSE_q5[,3]) # Values were being read as chars
cum_MSE_q5
# Plot
cum_MSE_plot_q5 <- ggplot(cum_MSE_q5, aes(x=Index, y=Value, group=Data, color=Data)) +
  geom_hline(yintercept=0) + geom_line(cex = 0.8) + 
  scale_colour_manual(values = c("#fc2e17", "#13a600","#0285c2")) +
  xlab("Month") + ylab("Cumulative MSE") + ggtitle("Cumulative Mean Square Errors")
cum_MSE_plot_q5
