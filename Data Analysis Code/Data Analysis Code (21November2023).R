# November 21, 2023

library(timetk)
library(tidyverse)
library(tvReg)
library(dplyr)
library(car)
library(lubridate)
library(broom)
library(knitr)
library(gt)
library(ggplot2)
library(segmented)
library(strucchange)
library(Ake)
library(nprobust)
library(ggpattern)
library(ggpubr)

#### Functions ####

## Functions that use cumulative counts

localpoly <- function(df, alpha = 0.1, method = "ll", min_lead = 5, max_lead = 21, CI = FALSE) {
  
  library(tidyverse)
  library(tvReg)
  library(dplyr)
  library(ggplot2)
  library(timetk)
  
  names(df) <- c("Date", "y", "cumulative_y", "x", "cumulative_x")
  # reset the cumulative y/cases using the start_date
  df$cumulative_y_reset <- df$cumulative_y - df$cumulative_y[1] 
  df$cumulative_x_reset <- df$cumulative_x - df$cumulative_x[1]
  
  # Return the optimal choice of lead for Tvreg
  tv <- function(df){
    
    #Create a variable to store Mean Squared Prediction Error for each possible lead value
    MSPE_List <- c()
    
    #Loop through all possible leads and record the MSPE
    for (n_ahead_tv in min_lead:max_lead) {
      
      df_tv <- df %>%
        tk_augment_lags(c(cumulative_y_reset, y), .lags = -n_ahead_tv, .names = c("y_best_lead_tv", "daily_y_lead"))
      
      names(df_tv) <- names(df_tv) %>% str_replace_all("lag-|lag", "lead")
      
      D1 <- max(df_tv$Date) -  2 * max_lead
      
      df_tv_insample <- df_tv %>% filter (Date <= D1)
      
      tvLM_insample <- tvLM(y_best_lead_tv ~ cumulative_x_reset, 
                            est = method, 
                            data = df_tv_insample,
                            #tkernel = "Gaussian",
                            singular.ok = FALSE)
      
      newdf1 <- df_tv %>% filter (Date > D1 & Date <= D1 + n_ahead_tv) %>% 
        dplyr::select (c("cumulative_x_reset", "y_best_lead_tv", "daily_y_lead"))
      
      #get predictions for cases outside the insample
      tvLM_df_pred1 <- forecast(tvLM_insample, newdata = as.matrix(newdf1[,1]), n.ahead = n_ahead_tv)
      
      #residual <- c(tvLM_insample$residuals, tvLM_df_pred1 - newdf1$y_best_lead_tv)
      #residual <- tvLM_df_pred1 - newdf1$y_best_lead_tv
      
      daily_pred <- tail(tvLM_df_pred1, -1) - head(tvLM_df_pred1, -1)
      daily_residual <- daily_pred - newdf1$daily_y_lead[-1]
      MSPE <- mean((daily_residual)^2)
      
      MSPE_List <- c(MSPE_List, MSPE)
    }
    
    #Return the lead with the smallest MSPE
    ans <- which.min(MSPE_List) + (min_lead - 1)
    print (paste("best lead for tv:", ans))
    return(ans)
  }
  
  n_ahead_tv <- tv(df)
  
  #predictions based on TvReg
  df_tv <- df %>%
    tk_augment_lags(cumulative_y_reset, .lags = -n_ahead_tv, .names = "y_best_lead_tv")
  names(df_tv) <- names(df_tv) %>% str_replace_all("lag-|lag", "lead")
  
  D1 <- max(df_tv$Date) - 2*n_ahead_tv 
  D2 <- max(df_tv$Date) - n_ahead_tv + 1
  
  df_tv_insample <- df_tv %>% filter (Date <= D1 )
  
  tvLM_insample <- tvLM(y_best_lead_tv ~ cumulative_x_reset, 
                        est = method, 
                        data = df_tv_insample,
                        #tkernel = "Epa",
                        singular.ok = FALSE)
  
  newdf1 <- df_tv %>% filter (Date > D1 & Date < D2) %>% dplyr::select (c("cumulative_x_reset", "y_best_lead_tv"))
  
  tvLM_df_pred1 <- forecast (tvLM_insample, newdata = as.matrix(newdf1[,1]), n.ahead = n_ahead_tv)
  
  
  pred_tv <- c(tvLM_insample$fitted, tvLM_df_pred1)
  
  n <- dim(df_tv)[1] # same as dim(df_tv)[1]
  opt_tv <- data.frame(pred_tv, head(df_tv$y_best_lead_tv, n - n_ahead_tv))
  names(opt_tv) <- c("predicted_tv", "observed_tv")
  n_tv <- length(pred_tv)
  opt_tv$pred_daily_tv <- c(NA, pred_tv[2:n_tv] - pred_tv[1:(n_tv - 1)])
  
  opt_tv$date <- head(df_tv$Date, n - n_ahead_tv) + n_ahead_tv ### date being moved forward
  
  opt_tv$observed_daily <- df_tv$y[df_tv$Date %in% opt_tv$date]
  
  mspe <- mean((tail(opt_tv$observed_daily, n_ahead_tv) - tail(opt_tv$pred_daily_tv, n_ahead_tv))^2, na.rm = TRUE)  
  
  if (CI == TRUE){
    # Pointwise confidence bands
    point_conf <- function (data, n_ahead, B = 200, alpha = 0.05) { # df after the best lead being selected
      #B = 1000 ## number of bootstrap draws
      data_insample <- head (data, dim(data)[1] - n_ahead)
      nobs <- dim(data_insample)[1]
      object <- tvLM(y_best_lead_tv ~ cumulative_x_reset,
                     data = data_insample,
                     est = method)

      residuals_raw <- object$residuals - mean(object$residuals)
      residuals_b <- matrix(sample (residuals_raw * rnorm(B * nobs), size = B * nobs, replace = TRUE),
                            nrow = B, ncol = nobs)
      y_b <- matrix(object$fitted, nrow = B, ncol = nobs, byrow = TRUE) + residuals_b ## synthetic y

      ## out-of-sample prediction 1 with death counts observable
      prediction1 <- matrix(NA, nrow = B, ncol = n_ahead)
      #prediction_daily <- matrix(NA, nrow = B, ncol = n_ahead - 1)
      totobs <- nobs + n_ahead

      newdf <- data %>% tail(n_ahead) %>% dplyr::select (c("cumulative_x_reset", "y_best_lead_tv"))
      newdata <- cbind(rep(1, n_ahead), as.matrix(newdf[, 1]))
      pred_raw <- forecast (object, newdata = newdata, n.ahead = n_ahead)

      for (k in 1: B) {
        tmp <- tvLM(y_b[k, ] ~ cumulative_x_reset,
                    data = data_insample,# z = (1:nobs)/nobs,
                    est = "ll")
        #prediction1[k, ] <- predict(tmp, newdata = newdata, newz = seq(nobs + 1, nobs + n_ahead)/(nobs + n_ahead))
        prediction1[k, ] <- forecast(tmp, newdata = newdata, n.ahead = n_ahead)
        #prediction_daily[k, ] <- tail(prediction1[k, ], -1) -  head(prediction1[k, ], -1)
      }# end of k loop


      sd.est <- apply(prediction1, 2, sd)
      Q <- (prediction1 - matrix(pred_raw, nrow=B, ncol=length(pred_raw), byrow=TRUE))/
        matrix(sd.est, nrow=B, ncol=length(sd.est), byrow=TRUE)
      calpha <- apply(Q, 2, function(x){quantile(x, 1-alpha/2, na.rm = TRUE)})

      # output
      ans <- list (lower = pred_raw - sd.est * calpha, upper = pred_raw + sd.est * calpha)
    }

    # Simultaneous confidence bands
    joint_conf <- function (data, n_ahead, B = 200, alpha = 0.05) { # df after the best lead being selected
      #B = 1000 ## number of bootstrap draws
      data_insample <- head (data, dim(data)[1] - n_ahead)
      nobs <- dim(data_insample)[1]
      object <- tvLM(y_best_lead_tv ~ cumulative_x_reset,
                     data = data_insample,
                     est = method)

      residuals_raw <- object$residuals - mean(object$residuals)
      residuals_b <- matrix(sample (residuals_raw * rnorm(2 * B * nobs), size = 2 * B * nobs, replace = TRUE),
                            nrow = 2 * B, ncol = nobs)
      y_b <- matrix(object$fitted, nrow = 2 * B, ncol = nobs, byrow = TRUE) + residuals_b ## synthetic y

      ## out-of-sample prediction 1 with death counts observable
      prediction1 <- matrix(NA, nrow = 2 * B, ncol = n_ahead)
      totobs <- nobs + n_ahead

      newdf <- data %>% tail(n_ahead) %>% dplyr::select (c("cumulative_x_reset", "y_best_lead_tv"))
      newdata <- cbind(rep(1, n_ahead), as.matrix(newdf[, 1]))
      pred_raw <- forecast (object, newdata = newdata, n.ahead = n_ahead)

      for (k in 1: (2 * B)) {
        tmp <- tvLM(y_b[k, ] ~ cumulative_x_reset,
                    data = data_insample, #z = (1:nobs)/nobs,
                    est = "ll")
        #prediction1[k, ] <- predict(tmp, newdata = newdata, newz = seq(nobs + 1, nobs + n_ahead)/(nobs + n_ahead))
        prediction1[k, ] <- forecast(tmp, newdata = newdata, n.ahead = n_ahead)
      }# end of k loop

      sd.est <- apply(prediction1[1:B,], 2, sd)

      # estimate common c_alpha
      prediction2 <- prediction1[-(1:B), ]
      Q <- abs(prediction2 -
                 matrix(pred_raw, nrow=B, ncol=length(pred_raw), byrow=TRUE))/
        matrix(sd.est, nrow=B, ncol=length(sd.est), byrow=TRUE)
      Qstar <- apply(Q, 2, max)
      calpha <- quantile(Qstar, 1 - alpha/2, na.rm = TRUE)
      # output
      ans <- list (lower = pred_raw - sd.est * calpha, upper = pred_raw + sd.est * calpha)
    }

    #Add cumulative prediction bands for TVReg
    tv_conf_point <- point_conf(df_tv %>% filter (Date < D2), n_ahead_tv)
    opt_tv$plower <- c(rep(NA, length(tvLM_insample$fitted)), tv_conf_point$lower)
    opt_tv$pupper <-  c(rep(NA, length(tvLM_insample$fitted)), tv_conf_point$upper)

    tv_conf_joint <- joint_conf(df_tv %>% filter (Date < D2), n_ahead_tv)
    opt_tv$jlower <- c(rep(NA, length(tvLM_insample$fitted)), tv_conf_joint$lower)
    opt_tv$jupper <-  c(rep(NA, length(tvLM_insample$fitted)), tv_conf_joint$upper)

    # Plots
    names(opt_tv) <- c(
      "Predicted cumulative counts",
      "Reported cumulative counts",
      "Predicted daily counts",
      "date",
      "Reported daily counts",
      "Lower Confidence Bound (Pointwise)",
      "Upper Confidence Bound (Pointwise)",
      "Lower Confidence Bound (Joint)",
      "Upper Confidence Bound (Joint)"
    )
    
    region <- "XXX"

    gg_cumulative <- opt_tv %>%
      pivot_longer(cols = where(is.numeric)) %>%
      filter(name %in% c("Lower Confidence Bound (Joint)",
                         "Lower Confidence Bound (Pointwise)",
                         "Predicted cumulative counts",
                         "Reported cumulative counts",
                         "Upper Confidence Bound (Joint)",
                         "Upper Confidence Bound (Pointwise)"
      )) %>%
      ggplot(aes(date, value, color = name, linetype = name)) +
      #ylim(range(c(pred, obs), na.rm = TRUE)) +
      geom_line(size = 0.7) +
      geom_vline(xintercept = D2, color = "orange") +
      scale_linetype_manual(values=c("dotted", "dashed", "dotdash", "solid", "dotted", "dashed")) +
      scale_color_manual(values=c('purple', "blue", "red", "black", "purple", "blue")) +
      labs(
        title = paste0("Reported vs. Predicted Cumulative Counts ", region),
        subtitle = "Inputs = Cumulative Cases",
        x = "Date",
        y = "Count",
        caption = paste(method, "method lag=", n_ahead_tv, sep = " ")
      )

  } else {

    names(opt_tv) <- c(
      "Predicted cumulative counts",
      "Reported cumulative counts",
      "Predicted daily counts",
      "date",
      "Reported daily counts"
    )

    region <- "XXX"
    
    gg_cumulative <- opt_tv %>%
      pivot_longer(cols = where(is.numeric)) %>%
      filter(name %in% c("Predicted cumulative counts",
                         "Reported cumulative counts")) %>%
      ggplot(aes(date, value, color = name, linetype = name)) +
      #ylim(range(c(pred, obs), na.rm = TRUE)) +
      geom_line(size = 0.7) +
      geom_vline(xintercept = D2, color = "orange") +
      scale_linetype_manual(values=c("dotted", "solid")) +
      scale_color_manual(values=c("red", "black")) +
      labs(
        title = paste0("Reported vs. Predicted Cumulative Counts ", region),
        subtitle = "Inputs = Cumulative Cases",
        x = "Date",
        y = "Count",
        caption = paste(method, "method lag=", n_ahead_tv, sep = " ")
      )

  } # else for CI

  gg_daily <- opt_tv %>%
    pivot_longer(cols = where(is.numeric)) %>%
    filter(name %in% c("Predicted daily counts",
                       "Reported daily counts")) %>%
    ggplot(aes(date, value, color = name, linetype = name, size = name)) +
    geom_point() + geom_line() +
    geom_vline(xintercept = D2, color = "orange") +
    scale_linetype_manual(values=c("dashed", "solid")) +
    scale_color_manual(values=c('Red', 'Black')) +
    scale_size_manual(values = c(1, 0.7)) +
    labs(
      title = paste("Reported vs. Predicted Daily Counts,", region, sep = " "),
      subtitle = "Inputs = Cumulative Cases",
      x = "Date",
      y = "Count",
      caption = paste(method, "method lag =", n_ahead_tv, "days", sep = " ")
    )
  
  return(list(daily = gg_daily, cumulative = gg_cumulative, mspe = mspe, out = opt_tv, lag = n_ahead_tv, D2 = D2, bw = tvLM_insample$bw))
  
}

localpoly_bandwidthgridsearch <- function(df, alpha = 0.1, method = "ll", min_lead = 5, max_lead = 21, n_bandwidths = 10) {
  
  library(tidyverse)
  library(tvReg)
  library(dplyr)
  library(ggplot2)
  library(timetk)
  
  names(df) <- c("Date", "y", "cumulative_y", "x", "cumulative_x")
  # reset the cumulative y/cases using the start_date
  df$cumulative_y_reset <- df$cumulative_y - df$cumulative_y[1] 
  df$cumulative_x_reset <- df$cumulative_x - df$cumulative_x[1]
  
  # Return the optimal choice of lead for Tvreg
  tv <- function(df){
    
    min_bandwidth <- 5/(nrow(df) - (2 * max_lead))
    bandwidth_testing_list <- c(seq(from = min_bandwidth, to = 1, length.out = n_bandwidths), 20)
    
    MSPE_Matrix <- matrix(data = rep(NA,length(min_lead:max_lead)*length(bandwidth_testing_list)), nrow = length(min_lead:max_lead), ncol = length(bandwidth_testing_list))
    
    for (n_ahead_index in 1:length(min_lead:max_lead)) {
      
      n_ahead_tv <- (min_lead:max_lead)[n_ahead_index]
      
      for (h_index in 1:length(bandwidth_testing_list)) {
        
        h_candidate <- bandwidth_testing_list[h_index]
        
        df_tv <- df %>%
          tk_augment_lags(c(cumulative_y_reset, y), .lags = -n_ahead_tv, .names = c("y_best_lead_tv", "daily_y_lead"))
        
        names(df_tv) <- names(df_tv) %>% str_replace_all("lag-|lag", "lead")
        
        D1 <- max(df_tv$Date) -  2 * max_lead
        
        df_tv_insample <- df_tv %>% filter (Date <= D1)
        
        singular_fit_test <- try(
          tvLM(y_best_lead_tv ~ cumulative_x_reset, 
               est = method, 
               data = df_tv_insample,
               #tkernel = "Gaussian",
               bw = h_candidate,
               singular.ok = FALSE),
          silent = TRUE
        )
        
        
        if (class(singular_fit_test) == "tvlm") {
          
          # No errors when fitting model
          
          newdf1 <- df_tv %>% filter (Date > D1 & Date <= D1 + n_ahead_tv) %>% 
            dplyr::select (c("cumulative_x_reset", "y_best_lead_tv", "daily_y_lead"))
          
          #get predictions for cases outside the insample
          tvLM_df_pred1 <- forecast(tvLM_insample, newdata = as.matrix(newdf1[,1]), n.ahead = n_ahead_tv)
          
          #residual <- c(tvLM_insample$residuals, tvLM_df_pred1 - newdf1$y_best_lead_tv)
          #residual <- tvLM_df_pred1 - newdf1$y_best_lead_tv
          
          daily_pred <- tail(tvLM_df_pred1, -1) - head(tvLM_df_pred1, -1)
          daily_residual <- daily_pred - newdf1$daily_y_lead[-1]
          MSPE <- mean((daily_residual)^2)
          
          MSPE_Matrix[n_ahead_index,h_index] <- MSPE
          
        } else if ( (class(singular_fit_test) == "try-error") & (grepl("singular fit encountered", singular_fit_test))) {
          
          # Singular fit error when fitting model
          MSPE_Matrix[n_ahead_index,h_index] <- NA
          
        } else {
          
          # Another type of error when fitting model
          MSPE_Matrix[n_ahead_index,h_index] <- NA
          
        }
        
      }
      
    }
    
    #which(MSPE_Matrix == min(MSPE_Matrix), arr.ind = TRUE)
    
    best_lead <- (min_lead:max_lead)[which(MSPE_Matrix == min(MSPE_Matrix, na.rm = TRUE), arr.ind = TRUE)[1]]
    best_bandwidth <- bandwidth_testing_list[which(MSPE_Matrix == min(MSPE_Matrix, na.rm = TRUE), arr.ind = TRUE)[2]]
    
    return(c(best_lead, best_bandwidth))
    
  }
  
  gridsearch_results <- tv(df)
  
  n_ahead_tv <- gridsearch_results[1]
  best_bandwidth <- gridsearch_results[2]
  
  #predictions based on TvReg
  df_tv <- df %>%
    tk_augment_lags(cumulative_y_reset, .lags = -n_ahead_tv, .names = "y_best_lead_tv")
  names(df_tv) <- names(df_tv) %>% str_replace_all("lag-|lag", "lead")
  
  D1 <- max(df_tv$Date) - 2*n_ahead_tv 
  D2 <- max(df_tv$Date) - n_ahead_tv + 1
  
  df_tv_insample <- df_tv %>% filter (Date <= D1 )
  
  tvLM_insample <- tvLM(y_best_lead_tv ~ cumulative_x_reset, 
                        est = method, 
                        data = df_tv_insample,
                        #tkernel = "Epa",
                        bw = best_bandwidth,
                        singular.ok = FALSE)
  
  newdf1 <- df_tv %>% filter (Date > D1 & Date < D2) %>% dplyr::select (c("cumulative_x_reset", "y_best_lead_tv"))
  
  tvLM_df_pred1 <- forecast (tvLM_insample, newdata = as.matrix(newdf1[,1]), n.ahead = n_ahead_tv)
  
  
  pred_tv <- c(tvLM_insample$fitted, tvLM_df_pred1)
  
  n <- dim(df_tv)[1] # same as dim(df_tv)[1]
  opt_tv <- data.frame(pred_tv, head(df_tv$y_best_lead_tv, n - n_ahead_tv))
  names(opt_tv) <- c("predicted_tv", "observed_tv")
  n_tv <- length(pred_tv)
  opt_tv$pred_daily_tv <- c(NA, pred_tv[2:n_tv] - pred_tv[1:(n_tv - 1)])
  
  opt_tv$date <- head(df_tv$Date, n - n_ahead_tv) + n_ahead_tv ### date being moved forward
  
  opt_tv$observed_daily <- df_tv$y[df_tv$Date %in% opt_tv$date]
  
  mspe <- mean((tail(opt_tv$observed_daily, n_ahead_tv) - tail(opt_tv$pred_daily_tv, n_ahead_tv))^2, na.rm = TRUE)  
  
  return(list(mspe = mspe, out = opt_tv, lag = n_ahead_tv, D2 = D2, bw = best_bandwidth))
  
}

lp_predict <- function(y_observed, x_observed, newdata, n.ahead = 1, p = 1, kernel = "tri", tau = "bc") {
  
  # Check arguments are valid
  if (n.ahead == 1) {
    newdata <- matrix(newdata, ncol = length(newdata))
  }
  
  if (NROW(newdata) != n.ahead) {
    stop("\nDimensions of 'newdata' are not compatible with 'n.ahead'.\n")
  }
  
  if (length(y_observed) != length(x_observed)) {
    stop("\nDimensions of x and y are identical.\n")
  }
  
  if (tau == "bc") {
    col = 6
  } else if (tau == "us") {
    col = 5
  } else if (tau != "bc" & tau != "us") {
    stop("\ntau must be either 'bc' or 'us'.\n")
  }
  
  
  
  
  # Fit a local polynomial model to the insample data. This is only for reference (this gets outputted at the end of the function) - we don't actually use this for the predictions
  ll_results_original <- lprobust(y = y_observed,
                                  x = x_observed,
                                  eval = x_observed,
                                  p = p, #Degree of polynomial to be fitted
                                  kernel = kernel
  )
  
  # Define some variables that will be used in the prediction loop  
  obs <- length(y_observed)
  prediction <- numeric(n.ahead)
  totobs <- obs + n.ahead
  X <- x_observed
  Y <- y_observed
  
  X_uncorrected <- x_observed
  Y_uncorrected <- y_observed
  
  # Make predictions. Predictions are calculated one-at-a-time, and after each prediction the next new X value and the corresponding prediction are added on to the X and Y variables so that they can be used in the next prediction
  for (t in 1:n.ahead) {
    i <- 1
    prediction_x <- X[i:(obs + t - 1)]
    prediction_y <- Y[i:(obs + t - 1)]
    
    ll_results_new <- lprobust(y = prediction_y,
                               x = prediction_x,
                               eval = newdata[t],
                               p = p, #Degree of polynomial to be fitted
                               kernel = kernel
    )
    
    # The prediction is the tau.bc value returned in the "Estimate" data frame
    prediction[t] <- ll_results_new$Estimate[,col]
    
    X <- c(X, newdata[t])
    Y <- c(Y, prediction[t])
  }
  
  return(list(OriginalModel = ll_results_original$Estimate, NewData = newdata, Predictions = tail(Y, n.ahead)))
}

lprobust_search <- function(df, method = "ll", min_lead = 5, max_lead = 21, tau = "bc") {
  
  names(df) <- c("Date", "y", "cumulative_y", "x", "cumulative_x")
  # reset the cumulative y/cases using the start_date
  df$cumulative_y_reset <- df$cumulative_y - df$cumulative_y[1] 
  df$cumulative_x_reset <- df$cumulative_x - df$cumulative_x[1]
  
  if (method == "lc") {
    p = 0
  } else if (method == "ll") {
    p = 1
  } else {
    p = 1
  }
  
  # Return the optimal choice of lead for Tvreg
  lpselection <- function(df, p){
    
    #Create a variable to store Mean Squared Prediction Error for each possible lead value
    MSPE_List <- c()
    
    #Loop through all possible leads and record the MSPE
    for (n_ahead_lp in min_lead:max_lead) {
      
      df_tv <- df %>%
        tk_augment_lags(c(cumulative_y_reset, y), .lags = -n_ahead_lp, .names = c("y_best_lead_tv", "daily_y_lead"))
      
      names(df_tv) <- names(df_tv) %>% str_replace_all("lag-|lag", "lead")
      
      D1 <- max(df_tv$Date) -  2 * max_lead
      
      df_tv_insample <- df_tv %>% filter (Date <= D1)
      
      newdf1 <- df_tv %>% filter (Date > D1 & Date <= D1 + n_ahead_lp) %>% 
        dplyr::select (c("cumulative_x_reset", "y_best_lead_tv", "daily_y_lead"))
      
      
      lp_forecast_results <- lp_predict(y_observed = df_tv_insample$y_best_lead_tv,
                                        x_observed = df_tv_insample$cumulative_x_reset,
                                        newdata = newdf1$cumulative_x_reset,
                                        p = p,
                                        n.ahead = length(newdf1$cumulative_x_reset),
                                        tau = tau)
      
      daily_pred <- tail(lp_forecast_results$Predictions, -1) - head(lp_forecast_results$Predictions, -1)
      daily_residual <- daily_pred - newdf1$daily_y_lead[-1]
      MSPE <- mean((daily_residual)^2)
      
      MSPE_List <- c(MSPE_List, MSPE)
    }
    
    #Return the lead with the smallest MSPE
    ans <- which.min(MSPE_List) + (min_lead - 1)
    print (paste("best lead for lp:", ans))
    return(ans)
  }
  
  n_ahead_lp <- lpselection(df, p)
  
  #predictions based on TvReg
  df_tv <- df %>%
    tk_augment_lags(cumulative_y_reset, .lags = -n_ahead_lp, .names = "y_best_lead_tv")
  names(df_tv) <- names(df_tv) %>% str_replace_all("lag-|lag", "lead")
  
  D1 <- max(df_tv$Date) - 2*n_ahead_lp 
  D2 <- max(df_tv$Date) - n_ahead_lp + 1
  
  df_tv_insample <- df_tv %>% filter (Date <= D1)
  
  newdf1 <- df_tv %>% filter (Date > D1 & Date < D2) %>% dplyr::select (c("cumulative_x_reset", "y_best_lead_tv"))
  
  lp_forecast_results <- lp_predict(y_observed = df_tv_insample$y_best_lead_tv,
                                    x_observed = df_tv_insample$cumulative_x_reset,
                                    newdata = newdf1$cumulative_x_reset,
                                    p = p,
                                    n.ahead = length(newdf1$cumulative_x_reset),
                                    tau = tau)
  
  
  
  
  pred_lp <- c(lp_forecast_results$OriginalModel[,6], lp_forecast_results$Predictions)
  
  n <- dim(df_tv)[1] # same as dim(df_tv)[1]
  opt_tv <- data.frame(pred_lp, head(df_tv$y_best_lead_tv, n - n_ahead_lp))
  names(opt_tv) <- c("predicted_tv", "observed_tv")
  n_tv <- length(pred_lp)
  opt_tv$pred_daily_tv <- c(NA, pred_lp[2:n_tv] - pred_lp[1:(n_tv - 1)])
  
  opt_tv$date <- head(df_tv$Date, n - n_ahead_lp) + n_ahead_lp ### date being moved forward
  
  opt_tv$observed_daily <- df_tv$y[df_tv$Date %in% opt_tv$date]
  
  mspe <- mean((tail(opt_tv$observed_daily, n_ahead_lp) - tail(opt_tv$pred_daily_tv, n_ahead_lp))^2, na.rm = TRUE)  
  
  return(list(mspe = mspe, out = opt_tv, lag = n_ahead_lp, D2 = D2))
  
}

## Functions that use cumulative x counts to predict raw y counts

localpoly_mixed <- function(df, alpha = 0.1, method = "ll", min_lead = 5, max_lead = 21) {
  
  library(tidyverse)
  library(tvReg)
  library(dplyr)
  library(ggplot2)
  library(timetk)
  
  names(df) <- c("Date", "y", "cumulative_y", "x", "cumulative_x")
  # reset the cumulative y/cases using the start_date
  df$cumulative_y_reset <- df$cumulative_y - df$cumulative_y[1] 
  df$cumulative_x_reset <- df$cumulative_x - df$cumulative_x[1]
  
  # Return the optimal choice of lead for Tvreg
  tv <- function(df){
    
    #Create a variable to store Mean Squared Prediction Error for each possible lead value
    MSPE_List <- c()
    
    #Loop through all possible leads and record the MSPE
    for (n_ahead_tv in min_lead:max_lead) {
      
      df_tv <- df %>%
        tk_augment_lags(c(cumulative_y_reset, y), .lags = -n_ahead_tv, .names = c("y_best_lead_tv", "daily_y_lead"))
      
      names(df_tv) <- names(df_tv) %>% str_replace_all("lag-|lag", "lead")
      
      D1 <- max(df_tv$Date) -  2 * max_lead
      
      df_tv_insample <- df_tv %>% filter (Date <= D1)
      
      tvLM_insample <- tvLM(daily_y_lead ~ cumulative_x_reset, 
                            est = method, 
                            data = df_tv_insample,
                            #tkernel = "Gaussian",
                            singular.ok = FALSE)
      
      newdf1 <- df_tv %>% filter (Date > D1 & Date <= D1 + n_ahead_tv) %>% 
        dplyr::select (c("cumulative_x_reset", "y_best_lead_tv", "daily_y_lead"))
      
      #get predictions for cases outside the insample
      tvLM_df_pred1 <- forecast(tvLM_insample, newdata = as.matrix(newdf1[,1]), n.ahead = n_ahead_tv)
      
      daily_residual <- tvLM_df_pred1 - newdf1$daily_y_lead
      MSPE <- mean((daily_residual)^2)
      
      MSPE_List <- c(MSPE_List, MSPE)
    }
    
    #Return the lead with the smallest MSPE
    ans <- which.min(MSPE_List) + (min_lead - 1)
    print (paste("best lead for tv:", ans))
    return(ans)
  }
  
  n_ahead_tv <- tv(df)
  
  #predictions based on TvReg
  df_tv <- df %>%
    tk_augment_lags(c(cumulative_y_reset, y), .lags = -n_ahead_tv, .names = c("y_best_lead_tv", "daily_y_lead"))
  names(df_tv) <- names(df_tv) %>% str_replace_all("lag-|lag", "lead")
  
  D1 <- max(df_tv$Date) - 2*n_ahead_tv 
  D2 <- max(df_tv$Date) - n_ahead_tv + 1
  
  df_tv_insample <- df_tv %>% filter (Date <= D1 )
  
  tvLM_insample <- tvLM(daily_y_lead ~ cumulative_x_reset, 
                        est = method, 
                        data = df_tv_insample,
                        #tkernel = "Epa",
                        singular.ok = FALSE)
  
  newdf1 <- df_tv %>% filter(Date > D1 & Date < D2) %>% dplyr::select (c("cumulative_x_reset", "y_best_lead_tv", "daily_y_lead"))
  
  tvLM_df_pred1 <- forecast(tvLM_insample, newdata = as.matrix(newdf1[,1]), n.ahead = n_ahead_tv)
  
  # pred_tv <- c(tvLM_insample$fitted, tvLM_df_pred1)
  # n <- dim(df_tv)[1] # same as dim(df_tv)[1]
  # opt_tv <- data.frame(pred_tv, head(df_tv$y_best_lead_tv, n - n_ahead_tv))
  # names(opt_tv) <- c("predicted_tv", "observed_tv")
  # n_tv <- length(pred_tv)
  # opt_tv$pred_daily_tv <- c(NA, pred_tv[2:n_tv] - pred_tv[1:(n_tv - 1)])
  # opt_tv$date <- head(df_tv$Date, n - n_ahead_tv) + n_ahead_tv ### date being moved forward
  # opt_tv$observed_daily <- df_tv$y[df_tv$Date %in% opt_tv$date]
  # mspe <- mean((tail(opt_tv$observed_daily, n_ahead_tv) - tail(opt_tv$pred_daily_tv, n_ahead_tv))^2, na.rm = TRUE)  
  
  daily_residual <- tvLM_df_pred1 - newdf1$daily_y_lead
  mspe <- mean((daily_residual)^2)
  
  #out = opt_tv
  return(list(mspe = mspe, lag = n_ahead_tv, D2 = D2, bw = tvLM_insample$bw))
  
}

lprobust_search_mixed <- function(df, method = "ll", min_lead = 5, max_lead = 21, tau = "bc") {
  
  names(df) <- c("Date", "y", "cumulative_y", "x", "cumulative_x")
  # reset the cumulative y/cases using the start_date
  df$cumulative_y_reset <- df$cumulative_y - df$cumulative_y[1] 
  df$cumulative_x_reset <- df$cumulative_x - df$cumulative_x[1]
  
  if (method == "lc") {
    p = 0
  } else if (method == "ll") {
    p = 1
  } else {
    p = 1
  }
  
  # Return the optimal choice of lead for Tvreg
  lpselection <- function(df, p){
    
    #Create a variable to store Mean Squared Prediction Error for each possible lead value
    MSPE_List <- c()
    
    #Loop through all possible leads and record the MSPE
    for (n_ahead_lp in min_lead:max_lead) {
      
      df_tv <- df %>%
        tk_augment_lags(c(cumulative_y_reset, y), .lags = -n_ahead_lp, .names = c("y_best_lead_tv", "daily_y_lead"))
      
      names(df_tv) <- names(df_tv) %>% str_replace_all("lag-|lag", "lead")
      
      D1 <- max(df_tv$Date) -  2 * max_lead
      
      df_tv_insample <- df_tv %>% filter (Date <= D1)
      
      newdf1 <- df_tv %>% filter (Date > D1 & Date <= D1 + n_ahead_lp) %>% 
        dplyr::select (c("cumulative_x_reset", "y_best_lead_tv", "daily_y_lead"))
      
      
      lp_forecast_results <- lp_predict(y_observed = df_tv_insample$daily_y_lead,
                                        x_observed = df_tv_insample$cumulative_x_reset,
                                        newdata = newdf1$cumulative_x_reset,
                                        p = p,
                                        n.ahead = length(newdf1$cumulative_x_reset),
                                        tau = tau)
      
      
      
      
      daily_residual <- lp_forecast_results$Predictions - newdf1$daily_y_lead
      MSPE <- mean((daily_residual)^2)
      
      MSPE_List <- c(MSPE_List, MSPE)
    }
    
    #Return the lead with the smallest MSPE
    ans <- which.min(MSPE_List) + (min_lead - 1)
    print (paste("best lead for lp:", ans))
    return(ans)
  }
  
  n_ahead_lp <- lpselection(df, p)
  
  #predictions based on TvReg
  df_tv <- df %>%
    tk_augment_lags(c(cumulative_y_reset, y), .lags = -n_ahead_lp, .names = c("y_best_lead_tv", "daily_y_lead"))
  names(df_tv) <- names(df_tv) %>% str_replace_all("lag-|lag", "lead")
  
  D1 <- max(df_tv$Date) - 2*n_ahead_lp 
  D2 <- max(df_tv$Date) - n_ahead_lp + 1
  
  df_tv_insample <- df_tv %>% filter (Date <= D1)
  
  newdf1 <- df_tv %>% filter (Date > D1 & Date < D2) %>% dplyr::select (c("cumulative_x_reset", "y_best_lead_tv", "daily_y_lead"))
  
  lp_forecast_results <- lp_predict(y_observed = df_tv_insample$daily_y_lead,
                                    x_observed = df_tv_insample$cumulative_x_reset,
                                    newdata = newdf1$cumulative_x_reset,
                                    p = p,
                                    n.ahead = length(newdf1$cumulative_x_reset),
                                    tau = tau)
  
  
  
  
  
  daily_residual <- lp_forecast_results$Predictions - newdf1$daily_y_lead
  mspe <- mean((daily_residual)^2)
  
  # pred_lp <- c(lp_forecast_results$OriginalModel[,6], lp_forecast_results$Predictions)
  # n <- dim(df_tv)[1] # same as dim(df_tv)[1]
  # opt_tv <- data.frame(pred_lp, head(df_tv$y_best_lead_tv, n - n_ahead_lp))
  # names(opt_tv) <- c("predicted_tv", "observed_tv")
  # n_tv <- length(pred_lp)
  # opt_tv$pred_daily_tv <- c(NA, pred_lp[2:n_tv] - pred_lp[1:(n_tv - 1)])
  # opt_tv$date <- head(df_tv$Date, n - n_ahead_lp) + n_ahead_lp ### date being moved forward
  # opt_tv$observed_daily <- df_tv$y[df_tv$Date %in% opt_tv$date]
  # mspe <- mean((tail(opt_tv$observed_daily, n_ahead_lp) - tail(opt_tv$pred_daily_tv, n_ahead_lp))^2, na.rm = TRUE)  
  # 
  return(list(mspe = mspe, lag = n_ahead_lp, D2 = D2))
  
}

#### Country Results - Setup ####

#Load in the country-level data and select the variables we are interested in
urlnewOWIDdata <- 'https://covid.ourworldindata.org/data/owid-covid-data.csv'
NewOWIDData <- read.csv(urlnewOWIDdata)
NewOWIDData <- NewOWIDData %>% dplyr::select(date, location, total_cases, new_cases, total_deaths, new_deaths, icu_patients, hosp_patients)

#Reformat the date values so that we can easily filter by them
NewOWIDData$Date <- as.Date(NewOWIDData$date, format="%Y-%m-%d") 


# #Exclude all data after a given end date so that results are reproducible
# cutoff_end <- "2022-12-31"
# NewOWIDData <- NewOWIDData %>% filter(Date <= cutoff_end)

#Filter to the countries of interest
LocationsList <- c("Canada", "United States", "Israel", "South Africa", "Japan", "South Korea", "Brazil", "Pakistan", "Thailand")

NewOWIDRegionalData <- NewOWIDData %>% filter(location %in% LocationsList)

#Omicron wave start dates for each country (found by visual inspection of the daily new cases plots)
Omicron_Start_Dates_Countries <- data.frame("Canada" = "2021-12-12", 
                                            "United States" = "2021-12-11", 
                                            "Israel" = "2021-12-30", 
                                            "Germany" = "2021-10-27", 
                                            "South Africa" = "2021-11-27", 
                                            "Japan" = "2022-01-10", 
                                            "South Korea" = "2022-01-30",
                                            "Brazil" = "2022-01-01",
                                            "France" = "2021-12-01",
                                            "Pakistan" = "2022-01-01",
                                            "Thailand" = "2022-01-01")

#Reporting end dates for each country (found by visual inspection of the daily new cases plots)
Reporting_End_Dates_Countries <- data.frame("Canada" = "2022-03-12", 
                                            "United States" = "2022-03-11", 
                                            "Israel" = "2022-03-30", 
                                            "Germany" = "2022-01-25", 
                                            "South Africa" = "2022-02-25", 
                                            "Japan" = "2022-04-10", 
                                            "South Korea" = "2022-04-30",
                                            "Brazil" = "2022-04-01",
                                            "France" = "2022-03-01",
                                            "Pakistan" = "2022-04-01",
                                            "Thailand" = "2022-04-01")



#### SECTION 3.2 ####
#### Country Results - Gridsearch Comparison -  Deaths (Cumulative X Predicting Cumulative Y) ####

#Set working directory
#setwd("C:/Users/gamem/Documents/R Files/Research/COVID-Modeling/Plots/Plots27July2023")

# Countries with Deaths data
LocationsList <- c("Canada", "United States", "Israel", "South Africa", "Japan", "South Korea", "Brazil", "Pakistan", "Thailand")

CountryDeathsGridMSPEs <- data.frame(Location = rep(LocationsList, each = 4),
                                     Model = rep(c("tvReg_LC","tvReg_LC_grid","tvReg_LL","tvReg_LL_grid"),length(LocationsList)),
                                     MSPE = rep(NA, length(LocationsList)*4),
                                     Lag = rep(NA, length(LocationsList)*4))

for (reg in 1:length(LocationsList)) {
  
  region <- LocationsList[reg]
  
  #Print the name of the current country to the console
  print(paste(LocationsList[reg], "Deaths"))
  
  #Filter the data to a single region
  NewOWIDFullSingleRegionData <- NewOWIDRegionalData %>% filter(location %in% LocationsList[reg])
  
  nDays <- 90
  
  #Filter the data to only use the correct number of days
  NewOWIDSingleRegionData <- NewOWIDFullSingleRegionData %>% filter(Date >= Omicron_Start_Dates_Countries[[sub(" ",".", LocationsList[reg])]])
  NewOWIDSingleRegionData <- NewOWIDSingleRegionData %>% filter(Date < as.Date(Omicron_Start_Dates_Countries[[sub(" ",".", LocationsList[reg])]]) + nDays)
  
  
  #Select the response variables
  Deaths_long <- NewOWIDSingleRegionData %>% dplyr::select(date, new_deaths)
  names(Deaths_long) <- c("Date", "Deaths")
  Deaths_long$Deaths[is.na(Deaths_long$Deaths)] <- 0
  Deaths_long$`Deaths Cumulative` <- cumsum(Deaths_long$Deaths)
  
  #Select the predictor variables
  cases_long <- NewOWIDSingleRegionData %>% dplyr::select(date, new_cases, total_cases)
  names(cases_long) <- c("Date","cases", "cumulative_cases")
  
  #Combine the response and predictor variables into a single data frame
  df <-  merge (Deaths_long[,c("Date", "Deaths", "Deaths Cumulative")], 
                cases_long[,c( "cases", "cumulative_cases", "Date" ) ], 
                by = "Date") 
  
  #Reformat the dates in the data frame
  df$Date <- as.Date(df$Date, format="%Y-%m-%d")
  
  debugonce(localpoly)
  
  lcresults <- localpoly(df, alpha = 0.1, method = "lc", min_lead = 5, max_lead = 21, CI = FALSE)
  lcresults_grid <- localpoly_bandwidthgridsearch(df, alpha = 0.1, method = "lc", min_lead = 5, max_lead = 21, n_bandwidths = 50)
  
  llresults <- localpoly(df, alpha = 0.1, method = "ll", min_lead = 5, max_lead = 21, CI = FALSE)
  llresults_grid <- localpoly_bandwidthgridsearch(df, alpha = 0.1, method = "ll", min_lead = 5, max_lead = 21, n_bandwidths = 50)
  
  CountryDeathsGridMSPEs[((reg-1)*4)+1,3] <- lcresults$mspe
  CountryDeathsGridMSPEs[((reg-1)*4)+2,3] <- lcresults_grid$mspe
  CountryDeathsGridMSPEs[((reg-1)*4)+3,3] <- llresults$mspe
  CountryDeathsGridMSPEs[((reg-1)*4)+4,3] <- llresults_grid$mspe
  
  CountryDeathsGridMSPEs[((reg-1)*4)+1,4] <- lcresults$lag
  CountryDeathsGridMSPEs[((reg-1)*4)+2,4] <- lcresults_grid$lag
  CountryDeathsGridMSPEs[((reg-1)*4)+3,4] <- llresults$lag
  CountryDeathsGridMSPEs[((reg-1)*4)+4,4] <- llresults_grid$lag
  
}

View(CountryDeathsGridMSPEs)


ggplot(data=CountryDeathsGridMSPEs, aes(fill=Model, y=log(MSPE), x=Location, pattern=Model)) +
  geom_col_pattern(position="dodge",
                   colour = "black",
                   #pattern_key_scale_factor = 0.5,
                   #pattern_density = 0.1,
                   #pattern_spacing = 0.01,
                   pattern_fill = "black"
  ) +
  scale_fill_manual(values=rep(c("#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65"),times=3)) +
  scale_pattern_manual(values = rep(c("none", "stripe", "circle"), each=6))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  ggtitle("Log of MSPE") +
  theme(plot.title = element_text(hjust = 0.5))







#### SECTION 3.3 ####
#### Country Results - Transformation Comparison - ICU Patients ####

#Set working directory
#setwd("C:/Users/gamem/Documents/R Files/Research/COVID-Modeling/Plots/Plots27July2023")

CountryICUMSPEs <- data.frame(Location = rep(LocationsList, each = 18),
                              Model = rep(c("tvReg_LC_CumulativeY","tvReg_LC_mixedY","tvReg_LC_rawY","tvReg_LL_CumulativeY","tvReg_LL_mixedY","tvReg_LL_rawY","LProbust_LC_bc_CumulativeY","LProbust_LC_bc_mixedY","LProbust_LC_bc_rawY","LProbust_LC_us_CumulativeY","LProbust_LC_us_mixedY","LProbust_LC_us_rawY","LProbust_LL_bc_CumulativeY","LProbust_LL_bc_mixedY","LProbust_LL_bc_rawY","LProbust_LL_us_CumulativeY","LProbust_LL_us_mixedY","LProbust_LL_us_rawY"),length(LocationsList)),
                              MSPE = rep(NA, length(LocationsList)*18),
                              Lag = rep(NA, length(LocationsList)*18))

#for (reg in c(1,2,3,4,6)) {
for (reg in 1:length(LocationsList)) {
  
  region <- LocationsList[reg]
  
  #Filter the data to a single region
  NewOWIDSingleRegionData <- NewOWIDRegionalData %>% filter(location %in% LocationsList[reg])
  
  #Filter the data to only use data after the start of the Omicron wave in the current region
  NewOWIDSingleRegionData <- NewOWIDSingleRegionData %>% filter(Date >= Omicron_Start_Dates_Countries[[sub(" ",".", LocationsList[reg])]])
  NewOWIDSingleRegionData <- NewOWIDSingleRegionData %>% filter(Date <= Reporting_End_Dates_Countries[[sub(" ",".", LocationsList[reg])]])
  
  
  #Select the response variables
  ICU_long <- NewOWIDSingleRegionData %>% dplyr::select(date, icu_patients)
  names(ICU_long) <- c("Date", "ICU Patients")
  ICU_long$`ICU Patients Cumulative` <- cumsum(ICU_long$`ICU Patients`)
  
  #Select the predictor variables
  cases_long <- NewOWIDSingleRegionData %>% dplyr::select(date, new_cases, total_cases)
  names(cases_long) <- c("Date","cases", "cumulative_cases")
  
  #Combine the response and predictor variables into a single data frame
  df <-  merge (ICU_long[,c("Date", "ICU Patients", "ICU Patients Cumulative")], 
                cases_long[,c( "cases", "cumulative_cases", "Date" ) ], 
                by = "Date") 
  
  #Reformat the dates in the data frame
  df$Date <- as.Date(df$Date, format="%Y-%m-%d")
  
  #Print the name of the current country to the console
  print(paste(LocationsList[reg], "ICU Patients"))
  
  lcresults <- localpoly(df, alpha = 0.1, method = "lc", min_lead = 5, max_lead = 21, CI = FALSE)
  lcresults_mixed <- localpoly_mixed(df, alpha = 0.1, method = "lc", min_lead = 5, max_lead = 21)
  lcresults_raw <- localpoly_raw(df, alpha = 0.1, method = "lc", min_lead = 5, max_lead = 21)
  
  llresults <- localpoly(df, alpha = 0.1, method = "ll", min_lead = 5, max_lead = 21, CI = FALSE)
  llresults_mixed <- localpoly_mixed(df, alpha = 0.1, method = "ll", min_lead = 5, max_lead = 21)
  llresults_raw <- localpoly_raw(df, alpha = 0.1, method = "ll", min_lead = 5, max_lead = 21)
  
  lprobust_results_lc_bc <- lprobust_search(df, method = "lc", min_lead = 5, max_lead = 21, tau = "bc")
  lprobust_results_lc_bc_mixed <- lprobust_search_mixed(df, method = "lc", min_lead = 5, max_lead = 21, tau = "bc")
  lprobust_results_lc_bc_raw <- lprobust_search_raw(df, method = "lc", min_lead = 5, max_lead = 21, tau = "bc")
  
  lprobust_results_lc_us <- lprobust_search(df, method = "lc", min_lead = 5, max_lead = 21, tau = "us")
  lprobust_results_lc_us_mixed <- lprobust_search_mixed(df, method = "lc", min_lead = 5, max_lead = 21, tau = "us")
  lprobust_results_lc_us_raw <- lprobust_search_raw(df, method = "lc", min_lead = 5, max_lead = 21, tau = "us")
  
  lprobust_results_ll_bc <- lprobust_search(df, method = "ll", min_lead = 5, max_lead = 21, tau = "bc")
  lprobust_results_ll_bc_mixed <- lprobust_search_mixed(df, method = "ll", min_lead = 5, max_lead = 21, tau = "bc")
  lprobust_results_ll_bc_raw <- lprobust_search_raw(df, method = "ll", min_lead = 5, max_lead = 21, tau = "bc")
  
  lprobust_results_ll_us <- lprobust_search(df, method = "ll", min_lead = 5, max_lead = 21, tau = "us")
  lprobust_results_ll_us_mixed <- lprobust_search_mixed(df, method = "ll", min_lead = 5, max_lead = 21, tau = "us")
  lprobust_results_ll_us_raw <- lprobust_search_raw(df, method = "ll", min_lead = 5, max_lead = 21, tau = "us")
  
  
  CountryICUMSPEs[((reg-1)*18)+1,3] <- lcresults$mspe
  CountryICUMSPEs[((reg-1)*18)+2,3] <- lcresults_mixed$mspe
  CountryICUMSPEs[((reg-1)*18)+3,3] <- lcresults_raw$mspe
  CountryICUMSPEs[((reg-1)*18)+4,3] <- llresults$mspe
  CountryICUMSPEs[((reg-1)*18)+5,3] <- llresults_mixed$mspe
  CountryICUMSPEs[((reg-1)*18)+6,3] <- llresults_raw$mspe
  CountryICUMSPEs[((reg-1)*18)+7,3] <- lprobust_results_lc_bc$mspe
  CountryICUMSPEs[((reg-1)*18)+8,3] <- lprobust_results_lc_bc_mixed$mspe
  CountryICUMSPEs[((reg-1)*18)+9,3] <- lprobust_results_lc_bc_raw$mspe
  CountryICUMSPEs[((reg-1)*18)+10,3] <- lprobust_results_lc_us$mspe
  CountryICUMSPEs[((reg-1)*18)+11,3] <- lprobust_results_lc_us_mixed$mspe
  CountryICUMSPEs[((reg-1)*18)+12,3] <- lprobust_results_lc_us_raw$mspe
  CountryICUMSPEs[((reg-1)*18)+13,3] <- lprobust_results_ll_bc$mspe
  CountryICUMSPEs[((reg-1)*18)+14,3] <- lprobust_results_ll_bc_mixed$mspe
  CountryICUMSPEs[((reg-1)*18)+15,3] <- lprobust_results_ll_bc_raw$mspe
  CountryICUMSPEs[((reg-1)*18)+16,3] <- lprobust_results_ll_us$mspe
  CountryICUMSPEs[((reg-1)*18)+17,3] <- lprobust_results_ll_us_mixed$mspe
  CountryICUMSPEs[((reg-1)*18)+18,3] <- lprobust_results_ll_us_raw$mspe
  
  CountryICUMSPEs[((reg-1)*18)+1,4] <- lcresults$lag
  CountryICUMSPEs[((reg-1)*18)+2,4] <- lcresults_mixed$lag
  CountryICUMSPEs[((reg-1)*18)+3,4] <- lcresults_raw$lag
  CountryICUMSPEs[((reg-1)*18)+4,4] <- llresults$lag
  CountryICUMSPEs[((reg-1)*18)+5,4] <- llresults_mixed$lag
  CountryICUMSPEs[((reg-1)*18)+6,4] <- llresults_raw$lag
  CountryICUMSPEs[((reg-1)*18)+7,4] <- lprobust_results_lc_bc$lag
  CountryICUMSPEs[((reg-1)*18)+8,4] <- lprobust_results_lc_bc_mixed$lag
  CountryICUMSPEs[((reg-1)*18)+9,4] <- lprobust_results_lc_bc_raw$lag
  CountryICUMSPEs[((reg-1)*18)+10,4] <- lprobust_results_lc_us$lag
  CountryICUMSPEs[((reg-1)*18)+11,4] <- lprobust_results_lc_us_mixed$lag
  CountryICUMSPEs[((reg-1)*18)+12,4] <- lprobust_results_lc_us_raw$lag
  CountryICUMSPEs[((reg-1)*18)+13,4] <- lprobust_results_ll_bc$lag
  CountryICUMSPEs[((reg-1)*18)+14,4] <- lprobust_results_ll_bc_mixed$lag
  CountryICUMSPEs[((reg-1)*18)+15,4] <- lprobust_results_ll_bc_raw$lag
  CountryICUMSPEs[((reg-1)*18)+16,4] <- lprobust_results_ll_us$lag
  CountryICUMSPEs[((reg-1)*18)+17,4] <- lprobust_results_ll_us_mixed$lag
  CountryICUMSPEs[((reg-1)*18)+18,4] <- lprobust_results_ll_us_raw$lag
}

View(CountryICUMSPEs)

# Nice colour palette: ["#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65", "#beb9db", "#fdcce5", "#8bd3c7"]

ggplot(data=CountryICUMSPEs, aes(fill=Model, y=log(MSPE), x=Location, pattern=Model)) +
  geom_col_pattern(position="dodge", 
                   colour = "black",
                   #pattern_key_scale_factor = 0.5,
                   #pattern_density = 0.1,
                   #pattern_spacing = 0.01,
                   pattern_fill = "black"
  ) + 
  scale_fill_manual(values=rep(c("#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65"),each=3)) +
  scale_pattern_manual(values = rep(c("none", "stripe", "circle"), 6))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  ggtitle("Log of MSPE") +
  theme(plot.title = element_text(hjust = 0.5))


# ggplot(data=CountryICUMSPEs, aes(fill=Model, y=log(MSPE), x=Location, pattern=Model)) +
#   geom_col_pattern(position="dodge", 
#                    colour = "black",
#                    pattern_fill = "black",
#                    pattern_density = 0.1,
#                    pattern_spacing = 0.01,
#                    pattern_key_scale_factor = 0.5) + 
#   scale_fill_manual(values=rep(c("#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65"),each=3)) +
#   scale_pattern_manual(values = rep(c("none", "stripe", "circle"), 6))+
#   scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
#   ggtitle("Log of MSPE") +
#   theme(plot.title = element_text(hjust = 0.5))


# Model3 <- lm(MSPE ~ Model, data = CountryICUMSPEs)
# anova(Model3)
# summary(Model3)



#### Country Results - Transformation Comparison - Hospitalizations ####

#Set working directory
#setwd("C:/Users/gamem/Documents/R Files/Research/COVID-Modeling/Plots/Plots27July2023")

CountryHospMSPEs <- data.frame(Location = rep(LocationsList, each = 18),
                               Model = rep(c("tvReg_LC_CumulativeY","tvReg_LC_mixedY","tvReg_LC_rawY","tvReg_LL_CumulativeY","tvReg_LL_mixedY","tvReg_LL_rawY","LProbust_LC_bc_CumulativeY","LProbust_LC_bc_mixedY","LProbust_LC_bc_rawY","LProbust_LC_us_CumulativeY","LProbust_LC_us_mixedY","LProbust_LC_us_rawY","LProbust_LL_bc_CumulativeY","LProbust_LL_bc_mixedY","LProbust_LL_bc_rawY","LProbust_LL_us_CumulativeY","LProbust_LL_us_mixedY","LProbust_LL_us_rawY"),length(LocationsList)),
                               MSPE = rep(NA, length(LocationsList)*18),
                               Lag = rep(NA, length(LocationsList)*18))


#for (reg in c(1,2,3,4,6)) {
for (reg in 1:length(LocationsList)) {
  
  region <- LocationsList[reg]
  
  #Filter the data to a single region
  NewOWIDSingleRegionData <- NewOWIDRegionalData %>% filter(location %in% LocationsList[reg])
  
  #Filter the data to only use data after the start of the Omicron wave in the current region
  NewOWIDSingleRegionData <- NewOWIDSingleRegionData %>% filter(Date >= Omicron_Start_Dates_Countries[[sub(" ",".", LocationsList[reg])]])
  NewOWIDSingleRegionData <- NewOWIDSingleRegionData %>% filter(Date <= Reporting_End_Dates_Countries[[sub(" ",".", LocationsList[reg])]])
  
  
  #Select the response variables
  Hosp_long <- NewOWIDSingleRegionData %>% dplyr::select(date, hosp_patients)
  names(Hosp_long) <- c("Date", "Hospitalizations")
  Hosp_long$`Hospitalizations Cumulative` <- cumsum(Hosp_long$Hospitalizations)
  
  #Select the predictor variables
  cases_long <- NewOWIDSingleRegionData %>% dplyr::select(date, new_cases, total_cases)
  names(cases_long) <- c("Date","cases", "cumulative_cases")
  
  #Combine the response and predictor variables into a single data frame
  df <-  merge (Hosp_long[,c("Date", "Hospitalizations", "Hospitalizations Cumulative")], 
                cases_long[,c( "cases", "cumulative_cases", "Date" ) ], 
                by = "Date") 
  
  #Reformat the dates in the data frame
  df$Date <- as.Date(df$Date, format="%Y-%m-%d")
  
  #Print the name of the current country to the console
  print(paste(LocationsList[reg], "ICU Patients"))
  
  lcresults <- localpoly(df, alpha = 0.1, method = "lc", min_lead = 5, max_lead = 21, CI = FALSE)
  lcresults_mixed <- localpoly_mixed(df, alpha = 0.1, method = "lc", min_lead = 5, max_lead = 21)
  lcresults_raw <- localpoly_raw(df, alpha = 0.1, method = "lc", min_lead = 5, max_lead = 21)
  
  llresults <- localpoly(df, alpha = 0.1, method = "ll", min_lead = 5, max_lead = 21, CI = FALSE)
  llresults_mixed <- localpoly_mixed(df, alpha = 0.1, method = "ll", min_lead = 5, max_lead = 21)
  llresults_raw <- localpoly_raw(df, alpha = 0.1, method = "ll", min_lead = 5, max_lead = 21)
  
  lprobust_results_lc_bc <- lprobust_search(df, method = "lc", min_lead = 5, max_lead = 21, tau = "bc")
  lprobust_results_lc_bc_mixed <- lprobust_search_mixed(df, method = "lc", min_lead = 5, max_lead = 21, tau = "bc")
  lprobust_results_lc_bc_raw <- lprobust_search_raw(df, method = "lc", min_lead = 5, max_lead = 21, tau = "bc")
  
  lprobust_results_lc_us <- lprobust_search(df, method = "lc", min_lead = 5, max_lead = 21, tau = "us")
  lprobust_results_lc_us_mixed <- lprobust_search_mixed(df, method = "lc", min_lead = 5, max_lead = 21, tau = "us")
  lprobust_results_lc_us_raw <- lprobust_search_raw(df, method = "lc", min_lead = 5, max_lead = 21, tau = "us")
  
  lprobust_results_ll_bc <- lprobust_search(df, method = "ll", min_lead = 5, max_lead = 21, tau = "bc")
  lprobust_results_ll_bc_mixed <- lprobust_search_mixed(df, method = "ll", min_lead = 5, max_lead = 21, tau = "bc")
  lprobust_results_ll_bc_raw <- lprobust_search_raw(df, method = "ll", min_lead = 5, max_lead = 21, tau = "bc")
  
  lprobust_results_ll_us <- lprobust_search(df, method = "ll", min_lead = 5, max_lead = 21, tau = "us")
  lprobust_results_ll_us_mixed <- lprobust_search_mixed(df, method = "ll", min_lead = 5, max_lead = 21, tau = "us")
  lprobust_results_ll_us_raw <- lprobust_search_raw(df, method = "ll", min_lead = 5, max_lead = 21, tau = "us")
  
  
  CountryHospMSPEs[((reg-1)*18)+1,3] <- lcresults$mspe
  CountryHospMSPEs[((reg-1)*18)+2,3] <- lcresults_mixed$mspe
  CountryHospMSPEs[((reg-1)*18)+3,3] <- lcresults_raw$mspe
  CountryHospMSPEs[((reg-1)*18)+4,3] <- llresults$mspe
  CountryHospMSPEs[((reg-1)*18)+5,3] <- llresults_mixed$mspe
  CountryHospMSPEs[((reg-1)*18)+6,3] <- llresults_raw$mspe
  CountryHospMSPEs[((reg-1)*18)+7,3] <- lprobust_results_lc_bc$mspe
  CountryHospMSPEs[((reg-1)*18)+8,3] <- lprobust_results_lc_bc_mixed$mspe
  CountryHospMSPEs[((reg-1)*18)+9,3] <- lprobust_results_lc_bc_raw$mspe
  CountryHospMSPEs[((reg-1)*18)+10,3] <- lprobust_results_lc_us$mspe
  CountryHospMSPEs[((reg-1)*18)+11,3] <- lprobust_results_lc_us_mixed$mspe
  CountryHospMSPEs[((reg-1)*18)+12,3] <- lprobust_results_lc_us_raw$mspe
  CountryHospMSPEs[((reg-1)*18)+13,3] <- lprobust_results_ll_bc$mspe
  CountryHospMSPEs[((reg-1)*18)+14,3] <- lprobust_results_ll_bc_mixed$mspe
  CountryHospMSPEs[((reg-1)*18)+15,3] <- lprobust_results_ll_bc_raw$mspe
  CountryHospMSPEs[((reg-1)*18)+16,3] <- lprobust_results_ll_us$mspe
  CountryHospMSPEs[((reg-1)*18)+17,3] <- lprobust_results_ll_us_mixed$mspe
  CountryHospMSPEs[((reg-1)*18)+18,3] <- lprobust_results_ll_us_raw$mspe
  
  CountryHospMSPEs[((reg-1)*18)+1,4] <- lcresults$lag
  CountryHospMSPEs[((reg-1)*18)+2,4] <- lcresults_mixed$lag
  CountryHospMSPEs[((reg-1)*18)+3,4] <- lcresults_raw$lag
  CountryHospMSPEs[((reg-1)*18)+4,4] <- llresults$lag
  CountryHospMSPEs[((reg-1)*18)+5,4] <- llresults_mixed$lag
  CountryHospMSPEs[((reg-1)*18)+6,4] <- llresults_raw$lag
  CountryHospMSPEs[((reg-1)*18)+7,4] <- lprobust_results_lc_bc$lag
  CountryHospMSPEs[((reg-1)*18)+8,4] <- lprobust_results_lc_bc_mixed$lag
  CountryHospMSPEs[((reg-1)*18)+9,4] <- lprobust_results_lc_bc_raw$lag
  CountryHospMSPEs[((reg-1)*18)+10,4] <- lprobust_results_lc_us$lag
  CountryHospMSPEs[((reg-1)*18)+11,4] <- lprobust_results_lc_us_mixed$lag
  CountryHospMSPEs[((reg-1)*18)+12,4] <- lprobust_results_lc_us_raw$lag
  CountryHospMSPEs[((reg-1)*18)+13,4] <- lprobust_results_ll_bc$lag
  CountryHospMSPEs[((reg-1)*18)+14,4] <- lprobust_results_ll_bc_mixed$lag
  CountryHospMSPEs[((reg-1)*18)+15,4] <- lprobust_results_ll_bc_raw$lag
  CountryHospMSPEs[((reg-1)*18)+16,4] <- lprobust_results_ll_us$lag
  CountryHospMSPEs[((reg-1)*18)+17,4] <- lprobust_results_ll_us_mixed$lag
  CountryHospMSPEs[((reg-1)*18)+18,4] <- lprobust_results_ll_us_raw$lag
}

View(CountryHospMSPEs)

# Nice colour palette: ["#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65", "#beb9db", "#fdcce5", "#8bd3c7"]

ggplot(data=CountryHospMSPEs, aes(fill=Model, y=log(MSPE), x=Location, pattern=Model)) +
  geom_col_pattern(position="dodge", 
                   colour = "black",
                   #pattern_key_scale_factor = 0.5,
                   #pattern_density = 0.1,
                   #pattern_spacing = 0.01,
                   pattern_fill = "black"
  ) + 
  scale_fill_manual(values=rep(c("#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65"),each=3)) +
  scale_pattern_manual(values = rep(c("none", "stripe", "circle"), 6))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  ggtitle("Log of MSPE") +
  theme(plot.title = element_text(hjust = 0.5))

# Model3 <- lm(MSPE ~ Model, data = CountryICUMSPEs)
# anova(Model3)
# summary(Model3)



#### SECTION 3.4 ####
#### Country Results - DF Length Comparison - Hospitalizations (Cumulative X Predicting Raw Y) ####

#Set working directory
#setwd("C:/Users/gamem/Documents/R Files/Research/COVID-Modeling/Plots/Plots27July2023")

# Countries with Hospitalizations Data
LocationsList <- c("United States", "Israel")

CountryHospMSPEs <- data.frame(Location = rep(LocationsList, each = 18),
                               Model = rep(c("tvReg_LC_70","tvReg_LC_140","tvReg_LC_210","tvReg_LL_70","tvReg_LL_140","tvReg_LL_210","LProbust_LC_bc_70","LProbust_LC_bc_140","LProbust_LC_bc_210","LProbust_LC_us_70","LProbust_LC_us_140","LProbust_LC_us_210","LProbust_LL_bc_70","LProbust_LL_bc_140","LProbust_LL_bc_210","LProbust_LL_us_70","LProbust_LL_us_140","LProbust_LL_us_210"),length(LocationsList)),
                               MSPE = rep(NA, length(LocationsList)*18),
                              Lag = rep(NA, length(LocationsList)*18))


for (reg in 1:length(LocationsList)) {
  
  region <- LocationsList[reg]
  
  #Print the name of the current country to the console
  print(paste(LocationsList[reg], "Hospitalizations"))
  
  #Filter the data to a single region
  NewOWIDFullSingleRegionData <- NewOWIDRegionalData %>% filter(location %in% LocationsList[reg])
  
  # List and loop over data frames of different lengths
  DayLengthList <- c(70,140,210)
  
  for (DayIndex in 1:length(DayLengthList)) {
    
    nDays <- DayLengthList[DayIndex]
    
    #Filter the data to only use the correct number of days
    NewOWIDSingleRegionData <- NewOWIDFullSingleRegionData %>% filter(Date >= Omicron_Start_Dates_Countries[[sub(" ",".", LocationsList[reg])]])
    NewOWIDSingleRegionData <- NewOWIDSingleRegionData %>% filter(Date < as.Date(Omicron_Start_Dates_Countries[[sub(" ",".", LocationsList[reg])]]) + nDays)
    
    
    #Select the response variables
    Hosp_long <- NewOWIDSingleRegionData %>% dplyr::select(date, hosp_patients)
    names(Hosp_long) <- c("Date", "Hospitalizations")
    Hosp_long$Hospitalizations[is.na(Hosp_long$Hospitalizations)] <- 0
    Hosp_long$`Hospitalizations Cumulative` <- cumsum(Hosp_long$Hospitalizations)
    
    #Select the predictor variables
    cases_long <- NewOWIDSingleRegionData %>% dplyr::select(date, new_cases, total_cases)
    names(cases_long) <- c("Date","cases", "cumulative_cases")
    
    #Combine the response and predictor variables into a single data frame
    df <-  merge (Hosp_long[,c("Date", "Hospitalizations", "Hospitalizations Cumulative")], 
                  cases_long[,c( "cases", "cumulative_cases", "Date" ) ], 
                  by = "Date") 
    
    #Reformat the dates in the data frame
    df$Date <- as.Date(df$Date, format="%Y-%m-%d")
    
    lcresults <- localpoly_mixed(df, alpha = 0.1, method = "lc", min_lead = 5, max_lead = 21)
    llresults <- localpoly_mixed(df, alpha = 0.1, method = "ll", min_lead = 5, max_lead = 21)
    
    lprobust_results_lc_bc <- lprobust_search_mixed(df, method = "lc", min_lead = 5, max_lead = 21, tau = "bc")
    lprobust_results_lc_us <- lprobust_search_mixed(df, method = "lc", min_lead = 5, max_lead = 21, tau = "us")
    
    lprobust_results_ll_bc <- lprobust_search_mixed(df, method = "ll", min_lead = 5, max_lead = 21, tau = "bc")
    lprobust_results_ll_us <- lprobust_search_mixed(df, method = "ll", min_lead = 5, max_lead = 21, tau = "us")
    
    CountryHospMSPEs[((reg-1)*18)+1+(DayIndex-1),3] <- lcresults$mspe
    CountryHospMSPEs[((reg-1)*18)+4+(DayIndex-1),3] <- llresults$mspe
    CountryHospMSPEs[((reg-1)*18)+7+(DayIndex-1),3] <- lprobust_results_lc_bc$mspe
    CountryHospMSPEs[((reg-1)*18)+10+(DayIndex-1),3] <- lprobust_results_lc_us$mspe
    CountryHospMSPEs[((reg-1)*18)+13+(DayIndex-1),3] <- lprobust_results_ll_bc$mspe
    CountryHospMSPEs[((reg-1)*18)+16+(DayIndex-1),3] <- lprobust_results_ll_us$mspe
    
    CountryHospMSPEs[((reg-1)*18)+1+(DayIndex-1),4] <- lcresults$lag
    CountryHospMSPEs[((reg-1)*18)+4+(DayIndex-1),4] <- llresults$lag
    CountryHospMSPEs[((reg-1)*18)+7+(DayIndex-1),4] <- lprobust_results_lc_bc$lag
    CountryHospMSPEs[((reg-1)*18)+10+(DayIndex-1),4] <- lprobust_results_lc_us$lag
    CountryHospMSPEs[((reg-1)*18)+13+(DayIndex-1),4] <- lprobust_results_ll_bc$lag
    CountryHospMSPEs[((reg-1)*18)+16+(DayIndex-1),4] <- lprobust_results_ll_us$lag
  }
  
}

View(CountryHospMSPEs)

#CountryHospMSPEs7090180

# Nice colour palette: ["#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65", "#beb9db", "#fdcce5", "#8bd3c7"]

#Reorder the columns for graphing
CountryHospMSPEsGraphingDF <- CountryHospMSPEs %>% 
  mutate(Model = fct_relevel(Model, 
                             "tvReg_LC_70","tvReg_LL_70","LProbust_LC_bc_70","LProbust_LC_us_70","LProbust_LL_bc_70","LProbust_LL_us_70",
                             "tvReg_LC_140","tvReg_LL_140","LProbust_LC_bc_140","LProbust_LC_us_140","LProbust_LL_bc_140","LProbust_LL_us_140",
                             "tvReg_LC_210","tvReg_LL_210","LProbust_LC_bc_210","LProbust_LC_us_210","LProbust_LL_bc_210","LProbust_LL_us_210"))


ggplot(data=CountryHospMSPEsGraphingDF, aes(fill=Model, y=log(MSPE), x=Location, pattern=Model)) +
  geom_col_pattern(position="dodge",
                   colour = "black",
                   #pattern_key_scale_factor = 0.5,
                   #pattern_density = 0.1,
                   #pattern_spacing = 0.01,
                   pattern_fill = "black"
  ) +
  scale_fill_manual(values=rep(c("#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65"),times=3)) +
  scale_pattern_manual(values = rep(c("none", "stripe", "circle"), each=6))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  ggtitle("Log of MSPE") +
  theme(plot.title = element_text(hjust = 0.5))



#### Country Results - DF Length Comparison - ICU Patients (Cumulative X Predicting Raw Y) ####

#Set working directory
#setwd("C:/Users/gamem/Documents/R Files/Research/COVID-Modeling/Plots/Plots27July2023")

# Countries with ICU Data
LocationsList <- c("United States", "Israel", "South Korea")

CountryICUMSPEs <- data.frame(Location = rep(LocationsList, each = 18),
                              Model = rep(c("tvReg_LC_70","tvReg_LC_140","tvReg_LC_210","tvReg_LL_70","tvReg_LL_140","tvReg_LL_210","LProbust_LC_bc_70","LProbust_LC_bc_140","LProbust_LC_bc_210","LProbust_LC_us_70","LProbust_LC_us_140","LProbust_LC_us_210","LProbust_LL_bc_70","LProbust_LL_bc_140","LProbust_LL_bc_210","LProbust_LL_us_70","LProbust_LL_us_140","LProbust_LL_us_210"),length(LocationsList)),
                              MSPE = rep(NA, length(LocationsList)*18),
                              Lag = rep(NA, length(LocationsList)*18))

for (reg in 1:length(LocationsList)) {
  
  region <- LocationsList[reg]
  
  #Print the name of the current country to the console
  print(paste(LocationsList[reg], "ICU Patients"))
  
  #Filter the data to a single region
  NewOWIDFullSingleRegionData <- NewOWIDRegionalData %>% filter(location %in% LocationsList[reg])
  
  # List and loop over data frames of different lengths
  DayLengthList <- c(70,140,210)
  
  for (DayIndex in 1:length(DayLengthList)) {
    
    nDays <- DayLengthList[DayIndex]
    
    #Filter the data to only use the correct number of days
    NewOWIDSingleRegionData <- NewOWIDFullSingleRegionData %>% filter(Date >= Omicron_Start_Dates_Countries[[sub(" ",".", LocationsList[reg])]])
    NewOWIDSingleRegionData <- NewOWIDSingleRegionData %>% filter(Date < as.Date(Omicron_Start_Dates_Countries[[sub(" ",".", LocationsList[reg])]]) + nDays)
    
    
    #Select the response variables
    ICU_long <- NewOWIDSingleRegionData %>% dplyr::select(date, icu_patients)
    names(ICU_long) <- c("Date", "ICU Patients")
    ICU_long$`ICU Patients`[is.na(ICU_long$`ICU Patients`)] <- 0
    ICU_long$`ICU Patients Cumulative` <- cumsum(ICU_long$`ICU Patients`)
    
    #Select the predictor variables
    cases_long <- NewOWIDSingleRegionData %>% dplyr::select(date, new_cases, total_cases)
    names(cases_long) <- c("Date","cases", "cumulative_cases")
    
    #Combine the response and predictor variables into a single data frame
    df <-  merge (ICU_long[,c("Date", "ICU Patients", "ICU Patients Cumulative")], 
                  cases_long[,c( "cases", "cumulative_cases", "Date" ) ], 
                  by = "Date") 
    
    #Reformat the dates in the data frame
    df$Date <- as.Date(df$Date, format="%Y-%m-%d")
    
    lcresults <- localpoly_mixed(df, alpha = 0.1, method = "lc", min_lead = 5, max_lead = 21)
    llresults <- localpoly_mixed(df, alpha = 0.1, method = "ll", min_lead = 5, max_lead = 21)
    
    lprobust_results_lc_bc <- lprobust_search_mixed(df, method = "lc", min_lead = 5, max_lead = 21, tau = "bc")
    lprobust_results_lc_us <- lprobust_search_mixed(df, method = "lc", min_lead = 5, max_lead = 21, tau = "us")
    
    lprobust_results_ll_bc <- lprobust_search_mixed(df, method = "ll", min_lead = 5, max_lead = 21, tau = "bc")
    lprobust_results_ll_us <- lprobust_search_mixed(df, method = "ll", min_lead = 5, max_lead = 21, tau = "us")
    
    CountryICUMSPEs[((reg-1)*18)+1+(DayIndex-1),3] <- lcresults$mspe
    CountryICUMSPEs[((reg-1)*18)+4+(DayIndex-1),3] <- llresults$mspe
    CountryICUMSPEs[((reg-1)*18)+7+(DayIndex-1),3] <- lprobust_results_lc_bc$mspe
    CountryICUMSPEs[((reg-1)*18)+10+(DayIndex-1),3] <- lprobust_results_lc_us$mspe
    CountryICUMSPEs[((reg-1)*18)+13+(DayIndex-1),3] <- lprobust_results_ll_bc$mspe
    CountryICUMSPEs[((reg-1)*18)+16+(DayIndex-1),3] <- lprobust_results_ll_us$mspe
    
    CountryICUMSPEs[((reg-1)*18)+1+(DayIndex-1),4] <- lcresults$lag
    CountryICUMSPEs[((reg-1)*18)+4+(DayIndex-1),4] <- llresults$lag
    CountryICUMSPEs[((reg-1)*18)+7+(DayIndex-1),4] <- lprobust_results_lc_bc$lag
    CountryICUMSPEs[((reg-1)*18)+10+(DayIndex-1),4] <- lprobust_results_lc_us$lag
    CountryICUMSPEs[((reg-1)*18)+13+(DayIndex-1),4] <- lprobust_results_ll_bc$lag
    CountryICUMSPEs[((reg-1)*18)+16+(DayIndex-1),4] <- lprobust_results_ll_us$lag
  }
  
}

#CountryICUMSPEs7090180

View(CountryICUMSPEs)

# Nice colour palette: ["#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65", "#beb9db", "#fdcce5", "#8bd3c7"]

#Reorder the columns for graphing
CountryICUMSPEsGraphingDF <- CountryICUMSPEs %>% 
  mutate(Model = fct_relevel(Model, 
                             "tvReg_LC_70","tvReg_LL_70","LProbust_LC_bc_70","LProbust_LC_us_70","LProbust_LL_bc_70","LProbust_LL_us_70",
                             "tvReg_LC_140","tvReg_LL_140","LProbust_LC_bc_140","LProbust_LC_us_140","LProbust_LL_bc_140","LProbust_LL_us_140",
                             "tvReg_LC_210","tvReg_LL_210","LProbust_LC_bc_210","LProbust_LC_us_210","LProbust_LL_bc_210","LProbust_LL_us_210"))


ggplot(data=CountryICUMSPEsGraphingDF, aes(fill=Model, y=log(MSPE), x=Location, pattern=Model)) +
  geom_col_pattern(position="dodge",
                   colour = "black",
                   #pattern_key_scale_factor = 0.5,
                   #pattern_density = 0.1,
                   #pattern_spacing = 0.01,
                   pattern_fill = "black"
  ) +
  scale_fill_manual(values=rep(c("#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65"),times=3)) +
  scale_pattern_manual(values = rep(c("none", "stripe", "circle"), each=6))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  ggtitle("Log of MSPE") +
  theme(plot.title = element_text(hjust = 0.5))



# ggplot(data=CountryICUMSPEs, aes(fill=Model, y=log(MSPE), x=Location, pattern=Model)) +
#   geom_col_pattern(position="dodge", 
#                    colour = "black",
#                    pattern_fill = "black",
#                    pattern_density = 0.1,
#                    pattern_spacing = 0.01,
#                    pattern_key_scale_factor = 0.5) + 
#   scale_fill_manual(values=rep(c("#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65"),each=3)) +
#   scale_pattern_manual(values = rep(c("none", "stripe", "circle"), 6))+
#   scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
#   ggtitle("Log of MSPE") +
#   theme(plot.title = element_text(hjust = 0.5))


# Model3 <- lm(MSPE ~ Model, data = CountryICUMSPEs)
# anova(Model3)
# summary(Model3)


#### EXTRA ####
#### Country Results - DF Length Comparison - Deaths (Cumulative X Predicting Cumulative Y) ####

#Set working directory
#setwd("C:/Users/gamem/Documents/R Files/Research/COVID-Modeling/Plots/Plots27July2023")

# Countries with Deaths data
LocationsList <- c("Canada", "United States", "Israel", "South Africa", "South Korea")

CountryDeathsMSPEs <- data.frame(Location = rep(LocationsList, each = 18),
                                 Model = rep(c("tvReg_LC_70","tvReg_LC_140","tvReg_LC_210","tvReg_LL_70","tvReg_LL_140","tvReg_LL_210","LProbust_LC_bc_70","LProbust_LC_bc_140","LProbust_LC_bc_210","LProbust_LC_us_70","LProbust_LC_us_140","LProbust_LC_us_210","LProbust_LL_bc_70","LProbust_LL_bc_140","LProbust_LL_bc_210","LProbust_LL_us_70","LProbust_LL_us_140","LProbust_LL_us_210"),length(LocationsList)),
                                 MSPE = rep(NA, length(LocationsList)*18),
                                 Lag = rep(NA, length(LocationsList)*18))

for (reg in 1:length(LocationsList)) {
  
  region <- LocationsList[reg]
  
  #Print the name of the current country to the console
  print(paste(LocationsList[reg], "Deaths"))
  
  #Filter the data to a single region
  NewOWIDFullSingleRegionData <- NewOWIDRegionalData %>% filter(location %in% LocationsList[reg])
  
  # List and loop over data frames of different lengths
  DayLengthList <- c(70,140,210)
  
  for (DayIndex in 1:length(DayLengthList)) {
    
    nDays <- DayLengthList[DayIndex]
    
    #Filter the data to only use the correct number of days
    NewOWIDSingleRegionData <- NewOWIDFullSingleRegionData %>% filter(Date >= Omicron_Start_Dates_Countries[[sub(" ",".", LocationsList[reg])]])
    NewOWIDSingleRegionData <- NewOWIDSingleRegionData %>% filter(Date < as.Date(Omicron_Start_Dates_Countries[[sub(" ",".", LocationsList[reg])]]) + nDays)
    
    
    #Select the response variables
    Deaths_long <- NewOWIDSingleRegionData %>% dplyr::select(date, new_deaths)
    names(Deaths_long) <- c("Date", "Deaths")
    Deaths_long$Deaths[is.na(Deaths_long$Deaths)] <- 0
    Deaths_long$`Deaths Cumulative` <- cumsum(Deaths_long$Deaths)
    
    #Select the predictor variables
    cases_long <- NewOWIDSingleRegionData %>% dplyr::select(date, new_cases, total_cases)
    names(cases_long) <- c("Date","cases", "cumulative_cases")
    
    #Combine the response and predictor variables into a single data frame
    df <-  merge (Deaths_long[,c("Date", "Deaths", "Deaths Cumulative")], 
                  cases_long[,c( "cases", "cumulative_cases", "Date" ) ], 
                  by = "Date") 
    
    #Reformat the dates in the data frame
    df$Date <- as.Date(df$Date, format="%Y-%m-%d")
    
    lcresults <- localpoly(df, alpha = 0.1, method = "lc", min_lead = 5, max_lead = 21, CI = FALSE)
    llresults <- localpoly(df, alpha = 0.1, method = "ll", min_lead = 5, max_lead = 21, CI = FALSE)
    
    lprobust_results_lc_bc <- lprobust_search(df, method = "lc", min_lead = 5, max_lead = 21, tau = "bc")
    lprobust_results_lc_us <- lprobust_search(df, method = "lc", min_lead = 5, max_lead = 21, tau = "us")
    
    lprobust_results_ll_bc <- lprobust_search(df, method = "ll", min_lead = 5, max_lead = 21, tau = "bc")
    lprobust_results_ll_us <- lprobust_search(df, method = "ll", min_lead = 5, max_lead = 21, tau = "us")
    
    CountryDeathsMSPEs[((reg-1)*18)+1+(DayIndex-1),3] <- lcresults$mspe
    CountryDeathsMSPEs[((reg-1)*18)+4+(DayIndex-1),3] <- llresults$mspe
    CountryDeathsMSPEs[((reg-1)*18)+7+(DayIndex-1),3] <- lprobust_results_lc_bc$mspe
    CountryDeathsMSPEs[((reg-1)*18)+10+(DayIndex-1),3] <- lprobust_results_lc_us$mspe
    CountryDeathsMSPEs[((reg-1)*18)+13+(DayIndex-1),3] <- lprobust_results_ll_bc$mspe
    CountryDeathsMSPEs[((reg-1)*18)+16+(DayIndex-1),3] <- lprobust_results_ll_us$mspe
    
    CountryDeathsMSPEs[((reg-1)*18)+1+(DayIndex-1),4] <- lcresults$lag
    CountryDeathsMSPEs[((reg-1)*18)+4+(DayIndex-1),4] <- llresults$lag
    CountryDeathsMSPEs[((reg-1)*18)+7+(DayIndex-1),4] <- lprobust_results_lc_bc$lag
    CountryDeathsMSPEs[((reg-1)*18)+10+(DayIndex-1),4] <- lprobust_results_lc_us$lag
    CountryDeathsMSPEs[((reg-1)*18)+13+(DayIndex-1),4] <- lprobust_results_ll_bc$lag
    CountryDeathsMSPEs[((reg-1)*18)+16+(DayIndex-1),4] <- lprobust_results_ll_us$lag
  }
  
}

View(CountryDeathsMSPEs)

# Nice colour palette: ["#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65", "#beb9db", "#fdcce5", "#8bd3c7"]


#Reorder the columns for graphing
CountryDeathsMSPEsGraphingDF <- CountryDeathsMSPEs %>% 
  mutate(Model = fct_relevel(Model, 
                             "tvReg_LC_70","tvReg_LL_70","LProbust_LC_bc_70","LProbust_LC_us_70","LProbust_LL_bc_70","LProbust_LL_us_70",
                             "tvReg_LC_140","tvReg_LL_140","LProbust_LC_bc_140","LProbust_LC_us_140","LProbust_LL_bc_140","LProbust_LL_us_140",
                             "tvReg_LC_210","tvReg_LL_210","LProbust_LC_bc_210","LProbust_LC_us_210","LProbust_LL_bc_210","LProbust_LL_us_210"))


ggplot(data=CountryDeathsMSPEsGraphingDF, aes(fill=Model, y=log(MSPE), x=Location, pattern=Model)) +
  geom_col_pattern(position="dodge",
                   colour = "black",
                   #pattern_key_scale_factor = 0.5,
                   #pattern_density = 0.1,
                   #pattern_spacing = 0.01,
                   pattern_fill = "black"
  ) +
  scale_fill_manual(values=rep(c("#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65"),times=3)) +
  scale_pattern_manual(values = rep(c("none", "stripe", "circle"), each=6))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  ggtitle("Log of MSPE") +
  theme(plot.title = element_text(hjust = 0.5))




