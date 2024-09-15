# Master's Thesis: Time Series Analysis with VAR Models

###########################
# Load necessary libraries
###########################

library(readxl)       # To read Excel files
library(tseries)      # For time series data handling
library(ggplot2)      # For data visualization
library(gridExtra)    # Arrange multiple ggplot2 plots in a grid layout
library(ggfortify)    # For visualizing time series objects with ggplot2
library(forecast)     # For forecasting time series
library(tidyr)        # For data tidying
library(urca)         # For unit root and cointegration tests
library(vars)         # For VAR models
library(svars)        # For structural VAR analysis
library(lubridate)    # For date and time manipulation
library(bvarsv)       # For time-varying parameter VAR (not used due to performance issues)
library(Matrix)       # Matrix computations
library(plotly)       # Interactive plotting
library(glue)         # String interpolation
library(scales)       # Scaling functions for plots
library(viridis)      # Color palettes for visualization
library(ggthemes)     # Additional themes for ggplot2

#######################
# Functions definition
#######################

# Helper Functions

# Function to create a styled plot

create_styled_plot <- function(data, x, y, title, y_label, y_suffix = "") {
  ggplot(data, aes(x = !!sym(x), y = !!sym(y))) +
    geom_line(color = "steelblue", linewidth = 1.2) +
    labs(title = title,
         x = "Year",
         y = y_label,
         caption = "Source: BDREMS") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "ghostwhite"),
      legend.position = "none",
      plot.margin = margin(20, 20, 20, 20)
    ) +
    scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
    scale_y_continuous(labels = function(x) paste0(comma(x), y_suffix))
}

theme_econ <- function(base_size = 12, base_family = "Helvetica") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_text(hjust = 0.5, size = rel(0.9), margin = margin(b = 10)),
      plot.caption = element_text(hjust = 1, size = rel(0.8), margin = margin(t = 10)),
      axis.title = element_text(face = "bold", size = rel(0.9)),
      axis.text = element_text(size = rel(0.8)),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = rel(0.8)),
      panel.grid.major = element_line(color = "gray90", size = 0.1),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 0.5),
      plot.margin = margin(20, 20, 20, 20)
    )
}

# Function to create ACF and PACF plots

create_acf_pacf_plot <- function(series, title) {
  acf_data <- acf(series, plot = FALSE)
  pacf_data <- pacf(series, plot = FALSE)
  
  acf_df <- data.frame(lag = acf_data$lag, acf = acf_data$acf)
  pacf_df <- data.frame(lag = pacf_data$lag, pacf = pacf_data$acf)
  
  acf_plot <- ggplot(acf_df, aes(x = lag, y = acf)) +
    geom_bar(stat = "identity", fill = "#1f77b4", width = 0.5) +
    geom_hline(yintercept = c(0, qnorm(c(0.025, 0.975)) / sqrt(length(series))), 
               linetype = "dashed", color = "red") +
    labs(title = "ACF", x = "Lag", y = "ACF") +
    theme_econ()
  
  pacf_plot <- ggplot(pacf_df, aes(x = lag, y = pacf)) +
    geom_bar(stat = "identity", fill = "#1f77b4", width = 0.5) +
    geom_hline(yintercept = c(0, qnorm(c(0.025, 0.975)) / sqrt(length(series))), 
               linetype = "dashed", color = "red") +
    labs(title = "PACF", x = "Lag", y = "PACF") +
    theme_econ()
  
  grid.arrange(acf_plot, pacf_plot, ncol = 2, top = title)
}

############################ SVAR FUNCTIONS ####################################

#' Generate Q matrices
#'
#' @param n Integer, dimension of the matrix
#' @return A randomly generated Q matrix of size n x n.

generate_Q <- function(n) {
  qr.Q(qr(matrix(rnorm(n * n), ncol = n))) # Generate Q matrix using QR decomposition
}

#' Compute IRFs for SVAR
#'
#' @param svar_model SVAR model object
#' @param b_matrix B matrix
#' @param horizon Integer, forecast horizon
#' @return IRF results as a list with the same structure as vars::irf outpu

compute_irfs_svar <- function(svar_model, B_matrix, horizon) {
  if (horizon == 0) {
    # Extract the varresult list
    varresult <- svar_model$var$varresult 
    
    # Number of endogenous variables
    k <- length(varresult)
    
    # Number of lags
    p <- svar_model$var$p
    
    # Initialize the sum of A matrices
    A_sum <- matrix(0, nrow = k, ncol = k)
    
    # Reconstruct the A matrices from the coefficients in varresult
    for (i in 1:p) {
      A_matrix <- matrix(0, nrow = k, ncol = k)
      for (j in 1:k) {
        # Extract the coefficients for the j-th variable
        coefficients <- as.vector(varresult[[j]]$coefficients)
        # Extract the coefficients for the i-th lag, skipping the constant term
        lag_indices <- (i - 1) * k + 1:k
        lag_coefficients <- coefficients[lag_indices]
        A_matrix[j, ] <- lag_coefficients
      }
      A_sum <- A_sum + A_matrix
    }
    
    # Compute the long-term impact matrix
    long_term_impact <- solve(diag(k) - A_sum) %*% B_matrix
    
    # Create a list with the same structure as vars::irf output
    irf_list <- list()
    var_names <- names(varresult)
    for (var in var_names) {
      irf_list[[var]] <- matrix(long_term_impact[, which(var_names == var)], nrow = 1, ncol = k)
      colnames(irf_list[[var]]) <- var_names
    }
    
    irf_result <- list(
      irf = irf_list,
      Lower = NULL,
      Upper = NULL,
      response = var_names,
      impulse = var_names,
      ortho = TRUE,
      cumulative = FALSE,
      runs = 1,
      ci = 0.95,
      boot = FALSE,
      model = "svarest"
    )
    
    class(irf_result) <- "varirf"
    
    return(irf_result)
  } else {
    # Compute IRFs using the provided B matrix
    svar_model$B <- B_matrix  # Update model with new B matrix
    return(vars::irf(svar_model, n.ahead = horizon, boot = FALSE))
  }
}

#' Check sign restrictions
#'
#' @param irf_result IRF result object from SVAR
#' @param sign_restrictions Matrix of sign restrictions (+1, -1, NA)
#' @param periods_to_check Integer, number of periods to check sign restrictions
#' @return Boolean, TRUE if restrictions are satisfied

check_sign_restrictions_svar <- function(irf_result, sign_restrictions, periods_to_check) {
  for (shock_index in seq_len(ncol(sign_restrictions))) {
    for (var_index in seq_len(nrow(sign_restrictions))) {
      irf_series <- irf_result$irf[[shock_index]][, var_index]
      expected_sign <- sign_restrictions[var_index, shock_index]
      
      # Ensure only +1, -1, and NA are allowed as sign restrictions
      if (!expected_sign %in% c(1, -1, NA)) {
        stop("Sign restrictions must be +1, -1, or NA")
      }
      
      if (!all(sign(irf_series[1:periods_to_check]) == expected_sign | is.na(expected_sign))) {
        return(FALSE)  
      }
    }
  }
  return(TRUE)  
}

#' Compute summary statistics of IRFs
#'
#' @param irf_list List of IRF results from multiple draws
#' @param percentiles Numeric vector, percentiles to compute (default is 16th and 84th percentiles)
#' @return A list containing mean, median, and percentile IRFs for each shock
#' @export

compute_irf_summary_svar <- function(irf_list, percentiles = c(16, 84)) {
  # Extract the number of draws, shocks, periods, and variables
  n_draws <- length(irf_list)
  n_shocks <- length(irf_list[[1]])
  n_periods <- nrow(irf_list[[1]][[1]])
  n_vars <- ncol(irf_list[[1]][[1]])
  
  # Initialize lists to store the results
  mean_irf <- vector("list", n_shocks)
  median_irf <- vector("list", n_shocks)
  percentile_irfs <- vector("list", n_shocks)
  
  # Iterate over each shock
  for (shock_idx in seq_len(n_shocks)) {
    mean_irf[[shock_idx]] <- matrix(0, n_periods, n_vars)
    median_irf[[shock_idx]] <- matrix(0, n_periods, n_vars)
    percentile_irfs[[shock_idx]] <- array(0, dim = c(n_periods, n_vars, length(percentiles)))
    
    # Convert each shock's IRFs across the list into a 3D array
    irf_array <- array(NA, dim = c(n_periods, n_vars, n_draws))
    for (draw_idx in seq_along(irf_list)) {
      irf_array[,,draw_idx] <- irf_list[[draw_idx]][[shock_idx]]
    }
    
    # Compute mean, median, and percentiles for the current shock
    mean_irf[[shock_idx]] <- apply(irf_array, c(1, 2), mean, na.rm = TRUE)
    median_irf[[shock_idx]] <- apply(irf_array, c(1, 2), median, na.rm = TRUE)
    
    # Compute percentiles
    for (p_idx in seq_along(percentiles)) {
      percentile_irfs[[shock_idx]][,,p_idx] <- apply(irf_array, c(1, 2), quantile, probs = percentiles[p_idx] / 100, na.rm = TRUE)
    }
  }
  
  return(list(mean = mean_irf, median = median_irf, percentiles = percentile_irfs))
}

#' Plot Scaled IRF with Professional Formatting and Custom Y-Axis Limits
#'
#' This function generates plots for Impulse Response Functions (IRFs) from VAR or SVAR models,
#' scaling the results by 100 and allowing customization of the summary statistic, confidence bounds, and y-axis limits.
#'
#' @param irf_summary Summary of IRF results from `compute_irf_summary_svar`.
#' @param var_index Integer, index of the variable to plot.
#' @param shock_index Integer, index of the shock to plot.
#' @param summary_type String, either "mean" or "median" to specify the summary statistic.
#' @param percentiles Numeric vector of two percentiles for confidence bounds.
#' @param plot_title String, title of the plot (optional).
#' @param x_axis_label String, label for the x-axis (default is "Time Horizon").
#' @param y_axis_label String, label for the y-axis (default is "IRF").
#' @param y_min Numeric, minimum value for the y-axis (optional).
#' @param y_max Numeric, maximum value for the y-axis (optional).
#' @return A ggplot object representing the IRF plot.
#' @export

plot_irf_svar <- function(irf_summary, var_index = 1, shock_index = 1, 
                          summary_type = "mean", percentiles = c(16, 84),
                          plot_title = NULL, x_axis_label = "Time Horizon", 
                          y_axis_label = "IRF", y_min = NULL, y_max = NULL) {
  
  # Ensure the variable index is within range
  if (var_index > ncol(irf_summary$mean[[shock_index]]) || var_index < 1) {
    stop("Variable index out of range")
  }
  
  # Ensure the shock index is within range
  if (shock_index > length(irf_summary$mean) || shock_index < 1) {
    stop("Shock index out of range")
  }
  
  # Extract time horizon length
  time_horizon_length <- nrow(irf_summary$mean[[shock_index]])
  time_horizon <- 1:time_horizon_length
  
  # Extract the central IRF and scale by 100
  if (summary_type == "mean") {
    central_irf <- irf_summary$mean[[shock_index]][, var_index] * 100
  } else if (summary_type == "median") {
    central_irf <- irf_summary$median[[shock_index]][, var_index] * 100
  } else {
    stop("Invalid summary_type. Use 'mean' or 'median'.")
  }
  
  # Extract and scale confidence bounds
  lower_percentile <- percentiles[1]
  upper_percentile <- percentiles[2]
  
  lower_bound <- irf_summary$percentiles[[shock_index]][, var_index, which(percentiles == lower_percentile)]*100
  upper_bound <- irf_summary$percentiles[[shock_index]][, var_index, which(percentiles == upper_percentile)]*100
  
  # Check if the lengths match
  if (length(time_horizon) != length(central_irf) || 
      length(central_irf) != length(lower_bound) || 
      length(lower_bound) != length(upper_bound)) {
    stop("Length mismatch between time horizon and IRF data")
  }
  
  # Create a data frame for plotting
  irf_df <- data.frame(
    Time = time_horizon,
    CentralIRF = central_irf,
    LowerBound = lower_bound,
    UpperBound = upper_bound
  )
  
  # Set a default title if none is provided
  if (is.null(plot_title)) {
    plot_title <- paste(summary_type, "IRF for Variable", var_index, "Shock", shock_index)
  }
  
  # Create the plot using ggplot2
  p <- ggplot(irf_df, aes(x = Time)) +
    geom_line(aes(y = CentralIRF), color = "#0072B2", size = 1.2) +  # Professional blue color, thicker line
    geom_ribbon(aes(ymin = LowerBound, ymax = UpperBound), fill = "#0072B2", alpha = 0.2) +  # Semi-transparent confidence interval
    labs(title = plot_title, y = y_axis_label, x = x_axis_label) +
    theme_minimal(base_size = 15) +  # Clean, professional theme
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),  # Centered and bold title
      axis.title.x = element_text(face = "bold", size = 12),  # Bold and slightly larger x-axis label
      axis.title.y = element_text(face = "bold", size = 12),  # Bold and slightly larger y-axis label
      axis.text = element_text(size = 10),  # Standardized axis text size
      panel.grid.major = element_line(color = "grey80", size = 0.5),  # Major grid lines
      panel.grid.minor = element_line(color = "grey90", size = 0.25)  # Minor grid lines
    ) +
    scale_y_continuous(labels = scales::percent_format(scale = 1))  # Format y-axis as percentage
  
  # Add custom y-axis limits if specified
  if (!is.null(y_min) & !is.null(y_max)) {
    p <- p + coord_cartesian(ylim = c(y_min, y_max))
  }
  
  return(p)
}

############################ TVP VAR FUNCTIONS #################################

# TVP-VAR Functions

#' Check Sign Restrictions for TVP-VAR B Matrix
#'
#' This function checks if a candidate B matrix satisfies the given sign restrictions. 
#' Sign restrictions are expected to be +1, -1, or NA, where +1 indicates a positive effect, 
#' -1 indicates a negative effect, and NA indicates no restriction.
#'
#' @param B_candidate A matrix of candidate B values to be checked against the sign restrictions.
#' @param sign_restrictions A matrix of sign restrictions where +1, -1, or NA indicate the expected sign of the elements in the B matrix.
#'
#' @return A boolean value, `TRUE` if all sign restrictions are satisfied, `FALSE` otherwise.
#'
#' @export

check_sign_restrictions_tvp <- function(B_candidate, sign_restrictions) {
  satisfied <- TRUE
  
  for (var_index in seq_len(nrow(sign_restrictions))) {
    for (shock_index in seq_len(ncol(sign_restrictions))) {
      expected_sign <- sign_restrictions[var_index, shock_index]
      
      # Ensure only +1, -1, and NA are allowed as sign restrictions
      if (!expected_sign %in% c(1, -1, NA)) {
        stop("Sign restrictions must be +1, -1, or NA")
      }
      
      contemporaneous_effect <- B_candidate[var_index, shock_index]
      
      # Check if the contemporaneous effect meets the expected sign
      if (!(sign(contemporaneous_effect) == expected_sign | is.na(expected_sign))) {
        satisfied <- FALSE
        break
      }
    }
    if (!satisfied) break
  }
  
  if (satisfied) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Calculate IRF for TVP-VAR
#'
#' This function calculates Impulse Response Functions (IRFs) for a Time-Varying Parameter VAR (TVP-VAR) model. 
#' It computes the IRFs over a specified forecast horizon and period, using provided draws for the Beta and B matrices.
#'
#' @param n_variables Integer, the number of variables in the VAR model.
#' @param p Integer, the number of lags in the VAR model.
#' @param n_periods Integer, the number of periods for which IRFs are calculated.
#' @param n_reps Integer, the number of repetitions or draws for the Beta and B matrices.
#' @param Beta_draws Array of Beta draws, with dimensions [number of variables + p * number of variables^2 + number of variables, number of periods, number of repetitions].
#' @param B_draws Array of B draws, with dimensions [number of variables, number of variables * number of periods, number of repetitions].
#' @param horizon Integer, the forecast horizon for the IRFs.
#' @param use_intercept Boolean, whether to include intercepts in the IRF calculation. Default is TRUE.
#'
#' @return An array of IRF results with dimensions [number of variables, number of variables, horizon + 1, number of periods, number of repetitions].
#'
#' @export

calculate_irf_tvp <- function(n_variables, p, n_periods, n_reps, Beta_draws, B_draws, horizon, use_intercept = TRUE) {
  irf_results <- array(0, dim = c(n_variables, n_variables, horizon + 1, n_periods, n_reps))
  
  for (t in 1:n_periods) {
    for (i in 1:n_reps) {
      B_current <- B_draws[, ((t - 1) * n_variables + 1):(t * n_variables), i]
      Beta_current <- Beta_draws[, t, i]
      
      intercept <- Beta_current[1:n_variables]
      
      coefficient_matrices <- array(Beta_current[(n_variables + 1):(n_variables + p * n_variables^2 + n_variables)], 
                                    dim = c(n_variables, n_variables, p))
      
      # Transpose each coefficient matrix for every lag
      for (j in 1:p) {
        coefficient_matrices[,,j] <- t(coefficient_matrices[,,j])
      }
      
      irf_current <- array(0, dim = c(n_variables, n_variables, horizon + 1))
      irf_current[,,1] <- B_current
      
      for (h in 2:(horizon + 1)) {
        temp_response <- array(0, dim = c(n_variables, n_variables))
        
        for (j in 1:min(h - 1, p)) {
          temp_response <- temp_response + coefficient_matrices[,,j] %*% irf_current[,,h - j]
        }
        
        if (use_intercept) {
          temp_response <- temp_response + matrix(intercept, nrow = n_variables, ncol = n_variables, byrow = TRUE)
        }
        
        irf_current[,,h] <- temp_response
      }
      
      irf_results[, , , t, i] <- irf_current
    }
  }
  
  return(irf_results)
}

#' Compute Summary Statistics for IRFs
#'
#' This function computes summary statistics for Impulse Response Functions (IRFs) from a TVP-VAR model. 
#' It calculates the mean, median, and specified percentiles of the IRFs across all periods and repetitions.
#'
#' @param irf_array Array of IRF results with dimensions [number of variables, number of shocks, horizon + 1, number of periods, number of repetitions].
#' @param percentiles Vector of percentiles (in percentage) for which to compute the IRF bounds. Default is c(16, 84).
#'
#' @return A list containing:
#' \item{mean}{Array of mean IRFs with dimensions [number of variables, number of shocks, horizon + 1, number of periods].}
#' \item{median}{Array of median IRFs with dimensions [number of variables, number of shocks, horizon + 1, number of periods].}
#' \item{percentiles}{A list of arrays, each containing IRFs at specified percentiles with dimensions [number of variables, number of shocks, horizon + 1, number of periods].}
#' \item{overall_mean}{Array of overall mean IRFs across all periods with dimensions [number of variables, number of shocks, horizon + 1].}
#' \item{overall_median}{Array of overall median IRFs across all periods with dimensions [number of variables, number of shocks, horizon + 1].}
#' \item{overall_percentiles}{A list of arrays, each containing overall IRFs at specified percentiles with dimensions [number of variables, number of shocks, horizon + 1].}
#'
#' @export

compute_irf_summary_tvp <- function(irf_array, percentiles = c(16, 84)) {
  n_vars <- dim(irf_array)[1]
  n_shocks <- dim(irf_array)[2]
  horizon <- dim(irf_array)[3] - 1
  n_periods <- dim(irf_array)[4]
  n_reps <- dim(irf_array)[5]
  
  mean_irf <- array(0, dim = c(n_vars, n_shocks, horizon + 1, n_periods))
  median_irf <- array(0, dim = c(n_vars, n_shocks, horizon + 1, n_periods))
  percentile_irfs <- vector("list", length(percentiles))
  
  for (p in 1:length(percentiles)) {
    percentile_irfs[[p]] <- array(0, dim = c(n_vars, n_shocks, horizon + 1, n_periods))
  }
  
  overall_mean_irf <- array(0, dim = c(n_vars, n_shocks, horizon + 1))
  overall_median_irf <- array(0, dim = c(n_vars, n_shocks, horizon + 1))
  overall_percentile_irfs <- vector("list", length(percentiles))
  
  for (p in 1:length(percentiles)) {
    overall_percentile_irfs[[p]] <- array(0, dim = c(n_vars, n_shocks, horizon + 1))
  }
  
  for (t in 1:n_periods) {
    for (v in 1:n_vars) {
      for (s in 1:n_shocks) {
        for (h in 1:(horizon + 1)) {
          irf_data <- irf_array[v, s, h, t, ]
          
          mean_irf[v, s, h, t] <- mean(irf_data, na.rm = TRUE)
          median_irf[v, s, h, t] <- median(irf_data, na.rm = TRUE)
          
          for (p in 1:length(percentiles)) {
            percentile_irfs[[p]][v, s, h, t] <- quantile(irf_data, probs = percentiles[p] / 100, na.rm = TRUE)
          }
        }
      }
    }
  }
  
  for (v in 1:n_vars) {
    for (s in 1:n_shocks) {
      for (h in 1:(horizon + 1)) {
        overall_irf_data <- as.vector(irf_array[v, s, h, , ])
        
        overall_mean_irf[v, s, h] <- mean(overall_irf_data, na.rm = TRUE)
        overall_median_irf[v, s, h] <- median(overall_irf_data, na.rm = TRUE)
        
        for (p in 1:length(percentiles)) {
          overall_percentile_irfs[[p]][v, s, h] <- quantile(overall_irf_data, probs = percentiles[p] / 100, na.rm = TRUE)
        }
      }
    }
  }
  
  return(list(
    mean = mean_irf,
    median = median_irf,
    percentiles = percentile_irfs,
    overall_mean = overall_mean_irf,
    overall_median = overall_median_irf,
    overall_percentiles = overall_percentile_irfs
  ))
}

#' Plot Scaled IRF for TVP-VAR with Professional Formatting
#'
#' This function generates plots for Impulse Response Functions (IRFs) from a Time-Varying Parameter VAR (TVP-VAR) model. 
#' It scales the IRF results by a specified factor (default 100) and provides customizable options for plotting, including variable, shock, and period indices, 
#' summary statistics, and confidence bounds.
#'
#' @param irf_summary Summary of IRF results, including mean, median, and percentiles.
#' @param var_index Integer, index of the variable to plot.
#' @param shock_index Integer, index of the shock to plot.
#' @param period_index Integer, index of the period to plot. Use 0 for overall IRF across periods.
#' @param summary_type String, type of summary statistic to plot ("mean" or "median").
#' @param percentiles Vector of two percentiles for confidence bounds (e.g., c(16, 84)).
#' @param plot_title String, title of the plot (optional).
#' @param x_axis_label String, label for the x-axis (default is "Time Horizon").
#' @param y_axis_label String, label for the y-axis (default is "IRF (%)").
#' @param scale_factor Numeric, factor by which to scale the IRF results (default is 100).
#'
#' @return A `ggplot` object of the IRF plot.
#'
#' @export

plot_irf_tvp <- function(irf_summary, var_index = 1, shock_index = 1, 
                         period_index = 1, summary_type = "mean", 
                         percentiles = c(16, 84), plot_title = NULL, 
                         x_axis_label = "Time Horizon", y_axis_label = "IRF (%)",
                         scale_factor = 100) {
  
  # Ensure the indices are within range
  if (var_index > dim(irf_summary$mean)[1] || var_index < 1) {
    stop("Variable index out of range")
  }
  if (shock_index > dim(irf_summary$mean)[2] || shock_index < 1) {
    stop("Shock index out of range")
  }
  if (period_index > dim(irf_summary$mean)[4] || period_index < 0) {
    stop("Period index out of range")
  }
  
  # Extract time horizon length
  time_horizon_length <- dim(irf_summary$mean)[3]
  time_horizon <- 1:time_horizon_length
  
  # Determine which IRF to plot (overall or specific period)
  if (period_index == 0) {
    if (summary_type == "mean") {
      central_irf <- irf_summary$overall_mean[var_index, shock_index, ] * scale_factor
    } else if (summary_type == "median") {
      central_irf <- irf_summary$overall_median[var_index, shock_index, ] * scale_factor
    } else {
      stop("Invalid summary_type. Use 'mean' or 'median'.")
    }
    
    lower_bound <- irf_summary$overall_percentiles[[which(percentiles == percentiles[1])]][var_index, shock_index, ] * scale_factor
    upper_bound <- irf_summary$overall_percentiles[[which(percentiles == percentiles[2])]][var_index, shock_index, ] * scale_factor
    
    if (is.null(plot_title)) {
      plot_title <- paste(summary_type, "Overall IRF for Variable", var_index, 
                          "Shock", shock_index)
    }
  } else {
    if (summary_type == "mean") {
      central_irf <- irf_summary$mean[var_index, shock_index, , period_index] * scale_factor
    } else if (summary_type == "median") {
      central_irf <- irf_summary$median[var_index, shock_index, , period_index] * scale_factor
    } else {
      stop("Invalid summary_type. Use 'mean' or 'median'.")
    }
    
    lower_bound <- irf_summary$percentiles[[which(percentiles == percentiles[1])]][var_index, shock_index, , period_index] * scale_factor
    upper_bound <- irf_summary$percentiles[[which(percentiles == percentiles[2])]][var_index, shock_index, , period_index] * scale_factor
    
    if (is.null(plot_title)) {
      plot_title <- paste(summary_type, "IRF for Variable", var_index, 
                          "Shock", shock_index, "at Period", period_index)
    }
  }
  
  # Check for length consistency
  if (length(time_horizon) != length(central_irf) || 
      length(central_irf) != length(lower_bound) || 
      length(lower_bound) != length(upper_bound)) {
    stop("Length mismatch between time horizon and IRF data")
  }
  
  # Create a data frame for plotting
  irf_df <- data.frame(
    Time = time_horizon,
    CentralIRF = central_irf,
    LowerBound = lower_bound,
    UpperBound = upper_bound
  )
  
  # Generate the plot using ggplot2
  p <- ggplot(irf_df, aes(x = Time)) +
    geom_line(aes(y = CentralIRF), color = "#0072B2", size = 1.2) +  # Professional blue color, thicker line
    geom_ribbon(aes(ymin = LowerBound, ymax = UpperBound), fill = "#0072B2", alpha = 0.2) +  # Confidence interval
    labs(title = plot_title, y = y_axis_label, x = x_axis_label) +
    theme_minimal(base_size = 15) +  # Clean, professional theme
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),  # Centered and bold title
      axis.title.x = element_text(face = "bold", size = 12),  # Bold and larger x-axis label
      axis.title.y = element_text(face = "bold", size = 12),  # Bold and larger y-axis label
      axis.text = element_text(size = 10),  # Standardized axis text size
      panel.grid.major = element_line(color = "grey80", size = 0.5),  # Major grid lines
      panel.grid.minor = element_line(color = "grey90", size = 0.25)  # Minor grid lines
    ) +
    scale_y_continuous(labels = scales::percent_format(scale = 1))  # Format y-axis as percentage
  
  return(p)
}

#' Plot IRF by Horizon with Custom X-Axis Label Separation
#'
#' This function generates a time series plot of Impulse Response Functions (IRFs) by a specific horizon from a Time-Varying Parameter VAR (TVP-VAR) model. 
#' It includes customizable options for x-axis labels, date frequency, and scaling of IRF values.
#'
#' @param irf_summary Summary of IRF results, as computed by the `compute_irf_summary_tvp` function.
#' @param var_index Integer, index of the variable for which to plot the IRF (e.g., 1 for GDP).
#' @param shock_index Integer, index of the shock to plot (e.g., 1 for E3CI Shock).
#' @param horizon_index Integer, index of the horizon to plot (e.g., 3 for Horizon 3).
#' @param summary_type String, type of summary statistic to plot: "mean" or "median".
#' @param percentiles Vector of percentiles for the confidence intervals (e.g., c(16, 84)).
#' @param plot_title String, custom title for the plot. If NULL, a default title will be generated.
#' @param x_axis_label String, label for the x-axis. Default is "Calendar Date".
#' @param y_axis_label String, label for the y-axis. Default is "IRF".
#' @param start_date Date string (e.g., "2000-01-01"), the start date for the x-axis.
#' @param freq String, frequency of the periods ("monthly", "quarterly", "yearly"). Default is "monthly".
#' @param x_breaks String, the separation between x-axis labels (e.g., "1 month", "1 quarter", "1 year"). Default is "1 year".
#' @param scale_factor Numeric, factor to scale the IRF values. Default is 100.
#'
#' @return A `ggplot` object representing the IRF plot for the specified horizon.
#'
#' @export

plot_irf_tvp_by_horizon <- function(irf_summary, var_index = 1, shock_index = 1, 
                                    horizon_index = 1, summary_type = "mean", 
                                    percentiles = c(16, 84), plot_title = NULL, 
                                    x_axis_label = "Calendar Date", y_axis_label = "IRF",
                                    start_date = "2000-01-01", freq = "monthly", 
                                    x_breaks = "1 year", scale_factor = 100) {
  
  # Ensure the indices are within range
  if (var_index > dim(irf_summary$mean)[1] || var_index < 1) {
    stop("Variable index out of range")
  }
  if (shock_index > dim(irf_summary$mean)[2] || shock_index < 1) {
    stop("Shock index out of range")
  }
  if (horizon_index > dim(irf_summary$mean)[3] || horizon_index < 1) {
    stop("Horizon index out of range")
  }
  
  # Extract the length of the time periods (number of periods to plot)
  time_period_length <- dim(irf_summary$mean)[4]
  
  # Generate calendar dates for the x-axis based on the start_date and frequency
  if (freq == "monthly") {
    time_periods <- seq.Date(from = as.Date(start_date), by = "month", length.out = time_period_length)
  } else if (freq == "quarterly") {
    time_periods <- seq.Date(from = as.Date(start_date), by = "quarter", length.out = time_period_length)
  } else if (freq == "yearly") {
    time_periods <- seq.Date(from = as.Date(start_date), by = "year", length.out = time_period_length)
  } else {
    stop("Unsupported frequency. Use 'monthly', 'quarterly', or 'yearly'.")
  }
  
  # Extract and scale IRF data based on summary type
  if (summary_type == "mean") {
    central_irf <- irf_summary$mean[var_index, shock_index, horizon_index, ] * scale_factor
  } else if (summary_type == "median") {
    central_irf <- irf_summary$median[var_index, shock_index, horizon_index, ] * scale_factor
  } else {
    stop("Invalid summary_type. Use 'mean' or 'median'.")
  }
  
  # Extract and scale confidence intervals
  lower_bound <- irf_summary$percentiles[[which(percentiles == percentiles[1])]][var_index, shock_index, horizon_index, ] * scale_factor
  upper_bound <- irf_summary$percentiles[[which(percentiles == percentiles[2])]][var_index, shock_index, horizon_index, ] * scale_factor
  
  # Generate default plot title if not provided
  if (is.null(plot_title)) {
    plot_title <- paste(summary_type, "IRF for Variable", var_index, 
                        "Shock", shock_index, "at Horizon", horizon_index)
  }
  
  # Create data frame for ggplot
  irf_df <- data.frame(
    Time = time_periods,
    CentralIRF = central_irf,
    LowerBound = lower_bound,
    UpperBound = upper_bound
  )
  
  # Decide the date label format based on the x_breaks
  if (x_breaks == "1 year") {
    date_format <- "%Y"  # Only show the year for yearly breaks
  } else {
    date_format <- "%Y-%m"  # Show year and month for shorter intervals
  }
  
  # Generate the plot using ggplot2
  p <- ggplot(irf_df, aes(x = Time)) +
    geom_line(aes(y = CentralIRF), color = "#0072B2", size = 1.2) +  # Professional blue color, thicker line
    geom_ribbon(aes(ymin = LowerBound, ymax = UpperBound), fill = "#0072B2", alpha = 0.2) +  # Confidence interval
    labs(title = plot_title, y = y_axis_label, x = x_axis_label) +
    theme_minimal(base_size = 15) +  # Clean, professional theme
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),  # Centered and bold title
      axis.title.x = element_text(face = "bold", size = 12),  # Bold and larger x-axis label
      axis.title.y = element_text(face = "bold", size = 12),  # Bold and larger y-axis label
      axis.text = element_text(size = 10),  # Standardized axis text size
      panel.grid.major = element_line(color = "grey80", size = 0.5),  # Major grid lines
      panel.grid.minor = element_line(color = "grey90", size = 0.25)  # Minor grid lines
    ) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Format y-axis as percentage
    scale_x_date(date_labels = date_format, date_breaks = x_breaks)  # Custom x-axis label breaks and format
  
  return(p)
}

#' Create 3D IRF plot
#'
#'
#' This function generates a 3D plot of Impulse Response Functions (IRFs) from a Time-Varying Parameter VAR (TVP-VAR) model. 
#' It visualizes the IRFs over different horizons and time periods using a 3D surface plot with Plotly.
#'
#' @param irf_summary Summary of IRF results, as computed by the `compute_irf_summary_tvp` function.
#' @param var_index Integer, index of the variable for which to plot the IRF.
#' @param shock_index Integer, index of the shock to plot.
#' @param summary_type String, type of summary statistic to use: "mean" or "median".
#' @param initial_start_year Integer, starting year for the time series. Default is 1983.
#' @param initial_start_quarter Integer, starting quarter for the time series. Default is 1.
#'
#' @return A `plotly` object representing the 3D IRF plot.
#'
#' @export

create_3d_irf_plot <- function(irf_summary, var_index, shock_index, summary_type = "median", 
                               initial_start_year = 1983, initial_start_quarter = 1, 
                               z_axis_label = "Response (%)", plot_title = NULL) {
  
  # Select the appropriate summary data (mean or median)
  if (summary_type == "mean") {
    irf_data <- irf_summary$mean[var_index, shock_index, , ]
  } else if (summary_type == "median") {
    irf_data <- irf_summary$median[var_index, shock_index, , ]
  } else {
    stop("Invalid summary_type. Use 'mean' or 'median'.")
  }
  
  # Scale the data for plotting
  irf_data <- irf_data * 100
  
  # Dimensions of the IRF data
  n_horizons <- dim(irf_data)[1]  # Number of horizons
  n_periods <- dim(irf_data)[2]   # Number of periods
  
  # Generate a sequence of dates corresponding to the periods
  start_year <- initial_start_year
  start_quarter <- initial_start_quarter
  
  dates <- seq(
    as.Date(paste0(start_year, "-", (start_quarter - 1) * 3 + 1, "-01")),
    by = "quarter", length.out = n_periods
  )
  
  # Create matrices for the 3D plot
  x <- matrix(rep(1:n_periods, each = n_horizons), nrow = n_horizons)  # Time Period
  y <- matrix(rep(1:n_horizons, n_periods), nrow = n_horizons)  # Horizon 
  z <- irf_data  # Magnitude of Impulse Response
  
  # Generate the 3D surface plot using plotly
  p <- plot_ly() %>%
    add_surface(
      x = ~x, y = ~y, z = ~z,
      colorscale = list(c(0, "#0072B2"), c(1, "#E69F00")),  # Color scale: blue (low) to orange (high)
      cmin = min(irf_data),
      cmax = max(irf_data),
      colorbar = list(title = z_axis_label)  # Use the provided z-axis label
    )
  
  # Add individual lines for each period (constant date) with a consistent color
  for (i in 1:n_periods) {
    p <- p %>%
      add_trace(
        x = rep(i, n_horizons),
        y = 1:n_horizons,  # y-axis directly matches the horizons
        z = irf_data[, i],
        type = "scatter3d",
        mode = "lines",
        line = list(color = "black", width = 4),  # Consistent blue color for lines
        showlegend = FALSE
      )
  }
  
  # Generate a default plot title if none is provided
  if (is.null(plot_title)) {
    plot_title <- paste("IRF for Variable", var_index, "to Shock", shock_index)
  }
  
  # Layout the plot with customizable labels
  p %>%
    layout(
      scene = list(
        xaxis = list(
          title = "Date", 
          tickvals = seq(1, n_periods, by = 16),  # Tick every 4 years (16 quarters)
          ticktext = format(dates[seq(1, n_periods, by = 16)], "%Y")
        ),
        yaxis = list(
          title = "Horizon period (in quarters)", 
          tickvals = seq(1, n_horizons),
          range = c(n_horizons, 1)  # Reverse y-axis to match z-axis orientation
        ),
        zaxis = list(
          title = z_axis_label  # Use the provided z-axis label
        ),
        camera = list(eye = list(x = -1.5, y = -1.5, z = 0.5))
      ),
      title = plot_title  # Use the provided plot title
    )
}

###################
# Data Preparation
###################

# Load the data from an Excel file
dataprep <-read_excel("C:/_casa/ROBER/MASTER/TFM/R files/BDREMS.xlsx", 
                      sheet = 3, col_names = TRUE)
# Convert to a data frame
data <- as.data.frame(dataprep)
head(data)

# UNIVARIANT ANALYSIS I
# Select necessary columns 
data1 <- data[, c("Fecha", "PIBpm","ppm", "Tasa de paro trimestral", "E3CI")]

# Convert the 'Fecha' column to Date type
data1$Fecha <- as.Date(data1$Fecha)

# Plot series
# Quarterly GDP [M€]
p1 <- create_styled_plot(data1, "Fecha", "PIBpm", 
                         "Quarterly GDP Over Time", 
                         "Gross Domestic Product (Millions €)", 
                         " M€")
# Display plots
print(p1)

# ACF and PACF plots
create_acf_pacf_plot(data1$PIBpm, "ACF and PACF for GDP")

# ADF TEST
adf.test(data1$PIBpm, alternative = "stationary")

# Plot series
# Quarterly GDP Deflator (2015 = 1)
p2 <- create_styled_plot(data1, "Fecha", "ppm", 
                         "Quarterly GDP Deflator Over Time", 
                         "GDP Deflator (2015 = 1)")

# Display plots
print(p2)

# ACF and PACF plots
create_acf_pacf_plot(data1$ppm, "ACF and PACF for GDP Deflator")

# ADF TEST
adf.test(data1$ppm, alternative = "stationary")

# Plot series
# Quarterly Unemployment Rate [%]
p3 <- create_styled_plot(data1, "Fecha", "Tasa de paro trimestral", 
                         "Quarterly Unemployment Rate Over Time", 
                         "Unemployment Rate (%)", 
                         "%") +
  scale_y_continuous(labels = function(x) paste0(comma(x * 100), "%"))

# Display plots
print(p3)

# ACF and PACF plots
create_acf_pacf_plot(data1$`Tasa de paro trimestral`, "ACF and PACF for Unemployment Rate")

#ADF TEST
adf.test(data1$`Tasa de paro trimestral`, alternative = "stationary")

#### UNIVARIANT ANALYSIS II
# Select relevant columns for analysis
data2 <- data[, c("Fecha", "Var anual PIB (cons.)", "Var anual Deflactor del PIB", 
                  "Tasa de paro media anual", "E3CI")]

# Convert the 'Fecha' column to Date type
data2$Fecha <- as.Date(data2$Fecha)

year<- data2$Fecha
gdp<- data2$`Var anual PIB (cons.)`
gdp_deflator<- data2$`Var anual Deflactor del PIB`
unemployment<- data2$`Tasa de paro media anual`
e3ci<- data2$`E3CI`

# PLOTS
# GDP Annual Growth Rate [%]
p4 <- create_styled_plot(data2, "Fecha", "Var anual PIB (cons.)", 
                         "GDP Annual Growth Rate Over Time", 
                         "GDP Growth Rate (%)", 
                         "%") +
  scale_y_continuous(labels = function(x) paste0(comma(x * 100), "%"))

# Display plots
print(p4)

# ACF and PACF plots
create_acf_pacf_plot(data2$`Var anual PIB (cons.)`, "ACF and PACF for GDP Growth Rate")

# ADF TEST
adf.test(gdp, alternative = "stationary")

# PLOTS
# GDP Deflator Annual Growth Rate [%]
p5 <- create_styled_plot(data2, "Fecha", "Var anual Deflactor del PIB", 
                         "GDP Deflator Annual Growth Rate Over Time", 
                         "GDP Deflator Growth Rate (%)", 
                         "%") +
  scale_y_continuous(labels = function(x) paste0(comma(x * 100), "%"))

# Display plots
print(p5)

# ACF and PACF plots
create_acf_pacf_plot(data2$`Var anual Deflactor del PIB`, "ACF and PACF for GDP Deflator Growth Rate")

# ADF test for GDP DEFLATOR
adf.test(gdp_deflator, alternative = "stationary")

# PLOTS
# Annual Average Unemployment Rate [%]
p6 <- create_styled_plot(data2, "Fecha", "Tasa de paro media anual", 
                         "Annual Average Unemployment Rate Over Time", 
                         "Unemployment Rate (%)", 
                         "%") +
  scale_y_continuous(labels = function(x) paste0(comma(x * 100), "%"))

# Display plots
print(p6)

# ACF and PACF plots
create_acf_pacf_plot(data2$`Tasa de paro media anual`, "ACF and PACF for Average Unemployment Rate")

# ADF test for UNEMPLOYMENY 
adf.test(unemployment, alternative = "stationary")

## PLOTS
# Quarterly Average E3CI
p7 <- create_styled_plot(data2, "Fecha", "E3CI", 
                         "Quarterly Average E3CI Over Time", 
                         "E3CI") +
  labs(caption = "Source: E3CI Dataclime")

# Display plots
print(p7)

# ACF and PACF plots
create_acf_pacf_plot(data2$E3CI, "ACF and PACF for Quarterly Average E3CI")      
# ADF test for E3CI 
adf.test(e3ci, alternative = "stationary")

# Create a list of all your plots
plot_list <- list(p1, p2, p3, p4, p5, p6, p7)

# Create a list of names for your plots
plot_names <- c("Quarterly_GDP", "Quarterly_GDP_Deflator", "Quarterly_Unemployment_Rate",
                "GDP_Annual_Growth_Rate", "GDP_Deflator_Annual_Growth_Rate",
                "Annual_Average_Unemployment_Rate", "Quarterly_Average_E3CI")

# Create a directory to save the plots (if it doesn't exist)
dir.create("economic_plots", showWarnings = FALSE)

# Loop through the plots and save them
for (i in 1:length(plot_list)) {
  ggsave(
    filename = paste0("economic_plots/", plot_names[i], ".png"),
    plot = plot_list[[i]],
    width = 10,
    height = 6,
    units = "in",
    dpi = 300
  )
}

print("All plots have been saved in the 'economic_plots' folder.")

## DATA TRANSFORMATIONS TO ACHIEVE STATIONARITY
### 1. First Differencing
# Compute the first difference for relevant time series
diff_gdp <- diff(gdp)
diff_gdp_deflator <- diff(gdp_deflator)
diff_unemployment <- diff(unemployment)

### 2. Stationarity Tests
# Perform Augmented Dickey-Fuller (ADF) tests on differenced series
adf.test(diff_gdp, alternative = "stationary")
adf.test(diff_gdp_deflator, alternative = "stationary")
adf.test(diff_unemployment, alternative = "stationary")

### 3. Seasonal Differencing
# Perform seasonal differencing on the unemployment rate
seasonal_diff_unemployment <- diff(diff_unemployment, lag = 4)

# ADF test for the seasonally differenced series
adf.test(seasonal_diff_unemployment, alternative = "stationary")

# Prepare the final dataset for VAR modeling
data3 <- data2[-(1:5), ]  # Remove the first 5 rows to align with seasonal_diff_paro
data3$`Var anual PIB (cons.)` <- diff_gdp[-(1:4)]
data3$`Var anual Deflactor del PIB` <- diff_gdp_deflator[-(1:4)]
data3$`Tasa de paro media anual` <- seasonal_diff_unemployment
data3$E3CI <- data2$E3CI[-(1:5)]

## Vector Autoregression (VAR) Modeling
# Define the start date for the time series
start_date3 <- c(year(min(year)) + 1, quarter(min(year)) + 1)

# Convert the data to a time series object
data3_ts <- ts(data3[, -which(names(data3) == "Fecha")],  # Exclude the 'Fecha' column if it exists
               start = start_date3,
               frequency = 4)  # Quarterly data

# Prepare the data for plotting
data3$Fecha <- seq(as.Date(paste0(year(min(data3$Fecha)), "-01-01")), 
                   by = "quarter", 
                   length.out = nrow(data3))

# Create a list to store all plots
plot_list <- list()

# Create plots for each variable with improved labels
for (col in names(data3)[names(data3) != "Fecha"]) {
  # Create more descriptive titles and labels
  if (col == "Var anual PIB (cons.)") {
    plot_title <- "GDP  Over Time"
    y_label <- "GDP (%)"
  } else if (col == "Var anual Deflactor del PIB") {
    plot_title <- "GDP Deflator Over Time"
    y_label <- "GDP Deflator (%)"
  } else if (col == "Tasa de paro media anual") {
    plot_title <- "Unemployment Rate Over Time"
    y_label <- "Unemployment Rate (%)"
  } else if (col == "E3CI") {
    plot_title <- "Transformed E3CI"
    y_label <- "Transformed E3CI"
  } else {
    plot_title <- paste("Transformed", col, "Over Time")
    y_label <- paste("Transformed", col)
  }
  
  plot_list[[col]] <- create_styled_plot(data3, "Fecha", col, plot_title, y_label, "%")
}

# Display all plots
for (p in plot_list) {
  print(p)
}

# Save all plots
dir.create("transformed_plots", showWarnings = FALSE)

for (i in seq_along(plot_list)) {
  ggsave(
    filename = paste0("transformed_plots/Transformed_", names(plot_list)[i], ".png"),
    plot = plot_list[[i]],
    width = 10,
    height = 6,
    units = "in",
    dpi = 300
  )
}

print("All transformed data plots have been saved in the 'transformed_plots' folder.")

#####################################
# Initial call to VAR and SVAR models
#####################################

# Combine the transformed series into a data frame
transformed_data <- data3[, -which(names(data3) == "Fecha")]

# Define the number of lags
p = 2 # Reference paper and BdE model

# Estimate the VAR model
var_model <- VAR(transformed_data, p, type = "const")

# Extract the covariance matrix of residuals (Sigma)
sigma <- summary(var_model)$covres

# Cholesky decomposition of the covariance matrix of residuals
P_chol <- t(chol(sigma))

# Define a B matrix with NA values for SVAR model initiliazation (on top of Cholesky)
Bmat <- P_chol
Bmat[upper.tri(Bmat)] <- NA
for (i in 1:4) {Bmat[i,i] <- NA}

# Identify the B matrix of the SVAR model, and set it back to Cholesky
svar_model <- SVAR(var_model, Bmat = Bmat, max.iter =1000)
svar_model$B <- P_chol

# ####################
# Scenario definition
# ####################

# Set seed for reproducibility
set.seed(123)

# Define the number of trials
n_Q <- 50000  # Number of Q matrices to generate

# Define sign restrictions matrix
sign_restrictions <- matrix(c(
  NA,  NA,  NA, -1,  # Variable 1: Positive sign expected for Shock 1, Negative for Shock 2, No restriction for Shock 3, Positive for Shock 4
  NA,  NA,  NA,  1,  # Variable 2: Negative sign expected for Shock 1, Positive for Shock 2, No restriction for Shock 3, Negative for Shock 4
  NA,  NA,  NA,  1,  # Variable 3: No sign restrictions for Shock 1 and 2, Positive for Shock 3, No restriction for Shock 4
  NA,  NA,  NA,  NA   # Variable 4: Positive sign expected for Shock 1, Negative for Shock 2, No restriction for Shock 3, Positive for Shock 4
), nrow = 4, byrow = TRUE)

# Set the number of periods to check in the IRF (1 = only contemporaneous effects)
periods_to_check <- 1 # You can adjust this value as needed >= 1 # Reference paper and BdE model = 1

# Set the number of periods in the IRF horizon (number of periods = horizon + 1)
horizon <- 12 # Ensure that horizon >= periods_to_check -1 # Reference paper and BdE model = 12

# Set the sign restrictions to contemparenous/short term (FALSE) or long-term (TRUE)
long_term = FALSE

if (long_term) {
  horizon_to_check <- 0
  periods_to_check_lt <- 1
} else {
  horizon_to_check <- periods_to_check
  periods_to_check_lt <- periods_to_check
}

# Set the percentiles for the upper and lower bounds in the plot
percentiles <- c(5, 95)  ## en la funci?n para el plot esta definido 14 y 86

# #########################
# Variables initialization
# #########################

valid_B_matrices <- list()
valid_irf_list <- list()

# #####################
# Scenario computation
# #####################

for (i in 1:n_Q) {
  Q <- generate_Q(n = ncol(var_model$y))
  B_candidate <- P_chol %*% Q
  
  irf_result_svar <- compute_irfs_svar(svar_model, B_candidate, horizon_to_check)
  
  if (check_sign_restrictions_svar(irf_result_svar, sign_restrictions, periods_to_check_lt)) {
    valid_B_matrices <- append(valid_B_matrices, list(B_candidate))
    irf_result_svar <- compute_irfs_svar(svar_model, B_candidate, horizon)
    valid_irf_list <- append(valid_irf_list, list(irf_result_svar$irf))
  }
}

cat("Number of valid B matrices:", length(valid_B_matrices), "\n")

#####################################################################
# Scenario statistics and plots preparation for the VAR and SVAR model
#####################################################################

# Calculate the IRFs of the VAR model
irf_var <- vars::irf(var_model, n.ahead = horizon, ortho = TRUE, boot = FALSE)$irf

# Compute the summary statistics for the valid IRFs
irf_summary_svar <- compute_irf_summary_svar(valid_irf_list, percentiles)

# Set the shock index and plot title and x-axis label
shock_index <- 4
plot_title <- "IRF for a shock on the E3CI variable using a VAR model or a SVAR model with sign restrictions"
x_axis_label <- "Horizon period (in quarters)"
summary_type <- "median"

#######################################
# Generate the plots for the 3 variables
#######################################

# Plot the IRFs with upper and lower bounds

time_horizon <- 1:(horizon + 1)

# Plot for the 1st variable
var_index <- 1
overlay_irf <- irf_var[[shock_index]][, var_index] * 100  # Scale overlay IRF
y_axis_label <- "GDP Annual Growth Rate (%)"
base_plot <- plot_irf_svar(irf_summary_svar, var_index, shock_index, summary_type, percentiles, plot_title, x_axis_label, y_axis_label, y_min=NULL, y_max=NULL)
final_plot <- base_plot +
  geom_line(aes(x = time_horizon, y = overlay_irf), color = "red", size = 1.5, linetype = "dashed") +  # Make VAR IRF more prominent
  labs(subtitle = "Dashed line represents VAR model IRF") +  # Add subtitle for clarity
  theme(legend.position = "none")  # Remove legend since it's not necessary

print(final_plot)

# Plot for the 2nd variable
var_index <- 2
overlay_irf <- irf_var[[shock_index]][, var_index] * 100  # Scale overlay IRF
y_axis_label <- "GDP Deflator Annual Growth Rate (%)"
base_plot <- plot_irf_svar(irf_summary_svar, var_index, shock_index, summary_type, percentiles, plot_title, x_axis_label, y_axis_label)
final_plot <- base_plot +
  geom_line(aes(x = time_horizon, y = overlay_irf), color = "red", size = 1.5, linetype = "dashed") +  # Make VAR IRF more prominent
  labs(subtitle = "Dashed line represents VAR model IRF") +
  theme(legend.position = "none")

print(final_plot)

# Plot for the 3rd variable
var_index <- 3
overlay_irf <- irf_var[[shock_index]][, var_index] * 100  # Scale overlay IRF
y_axis_label <- "Annual Average Unemployment Rate (%)"
base_plot <- plot_irf_svar(irf_summary_svar, var_index, shock_index, summary_type, percentiles, plot_title, x_axis_label, y_axis_label)
final_plot <- base_plot +
  geom_line(aes(x = time_horizon, y = overlay_irf), color = "red", size = 1.5, linetype = "dashed") +  # Make VAR IRF more prominent
  labs(subtitle = "Dashed line represents VAR model IRF") +
  theme(legend.position = "none")

print(final_plot)

####################################
# TVP BVAR without sign restrictions
####################################

# Parameters
thinfac <- 1
tau <- 40

n_hor<- horizon

nburn <- 5000
nrep <- 20000

t_comp <- 1
variable <- 1

use_intercept <- FALSE

# Model 
fit <- bvar.sv.tvp(data3_ts, p = p , nburn = nburn, nrep = nrep, thinfac = thinfac)

# Parameters Initialization 
Beta_draws <- fit$Beta.draws[,,]
H_draws <- fit$H.draws[,,]
dim_H <- dim(H_draws)
B_draws <- array(0, dim = dim_H)

n_periods <- nrow(data3_ts) - tau - p

for (t in 1:n_periods) {
  for (i in 1:nrep) {
    start_col <- (t - 1) * 4 + 1  # Columna inicial del bloque
    end_col <- t * 4              # Columna final del bloque
    B_draws[, start_col:end_col, i] <-  t(chol(H_draws[, start_col:end_col, i]))
  }
}

#################################
# TVP BVAR with sign restrictions
#################################

max_trials<- 10000

# Nested loops to update H_draws with sign-restricted Sigma
for (t in 1:n_periods) {
  for (i in 1:nrep) {
    # Calcula las columnas que corresponden al instante de tiempo t
    start_col <- (t - 1) * 4 + 1  # Columna inicial del bloque
    end_col <- t * 4              # Columna final del bloque
    
    # Extrae la matriz Sigma del bloque correspondiente
    Sigma <- H_draws[, start_col:end_col, i]
    
    # Generate candidate B_candidate and check sign restrictions
    candidate<- FALSE
    for (n in 1:max_trials) {
      Q <- generate_Q(4)
      B_candidate <- t(chol(Sigma)) %*% Q
      
      if (check_sign_restrictions_tvp(B_candidate, sign_restrictions)) {
        # Calcula la nueva Sigma y actualiza el bloque en H_draws
        B_draws[, start_col:end_col, i] <-  B_candidate
        
        candidate<- TRUE
        break
      }
    }
    if (candidate == FALSE) {
      print( "NO B valid matrix found")}
  }
}

time_irf<- n_periods- n_hor # para calcualar irf en instantes de tiempo en los que tenemos datos

irf_results_tvp<-calculate_irf_tvp(n_variables=4, p= p, n_periods= time_irf, n_reps= nrep, Beta_draws= Beta_draws, B_draws= B_draws, horizon= n_hor, use_intercept = use_intercept)

irf_summary_tvp <- compute_irf_summary_tvp(irf_results_tvp, percentiles = percentiles)

######################################### 
#IRF PLOT FOR A SPECIPIC PERIOD OF TIME
#########################################

plot_irf_tvp(
  irf_summary = irf_summary_tvp, 
  var_index = 1, 
  shock_index = 4, 
  period_index = 0, 
  summary_type = "median", 
  percentiles = percentiles,
  plot_title = " IRF for a shock on the E3CI variable using a Bayesian TVP-VAR model with sign restrictions", 
  x_axis_label = "Horizon period (in quarters)", 
  y_axis_label = "GDP Annual Growth Rate (%)",
  scale_factor = 100
)

plot_irf_tvp(
  irf_summary = irf_summary_tvp, 
  var_index = 2, 
  shock_index = 4, 
  period_index = 0, 
  summary_type = "median", 
  percentiles = percentiles,
  plot_title = "IRF for a shock on the E3CI variable using a Bayesian TVP-VAR model with sign restrictions", 
  x_axis_label = "Horizon period (in quarters)", 
  y_axis_label = "GDP Deflator Annual Growth Rate (%)",
  scale_factor = 100
)

plot_irf_tvp(
  irf_summary = irf_summary_tvp, 
  var_index = 3, 
  shock_index = 4, 
  period_index = 0, 
  summary_type = "median", 
  percentiles = percentiles,
  plot_title = "IRF for a shock on the E3CI variable using a Bayesian TVP-VAR model with sign restrictions", 
  x_axis_label = "Horizon period (in quarters)", 
  y_axis_label = "Annual Average Unemployment Rate (%)",
  scale_factor = 100
)

######################################### 
#IRF PLOT FOR A SPECIFIC HORIZON ANALYSIS
#########################################

percentiles <- c(32, 68)  ## en la funci?n para el plot esta definido 14 y 86

irf_summary_tvp <- compute_irf_summary_tvp(irf_results_tvp, percentiles = percentiles)

plot_irf_tvp_by_horizon(
  irf_summary = irf_summary_tvp,  # Replace with your actual IRF summary object
  var_index = 1,                  # Index of the variable (e.g., 1 for GDP)
  shock_index = 4,                # Index of the shock (e.g., 4 for E3CI Shock)
  horizon_index = 3,              # Horizon index (e.g., 3 for Horizon 3)
  summary_type = "median",        # Summary type ('mean' or 'median')
  percentiles = percentiles,      # Percentiles for confidence intervals
  plot_title = "IRF for a shock on the E3CI variable using a Bayesian TVP-VAR model with sign restrictions",  # Custom plot title
  x_axis_label = "Date",   # Label for x-axis
  y_axis_label = "GDP Annual Growth Rate (%)",           # Label for y-axis
  start_date = "1993-01-01",
  freq = "quarterly",
  x_breaks = "1 year",
  scale_factor = 100              # Scaling factor
)

plot_irf_tvp_by_horizon(
  irf_summary = irf_summary_tvp,  # Replace with your actual IRF summary object
  var_index = 2,                  # Index of the variable (e.g., 1 for GDP)
  shock_index = 4,                # Index of the shock (e.g., 4 for E3CI Shock)
  horizon_index = 3,              # Horizon index (e.g., 3 for Horizon 3)
  summary_type = "median",        # Summary type ('mean' or 'median')
  percentiles = percentiles,      # Percentiles for confidence intervals
  plot_title = "IRF for a shock on the E3CI variable using a Bayesian TVP-VAR model with sign restrictions",  # Custom plot title
  x_axis_label = "Date",   # Label for x-axis
  y_axis_label = "GDP Deflator Annual Growth Rate (%)" ,          # Label for y-axis
  start_date = "1993-01-01",
  freq = "quarterly",
  x_breaks = "1 year",
  scale_factor = 100  
)

plot_irf_tvp_by_horizon(
  irf_summary = irf_summary_tvp,  # Replace with your actual IRF summary object
  var_index = 3,                  # Index of the variable (e.g., 1 for GDP)
  shock_index = 4,                # Index of the shock (e.g., 4 for E3CI Shock)
  horizon_index = 3,              # Horizon index (e.g., 3 for Horizon 3)
  summary_type = "median",        # Summary type ('mean' or 'median')
  percentiles = percentiles,      # Percentiles for confidence intervals
  plot_title = "IRF for a shock on the E3CI variable using a Bayesian TVP-VAR model with sign restrictions",  # Custom plot title
  x_axis_label = "Date",   # Label for x-axis
  y_axis_label = "Annual Average Unemployment Rate (%)"  ,         # Label for y-axis
  start_date = "1993-01-01",
  freq = "quarterly",
  x_breaks = "1 year",
  scale_factor = 100  
)

####################### 
#3D PLOTS
#######################

# Plot IRF for GDP
plot_gdp <- create_3d_irf_plot(irf_summary_tvp, var_index = 1, shock_index = 4, summary_type = "median", initial_start_year = 1993, initial_start_quarter = 1, z_axis_label = "GDP Annual Growth Rate (%)", plot_title = "IRF for a shock on the E3CI variable using a Bayesian TVP-VAR model with sign restrictions")
plot_gdp

# Plot IRF for GDP Deflator
plot_deflator <- create_3d_irf_plot(irf_summary_tvp, var_index = 2, shock_index = 4, summary_type = "median", initial_start_year = 1993, initial_start_quarter = 1, z_axis_label = "GDP Deflator Annual Growth Rate (%)", plot_title = "IRF for a shock on the E3CI variable using a Bayesian TVP-VAR model with sign restrictions")
plot_deflator

# Plot IRF for Unemployment
plot_unemployment <- create_3d_irf_plot(irf_summary_tvp, var_index = 3, shock_index = 4, summary_type = "median", initial_start_year = 1993, initial_start_quarter = 1, z_axis_label = "Annual Average Unemployment Rate (%)", plot_title = "IRF for a shock on the E3CI variable using a Bayesian TVP-VAR model with sign restrictions")
plot_unemployment

##############################################
#ADDITIONAL ANALYSIS
##############################################

# #####################
# NUMBER OF LAGS
# #####################

# Define the number of lags
p = 6 # Reference paper and BdE model

# Estimate the VAR model
var_model <- VAR(transformed_data, p, type = "const")

# Extract the covariance matrix of residuals (Sigma)
sigma <- summary(var_model)$covres

# Cholesky decomposition of the covariance matrix of residuals
P_chol <- t(chol(sigma))

# Define a B matrix with NA values for SVAR model initiliazation (on top of Cholesky)
Bmat <- P_chol
Bmat[upper.tri(Bmat)] <- NA
for (i in 1:4) {Bmat[i,i] <- NA}

# Identify the B matrix of the SVAR model, and set it back to Cholesky
svar_model <- SVAR(var_model, Bmat = Bmat, max.iter =1000)
svar_model$B <- P_chol

# Set seed for reproducibility
set.seed(123)

# Define the number of trials
n_Q <- 50000  # Number of Q matrices to generate

# Define sign restrictions matrix
sign_restrictions <- matrix(c(
  NA,  NA,  NA, -1,  # Variable 1: Positive sign expected for Shock 1, Negative for Shock 2, No restriction for Shock 3, Positive for Shock 4
  NA,  NA,  NA,  1,  # Variable 2: Negative sign expected for Shock 1, Positive for Shock 2, No restriction for Shock 3, Negative for Shock 4
  NA,  NA,  NA,  1,  # Variable 3: No sign restrictions for Shock 1 and 2, Positive for Shock 3, No restriction for Shock 4
  NA,  NA,  NA,  NA   # Variable 4: Positive sign expected for Shock 1, Negative for Shock 2, No restriction for Shock 3, Positive for Shock 4
), nrow = 4, byrow = TRUE)

# Set the number of periods to check in the IRF (1 = only contemporaneous effects)
periods_to_check <- 1 # You can adjust this value as needed >= 1 # Reference paper and BdE model = 1

# Set the number of periods in the IRF horizon (number of periods = horizon + 1)
horizon <- 12 # Ensure that horizon >= periods_to_check -1 # Reference paper and BdE model = 12

# Set the sign restrictions to contemparenous/short term (FALSE) or long-term (TRUE)
long_term = FALSE

if (long_term) {
  horizon_to_check <- 0
  periods_to_check_lt <- 1
} else {
  horizon_to_check <- periods_to_check
  periods_to_check_lt <- periods_to_check
}

# Set the percentiles for the upper and lower bounds in the plot
percentiles <- c(5, 95)  ## en la funci?n para el plot esta definido 14 y 86

valid_B_matrices <- list()
valid_irf_list <- list()

for (i in 1:n_Q) {
  Q <- generate_Q(n = ncol(var_model$y))
  B_candidate <- P_chol %*% Q
  
  irf_result_svar <- compute_irfs_svar(svar_model, B_candidate, horizon_to_check)
  
  if (check_sign_restrictions_svar(irf_result_svar, sign_restrictions, periods_to_check_lt)) {
    valid_B_matrices <- append(valid_B_matrices, list(B_candidate))
    irf_result_svar <- compute_irfs_svar(svar_model, B_candidate, horizon)
    valid_irf_list <- append(valid_irf_list, list(irf_result_svar$irf))
  }
}

cat("Number of valid B matrices:", length(valid_B_matrices), "\n")

# Calculate the IRFs of the VAR model
irf_var <- vars::irf(var_model, n.ahead = horizon, ortho = TRUE, boot = FALSE)$irf

# Compute the summary statistics for the valid IRFs
irf_summary_svar <- compute_irf_summary_svar(valid_irf_list, percentiles)

# Set the shock index and plot title and x-axis label
shock_index <- 4
plot_title <- "IRF for a shock on the E3CI variable using a VAR model or a SVAR model with sign restrictions"
x_axis_label <- "Horizon period (in quarters)"
summary_type <- "median"

# Plot the IRFs with upper and lower bounds

time_horizon <- 1:(horizon + 1)

# Plot for the 1st variable
var_index <- 1
overlay_irf <- irf_var[[shock_index]][, var_index] * 100  # Scale overlay IRF
y_axis_label <- "GDP Annual Growth Rate (%)"
base_plot <- plot_irf_svar(irf_summary_svar, var_index, shock_index, summary_type, percentiles, plot_title, x_axis_label, y_axis_label, y_min=NULL, y_max=NULL)
final_plot <- base_plot +
  geom_line(aes(x = time_horizon, y = overlay_irf), color = "red", size = 1.5, linetype = "dashed") +  # Make VAR IRF more prominent
  labs(subtitle = "Dashed line represents VAR model IRF") +  # Add subtitle for clarity
  theme(legend.position = "none")  # Remove legend since it's not necessary

print(final_plot)

# Define the number of lags
p = 9 # Reference paper and BdE model

# Estimate the VAR model
var_model <- VAR(transformed_data, p, type = "const")

# Extract the covariance matrix of residuals (Sigma)
sigma <- summary(var_model)$covres

# Cholesky decomposition of the covariance matrix of residuals
P_chol <- t(chol(sigma))

# Define a B matrix with NA values for SVAR model initiliazation (on top of Cholesky)
Bmat <- P_chol
Bmat[upper.tri(Bmat)] <- NA
for (i in 1:4) {Bmat[i,i] <- NA}

# Identify the B matrix of the SVAR model, and set it back to Cholesky
svar_model <- SVAR(var_model, Bmat = Bmat, max.iter =1000)
svar_model$B <- P_chol

# Set seed for reproducibility
set.seed(123)

# Define the number of trials
n_Q <- 50000  # Number of Q matrices to generate

# Define sign restrictions matrix
sign_restrictions <- matrix(c(
  NA,  NA,  NA, -1,  # Variable 1: Positive sign expected for Shock 1, Negative for Shock 2, No restriction for Shock 3, Positive for Shock 4
  NA,  NA,  NA,  1,  # Variable 2: Negative sign expected for Shock 1, Positive for Shock 2, No restriction for Shock 3, Negative for Shock 4
  NA,  NA,  NA,  1,  # Variable 3: No sign restrictions for Shock 1 and 2, Positive for Shock 3, No restriction for Shock 4
  NA,  NA,  NA,  NA   # Variable 4: Positive sign expected for Shock 1, Negative for Shock 2, No restriction for Shock 3, Positive for Shock 4
), nrow = 4, byrow = TRUE)

# Set the number of periods to check in the IRF (1 = only contemporaneous effects)
periods_to_check <- 1 # You can adjust this value as needed >= 1 # Reference paper and BdE model = 1

# Set the number of periods in the IRF horizon (number of periods = horizon + 1)
horizon <- 12 # Ensure that horizon >= periods_to_check -1 # Reference paper and BdE model = 12

# Set the sign restrictions to contemparenous/short term (FALSE) or long-term (TRUE)
long_term = FALSE

if (long_term) {
  horizon_to_check <- 0
  periods_to_check_lt <- 1
} else {
  horizon_to_check <- periods_to_check
  periods_to_check_lt <- periods_to_check
}

# Set the percentiles for the upper and lower bounds in the plot
percentiles <- c(5, 95)  ## en la funci?n para el plot esta definido 14 y 86

valid_B_matrices <- list()
valid_irf_list <- list()

for (i in 1:n_Q) {
  Q <- generate_Q(n = ncol(var_model$y))
  B_candidate <- P_chol %*% Q
  
  irf_result_svar <- compute_irfs_svar(svar_model, B_candidate, horizon_to_check)
  
  if (check_sign_restrictions_svar(irf_result_svar, sign_restrictions, periods_to_check_lt)) {
    valid_B_matrices <- append(valid_B_matrices, list(B_candidate))
    irf_result_svar <- compute_irfs_svar(svar_model, B_candidate, horizon)
    valid_irf_list <- append(valid_irf_list, list(irf_result_svar$irf))
  }
}

cat("Number of valid B matrices:", length(valid_B_matrices), "\n")

# Calculate the IRFs of the VAR model
irf_var <- vars::irf(var_model, n.ahead = horizon, ortho = TRUE, boot = FALSE)$irf

# Compute the summary statistics for the valid IRFs
irf_summary_svar <- compute_irf_summary_svar(valid_irf_list, percentiles)

# Set the shock index and plot title and x-axis label
shock_index <- 4
plot_title <- "IRF for a shock on the E3CI variable using a VAR model or a SVAR model with sign restrictions"
x_axis_label <- "Horizon period (in quarters)"
summary_type <- "median"

# Plot the IRFs with upper and lower bounds

time_horizon <- 1:(horizon + 1)

# Plot for the 1st variable
var_index <- 1
overlay_irf <- irf_var[[shock_index]][, var_index] * 100  # Scale overlay IRF
y_axis_label <- "GDP Annual Growth Rate (%)"
base_plot <- plot_irf_svar(irf_summary_svar, var_index, shock_index, summary_type, percentiles, plot_title, x_axis_label, y_axis_label, y_min=NULL, y_max=NULL)
final_plot <- base_plot +
  geom_line(aes(x = time_horizon, y = overlay_irf), color = "red", size = 1.5, linetype = "dashed") +  # Make VAR IRF more prominent
  labs(subtitle = "Dashed line represents VAR model IRF") +  # Add subtitle for clarity
  theme(legend.position = "none")  # Remove legend since it's not necessary

print(final_plot)

# #####################
# TYPE AND HORIZON OF SIGN RESTRICTIONS
# #####################

# Define the number of lags
p = 2 # Reference paper and BdE model

# Estimate the VAR model
var_model <- VAR(transformed_data, p, type = "const")

# Extract the covariance matrix of residuals (Sigma)
sigma <- summary(var_model)$covres

# Cholesky decomposition of the covariance matrix of residuals
P_chol <- t(chol(sigma))

# Define a B matrix with NA values for SVAR model initiliazation (on top of Cholesky)
Bmat <- P_chol
Bmat[upper.tri(Bmat)] <- NA
for (i in 1:4) {Bmat[i,i] <- NA}

# Identify the B matrix of the SVAR model, and set it back to Cholesky
svar_model <- SVAR(var_model, Bmat = Bmat, max.iter =1000)
svar_model$B <- P_chol

# Set seed for reproducibility
set.seed(123)

# Define the number of trials
n_Q <- 50000  # Number of Q matrices to generate

# Define sign restrictions matrix
sign_restrictions <- matrix(c(
  NA,  NA,  NA, -1,  # Variable 1: Positive sign expected for Shock 1, Negative for Shock 2, No restriction for Shock 3, Positive for Shock 4
  NA,  NA,  NA,  1,  # Variable 2: Negative sign expected for Shock 1, Positive for Shock 2, No restriction for Shock 3, Negative for Shock 4
  NA,  NA,  NA,  1,  # Variable 3: No sign restrictions for Shock 1 and 2, Positive for Shock 3, No restriction for Shock 4
  NA,  NA,  NA,  NA   # Variable 4: Positive sign expected for Shock 1, Negative for Shock 2, No restriction for Shock 3, Positive for Shock 4
), nrow = 4, byrow = TRUE)

# Set the number of periods to check in the IRF (1 = only contemporaneous effects)
periods_to_check <- 1 # You can adjust this value as needed >= 1 # Reference paper and BdE model = 1

# Set the number of periods in the IRF horizon (number of periods = horizon + 1)
horizon <- 20 # Ensure that horizon >= periods_to_check -1 # Reference paper and BdE model = 12

# Set the sign restrictions to contemparenous/short term (FALSE) or long-term (TRUE)
long_term = FALSE

if (long_term) {
  horizon_to_check <- 0
  periods_to_check_lt <- 1
} else {
  horizon_to_check <- periods_to_check
  periods_to_check_lt <- periods_to_check
}

# Set the percentiles for the upper and lower bounds in the plot
percentiles <- c(5, 95)  ## en la funci?n para el plot esta definido 14 y 86

valid_B_matrices <- list()
valid_irf_list <- list()

for (i in 1:n_Q) {
  Q <- generate_Q(n = ncol(var_model$y))
  B_candidate <- P_chol %*% Q
  
  irf_result_svar <- compute_irfs_svar(svar_model, B_candidate, horizon_to_check)
  
  if (check_sign_restrictions_svar(irf_result_svar, sign_restrictions, periods_to_check_lt)) {
    valid_B_matrices <- append(valid_B_matrices, list(B_candidate))
    irf_result_svar <- compute_irfs_svar(svar_model, B_candidate, horizon)
    valid_irf_list <- append(valid_irf_list, list(irf_result_svar$irf))
  }
}

cat("Number of valid B matrices:", length(valid_B_matrices), "\n")

# Calculate the IRFs of the VAR model
irf_var <- vars::irf(var_model, n.ahead = horizon, ortho = TRUE, boot = FALSE)$irf

# Compute the summary statistics for the valid IRFs
irf_summary_svar <- compute_irf_summary_svar(valid_irf_list, percentiles)

# Set the shock index and plot title and x-axis label
shock_index <- 4
plot_title <- "IRF for a shock on the E3CI variable using a VAR model or a SVAR model with sign restrictions"
x_axis_label <- "Horizon period (in quarters)"
summary_type <- "median"

# Plot the IRFs with upper and lower bounds

time_horizon <- 1:(horizon + 1)

# Plot for the 1st variable
var_index <- 1
overlay_irf <- irf_var[[shock_index]][, var_index] * 100  # Scale overlay IRF
y_axis_label <- "GDP Annual Growth Rate (%)"
base_plot <- plot_irf_svar(irf_summary_svar, var_index, shock_index, summary_type, percentiles, plot_title, x_axis_label, y_axis_label, y_min=NULL, y_max=NULL)
final_plot <- base_plot +
  geom_line(aes(x = time_horizon, y = overlay_irf), color = "red", size = 1.5, linetype = "dashed") +  # Make VAR IRF more prominent
  labs(subtitle = "Dashed line represents VAR model IRF") +  # Add subtitle for clarity
  theme(legend.position = "none")  # Remove legend since it's not necessary

print(final_plot)

# Define the number of lags
p = 2 # Reference paper and BdE model

# Estimate the VAR model
var_model <- VAR(transformed_data, p, type = "const")

# Extract the covariance matrix of residuals (Sigma)
sigma <- summary(var_model)$covres

# Cholesky decomposition of the covariance matrix of residuals
P_chol <- t(chol(sigma))

# Define a B matrix with NA values for SVAR model initiliazation (on top of Cholesky)
Bmat <- P_chol
Bmat[upper.tri(Bmat)] <- NA
for (i in 1:4) {Bmat[i,i] <- NA}

# Identify the B matrix of the SVAR model, and set it back to Cholesky
svar_model <- SVAR(var_model, Bmat = Bmat, max.iter =1000)
svar_model$B <- P_chol

# Set seed for reproducibility
set.seed(123)

# Define the number of trials
n_Q <- 50000  # Number of Q matrices to generate

# Define sign restrictions matrix
sign_restrictions <- matrix(c(
  NA,  NA,  NA, -1,  # Variable 1: Positive sign expected for Shock 1, Negative for Shock 2, No restriction for Shock 3, Positive for Shock 4
  NA,  NA,  NA,  1,  # Variable 2: Negative sign expected for Shock 1, Positive for Shock 2, No restriction for Shock 3, Negative for Shock 4
  NA,  NA,  NA,  1,  # Variable 3: No sign restrictions for Shock 1 and 2, Positive for Shock 3, No restriction for Shock 4
  NA,  NA,  NA,  NA   # Variable 4: Positive sign expected for Shock 1, Negative for Shock 2, No restriction for Shock 3, Positive for Shock 4
), nrow = 4, byrow = TRUE)

# Set the number of periods to check in the IRF (1 = only contemporaneous effects)
periods_to_check <- 8 # You can adjust this value as needed >= 1 # Reference paper and BdE model = 1

# Set the number of periods in the IRF horizon (number of periods = horizon + 1)
horizon <- 20 # Ensure that horizon >= periods_to_check -1 # Reference paper and BdE model = 12

# Set the sign restrictions to contemparenous/short term (FALSE) or long-term (TRUE)
long_term = FALSE

if (long_term) {
  horizon_to_check <- 0
  periods_to_check_lt <- 1
} else {
  horizon_to_check <- periods_to_check
  periods_to_check_lt <- periods_to_check
}

# Set the percentiles for the upper and lower bounds in the plot
percentiles <- c(5, 95)  ## en la funci?n para el plot esta definido 14 y 86

valid_B_matrices <- list()
valid_irf_list <- list()

for (i in 1:n_Q) {
  Q <- generate_Q(n = ncol(var_model$y))
  B_candidate <- P_chol %*% Q
  
  irf_result_svar <- compute_irfs_svar(svar_model, B_candidate, horizon_to_check)
  
  if (check_sign_restrictions_svar(irf_result_svar, sign_restrictions, periods_to_check_lt)) {
    valid_B_matrices <- append(valid_B_matrices, list(B_candidate))
    irf_result_svar <- compute_irfs_svar(svar_model, B_candidate, horizon)
    valid_irf_list <- append(valid_irf_list, list(irf_result_svar$irf))
  }
}

cat("Number of valid B matrices:", length(valid_B_matrices), "\n")

# Calculate the IRFs of the VAR model
irf_var <- vars::irf(var_model, n.ahead = horizon, ortho = TRUE, boot = FALSE)$irf

# Compute the summary statistics for the valid IRFs
irf_summary_svar <- compute_irf_summary_svar(valid_irf_list, percentiles)

# Set the shock index and plot title and x-axis label
shock_index <- 4
plot_title <- "IRF for a shock on the E3CI variable using a VAR model or a SVAR model with sign restrictions"
x_axis_label <- "Horizon period (in quarters)"
summary_type <- "median"

# Plot the IRFs with upper and lower bounds

time_horizon <- 1:(horizon + 1)

# Plot for the 1st variable
var_index <- 1
overlay_irf <- irf_var[[shock_index]][, var_index] * 100  # Scale overlay IRF
y_axis_label <- "GDP Annual Growth Rate (%)"
base_plot <- plot_irf_svar(irf_summary_svar, var_index, shock_index, summary_type, percentiles, plot_title, x_axis_label, y_axis_label, y_min=NULL, y_max=NULL)
final_plot <- base_plot +
  geom_line(aes(x = time_horizon, y = overlay_irf), color = "red", size = 1.5, linetype = "dashed") +  # Make VAR IRF more prominent
  labs(subtitle = "Dashed line represents VAR model IRF") +  # Add subtitle for clarity
  theme(legend.position = "none")  # Remove legend since it's not necessary

print(final_plot)

# Define the number of lags
p = 2 # Reference paper and BdE model

# Estimate the VAR model
var_model <- VAR(transformed_data, p, type = "const")

# Extract the covariance matrix of residuals (Sigma)
sigma <- summary(var_model)$covres

# Cholesky decomposition of the covariance matrix of residuals
P_chol <- t(chol(sigma))

# Define a B matrix with NA values for SVAR model initiliazation (on top of Cholesky)
Bmat <- P_chol
Bmat[upper.tri(Bmat)] <- NA
for (i in 1:4) {Bmat[i,i] <- NA}

# Identify the B matrix of the SVAR model, and set it back to Cholesky
svar_model <- SVAR(var_model, Bmat = Bmat, max.iter =1000)
svar_model$B <- P_chol

# Set seed for reproducibility
set.seed(123)

# Define the number of trials
n_Q <- 50000  # Number of Q matrices to generate

# Define sign restrictions matrix
sign_restrictions <- matrix(c(
  NA,  NA,  NA, -1,  # Variable 1: Positive sign expected for Shock 1, Negative for Shock 2, No restriction for Shock 3, Positive for Shock 4
  NA,  NA,  NA,  1,  # Variable 2: Negative sign expected for Shock 1, Positive for Shock 2, No restriction for Shock 3, Negative for Shock 4
  NA,  NA,  NA,  1,  # Variable 3: No sign restrictions for Shock 1 and 2, Positive for Shock 3, No restriction for Shock 4
  NA,  NA,  NA,  NA   # Variable 4: Positive sign expected for Shock 1, Negative for Shock 2, No restriction for Shock 3, Positive for Shock 4
), nrow = 4, byrow = TRUE)

# Set the number of periods to check in the IRF (1 = only contemporaneous effects)
periods_to_check <- 1 # You can adjust this value as needed >= 1 # Reference paper and BdE model = 1

# Set the number of periods in the IRF horizon (number of periods = horizon + 1)
horizon <- 20 # Ensure that horizon >= periods_to_check -1 # Reference paper and BdE model = 12

# Set the sign restrictions to contemparenous/short term (FALSE) or long-term (TRUE)
long_term = TRUE

if (long_term) {
  horizon_to_check <- 0
  periods_to_check_lt <- 1
} else {
  horizon_to_check <- periods_to_check
  periods_to_check_lt <- periods_to_check
}

# Set the percentiles for the upper and lower bounds in the plot
percentiles <- c(5, 95)  ## en la funci?n para el plot esta definido 14 y 86

valid_B_matrices <- list()
valid_irf_list <- list()

for (i in 1:n_Q) {
  Q <- generate_Q(n = ncol(var_model$y))
  B_candidate <- P_chol %*% Q
  
  irf_result_svar <- compute_irfs_svar(svar_model, B_candidate, horizon_to_check)
  
  if (check_sign_restrictions_svar(irf_result_svar, sign_restrictions, periods_to_check_lt)) {
    valid_B_matrices <- append(valid_B_matrices, list(B_candidate))
    irf_result_svar <- compute_irfs_svar(svar_model, B_candidate, horizon)
    valid_irf_list <- append(valid_irf_list, list(irf_result_svar$irf))
  }
}

cat("Number of valid B matrices:", length(valid_B_matrices), "\n")

# Calculate the IRFs of the VAR model
irf_var <- vars::irf(var_model, n.ahead = horizon, ortho = TRUE, boot = FALSE)$irf

# Compute the summary statistics for the valid IRFs
irf_summary_svar <- compute_irf_summary_svar(valid_irf_list, percentiles)

# Set the shock index and plot title and x-axis label
shock_index <- 4
plot_title <- "IRF for a shock on the E3CI variable using a VAR model or a SVAR model with sign restrictions"
x_axis_label <- "Horizon period (in quarters)"
summary_type <- "median"

# Plot the IRFs with upper and lower bounds

time_horizon <- 1:(horizon + 1)

# Plot for the 1st variable
var_index <- 1
overlay_irf <- irf_var[[shock_index]][, var_index] * 100  # Scale overlay IRF
y_axis_label <- "GDP Annual Growth Rate (%)"
base_plot <- plot_irf_svar(irf_summary_svar, var_index, shock_index, summary_type, percentiles, plot_title, x_axis_label, y_axis_label, y_min=NULL, y_max=NULL)
final_plot <- base_plot +
  geom_line(aes(x = time_horizon, y = overlay_irf), color = "red", size = 1.5, linetype = "dashed") +  # Make VAR IRF more prominent
  labs(subtitle = "Dashed line represents VAR model IRF") +  # Add subtitle for clarity
  theme(legend.position = "none")  # Remove legend since it's not necessary

print(final_plot)

# #####################
# ORDER OF VARIABLES
# #####################

# To check the order of data
transformed_data <- transformed_data [c(3, 1, 2, 4)]

# Define the number of lags
p = 2 # Reference paper and BdE model

# Estimate the VAR model
var_model <- VAR(transformed_data, p, type = "const")

# Extract the covariance matrix of residuals (Sigma)
sigma <- summary(var_model)$covres

# Cholesky decomposition of the covariance matrix of residuals
P_chol <- t(chol(sigma))

# Define a B matrix with NA values for SVAR model initiliazation (on top of Cholesky)
Bmat <- P_chol
Bmat[upper.tri(Bmat)] <- NA
for (i in 1:4) {Bmat[i,i] <- NA}

# Identify the B matrix of the SVAR model, and set it back to Cholesky
svar_model <- SVAR(var_model, Bmat = Bmat, max.iter =1000)
svar_model$B <- P_chol

# Set seed for reproducibility
set.seed(123)

# Define the number of trials
n_Q <- 50000  # Number of Q matrices to generate

# Define sign restrictions matrix (updated according to the new order of variables)
sign_restrictions <- matrix(c(
  NA,  NA,  NA,  1,  
  NA,  NA,  NA,  -1,  
  NA,  NA,  NA,  1,  
  NA,  NA,  NA,  NA   
), nrow = 4, byrow = TRUE)

# Set the number of periods to check in the IRF (1 = only contemporaneous effects)
periods_to_check <- 1 # You can adjust this value as needed >= 1 # Reference paper and BdE model = 1

# Set the number of periods in the IRF horizon (number of periods = horizon + 1)
horizon <- 12 # Ensure that horizon >= periods_to_check -1 # Reference paper and BdE model = 12

# Set the sign restrictions to contemparenous/short term (FALSE) or long-term (TRUE)
long_term = FALSE

if (long_term) {
  horizon_to_check <- 0
  periods_to_check_lt <- 1
} else {
  horizon_to_check <- periods_to_check
  periods_to_check_lt <- periods_to_check
}

# Set the percentiles for the upper and lower bounds in the plot
percentiles <- c(5, 95)  ## en la funci?n para el plot esta definido 14 y 86

valid_B_matrices <- list()
valid_irf_list <- list()

for (i in 1:n_Q) {
  Q <- generate_Q(n = ncol(var_model$y))
  B_candidate <- P_chol %*% Q
  
  irf_result_svar <- compute_irfs_svar(svar_model, B_candidate, horizon_to_check)
  
  if (check_sign_restrictions_svar(irf_result_svar, sign_restrictions, periods_to_check_lt)) {
    valid_B_matrices <- append(valid_B_matrices, list(B_candidate))
    irf_result_svar <- compute_irfs_svar(svar_model, B_candidate, horizon)
    valid_irf_list <- append(valid_irf_list, list(irf_result_svar$irf))
  }
}

cat("Number of valid B matrices:", length(valid_B_matrices), "\n")

# Calculate the IRFs of the VAR model
irf_var <- vars::irf(var_model, n.ahead = horizon, ortho = TRUE, boot = FALSE)$irf

# Compute the summary statistics for the valid IRFs
irf_summary_svar <- compute_irf_summary_svar(valid_irf_list, percentiles)

# Set the shock index and plot title and x-axis label
shock_index <- 4
plot_title <- "IRF for a shock on the E3CI variable using a VAR model or a SVAR model with sign restrictions"
x_axis_label <- "Horizon period (in quarters)"
summary_type <- "median"

# Plot the IRFs with upper and lower bounds
time_horizon <- 1:(horizon + 1)

# Plot for the 2nd variable (due to the new order of variables)
var_index <- 2
overlay_irf <- irf_var[[shock_index]][, var_index] * 100  # Scale overlay IRF
y_axis_label <- "GDP Annual Growth Rate (%)"
base_plot <- plot_irf_svar(irf_summary_svar, var_index, shock_index, summary_type, percentiles, plot_title, x_axis_label, y_axis_label, y_min=NULL, y_max=NULL)
final_plot <- base_plot +
  geom_line(aes(x = time_horizon, y = overlay_irf), color = "red", size = 1.5, linetype = "dashed") +  # Make VAR IRF more prominent
  labs(subtitle = "Dashed line represents VAR model IRF") +  # Add subtitle for clarity
  theme(legend.position = "none")  # Remove legend since it's not necessary

print(final_plot)
