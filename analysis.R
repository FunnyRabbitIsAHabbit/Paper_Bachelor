# Plain R ----------
# setwd(getSrcDirectory()[1])
# RStudio ----------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Libraries ----------

library(dplyr)
library(openxlsx)
library(stargazer)
library(pastecs)
library(ggplot2)
library(forecast)
library(ecm)
library(gtools)
library(stringr)
library(lmtest)
library(urca)
library(vars)
library(tsDyn)


# Define string formatting operator ----------

`%format%` <- function(x, y) {
    
    do.call(sprintf, c(list(x), y))
}


# Modeling functions ----------

do_granger <- function(data_to_test,
                       num_of_possible_lags,
                       alpha_value = 0.1) {
    vars_to_test <- names(data_to_test)
    x <- permutations(length(vars_to_test), 2)
    
    for (i in 1:NROW(x)) {
        x1 <- x[i, 1]
        x2 <- x[i, 2]
        names <- c(vars_to_test[x1],
                   vars_to_test[x2])
        for (j in 1:length(num_of_possible_lags)) {
            test <- grangertest(formula = as.formula("%s ~ %s" %format% names),
                                order = num_of_possible_lags[j],
                                data = data_to_test)
            if (test$`Pr(>F)`[2] < alpha_value) {
                print(test)
            }
        }
    }
}

better_to_formula <- function(var_names,
                              max_num_lags,
                              dep_var_name,
                              other_no_lag_varsnames) {
    
    
    len <- length(var_names)
    possible_lags_matrix <- matrix(nrow = max_num_lags,
                                   ncol = 1)
    for (i in 1:max_num_lags) {
        possible_lags_matrix[i, ] <- paste(c(rep(1:max_num_lags, len = i),
                                             rep(0, len = max_num_lags - i)), collapse = "")
    }
    
    x <- permutations(max_num_lags, len, possible_lags_matrix, repeats.allowed = TRUE)
    N <- NROW(x)
    to_return_df <- data.frame(N = 1:N, formula = NA)
    for (i in 1:N) {
        mid_str_out <- ""
        for (j in 1:len) {
            lag_vector <- strsplit(x[i, j], NULL)[[1]]
            for (char in lag_vector) {
                if (char != "0") {
                    mid_str_out <- paste(mid_str_out, paste(var_names[j], "lag_", char, sep = ""), sep = "+")
                }
            }
        }
        if (!NA %in% other_no_lag_varsnames) {
            to_return_df$formula[N = i] <- paste(paste(dep_var_name, "~",
                                                       sep = " "),
                                                 paste(mid_str_out,
                                                       paste(other_no_lag_varsnames,
                                                             collapse = " + "),
                                                       sep = " + "))
        } else {
            to_return_df$formula[N = i] <- paste(paste(dep_var_name, "~",
                                                       sep = " "),
                                                 mid_str_out)
        }
        
    }
    
    return(to_return_df$formula)
}

to_formula <- function(name, nums, names, factor_name) {
    str_out <- NA
    for (i in 1:stringr::str_length(nums)) {
        current <- ifelse(substr(nums, i, i) == '0',
                          NA,
                          names[i])
        if (!is.na(current)) {
            if (!is.na(str_out)) {
                str_out <- paste(str_out, current, sep = '+')
            } else {
                str_out <- current
            }
            
        }
        
    }
    form <- paste(name, '~', str_out, '+', factor_name)
    return(as.formula(form))
}

give_num_combination_str_df <- function(list_to_replicate,
                                        num_reps_for_list_item) {
    
    num_combs <- NROW(num_reps_for_list_item)
    str_out_df <- data.frame(N = 1:num_combs, str_var = NA, str_res = NA)
    all_str_df <- data.frame(N = 1:num_combs ^ 2, str = NA)
    
    for (i in 1:num_combs) {
        str_out <- ''
        for (j in 1:length(list_to_replicate)) {
            new_str <- paste(rep(paste(list_to_replicate[, j],
                                       collapse = ''),
                                 num_reps_for_list_item[i, j]),
                             collapse = '')
            str_out <- paste(str_out, new_str, sep = '')
        }
        
        str_out_df$str_var[N = i] <- str_out
        str_out_df$str_res[N = i] <- str_out
    }
    
    perms <- gtools::permutations(num_combs, 2, repeats.allowed = TRUE)
    
    for (i in 1:NROW(perms)) {
        all_str_df$str[N = i] <- paste(str_out_df$str_var[perms[i, 1]],
                                       str_out_df$str_res[perms[i, 2]], sep = '')
    }
    
    return(all_str_df)
}

give_num_combination_str_df_intermediary <- function(num_of_vars_with_res,
                                                     num_of_possible_lags,
                                                     num_of_possible_res_lags) {
    lags_str_var_res <- data.frame(c(1, 0, 0),
                                   c(1, 2, 0),
                                   c(1, 2, 3))
    
    xyz <- gtools::permutations(num_of_vars_with_res + 1, num_of_vars_with_res, 0:num_of_vars_with_res,
                                repeats.allowed = TRUE)
    xyz <- as.data.frame.matrix(xyz)
    condition_new <- rowSums(xyz) == num_of_vars_with_res
    num_repeats <- data.frame(row.names = 1:sum(condition_new))
    
    for (i in 1:num_of_vars_with_res) {
        num_repeats <- cbind(num_repeats, xyz[, i][condition_new])
    }
    
    names(num_repeats) <- 1:NCOL(num_repeats)
    df_choose_lags <- give_num_combination_str_df(lags_str_var_res, num_repeats)
    num_perms <- NROW(df_choose_lags) * length(lags_str_var_res)
    all_str_df <- data.frame(row.names = 1:num_perms)
    
    i <- 0
    for (j in 1:length(lags_str_var_res)) {
        for (k in 1:NROW(df_choose_lags)) {
            i <- i + 1
            all_str_df$str[i] <- paste(paste(lags_str_var_res[, j], collapse = ''),
                                       df_choose_lags$str[k], collapse = '', sep = '')
        }
    }
    
    return(all_str_df)
}

ecm_build <- function(all_str_df,
                      acf_plot_title_ecm,
                      fitted_plot_title_ecm,
                      fits,
                      X) {
    min_aic <- 0
    
    for (j in 1:NROW(all_str_df)) {
        
        explaining_var_names <- names(X)[-c(1)]
        numbers <- all_str_df$str[j]
        formula <- to_formula(name = names(X)[c(1)],
                              nums = numbers,
                              names = explaining_var_names,
                              factor_name = 'crisis')
        print(formula)
        ecm_fit <- lm(data = X, formula = formula)
        
        if (!NA %in% ecm_fit$coefficients) {
            
            fits$aic[j] <- AIC(ecm_fit)
            min_aic <- min(fits$aic[j], min_aic)
            fits$number[j] <- j
            fits$coeffs[j] <- ecm_fit$coefficients
            fits$res[j] <- ecm_fit$residuals
            
            if (fits$aic[j] == min_aic) {
                aic_best <- ecm_fit
            }
            
        } else {
            print(paste("NA here", j, "/", NROW(all_str_df)))
        }
    }
    
    if (exists("aic_best")) {
        
        print(summary(aic_best))
        
        a <- acf(aic_best$residuals, plot = FALSE, lag.max = lag_max)
        plot_ict_ecm_res <- ggplot(data = data.frame(Lag = a$lag, ACF = a$acf),
                                   aes(x = Lag, y = ACF)) +
            geom_line() +
            ggtitle(acf_graph_title, subtitle = "ECM fit residuals")
        filename <- paste(country, acf_plot_title_ecm, sep = "_")
        pdf(filename)
        print(plot_ict_ecm_res)
        dev.off()
        
        plot_ecm_fit <- ggplot(data = aic_best,
                               aes(x = years[5]:years[length(years)],
                                   y = aic_best$model$delta_log_gdp, colour = '1')) +
            geom_line() +
            geom_point(aes(y = aic_best$fitted.values,
                           colour = '2')) +
            xlab('Year') + ylab('Value') +
            scale_colour_manual(name = 'Series',
                                labels = c('observed', 'fitted'),
                                values = c("red", "blue")) +
            ggtitle("Observed VS Fitted", subtitle = "%s. First difference of logGDP. ECM." %format% c(country))
        filename <- paste(country, fitted_plot_title_ecm, sep = "_")
        pdf(filename)
        print(plot_ecm_fit)
        dev.off()
    }
}

ecm_build_by_formula <- function(formula_str_vector,
                                 acf_plot_title_ecm,
                                 fitted_plot_title_ecm,
                                 df,
                                 fits) {
    N <- NROW(formula_str_vector)
    min_aic <- 0
    
    for (j in 1:N) {
        ecm_fit <- lm(data = df, formula = formula_str_vector[j])
        
        if (!NA %in% ecm_fit$coefficients) {
            
            fits$aic[j] <- AIC(ecm_fit)
            min_aic <- min(fits$aic[j], min_aic)
            fits$number[j] <- j
            fits$coeffs[j] <- ecm_fit$coefficients
            fits$res[j] <- ecm_fit$residuals
            
            if (fits$aic[j] == min_aic) {
                aic_best <- ecm_fit
            }
            
        } else {
            print(paste("NA here", j, "/", N))
        }
    }
    
    if (exists("aic_best")) {
        
        print(summary(aic_best))
        
        a <- acf(aic_best$residuals, plot = FALSE, lag.max = lag_max)
        plot_ict_ecm_res <- ggplot(data = data.frame(Lag = a$lag, ACF = a$acf),
                                   aes(x = Lag, y = ACF)) +
            geom_line() +
            ggtitle(acf_graph_title, subtitle = "ECM fit residuals")
        filename <- paste(country, acf_plot_title_ecm, sep = "_")
        pdf(filename)
        print(plot_ict_ecm_res)
        dev.off()
        
        plot_ecm_fit <- ggplot(data = aic_best,
                               aes(x = years[5]:years[length(years)],
                                   y = aic_best$model$delta_log_gdp, colour = '1')) +
            geom_line() +
            geom_point(aes(y = aic_best$fitted.values,
                           colour = '2')) +
            xlab('Year') + ylab('Value') +
            scale_colour_manual(name = 'Series',
                                labels = c('observed', 'fitted'),
                                values = c("red", "blue")) +
            ggtitle("Observed VS Fitted", subtitle = "%s. First difference of logGDP. ECM." %format% c(country))
        filename <- paste(country, fitted_plot_title_ecm, sep = "_")
        pdf(filename)
        print(plot_ecm_fit)
        dev.off()
    }
    
}

# Set input parameters ----------

countries <- c("UScuip", "USrd", "USse")
years <- 1990:2019
lag_max <- years[length(years)] - years[1]
acf_graph_title <- "Auto-Correlation Function"

# Main function ----------

main <- function(country, years, lag_max) {
    # Get data ----------
    
    data_file_country <- "%s_data.xlsx" %format% c(country)
    sheet_names <- openxlsx::getSheetNames(data_file_country)
    analysis_data <- list()
    
    for (name in sheet_names) {
        N <- lag_max + 2
        analysis_data[[name]] <- t(openxlsx::read.xlsx(data_file_country, sheet = name)[1, 2:N])
    }
    
    new_data_country <- data.frame(year = years)
    new_data_country$ict <- analysis_data$ict * analysis_data$gdp_defl_lcu
    new_data_country$gdp <- analysis_data$gdp_lcu * analysis_data$gdp_defl_lcu
    new_data_country$l <- analysis_data$labor_total * analysis_data$wages_lcu * analysis_data$gdp_defl_lcu
    new_data_country$k <- analysis_data$gross_cap_form_lcu * analysis_data$gdp_defl_lcu
    
    desc <- pastecs::stat.desc(new_data_country)
    # View(new_data_country)
    openxlsx::write.xlsx(desc, file = '%s_desc.xlsx' %format% c(country), row.names = TRUE)
    
    d <- data.frame(year = new_data_country$year, value = new_data_country$gdp, name = 'gdp')
    for (column in names(new_data_country)[-c(1)]) {
        d <- rbind(d, data.frame(year = new_data_country$year,
                                 value = new_data_country[[column]],
                                 name = column))
    }
    names(d) <- c('year', 'value', 'name')
    
    # plot_data <- ggplot(data = d, aes(x = year, y = value)) +
    #     geom_point() +
    #     xlab('Year') + ylab('USD') +
    #     facet_wrap(~name, scales = 'free_y') +
    #     ggtitle("Data", subtitle = country)
    # filename <- paste(country, "plot_data.pdf", sep = "_")
    # pdf(filename)
    # print(plot_data)
    # dev.off()
    
    # Modeling ----------
    
    fit <- lm(data = new_data_country, log(gdp) ~ log(k) + log(l) + ict) # Target model (6)
    num_of_possible_lags <- 1:3
    # print(summary(fit))
    
    
    # fit_plot <- ggplot(data = new_data_country,
    #                    aes(x = year)) +
    #     geom_line(aes(y = log(gdp), colour = '1')) +
    #     geom_point(aes(y = fit$fitted.values, colour = '2')) +
    #     scale_colour_manual(name = 'Series',
    #                         labels = c('observed', 'fitted'),
    #                         values = c("red", "blue")) +
    #     ggtitle("Observed VS Fitted", subtitle = "%s GDP model" %format% c(country))
    # filename <- paste(country, "fit_plot.pdf", sep = "_")
    # pdf(filename)
    # print(fit_plot)
    # dev.off()
    
    # Plot ACF fit ----------
    # a <- acf(fit$residuals, plot = FALSE, lag.max = lag_max)
    # plot <- ggplot(data = data.frame(Lag = a$lag, ACF = a$acf), aes(x = Lag, y = ACF)) +
    #     geom_line() +
    #     ggtitle(acf_graph_title, subtitle = "%s GDP model residuals" %format% c(country))
    # filename <- paste(country, "plot.pdf", sep = "_")
    # pdf(filename)
    # print(plot)
    # dev.off()
    
    # Plot ACF for input data ----------
    
    # a <- acf(log(new_data_country$gdp), plot = FALSE, lag.max = lag_max)
    # plot_gdp <- ggplot(data = data.frame(Lag = a$lag, ACF = a$acf), aes(x = Lag, y = ACF)) +
    #     geom_line() +
    #     ggtitle(acf_graph_title, subtitle = "%s logGDP" %format% c(country))
    # filename <- paste(country, "plot_log_gdp.pdf", sep = "_")
    # pdf(filename)
    # print(plot_gdp)
    # dev.off()
    # 
    # a <- acf(log(new_data_country$k), plot = FALSE, lag.max = lag_max)
    # plot_k<- ggplot(data = data.frame(Lag = a$lag, ACF = a$acf), aes(x = Lag, y = ACF)) +
    #     geom_line() +
    #     ggtitle(acf_graph_title, subtitle = "%s logCapital" %format% c(country))
    # filename <- paste(country, "plot_log_k.pdf", sep = "_")
    # pdf(filename)
    # print(plot_k)
    # dev.off()
    # 
    # a <- acf(log(new_data_country$l), lag = 1, plot = FALSE, lag.max = lag_max)
    # plot_l <- ggplot(data = data.frame(Lag = a$lag, ACF = a$acf), aes(x = Lag, y = ACF)) +
    #     geom_line() +
    #     ggtitle(acf_graph_title, subtitle = "%s logLabor" %format% c(country))
    # filename <- paste(country, "plot_log_l.pdf", sep = "_")
    # pdf(filename)
    # print(plot_l)
    # dev.off()
    # 
    # a <- acf(diff(log(new_data_country$gdp), lag = 1), plot = FALSE, lag.max = lag_max)
    # plot_gdp <- ggplot(data = data.frame(Lag = a$lag, ACF = a$acf), aes(x = Lag, y = ACF)) +
    #     geom_line() +
    #     ggtitle(acf_graph_title, subtitle = "%s d.logGDP" %format% c(country))
    # filename <- paste(country, "plot_d_log_gdp.pdf", sep = "_")
    # pdf(filename)
    # print(plot_gdp)
    # dev.off()
    # 
    # a <- acf(diff(log(new_data_country$k), lag = 1), plot = FALSE, lag.max = lag_max)
    # plot_k<- ggplot(data = data.frame(Lag = a$lag, ACF = a$acf), aes(x = Lag, y = ACF)) +
    #     geom_line() +
    #     ggtitle(acf_graph_title, subtitle = "%s d.logCapital" %format% c(country))
    # filename <- paste(country, "plot_d_log_k.pdf", sep = "_")
    # pdf(filename)
    # print(plot_k)
    # dev.off()
    # 
    # a <- acf(diff(log(new_data_country$l), lag = 1), plot = FALSE, lag.max = lag_max)
    # plot_l <- ggplot(data = data.frame(Lag = a$lag, ACF = a$acf), aes(x = Lag, y = ACF)) +
    #     geom_line() +
    #     ggtitle(acf_graph_title, subtitle = "%s d.logLabor" %format% c(country))
    # filename <- paste(country, "plot_d_log_l.pdf", sep = "_")
    # pdf(filename)
    # print(plot_l)
    # dev.off()
    # 
    # a <- acf(diff(new_data_country$ict, lag = 1), plot = FALSE, lag.max = lag_max)
    # plot_ict <- ggplot(data = data.frame(Lag = a$lag, ACF = a$acf), aes(x = Lag, y = ACF)) +
    #     geom_line() +
    #     ggtitle(acf_graph_title, subtitle = "%s d.ICT" %format% c(country))
    # filename <- paste(country, "plot_d_ict.pdf", sep = "_")
    # pdf(filename)
    # print(plot_ict)
    # dev.off()
    # 
    # a <- acf(diff(log(new_data_country$ict), lag = 1), plot = FALSE, lag.max = lag_max)
    # plot_ict <- ggplot(data = data.frame(Lag = a$lag, ACF = a$acf), aes(x = Lag, y = ACF)) +
    #     geom_line() +
    #     ggtitle(acf_graph_title, subtitle = "%s d.logICT" %format% c(country))
    # filename <- paste(country, "plot_d_log_ict.pdf", sep = "_")
    # pdf(filename)
    # print(plot_ict)
    # dev.off()
    # 
    # a <- acf(log(new_data_country$ict), plot = FALSE, lag.max = lag_max)
    # plot_ict <- ggplot(data = data.frame(Lag = a$lag, ACF = a$acf), aes(x = Lag, y = ACF)) +
    #     geom_line() +
    #     ggtitle(acf_graph_title, subtitle = "%s logICT" %format% c(country))
    # filename <- paste(country, "plot_log_ict.pdf", sep = "_")
    # pdf(filename)
    # print(plot_ict)
    # dev.off()
    
    # Pair co-integration ----------
    # coint_fit_k <- lm(data = new_data_country, log(gdp) ~ log(k))
    # coint_fit_l <- lm(data = new_data_country, log(gdp) ~ log(l))
    # coint_fit_ict <- lm(data = new_data_country, log(gdp) ~ ict)
    # u <- data.frame(k = coint_fit_k$residuals,
    #                 l = coint_fit_l$residuals,
    #                 ict = coint_fit_ict$residuals)
    # z <- data.frame(k = log(new_data_country$k),
    #                 l = log(new_data_country$l),
    #                 ict = new_data_country$ict)
    # coint_a <- data.frame(k = acf(u$k, plot = FALSE, lag.max = lag_max)$acf,
    #                       l = acf(u$l, plot = FALSE, lag.max = lag_max)$acf,
    #                       ict = acf(u$ict, plot = FALSE, lag.max = lag_max)$acf,
    #                       lag = acf(u$ict, plot = FALSE, lag.max = lag_max)$lag)
    # plot_ict_coint <- ggplot(data = coint_a, aes(x = lag)) +
    #     geom_line(aes(y = k, colour = '1')) +
    #     geom_line(aes(y = l, colour = '2')) +
    #     geom_line(aes(y = ict, colour = '3')) +
    #     scale_colour_manual(name = 'Series',
    #                         labels = c('logCapital', 'logLabor', 'ICT'),
    #                         values = c("red", "blue", "purple")) +
    #     ggtitle(acf_graph_title,
    #             subtitle = 'USA logGDP on ... Co-integraion regression residuals')
    # filename <- paste(country, "plot_ict_coint.pdf", sep = "_")
    # pdf(filename)
    # print(plot_ict_coint)
    # dev.off()
    
    # Total ECM build ----------
    
    
    # X <- data.frame(diff(log(new_data_country$gdp), lag = 1),
    #                 diff(log(new_data_country$k), lag = 1),
    #                 diff(log(new_data_country$l), lag = 1),
    #                 diff(fit$residuals, lag = 1))
    # cols <- c('delta_log_gdp', 'delta_log_k', 'delta_log_l',
    #           'coint_res')
    # names(X) <- cols
    

    # for (col in cols) {
    #     for (i in num_of_possible_lags) {
    #         new_col <- paste(col, 'lag', i, sep = '_')
    #         X[[new_col]] <- ecm::lagpad(X[[col]], k = i)
    #     }
    # 
    # }
    # 
    # X$delta_ict_lag_1 <- ecm::lagpad(diff(new_data_country$ict, lag = 1))
    # X$delta_ict_lag_2 <- ecm::lagpad(X$delta_ict_lag_1)
    # X$delta_ict_lag_3 <- ecm::lagpad(X$delta_ict_lag_2)
    # X <- na.omit(X[-c(2:4)])
    # X$crisis <- as.factor(ifelse(row.names(X) == '2008' | row.names(X) == '2009', 1, 0))
    # View(X)

    # ecm_fits <- list()
    # acf_plot_title_ecm <- "plot_ict_total_ecm_res.pdf"
    # fitted_plot_title_ecm <- "plot_total_ecm_fit.pdf"
    # all_str_formula <- better_to_formula(c("delta_log_gdp_",
    #                                        "delta_log_k_",
    #                                        "delta_log_l_",
    #                                        "coint_res_",
    #                                        "delta_ict_"),
    #                                      3,
    #                                      "delta_log_gdp",
    #                                      c("crisis"))
    # ecm_build_by_formula(all_str_formula,
    #                      acf_plot_title_ecm,
    #                      fitted_plot_title_ecm,
    #                      X,
    #                      ecm_fits)
    
    # Granger causality -----------
    
    # data_to_test <- data.frame(log_gdp_d = c(NA, diff(log(new_data_country$gdp), lag = 1)),
    #                            log_capital_d = c(NA, diff(log(new_data_country$k), lag = 1)),
    #                            log_labor_d = c(NA, diff(log(new_data_country$l), lag = 1)),
    #                            ict_d = c(NA, diff(new_data_country$ict, lag = 1)))
    # data_to_test <- data_to_test[-c(1), ]
    # do_granger(data_to_test = data_to_test,
    #            num_of_possible_lags = num_of_possible_lags)
    
    # another_data_to_test <- data.frame(gdp_d = c(NA, diff(new_data_country$gdp, lag = 1)),
    #                                    capital_d = c(NA, diff(new_data_country$k, lag = 1)),
    #                                    labor_d = c(NA, diff(new_data_country$l, lag = 1)),
    #                                    ict_d = c(NA, diff(new_data_country$ict, lag = 1)))
    # another_data_to_test <- another_data_to_test[-c(1), ]
    # do_granger(data_to_test = another_data_to_test,
    #            num_of_possible_lags = num_of_possible_lags)
    
    # data_to_test <- data.frame(log_gdp_d = c(NA, diff(log(new_data_country$gdp), lag = 1)),
    #                            log_capital_d = c(NA, diff(log(new_data_country$k), lag = 1)),
    #                            log_labor_d = c(NA, diff(log(new_data_country$l), lag = 1)),
    #                            ict_d = c(NA, diff(log(new_data_country$ict), lag = 1)))
    # data_to_test <- data_to_test[-c(1), ]
    # do_granger(data_to_test = data_to_test,
    #            num_of_possible_lags = num_of_possible_lags)
    
    # Smaller ECMs build ----------
    
    # X_spec <- data.frame(diff(log(new_data_country$gdp), lag = 1),
    #                      diff(log(new_data_country$k), lag = 1),
    #                      diff(log(new_data_country$l), lag = 1),
    #                      diff(u$k, lag = 1),
    #                      diff(u$l, lag = 1),
    #                      diff(u$ict, lag = 1))
    # cols <- c('delta_log_gdp', 'delta_log_k', 'delta_log_l',
    #           'coint_res_k', 'coint_res_l', 'coint_res_ict')
    # names(X_spec) <- cols
    # 
    # for (col in cols) {
    #     for (i in 1:3) {
    #         new_col <- paste(col, 'lag', i, sep = '_')
    #         X_spec[[new_col]] <- ecm::lagpad(X_spec[[col]], k = i)
    #     }
    # 
    # }
    # 
    # X_spec$delta_ict_lag_1 <- ecm::lagpad(diff(new_data_country$ict, lag = 1))
    # X_spec$delta_ict_lag_2 <- ecm::lagpad(X_spec$delta_ict_lag_1)
    # X_spec$delta_ict_lag_3 <- ecm::lagpad(X_spec$delta_ict_lag_2)
    # X_spec <- na.omit(X_spec[-c(2:6)])
    # X_spec$crisis <- as.factor(ifelse(row.names(X_spec) == '2008' | row.names(X_spec) == '2009', 1, 0))
    # 
    # num_of_vars_with_res <- 3
    # fits <- list()
    # 
    # all_str_df <- give_num_combination_str_df_intermediary(num_of_vars_with_res,
    #                                                        num_of_possible_lags,
    #                                                        num_of_possible_res_lags)
    # acf_plot_title_ecm <- "plot_ict_ecm_res.pdf"
    # fitted_plot_title_ecm <- "plot_ecm_fit.pdf"
    # ecm_build(all_str_df, acf_plot_title_ecm, fitted_plot_title_ecm, fits, X_spec)
    
    # VECM ----------
    
    vecm_data <- data.frame(row.names = new_data_country$year)
    vecm_data$gdp <- new_data_country$gdp
    vecm_data$ict <- new_data_country$ict
    vecm_data$capital <- new_data_country$k
    vecm_data$labor <- new_data_country$l
    vecm_data_log <- log(vecm_data)
    optimal_VAR <- VARselect(vecm_data_log,
                             lag.max = num_of_possible_lags[length(num_of_possible_lags)],
                             type = "const")
    optimal_VAR_lag <- optimal_VAR$selection[1] - 1
    print(optimal_VAR_lag)
    coint_relations <- ca.jo(vecm_data_log,
                             type = "eigen",
                             ecdet = "const",
                             K = optimal_VAR_lag,
                             spec = "longrun")
    print(summary(coint_relations))
    vecm_fit <- cajorls(coint_relations,
                        r = 1)
    print(vecm_fit)
    print(summary(vecm_fit$rlm))
    
    
    to_var <- vec2var(coint_relations)
    
    # Plot shock tests ----------
    
    filename <- paste(country, "plot_irf_10_gdp_ict.pdf", sep = "_")
    pdf(filename)
    print(plot(irf(to_var, impulse = "gdp", response = "ict",
                   n.ahead = 10),
               main = "IRF", ylab = "lnICT", xlab = "lnGDP"))
    dev.off()
    
    filename <- paste(country, "plot_irf_100_gdp_ict.pdf", sep = "_")
    pdf(filename)
    print(plot(irf(to_var, impulse = "gdp", response = "ict",
                   n.ahead = 100),
               main = "IRF", ylab = "lnICT", xlab = "lnGDP"))
    dev.off()
    
    filename <- paste(country, "plot_irf_10_ict_gdp.pdf", sep = "_")
    pdf(filename)
    print(plot(irf(to_var, impulse = "ict", response = "gdp",
                   n.ahead = 10),
               main = "IRF", ylab = "lnGDP", xlab = "lnICT"))
    dev.off()
    
    filename <- paste(country, "plot_irf_100_ict_gdp.pdf", sep = "_")
    pdf(filename)
    print(plot(irf(to_var, impulse = "ict", response = "gdp",
                   n.ahead = 100),
               main = "IRF", ylab = "lnGDP", xlab = "lnICT"))
    dev.off()
    
    filename <- paste(country, "plot_irf_10_ict_ict.pdf", sep = "_")
    pdf(filename)
    print(plot(irf(to_var, impulse = "ict", response = "ict",
                   n.ahead = 10),
               main = "IRF", ylab = "lnICT", xlab = "lnICT"))
    dev.off()
    
    filename <- paste(country, "plot_irf_100_ict_ict.pdf", sep = "_")
    pdf(filename)
    print(plot(irf(to_var, impulse = "ict", response = "ict",
                   n.ahead = 100),
               main = "IRF", ylab = "lnICT", xlab = "lnICT"))
    dev.off()
    
    filename <- paste(country, "plot_fevd_10.pdf", sep = "_")
    pdf(filename)
    print(plot(fevd(to_var)))
    dev.off()
    
    plotres(coint_relations)
    
    # Tests -----------
    
    
    print(serial.test(to_var, type = c("PT.asymptotic")))
    print(arch.test(to_var))
    print(normality.test(to_var))

}



# Run ----------

main(country = countries[1], years = years, lag_max = lag_max)
# main(country = countries[2], years = years, lag_max = lag_max)
# main(country = countries[3], years = years, lag_max = lag_max)


