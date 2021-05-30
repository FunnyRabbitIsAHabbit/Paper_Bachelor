# Plain R ----------
# setwd(getSrcDirectory()[1])
# RStudio ----------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Libraries ----------

library(dplyr)
library(openxlsx)
library(ggplot2)
library(stringr)
library(zoo)
library(urca)
library(vars)
library(tsDyn)

`%format%` <- function(x, y) {
    
    do.call(sprintf, c(list(x), y))
}

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
    
    return(new_data_country)
    
}

main_non_currency_ict <- function(country, years, lag_max) {
    # Get data ----------
    
    data_file_country <- "%s_data.xlsx" %format% c(country)
    sheet_names <- openxlsx::getSheetNames(data_file_country)
    analysis_data <- list()
    
    for (name in sheet_names) {
        N <- lag_max + 2
        analysis_data[[name]] <- t(openxlsx::read.xlsx(data_file_country, sheet = name)[1, 2:N])
    }
    
    new_data_country <- data.frame(year = years)
    new_data_country$ict <- analysis_data$ict
    new_data_country$gdp <- analysis_data$gdp_lcu * analysis_data$gdp_defl_lcu
    new_data_country$l <- analysis_data$labor_total * analysis_data$wages_lcu * analysis_data$gdp_defl_lcu
    new_data_country$k <- analysis_data$gross_cap_form_lcu * analysis_data$gdp_defl_lcu
    
    return(new_data_country)
    
}



# Run ----------

countries <- c("UScuip", "USrd", "USse", "USfix")
years <- 1990:2019
lag_max <- years[length(years)] - years[1]
acf_graph_title <- "Auto-Correlation Function"
num_of_possible_lags <- 1:3

df1 <- main(country = countries[1], years = years, lag_max = lag_max)
df2 <- main(country = countries[2], years = years, lag_max = lag_max)
df3 <- main(country = countries[3], years = years, lag_max = lag_max)
df4 <- main_non_currency_ict(country = countries[4], years = years, lag_max = lag_max)
# View(df1)
# View(df2)
# View(df3)
# View(df4)

df1$ict <- ifelse(na.spline.default(df1$ict) < 0, NA, na.spline.default(df1$ict))
df2$ict <- ifelse(na.spline.default(df2$ict) < 0, NA, na.spline.default(df2$ict))
df3$ict <- ifelse(na.spline.default(df3$ict) < 0, NA, na.spline.default(df3$ict))
df4$ict <- ifelse(na.spline.default(df4$ict) < 0, NA, na.spline.default(df4$ict))

icts <- data.frame(df1$ict,
                   df2$ict,
                   df3$ict,
                   df4$ict)
for (col in names(icts)) {
    for (i in 1:NROW(icts[[col]])) {
        if (is.na(icts[[col]][i])) {
            icts[[col]][i] <- icts[[col]][!is.na(icts[[col]])][1]
        }
    }
    
}

icts0 <- as.data.frame(log(as.matrix(icts)))
names(icts0) <- countries
icts <- diff(icts0, lag = 1)
icts <- as.data.frame(icts)
names(icts) <- countries
# test_variance_df <- data.frame(value = unlist(icts),
#                                group = as.factor(rep(1:4, each = NROW(icts))))
# View(test_variance_df)
# 
# anova_model <- aov(data = test_variance_df, value ~ group)
# summary(anova_model)
# oneway.test(data = test_variance_df, value ~ group, var.equal = TRUE)
# cor(icts)
# cor.test(icts$UScuip, icts$USrd, method = "pearson")
# cor.test(icts$UScuip, icts$USse, method = "pearson")
# cor.test(icts$UScuip, icts$USfix, method = "pearson")
# cor.test(icts$USrd, icts$USse, method = "pearson")
# cor.test(icts$USrd, icts$USfix, method = "pearson")

other <- as.data.frame(log(as.matrix(data.frame(df1$gdp,
                                                df1$l,
                                                df1$k))))
names(other) <- c("USgdp", "USlabor", "UScapital")
df <- cbind(other, icts0)

optimal_VAR <- VARselect(df,
                         lag.max = num_of_possible_lags[length(num_of_possible_lags)],
                         type = "const")
optimal_VAR_lag <- optimal_VAR$selection[1] - 1
coint_relations <- ca.jo(df,
                         type = "eigen",
                         ecdet = "const",
                         K = optimal_VAR_lag,
                         spec = "longrun")
print(summary(coint_relations))
# vecm_fit <- cajorls(coint_relations,
#                     r = 3)
# print(vecm_fit)
# print(summary(vecm_fit$rlm))
# 
# to_var <- vec2var(coint_relations)
# print(serial.test(to_var, type = c("PT.asymptotic")))
# print(arch.test(to_var))
# print(normality.test(to_var))
# # plotres(coint_relations)
# plot(fevd(to_var))






