# adapted from https://otexts.com/fpp2/ts-objects.html

# objective: find a statistical confirmation of seasonality
# performed with either  stl(ts_data, s.window = "periodic")
# or a sarima mocdel auto.arima(ts_data)


# padding data ####
# to account for unequal sampling each year and missing samples
library(lubridate)

# Generate a sequence of dates from January 2020 to December 2022
date_sequence <- seq(ym("2020-01"), ym("2023-12"), by = "1 month") %>% 
  as.data.frame() %>% dplyr::rename("Date" = ".")
date_sequence$Date <- format(date_sequence$Date, "%Y-%m")

# add env data
lobau_data_combined <- read_excel("~/Desktop/R/Christina/poster_fungi/lobau_data_combined.xlsx")
lobau_data_combined$date_yyyy_mm_dd <- as.Date(lobau_data_combined$date_yyyy_mm_dd)
lobau_data_combined$date_yyyy_mm_dd <- format(lobau_data_combined$date_yyyy_mm_dd, "%Y-%m")

# join together
left_join(date_sequence, 
          lobau_data_combined[which(lobau_data_combined$well_id == "D15" & 
                                      lobau_data_combined$sample_type == "pumped_groundwater"),], 
          by = c("Date"= "date_yyyy_mm_dd") ) -> time_data

ts(time_data, start = 2020, frequency=12) -> time_data

library(ggfortify) # needed for autoplot
ggplot2::autoplot(time_data[, "oxygen_mg_L"]) +
  ggtitle("Seasonal changes in oxygen in D15") +
  ylab("Oxygen [mg/L]") +
  xlab("Year") +
  theme_minimal()


# ARIMA ####
library(forecast)
library(imputeTS)

imputeTS::na_ma(time_data[, "oxygen_mg_L"]) -> ts_no_na
forecast::auto.arima(ts_no_na, seasonal = TRUE) # output: ARIMA(1,1,0)(0,1,0)[12] 

# # interpretation # # 
# p = 0, d = 0, q = 0: These are the non-seasonal ARIMA parameters.
# P = 0, D = 1, Q = 0: These are the seasonal ARIMA parameters.
# [12] indicates a seasonal period of 12 (monthly data).

plot(ts_no_na)


# start from scratch with model
ts_no_na %>%
  decompose() %>%
  autoplot() 


forecast::forecast(ts_no_na, h=12) %>% plot()
 
# testing
library(seastests)
seastests::isSeasonal(ts_no_na, test="kw", freq=12)
seastests::kw(ts_no_na, freq = 12, diff = T, residuals = F, autoarima = T)
                      

# temperature ####
ggplot2::autoplot(time_data[, "temp_C"]) +
  ggtitle("Seasonal changes in temperature in D15") +
  ylab("Temperature [C]") +
  xlab("Year") +
  theme_minimal()

imputeTS::na_ma(time_data[, "temp_C"]) -> ts_no_na_temp
plot(ts_no_na_temp)

seastests::isSeasonal(ts_no_na_temp, test="kw", freq=12)
seastests::kw(ts_no_na_temp, freq = 12, diff = T, residuals = F, autoarima = T)

