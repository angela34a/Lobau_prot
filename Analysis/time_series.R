# adapted from https://otexts.com/fpp2/ts-objects.html



# input is a data frame with two columns
# col1 is Date, col2 is values in question
env_lobau_final_r$Date <- as.Date(env_lobau_final_r$Date)

# temp ####
env_lobau_final_r %>% 
  dplyr::select("Date", "temp_C") %>% 
  as.data.frame() -> data1

data
data1 %>%  as.data.frame()
class(data1$temp_C)


library(forecast)
AutoArimaModel=auto.arima( data1 )
AutoArimaModel