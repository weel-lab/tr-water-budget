---
title: "Tres Rios Water Budget Data Processing"
output: html_notebook
---

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Cmd+Shift+Enter*. 

```{r}
#HOUSEKEEPING

# libraries ----
library('plyr') # always load plyr before dplyr
library("dplyr")
library("tidyr")
library("stringr")
library("lubridate")

#tr met station dat import and formatting
tr_met <- read.csv('/Volumes/GoogleDrive/My Drive/weel/tres rios/datasets/tr-water-budget/data/tr_metstation_data.csv', stringsAsFactors = FALSE)
tr_met[tr_met == ''] <- NA #change all missing values to NA
tr_met[,3:8] <- sapply(tr_met[,3:8], function(x) as.numeric(x))
tr_met <- tr_met %>%
  mutate(date = paste(date, hour, sep = " ")) %>% #combine date and hour columns
  mutate(date = mdy_h(date)) #read date column as date formatted via lubridate

#tr_met quality control
tr_met <- tr_met %>%
  mutate(air_temp_f = ifelse(air_temp_f >= 125, NA, air_temp_f)) %>%
  mutate(par_w.k2 = ifelse(par_w.k2 >= 1000, NA, par_w.k2)) %>%
  mutate(rh_percent = ifelse(rh_percent >= 101, NA, rh_percent)) %>%
  mutate(wind_speed_m.s = ifelse(wind_speed_m.s >= 40, NA, wind_speed_m.s)) %>%
  mutate(bar_press_hg = ifelse(bar_press_hg >= 35, NA, bar_press_hg))

#tr_met calculating additional variables and formats
tr_met <- tr_met %>%
  mutate(air_temp_c = (air_temp_f-32)*(5/9)) %>%
  mutate(rh_decimal = rh_percent/100)

#goodyear precip data import
daily_precip <- read.csv('/Volumes/GoogleDrive/My Drive/weel/tres rios/datasets/tr-water-budget/data/tr_precip_timeseries.csv')
daily_precip <- mutate(daily_precip, date = mdy(date)) #read date column as date formatted via lubridate

#scale daily goodyear precip data to month
goodyear_monthly_precip <- daily_precip %>%
  group_by(month(date), year(date)) %>%
  summarize(
    monthly_precip_mm = sum(prcp),
    monthly_precip_m3 = ((21*10000)*(0.001*monthly_precip_mm))
  ) %>%
  ungroup()
names(goodyear_monthly_precip)[1:2] <- c("month", "year") #rename date columns
goodyear_monthly_precip <- arrange(goodyear_monthly_precip, year, month)

#frw daily flow data import and edits
tr_daily_flows <- read.csv('/Volumes/GoogleDrive/My Drive/weel/tres rios/datasets/tr-water-budget/data/tr_waterflow_data.csv', stringsAsFactors = FALSE)
tr_daily_flows[tr_daily_flows == ''] <- NA #change all missing values to NA
tr_daily_flows[is.na(tr_daily_flows)] <- 0 #change all NA vla values to NA
tr_daily_flows <- mutate(tr_daily_flows, date = mdy(date)) #read combined date-hour column as date formatted via lubridate
names(tr_daily_flows)[2:5] <- c("influent_flow_mgd", "frw1_outflow_mgd", "frw2_outflow_mgd", "frw3_outflow_mgd") #uses name to append units to column names 
tr_daily_flows <- tr_daily_flows %>% #calculate additional variables and unit conversions from original CoP flow data
  mutate(total_outflow_mgd = frw1_outflow_mgd + frw2_outflow_mgd +frw3_outflow_mgd) %>%
  mutate(frw1_percent_outflow = frw1_outflow_mgd/total_outflow_mgd) %>%
  mutate(frw1_inflow_mgd = frw1_percent_outflow*influent_flow_mgd) %>%
  mutate(frw1_outflow_m3 = (frw1_outflow_mgd*1000000)*0.00379) %>% #unit conversion mgd -> m3
  mutate(frw1_inflow_m3 = (frw1_inflow_mgd*1000000)*0.00379) %>% #unit conversion mgd -> m3
  mutate(frw1_deficit_m3 =  frw1_inflow_m3 - frw1_outflow_m3)

#create a data frame to sum and store monthly flow data based on daily flow data (above)
tr_monthly_flows <- tr_daily_flows %>%
  group_by(month(date), year(date)) %>%
  summarize(
    monthly_frw1_outflow_m3 = sum(frw1_outflow_m3),
    monthly_frw1_inflow_m3 = sum(frw1_inflow_m3),
    monthly_frw1_deficit_m3 = sum(frw1_deficit_m3)
  ) %>%
  ungroup()
names(tr_monthly_flows)[1:2] <- c("month", "year") #rename date columns
tr_monthly_flows <- arrange(tr_monthly_flows, year, month)

#et met data import (bring in CoP met data to eventually use to calculate IRGApar and IRGAtemp)
hourly_ET <- read.csv('/Volumes/GoogleDrive/My Drive/weel/tres rios/datasets/tr-water-budget/data/tr_metstation_data.csv', stringsAsFactors = FALSE)
hourly_ET[hourly_ET == ''] <- NA # change all missing values to NA

#converting data classes & trimming dataset for ET calcs
hourly_ET[,3:8] <- sapply(hourly_ET[,3:8], function(x) as.numeric(as.character(x))) #convert all columns except date and hour to numeric
hourly_ET$hour <- as.character(hourly_ET$hour) #convert hour column to character (initially integer)
hourly_ET <- hourly_ET %>%
  mutate(date = paste(date, hour, sep = " ")) %>% #combine date and hour columns
  mutate(date = mdy_h(date)) #read combined date-hour column as date formatted via lubridate
#use select to remove all columns except those necessary for ET calculations
hourly_ET <- select(hourly_ET, date, hour, air_temp_f, par_w.k2)

#hourly biomass data import (stop-gap until I can do this part in R)
hourly_biomass <- read.csv('/Volumes/GoogleDrive/My Drive/weel/tres rios/datasets/tr-water-budget/data/tr_hourly_biomass.csv', stringsAsFactors = FALSE)
#converting data classes
hourly_biomass <- hourly_biomass %>%
  mutate(date = paste(date, hour, sep = " ")) %>% #combine date and hour columns
  mutate(date = mdy_h(date)) #read combined date-hour column as date formatted via lubridate


```

```{r ET models by species}
#models
irga_temp_model <- function(x) 
  ((0.819*hourly_ET$air_temp_f)+9.315)
irga_par_model <- function(x) 
  (0.769*(hourly_ET$par_w.k2)+88.55)
sam_et_model <- function(x)
  ((0.00089*hourly_ET$irga_temp)+(0.0000189*hourly_ET$irga_par)-0.017)
tlat_et_model <- function(x)
  ((0.000305*hourly_ET$irga_temp)+(0.00000616*hourly_ET$irga_par))
tdom_et_model <- function(x)
  ((0.000834*hourly_ET$irga_temp)+(0.0000177*hourly_ET$irga_par))
typha_et_model <- function(x)
  ((0.00032*hourly_ET$irga_temp)+(0.00000613*hourly_ET$irga_par)-0.0051)
stab_et_model <- function(x)
  ((0.00058*hourly_ET$irga_temp)+(0.00002*hourly_ET$irga_par))
sac_et_model <- function(x)
  ((0.00069*hourly_ET$irga_temp)+(0.00002*hourly_ET$irga_par))
sacstab_et_model <- function(x)
  ((0.00089*hourly_ET$irga_temp)+(0.0000189*hourly_ET$irga_par)-0.017)
scal_et_model <- function(x)
  ((0.000634*hourly_ET$irga_temp)+(0.000011*hourly_ET$irga_par)-0.0076)
```

```{r TIME SERIES ET}
#TIME SERIES ET

#calculate IRGA temp and PAR
hourly_ET <- hourly_ET %>%
  mutate(irga_temp = irga_temp_model()) %>%
  mutate(irga_par = ifelse(par_w.k2 > 1, irga_par_model(), 0))

#calculate ET for each species, as well as two combinations of species (typha and sacstab), in mmol H2O/g*s
hourly_ET <- hourly_ET %>%
  mutate(sam_et_mmolH2O.g.s = sam_et_model()) %>%
  mutate(tlat_et_mmolH2O.g.s = tlat_et_model()) %>%
  mutate(tdom_et_mmolH2O.g.s = tdom_et_model()) %>%
  mutate(typha_et_mmolH2O.g.s = typha_et_model()) %>%
  mutate(stab_et_mmolH2O.g.s = stab_et_model()) %>%
  mutate(sac_et_mmolH2O.g.s = sac_et_model()) %>%
  mutate(sacstab_et_mmolH2O.g.s = sacstab_et_model()) %>%
  mutate(scal_et_mmolH2O.g.s = scal_et_model())

#join tables to perform calcs. maybe not necessary but for now this avoids the issues with different total rows in each spreadsheet by joining by "date"
hourly_ET <- merge(hourly_ET, hourly_biomass, by = c("date")) %>%
  arrange(date)

#calculate ET for each species per second (multiply previous values by hourly biomass)
hourly_ET <- hourly_ET %>%
  mutate(sam_et_mmolH2O.s = (sam_et_mmolH2O.g.s * 1000) * sam_bio.kg) %>%
  mutate(typha_et_mmolH2O.s = (typha_et_mmolH2O.g.s * 1000) * typha_bio.kg) %>%
  mutate(sacstab_et_mmolH2O.s = (sacstab_et_mmolH2O.g.s * 1000) * sacstab_bio.kg) %>%
  mutate(scal_et_mmolH2O.s = (scal_et_mmolH2O.g.s * 1000) * scal_bio.kg)

#perform dimensional analysis to convert from mmolH2O/s to m3 H2O/hr
hourly_ET <- hourly_ET %>%
  mutate(sam_et_m3H2O.hr = (sam_et_mmolH2O.s * (18.015*3600)/1000000000)) %>%
  mutate(typha_et_m3H2O.hr = (typha_et_mmolH2O.s * (18.015*3600)/1000000000)) %>%
  mutate(sacstab_et_m3H2O.hr = (sacstab_et_mmolH2O.s * (18.015*3600)/1000000000)) %>%
  mutate(scal_et_m3H2O.hr = (scal_et_mmolH2O.s * (18.015*3600)/1000000000))

#convert all NA in ET values as a result of missing met data to 0
#hourly_ET[,7:32][is.na(hourly_ET[,7:32])] <- 0

#create a data frame to calculate and store daily_ET using the hourly_ET data
daily_ET <- hourly_ET %>%
  group_by(month(date), day(date), year(date)) %>%
  summarize(
    sam_et_m3H2O = sum(sam_et_m3H2O.hr),
    typha_et_m3H2O = sum(typha_et_m3H2O.hr),
    sacstab_et_m3H2O = sum(sacstab_et_m3H2O.hr),
    scal_et_m3H2O = sum(scal_et_m3H2O.hr)
  ) %>%
  ungroup() %>%
  mutate(allspp_et_m3H2O = sam_et_m3H2O + typha_et_m3H2O + sacstab_et_m3H2O + scal_et_m3H2O) #add new column that sums ET across all spp
names(daily_ET)[1:3] <- c("month", "day", "year") #rename first three columns

#date formatting for daily_ET
daily_ET <- daily_ET %>%
  mutate(date = paste(month, day, year, sep = "/")) %>%
  mutate(date = mdy(date)) %>%
  arrange(date)

#QA for daily_ET
daily_ET$sam_et_m3H2O <- ifelse(daily_ET$sam_et_m3H2O < 0, NA, daily_ET$sam_et_m3H2O)
daily_ET$typha_et_m3H2O <- ifelse(daily_ET$typha_et_m3H2O < 0, NA, daily_ET$typha_et_m3H2O)
daily_ET$sacstab_et_m3H2O <- ifelse(daily_ET$sacstab_et_m3H2O < 0, NA, daily_ET$sacstab_et_m3H2O)
daily_ET$scal_et_m3H2O <- ifelse(daily_ET$scal_et_m3H2O < 0, NA, daily_ET$scal_et_m3H2O)
daily_ET$allspp_et_m3H2O <- ifelse(daily_ET$allspp_et_m3H2O < 0, NA, daily_ET$allspp_et_m3H2O)

monthly_ET <- daily_ET %>%
  group_by(month, year) %>%
  summarize(
    sam_et_m3H2O = sum(sam_et_m3H2O, na.rm = TRUE),
    typha_et_m3H2O = sum(typha_et_m3H2O, na.rm = TRUE),
    sacstab_et_m3H2O = sum(sacstab_et_m3H2O, na.rm = TRUE),
    scal_et_m3H2O = sum(scal_et_m3H2O, na.rm = TRUE),
    allspp_et_m3H2O = sum(allspp_et_m3H2O, na.rm = TRUE)
  ) %>%
  arrange(year, month) %>%
  ungroup()

```

```{r}

#open water evap calculations

#calculate additional met variables that are used as inputs to penmann-monteith equation
hourly_water_evap <- tr_met %>%
  mutate(air_temp_K = air_temp_f + 273) %>%
  mutate(par_mj.m2.hr = par_w.k2*0.0036) %>%
  mutate(bar_press_kPa = bar_press_hg*3.39) %>%
  mutate(slope_vpsatcurve_kPa.K = 0.00000022*exp(0.062*air_temp_K)) %>%
  mutate(psych_constant_kpa.kg = (0.0016286*bar_press_kPa)/2.39) %>%
  mutate(vapor_press_deficit_k = ((1-rh_decimal)*exp(21.07-(5336/air_temp_K))))

#use penmann-monteith equation to calculate hourly water evap in mm.hr
hourly_water_evap <- hourly_water_evap %>%
  mutate(hourly_evap_mm = ((slope_vpsatcurve_kPa.K*par_mj.m2.hr)+(psych_constant_kpa.kg*(6.43*(1+0.536*wind_speed_m.s))*vapor_press_deficit_k))/(2.39*(slope_vpsatcurve_kPa.K+psych_constant_kpa.kg)))
hourly_water_evap[is.na(hourly_water_evap)] <- 0

#scale hourly evap calculations to days
daily_water_evap <- hourly_water_evap %>%
  group_by(month(date), day(date), year(date)) %>%
  summarize(
    daily_evap_mm = sum(hourly_evap_mm),
    daily_evap_m3 = (sum(hourly_evap_mm)*210)
  ) %>%
  ungroup()
names(daily_water_evap)[1:3] <- c("month", "day", "year") #rename first three columns
daily_water_evap <- arrange(daily_water_evap, year, month, day)

#scale daily evap sums to monthly sums
monthly_water_evap <- daily_water_evap %>%
  group_by(month, year) %>%
  summarize(
    monthly_evap_m3 = sum(daily_evap_m3)
  ) %>%
  arrange(year, month) %>%
  ungroup()



```

```{r}
#WATER BUDGET
#This chunk will pull together all of the various monthly data from above into one data frame for comparative graphing and output

water_budget_inputs <- list(monthly_ET, tr_monthly_flows, goodyear_monthly_precip, monthly_water_evap) #creates a list that stores all of the monthly data frames to use as inputes to water_budget

water_budget <- water_budget_inputs %>%
    Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by=c("month", "year")), .) %>%
  mutate(date = paste(month, "01", year, sep = "/")) %>% #make a new date column
  mutate(date = mdy(date)) %>% #reads new date column as date via lubridate package
  select(date, everything()) #moves date to the first column position

```

*EXPORT DATA*

```{r RESULTS EXPORT}

#write certain datasheets as csv to the current working directory

write.csv(water_budget, file = "/Volumes/GoogleDrive/My Drive/weel/tres rios/datasets/tr-water-budget/results/tr_water_budget.csv")


```



1) As of 8/23, hourly_ET only goes through "2016-07-31 01:00:00". So select for all dates before that
2) See if there is a way to get all dates consistently copied over from initial CoP data. In most cases the dates from those spreadsheets already contain time info embedded (though may not be consistent pre-2014). If so, then you cna eliminate manually updating in hours column in excel and edit code lines 24-38 to just extract mdyh from date column
 -building off of this, is there a way to run this just starting from the tr_weather_alldata.csv file?
3) need to see if there's a way to resolve issue in tr_metstation_data.csv - how to get hourly data pre-2012? this would allow me to just use tr_metstation_data as the starting point for hourly_et process

20101029 - with respect to item 3 above, this seems to cause me to need to bring in this random old timeseries_et file. can I just take the CoP data from there and paste it into the CoP met data?
need to notate new data source for precip - used to be godo year but now 
20101030 - THIS HAS BEEN FIXED! Eliminated the need for "timeseries_et", can now just use met station file for all ET calcs
-also, air_temp in tr_met is already in C, need to remove code that does this conversion and fix the labeling of this in the Evap calculations
-also need to find a way to generate a QC'd version of TR_metstation_data to use for ET calcs...otherwise this is causing huge spikes where there is a QC issue in hourly temperature data.

Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Cmd+Option+I*.

When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Cmd+Shift+K* to preview the HTML file).
