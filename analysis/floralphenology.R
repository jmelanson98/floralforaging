##### Get blueberry phenology
### April 14, 2026
### J Melanson


### Load in packages
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(ggplot2)

bombus_path = "/Users/jenna1/Documents/UBC/bombus_project/"

### Load in floral survey data
specimenData2022 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2022specimendata.csv"), sep = ",", header = T))
specimenData2023 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2023specimendata.csv"), sep = ",", header = T))
veg2022 = read.csv(paste0(bombus_path, "raw_data/2022vegetationdata.csv"))
veg2023 = read.csv(paste0(bombus_path, "raw_data/2023vegetationdata.csv"))
sample2022 = read.csv(paste0(bombus_path, "raw_data/2022sampledata.csv"))
sample2023 = read.csv(paste0(bombus_path, "raw_data/2023sampledata.csv"))

# get sites where we observed blueberry
blueberry2022 = veg2022[c("site", "round", "sample_point", "sample_id", "VACO")]
blueberry2023 = veg2023[c("site", "round", "sample_point", "sample_id", "VACO")]

filtered2022 = blueberry2022 %>% 
  group_by(sample_point) %>%
  filter(sum(VACO, na.rm = TRUE) > 0)
filtered2022$VACO[is.na(filtered2022$VACO)] = 0
filtered2023 = blueberry2023 %>%
  group_by(sample_point) %>%
  filter(sum(VACO, na.rm = TRUE) > 0)
filtered2023$VACO[is.na(filtered2023$VACO)] = 0

#add julian date to sample effort data frame
sample2022$date = paste(sample2022$day, sample2022$month, sample2022$year)
sample2022$date = gsub(" ", "", sample2022$date, fixed = TRUE)
sample2022$date <- as.POSIXlt(sample2022$date, format = "%d%b%y")
sample2022$julian_date = sample2022$date$yday
dates2022 = sample2022[c("round", "julian_date", "sample_point")]
filtered2022 = left_join(filtered2022, dates2022, by = c("sample_point", "round"))

sample2023$date = paste(sample2023$day, sample2023$month, sample2023$year)
sample2023$date = gsub(" ", "", sample2023$date, fixed = TRUE)
sample2023$date <- as.POSIXlt(sample2023$date, format = "%d%b%y")
sample2023$julian_date = sample2023$date$yday
dates2023 = sample2023[c("round", "julian_date", "sample_point")]
filtered2023 = left_join(filtered2023, dates2023, by = c("sample_point", "round"))


# plot by round
ggplot(data = filtered2023, aes(x = julian_date, y = VACO, colour = sample_point)) +
  geom_point() +
  geom_line() + 
  theme_bw() 


# get min and max
filtered2022_max = filtered2022 %>%
  filter(VACO>0)
print(c(min(filtered2022_max$julian_date), max(filtered2022_max$julian_date)))
print(c(min(filtered2022_max$round), max(filtered2022_max$round)))


filtered2023_max = filtered2023 %>%
  filter(VACO>0)
print(c(min(filtered2023_max$julian_date), max(filtered2023_max$julian_date)))
print(c(min(filtered2023_max$round), max(filtered2023_max$round)))



# get sites where we observed blackberry
blackberry2022 = veg2022[c("site", "round", "sample_point", "sample_id", "RULA", "RUAR")]
blackberry2023 = veg2023[c("site", "round", "sample_point", "sample_id", "RULA", "RUAR")]

filtered2022 = blackberry2022 %>% 
  group_by(sample_point) %>%
  filter(sum(RULA, na.rm = TRUE) + sum(RUAR, na.rm = TRUE) > 0)
filtered2022$RUAR[is.na(filtered2022$RUAR)] = 0
filtered2022$RULA[is.na(filtered2022$RULA)] = 0
filtered2023 = blackberry2023 %>% 
  group_by(sample_point) %>%
  filter(sum(RULA, na.rm = TRUE) + sum(RUAR, na.rm = TRUE) > 0)
filtered2023$RUAR[is.na(filtered2023$RUAR)] = 0
filtered2023$RULA[is.na(filtered2023$RULA)] = 0

#add julian date
filtered2022 = left_join(filtered2022, dates2022, by = c("sample_point", "round"))
filtered2023 = left_join(filtered2023, dates2023, by = c("sample_point", "round"))


# plot by round
ggplot(data = filtered2023, aes(x = julian_date, y = RUAR, colour = sample_point)) +
  geom_point() +
  geom_line() + 
  theme_bw() +
  theme(legend.position = "none")


# get min and max
filtered2022_max = filtered2022 %>%
  filter(RUAR>0 | RULA >0)
print(c(min(filtered2022_max$julian_date), max(filtered2022_max$julian_date)))
print(c(min(filtered2022_max$round), max(filtered2022_max$round)))


filtered2023_max = filtered2023 %>%
  filter(RUAR>0 | RULA>0)
print(c(min(filtered2023_max$julian_date), max(filtered2023_max$julian_date)))
print(c(min(filtered2023_max$round), max(filtered2023_max$round)))


