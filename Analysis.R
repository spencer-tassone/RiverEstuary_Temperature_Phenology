library(ggplot2)
library(lubridate)
library(dplyr)
library(tidyverse)
library(dataRetrieval)

rm(list = ls())
dev.off()

setwd("D:/School/USGSdata/GitHub")

Wtemp_daily_dat <- read.csv('Wtemp_daily_dat.csv', colClasses = c(site_no = "character",
                                                                  Date = "Date"))

dat <- Wtemp_daily_dat[,c(1,2,24)]
colnames(dat)[3] <- "temperature"

dat_table10 <- dat %>%
  mutate(year = year(Date),
         doy = yday(Date)) %>%
  group_by(site_no, year) %>%
  mutate(gt_10 = temperature >= 10,
         lt_10 = temperature <= 10,
         peak_doy = doy[which.max(temperature)],
         below_peak = doy < peak_doy,
         after_peak = doy > peak_doy,
         run = data.table::rleid(lt_10),
         roll_temp5day = zoo::rollmean(temperature, k = 5, fill = NA, align = 'center')) %>%
  summarise(first_above10 = doy[below_peak & gt_10][1],
            sustain_above10 = first(doy[run == max(run[below_peak])]),
            first_below10 = doy[after_peak & lt_10][1],
            sustain_below10 = first(doy[run == max(run[after_peak])]), .groups = 'drop',
            mid_above10 = mean(first_above10, sustain_above10),
            mid_below10 = mean(first_below10, sustain_below10),
            max_temp = max(temperature),
            min_temp = min(temperature),
            max_temp_doy = doy[which.max(temperature)],
            min_temp_doy = doy[which.min(temperature)],
            max_temp_5day = round(max(roll_temp5day, na.rm = T),1),
            min_temp_5day = round(min(roll_temp5day, na.rm = T),1),
            max_temp_doy_5day = doy[which.max(roll_temp5day)],
            min_temp_doy_5day = doy[which.min(roll_temp5day)])

station_details <- read.csv('Station_Details.csv', colClasses = c(site_no = "character"))
dat_table10 <- merge(dat_table10,station_details, by ="site_no")

ggplot(dat_table10, aes(x = year, y = first_above10)) +
  stat_smooth(method = 'lm', se = F) +
  geom_point() +
  xlab("") +
  ylab(expression(paste("First Day of Year above 10 ", degree, "C"))) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black", angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(~site_no)

# setwd("D:/School/USGSdata/TempPhenology")
# write.csv(dat,'usgs_wtemp_phenology.csv', row.names = FALSE)

### middle doy below 10 C

dat_mid_below10 <- dat_table10 %>%
  drop_na(mid_below10)

mid_below_fit10 <- dat_mid_below10 %>%
  group_by(site_no) %>%
  do(model = tidy(lm(mid_below10~year, na.action=na.exclude, data = .))) %>%
  unnest(model)
toDelete <- seq(0,nrow(mid_below_fit10),2)
mid_below_fit10 <- mid_below_fit10[toDelete,]
mid_below_fit10$slope_perdecade = round(mid_below_fit10$estimate * 10)
mid_below_fit10$p.value <- round(mid_below_fit10$p.value, 3)
sum(mid_below_fit10$p.value <= 0.05)
mid_below10_sig_sites <- mid_below_fit10 %>%
  group_by(site_no) %>%
  filter(p.value <= 0.05)
mid_below10_sig_sites <- merge(mid_below10_sig_sites, station_details, by = "site_no")

### middle doy above 10 C

dat_mid_above10 <- dat_table10 %>%
  drop_na(mid_above10)

mid_above_fit10 <- dat_mid_above10 %>%
  group_by(site_no) %>%
  do(model = tidy(lm(mid_above10~year, na.action=na.exclude, data = .))) %>%
  unnest(model)
toDelete <- seq(0,nrow(mid_above_fit10),2)
mid_above_fit10 <- mid_above_fit10[toDelete,]
mid_above_fit10$slope_perdecade = round(mid_above_fit10$estimate * 10)
mid_above_fit10$p.value <- round(mid_above_fit10$p.value, 3)
sum(mid_above_fit10$p.value <= 0.05)
mid_above10_sig_sites <- mid_above_fit10 %>%
  group_by(site_no) %>%
  filter(p.value <= 0.05)
mid_above10_sig_sites <- merge(mid_above10_sig_sites, station_details, by = "site_no")

### magnitude of highest 5 day average

dat_5d_max <- dat_table10 %>%
  drop_na(max_temp_5day)

dat_5d_max_fit <- dat_5d_max %>%
  group_by(site_no) %>%
  do(model = tidy(lm(max_temp_5day~year, na.action=na.exclude, data = .))) %>%
  unnest(model)
toDelete <- seq(0,nrow(dat_5d_max_fit),2)
dat_5d_max_fit <- dat_5d_max_fit[toDelete,]
dat_5d_max_fit$slope_perdecade = round(dat_5d_max_fit$estimate * 10)
dat_5d_max_fit$p.value <- round(dat_5d_max_fit$p.value, 3)
sum(dat_5d_max_fit$p.value <= 0.05)
dat_5d_max_sig_sites <- dat_5d_max_fit %>%
  group_by(site_no) %>%
  filter(p.value <= 0.05)
dat_5d_max_sig_sites <- merge(dat_5d_max_sig_sites, station_details, by = "site_no")

### timing (doy) of highest 5 day average

dat_5day_max_doy <- dat_table10 %>%
  drop_na(max_temp_doy_5day)

dat_5d_max_doy_fit <- dat_5day_max_doy %>%
  group_by(site_no) %>%
  do(model = tidy(lm(max_temp_doy_5day~year, na.action=na.exclude, data = .))) %>%
  unnest(model)
toDelete <- seq(0,nrow(dat_5d_max_doy_fit),2)
dat_5d_max_doy_fit <- dat_5d_max_doy_fit[toDelete,]
dat_5d_max_doy_fit$slope_perdecade = round(dat_5d_max_doy_fit$estimate * 10)
dat_5d_max_doy_fit$p.value <- round(dat_5d_max_doy_fit$p.value, 3)
sum(dat_5d_max_doy_fit$p.value <= 0.05)
dat_5d_max_doy_sig_sites <- dat_5d_max_doy_fit %>%
  group_by(site_no) %>%
  filter(p.value <= 0.05)
dat_5d_max_doy_sig_sites <- merge(dat_5d_max_doy_sig_sites, station_details, by = "site_no")

dat_table20 <- dat %>%
  mutate(year = year(Date),
         doy = yday(Date)) %>%
  group_by(site_no, year) %>%
  mutate(gt_20 = temperature >= 20,
         lt_20 = temperature <= 20,
         peak_doy = doy[which.max(temperature)],
         below_peak = doy < peak_doy,
         after_peak = doy > peak_doy,
         run = data.table::rleid(lt_20)) %>%
  summarise(first_above20 = doy[below_peak & gt_20][1],
            sustain_above20 = first(doy[run == max(run[below_peak])]),
            first_below20 = doy[after_peak & lt_20][1],
            sustain_below20 = first(doy[run == max(run[after_peak])]), .groups = 'drop',
            mid_above20 = mean(first_above20, sustain_above20),
            mid_below20 = mean(first_below20, sustain_below20))

### middle doy below 20 C

dat_mid_below20 <- dat_table20 %>%
  drop_na(mid_below20)

mid_below20_fit <- dat_mid_below20 %>%
  group_by(site_no) %>%
  do(model = tidy(lm(mid_below20~year, na.action=na.exclude, data = .))) %>%
  unnest(model)
toDelete <- seq(0,nrow(mid_below20_fit),2)
mid_below20_fit <- mid_below20_fit[toDelete,]
mid_below20_fit$slope_perdecade = round(mid_below20_fit$estimate * 20)
mid_below20_fit$p.value <- round(mid_below20_fit$p.value, 3)
sum(mid_below20_fit$p.value <= 0.05)
mid_below20_sig_sites <- mid_below20_fit %>%
  group_by(site_no) %>%
  filter(p.value <= 0.05)
mid_below20_sig_sites <- merge(mid_below20_sig_sites, station_details, by = "site_no")

### middle doy above 20 C

dat_mid_above20 <- dat_table20 %>%
  drop_na(mid_above20)

mid_above20_fit <- dat_mid_above20 %>%
  group_by(site_no) %>%
  do(model = tidy(lm(mid_above20~year, na.action=na.exclude, data = .))) %>%
  unnest(model)
toDelete <- seq(0,nrow(mid_above20_fit),2)
mid_above20_fit <- mid_above20_fit[toDelete,]
mid_above20_fit$slope_perdecade = round(mid_above20_fit$estimate * 20)
mid_above20_fit$p.value <- round(mid_above20_fit$p.value, 3)
sum(mid_above20_fit$p.value <= 0.05, na.rm = T)
mid_above20_sig_sites <- mid_above20_fit %>%
  group_by(site_no) %>%
  filter(p.value <= 0.05)
mid_above20_sig_sites <- merge(mid_above20_sig_sites, station_details, by = "site_no")

NROW(mid_below10_sig_sites) # 14
NROW(mid_above10_sig_sites) # 3
NROW(dat_5d_max_sig_sites) # 3
NROW(dat_5d_max_doy_sig_sites) # 14
NROW(mid_below20_sig_sites) # 9
NROW(mid_above20_sig_sites) # 6

mid_below10_sig_sites$type = "mid below 10C"
mid_above10_sig_sites$type = "mid above 10C"
dat_5d_max_sig_sites$type = "5day max temp"
dat_5d_max_doy_sig_sites$type = "5day max temp DoY"
mid_below20_sig_sites$type = "mid below 20C"
mid_above20_sig_sites$type = "mid above 20C"

output_table <- rbind(mid_below10_sig_sites,mid_above10_sig_sites,
                      dat_5d_max_sig_sites,dat_5d_max_doy_sig_sites,
                      mid_below20_sig_sites,mid_above20_sig_sites)
output_table <- output_table[,c(8,1,9:14,3:6)]
