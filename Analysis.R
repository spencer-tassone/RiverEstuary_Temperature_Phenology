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

dat <- dat %>%
  mutate(year = year(Date),
         doy = yday(Date)) %>%
  group_by(site_no, year) %>%
  mutate(gt_10 = temperature >= 10,
         lt_10 = temperature <= 10,
         # max_t = which.max(temperature),
         # min_t = which.min(temperature),
         peak_doy = doy[which.max(temperature)],
         below_peak = doy < peak_doy,
         after_peak = doy > peak_doy,
         run = data.table::rleid(lt_10)) %>%
  summarise(first_above = doy[below_peak & gt_10][1],
            sustain_above = first(doy[run == max(run[below_peak])]),
            first_below = doy[after_peak & lt_10][1],
            sustain_below = first(doy[run == max(run[after_peak])]), .groups = 'drop',
            mid_above = mean(first_above, sustain_above),
            mid_below = mean(first_below, sustain_below))
View(dat)

station_details <- read.csv('Station_Details.csv', colClasses = c(site_no = "character"))
dat <- merge(dat,station_details, by ="site_no")

ggplot(dat, aes(x = year, y = first_above)) +
  stat_smooth(method = 'lm', se = F) +
  geom_point() +
  xlab("") +
  ylab(expression(paste("First Day of Year above 10 ", degree, "C"))) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black", angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(~site_no)

setwd("D:/School/USGSdata/TempPhenology")
# write.csv(dat,'usgs_wtemp_phenology.csv', row.names = FALSE)

dat_mid_below <- dat %>%
  drop_na(mid_below)
dat_mid_above <- dat %>%
  drop_na(mid_above)

mid_below_fit <- dat_mid_below %>%
  group_by(site_no) %>%
  do(model = tidy(lm(mid_below~year, na.action=na.exclude, data = .))) %>%
  unnest(model)
toDelete <- seq(0,nrow(mid_below_fit),2)
mid_below_fit <- mid_below_fit[toDelete,]
mid_below_fit$slope_perdecade = round(mid_below_fit$estimate * 10)
mid_below_fit$p.value <- round(mid_below_fit$p.value, 3)
sum(mid_below_fit$p.value <= 0.05)
mid_below_sig_sites <- mid_below_fit %>%
  group_by(site_no) %>%
  filter(p.value <= 0.05)
mid_below_sig_sites <- merge(mid_below_sig_sites, station_details, by = "site_no")
View(mid_below_sig_sites)

mid_above_fit <- dat_mid_above %>%
  group_by(site_no) %>%
  do(model = tidy(lm(mid_above~year, na.action=na.exclude, data = .))) %>%
  unnest(model)
toDelete <- seq(0,nrow(mid_above_fit),2)
mid_above_fit <- mid_above_fit[toDelete,]
mid_above_fit$slope_perdecade = round(mid_above_fit$estimate * 10)
mid_above_fit$p.value <- round(mid_above_fit$p.value, 3)
sum(mid_above_fit$p.value <= 0.05)
mid_above_sig_sites <- mid_above_fit %>%
  group_by(site_no) %>%
  filter(p.value <= 0.05)
mid_above_sig_sites <- merge(mid_above_sig_sites, station_details, by = "site_no")
View(mid_above_sig_sites)
