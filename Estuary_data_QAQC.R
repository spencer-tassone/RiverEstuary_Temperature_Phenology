library(ggplot2)
library(dplyr)
library(lubridate)
library(SWMPr)

rm(list = ls())
dev.off()

path <- "D:/School/USGSdata/TempPhenology/NERRSdat/217117.zip"

acesp <- import_local(path,'acespwq')
apaeb <- import_local(path,'apaebwq')
apaes <- import_local(path,'apaeswq')
cbvtc <- import_local(path,'cbvtcwq')
delbl <- import_local(path,'delblwq')
delsl <- import_local(path,'delslwq')
elkap <- import_local(path,'elkapwq')
elksm <- import_local(path,'elksmwq')
grbgb <- import_local(path,'grbgbwq')
hudtn <- import_local(path,'hudtnwq')
hudts <- import_local(path,'hudtswq')
jacb6 <- import_local(path,'jacb6wq')
jacb9 <- import_local(path,'jacb9wq')
jacba <- import_local(path,'jacbawq')
jacne <- import_local(path,'jacnewq')
job09 <- import_local(path,'job09wq')
narpc <- import_local(path,'narpcwq')
niwol <- import_local(path,'niwolwq')
niwta <- import_local(path,'niwtawq')
nocec <- import_local(path,'nocecwq')
nocrc <- import_local(path,'nocrcwq')
owcwm <- import_local(path,'owcwmwq')
pdbby <- import_local(path,'pdbbywq')
soswi <- import_local(path,'soswiwq')
tjros <- import_local(path,'tjroswq')
welht <- import_local(path,'welhtwq')
welin <- import_local(path,'welinwq')
wkbfr <- import_local(path,'wkbfrwq')
wkbwb <- import_local(path,'wkbwbwq')

acesp <- qaqc(acesp,qaqc_keep = c("-2","0","4","5"))
apaeb <- qaqc(apaeb,qaqc_keep = c("-2","0","4","5"))
apaes <- qaqc(apaes,qaqc_keep = c("-2","0","4","5"))
cbvtc <- qaqc(cbvtc,qaqc_keep = c("-2","0","4","5"))
delbl <- qaqc(delbl,qaqc_keep = c("-2","0","4","5"))
delsl <- qaqc(delsl,qaqc_keep = c("-2","0","4","5"))
elkap <- qaqc(elkap,qaqc_keep = c("-2","0","4","5"))
elksm <- qaqc(elksm,qaqc_keep = c("-2","0","4","5"))
grbgb <- qaqc(grbgb,qaqc_keep = c("-2","0","4","5"))
hudtn <- qaqc(hudtn,qaqc_keep = c("-2","0","4","5"))
hudts <- qaqc(hudts,qaqc_keep = c("-2","0","4","5"))
jacb6 <- qaqc(jacb6,qaqc_keep = c("-2","0","4","5"))
jacb9 <- qaqc(jacb9,qaqc_keep = c("-2","0","4","5"))
jacba <- qaqc(jacba,qaqc_keep = c("-2","0","4","5"))
jacne <- qaqc(jacne,qaqc_keep = c("-2","0","4","5"))
job09 <- qaqc(job09,qaqc_keep = c("-2","0","4","5"))
narpc <- qaqc(narpc,qaqc_keep = c("-2","0","4","5"))
niwol <- qaqc(niwol,qaqc_keep = c("-2","0","4","5"))
niwta <- qaqc(niwta,qaqc_keep = c("-2","0","4","5"))
nocec <- qaqc(nocec,qaqc_keep = c("-2","0","4","5"))
nocrc <- qaqc(nocrc,qaqc_keep = c("-2","0","4","5"))
owcwm <- qaqc(owcwm,qaqc_keep = c("-2","0","4","5"))
pdbby <- qaqc(pdbby,qaqc_keep = c("-2","0","4","5"))
soswi <- qaqc(soswi,qaqc_keep = c("-2","0","4","5"))
tjros <- qaqc(tjros,qaqc_keep = c("-2","0","4","5"))
welht <- qaqc(welht,qaqc_keep = c("-2","0","4","5"))
welin <- qaqc(welin,qaqc_keep = c("-2","0","4","5"))
wkbfr <- qaqc(wkbfr,qaqc_keep = c("-2","0","4","5"))
wkbwb <- qaqc(wkbwb,qaqc_keep = c("-2","0","4","5"))

acesp <- acesp[,1:2]
apaeb <- apaeb[,1:2]
apaes <- apaes[,1:2]
cbvtc <- cbvtc[,1:2]
delbl <- delbl[,1:2]
delsl <- delsl[,1:2]
elkap <- elkap[,1:2]
elksm <- elksm[,1:2]
grbgb <- grbgb[,1:2]
hudtn <- hudtn[,1:2]
hudts <- hudts[,1:2]
jacb6 <- jacb6[,1:2]
jacb9 <- jacb9[,1:2]
jacba <- jacba[,1:2]
jacne <- jacne[,1:2]
job09 <- job09[,1:2]
narpc <- narpc[,1:2]
niwol <- niwol[,1:2]
niwta <- niwta[,1:2]
nocec <- nocec[,1:2]
nocrc <- nocrc[,1:2]
owcwm <- owcwm[,1:2]
pdbby <- pdbby[,1:2]
soswi <- soswi[,1:2]
tjros <- tjros[,1:2]
welht <- welht[,1:2]
welin <- welin[,1:2]
wkbfr <- wkbfr[,1:2]
wkbwb <- wkbwb[,1:2]

acesp$site_nm <- "acesp"
apaeb$site_nm <- "apaeb"
apaes$site_nm <- "apaes"
cbvtc$site_nm <- "cbvtc"
delbl$site_nm <- "delbl"
delsl$site_nm <- "delsl"
elkap$site_nm <- "elkap"
elksm$site_nm <- "elksm"
grbgb$site_nm <- "grbgb"
hudtn$site_nm <- "hudtn"
hudts$site_nm <- "hudts"
jacb6$site_nm <- "jacb6"
jacb9$site_nm <- "jacb9"
jacba$site_nm <- "jacba"
jacne$site_nm <- "jacne"
job09$site_nm <- "job09"
narpc$site_nm <- "narpc"
niwol$site_nm <- "niwol"
niwta$site_nm <- "niwta"
nocec$site_nm <- "nocec"
nocrc$site_nm <- "nocrc"
owcwm$site_nm <- "owcwm"
pdbby$site_nm <- "pdbby"
soswi$site_nm <- "soswi"
tjros$site_nm <- "tjros"
welht$site_nm <- "welht"
welin$site_nm <- "welin"
wkbfr$site_nm <- "wkbfr"
wkbwb$site_nm <- "wkbwb"

# unexplained and unflagged doubling in water temperature in < 24 hours, therefore removed

nocrc[720054:720395,2] <- NA

# unflagged error where sensor came back online ~10C warmer than what was likely

elkap[11194:11232,2] <- NA

all_sites_hf <- rbind(acesp,apaeb,apaes,cbvtc,delbl,delsl,elkap,elksm,grbgb,hudtn,
                      hudts,jacb6,jacb9,jacba,jacne,job09,narpc,niwol,niwta,nocec,
                      nocrc,owcwm,pdbby,soswi,tjros,welht,welin,wkbfr,wkbwb)

rm(list=setdiff(ls(), "all_sites_hf"))

all_sites_daily <- all_sites_hf %>%
  mutate(year = year(datetimestamp),
         date = date(datetimestamp)) %>%
  group_by(site_nm,date) %>%
  summarise(meanTemp = round(mean(temp, na.rm = T),2)) %>%
  mutate_at("meanTemp", ~ifelse(is.nan(.), NA, .))

st <- as.Date("1996-01-01")
en <- as.Date("2021-12-31")
ll <- as.data.frame(rep(seq(st,en, by = "1 day"),29))
colnames(ll)[1] <- "date"
zz <- as.data.frame(sort(rep(c('acesp','apaeb','apaes','cbvtc','delbl','delsl','elkap','elksm','grbgb','hudtn',
            'hudts','jacb6','jacb9','jacba','jacne','job09','narpc','niwol','niwta','nocec',
            'nocrc','owcwm','pdbby','soswi','tjros','welht','welin','wkbfr','wkbwb'),9497)))
colnames(zz)[1] <- "site_nm"
aa <- cbind(zz,ll)
all_sites_daily <- merge(all_sites_daily,aa,by = c('site_nm','date'), all = TRUE)

rm(list=setdiff(ls(), c("all_sites_hf", "all_sites_daily")))

### Linear interpolate for data gaps less than or equal to 3 days

library(zoo)

interpolate <- function(df){
  df$Wtemp_int <- na.approx(df$meanTemp, maxgap = 3, na.rm = T)
  return(df)
}

all_sites_daily <- interpolate(all_sites_daily)
all_sites_daily$meanTemp <- ifelse(is.na(all_sites_daily$meanTemp),all_sites_daily$Wtemp_int, all_sites_daily$meanTemp)

# jacb9 given lat and long ~900 m south due to daymet coverage
# pdbby given longitude ~900 m east due to daymet coverage
lat_long <- data.frame(matrix(ncol = 3, nrow = 29))
colnames(lat_long)[1:3] <- c("site","lat","long")
lat_long$site <- c('acesp','apaeb','apaes','cbvtc','delbl','delsl','elkap','elksm','grbgb','hudtn',
                      'hudts','jacb6','jacb9','jacba','jacne','job09','narpc','niwol','niwta','nocec',
                      'nocrc','owcwm','pdbby','soswi','tjros','welht','welin','wkbfr','wkbwb')
lat_long$lat <- c(32.528,29.787222,29.787222,37.414986,39.388755,39.084977,
                 36.845733,36.818056,43.072222,42.036544,42.027039,39.5079,
                 39.49,39.593661,39.547881,17.943056,41.640361,33.349361,
                 33.299175,33.939861,34.156028,41.3825,48.496139,43.282464,
                 32.568333,43.298347,43.320089,30.416667,30.383333)
lat_long$long <- c(-80.36144,-84.881111,-84.881111,-76.71442,-75.635997,-75.460585,
                   -121.753678,-121.739397,-70.869444,-73.925325,-73.925958,
                   -74.338519,-74.38,-74.551511,-74.460769,-66.238583,
                   -71.340879,-79.188881,-79.256042,-77.941139,-77.849972,-82.514722,
                   -122.495,-124.320424,-117.131278,-70.587094,-70.563442,-87.816667,-87.833333)

# Grab meteorological data from Daymet web services to build water temperature multiple linear regressions
# https://daac.ornl.gov/
# https://www.nature.com/articles/s41597-021-00973-0#code-availability
library(daymetr)

write.table(lat_long, paste0(tempdir(),"/lat_long.csv"),
            sep = ",",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

met_dat <- download_daymet_batch(file_location = paste0(tempdir(),
                                                        "/lat_long.csv"),
                                 start = 1996,
                                 end = 2021,
                                 internal = TRUE,
                                 simplify = TRUE)
library(tidyr)
met_dat_wide <- spread(met_dat, measurement, value)
met_dat_wide <- met_dat_wide %>%
  mutate(AirTempMean = (tmax..deg.c. + tmin..deg.c.)/2,
         date = as.Date(paste(year, yday, sep = "-"), "%Y-%j"),
         totalRadiation = (met_dat_wide$srad..W.m.2.*met_dat_wide$dayl..s.)/1000000) # calculation based on daymetr website https://daymet.ornl.gov/overview
met_dat_wide <- met_dat_wide[,c(1,16,9,11:15,17)]
colnames(met_dat_wide)[1] <- "site_nm"
wmet <- merge(met_dat_wide, all_sites_daily, by = c("site_nm","date"), all = T)
wmet <- wmet[,1:10]
colnames(wmet)[c(8,10)] <- c("atemp_mean","wtemp_mean")

# remove leap days
remove_leap <- as.Date(c("1996-02-29","2000-02-29","2004-02-29",
                         "2008-02-29","2012-02-29","2016-02-29","2020-02-29"))
wmet <- wmet[!wmet$date %in% remove_leap,]

# day of year that does not recognize leap day
wmet <- wmet %>% 
  mutate(doy = day(date),
         month = month(date),
         year = year(date)) %>% 
  group_by(year, month, site_nm) %>%
  mutate(doy = doy - lag(doy, default = 0)) %>%
  group_by(year,site_nm) %>%
  mutate(DoY = cumsum(doy)) %>%
  select(-month)
wmet <- as.data.frame(wmet)

library(broom)

models_fit <- wmet %>%
  group_by(site_nm) %>%
  do(model = glance(lm(wtemp_mean~prcp..mm.day.+vp..Pa.+atemp_mean+tmax..deg.c.+tmin..deg.c.+totalRadiation+doy, data = .))) %>%
  unnest(model)
models_fit$adj.r.squared <- round(models_fit$adj.r.squared, digits = 2)
good_model_fits <- models_fit[models_fit$adj.r.squared >= 0.80,] # 28 sites with r-square >= 0.80
round(mean(good_model_fits$adj.r.squared),digits = 2) # answer is 0.90
round(sd(good_model_fits$adj.r.squared),digits = 2) # answer is 0.04
wmet <- wmet[wmet$site_nm %in% good_model_fits$site_nm,]

library(purrr)

f <- function (.fit, .new_data) {
  predict(.fit, newdata = .new_data)
}

set.seed(8992)
wmet <- wmet %>%
  nest(data = -site_nm) %>% 
  mutate(fit  = map(data, ~ lm(wtemp_mean~prcp..mm.day.+vp..Pa.+atemp_mean+tmax..deg.c.+tmin..deg.c.+totalRadiation+doy, data = .x)),
         yhat = map2(.x = fit, .y = data, f)) %>% 
  unnest(cols = c(data, yhat)) %>% 
  select(-fit)

wmet$predWtemp <- ifelse(is.na(wmet$wtemp_mean),wmet$yhat,wmet$wtemp_mean)
wmet$predWtemp <- round(wmet$predWtemp, digits = 2)
wmet$corWtemp <- ifelse(wmet$predWtemp < 0,0,wmet$predWtemp)

remove_extra <- as.Date("2022-01-01")
wmet <- wmet[!wmet$date %in% remove_extra,]

ggplot(wmet, aes(x = date, y = corWtemp)) +
  geom_line() +
  facet_wrap(~site_nm)

#########################################################

dat <- wmet[,c(1,2,16)]
colnames(dat)[3] <- "temperature"

dat_table10 <- dat %>%
  mutate(year = year(date),
         doy = yday(date)) %>%
  group_by(site_nm, year) %>%
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

ggplot(dat_table10, aes(x = year, y = first_above10)) +
  stat_smooth(method = 'lm', se = F) +
  geom_point() +
  xlab("") +
  ylab(expression(paste("First Day of Year above 10 ", degree, "C"))) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black", angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(~site_nm)

setwd("D:/School/USGSdata/TempPhenology")
# write.csv(dat,'usgs_wtemp_phenology.csv', row.names = FALSE)

### middle doy below 10 C

dat_mid_below10 <- dat_table10 %>%
  drop_na(mid_below10)

mid_below_fit10 <- dat_mid_below10 %>%
  group_by(site_nm) %>%
  do(model = tidy(lm(mid_below10~year, na.action=na.exclude, data = .))) %>%
  unnest(model)
toDelete <- seq(0,nrow(mid_below_fit10),2)
mid_below_fit10 <- mid_below_fit10[toDelete,]
mid_below_fit10$slope_perdecade = round(mid_below_fit10$estimate * 10)
mid_below_fit10$p.value <- round(mid_below_fit10$p.value, 3)
sum(mid_below_fit10$p.value <= 0.05, na.rm = T)
mid_below10_sig_sites <- mid_below_fit10 %>%
  group_by(site_nm) %>%
  filter(p.value <= 0.05)


### middle doy above 10 C

dat_mid_above10 <- dat_table10 %>%
  drop_na(mid_above10)

mid_above_fit10 <- dat_mid_above10 %>%
  group_by(site_nm) %>%
  do(model = tidy(lm(mid_above10~year, na.action=na.exclude, data = .))) %>%
  unnest(model)
toDelete <- seq(0,nrow(mid_above_fit10),2)
mid_above_fit10 <- mid_above_fit10[toDelete,]
mid_above_fit10$slope_perdecade = round(mid_above_fit10$estimate * 10)
mid_above_fit10$p.value <- round(mid_above_fit10$p.value, 3)
sum(mid_above_fit10$p.value <= 0.05)
mid_above10_sig_sites <- mid_above_fit10 %>%
  group_by(site_nm) %>%
  filter(p.value <= 0.05)

### magnitude of highest 5 day average

dat_5d_max <- dat_table10 %>%
  drop_na(max_temp_5day)

dat_5d_max_fit <- dat_5d_max %>%
  group_by(site_nm) %>%
  do(model = tidy(lm(max_temp_5day~year, na.action=na.exclude, data = .))) %>%
  unnest(model)
toDelete <- seq(0,nrow(dat_5d_max_fit),2)
dat_5d_max_fit <- dat_5d_max_fit[toDelete,]
dat_5d_max_fit$slope_perdecade = round(dat_5d_max_fit$estimate * 10)
dat_5d_max_fit$p.value <- round(dat_5d_max_fit$p.value, 3)
sum(dat_5d_max_fit$p.value <= 0.05)
dat_5d_max_sig_sites <- dat_5d_max_fit %>%
  group_by(site_nm) %>%
  filter(p.value <= 0.05)

### timing (doy) of highest 5 day average

dat_5day_max_doy <- dat_table10 %>%
  drop_na(max_temp_doy_5day)

dat_5d_max_doy_fit <- dat_5day_max_doy %>%
  group_by(site_nm) %>%
  do(model = tidy(lm(max_temp_doy_5day~year, na.action=na.exclude, data = .))) %>%
  unnest(model)
toDelete <- seq(0,nrow(dat_5d_max_doy_fit),2)
dat_5d_max_doy_fit <- dat_5d_max_doy_fit[toDelete,]
dat_5d_max_doy_fit$slope_perdecade = round(dat_5d_max_doy_fit$estimate * 10)
dat_5d_max_doy_fit$p.value <- round(dat_5d_max_doy_fit$p.value, 3)
sum(dat_5d_max_doy_fit$p.value <= 0.05)
dat_5d_max_doy_sig_sites <- dat_5d_max_doy_fit %>%
  group_by(site_nm) %>%
  filter(p.value <= 0.05)

dat_table20 <- dat %>%
  mutate(year = year(date),
         doy = yday(date)) %>%
  group_by(site_nm, year) %>%
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
  group_by(site_nm) %>%
  do(model = tidy(lm(mid_below20~year, na.action=na.exclude, data = .))) %>%
  unnest(model)
toDelete <- seq(0,nrow(mid_below20_fit),2)
mid_below20_fit <- mid_below20_fit[toDelete,]
mid_below20_fit$slope_perdecade = round(mid_below20_fit$estimate * 20)
mid_below20_fit$p.value <- round(mid_below20_fit$p.value, 3)
sum(mid_below20_fit$p.value <= 0.05)
mid_below20_sig_sites <- mid_below20_fit %>%
  group_by(site_nm) %>%
  filter(p.value <= 0.05)

### middle doy above 20 C

dat_mid_above20 <- dat_table20 %>%
  drop_na(mid_above20)

mid_above20_fit <- dat_mid_above20 %>%
  group_by(site_nm) %>%
  do(model = tidy(lm(mid_above20~year, na.action=na.exclude, data = .))) %>%
  unnest(model)
toDelete <- seq(0,nrow(mid_above20_fit),2)
mid_above20_fit <- mid_above20_fit[toDelete,]
mid_above20_fit$slope_perdecade = round(mid_above20_fit$estimate * 20)
mid_above20_fit$p.value <- round(mid_above20_fit$p.value, 3)
sum(mid_above20_fit$p.value <= 0.05, na.rm = T)
mid_above20_sig_sites <- mid_above20_fit %>%
  group_by(site_nm) %>%
  filter(p.value <= 0.05)

NROW(mid_below10_sig_sites) # 7
NROW(mid_above10_sig_sites) # 0
NROW(dat_5d_max_sig_sites) # 12
NROW(dat_5d_max_doy_sig_sites) # 1
NROW(mid_below20_sig_sites) # 9
NROW(mid_above20_sig_sites) # 4

mid_below10_sig_sites$type = "mid below 10C"
mid_above10_sig_sites$type = "mid above 10C"
dat_5d_max_sig_sites$type = "5day max temp"
dat_5d_max_doy_sig_sites$type = "5day max temp DoY"
mid_below20_sig_sites$type = "mid below 20C"
mid_above20_sig_sites$type = "mid above 20C"

output_table <- rbind(mid_below10_sig_sites,mid_above10_sig_sites,
                      dat_5d_max_sig_sites,dat_5d_max_doy_sig_sites,
                      mid_below20_sig_sites,mid_above20_sig_sites)
output_table <- output_table[,c(1,8,3:6)]
colnames(output_table)[3:5] <- c("slope","slope_se","t_stat")
write.csv(output_table,"nerrs_phenology_sigOutput.csv", row.names = FALSE)

test <- merge(dat_table10,dat_table20, by = c("site_nm","year"))
write.csv(test,"nerrs_phenology_metrics.csv", row.names = FALSE)
