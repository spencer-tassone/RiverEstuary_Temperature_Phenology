library(ggplot2)
library(dplyr)
library(lubridate)
library(SWMPr)

rm(list = ls())

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

all_sites_hf <- rbind(acesp,apaeb,apaes,cbvtc,delbl,delsl,elkap,elksm,grbgb,hudtn,
                      hudts,jacb6,jacb9,jacba,jacne,job09,narpc,niwol,niwta,nocec,
                      nocrc,owcwm,pdbby,soswi,tjros,welht,welin,wkbfr,wkbwb)

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

### Linear interpolate for data gaps less than or equal to 2 days

library(zoo)

interpolate <- function(df){
  df$Wtemp_int <- na.approx(df$meanTemp, maxgap = 2, na.rm = T)
  return(df)
}

all_sites_daily <- interpolate(all_sites_daily)
all_sites_daily$meanTemp <- ifelse(is.na(all_sites_daily$meanTemp),all_sites_daily$Wtemp_int, all_sites_daily$meanTemp)

ggplot(all_sites_daily, aes(x = date, y = meanTemp)) +
  geom_point() +
  facet_wrap(~site_nm)

test1 <- all_sites_daily %>%
  group_by(site_nm) %>%
  summarise(count_na = sum(is.na(meanTemp)))
View(test1)

test2 <- all_sites_daily %>%
  mutate(year = year(date)) %>%
  group_by(site_nm,year) %>%
  summarise(count_na = sum(is.na(meanTemp)))
View(test2)
