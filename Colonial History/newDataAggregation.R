#Librarires
library(dplyr)

#load data
folder <- "New colonial data/"

#settler mortality
setMort <- read.csv(paste0(folder,"SettlerMortality.csv"), sep=";", dec=",")%>%
  select(cname, ccodealp, ajr_settmort)%>%
  unique()%>%
  rename(Country = cname,
         iso3 = ccodealp)#%>%
  select(iso3,setMort =ajr_settmort)

#colonizer
col <- read.csv(paste0(folder,"colonies.csv"))%>%
  filter(!(European.colonial.power..grouped. %in% c("z. Multiple colonizers",
                                                 "zz. Colonizer",
                                                 "zzzz. No longer colonized")))%>%
  group_by(Entity)%>%
  slice_max(Year)%>%
  ungroup()%>%
  rename(Country = Entity,
         iso3 = Code)#%>%
  select(iso3,colonizer = European.colonial.power..grouped.)

#primary school enrollment 1900
prienr1900 <- read.csv(paste0(folder,"prienr.csv"))%>%
  filter(Year <= 1900)%>%
  filter(!(Code %in% c("","OWID_WRL")))%>%
  group_by(Code)%>%
  slice_max(Year)%>%
  ungroup()%>%
  rename(Country = Entity,
         iso3 = Code)#%>%
  select(iso3,prienr1900 = Combined.total.net.enrolment.rate..primary..both.sexes )


MergedDF <- merge(setMort, col, by=c("iso3"), all=TRUE) %>%
  merge(prienr1900, by=c("iso3"), all=TRUE)%>%
  select(iso3,Country,Country.x, Country.y)
           