library(readxl)
library(dplyr)
library(tidyr)


#Load Data ---------------------------------------
# Define file paths
data_dir <- "Data/"
subfolder <- "Justus and william/"

files <- list(
  SPI_Time_Series = paste0(data_dir, "SPI Time Series 1990-2020 data file.xlsx"),
  PrimSchoolEnroll = paste0(data_dir, "PrimarySchoolEnrollment.csv"),
  SafelyManagedSaniation = paste0(data_dir, "SafelyManagedSanitation.csv"),
  gdp = paste0(data_dir,subfolder, "gdp_data.csv"),
  caloric_supply = paste0(data_dir,subfolder, "daily-per-capita-caloric-supply.csv"),
  healthy_life_expectancy = paste0(data_dir,subfolder, "Health Life Expectance (HALE) at birth.csv"),
  protein_supply = paste0(data_dir,subfolder, "protein_supply.csv"),
  child_mortality = paste0(data_dir,subfolder, "Under Five child mortality Rate (WHO).csv"),
  physicians = paste0(data_dir,subfolder, "Number of Physicians per 1000 people.csv"),
  hospital_beds = paste0(data_dir,subfolder, "Number of Hospital Beds per 1000 people.csv"),
  essential_health_services = paste0(data_dir,subfolder, "Coverage of Essential Health Services.xlsx"),
  basic_drinking_water = paste0(data_dir,subfolder, "People using at least basic drinking water services (% of population).csv"),
  safely_managed_drinking_water = paste0(data_dir,subfolder, "People using safely managed drinking water services (% of population).csv"),
  basic_handwashing = paste0(data_dir,subfolder, "People with basic handwashing facilities including soap and water (% of population).csv"),
  clean_fuels = paste0(data_dir,subfolder, "Access to clean fuels and technologies for cooking (% of population).csv"),
  house_price_income_ratio = paste0(data_dir,subfolder, "House price-to-income ratio.csv"),
  share_slums = paste0(data_dir,subfolder, "share_of_urban_population_living_in_slums.csv"),
  share_overcrowding = paste0(data_dir,subfolder, "percentage_population_living in overcrwoding conditions_OECD.csv"),
  intentional_homicide = paste0(data_dir,subfolder, "data_cts_intentional_homicide.xlsx"),
  ieaworldind = paste0(data_dir,subfolder, "WORLDIND.xlsx"),
  GDPUN = paste0(data_dir,subfolder, "GDPUN.xlsx")
)

#Load and prepare SPI Data --------------------------------
SPI_Time_Series <- read_excel(files$SPI_Time_Series, sheet = "1990-2020 Time-series data")[, 1:132] %>%
    select(-c(7:24,77:132)) %>% #Removes aggregated scores (SPI score, dimensions and component scores)
    setNames(.[1, ]) %>% # Rename columns using the first row
    slice(-1) %>%
    rename("SPI_year" = "SPI \r\nyear",
           "SPI_rank" = "SPI\r\nRank",
           "SPI_countrycode" = "SPI \r\ncountry \r\ncode") %>%
    filter(!is.na(Region))%>% #Keep countries only
    select(-c(SPI_rank,Status))%>% #Remove rank and status
    mutate(SPI_year = as.integer(SPI_year))

nIndicators <- ncol(SPI_Time_Series)-4
nYears <- length(unique(SPI_Time_Series$SPI_year))

#Filter out countries with high missingness
na <- SPI_Time_Series %>%
  group_by(Country) %>%
  summarise_all(~sum(is.na(.)))%>%
  mutate(total_NA = rowSums(across(where(is.numeric))))%>%
  mutate(frac_NA = total_NA/(nIndicators*nYears)) #Calculate fraction of NAs per country

#Inspect distribution
hist(na$frac_NA,breaks=100, main = "Distribution of missing values across countries",
     xlab="Fraction of missing values", ylab="Number of countries")
abline(v = 0.3, col = "red", lwd = 2, lty = 2)
legend("topright",
       legend = "Cutoff",
       col = "red",
       lwd = 2,
       lty = 2,
       bty = "n")  

#Remove countries 
naCountries <- as.matrix(na[na$frac_NA>0.3,1])
SPI_Time_Series <- SPI_Time_Series%>%
  filter(!(Country %in% naCountries))

nCountries <- length(unique(SPI_Time_Series$Country))

#Change names for Mice
names(SPI_Time_Series) <- gsub(" ", "_", names(SPI_Time_Series))
names(SPI_Time_Series) <- gsub("%", "p", names(SPI_Time_Series))
names(SPI_Time_Series) <- gsub("/", "per", names(SPI_Time_Series))
names(SPI_Time_Series) <- gsub("[\\+ | . | ( | ) | = | , | ;]", "", names(SPI_Time_Series))
names(SPI_Time_Series) <- gsub("__", "_", names(SPI_Time_Series))
names(SPI_Time_Series) <- gsub("_score_0-100", "", names(SPI_Time_Series))


#Load and prepare other sources ------------------------------------------
datasets <- list(
  PrimSchoolEnroll <- read.csv(files$PrimSchoolEnroll,skip=3)%>%
    pivot_longer(
      cols = starts_with("X"), 
      names_to = "Year",
      values_to = "Prim_School_Enroll",
      names_transform = list(Year = ~sub("^X", "", .x)))%>% #Removes X prefix
    filter((Year >= 1990) & (Year <= 2020))%>%
    select(-c(Indicator.Name, Indicator.Code,Country.Code ))%>%
    rename(Country = Country.Name,
           SPI_year = Year)%>%
    mutate(SPI_year=as.integer(SPI_year),
           Country = as.factor(Country)),
  
  safely_managed_sanitation <- read.csv(files$SafelyManagedSaniation,skip=3)%>%
    pivot_longer(
      cols = starts_with("X"), 
      names_to = "Year",
      values_to = "Safely_managed_saniation",
      names_transform = list(Year = ~sub("^X", "", .x)))%>% #Removes X prefix
    filter((Year >= 1990) & (Year <= 2020))%>%
    select(-c(Indicator.Name, Indicator.Code,Country.Code ))%>%
    rename(Country = Country.Name,
           SPI_year = Year)%>%
    mutate(SPI_year=as.integer(SPI_year),
           Country = as.factor(Country)),
  
  gdp <- read.csv(files$gdp) %>%
    select(Country = Country, SPI_year = Time, GDP) %>%
    mutate(SPI_year = as.integer(SPI_year), GDP = as.numeric(GDP)),
  
  caloric_supply <- read.csv(files$caloric_supply) %>%
    filter(Code != "", Year >= 1990 & Year <= 2020) %>%
    select(Country = Entity, SPI_year = Year, Daily_calorie_supply_pc = Daily.calorie.supply.per.person)%>%
    mutate(SPI_year = as.integer(SPI_year), Daily_calorie_supply_pc = as.numeric(Daily_calorie_supply_pc)),
  # select(all_of(c("Country", "SPI_year", "Daily_calorie_supply_per_person"))),
  
  healthy_life_expectancy <- read.csv(files$healthy_life_expectancy) %>%
    filter(Indicator == "Healthy life expectancy (HALE) at birth (years)" & Dim1 == "Both sexes") %>%
    select(Country = Location, SPI_year = Period, Healthy_Life_Expectancy = FactValueNumeric) %>%
    mutate(SPI_year = as.integer(SPI_year), Healthy_Life_Expectancy = as.numeric(Healthy_Life_Expectancy)),
  # select(all_of(c("Country", "SPI_year", "Healthy_Life_Expectancy"))),
  
  protein_supply <- read.csv(files$protein_supply) %>%
    select(Country = Country, SPI_year = Year, Protein_supply = Food.supply..Protein.g.per.capita.per.day.) %>%
    mutate(SPI_year = as.integer(SPI_year), Protein_supply = as.numeric(Protein_supply)),
  # select(all_of(c("Country", "SPI_year", "Protein_supply"))),
  
  child_mortality <- read.csv(files$child_mortality) %>%
    filter(Dim1 == "Both sexes", Period >= 1990, Period <= 2020) %>%
    select(Country = Location, SPI_year = Period, Child_Mortality = FactValueNumeric) %>%
    mutate(SPI_year = as.integer(SPI_year), Child_Mortality = as.numeric(Child_Mortality)),
  # select(all_of(c("Country", "SPI_year", "Child_Mortality"))),
  
  physicians <- read.csv(files$physicians)[-c(8247:8251), ] %>%
    select(Country = Country.Name, SPI_year = Year, Physicians_per_1000 = Physicians..per.1.000.people...SH.MED.PHYS.ZS.) %>%
    mutate(SPI_year = as.integer(SPI_year), Physicians_per_1000 = as.numeric(Physicians_per_1000)),
  # select(all_of(c("Country", "SPI_year", "Physicians_per_1000"))),
  
  hospital_beds <- read.csv(files$hospital_beds)[-c(8247:8251), ] %>%
    select(Country = Country.Name, SPI_year = Time, Hospital_Beds_per_1000 = Hospital.beds..per.1.000.people...SH.MED.BEDS.ZS.) %>%
    mutate(SPI_year = as.integer(SPI_year), Hospital_Beds_per_1000 = as.numeric(Hospital_Beds_per_1000)),
  # select(all_of(c("Country", "SPI_year", "Hospital_Beds_per_1000"))),
  
  essential_health_services <- read_excel(files$essential_health_services, sheet = "Record format") %>%
    select(Country = GeoAreaName, SPI_year = TimePeriod, Essential_Health_Coverage = Value) %>%
    mutate(SPI_year = as.integer(SPI_year), Essential_Health_Coverage = as.numeric(Essential_Health_Coverage)),
  
  basic_drinking_water <- read.csv(files$basic_drinking_water)[-c(8248:9050), ] %>%
    select(Country = Country.Name, SPI_year = Year, Basic_Drinking_Water = People.using.at.least.basic.drinking.water.services....of.population...SH.H2O.BASW.ZS.) %>%
    mutate(SPI_year = as.integer(SPI_year), Basic_Drinking_Water = as.numeric(Basic_Drinking_Water)),
  
  safely_managed_drinking_water <- read.csv(files$safely_managed_drinking_water)[-c(8248:9050), ] %>%
    select(Country = Country.Name, SPI_year = Year, Safely_Managed_Drinking_Water = People.using.safely.managed.drinking.water.services....of.population...SH.H2O.SMDW.ZS.) %>%
    mutate(SPI_year = as.integer(SPI_year), Safely_Managed_Drinking_Water = as.numeric(Safely_Managed_Drinking_Water)),
  
  basic_handwashing <- read.csv(files$basic_handwashing)[-c(8248:9050), ] %>%
    select(Country = Country.Name, SPI_year = Year, Basic_Handwashing = People.with.basic.handwashing.facilities.including.soap.and.water....of.population...SH.STA.HYGN.ZS.) %>%
    mutate(SPI_year = as.integer(SPI_year), Basic_Handwashing = as.numeric(Basic_Handwashing)),
  
  clean_fuels <- read.csv(files$clean_fuels)[-c(8248:8252), ] %>%
    select(Country = Country.Name, SPI_year = Time, Clean_Fuels = Access.to.clean.fuels.and.technologies.for.cooking....of.population...EG.CFT.ACCS.ZS.) %>%
    mutate(SPI_year = as.integer(SPI_year), Clean_Fuels = as.numeric(Clean_Fuels)),
  
  house_price_income_ratio <- read.csv(files$house_price_income_ratio) %>%
    select(Country = Reference.area, SPI_year = TIME_PERIOD, House_Price_Income_Ratio = OBS_VALUE) %>%
    mutate(SPI_year = as.integer(SPI_year), House_Price_Income_Ratio = as.numeric(House_Price_Income_Ratio)),
  
  share_slums <- read.csv(files$share_slums) %>%
    select(Country = Entity, SPI_year = Year, Share_Slums = X11.1.1...Proportion.of.urban.population.living.in.slums.......EN_LND_SLUM) %>%
    mutate(SPI_year = as.integer(SPI_year), Share_Slums = as.numeric(Share_Slums)),
  
  share_overcrowding <- read.csv(files$share_overcrowding) %>%
    select(Country =        Reference.area, SPI_year = TIME_PERIOD, Share_Overcrowding = OBS_VALUE) %>%
    mutate(SPI_year = as.integer(SPI_year), Share_Overcrowding = as.numeric(Share_Overcrowding)),
  
  intentional_homicide <- read_excel(files$intentional_homicide, sheet = "data_cts_intentional_homicide") %>%
    filter(
      Indicator == "Victims of intentional homicide",
      Dimension == "Total",
      Category == "Total",
      Sex == "Total",
      Age == "Total",
      `Unit of measurement` == "Rate per 100,000 population") %>%
    select(Country = Country, SPI_year = Year, Intentional_Homicide_Victims = VALUE) %>%
    mutate(SPI_year = as.integer(SPI_year), Intentional_Homicide_Victims = as.numeric(Intentional_Homicide_Victims))
)


#Join datasets
extendedData <- SPI_Time_Series
for (dataset in datasets) {
  extendedData <- left_join(extendedData, dataset, by = c("Country", "SPI_year"))
}



#Filter data-----------------------------------------------
naMap <- extendedData %>%
  group_by(SPI_year)%>%
  summarise_all(~sum(is.na(.)))%>%
  select(-SPI_year)
library(lattice)
levelplot(as.matrix(naMap))

#Removal due to high missingness
extendedData <- extendedData%>%
  select(-c(Share_Overcrowding,
            House_Price_Income_Ratio))

#Removal of unfit indicators
extendedData <- extendedData%>%
  select(-c(GDP, #GDP doesn't measure well being outcomes
            Basic_Handwashing, #High amounts of missingness and similar measure exists in SPI
            Basic_Drinking_Water, #Not sufficient for DLS
            Essential_Health_Coverage, #Measures UHC, which is already in SPI
            Child_Mortality)) #Already measured in SPI


naMap <- extendedData %>%
  group_by(SPI_year)%>%
  summarise_all(~sum(is.na(.)))%>%
  select(-SPI_year)

levelplot(as.matrix(naMap))

#Remove early years due to important indicators missing data
extendedData <- extendedData%>%
  filter(SPI_year>=2000)

summary(extendedData)

#write.csv(extendedData,"extendedData.csv")

