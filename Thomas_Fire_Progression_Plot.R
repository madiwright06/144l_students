#### Thomas Fire Progression Project ####
#Madison Wright
#10/23/20

##### Load Data #####

library(tidyverse)
library(readxl)

excel_sheets("Input_Data/week1/Thomas_Fire_Progression.xlsx")


data <- read_excel("Input_Data/week1/Thomas_Fire_Progression.xlsx", sheet = "Data")

metadata <- read_excel("Input_Data/week1/Thomas_Fire_Progression.xlsx", sheet = 2)

#### Data Exploration ####
names(data)
names(metadata)
dim(data)
dim(metadata)

#### data wrangling ####

df1 <- select(data, Date, Acres_Burned, Containment)
unique(df1$Acres_Burned)
df4 <- filter(df1, Acres_Burned %in% c("50000", "230000", "256000", "272600", "281893"))

df2 <- select(data, Acres_Burned, PM10)

df3 <- select(data, Acres_Burned, PM25)

#### Piping ####

library(lubridate)

thomas.fires <- data %>%
  select(Date, Acres_Burned, Containment) %>% 
  arrange(desc(Acres_Burned))

##### ggplot #####

ggplot(thomas.fires, aes(x = Date, y = Acres_Burned)) +
  geom_point(aes(color = Containment)) +
  ggtitle("2017 - 2018 Thomas Fire Progression") +
  labs(x = "Date", y = "Total Acres Burned") +
  theme_dark()

ggplot(df2, aes(x = Acres_Burned, y = PM25)) +
  geom_point() +
  ggtitle ("Acres Burned Associated with Air Quality") +
  labs(x = "Amount of Acres Burned", y = "Particulate Matter Smaller Than 10 Microns In Diameter")
  