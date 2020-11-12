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

progression.plot <- thomas.fires %>% 
  ggplot(aes(x = Date, y = Acres_Burned)) +
  geom_point(aes(color = Containment)) +
  ggtitle("2017 - 2018 Thomas Fire Progression") +
  labs(x = "Date", y = "Total Acres Burned") +
  theme_dark()


air_quality.plot <- df2 %>% 
  ggplot(aes(x = Acres_Burned, y = PM25)) +
  geom_point() +
  ggtitle ("Acres Burned Associated with Air Quality") +
  labs(x = "Amount of Acres Burned", y = "Particulate Matter Smaller Than 10 Microns In Diameter")

#### Save Data and Plots ####

saveRDS(thomas.fires, file = "Output_Data/week 1/Thomas_Fires_Plots.rds")
write_csv(thomas.fires, "Output_Data/week 1/Thomas_Fires_Plots.csv")

ggsave(filename = "Thomas_Fire_Data", progression.plot, device = "jpeg", "Output_Data/week 1/")
  
