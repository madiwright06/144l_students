##### EEMB 144L Intro to R with CAL FIRE data #####
#Madison Wright
#10/2/2020

# This is comment, below this will be a command

#tidyverse is a collection of many useful packages
install.packages("tidyverse")
library(tidyverse) 

# readxl is used to read excel files
install.packages("readxl")
library(readxl)

##### Load Data #####

# unlike .csv files .xlsx files can have multiple sheets

excel_sheets("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx")

calfire.metadata <- read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx", sheet = "Metadata")
calfire.data <- read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx", sheet = 2)

#### Initial data exploration ####

names(calfire.data) #shows variable (column) names
dim(calfire.data) #dimensions of dataset (rows, columns)
class(calfire.data) #data class
head(calfire.data) #shows first 6 lines of dataset
tail(calfire.data) #shows last 6 lines of data set

# Want to know how a function works?

?names # single ? brings up R doc for that function
??names #double ?? brings up every function that might contain "names"

# single columns can be reffered to using a '$'
county <- calfire.data$County_Unit
county

max_acres <- max(calfire.data$Total_Acres_Burned, na.rm = T)
max(calfire.data$Structures_Destroyed)
max(calfire.data$Structures_Destroyed, na.rm = T)

##### Basic data wrangling (dyplr functions) #####

df1 <- select(calfire.data, County_Unit:Controlled_Date, Total_Acres_Burned:Civil_Fatalities )
view(dfl)

df2 <- filter(dfl, County_Unit %in% c("SANTA BARBARA", "VENTURA", "LOS ANGELES", "ORANGE", "SAN DIEGO") & Total_Acres_Burned >= 500 | Fire_Name== "THOMAS")

df3 <- arrange(df2, desc(Start_Date), Total_Acres_Burned)
view(df3)

df4 <- mutate_at(df3, vars(Structures_Destroyed:Civil_Fatalities), replace_na, 0)

df5 <-  mutate(df4, Fatalities = Fire_Fatalities + Civil_Fatalities)
view(df5)

#mess with time
install.packages("lubridate")
library(lubridate)

df6 <-  mutate(df5, | | | | | | interv = interval(Start_Date, Controlled_Date), + | | | | | | dur = as.duration(interv))

## We used 15 lines to do all of that! Now we have 5 different dataframes! This seems a little inefficient. There is a better way - it's called "piping"

#### Intro to piping ####

# We want to restrict our data to the SoCal coast, exlude fires that burned less than 500 acres, add column that sums the number of fatalities, change NAs to 0s, arrange data

# the MAGICAL pipe operator: %>% (command + shift + m on a mac)

# Think of this as a code way of saying "and them..."

socal_fires <- calfire.data %>% filter(County_Unit %in% c("SANTA BARBARA", "VENTURA", "LOS ANGELES", "ORANGE", "SAN DIEGO") & Total_Acres_Burned >= 500 | Fire_Name== "THOMAS") mutate_at(df3, vars(Structures_Destroyed:Civil_Fatalities), replace_na, 0) mutate(Fatalities = Fire_Fatalities + Civil_Fatalities

#### Our first graphs in ggplot ####

# Three things you must tell R to make a graph in ggplot
# (1) That you're using ggplot
# (2) What data you're using (including what should be x and what should be y)
# (3) Type of graph that you want to create



