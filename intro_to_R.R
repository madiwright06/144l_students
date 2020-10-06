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
