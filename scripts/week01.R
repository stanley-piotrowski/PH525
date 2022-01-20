# week01.R ----------------------------------------------------------------

# Code for week 01 

# Setup -------------------------------------------------------------------

library(pacman)
p_load(tidyverse, UsingR, patchwork)

# Download files ----------------------------------------------------------

url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleMiceWeights.csv"
download.file(url, "./data/female_mice_weights.csv")

# Mouse data exercises ----------------------------------------------------

mice_data <- read_csv("./data/female_mice_weights.csv")
mice_data[12, 2] # 26.2
mice_data$Bodyweight[11] # 26.91
length(mice_data$Diet) # 24
mean(mice_data[mice_data$Diet == "hf", 2]$Bodyweight) # 26.83417

set.seed(1)
s <- sample(13:24, 1)
mice_data$Bodyweight[s]

# dplyr exercises ---------------------------------------------------------

url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/msleep_ggplot2.csv"
download.file(url, "./data/sleep_data.csv")

sleep_data <- read_csv("./data/sleep_data.csv")
class(sleep_data)
nrow(sleep_data[sleep_data$order == "Primates", ]) # 12 -- base R way
primates <- filter(sleep_data, order == "Primates") # tidyverse way
class(primates)
primates$sleep_total
class(filter(sleep_data, order == "Primates") %>% select(sleep_total))

filter(sleep_data, order == "Primates") %>% 
  select(sleep_total) %>% 
  unlist() %>% 
  mean()

filter(sleep_data, order == "Primates") %>% 
  select(sleep_total) %>% 
  summarise()

# qq-plot exercise --------------------------------------------------------

# Main takeaway here-- when data are not normally distributed, the mean and standard deviation are not good summaries of the data
# The normal distribution has many desirable properties, including being able to describe the distribution of the data using the mean and standard deviation

# Download data
url <- "http://courses.edx.org/c4x/HarvardX/PH525.1x/asset/skew.RData"
download.file(url, "./data/skew.RData")
load("./data/skew.RData")
dim(dat) # matches what we expect from the course description

# Create series of plots to see which ones are skewed from the normal distribution
dat_long <- data.frame(dat) %>% 
  pivot_longer(everything(), names_to = "col", values_to = "val") 
  
dat_long %>% 
  ggplot(aes(sample = val)) + 
  stat_qq() + 
  geom_qq_line() + 
  facet_wrap(~ col, scales = "free")

# Columns 4 and 9 are skewed-- pull these and create histograms
filter(dat_long, col %in% c("X4", "X9")) %>% 
  ggplot(aes(val)) + 
  geom_histogram() + 
  facet_wrap(~col, scales = "free") 

# column 4 is skewed to the right; column 9 is skewed to the left

# Insect spray data -------------------------------------------------------
InsectSprays %>% 
  ggplot(aes(fct_reorder(spray, count), count)) + 
  geom_boxplot() + 
  labs(x = "Spray")

# Spray C has the lowest median count-- most effective

# UsingR ------------------------------------------------------------------

# Load the nym.2002 data--NYC marathon data
data(nym.2002, package = "UsingR")

# Use boxplots and histograms to compare finishing times of males and females
nym.2002

# Boxplots
boxplots <- nym.2002 %>% 
  ggplot(aes(fct_reorder(gender, time), time)) + 
  geom_boxplot() + 
  labs(x = "Sex")

# Histograms
histograms <- nym.2002 %>% 
  ggplot(aes(time, fill = gender)) + 
  geom_histogram() + 
  facet_wrap(~gender, scales = "free") + 
  theme(legend.position = "none")

# Wrap
wrap_plots(boxplots, histograms)

# Observations-- pretty similar, but females had a higher median finishing time relative to males
# Both groups have slightly right-tailed distributions
median(nym.2002[nym.2002$gender == "Female", 5]) # 277.725
median(nym.2002[nym.2002$gender == "Male", 5]) # 256.0167

# So, both distributions are relatively similar, but the distribution for the males is shifted to the left (faster) by about 20 minutes


