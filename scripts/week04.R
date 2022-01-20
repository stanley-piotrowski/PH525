# week04.R ----------------------------------------------------------------

# This week's material focuses on exploratory data analysis techniques

# Setup -------------------------------------------------------------------

# Libraries
library(pacman)
p_load(tidyverse, rafalib, ggbeeswarm, UsingR, patchwork, dslabs)

# Quantile-quantile plots -------------------------------------------------

# Define qq-plot helper function
geom_qq_plot <- function(base_plot) {
  base_plot + 
    geom_qq(alpha = 0.5) + 
    geom_qq_line(linetype = "dashed") + 
    labs(x = "Theoretical quantiles", 
         y = "Sample quantiles")
}
# QQ-plot of the father son data against a standard normal distribution
ggplot(father.son, aes(sample = fheight)) %>% 
  geom_qq_plot()

# There is some slightly non-normal behavior in the tails-- however, we can see this for data sampled from a normal distribution as well using Monte Carlo simulations
data.frame(x = rnorm(1000)) %>% 
  ggplot(aes(sample = x)) %>% 
  geom_qq_plot()

# Now we can simulate drawing random variables from a t-distribution with different degrees of freedom
# We'll look at the differences between the histograms and the QQ-plots to see how the proportion of the data in the tails changes with the degrees of freedom
sim_df <- map(c(3, 6, 12, 30), ~ rt(1000, .x)) %>% 
  set_names(paste0("df_", c(3, 6, 12, 30))) %>% 
  bind_rows(.id = "df") %>% 
  pivot_longer(everything(), names_to = "df", values_to = "random_var") %>% 
  mutate(df = factor(gsub("df_", "df = ", df), 
                     levels = paste0("df = ", c(3, 6, 12, 30))))

# Build histograms and QQ-plots against standard normal distribution
histograms <- ggplot(sim_df, aes(random_var, fill = df)) + 
  geom_histogram() + 
  facet_wrap(~df, scales = "free", nrow = 1) + 
  labs(x = "Random variable", y = "Count") +
  theme(legend.position = "none")

qqplots <- ggplot(sim_df, aes(sample = random_var, color = df)) %>% 
  geom_qq_plot() + 
  facet_wrap(~df, scales = "free", nrow = 1) + 
  theme(legend.position = "none")

# Wrap together to see the differences in tail behavior with degrees of freedom
wrap_plots(histograms, qqplots, nrow = 2)

# Looking at the plots above, we can see that with smaller degrees of freedom, a larger proportion of the data lie in the tails of the distribution-- in the QQ-plots specifically, we can see that this deviates substantially from the theoretical expectations
# The degrees of freedom are the number of independent observations in the dataset-- this is typically the sample size minus 1
# For example, if we compute a statistic like the mean with a sample size of 8, the degrees of freedom is 7, because we only need to know the identity of 7 of those values-- if we know the 7 values and the mean, we can compute the value of the 8th, so there are 7 independent values that are free to vary
# Extending this to the shapes of the distributions, as the degrees of freedom, or number of independent observations increases, the distribution becomes better approximated by the standard normal distribution (hence, why the sample quantiles track the line of the theoretical quantiles in the QQ-plot)

# Boxplots ----------------------------------------------------------------

# Using the executive pay dataset from the UsingR package, which tabulate the compensation for 199 US CEOs in 2000
# These data are clearly not well approximated by the normal distribution
exec_pay <- data.frame(pay = exec.pay)
exec_pay_hist <- ggplot(exec_pay, aes(pay)) + 
  geom_histogram()

exec_pay_qq <- ggplot(exec_pay, aes(sample = pay)) %>% 
  geom_qq_plot()

# Include a boxplot as well
exec_pay_boxplot <- ggplot(exec_pay, aes(y = pay + 0.5)) + 
  geom_boxplot() + 
  scale_x_discrete() + 
  scale_y_log10() + 
  labs(y = "log10 (pay)")

exec_pay_hist + exec_pay_qq + exec_pay_boxplot

# Scatterplots and correlation --------------------------------------------

# Scatterplot of the height of fathers and sons to show correlation
ggplot(father.son, aes(fheight, sheight)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(se = FALSE) + 
  labs(x = "Father's height (inches)", 
       y = "Son's height (inches)")

# What's the correlation?
cor(father.son$fheight, father.son$sheight) # 0.5

 # Stratification ----------------------------------------------------------

# Look at the height of sons where father's heights are about 72 inches
# First, round everything in the data
father_son <- father.son %>% 
  mutate(across(everything(), ~ round(.x)))

# Look at the distribution of son's heights
ggplot(father_son, aes(fheight, sheight, group = fheight)) + 
  geom_boxplot() + 
  labs(x = "Father's height (inches)", 
       y = "Son's height (inches)")

# Scatterplot exercises ---------------------------------------------------

# Using NYC marathon data from 2002
males <- filter(nym.2002, gender == "Male")
females <- filter(nym.2002, gender == "Female")

# Correlation between age and time in males and females
cor(males$age, males$time) # 0.2432273
cor(females$age, females$time) # 0.2443156

# Look at scatterplots and boxplots stratified by age, binned into 5 years
range(nym.2002$age) # 5-81

# Scatterplots for times and age, by gender
age_time_scatter <- ggplot(nym.2002, aes(age, time, color = gender)) + 
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE) + 
  facet_wrap(~gender, scales = "free") + 
  theme(legend.position = "none")

# Make 17 groups-- 5-85 in groups of 5
age <- seq(5, 85, 1)
group <- unlist(map(seq(1, 16, 1), ~ rep(.x, 5)))

# Make lookup table-- note, add one group_17 to deal with the extra age
group_lookup <- data.frame(
  age = age,
  group = c(group, 17)
)

# Join with marathon data
nyc_marathon <- nym.2002 %>% 
  left_join(group_lookup, by = "age")

# Look at boxplots stratified by age and split by gender
age_time_box <- ggplot(nyc_marathon, aes(group, time, group = group, fill = gender)) + 
  geom_boxplot() + 
  facet_wrap(~gender, scales = "free") +
  scale_x_continuous(breaks = seq(1, 17, 1), 
                     labels = seq(1, 17, 1)) + 
  theme(legend.position = "none")

wrap_plots(age_time_scatter, age_time_box, nrow = 2)

# The course answer for question 3 was that the finish times are constant through 50-60, and then they increase-- I suppose that's a reasonable conclusion for each gender, although there are clearly differences in median times (see boxplots)

# Bi-variate normal distribution ------------------------------------------

# The basic idea of a bivariate normal data is that when looking at paired data (e.g., X and Y), if we take a random value of x, all pairs of X and Y with have Y following a normal distribution
# Let's explore this with the father and son height data-- using the rounded data frame we computed earlier
father_son %>% 
  ggplot(aes(sample = sheight)) %>%  
  geom_qq_plot() +
  facet_wrap(~fheight, scales = "free")

# We can see that most of the son's heights are normally distributed for each value of father's heights, except for in the tails of the distribution
# Let's also explore this with histograms
father_son %>% 
  ggplot(aes(sheight)) +
  geom_histogram(bins = 20) +
  facet_wrap(~fheight, scales = "free")

# Here, we see the same sort of behavior

# Robust summaries --------------------------------------------------------

# Symmetry of log ratios
# Simulate ratios of two random numbers
(x <- 2^(rnorm(100)))
(y <- 2^(rnorm(100)))
(ratios <- x / y)

# Create histogram of the ratios
data.frame(ratios) %>% 
  ggplot(aes(log2(ratios))) + 
  geom_histogram(bins = 30) + 
  geom_vline()

# Taking the log of the ratios means that they will be symmetrical around 1
# Simulate 1/32 to 32, then plot on the original scale and the log2 scale--reproduce the figure from the book
nums <- c(1/32, 1/16, 1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)
data.frame(
  nums = c(nums, log2(nums)), 
  type = c(rep("normal_scale", 11), rep("log_scale", 11)), 
  horiz_line = c(rep(1, 11), rep(log2(1), 11))
) %>% 
  ggplot(aes(nums, nums)) + 
  geom_point() + 
  geom_hline(aes(yintercept = horiz_line), linetype = "dashed") +
  facet_wrap(~type, scales = "free") + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  labs(x = " ", y = "Numbers")

# The plot re-created from the textbook shows that numbers on the log scale are symmetrical about log(1)-- the same number of observations are < 1 and > 1

# Symmetry of log ratios exercises ----------------------------------------

# Using the NYC marathon data
time <- sort(nym.2002$time)

# Fastest time divided by the median time
time[1] / median(time)

# Slowest time divided by the median time
time[length(time)] / median(time)

# Plots to avoid exercises ------------------------------------------------

# Using the divore_margarine data set from the dslabs package-- data on the US consumption of margarine and the divorces per 1,000 people in Maine from 2000 to 2009
# Plot margarine consumption and the divorce rate
ggplot(divorce_margarine, aes(margarine_consumption_per_capita, divorce_rate_maine)) +
  geom_point()

cor(divorce_margarine$divorce_rate_maine, divorce_margarine$margarine_consumption_per_capita) # 0.99 correlation-- correlation does not mean causation

# Robust summaries exercises ----------------------------------------------

# Using the ChickWeight data set-- weight of chicks (in grams) from day 0 to 21 varying by diet (coded as 1-4)
# Scatterplot
ChickWeight %>% 
  ggplot(aes(Time, weight, color = Diet)) + 
  geom_point()

# Pivot data to wider format-- use tidyr, not reshape
# Pivot so that the chick and diet fields are the IDs and the times are now in wide format
chick <- pivot_wider(
  ChickWeight, 
  id_cols = c("Chick", "Diet"), 
  names_from = "Time", 
  values_from = "weight", 
  names_prefix = "day_") %>% 
  na.omit()

# Compute the average weight of day 4 chicks with an outlier of 3000 grams, divided by the weight of the day 4 chicks without the outlier
mean(c(chick$day_4, 3000)) / mean(chick$day_4) # 2.062407

# Compute the same ratio with the median
median(c(chick$day_4, 3000)) / median(chick$day_4) # 1

# Compute the same ratio with the standard deviation
sd(c(chick$day_4, 3000)) / sd(chick$day_4) # 101.2859

# Compare the ratio with the sd() to the mad()
mad(c(chick$day_4, 3000)) / mad(c(chick$day_4)) # 1

# Compare the Pearson and Spearman correlations-- the former is influenced the presence of the outlier
# The scatterplot is colored by the chick ID-- in general, chicks with relatively low weight on day 4 also have relatively low weight on day 21
ggplot(chick, aes(day_4, day_21, color = Chick)) + 
  geom_point() + 
  theme(legend.position = "none")

# Compute the correlations with different methods
cor_outlier <- cor(c(chick$day_4, 3000), c(chick$day_21, 3000))
cor_normal <- cor(chick$day_4, chick$day_21)
cor_outlier / cor_normal

# Wilcoxon rank-sum test --------------------------------------------------

# Both the mean and sd are influenced by outliers, and since the t-test relies on the normal distribution controlled by these two parameters, outliers can seriously influence results
# We'll compare the results of a t-test and Wilcoxon rank sum test for two sets of numbers-- in one set, we'll use everything as is, where data are normal and the null hypothesis is true
# In the second set, we'll change two numbers in one set such that they are outliers
set.seed(779) # same as the book
x <- rnorm(25)
y <- rnorm(25)  

# Perform t-test and Wilcoxon rank sum test with the original data
t.test(x, y)$p.value # 0.159
wilcox.test(x, y)$p.value # 0.306

# Add outliers and redo tests
x[1] <- 5
x[2] <- 7
t.test(x, y)$p.value # 0.044
wilcox.test(x, y)$p.value # 0.131 -- still not significant

# Plot the original data and the ranks
rank_x <- rank(x, )
rank_y <- rank(y)

plot_input <- data.frame(value = c(x, rank_x, y, rank_y), 
                         type = c(rep("num", 25), rep("rank", 25), rep("num", 25), rep("rank", 25)),
                         group = c(rep(1, 50), rep(2, 50)))

# This is the basic idea-- although may need to change this
plot_input %>% 
  ggplot(aes(group, value)) + 
  geom_point() + 
  facet_wrap(~type, scales = "free")

# Wilcoxon rank-sum test exercises ----------------------------------------

# Using chick weight data
x <- filter(chick, Diet == 1) %>% pull("day_4")
y <- filter(chick, Diet == 4) %>% pull("day_4")

# t-test and Wilcoxon on regular data
t.test(x, y)$p.value # 7.320259e-06
wilcox.test(x, y) # 0.0002012

# Perform tests again after adding an outlier
t.test(c(x, 200), y)$p.value # 0.9380347
wilcox.test(c(x, 200), y)$p.value # 0.0009840921

# Generate three boxplots with different outliers added
df1 <- data.frame(weight = c(x, y), diet = c(rep(1, length(x)), rep(4, length(y))), type = "y")
df2 <- data.frame(weight = c(x, y + 10), diet = c(rep(1, length(x)), rep(4, length(y))), type = "y + 10")
df3 <- data.frame(weight = c(x, y + 100), diet = c(rep(1, length(x)), rep(4, length(y))), type = "y + 100")
master_df <- bind_rows(df1, df2, df3)

master_df %>% 
  mutate(diet = factor(diet)) %>% 
  ggplot(aes(diet, weight)) + 
  geom_boxplot() + 
  facet_wrap(~type, scales = "free")

# Look at the difference in the t-statistic
(t.test(x, y + 10)$statistic) - (t.test(x, y + 100)$statistic)