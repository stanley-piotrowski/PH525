# week02.R ----------------------------------------------------------------

# Week 02: random variables and probability distributions

# Setup -------------------------------------------------------------------

library(pacman)
p_load(tidyverse, gapminder, rafalib, patchwork)

# Control population data -------------------------------------------------

# Download data
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
download.file(url, "./data/female_controls_population.csv")
population <- read_csv("./data/female_controls_population.csv")
female_mice <- read_csv("./data/female_mice_weights.csv")

# Questions
mean(population$Bodyweight)

set.seed(5)
avg_sample <- sample(population$Bodyweight, 5)
abs(mean(avg_sample) - mean(population$Bodyweight))

# Simulate the null hypothesis by generating random samples from the control population
# the null hypothesis is that there's no difference between the treatment and control
# Put this into a loop to do it over and over again, say 10000 times
diff_means <- vector(mode = "numeric", length = 10000L)

for (i in 1:10000) {
  
  # Sample control and treatment mice to simulate null hypothesis
  control <- unlist(sample(population$Bodyweight, 12))
  treatment <- unlist(sample(population$Bodyweight, 12))
  
  # Populate the diff means vector
  diff_means[i] <- mean(treatment) - mean(control)
  
}

# We can create a histogram of the difference in means with the null hypothesis, then plot the actual observed difference as a vertical line
treatment <- unlist(female_mice[female_mice$Diet == "hf", 2])
control <- unlist(female_mice[female_mice$Diet == "chow", 2])
actual_diff <- mean(treatment) - mean(control) # about 3

data.frame(diff_means = diff_means) %>% 
  ggplot(aes(diff_means)) + 
  geom_histogram(bins = 30) + 
  geom_vline(xintercept = actual_diff, linetype = "dashed", color = "red")

# We can see that the actual difference in means is way out in the tail of the distribution of means

# Now to calcualte the p-value
# We need to compute the number of differences in means that are as extreme or more extreme than the one we observed
length(which((diff_means >= actual_diff))) / length(diff_means) # 0.019

# Another way to do it
sum(diff_means > actual_diff) / length(diff_means)  # 0.019 -- one-sided

# To get the two-sided p-value, we can multiply the one-sided p-value by 2
0.019 * 2

# Null distribution questions ---------------------------------------------

# Sample 5 mice 1,000 times and compute the averages
# Then compute the proportions of these averages that are more than 1 gram away from the actual population average
set.seed(1)
sample_means <- vector("numeric", 1000)
for (i in 1:1000) {
  x <- sample(population$Bodyweight, 5) 
  sample_means[i] <- mean(x)
}

pop_mean <- mean(population$Bodyweight)
diffs <- sample_means - pop_mean
length(which(abs(diffs) > 1)) / 1000 # 0.503

# Repeat the exercise with 10000 samples
n <- 10000
set.seed(1)
sample_means <- vector("numeric", n)
for (i in 1:n) {
  x <- sample(population$Bodyweight, 5)
  sample_means[i] <- mean(x)
}

# Repeat this using apply
sample_means_v2 <- vector("numeric", n)
sample_means_v2 <- sapply(1:n, function(x) {
  x <- sample(population$Bodyweight, 5)
  sample_means[i] <- mean(x)
})

# And one more way with purrr-- this is much faster
sample_means_v3 <- map_dbl(1:n, function(x) {
  mean(sample(population$Bodyweight, 5))
})


pop_mean <- mean(population$Bodyweight)
diffs <- sample_means - pop_mean
length(which(abs(diffs) > 1)) / n # 0.5084

# Gapminder exercises -----------------------------------------------------

# Get life expectancy for the year 1952
data(gapminder)
x <- dplyr::filter(gapminder, year == "1952")
ggplot(aes(lifeExp), data = x) + 
  geom_histogram()

# Proportion of countries in 1952 with life Expectancy less than or equal to 40
# We can use the ecdf() function to compute the empirical cumulative distribution function-- this will return the function that can be used to compute the probability of a random variable
func_x <- ecdf(x$lifeExp)
func_x(40) # 0.2887324

# Or, using the following
sum(x$lifeExp <= 40) / length(x$lifeExp) # same answer-- 0.2887324
mean(x$lifeExp <= 40) # same answer-- 0.2887324

# Height data and the normal distribution ---------------------------------

# The normal distribution can be described using just two numbers-- the mean and the standard deviation
# So, if data can be approximated by the normal distribution, we can describe the proportion of data at any random value-- further, we know that if data are approximately normal, about 68% of observations will be within 1 standard deviation of the mean

# We've already constructed a histogram-- now let's build a QQ-plot to check for normality in the life expectancy data from 1952
ggplot(aes(sample = lifeExp), data = x) + 
  stat_qq(alpha = 0.8) + 
  stat_qq_line(linetype = "dashed") + 
  labs(x = "Sample quantiles", 
       y = "Theoretical quantiles", 
       title = "QQ plot of life expectancies from 1952")

# Looking at the QQ-plot, the data are non-normal for much of the distribution
# If the data were normally distributed, the sample and theoretical quantiles (the latter are computed from a normal distribution) should be the same, and ideally all of the points would fall on the 1:1 line

# If our data follow the normal distribution, we can also calculate the Z-scores, or standardized units, which are named because they follow a standard normal distribution with mean 0 and standard deviation of 1
# Z-scores are computed by subtracting the mean from each observation and dividing by the standard deviation-- now, we're talking about the number of standard deviations from the mean 

# For example, we can look at the mice population data for this
# First, let's check to see if the data are normally distributed using a quantile-quantile plot
population %>% 
  ggplot(aes(sample = Bodyweight)) + 
  stat_qq(alpha = 0.5) + 
  stat_qq_line()

# In general, this looks pretty good
# Let's transform the data into z-scores
standardize_z <- function(x) {
  z <- (x - mean(x)) / sd(x)
  return(z)
}

population$z_score <- standardize_z(population$Bodyweight)

# Plot the z-scores
population %>% 
  ggplot(aes(z_score)) + 
  geom_histogram() + 
  geom_vline(xintercept = mean(population$z_score), linetype = "dashed")

# Now we can see that the z-scores are centered at a mean of 0

# Normal distribution exercises -------------------------------------------

# Using the female mice population data
x <- population$Bodyweight

# Take the mean of 5 mice weights 1,000 times
set.seed(1)
sample_means_5 <- vector("numeric", 1000L)
for (i in 1:1000) {
  y <- sample(x, 5)
  sample_means_5[i] <- mean(y)
}

# Take the mean of 50 mice weights 1,000 times
set.seed(1)
sample_means_50 <- vector("numeric", 1000L)
for (i in 1:1000) {
  y <- sample(x, 50) 
  sample_means_50[i] <- mean(y)
}

# Look at the distributions for the sample averages taken from 5 mice vs 50 mice
hist_df <- data.frame(
  means = c(sample_means_5, sample_means_50), 
  n_samples = c(rep(5, 1000), rep(50, 1000))
)

hist_df %>% 
  ggplot(aes(means, fill = n_samples)) + 
  geom_histogram() + 
  facet_wrap(~n_samples) + 
  theme(legend.position = "none")

# In general, the distributions are fairly similar
mean(sample_means_5) # 23.82128
mean(sample_means_50) # 23.90436

# The means are almost the same, but the tails are quite different-- the means taken from 5 samples have much longer tails relative to means taken from 50 samples
# We can also look at QQ-plots to see this
hist_df %>% 
  ggplot(aes(sample = means)) + 
  stat_qq() + 
  stat_qq_line() + 
  facet_wrap(~n_samples, scales = "free")

# Both are normally distributed from the QQ-plots though

# What proportions of the data from the sample size of 50 are between 23 and 25?
mean(sample_means_50 < 25 & sample_means > 23) # 0.706

# Using the pnorm with average of 23.9 and sd of 0.43
pnorm(25, mean = 23.9, sd = 0.43) - pnorm(23, mean = 23.9, sd = 0.43) # 0.9765648

# Recall- random variables will change if we take another sample from the population

# Population, samples, and estimates exercises ----------------------------

# Download the mice phenotype data
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/mice_pheno.csv"
download.file(url, destfile = "./data/mice_pheno.csv")

# Read data and remove missing values
mice_pheno <- read_csv("./data/mice_pheno.csv") %>% 
  na.omit() 

# Find means of all bodyweights for each sex and diet type
mice_pheno %>% 
  group_by(Sex, Diet) %>% 
  summarise(mean_bodyweight = mean(Bodyweight))

male_controls <- filter(mice_pheno, Sex == "M", Diet == "chow")
mean(male_controls$Bodyweight)

# Compute population standard deviation with rafalib package
# The function takes the variance and multiples the result by (n-1)/n to control for sample size
rafalib::popsd(male_controls$Bodyweight) # 4.420501

# Compare to normal population standard deviation function
abs(sd(male_controls$Bodyweight) - rafalib::popsd(male_controls$Bodyweight)) # 0.009944896

# Take a random sample of 25 bodyweights from male mice on the chow diet and compute the mean
set.seed(1)
mean(sample(male_controls$Bodyweight, 25))

# Compute the population mean of all males on the high fat diet
y <- filter(mice_pheno, Sex == "M", Diet == "hf") %>% 
  pull("Bodyweight")

mean(y) # 34.84793

# Compute the population standard deviation
rafalib::popsd(y)

# Take a sample mean from 25 mice
set.seed(1)
mean(sample(y, 25))

# What's the difference in absolute value between the sample averages of high fat vs control diets in males and the population averages?
x <- filter(mice_pheno, Sex == "M", Diet == "chow") %>% 
  pull("Bodyweight")

set.seed(1)
xs <- mean(sample(x, 25))

set.seed(1)
ys <- mean(sample(y, 25))
diff_s <- ys - xs # 5.284

set.seed(1)
xp <- mean(x)

set.seed(1)
yp <- mean(y)
diff_p <- yp - xp # 3.884116

abs(diff_p - diff_s) # 1.399884

# Now do the same thing for females
set.seed(2)
x <- filter(mice_pheno, Sex == "F", Diet == "chow") %>% 
  pull("Bodyweight")

y <- filter(mice_pheno, Sex == "F", Diet == "hf") %>% 
  pull("Bodyweight")


set.seed(2)
xs <- mean(sample(x, 25))

set.seed(2)
ys <- mean(sample(y, 25))
diff_s <- diff_ys - diff_xs

abs(diff_p - diff_s)

# One important thing to recall about the central limit theorem and the standard deviation-- the spread of the distribution, or the standard deviation, is partially controlled by the sample size, so if we take larger samples, we would expect the standard deviation to be smaller-- in fact, that is exactly what we saw when we compared the histograms from sample sizes of 10 vs sample sizes of 50

# Central limit theorem exercises -----------------------------------------

# Using mouse phenotype data
y <- filter(mice_pheno, Sex == "M", Diet == "chow") %>% 
  pull("Bodyweight")

y_sd <- rafalib::popsd(y)

# Standardize to z-scores
z <- (y - mean(y)) / rafalib::popsd(y)
mean(abs(z) <= 1) # 0.6950673

mean(abs(z) <= 2) # 0.9461883
mean(abs(z) <= 3) # 0.9910314

# Central limit theorem states that the distribution of sample means will be normally distributed, even if the data the means are computed from are not normally distributed

# Take sample size of 25 from male control diet and generate 10,000 averages
n <- 10000
means <- vector("numeric", n)
set.seed(1)
for (i in 1:n) {
  x <- sample(y, 25)
  means[i] <- mean(x)
}

data.frame(means = means) %>% 
  ggplot(aes(sample = means)) + 
  stat_qq() + 
  stat_qq_line()

# Average and standard deviation of the sample distribution of sample averages
mean(means)
rafalib::popsd(means)

# Compare to the actual population
abs(mean(y) - mean(means)) # 0.004744861

# t-test construction -----------------------------------------------------

# The beauty of the central limit theorem is that if our data are well approximated by the normal distribution, we can compute a t-statistic from the sample data and compare its value to a normal distribution with mean = 0 and standard deviation of 1-- in other words, we can make inferences about the population by only using the sample data

# Now we can compute the t-test by hand using the female mice data
control <- filter(female_mice, Diet == "chow") %>% 
  pull("Bodyweight")

treatment <- filter(female_mice, Diet == "hf") %>% 
  pull("Bodyweight")

# Construct the t-statistic
obs <- mean(treatment) - mean(control)
se <- sqrt(
  (var(control) / length(control)) +
  (var(treatment) / length(treatment))
)

t_stat <- obs / se # 2.055174

# Compute the p-value assuming the normal distribution
2 * (1 - pnorm(t_stat, mean = 0, sd = 1)) # 0.0398622-- double-sided p-value

# Simulate the normal distribution, sampling 12 individuals from the population
# We're not making a distinction between control and treatment because we're assuming that under the null, they're exactly the same

N <- 12
null <- vector("numeric", 10000L)
for (i in 1:10000) {
  control_sample <- sample(population$Bodyweight, 12)
  treatment_sample <- sample(population$Bodyweight, 12)
  
  obs <- mean(treatment_sample) - mean(control_sample)
  
  se <- sqrt(
    (var(treatment) / N) + 
    (var(control_sample) / N)
  )
  
  null[i] <- obs / se
  
}

# Plot histogram and QQ plot-- this looks good
data.frame(null) %>% 
  ggplot(aes(sample = null)) + 
  stat_qq() + 
  stat_qq_line()

# Show histogram with the observed t-statistic as a dashed black line
data.frame(null) %>% 
  ggplot(aes(null)) + 
  geom_histogram() + 
  geom_vline(xintercept = t_stat, linetype = "dashed")

# Running a t-test --------------------------------------------------------

# We can also use the built-in functions for performing t-tests
t.test(treatment, control) # we get roughly the same answer, although the p-value here is slightly larger because we're using a t-distribution, not assuming we can follow a normal distribution
# The major assumption of the t-distribution is that the input data need to be normally distributed-- which we can check with QQ-plots

data.frame(bodyweight = c(treatment, control), 
           groups = c(rep("treatment", 12), rep("control", 12))) %>% 
  ggplot(aes(sample = bodyweight)) + 
  stat_qq() + 
  stat_qq_line() + 
  facet_wrap(~groups, scales = "free")

# CLT and t-distribution exericses ----------------------------------------

# Roll 100 dice and then compute the proportion of 6's rolled-- then repeat this 10,000 times
# The resulting vector is a proportion of 6's rolled per 100 dice, repeated 10,000 times
N <- 10000
prop_rolls <- vector("numeric", N)
set.seed(1)
for (i in 1:N) {
  roll <- sample(1:6, 100, replace = TRUE)
  prop <- mean(roll == 6)
  prop_rolls[i] <- prop
}

# Calculate probability and variance, then calculate Z-scores
prob <- 1/6
variance <- (prob * (1 - prob)) / 100
z_scores <- sapply(prop_rolls, function(x) {
  (x - prob) / sqrt(variance)
})

# Find the proportion of times z was larger than 2 standard deviations
mean(abs(z_scores) > 2)

# Confirm that the z-scores are normally distributed
data.frame(z_scores) %>% 
  ggplot(aes(sample = z_scores)) + 
  stat_qq() + 
  stat_qq_line()

data.frame(z_scores) %>% 
  ggplot(aes(z_scores)) + 
  geom_histogram()

# Using different probability (e.g., weighted dice) and sample sizes (number of rolls per round)
p <- c(0.5, 0.5, 0.01, 0.01)
n <- c(5, 30, 30, 100)

prop_rolls <- function(p, n, N = 10000) {
  # Create vector of proportionss
  x <- vector("numeric", N)
  sides <- 1/p # make each side of the dice equal probability
  for (i in 1:N) {
    roll <- sample(1:sides, n, replace = TRUE)
    prop <- mean(roll == 1)
    x[i] <- prop
  }
  
  # Calculate z-scores
  variance <- (p * (1-p)) / n
  z_scores <- sapply(x, function(y) {
    (y - p) / sqrt(variance)
  })
  
}

z_scores_list <- map2(p, n, ~ prop_rolls(p = .x, n = .y))
qq_plots <- pmap(
  .l = list(z_scores_list, p, n),
  .f = ~ data.frame(z_scores = ..1) %>% 
    ggplot(aes(sample = z_scores)) + 
    stat_qq() + 
    stat_qq_line() + 
    labs(title = paste0("p = ", ..2, "; n = ", ..3))
)

wrap_plots(qq_plots)

# Now let's build the histograms
hists <- pmap(
  .l = list(z_scores_list, p, n), 
  .f = ~ data.frame(z_scores = ..1) %>% 
    ggplot(aes(z_scores)) + 
    geom_histogram() + 
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") + 
    labs(x = paste0("p = ", ..2, "; n = ", ..3))
)

wrap_plots(hists)

# From both sets of plots, p = 0.5 and n = 30 is best
# The QQ-plot for this combination looks the most normal, as does the histogram

# Questions related to the sample of female mice weights
x <- filter(female_mice, Diet == "chow") %>% 
  pull("Bodyweight")

y <- filter(female_mice, Diet == "hf") %>% 
  pull("Bodyweight")

mean(x)

# Comparing sample to population
control_pop <- filter(mice_pheno, Sex == "F", Diet == "chow") %>% 
  pull("Bodyweight")

treatment_pop <- filter(mice_pheno, Sex == "F", Diet == "hf") %>% 
  pull("Bodyweight")

data.frame(bw = control_pop) %>% 
  ggplot(aes(bw)) + 
  geom_histogram() + 
  geom_vline(xintercept = c(mean(x), mean(control_pop)), 
             linetype = "dashed", color = c("red", "black"))

# Both means are almost exactly the same
popsd(control_pop)/ 12
sd(x)

# Use the CLT to estimate the probability that the estimate of the sample mean is off by more than 2 grams from the population mean
2 * (1 - pnorm(2/sd(x) * sqrt(12)))

# What is the standard error of the difference in sample means?
# It's the square root of the variance of one group + variance of the other
var_x <- var(x) / 12
var_y <- var(y) / 12
se <- sqrt(var_x + var_y)

# Construct the t-statistic by hand
t_stat <- (mean(y) - mean(x)) / se

# Find the number of degrees of freedom-- use the t.test() function
ttest_res <- t.test(y, x)
ttest_res

# Calculate p-value using the central limit theorem-- mean of 0 and sd of 1
2 * (1 - pnorm(t_stat, mean = 0, sd = 1))

# p-value under the t-distribution
ttest_res$p.value

# More CLT and t-distribution practice ------------------------------------

# Calculate difference in means and the variance
diff <- mean(treatment) - mean(control)
se <- sqrt(
  var(treatment) / length(treatment) + 
    var(control) / length(control)
)

# Compute t-statistic
t_stat <- diff / se # 2.055174
 
# Simulate a null hypothesis by drawing samples from the population data
# Build a function to vary the sample size
simulate_null <- function(s) {
  
  # Build vector of t-statistics
  n <- 10000
  null <- vector("numeric", n)
  for (i in 1:n) {
    x <- sample(population$Bodyweight, s)
    y <- sample(population$Bodyweight, s)
    
    # Calculate t-statistic
    diff <- mean(y) - mean(x) 
    se <- sqrt(
      (var(x) / length(x)) + 
        (var(y) / length(y))
    )
    
    null[i] <- diff / se
  }
  
  # PLot histogram and QQ-plot
  hist_null <- data.frame(null) %>% 
    ggplot(aes(null)) + 
    geom_histogram() + 
    geom_vline(xintercept = mean(null))
  
  qqplot_null <- data.frame(null) %>% 
    ggplot(aes(sample = null)) + 
    stat_qq() + 
    stat_qq_line(linetype = "dashed")
  
  wrap_plots(hist_null, qqplot_null) + 
    plot_annotation(title = paste0("n = ", s))
}

# Vary sample sizes from 3, 12, 50 and 100
null_plots <- map(c(3, 12, 50, 100), ~ simulate_null(.x))
wrap_plots(null_plots)

# From these plots, we can see that a sample size of 12 is reasonably approximated by the normal distribution, but it's quite a bit off using a sample size of 3
