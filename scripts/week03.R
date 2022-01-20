# week03.R ----------------------------------------------------------------

# This week's materials are focused around statistical inference

# Setup -------------------------------------------------------------------

# Libraries
library(pacman)
p_load(tidyverse, rafalib, ggbeeswarm)

# Read input files 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
download.file(url, "./data/babies.txt")

url <- "https://courses.edx.org/assets/courseware/v1/f3b6df96b94a01d80e35e3cecf3d83f0/asset-v1:HarvardX+PH525.1x+3T2020+type@asset+block/assoctest.csv"
download.file(url, "./data/association_test.csv")

babies <- read_table("./data/babies.txt")
mice_pheno <- read_csv("./data/mice_pheno.csv", show_col_types = FALSE)
female_weights <- read_csv("./data/female_mice_weights.csv", show_col_types = FALSE)
association_test <- read_csv("./data/association_test.csv", show_col_types = FALSE)

# T-test exercises --------------------------------------------------------

# Split babies data set to get the weight of babies from non-smokers and smokers
bw_smokers <- filter(babies, smoke == "1") %>% 
  pull("bwt")

bw_nonsmokers <- filter(babies, smoke == "0") %>% 
  pull("bwt")

# Get the difference in means and the population standard deviations for each
diff <- mean(bw_nonsmokers) - mean(bw_smokers) # about 8.9
popsd(bw_nonsmokers) # 17.3
popsd(bw_smokers) # 18

# We're interested in whether or not the bodyweights of babies born from non-smokers are significantly different than those from smokers

# Take sample of 25 bodyweights from each group
set.seed(1)
dat_ns <- sample(bw_nonsmokers, 25)
dat_s <- sample(bw_smokers, 25)

# Compute the t-statistic
diff <- mean(dat_ns) - mean(dat_s)
se <- sqrt(
  (var(dat_ns) / length(dat_ns)) + 
    (var(dat_s) / length(dat_s))
)

tval <- diff / se # 1.659325

# Compute the p-value by comparing the t-statistic to a standard normal distribution
2 * (1 - pnorm(tval)) # 0.09705034

# Recall that when applying the CLT, the standard deviation of the random variable will be equal to the population standard deviation divided by the square root of the sample size
# So, if we subtract the mean of the population from the mean of the sample, all divided by the standard error, that quantity will follow a standard normal distribution centered at mean 0 and standard deviation 1-- additionally, we know that about 95% of observations will fall within 2 standard deviations from the mean

# Simulating confidence intervals -----------------------------------------

# The confidence intervals describe the range of credible values a test statistic could take 95% of the time, with the true population parameter being somewhere in that range
# For example, if we wanted to know the weight of babies from non-smokers, we'd first need to take a sample 
samp_ns <- sample(bw_nonsmokers, 25)
mean_samp <- mean(samp_ns)
se_ns <- sd(samp_ns) / sqrt(length(samp_ns))
q <- qnorm(1 - 0.05/2)
ll <- mean_samp - q * se
ul <- mean_samp + q * se
c(ll, ul) # 95% CI-- 116 to 136

# Now we can do this many more times, then plot this on top of the true population mean-- which we'll assume we know
# Build a function to simulate the confidence intervals

sim_conf_int <- function(n_intervals, n_samples, distribution) {
  
  distribution <- match.arg(
    distribution, 
    choices = c("normal", "t")
  )
  
  # Initialize an empty data frame
  df <- data.frame(
    sample_means = vector("numeric", n_intervals),
    sample_se = vector("numeric", n_intervals),
    q_vec = vector("numeric", n_intervals),
    ll_vec = vector("numeric", n_intervals),
    ul_vec = vector("numeric", n_intervals),
    sim = 1:n_intervals, 
    width = vector("numeric", n_intervals), 
    within_ci = vector("logical", n_intervals)
  )
  
  for (i in 1:n_intervals) {
    
    s <- sample(bw_nonsmokers, n_samples)
    se <- sd(s) / sqrt(length(s))
    
    # Ideally, if sample size is <= 5, use t-distribution instead of normal
    # If using the t-distribution, the degrees of freedom are equal to the sample size - 2
    if (distribution == "t") {
      q <- qt(1 - 0.05/2, df = n_samples - 2)
    } else if (distribution == "normal") {
      q <- qnorm(1 - 0.05/2)
    }
    
    df$sample_means[i] <- mean(s)
    df$sample_se[i] <- se
    df$q_vec[i] <- q
    df$ll_vec[i] <- mean(s) - q * se
    df$ul_vec[i] <- mean(s) + q * se
    df$width[i] <- df$ul_vec[i] - df$ll_vec[i] # add interval width
    df$within_ci[i] <- between(mean(bw_nonsmokers), df$ll_vec[i], df$ul_vec[i])
  }
  
  return(df)
  
}

# Call with 100 simulated confidence intervals with sample size == 25
df <- sim_conf_int(100, 25, "normal") 

# Color lines based on whether or not the CI overlaps the population parameter
# Define plotting function
plot_sim_ci <- function(df) {
  
  ci_plot <- df %>% 
    ggplot(aes(sample_means, sim, group = sim, color = within_ci)) + 
    geom_errorbar(aes(xmin = ll_vec, xmax = ul_vec), alpha = 0.75) +
    geom_vline(xintercept = mean(bw_nonsmokers), linetype = "dashed") +
    labs(x = "Sample means", y = "Simulated confidence intervals", 
         title = "Distribution of confidence intervals from sampled data")
  
  return(ci_plot)
  
}

plot_sim_ci(df)

# Now we should compute the proportion of the intervals that contain the population parameter
# This should be about 95%
mean(map2_chr(df$ll_vec, df$ul_vec, ~ between(mean(bw_nonsmokers), .x, .y)) == TRUE) # 0.96

# Now repeat with 1,000 confidence intervals and compute the proportion that contain the population parameter
df2 <- sim_conf_int(1000, 25, "normal")
plot_sim_ci(df2)

mean(map2_chr(df2$ll_vec, df2$ul_vec, ~ between(mean(bw_nonsmokers), .x, .y)) == TRUE) # 0.945

# We can also use these two functions to show that confidence intervals don't perform all that well due to the impact of the standard error in the CLT with low sample sizes (n = 3)
df_small_n <- sim_conf_int(100, 3, "normal")
mean(map2_chr(df_small_n$ll_vec, df_small_n$ul_vec, ~ between(mean(bw_nonsmokers), .x, .y)) == TRUE) # 0.84

plot_sim_ci(df_small_n)

# We can see that using a smaller sample size, fewer confidence intervals overlap the true population parameter
# We can also see that the widths of the confidence intervals appear to be much wider-- let's look at the distribution of confidence interval widths for a few different scenarios-- sampling 3, 25, 50, and 100 mice from the bw non-smokers
sim_list <- map(c(3, 25, 50, 100), ~ sim_conf_int(1000, .x, "normal"))
names(sim_list) <- paste0("n_", c(3, 25, 50, 100))

bind_rows(sim_list, .id = "sample_size") %>% 
  ggplot(aes(width, fill = sample_size)) + 
  geom_histogram() + 
  facet_wrap(~sample_size, scales = "free") + 
  theme(legend.position = "none") +
  labs(x = "Confidence interval width", 
       title = "Distribution of confidence interval widths with sample size")

# These plots show that as the sample size increases, the distribution of the confidence intervals become more narrow
# Specifically, looking at the distribution with n = 3, the histogram is right-skewed

# Now let's look at how the confidence intervals change with small sample sizes if we use a t-distribution with fatter tails
df4 <- sim_conf_int(100, 3, "t")
plot_sim_ci(df4)
mean(map2_chr(df4$ll_vec, df4$ul_vec, ~ between(mean(bw_nonsmokers), .x, .y)) == TRUE) # 0.99

# Now, the tails are much wider, so a larger proportion of intervals contain the true population parameter

# The test statistics are compared to either a standard normal distribution (mean = 0; sd = 1) or a t-distribution (controlled by the mean and the degrees of freedom)
# We know that with a standard normal distribution, 95% of observations are within 2 standard deviations from the mean-- thus, if the t-statistic is greater than absolute value of 2, then it must be in the remaining 5% of the distribution (corresponding to p-value < 0.05)

# Using t.test function ---------------------------------------------------

# Now we'll look at the difference in means between treatment and control using the t.test() function
t.test(bw_smokers, bw_nonsmokers)

# Confidence interval exercises -------------------------------------------

# Use the babies bw vectors 
set.seed(1)
dat_ns <- sample(bw_nonsmokers, 25)
dat_s <- sample(bw_smokers, 25)

# What do we add and subtract to obtain a 99% confidence interval with 2*N-2 degrees of freedom?
# We need the quantile, as well as the standard error for both vectors
qt(1 - 0.01/2, 48) * sqrt(
  var(dat_ns) / 25 + 
    var(dat_s) / 25
)

# Power analyses are often used to determine the minimum sample size required to detect a result, if one exists
# Power is calculated as 1 - type II error rate (rate of false negatives)-- probability that you'll reject the null when the alternative is true
set.seed(1)
dat_ns <- sample(bw_nonsmokers, 5)
dat_s <- sample(bw_smokers, 5)
t.test(dat_ns, dat_s)

# Power calculations tutorial ---------------------------------------------

# Using the female mice phenotype data
control_pop <- filter(mice_pheno, Diet == "chow", Sex == "F") %>% pull("Bodyweight")
hf_pop <- filter(mice_pheno, Diet == "hf", Sex == "F") %>% pull("Bodyweight")

# Look at the actual difference between populations
abs(mean(control) - mean(hf)) # 2.375517

# Get percent increase 
(mean(hf) - mean(control)) / mean(control) * 100 # 10% increase with hf diet

# Sample 5 mice and perform the t-test
set.seed(1)
n <- 5
hf_sample <- sample(hf, n)
control_sample <- sample(control, n)
t.test(hf_sample, control_sample) # not significant

# Just because the p-value isn't < 0.05, even though we know there's an actual difference, simply means we don't have enough statistical power 
# Importantly, statistical power depends on the actual difference in the parameter (i.e., the true effect size), sample sizes, and population standard deviations

# Conduct the simulations for power
reject <- function(n, alpha = 0.05) {
  hf <- sample(hf_pop, n)
  control <- sample(control_pop, n)
  
  pval <- t.test(hf, control)$p.value
  return(pval < alpha)
}

# Conduct 2000 simulations
rejections <- replicate(2000, reject(n = 12))
mean(rejections) # 21% of the time we rejected the null

# We can see there's an issue with the sample size
# We can also simulate different sample sizes and see how the power to correctly reject the null increases
rejections <- map_dbl(seq(5, 50, 5), ~ mean(replicate(2000, reject(n = .x))))
data.frame(n = seq(5, 50, 5), power = rejections) %>% 
  ggplot(aes(n, power)) + 
  geom_point() + 
  geom_line() + 
  labs(
    x = "Sample size", 
    y = "Statistical power"
  )

# We can also vary alpha as well
ns <- seq(5, 50, 5)
alphas <- c(0.25, 0.1, 0.05, 0.01, 0.001, 0.001)
combos <- cross2(ns, alphas)
names(combos) <- paste0("combo_", 1:length(combos))

# Simulate rejections
rejection_sims <- map(combos, ~ replicate(2000, reject(.x[[1]], .x[[2]])))

# Plot power for all combinations
power_df <- data.frame(
  n = map_dbl(combos, ~ pluck(.x, 1)), 
  alpha = map_dbl(combos, ~ pluck(.x, 2)), 
  power = map_dbl(rejection_sims, mean)
)

power_df %>% 
  ggplot(aes(n, power, color = alpha)) + 
  geom_line() + 
  geom_point() + 
  facet_grid(~alpha) + 
  theme(legend.position = "none") + 
  labs(
    x = "Sample size",
    y = "Power", 
    title = "Power analysis with simulated data and different levels of alpha"
  )

# From these plots, we can see that when we are very strict with alpha, or type I error rate, at 0.001, we need our sample size to be about 50 before we can reject the null hypothesis about 25% of the time
# At the other extreme end of the spectrum, when we use a type I error rate of 0.25, even with only 10 samples, we rejected the null a little over 50% of the time

# Another important consideration from these simulations-- as you increase alpha, or the type I error rate, in general, holding all other factors constant, power also increases
# It's important to note that power is simply the average proportion of rejections-- it doesn't necessarily say anything about how well the type I errors are controlled
# As another exercise, say we wanted to know the sample size required to reach 90% power and control the type I error rate at 0.001--
error_control_rejections <- map(seq(20, 200, 20), ~ replicate(2000, reject(.x, alpha = 0.001)))
error_control_df <- data.frame(
  sample_size = seq(20, 200, 20), 
  power = map_dbl(error_control_rejections, mean)
)

error_control_df %>% 
  ggplot(aes(sample_size, power)) + 
  geom_line() + 
  geom_point() + 
  geom_hline(yintercept = 0.9, linetype = "dashed") + 
  labs(
    x = "Sample size", 
    y = "Power", 
    title = "Power analysis to find minimum sample size control alpha at 0.001"
  )

# This simulation shows that in order to achieve about 90% statistical power and control the type I error rate at about 0.001, you'd need to have a minimum of about 120 samples-- 90% power is shown as the dashed horizontal line

# P-values are arbitrary when the null is false-- we can simulate drawing varying sample sizes
calc_pvals <- function(n) {
  
  hf <- sample(hf_pop, n)
  control <- sample(control_pop, n)
  pval <- t.test(hf, control)$p.value
  
  return(pval)
}

# For each value of n, repeat 100 times to get a distribution of p-values
ns <- seq(10, 200, 10)
pval_list <- map(ns, ~ replicate(100, calc_pvals(.x)))
names(pval_list) <- paste0("n_", ns)

# Plot
data.frame(pval_list) %>% 
  tidyr::pivot_longer(everything(), names_to = "sample_size", values_to = "pval") %>% 
  mutate(sample_size = gsub("n_", "", sample_size)) %>% 
  ggplot(aes(fct_reorder(sample_size, as.numeric(sample_size)), log10(pval), fill = sample_size)) + 
  geom_violin(draw_quantiles = 0.5) +
  geom_hline(yintercept = c(log10(0.01), log10(0.05)), linetype = "dashed") +
  theme(legend.position = "none") + 
  labs(
    x = "Sample size", 
    y = "log10(p-value)"
  )

# From this plot, we can see that as the sample size increases, the p-values become smaller and smaller while holding alpha constant (note the y-axis is on the log10 scale)-- this is important, because it shows how the sample size influences the standard error 
# Additionally, p-values of 0.01 (top) and 0.05 (bottom), both on the log10 scale, are depicted as dashed horizontal lines to show where the simulated p-values fall in the distribution for each sample size

# We can also report effect sizes as percent increases, in addition to reporting the raw confidence intervals
# For example, simulating 12 samples drawn from each mouse population
calc_perc_increase <- function(n) {
  
  hf <- sample(hf_pop, n)
  control <- sample(control_pop, n)
  diff <- mean(hf) - mean(control)
  perc_increase <- diff / mean(control) * 100
  
  return(perc_increase)

}

# Population increase
diff_pop <- mean(hf_pop) - mean(control_pop)
perc_increase_pop <- diff_pop / mean(control_pop) * 100 # 9.942157 actual percent increase

# Conduct simulations and plot the distribution of percent increases
df <- data.frame(perc_increase = replicate(5000, calc_perc_increase(n = 12)))
df %>% 
  ggplot(aes(perc_increase)) + 
  geom_histogram() + 
  geom_vline(xintercept = mean(df$perc_increase), linetype = "dashed")

# Check normality
df %>% 
  ggplot(aes(sample = perc_increase)) + 
  stat_qq() + 
  stat_qq_line(color = "red", linetype = "dashed") + 
  labs(x = "Theoretical quantiles", y = "Sample quantiles")

# Shapiro-Wilk test for normality -----------------------------------------

# We can use the Shapiro-Wilk test for normality-- let's look at an example with the percent increases 
# The null hypothesis is that the data are normally distributed-- the alternative is the data are NOT normally distributed
shapiro.test(df$perc_increase) # p-value < 0.05

# This indicates the data are not normally distributed, which we can also see from the QQ-plot and histogram above in some areas of the tails-- but overall, it looks okay

# Power analysis exercises ------------------------------------------------

# Look at the true effect size
diff <- mean(bw_nonsmokers) - mean(bw_smokers) # 8.93

# Using the babies data
set.seed(1)
bw_ns <- sample(bw_nonsmokers, 5)
bw_s <- sample(bw_smokers, 5)
t.test(bw_ns, bw_s)$p.value # 0.1843195-- type II error, because we know the true effect size is non-zero

# Simulations
calc_rejects <- function(sample_size, alpha = 0.05, seed = 1) {
  
  set.seed(seed)
  
  rejects <- replicate(n = 10000, {
    bw_ns <- sample(bw_nonsmokers, sample_size)
    bw_s <- sample(bw_smokers, sample_size)
    
    t.test(bw_ns, bw_s)$p.value < alpha
  })
  
  return(mean(rejects))

}

calc_rejects(5, 0.05, 1) # 0.06

# Simulate with different sample sizes and find what sample size gets you 80% power
map_dbl(seq(30, 120, 30), ~ calc_rejects(.x)) # 60 samples 

# Do the same thing, now change alpha to 0.01
map_dbl(seq(30, 120, 30), ~ calc_rejects(.x, alpha = 0.01)) # 90 samples

# Power == 1 - type II error

# Monte Carlo simulations -------------------------------------------------

# Using the mice phenotype data
control <- filter(mice_pheno, Sex == "F" & Diet == "chow") %>%
  pull("Bodyweight")

# First, build a function to generate a t-statistic under the null hypothesis
tstat_generator <- function(x, n) {
  
  # Draw two samples to form the null
  cases <- sample(x, n)
  controls <- sample(x, n)
  
  # Generate t-stat
  tstat <- (mean(cases) - mean(controls)) / 
    sqrt(
      var(cases) / n + 
      var(controls) / n
    )
  
  return(tstat)
}

# Simulate t-statistic with 1,000 monte carlo simulations-- here, we expect that the distribution of test statistics should be centered around a mean of zero
# This is a Monte Carlo simulation
ttests <- replicate(1000, tstat_generator(10))
data.frame(ttests) %>% 
  ggplot(aes(ttests)) + 
  geom_histogram()

# Generating population data with parametric Monte Carlo simulation
# If we didn't have access to the population data as we do here, but we knew the mean and standard deviation of a sample and assumed the population that sample was drawn from was normally-distributed, we could generate population data
mean(control) # 24
popsd(control) # 3.5
sim_pop_control <- rnorm(
  n = 5000, 
  mean = mean(control), 
  sd = popsd(control)
)

# Re-write the function to generate random variables from a normal distribution
# This is the parametric simulation-- taking sample data and then using a model to simulate the population data
# Again, this is assuming that the data are normally distributed-- we can check a QQ-plot
tstat_generator2 <- function(n, mean, sd) {
  
  cases <- rnorm(n, mean, sd)
  controls <- rnorm(n, mean, sd)
  
  tstat <- (mean(cases) - mean(controls)) / 
    sqrt(
      var(cases) / n + 
        var(controls) / n
    )
  
  return(tstat)
}

tstats <- replicate(1000, tstat_generator2(3, mean(control), popsd(control)))

data.frame(tstats) %>% 
  ggplot(aes(sample = tstats)) + 
  geom_qq() + 
  geom_qq_line()

# The above plot is not very well approximated by the normal distribution-- let's check the t-distribution

# This section needs to be fixed
tstats_df <- data.frame(tstats)  
ps <- (seq(0, 999) + 0.05) / 1000
params <- qt(ps, 2 *3-2)
ggplot(tstats_df, aes(sample = tstats)) +
  geom_qq(distribution = qt, dparams = params) + 
  geom_qq_line(distribution = qt, dparams = params)

# Monte Carlo exercises ---------------------------------------------------

# Generate a t-statistic for a series of random variables drawn from the standard normal dist
set.seed(1)
x <- rnorm(5)
tstat <- (sqrt(length(x)) * mean(x)) / sd(x) # 0.3007746

# What proportion of the t-statistics is larger than 2?
tstat_generator <- function(n, reps, seed = 1) {
  
  set.seed(seed)
  
  tstats <- replicate(reps, {
    x <- rnorm(n)
    tstat <- sqrt(length(x)) * mean(x) / sd(x)
  })
  
  return(tstats)
}

tstats <- tstat_generator(5, reps = 1000)
mean(tstats > 2) # 0.068

# Simulate 
sim_tstats <- map(seq(5, 30, 5), ~ tstat_generator(.x, 1000))

sample_sizes <- paste0("n_", seq(5, 30, 5))
names(sim_tstats) <- sample_sizes

# Use Monte Carlo simulation to confirm that the t-test between means computed from normal distributions follows a t-distribution
mc_ttest <- function(n, reps) {
  
  tstats <- replicate(reps, {
    x <- rnorm(n)
    y <- rnorm(n)
    tstat <- t.test(x, y, var.equal = TRUE)$statistic
  })
  
  return(tstats)
}

# All of these look pretty spot on, regardless of the sample size
# Let's look at the QQ-plots
b <- 1000
ps <- seq(1/b, 1-1/b, len = b)
df <- 2*5-2
mc_stat_list <- map(seq(5, 30, 5), ~ data.frame(tstat = mc_ttest(.x, 1000), 
                                                tquantile = qt(ps, df)))

names(mc_stat_list) <- paste0("n_", seq(5, 30, 5))

# Bind to a single data frame
mc_df <- bind_rows(mc_stat_list, .id = "sample_size") %>% 
  mutate(sample_size = as.numeric(gsub("n_", "", sample_size))) 
  
# Histograms
mc_df %>% 
  ggplot(aes(tstat)) +
  geom_histogram() +
  facet_wrap(~sample_size, scales = "free") + 
  labs(
    x = "t-statistic", 
    y = "Count", 
    title = "Distributions of t-statistics drawn from Monte Carlo simulations"
  )

# QQ plots-- need to fix these
mc_df %>% 
  filter(sample_size == 5) %>% 
  ggplot(aes()) +
  geom_point()
  stat_qq() + 
  stat_qq_line(
    distribution = stats::qt, 
    dparams = list(p = ps, df = df), 
    linetype = "dashed", color = "red"
  ) +
  labs(x = "Theoretical quantiles", y = "Sample quantiles") + 
  facet_wrap(~sample_size, scales = "free")

y <- filter(mc_df, sample_size == "5") %>% pull("tstat")
qqplot(qt(ps, df), y)
abline(0, 1)

# Sampling binary data
binary_tstat <- replicate(1000, {
  x <- sample(c(-1, 1), 15, replace = TRUE)
  tstat <- sqrt(15)*mean(x) / sd(x)}
)

data.frame(binary_tstat) %>% 
  ggplot(aes(binary_tstat)) + 
  geom_histogram()

# Now repeat with 1000 samples instead of 15
binary_tstat <- replicate(1000, {
  x <- sample(c(-1, 1), 1000, replace = TRUE)
  tstat <- sqrt(1000)*mean(x) / sd(x)}
)

# Compare this histogram with a standard normal distribution
data.frame(
  binary_tstat = binary_tstat,
  normal_dist = rnorm(1000)
) %>% 
  pivot_longer(everything(), names_to = "type", values_to = "value") %>% 
  ggplot(aes(value, fill = type)) + 
  geom_histogram(position = "stack")

# The distributions essentially stack on one another

# Use Monte Carlo sampling to compute the median from a standard normal distribution
medians <- replicate(10000, {median(rnorm(1000))})
data.frame(medians) %>% 
  ggplot(aes(medians)) + 
  geom_histogram()

# QQ-plot-- approximately normally distributed
data.frame(medians) %>% 
  ggplot(aes(sample = medians)) + 
  geom_qq(alpha = 0.7) + 
  geom_qq_line(linetype = "dashed", color = "red")

# Check the mean and standard deviation
mean(medians) # 0.0003675629
1/sqrt(10000) < sd(medians) # standard deviation is slightly larger

# The last option is true-- sample median is approx normal with mean nearly 0 and standard deviation slightly larger than 1/sqrt(N)

# Permutation tests -------------------------------------------------------

# Observed difference between mouse weights
control <- filter(female_weights, Diet == "chow") %>% pull("Bodyweight")
treatment <- filter(female_weights, Diet == "hf") %>% pull("Bodyweight")
obs_diff <- mean(treatment) - mean(control)

# Create a function to permute the data-- essentially, randomly sample bodyweights from control and treatment
# Then, assign the first half to the new, permuted control group
# The second half gets assigned to the new, permuted treatment group
# finally, return the difference between the mean of the new permuted variables
shuffle <- function(n) {
  
  all <- sample(c(control, treatment))
  permute_control <- all[1:n]
  permute_treatment <- all[(n + 1):(2 * n)]
  
  return(mean(permute_treatment) - mean(permute_control))
}

# Repeat the permutation 1 million times
set.seed(1)
mc_diffs <- replicate(1e6, {shuffle(12)})

# Plot histogram of the null hypothesis and the observed difference
data.frame(mc_diffs) %>% 
  ggplot(aes(mc_diffs)) + 
  geom_histogram() + 
  geom_vline(xintercept = obs_diff, linetype = "dashed", color = "red") + 
  labs(x = "Difference in means", 
       y = "Count", 
       title = "Null distribution of difference in means with 1 million Monte Carlo simulations")

# Calculate the proportion of permuted means that are greater than or equal to the observed
# Note, add 1 to avoid p-value misestimation-- see Phipson and Smyth 2010 for more discussion on the topic of estimating p-values using randomization tests
# Essentially, when we perform the Monte Carlo simulations, we are generating estimates of the p-values, not enumerating all possible combinations to report exact p-values (see Fisher's exact test)
# In Phipson and Smyth 2010, the authors state that on average, permuatation tests p-value estimates are off by about 1/m, where m = the number of permutations
# So, in this example, if we used 1000 permutations, our p-value estimate would be off by about 0.001

# To compute the exact Monte Carlo p-value, compute (b + 1)/(m + 1), where b = proportion of simulated test statistics >= the observed, and m = the number of Monte Carlo simulations
est_mc_pval <- function(obs_diff, mc_diffs) {
  
  x <- sum(abs(mc_diffs) > abs(obs_diff)) + 1
  y <- length(mc_diffs) + 1
  return(x / y)
  
}

est_mc_pval(obs_diff, mean_diff) # 0.05794206

# Now repeat with a smaller sample size
est_mc_pval(obs_diff, replicate(1e5, {shuffle(5)})) # 0.2226678

# Permutation exercises ---------------------------------------------------

# Bodyweights of babies born from smokers and nonsmokers
bwt_nonsmokers <- filter(babies, smoke == 0) %>% pull("bwt")
bwt_smokers <- filter(babies, smoke == 1) %>% pull("bwt")

# Generate random variabble based on sample size of 10
n <- 10
set.seed(1)
nonsmokers <- sample(bwt_nonsmokers, n)
smokers <- sample(bwt_smokers, n)
obs <- mean(smokers) - mean(nonsmokers) # observed difference from random sample

# Simulate null distribution with 1000 Monte Carlo samples
set.seed(1)
null_diffs <- replicate(
  n = 1000, 
  expr = {
    sample_size <- 10
    dat <- c(smokers, nonsmokers)
    shuffle <- sample(dat)
    perm_smokers <- shuffle[1:sample_size]
    perm_nonsmokers <- shuffle[(sample_size + 1):(2*sample_size)]
    perm_diff <- mean(perm_smokers) - mean(perm_nonsmokers)
  }
)

# Calculate permutation derived p-value-- add 1 to the numerator and denominator so the p-value won't be zero-- this can't be possible with complete enumeration, but may happen when random samples are drawn
(sum(abs(null_diffs) >= abs(obs)) + 1) / (length(null_diffs) + 1) # 0.1108891

# Repeat process for the median instead of the mean
obs <- median(smokers) - median(nonsmokers) # -7.5

# Simulate with 1000 MC samples
set.seed(1)
null_medians <- replicate(
  n = 1000,
  expr = {
    sample_size <- 10
    shuffle <- sample(c(smokers, nonsmokers)) 
    perm_smokers <- shuffle[1:sample_size]
    perm_nonsmokers <- shuffle[(sample_size + 1):(2*sample_size)]
    perm_diff <- median(perm_smokers) - median(perm_nonsmokers)
  }
)

(sum(abs(null_medians) >= abs(obs)) + 1) / (length(null_medians) + 1) # 0.2877123

# Association tests -------------------------------------------------------

# Hypothesis test-- uses a statistic, which is a summary of the data, and its probability distribution, with the idea being to compute the probability that a particular result could be obtained simply by chance
# The key behind statistical inference is making summary conclusions about a population from a sample

# Simulate the tea tasting experiment
tab <- matrix(c(3, 1, 1, 3), 2, 2)
rownames(tab) <- c("poured before", "poured after")
colnames(tab) <- c("guessed before", "guessed after")

# Perform the Fisher exact test
fisher.test(tab, alternative = "greater")$p.value # 0.24
(choose(4, 3) * choose(4, 1)) / choose(8, 4) # 0.2285714

# Chi-square test ---------------------------------------------------------

# When sample sizes are large enough, we may be able to use the chi-square approximation (generally if there are any cells < 5, we need to use the exact test)

# Simulate disease and control genotype data
df <- data.frame(
  disease = factor(
    c(rep(0, 180), rep(1, 20), rep(0, 40), rep(1, 10)), 
    labels = c("control", "cases")), 
  genotype = factor(
    c(rep("AA/Aa", 200), rep("aa", 50)), 
    levels = c("AA/Aa", "aa"))
)

# Odds-ratio-- in this case, the probability that you have the disease if you are genotype "aa" divided by the probability that you have the disease if you are "AA/Aa"
(tab <- table(df))
(disease_minor <- tab[2, 2] / tab[1, 2]) # 0.25, 10/40
(disease_major <- tab[2, 1] / tab[1, 1]) # 0.111111, 20/180
(odds <- disease_minor / disease_major) # 2.25

# With an odds-ratio of 2.25, that means you are more likely to have the disease if you have the genotype 'aa'
# Now calculate the expected numbers in each cell-- that is, if the groups with different genotypes AA/Aa (200 total) and aa (50 total) were randommly assigned "case" or "control" status
(p <- mean(df$disease == "cases")) # 0.12-- proportion of cases in general

# Build the expected table
x <- c(1-p, p) * sum(df$genotype == "AA/Aa") # calculate expected numbers of cases for this genotype
y <- c(1-p, p) * sum(df$genotype == "aa")
expected <- rbind(x, y)
dimnames(expected) <- dimnames(tab)
expected

# Now we want to ask-- what's the probability that we would see a table that deviates as much or more from the expected as ours does
chisq.test(tab) # p-value == 0.08857
chisq.test(tab*10) # p-value == 1.22e-09-- not a 1:1 relationship between sample size and odds ratio, meaning that as the sample size increases, the probability that a particular result will be in the tails of the distribution becomes lower and lower, but the odds ratio stays the same

# CIs for the odds ratio --------------------------------------------------

# The odds-ratio is actually a ratio of ratios, and there is no simple way to derive a distribution to compute confidence intervals-- but, we can use the log odds-ratio, which is asymptotically normal
# To compute the log odds-ratio, we can use a generalized linear model
fit <- glm(disease ~ genotype, family = "binomial", data = df)
summary(fit)$coefficients

# This shows the estimate for the log odds of disease given genotype aa-- 0.81
# To form the confidence intervals, we need to multiple -2 and 2 by the standard error of the estimate
(aa <- summary(fit)$coefficients[2, 1])
(se_aa <- summary(fit)$coefficients[2, 2])
(ll <- aa -2*se_aa)
(ul <- aa + 2*se_aa)

# To get the odds back on the original scale, we need to exponentiate
exp(c(ll, ul)) # includes 1, meaning a plausible value for the odds ratio is 1

# Association test exercises ----------------------------------------------

# Using the association test data
tab <- table(association_test)
chi_tab <- chisq.test(tab) # x-square statistic 3.3437

# Plot against chi-square distribution with same degrees of freedom
data.frame(chi =rchisq(1e6, 3)) %>% 
  ggplot(aes(chi)) + 
  geom_histogram() + 
  geom_vline(xintercept = chi_tab$statistic, linetype = "dashed")

# From this, we can see that this probably isn't statistically significant-- p-value 0.06746

# Compute the Fisher's exact p-value on the same table
fisher.test(tab)

