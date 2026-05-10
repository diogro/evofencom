if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(wesanderson)) { install.packages("wesanderson"); library(wesanderson) }
if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(ggthemes)){install.packages("ggthemes"); library(ggthemes)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}

plot_pop = function(p, mean, sd, multi, ...){
  x <- seq(-1,11,0.01)
  MyDF <- data.frame(x = x, y = dnorm(x, mean, sd)*multi)
  shade <- rbind(c(0.1, 0), MyDF, c(MyDF[nrow(MyDF), "X"], 0))
  pop = p + geom_line(data = MyDF, aes(x, y)) + geom_polygon(data = shade, aes(x, y), ...)
  return(pop)
}

p1 <- ggplot() + labs(x = "Um locus", y= "") + theme_tufte() + theme(panel.grid.major = element_blank(),
                                                                               panel.grid.minor = element_blank(),
                                                                               axis.text = element_blank())+
  scale_y_continuous(breaks = NULL) + scale_x_continuous(breaks = NULL)
p1 = plot_pop(p1, 5.0, 0.8, 0.9, color = "lightblue", alpha = 0)
p1 = plot_pop(p1, 6.5, 1, 0.5, color = "lightblue", alpha = 0)
p1 = plot_pop(p1, 3.5, 1, 0.5, color = "lightblue", alpha = 0)
p1 = plot_pop(p1, 5, 1.57, 1.8, color = "blue", alpha = 0) + geom_hline(yintercept = 0)

p2 <- ggplot() + labs(x = "Dois loci", y= "") + theme_tufte() + theme(panel.grid.major = element_blank(),
                                                                     panel.grid.minor = element_blank(),
                                                                     axis.text = element_blank()) +
  scale_y_continuous(breaks = NULL) + scale_x_continuous(breaks = NULL)
p2 = plot_pop(p2, 5.0, 0.8, 0.7, color = "lightblue", alpha = 0)
p2 = plot_pop(p2, 6, 1, 0.5, color = "lightblue", alpha = 0)
p2 = plot_pop(p2, 4, 1, 0.5, color = "lightblue", alpha = 0)
p2 = plot_pop(p2, 7.5, 1, 0.1, color = "lightblue", alpha = 0)
p2 = plot_pop(p2, 2.5, 1, 0.1, color = "lightblue", alpha = 0)
p2 = plot_pop(p2, 5, 1.57, 1.8, color = "blue", alpha = 0) + geom_hline(yintercept = 0)
p2
p1p2 = plot_grid(p1, p2, ncol = 2, labels = c("A", "B"))
save_plot("~/Dropbox/Cursos/ModCurso2019/Aulas/IntroQuantGen/figures/discrete_gaussian.png", p1p2, base_height = 5, ncol = 2, base_aspect_ratio = 1.2)

## 

# Load required library
library(ggplot2)
library(patchwork)
library(cowplot)

# --- Parameters ---
p <- 0.1      # Population frequency of the focal allele (0 to 1)
a <- 1.0      # Additive effect (mean phenotype for genotype 1)
d <- 0.0      # Dominance effect (mean phenotype for genotype 1/2)
sigma <- 0.25 # Variation (standard deviation) around the mean phenotype
N <- 1000     # Total number of individuals to simulate

plot_average_effect <- function(p, a, d, sigma = 0.25, N = 1000){
# --- Simulation ---
# 1. Calculate expected genotype frequencies (Hardy-Weinberg proportions)
freq_0    <- (1 - p)^2
freq_half <- 2 * p * (1 - p)
freq_1    <- p^2

# 2. Calculate the number of individuals per genotype
n_0    <- round(N * freq_0)
n_half <- round(N * freq_half)
n_1    <- N - n_0 - n_half # Ensure the sum is exactly N

# 3. Generate genotype data (Within-individual allele frequency)
genotypes <- c(rep(0, n_0), rep(0.5, n_half), rep(1, n_1))

# 4. Generate corresponding mean genotypic values
mean_values <- c(rep(-a, n_0), rep(d, n_half), rep(a, n_1))

# 5. Simulate phenotypes by adding normally distributed noise (sigma)
phenotypes <- rnorm(N, mean = mean_values, sd = sigma)

# 6. Create a dataframe for plotting
df <- data.frame(
  Genotype = genotypes,
  Phenotype = phenotypes
)
up = max(a, d)
lw = min(-a, d)
# --- Plotting ---
ggplot(df, aes(x = Genotype, y = Phenotype)) +
  # Add points with horizontal jitter to prevent overplotting
  geom_jitter(width = 0.02, alpha = 0.4, color = "steelblue") +
  
  # Add the true genotypic means as large points
  stat_summary(fun = mean, geom = "point", color = "darkred", size = 4, shape = 18) +
  
  # Add the linear regression line (Average Effect of Substitution)
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 1) +
  
  # Format the X-axis exactly as requested
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "1/2", "1")) +
  scale_y_continuous(, limits = c(lw - 3*sigma, up + 3*sigma)) + 
  
  # Labels and theming
  labs(
    title = "Average Effect of Allele Substitution",
    subtitle = paste0("Parameters: p = ", p, ", a = ", a, ", d = ", d, ", sigma = ", sigma),
    x = "Within-Individual Allele Frequency",
    y = "Phenotype Value"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )
}
p1 = plot_average_effect(0.1, 1, 0.6)
p2 = plot_average_effect(0.5, 1, 0.6)
p3 = plot_average_effect(0.9, 1, 0.6)


save_plot("average_excess1.png", p1, base_height = 5, base_asp = 1.3)
save_plot("average_excess2.png", p2, base_height = 5, base_asp = 1.3)
save_plot("average_excess3.png", p3, base_height = 5, base_asp = 1.3)

p1 = plot_average_effect(0.9, 0.5, 1)
p2 = plot_average_effect(0.6, 0.5, 1)
p3 = plot_average_effect(0.1, 0.5, 1)

save_plot("average_excess_dom1.png", p1, base_height = 5, base_asp = 1.3)
save_plot("average_excess_dom2.png", p2, base_height = 5, base_asp = 1.3)
save_plot("average_excess_dom3.png", p3, base_height = 5, base_asp = 1.3)


library(ggplot2)

# --- Parameters ---
p <- 0.25     # Population frequency of the focal allele (0 to 1)
a <- 1.25     # Additive effect (mean for genotype 1)
d <- 0.125      # Dominance effect (mean for genotype 1/2)
sigma <- 0.01 # Phenotypic variation (sd) around the genotypic mean
N <- 5000     # Increased N for a smoother histogram

plot_allele_contribution <- function(p, a, d, sigma = 0.25, N = 5000){
  
# --- Simulation ---
freq_0    <- (1 - p)^2
freq_half <- 2 * p * (1 - p)
freq_1    <- p^2

n_0    <- round(N * freq_0)
n_half <- round(N * freq_half)
n_1    <- N - n_0 - n_half 

# Generate simulated data
genotypes <- c(rep("0", n_0), rep("1/2", n_half), rep("1", n_1))
mean_values <- c(rep(-a, n_0), rep(d, n_half), rep(a, n_1))
phenotypes <- rnorm(N, mean = mean_values, sd = sigma)

df <- data.frame(Genotype = factor(genotypes, levels = c("0", "1/2", "1")), 
                 Phenotype = phenotypes)

# --- Theoretical Density Functions for Overlays ---
# Scaled by their HW population frequencies
dist_0    <- function(x) dnorm(x, mean = -a, sd = sigma) * freq_0
dist_half <- function(x) dnorm(x, mean = d,  sd = sigma) * freq_half
dist_1    <- function(x) dnorm(x, mean = a,  sd = sigma) * freq_1

# Total population density is the sum of the components
dist_total <- function(x) dist_0(x) + dist_half(x) + dist_1(x)

up = max(a, d)
lw = min(-a, d)

# --- Plotting ---
ggplot(df, aes(x = Phenotype)) +
  # 1. Histogram of the simulated population data
  geom_histogram(aes(y = after_stat(density)), bins = 40, 
                 fill = "grey85", color = "grey60", alpha = 0.7) +
  
  # 2. Overlay individual genotype contributions (dashed lines)
  stat_function(fun = dist_0, color = "#E41A1C", linewidth = 1.2, alpha = 0.5) +
  stat_function(fun = dist_half, color = "#4DAF4A", linewidth = 1.2, alpha = 0.5) +
  stat_function(fun = dist_1, color = "#377EB8", linewidth = 1.2, alpha = 0.5) +
  
  scale_x_continuous(, limits = c(lw - 4*sigma, up + 4*sigma)) + 
  
  # 3. Overlay total theoretical population distribution (thick solid line)
  #stat_function(fun = dist_total, color = "black", linewidth = 1.5) +
  
  # 4. Labels and theming
  labs(
    title = "Locus Contribution to Phenotypic Variance",
    subtitle = bquote("Parameters:" ~ p == .(p) * "," ~ a == .(a) * "," ~ d == .(d) * "," ~ sigma == .(sigma)),
    x = "Phenotypic Value",
    y = "Probability Density"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )
}
plot_allele_contribution(0.5, a=1, d=0, sigma = 1, 1000)


library(ggplot2)

# --- Parameters ---
p_A <- 0.5   # Allele frequency at locus A
a_A <- 1.0   # Additive effect of locus A
d_A <- 0.   # Dominance effect of locus A (shift of the heterozygote)

p_B <- 0.5   # Allele frequency at locus B
a_B <- 1.0   # Additive effect of locus B
d_B <- 0.  # Dominance effect of locus B

sigma <- .5 # Phenotypic variation (sd) around the genotypic mean
N <- 10000   # Population size

# --- Simulation: Calculate the 9 genotypic classes ---
# Locus A probabilities and effects (now including dominance)
freqs_A <- c((1 - p_A)^2, 2 * p_A * (1 - p_A), p_A^2)
effs_A  <- c(-a_A, d_A, a_A) 

# Locus B probabilities and effects (now including dominance)
freqs_B <- c((1 - p_B)^2, 2 * p_B * (1 - p_B), p_B^2)
effs_B  <- c(-a_B, d_B, a_B)

# Create the 9 combinations assuming linkage equilibrium (no epistasis)
classes <- expand.grid(eff_A = effs_A, eff_B = effs_B)
classes$mean_val <- classes$eff_A + classes$eff_B

freq_grid <- expand.grid(freq_A = freqs_A, freq_B = freqs_B)
classes$freq <- freq_grid$freq_A * freq_grid$freq_B

# --- Simulate the Population ---
# Sample from the 9 classes based on their combined frequencies
sampled_classes <- sample(1:9, size = N, replace = TRUE, prob = classes$freq)
phenotypes <- rnorm(N, mean = classes$mean_val[sampled_classes], sd = sigma)

df <- data.frame(Phenotype = phenotypes)

# --- Theoretical Density Function ---
dist_total <- function(x) {
  y <- numeric(length(x))
  for(i in 1:nrow(classes)) {
    y <- y + dnorm(x, mean = classes$mean_val[i], sd = sigma) * classes$freq[i]
  }
  return(y)
}

# --- Plotting ---
p_plot <- ggplot(df, aes(x = Phenotype)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, 
                 fill = "grey85", color = "grey60", alpha = 0.7)

# Add the 9 individual sub-population densities
for(i in 1:nrow(classes)) {
  p_plot <- p_plot + stat_function(
    fun = function(x, mean, freq) dnorm(x, mean = mean, sd = sigma) * freq,
    args = list(mean = classes$mean_val[i], freq = classes$freq[i]),
    color = "steelblue", linewidth = 0.8, alpha = 0.6
  )
}

# Add the total theoretical population distribution
p_plot <- p_plot + 
  labs(
    title = "Two Loci: Contribution to Variance",
    subtitle = paste0("A: p=", p_A, ", a=", a_A, ", d=", d_A, " | B: p=", p_B, ", a=", a_B, ", d=", d_B),
    x = "Phenotypic Value",
    y = "Probability Density"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )

print(p_plot)

