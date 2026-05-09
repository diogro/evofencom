ratones <- read.table("https://raw.githubusercontent.com/diogro/evofencom/refs/heads/main/Tutoriais/ratones.tsv", header = TRUE)

pak::pkg_install("diogro/ratones")
rgl.useNULL = TRUE
library(ratones)

traits = names(ratones)[12:46]
f = paste0("cbind(", paste(traits, collapse = ","), ") ~ SEX + line") 
cov_mat = lm(as.formula(f), data = ratones) |> CalculateMatrix()

cor_mat = cov2cor(cov_mat)

cor_df = as.data.frame(cor_mat)
library(tidyr)
library(dplyr)

cor_raw = cor(ratones[, traits])

as.data.frame(cor_mat - cor_raw) |>
  mutate(t = colnames(cor_mat)) |>
pivot_longer(, col = traits) |>
  filter(t == "NSL_NA") |>
  arrange(desc(value)) |> print(n = 35)
  dat = ratones |>
  select(SEX, line, NA_BR, LD_AS)

cor_df["NSL_NA", "IS_PNS"]
cor_raw["NSL_NA", "IS_PNS"]

library(ggplot2)

ggplot(ratones, aes(NSL_NA, IS_PNS, color = SEX, shape = line, group = line)) + 
  geom_point() + stat_ellipse() + scale_color_manual(values = 1:2) + theme_classic()

ratones_residuals = ratones
ratones_residuals[,traits] = residuals(lm(as.formula(f), data = ratones))

ggplot(ratones_residuals, aes(NSL_NA, IS_PNS, color = SEX, shape = line)) + 
  geom_point() + scale_color_manual(values = 1:2) + theme_classic()

