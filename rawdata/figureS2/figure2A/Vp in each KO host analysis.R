library(tidyverse)
library(DescTools)
library(cowplot)

setwd("C:/Users/cf466/OneDrive - University of Exeter/Documents/Documents - CF laptop/PhD/Phage Cocktails/272-950 work/EllieMutant Plaque Assay Data")
okabe <- c("#8dd3c7", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#a35eb5","#80c683","#8c8cd4", "#bc80bd")

reformat <- function(df) {
  colnames(df)[1] <- "Well"
  new_names <- sapply(colnames(df)[-1], function(x) substr(x, 2, nchar(x) - 1))
  colnames(df)[-1] <- new_names
  df_long <- pivot_longer(df, cols = -Well, names_to = "Time_Hours", values_to = "values")
  df_long$Time_Hours <- as.numeric(as.character(df_long$Time_Hours))
  df_long <- mutate(df_long, Time_Hours = round(Time_Hours / 3600 / 0.25) * 0.25)
  df_long <- left_join(df_long, dictionary)
  df_long <- df_long %>% na.omit()
  return(df_long)
}

dictionary <- read.csv("Dictionary.csv")
df <- read.csv("Vp in each KO host data 9.9.24.csv") %>%
  reformat()

df_blks <- df %>%
  filter(Bacteria == "None") %>%
  filter(Phage == "BLK") %>% 
  group_by(Time_Hours) %>%
  reframe(mean_blk_od = mean(values))

df <- df %>%
  left_join(df_blks) %>%
  mutate(norm_od = values - mean_blk_od)

# ggplot(df, aes(x = Time_Hours, y = norm_od, group = Well)) +
#   geom_line() +
#   facet_wrap(Bacteria~Phage, ncol = 4) +
#   theme_cowplot(14)

df_control_auc <- df %>%
  filter(Phage == "Control") %>%
  group_by(Well, Bacteria) %>%
  reframe(cntrl_auc = AUC(x=Time_Hours, y=values, method='trapezoid')) %>%
  group_by(Bacteria) %>%
  reframe(mean_cntrl_auc = mean(cntrl_auc))

df_viurlence <- df %>%
  group_by(Well, Bacteria, Phage) %>%
  reframe(test_auc = AUC(x=Time_Hours, y=values, method='trapezoid')) %>%
  left_join(df_control_auc) %>%
  na.omit() %>%
  mutate(Vp = 1 - (test_auc/mean_cntrl_auc)) %>%
  mutate(Bacteria = factor(Bacteria, levels = c("PAO1.HsdR", "PAO1.PilA", "PAO1.GalU", "PAO1.PilA.GalU"))) %>%
  mutate(Phage = factor(Phage, levels = c("CPL00272", "CPL00950", "Cocktail", "Control")))

df_mean_vp <- df_viurlence %>%
  group_by(Bacteria, Phage) %>%
  reframe(mean_vp = mean(Vp))

Vp_plot <- ggplot(df_viurlence, aes(x=Bacteria, y=Vp, fill=Phage, colour = Phage)) + 
  geom_col(data = df_mean_vp, aes(x=Bacteria, y=mean_vp), position = "dodge", alpha = 0.5)+
  geom_point(position = position_dodge(width = 0.8)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), name = "Virulence Index (Vp)", breaks = seq(-1, 1.5, 0.2), limits = c(-1,1)) +
  scale_x_discrete(labels = c("Parent Strain", "??pilA", "??galU", "??pilA ??galU")) +
  #geom_point(data = df_stats, aes(x=bacterial_strain, y=EOP),
  #           position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), size = 2) +
  ggtitle("Vp of Cocktail Phage on each receptor knockout") +
  scale_fill_manual(values = c("CPL00272"="#E69F00", "CPL00950" = "#56B4E9","Cocktail" = "#009E73", "Control" = "black"))+
  scale_colour_manual(values = c("CPL00272"="#E69F00", "CPL00950" = "#56B4E9","Cocktail" = "#009E73", "Control" = "black"))+
  theme_cowplot() +
  xlab("PAO1 Receptor Deletion Mutant") +
  theme(panel.border = element_rect(colour = "black", fill=NA),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16, angle=80, hjust=1),
        plot.title = element_text(size=12, hjust=0.5))
Vp_plot

df_stats <- df_viurlence %>%
  mutate(Phage = factor(Phage, levels = c("Control", "CPL00272", "CPL00950", "Cocktail")))

# model <- lm(Vp ~ Bacteria * Phage, data = df_stats)
# summary(model)
# plot(model)

df_stats_2 <- df_stats %>%
  filter(Phage == "CPL00272") %>%
  mutate(Bacteria = factor(Bacteria, levels = c("PAO1.PilA.GalU", "PAO1.HsdR", "PAO1.PilA", "PAO1.GalU")))

model_2 <- lm(Vp ~ Bacteria, data = df_stats_2)
summary(model_2)
plot(model)

