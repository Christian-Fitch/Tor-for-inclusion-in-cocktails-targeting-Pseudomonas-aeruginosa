library(tidyverse)
library(DescTools)
library(cowplot)
library(car)

setwd("C:/Users/cf466/OneDrive - University of Exeter/Documents/Documents - CF laptop/PhD/Phage Cocktails/272-950 work/EllieMutant Plaque Assay Data")
df = read.csv("christian cocktail phage pfus.csv")

df <- df %>%
  mutate(phage = factor(phage, levels = c("P272","P950","cocktail")))

df_average <- df %>% 
  mutate(bacterial_strain = factor(bacterial_strain, levels = c("PAO1.hsdr", "PAO1.oprM","PAO1.pilA", "PAO1.wbpL", "PAO1.galU", "PAO1.algC", "PAO1.pilA.galU", "PAO1.pilA.algC"))) %>%
  mutate(phage = factor(phage, levels = c("P272","P950","cocktail"))) %>%
  group_by(bacterial_strain, phage) %>%
  reframe(av_PFU = mean(pfu), stdev = sd(pfu))


#log10 version

ggplot(df_average, aes(x=bacterial_strain, y=log10(av_PFU), fill=phage)) + 
  geom_bar(stat = "identity", position="dodge", alpha = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), name = "Log10 PFU/mL", breaks = seq(0, 13, 1)) +
  scale_x_discrete(labels = c("Parent Strain", "ΔoprM", "ΔpilA", "ΔwbpL", "ΔgalU", "ΔalgC", "ΔpilA ΔgalU", "ΔpilA ΔalgC")) +
  #geom_errorbar(aes(ymin = av_PFU, ymax = av_PFU + stdev), width = 0.2) +
  geom_point(data = df, aes(x=bacterial_strain, y=log10(pfu)),
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), size = 2) +
  scale_fill_manual(values = c("P272"="#E69F00", "P950" = "#56B4E9","cocktail" = "#009E73"))+
  theme_cowplot() +
  xlab("PAO1 Receptor Deletion Mutant") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        plot.title = element_text(size=12, hjust=0.5))




#  not log10 version

ggplot(df_average, aes(x=bacterial_strain, y=av_PFU, fill=phage)) + 
  geom_bar(stat = "identity", position="dodge") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(labels = c("∆hsdR", "ΔoprM", "ΔpilA", "ΔwbpL", "ΔgalU", "ΔalgC", "ΔpilA ΔgalU", "ΔpilA ΔalgC")) +
  #geom_errorbar(aes(ymin = av_PFU, ymax = av_PFU + stdev), width = 0.2) +
  geom_point(data = df, aes(x=bacterial_strain, y=pfu),
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), size = 2) +
  scale_fill_manual(values = c("P272"="#E69F00", "P950" = "#56B4E9","cocktail" = "#009E73"))+
  theme_bw() +
  ylab("Average PFU/ml") +
  xlab("Bacterial Strain") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        plot.title = element_text(size=12, hjust=0.5))

# stats
setwd("C:/Users/cf466/OneDrive - University of Exeter/Documents/Documents - CF laptop/PhD/Phage Cocktails/272-950 work/EllieMutant Plaque Assay Data")
df = read.csv("christian cocktail phage pfus.csv")

df <- df %>%
  mutate(phage = factor(phage, levels = c("P272","P950","cocktail")))

df_stats <- df %>% 
  mutate(bacterial_strain = factor(bacterial_strain, levels = c("PAO1.hsdr", "PAO1.oprM","PAO1.pilA", "PAO1.wbpL", "PAO1.galU", "PAO1.algC", "PAO1.pilA.galU", "PAO1.pilA.algC"))) %>%
  mutate(phage = factor(phage, levels = c("P272","P950","cocktail")))

df_control <- df_stats %>%
  filter(bacterial_strain == "PAO1.hsdr") %>%
  group_by(phage) %>%
  reframe(mean_control = mean(pfu))

df_stats <- df_stats %>% 
  left_join(df_control) %>% 
  mutate(EOP = pfu/mean_control)

#want differences in EOP between bacterial strain (compared to Control) to be tested for each phage.
# this will look something like EOP ~ phage*strain

#model <- lm(EOP ~ phage*bacterial_strain, data = df_stats)

#no sig difference between 272 and 950 or cocktail, so no difference, as expected.
#no sig difference in oprM, pilA, wbpL mutants, indicating these receptors are not
#involved in phage infection in any way.
#sig difference in galU, algC, pilA/algC, pilA/galU.
#no significant interaction terms between the strains and their phages, so the model
# should look like: EOP ~ strain

#model_updated <- lm(log10(pfu+0.0000001) ~ bacterial_strain, data = df_stats)
# summary(model_updated)
# plot(model_updated)

#not great on the plots, need to refine this. As technically the double mutants are 
#showing the same things, going to log the EOP values to limit the difference
#between the ones that infect and the ones that don't?

#also not sure if I can have the double mutants... they are technically not
#independent no?

df_stats_2 <- df_stats %>%
  filter(bacterial_strain != "PAO1.pilA.galU") %>%
  filter(bacterial_strain != "PAO1.pilA.algC")

model_updated <- lm(EOP ~ bacterial_strain, data = df_stats_2)
summary(model_updated)
plot(model_updated)

#not normal, going to bootstrap the graph!

betahat.boot <- Boot(model_updated, R=1000) # 1000 bootstrap samples
summary(betahat.boot)  # default summary
confint(betahat.boot)
hist(betahat.boot)

df_stats_3 <- df_stats_2 %>%
  group_by(bacterial_strain) %>%
  reframe(mean_EOP = mean(EOP))





















