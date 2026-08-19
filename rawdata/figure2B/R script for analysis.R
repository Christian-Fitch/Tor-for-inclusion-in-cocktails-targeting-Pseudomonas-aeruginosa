library(tidyverse)          #load in the package
library(DescTools)          #more utility tools, particularly area under curves
library(cowplot)            #for plot_grid function
library(scales)             #for the superscirpt in ggplots

#Boring beginning things to get my data in 
setwd("C:/Users/cf466/OneDrive - University of Exeter/Documents/Documents - CF laptop/PhD/Phage Cocktails/272-950 work/Phage breadth 27.08.24")
okabe <- c("#8dd3c7", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#a35eb5","#80c683","#8c8cd4", "#bc80bd")
dictionary <- read.csv("Dictionary.csv")
strains <- read.csv("strains.csv")

#read in to test 1 dataframe before writing a function for it
setwd("C:/Users/cf466/OneDrive - University of Exeter/Documents/Documents - CF laptop/PhD/Phage Cocktails/272-950 work/Phage breadth 27.08.24/0 hour")

#test <- read.csv("Wizard script `PhageCocktails`, Prot.1, Plate 1, 22 August 2024 150758.CSV")
#want to remove the columns "X" and "Content" and rename " Raw.Data..600."
#Also want to remove the first row of the df
#Will also need to create and left_join a dictionary to this!

# test <- read.csv("Wizard script `PhageCocktails`, Prot.1, Plate 1, 22 August 2024 150758.CSV", header = TRUE) %>%
#   select(-c("X", "Content")) %>%
#   rename(OD600 =  Raw.Data..600.) %>%
#   slice(-1) %>%
#   left_join(dictionary) %>%
#   left_join(strains)

#okay, so now creating the reformat function to do these things

reformat <- function(df){
  df <- df %>%
    select(-Content) %>%
    slice(-1) %>%
    left_join(dictionary)
  colnames(df)[2] <- "OD600"
  return(df)
}

#now need to make the function that reads in the files in the current working directory
#and adds the plate number, and the time.

produce_df <- function(file_name, time){
  df <- read.csv(file_name, header = TRUE, row.names = NULL)
  df <- reformat(df)
  file_base_name <- sub(".*(Plate \\d+).*", "\\1", basename(file_name))
  df$plate <- file_base_name
  df$time_hours <- time
  df <- df %>% left_join(strains)
  return(df)
}

#Now to go through this for timepoint 0 to make sure I can do it!

setwd("C:/Users/cf466/OneDrive - University of Exeter/Documents/Documents - CF laptop/PhD/Phage Cocktails/272-950 work/Phage breadth 27.08.24/0 hour")
csv_files <- list.files(pattern = "\\.CSV$")

#use purrr::map_df to apply the reformat_csv function to each CSV file and bind them together
df_0 <- map2(csv_files, c(0), produce_df)

#repeat for the other timepoints and then cbind them together
#t19

setwd("C:/Users/cf466/OneDrive - University of Exeter/Documents/Documents - CF laptop/PhD/Phage Cocktails/272-950 work/Phage breadth 27.08.24/19 hour")
csv_files <- list.files(pattern = "\\.CSV$")

#use purrr::map_df to apply the reformat_csv function to each CSV file and bind them together
df_19 <- map2(csv_files, c(19), produce_df)

#t24

setwd("C:/Users/cf466/OneDrive - University of Exeter/Documents/Documents - CF laptop/PhD/Phage Cocktails/272-950 work/Phage breadth 27.08.24/24 hour")
csv_files <- list.files(pattern = "\\.CSV$")

#use purrr::map_df to apply the reformat_csv function to each CSV file and bind them together
df_24 <- map2(csv_files, c(24), produce_df)

#t43

setwd("C:/Users/cf466/OneDrive - University of Exeter/Documents/Documents - CF laptop/PhD/Phage Cocktails/272-950 work/Phage breadth 27.08.24/43 hour")
csv_files <- list.files(pattern = "\\.CSV$")

#use purrr::map_df to apply the reformat_csv function to each CSV file and bind them together
df_43 <- map2(csv_files, c(43), produce_df)

#t48

setwd("C:/Users/cf466/OneDrive - University of Exeter/Documents/Documents - CF laptop/PhD/Phage Cocktails/272-950 work/Phage breadth 27.08.24/48 hour")
csv_files <- list.files(pattern = "\\.CSV$")

#use purrr::map_df to apply the reformat_csv function to each CSV file and bind them together
df_48 <- map2(csv_files, c(48), produce_df)

#produce a df_all for all the data
df_produce <- function(df_list){
  final_df <- Reduce(full_join, df_list)
  return(final_df)
}


#bind them all together
df_0 <- df_produce(df_0)
df_19 <- df_produce(df_19)
df_24 <- df_produce(df_24)
df_43 <- df_produce(df_43)
df_48 <- df_produce(df_48)
df_all <- rbind(df_0,df_19,df_24,df_43,df_48)

#remove the annoying X column
df_all <- df_all %>% select(-X)

#now as I messed up and added phage to B07 to B11, remove these! It also 
#turns out the OD600 row is a character, changing that to numeric!
dodgy <- c("B06", "B07", "B08", "B09", "B10", "B11")
df_all <- df_all %>%
  filter(!(Well %in% dodgy)) %>%
  mutate(OD600 = as.numeric(OD600))

#for each time point, need to normalise to the blanks. 
#There was a lot of variation it seems across the plates, so want to explore how bad
#That is first.

# df_blks <- df_all %>%
#   filter(bacteria == "LB")

# ggplot(df_blks, aes(x = time_hours, y = OD600, group = time_hours)) +
#   geom_point() +
#   geom_boxplot()

#It seems that generally, the values here are very similar for all wells, except for 
#timepoint 48. This has some wells where the values are way too high, so need to take an average
#for each plate at each timepoint and use this.

# df_blks <- df_blks %>%
#   group_by(time_hours) %>%
#   reframe(meanblkod = mean(OD600), stdev = sd(OD600))

#as expected the SD for those final wells high, might do this again, but filter out 
#any values over 0.75 as these seem very wrong

df_blks <- df_all %>%
  filter(bacteria == "LB") %>%
  filter(OD600 < 0.75) %>%
  group_by(time_hours) %>%
  reframe(meanblkod = mean(OD600))

# this looks a lot better. Now to leftjoin this back into the main dataframe
#while im there, can also normalise ODs!
#got quite a few negative values. Normalise these to 0 and remove the LB blanks

df_all <- df_all %>%
  left_join(df_blks) %>%
  mutate(NormOD = OD600 - meanblkod) %>%
  mutate(NormOD = ifelse(NormOD < 0, 0, NormOD)) %>%
  filter(bacteria != "LB")

#now for the normalisation purposes, want to calculate reduction percentage for each bacteria
#at each timepoint

df_controls <- df_all %>%
  filter(Phage == "Control") %>% 
  select(Well, NormOD, Phage, bacteria, time_hours) %>%
  group_by(bacteria, time_hours) %>%
  reframe(meancontrolod = mean(NormOD))

#now to do the maths. Also want to remove the timepoint 0s as theyre not really showing anything,
#and select only relevant columns as it is getting a tad messy. If the bacteria hasnt grown yet, then
#there will be a non-number (from 0/0 or if there is x/0 then -Inf). These only occurred at timepoint
# 19 in bacteria that dont grow very quickly, I am going to remove them.

df_all <- df_all %>%
  left_join(df_controls) %>%
  mutate(red_perc = ((meancontrolod-NormOD)/meancontrolod)*100) %>%
  filter(time_hours != 0) %>%
  na.omit() %>%
  filter(red_perc != "-Inf")

#now this is sorted, can create some preliminary plots.
#I need something to show a range of concentrations, over diferent times, across a number of
#phage and bacteria

#as this will be a lot of data, want a variation of the dataframe that has means calculated
#timepoint 43 looks too varied, going to remove it

df_all_means <- df_all %>%
  group_by(time_hours, bacteria, Phage, Concentration) %>%
  reframe(mean_red = mean(red_perc), sd = sd(red_perc)) %>%
  filter(Phage != "Control") %>%
  filter(time_hours != 43)

# ggplot(df_all_means, aes(x = Concentration, y = mean_red, colour = Phage)) +
#   geom_col() +
#   scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
#   geom_hline(yintercept = -100) +
#   geom_hline(yintercept = 100) +
#   geom_errorbar(aes(ymin = mean_red-sd, ymax = mean_red+sd), width = 0.2) + 
#   facet_grid(bacteria ~ Phage + time_hours) +
#   theme_cowplot(14)

#Not amazing in the sense that some are highly varied, and others are very much outside the
#bounds of okay. But, lets assess this in the heatmap like I would like to!
#This completely skews the plot, will need to do some manipulations.

#If reduction percentage is less than 0, it is being set to 0.

df_all_means_plot <- df_all_means %>%
  mutate(mean_red = ifelse(mean_red < 0, 0, mean_red)) %>%
  mutate(time_hours = ifelse(time_hours == 19, "19 Hours",
                             ifelse(time_hours == 24, "24 Hours",
                                    ifelse(time_hours == 48, "48 Hours", NA)))) %>%
  mutate(Phage = factor(Phage, levels = c("CPL00272", "CPL00950", "Cocktail"))) %>%
  mutate(Concentration = factor(Concentration))

custom_labels <- c("10\u00B2", "10\u2074", "10\u2076", "10\u2077", "10\u2078")
levels(df_all_means_plot$Concentration) <- custom_labels

ggplot(df_all_means_plot, aes(x = Phage, y = Concentration, fill = mean_red)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#968be0", name = "     Mean\n Reduction\nPercentage") +
  scale_x_discrete("Phage") +
  scale_y_discrete("Phage Concentration (PFU/mL)", labels = levels(df_all_means_plot$Concentration)) +
  facet_grid(bacteria ~ time_hours) +
  theme_cowplot(14) +
  theme(strip.text.y = element_text(angle = 0),
        strip.background = element_blank())


#so my question is, can we calculate the virulence of each phage on each host? 

df_test <- df_all_means_plot %>%
  mutate(time = as.numeric(str_extract(time_hours, "\\d+"))) %>%
  filter(Phage != "Cocktail")

ggplot(df_test, aes(x = time, y = mean_red)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean_red - sd, ymax = mean_red + sd), width = 0.2) +
  scale_y_continuous(limits = c(-150,150)) +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_hline(yintercept = 100, colour = "grey") +
  facet_grid(bacteria~Phage+Concentration) +
  theme_bw(14)

Pbreadth_analysis <- df_all %>%
  mutate(time = as.numeric(str_extract(time_hours, "\\d+"))) %>%
  filter(Phage != "Cocktail") %>%
  filter(Phage != "Control") %>%
  filter(time_hours != 43) %>%
  select(bacteria, Phage, time, Concentration, red_perc) %>%
  mutate(reduction = ifelse(red_perc > 50, 1, 0)) %>%
  select(bacteria, Phage, time, Concentration, reduction) %>%
  group_by(bacteria, Phage, time, Concentration) %>%
  reframe(percentage = 100*(sum(reduction)/3))

Pbreadth_analysis_beyondplot <- Pbreadth_analysis %>%
  select(bacteria, Phage, Concentration, time, percentage) %>%
  mutate(time = ifelse(time == 19, "h19", 
                       ifelse(time == 24, "h24",
                              ifelse(time == 48, "h48", NA)))) %>%
  pivot_wider(names_from = time, values_from = percentage) %>%
  mutate(at_least_one_timepoint = ifelse(h19 > 70|h24 > 70|h48 > 70, 1, 0)) %>%
  filter(at_least_one_timepoint == 1) %>%
  group_by(bacteria, Phage, Concentration) %>%
  reframe(total = sum(at_least_one_timepoint)) %>%
  group_by(Phage, Concentration) %>%
  reframe(count = sum(total))



  
ggplot(Pbreadth_analysis, aes(x = time, y = percentage, colour = as.factor(Concentration))) +
  geom_point() +
  geom_line() +
  facet_grid(bacteria~Phage+Concentration) +
  theme_bw(14)
















#now to run the stats. The question is here: is the cocktail better than 272 on its own? 
#to be more specific... is the reduction percentage at each of three timepoints higher 
#for the cocktail of CPL00950 vs CPL00272. I want to account for concentration as a
#continuous variable, but as this is meant to be representative of the patient population,
#I do not want to account for host strain. Perhaps that can be built in, but I'd possibly not
#Include each host as a separate thing.

#therefore, the model will be something like ReductionPercentage ~ Phage + Conc where
#CPL00272 is the baseline, CPL00950 is removed, and the cocktail is what its being compared to.
#Will also run a first model with the control as the baseline to make sure the phages are doing something

df_stats <- df_all %>%
  filter(Phage != "CPL00950") %>%
  mutate(Phage = factor(Phage, levels = c("Control", "CPL00272", "Cocktail")))

#null = CPL00272 and Cocktail provide no evidence of reduction
model_blk <- lm(red_perc ~ Phage + Concentration, df_stats)
plot(model_blk)
summary(model_blk)

#geninely horrendous model, think there are a lot of extremes, need to set min and max to 0 and 100 I think
df_stats <- df_all %>%
  filter(Phage != "CPL00950") %>%
  mutate(Phage = factor(Phage, levels = c("Control", "CPL00272", "Cocktail"))) %>%
  mutate(red_perc = ifelse(red_perc > 100, 100, 
                           ifelse(red_perc < 0, 0, red_perc)))

#null = CPL00272 and Cocktail provide no evidence of reduction
model_blk <- lm(red_perc ~ Phage + Concentration, df_stats)
plot(model_blk)
summary(model_blk)

#is a lot better, not hugely normal, but a lot better. 
#maybe 2 transforms to see if there is a better looking model
#not really the one above was the best I think. But It shows that phage were doing something.

#now to directly compare cocktail and cpl00272
df_stats <- df_all %>%
  filter(Phage != "CPL00950") %>%
  filter(Phage != "Control") %>%
  mutate(Phage = factor(Phage, levels = c("CPL00272", "Cocktail"))) %>%
  mutate(red_perc = ifelse(red_perc > 100, 100, 
                           ifelse(red_perc < 0, 0, red_perc)))

#null = no difference between CPL00272 and Cocktail at each of 19, 24, 48
df_stats_19 <- df_stats %>%
  filter(time_hours == 19)

model_cocktail_19 <- lm(red_perc ~ Phage + Concentration, df_stats_19)
plot(model_cocktail_19)
summary(model_cocktail_19)

df_stats_24 <- df_stats %>%
  filter(time_hours == 24)

model_cocktail_24 <- lm(red_perc ~ Phage + Concentration, df_stats_24)
plot(model_cocktail_24)
summary(model_cocktail_24)

df_stats_48 <- df_stats %>%
  filter(time_hours == 48)

model_cocktail_48 <- lm(red_perc ~ Phage + Concentration, df_stats_48)
plot(model_cocktail_48)
summary(model_cocktail_48)


#null = no difference between CPL00272 and Cocktail, trying to consider host
model_cocktail <- lm(red_perc ~ Phage + Concentration, df_stats)
plot(model_cocktail)
summary(model_cocktail)












