library(tidyverse)
library(cowplot)
library(DescTools)

#read in the synogram data
row_content <- read.csv("row_data.csv")
col_content <- read.csv("column_data.csv")

#This function will reformat an individual .CSV into the correct format by modifying the columns, reformatting the time
#adding in the dictionary that tells us what is in each well, and removes NAs for dodgy values
reformat <- function(df_test) {
  df_test <- df_test %>% select(-Well)
  df_test[1,1] <- "Well"
  #names(df_test) <- lapply(df_test[1, ], as.character)
  df_test <- df_test[-1,] 
  df_test <- rename(df_test, "Time_Hours" = "X")
  
  df_test <- pivot_longer(df_test, cols = -Time_Hours, names_to = "Well", values_to = "values") %>%
    filter(Time_Hours != "Time")
  
  df_test <- df_test %>%
    mutate(
      hours = as.numeric(str_extract(Time_Hours, "\\d+(?=\\s*h)")),
      minutes = as.numeric(str_extract(Time_Hours, "\\d+(?=\\s*min)")),
      hours = replace_na(hours, 0),
      minutes = replace_na(minutes, 0),
      Time_seconds = (hours * 3600) + (minutes * 60)
    ) %>%
    select(-hours, -minutes)  # Remove temporary columns
  
  df_test <- df_test %>%
    mutate(
      Row = str_extract(Well, "^[A-Z]"),
      Column = as.numeric(str_extract(Well, "\\d+$"))
    )
}


#This function now will take each of the filenames in the folder, and for each it'll take in
#the file name, run the reformatting function above, add in a column with the plate name
#(taken from the file name) 

produce_df <- function(file_name){
  df <- read.csv(file_name, header = TRUE, row.names = NULL)
  df <- reformat(df)
  file_base_name <- sub("\\.\\w+$", "", basename(file_name))
  df$plate <- file_base_name
  return(df)
}

#get the list of all the csv files in the folder
csv_files <- list.files(
  pattern = "\\.CSV$"
)

#use purrr::map_df to apply the produce_df function to each CSV file and bind them together
df_all <- map(csv_files, produce_df)

#df in list form, this isnt ideal, need a concatenated one
df_list <- function(df_list){
  final_df <- Reduce(full_join, df_list)
  return(final_df)
}
df_all <- df_list(df_all)

df_all <- left_join(df_all, row_content) %>% left_join(col_content) %>% na.omit() %>%
  mutate(Row = factor(Row, levels = c("A","B","C","D","E","F","G"))) %>%
  mutate(Column = factor(Column, levels = seq(2,10,1)))

temporary_plot <- ggplot(df_all, aes(x = Time_seconds, y = values, group = Well)) +
  geom_line() +
  facet_grid(plate+Row~Column)

#plate 7 definitely contaminated. Will have to re-run. The rest are okay though
df_all_analysis <- df_all %>%
  filter(plate != "Plate7") %>%
  mutate(values = as.numeric(values)) %>%
  filter(Time_Hours != "")

df_blanks <- df_all_analysis %>%
  filter(bacteria == "None") %>%
  filter(phage == "No Phage") %>%
  group_by(Time_seconds, plate) %>%
  mutate(mean_blk_od = mean(values)) %>%
  select(Time_seconds, plate, mean_blk_od)

df_all_analysis <- df_all_analysis %>%
  left_join(df_blanks) %>%
  mutate(norm_od = values - mean_blk_od) %>%
  filter(bacteria != "None") %>%
  mutate(Time_Hours = Time_seconds/(60*60)) %>%
  mutate(phage = factor(phage, levels = c("No Phage", "CPL00272", "PNM", "Cocktail"))) 

gc_plot <- ggplot(df_all_analysis, aes(x = Time_Hours, y = abs(norm_od), group = Well, colour = phage)) +
  geom_line() +
  facet_wrap(~bacteria, nrow = 3) +
  scale_x_continuous("Time (Hours)") +
  scale_y_continuous("Normalised OD600") +
  scale_colour_manual("Phage", values = okabe) +
  theme_cowplot(16)

control_auc_df <- df_all_analysis %>%
  filter(phage == "No Phage") %>%
  group_by(phage, bacteria, Well) %>%
  mutate(control_auc = DescTools::AUC(Time_Hours, norm_od, method = "trapezoid")) %>%
  group_by(bacteria) %>%
  reframe(mean_cntrl_auc = mean(control_auc)) %>%
  select(bacteria, mean_cntrl_auc)

virulence_df <- df_all_analysis %>%
  group_by(phage, bacteria, Well) %>%
  mutate(auc = DescTools::AUC(Time_Hours, norm_od, method = "trapezoid")) %>%
  select(Well, phage, bacteria, auc) %>%
  unique() %>%
  left_join(control_auc_df) %>%
  mutate(vp = 1-(auc/mean_cntrl_auc)) %>%
  mutate(phage = factor(phage, levels = c("No Phage", "CPL00272", "PNM", "Cocktail")))

virulence_means <- virulence_df %>%
  group_by(phage, bacteria) %>%
  reframe(mean_vp = mean(vp), stdev = sd(vp))
  
virulence_means_plotting <- virulence_means %>%
  filter(phage != "No Phage") %>%
  mutate(does_infect = ifelse(mean_vp > 0.499, "Inhibition", "No Inhibition")) %>%
  mutate(does_infect = factor(does_infect, levels = c("Inhibition", "No Inhibition"))) %>%
  mutate(mean_vp = ifelse(mean_vp < 0, 0, 
                          ifelse(mean_vp > 1, 1, mean_vp)))

host_range_pnm_272 <- ggplot(virulence_means_plotting, aes(x = phage, y = bacteria, fill = mean_vp)) +
  geom_tile() +
  scale_fill_gradient2("Mean Virulence", high = "#a002a2", mid = "#e28dbe", low = "white", midpoint = 0.5) +
  scale_y_discrete('Pseudomonas aeruginosa Strain') +
  scale_x_discrete("Phage") +
  labs(fill = "Inhibition") +
  theme_classic(base_size=16) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black")) +
  theme(axis.text.x = element_text(angle = 90))
host_range_pnm_272

only_infect <- ggplot(filter(virulence_means_plotting, does_infect == "Inhibition"), aes(x = phage, y = mean_vp)) +
  geom_boxplot() +
  geom_point(position = position_dodge2(width = 0.4)) +
  scale_y_continuous("Mean Virulence Index", limits = c(0,1.5), breaks = seq(0,1.5,0.2)) +
  scale_x_discrete("Phage") +
  scale_colour_manual("Inhibtion\nover 48h", values = okabe) +
  theme_cowplot(16)

quick_test <- virulence_means_plotting %>%
  filter(does_infect == "Inhibition") %>%
  select(bacteria, phage, mean_vp) %>% 
  pivot_wider(names_from = phage, values_from = mean_vp)

t.test(quick_test$PNM, quick_test$Cocktail)
  
right <- plot_grid(host_range_pnm_272, only_infect, nrow = 2, rel_heights = c(2,1), labels = c("B", "C"))

plot_grid(gc_plot, right, nrow = 1, rel_widths = c(3,1), labels = c("A", ""))  
  
  
  
  
  















