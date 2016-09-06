library(devtools)
install_github("silentspringinstitute/RNHANES")

library(RNHANES)
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)

##### BPA Exposure ######
bpa_data_2011 <- nhanes_load_data("EPH_G", "2011-2012", demographics = TRUE,
                                  destination = "./nhanes_data")

nhanes_quantile(bpa_data_2011, "URXBPH", "URDBPHLC", quantiles = c(0.5, 0.95))
nhanes_detection_frequency(bpa_data_2011, "URXBPH", "URDBPHLC")

##### Download PFCs dataset #####
pfc_data_2011 <- nhanes_load_data("PFC_G", "2011-2012", demographics = TRUE,
                                  destination = "./nhanes_data")

##### Columns we would like to analyze #####
columns <- c(
  "LBXPFOA",
  "LBXPFOS",
  "LBXPFHS",
  "LBXEPAH",
  "LBXMPAH",
  "LBXPFDE",
  "LBXPFBS",
  "LBXPFHP",
  "LBXPFNA",
  "LBXPFSA",
  "LBXPFUA",
  "LBXPFDO"
)

comment_columns <- c(
  "LBDPFOAL",
  "LBDPFOSL",
  "LBDPFHSL",
  "LBDEPAHL",
  "LBDMPAHL",
  "LBDPFDEL",
  "LBDPFBSL",
  "LBDPFHPL",
  "LBDPFNAL",
  "LBDPFSAL",
  "LBDPFUAL",
  "LBDPFDOL"
)

##### Computing statistics #####
nhanes_sample_size(pfc_data_2011, columns, comment_columns, weights_column = "WTSA2YR")
nhanes_detection_frequency(pfc_data_2011, columns, comment_columns, weights_column = "WTSA2YR")
nhanes_quantile(pfc_data_2011, columns, comment_columns, weights_column = "WTSA2YR", quantiles = c(0.5))

##### Barplot of medians #####
medians <- nhanes_quantile(pfc_data_2011, columns, comment_columns, weights_column = "WTSA2YR", quantiles = c(0.5))

ggplot(medians, aes(x = column, y = value)) + geom_bar(stat = "identity")

##### Comparing medians between groups #####

#
# Let's try comparing two groups to see how their median exposure to PFCs differ.
# nhanes_quantile lets you supply a filter for calculating quantiles for a subset of the data.
# RIAGENDR is the column that indicates the gender of a participant, with 1 coded
# as male and 2 as female.
#

medians_group_a <- nhanes_quantile(pfc_data_2011,
                                columns,
                                comment_columns,
                                weights_column = "WTSA2YR",
                                quantiles = c(0.5),
                                filter = RIAGENDR == 1)
medians_group_a$group = "Males" # Record which group these medians belong to

medians_group_b <- nhanes_quantile(pfc_data_2011,
                                  columns,
                                  comment_columns,
                                  weights_column = "WTSA2YR",
                                  quantiles = c(0.5),
                                  filter = RIAGENDR == 2)
medians_group_b$group = "Females"

# Combine the data from the males and females
medians <- rbind(medians_group_a, medians_group_b)

# Comparison barchart
ggplot(medians, aes(x = column, y = value, fill = group)) +
  geom_bar(stat = "identity", position = "dodge")


##### Box Plot #####

#
# Plotting medians is a good start, but it only shows us the middle of the distribution.
# We can get a better picture of the whole distribution through a boxplot.
#
quantiles <- nhanes_quantile(pfc_data_2011,
                             columns,
                             comment_columns,
                             weights_column = "WTSA2YR",
                             quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95),
                             destination = "./nhanes_data")

# 
# In order to plot this using ggplot, we need to reformat the data returned from nhanes_quantile.
# We need one row for each chemical, with columns for the 0th, 25th, 50th, 75th, and 100th percentiles.
# This code uses dplyr and tidyr to reformat the data. There are good introductions to these packages here:
# https://cran.rstudio.com/web/packages/dplyr/vignettes/introduction.html
# and here: https://blog.rstudio.org/2014/07/22/introducing-tidyr/
#
quantiles <- quantiles %>% select(-below_lod) %>% spread(quantile, value)

ggplot(quantiles, aes(x = column, ymin = `5%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `95%`)) +
  geom_boxplot(stat = "identity")


##### Multiple Cycles #####

# Let's look at one chemical over several cycles to see how it's median level changes.

file_names <- c("L24PFC_C", "PFC_D", "PFC_E", "PFC_F", "PFC_G")
cycles <- c("2003-2004", "2005-2006", "2007-2008", "2009-2010", "2011-2012")

pfcs <- nhanes_load_data(file_names,
                         cycles,
                         demographics = TRUE,
                         destination = "./nhanes_data")

analysis <- data.frame(
  file_name = file_names,
  cycle = cycles,
  column = "LBXPFOS",
  comment_column = "LBDPFOSL",
  stringsAsFactors = FALSE
)

medians <- nhanes_quantile(pfcs, analysis, quantiles = c(0.5))

ggplot(medians, aes(x = begin_year, y = value)) + geom_bar(stat = "identity")

# Let's kick it up a notch by looking at every PFC over these 4 cycles. The analysis data frame gets more complicated:
analysis <- data.frame(
  file_name = rep(file_names, length(columns)),
  cycle = rep(cycles, length(columns)),
  column = rep(columns, each = length(cycles)),
  comment_column = rep(comment_columns, each = length(cycles)),
  stringsAsFactors = FALSE
)

medians <- nhanes_quantile(pfcs, analysis, quantiles = c(0.5))

# Problem: some of these are below the limit of detection.
ggplot(medians, aes(x = begin_year, y = value)) + geom_bar(stat = "identity") +
  facet_wrap(~column, scales = "free")


##### Correlations #####

cors <- nhanes_vcov(pfc_data_2011, column = columns) %>% cov2cor %>% melt

#
# Heatmap
#
ggplot(cors, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2()

#
# Chord diagram
#
library(circlize)
chordDiagram(cors %>% filter(Var1 != Var2, abs(value) >= 0.5))

##### Food vs. PFCs #####

food <- nhanes_load_data("DR1IFF_G", "2011-2012", demographics = TRUE, destination = "./nhanes_data")

food_summary <- food %>% group_by(SEQN) %>% summarise(calories = sum(DR1IKCAL))

pfc_data_2011_with_food <- left_join(pfc_data_2011, food_summary, by = "SEQN") %>% as.data.frame()
pfc_data_2011_with_food <- pfc_data_2011_with_food %>% filter(!is.na(calories))

nhanes_quantile(pfc_data_2011_with_food,
                column = "calories",
                comment_column = FALSE,
                quantiles = c(0.5))


nhanes_quantile(pfc_data_2011_with_food,
                column = "LBXPFOS",
                comment_column = "LBDPFOSL",
                quantiles = c(0.5),
                filter = calories >= 1999)

# hmm, there might be something interesting here. Let's make a boxplot

low_group <- nhanes_quantile(pfc_data_2011_with_food,
                             column = "LBXPFOS",
                             comment_column = "LBDPFOSL",
                             quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95),
                             filter = calories < 1999)
low_group$group <- "Low"

high_group <- nhanes_quantile(pfc_data_2011_with_food,
                             column = "LBXPFOS",
                             comment_column = "LBDPFOSL",
                             quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95),
                             filter = calories >= 1999)
high_group$group <- "High"
combined <- rbind(low_group, high_group)

combined <- combined %>% spread(quantile, value)

ggplot(combined, aes(x = group, ymin = `5%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `95%`)) +
  geom_boxplot(stat = "identity")
