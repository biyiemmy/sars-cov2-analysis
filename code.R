library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
install.packages("ggpubr")
install.packages("RColorBrewer")

library(ggpubr)
library(RColorBrewer)

var <- read.csv("variants_long_table.csv")
head(var)

# Check the dimension of the data
dim(var)

# Check number of rows
nrow(var)

# Check number of columns
ncol(var)

# Display the structure of your R object
str(var)

# Summary statistics of the whole data or specified columns
## For the whole table 
summary(var)

## For non-numerical data
summary(var$SAMPLE)

## For numerical data
summary(var$DP)

# Check the class of your data
class(var)

# Check the class of an object
class(var$CHROM)

typeof(var$CHROM)

# Preview the data using a spreadsheet-style data viewer in RStudio
View(var)


# Check columns name
colnames(var)

# Select columns 1, 4, and 5
var[, c(1, 4, 5)]

# Select columns 1, 4, and 5 with default display
select(var, SAMPLE, REF, ALT)

# Select columns 1, 4, and 5 with selected display
select(var, SAMPLE, REF, ALT) %>% head(3)

# Select all columns except the column “CALLER” with selected display
select(var, -CALLER) %>% head(3)

# Transform the data frame into a tibble
var_tb <- as_tibble(var)
select(var_tb, SAMPLE, REF, ALT) %>% head(3)

# Select rows with selected display using base R code
var_tb[var_tb$SAMPLE == "SRR13500958",]

# Select rows with selected display using dplyr functions
filter(var_tb, SAMPLE == "SRR13500958") %>% head(3)

# Select sample type (rows) and variables (columns) with selected display
var_tb %>% filter(SAMPLE == "SRR13500958") %>% select(CHROM, POS, REF, ALT) %>% head(3)

# To select all data related to the sample specified
var_tb %>% filter(SAMPLE == "SRR13500958") %>% select(CHROM, POS, REF, ALT, DP)

# To select only values for which DP>=500 for the same sample
var_tb %>% filter(SAMPLE == "SRR13500958" & DP>=500) %>% select(CHROM, POS, REF, ALT, DP)

# To select only values for which DP>=1000 for the same sample
var_tb %>% filter(SAMPLE == "SRR13500958" & DP>=1000) %>% select(CHROM, POS, REF, ALT, DP)

# Count how many rows are associated with each sample in the data 
var_tb %>% count(SAMPLE)

# Sorting the counts 
var_tb %>% count(SAMPLE, sort = TRUE)

# Distribution of genes per sample and counts 
var_tb %>% count(SAMPLE, GENE, sort = TRUE) %>% head()

# Maximum value of column DP
max(var_tb$DP)

# Minimum value of column DP
min(var_tb$DP)

# Mean value of column DP
mean(var_tb$DP)

# Compute a LOG2 transformation on the DP values
var_tb_log <- var_tb %>% mutate(DP_log2 = log2(DP))

# View the table columns with the DP_log2 new column appended at the end
head(var_tb_log)

# View a selected content including the new column
select(var_tb_log, SAMPLE, REF, ALT, DP, DP_log2) %>% head()

# Show the maximum value of DP for each sample
var_tb %>% group_by(SAMPLE) %>% summarize(max(DP))

# Show the minimum value of DP for each sample
var_tb %>% group_by(SAMPLE) %>% summarize(min(DP))

# Link ggplot2 to a specific data frame
ggplot(data = var_tb)

# Link ggplot2 to specific variables using aesthetics
ggplot(data = var_tb, aes(x=SAMPLE, y=DP))


# Points (left-hand plot)
ggplot(data = var_tb, aes(x=SAMPLE, y=DP)) + geom_point()

# Boxplot (right-hand plot)
ggplot(data = var_tb, aes(x=SAMPLE, y=DP)) + geom_boxplot()


# Points (left-hand plot)
ggplot(data = var_tb, aes(x=SAMPLE, y=DP)) + geom_point() + ylim(0,10000)

# Boxplot (right-hand plot)
ggplot(data = var_tb, aes(x=SAMPLE, y=DP)) + geom_boxplot() + ylim(0,10000)

# Points (left-hand plot)
ggplot(data = var_tb, aes(x=SAMPLE, y=DP)) + geom_point() + scale_y_continuous(name="dp", limits=c(0, 10000))

# Boxplot (right-hand plot)
ggplot(data = var_tb, aes(x=SAMPLE, y=DP)) + geom_boxplot() + scale_y_continuous(name="dp", limits=c(0, 10000))


# Points (left-hand plot)
ggplot(data = var_tb, aes(x=SAMPLE, y=DP)) + geom_point() + scale_y_continuous(trans='log10')
ggplot(data = var_tb, aes(x=SAMPLE, y=DP)) + geom_point() + scale_y_log10()

# Boxplot (right-hand plot)
ggplot(data = var_tb, aes(x=SAMPLE, y=DP)) + geom_boxplot() + scale_y_continuous(trans='log10')
ggplot(data = var_tb, aes(x=SAMPLE, y=DP)) + geom_boxplot() + scale_y_log10()

# Colours of shapes
ggplot(data = var_tb, aes(x=SAMPLE, y=DP, colour = SAMPLE)) + geom_boxplot() + ylim(0,10000)

# Colours for filling options
ggplot(data = var_tb, aes(x=SAMPLE, y=DP, fill= SAMPLE)) + geom_boxplot() + ylim(0,10000)

# Colours for filling options with manual colors
ggplot(data = var_tb, aes(x=SAMPLE, y=DP, fill= SAMPLE)) + geom_boxplot() + ylim(0,10000) + scale_fill_manual(values=c("#cb6015", "#e1ad01", "#6d0016", "#808000", "#4e3524"))

# Colours for filling options with preset palettes
install.packages("RColorBrewer")
library(RColorBrewer)
ggplot(data = var_tb, aes(x=SAMPLE, y=DP, fill= SAMPLE)) + geom_boxplot() + ylim(0,10000) + scale_fill_brewer(palette="RdYlBu")

# All possible palettes can be displayed using:
display.brewer.all()

# Change legend position
ggplot(data = var_tb, aes(x=SAMPLE, y=DP, fill= SAMPLE)) + geom_boxplot() + ylim(0,10000) + scale_fill_brewer(palette="RdYlBu") + theme(legend.position="top")
ggplot(data = var_tb, aes(x=SAMPLE, y=DP, fill= SAMPLE)) + geom_boxplot() + ylim(0,10000) + scale_fill_brewer(palette="RdYlBu") + theme(legend.position="none")

# Change shapes, colors, and sizes
ggplot(data = var_tb, aes(x=SAMPLE, y=DP)) +  geom_point(shape = 21, fill = "#e4dbc1", color = "#b92e17", size = 6) + ylim(0,10000)
ggplot(data = var_tb, aes(x=SAMPLE, y=DP)) +  geom_point(shape = 23, color = "#e4dbc1", fill = "#b92e17", size = 5, alpha=0.5) + ylim(0,10000)

# All possible points can be displayed using:
ggpubr::show_point_shapes()


# Progressing in variants exploration
# View the data 
View(var_tb)

# Distribution of DP values per chromosome and per sample
ggplot(data = var_tb, aes(x=CHROM, y=DP, fill= SAMPLE)) + geom_boxplot() + ylim(0,10000) + scale_fill_brewer(palette="RdYlBu") + labs(title="DP_per_Chromosome") + facet_grid(. ~ SAMPLE)

# Define a variable with plotting options
p_DP_CHROM <- ggplot(data = var_tb, aes(x=CHROM, y=DP, fill= SAMPLE)) + ylim(0,10000) + scale_fill_brewer(palette="RdYlBu") + labs(title="DP_per_Chromosome") + theme(legend.position="bottom")

# Test boxplots with faceting 
p_DP_CHROM + geom_boxplot() + facet_grid(. ~ SAMPLE)

# Combine violin plots and boxplots with faceting
p_DP_CHROM + geom_violin(trim=FALSE) + facet_grid(. ~ SAMPLE) + geom_boxplot(width=0.1)

# Variants effects per sample
# 1. Plotting the variants effects

# Count number of different effects per sample
p_EFFECT <- ggplot(data = var_tb, aes(x=EFFECT, fill= SAMPLE)) + scale_fill_brewer(palette="RdBu") + labs(title="Effect_per_Sample") + theme(legend.position="bottom")
p_EFFECT + geom_bar()

# Flip orientation
p_EFFECT_flip <- ggplot(data = var_tb, aes(y=EFFECT, fill= SAMPLE)) + scale_fill_brewer(palette="RdBu") + labs(title="Effect_per_Sample") + theme(legend.position="bottom")
p_EFFECT_flip + geom_bar()

# 2. Counting the variants effects
# Count the number of different effects
var_tb %>% count(EFFECT)

# Count the number of different effects and link them to sample information
var_tb %>% count(EFFECT, SAMPLE, sort = TRUE)

# The genes showing variants effects
# 1. Counting and extracting all effects for all genes

# Counting the effects per gene
var_tb %>% count(EFFECT, GENE, sort = TRUE)

# 2. Counting and extracting specific effects for all genes
# Filtering option 1 to select for effect on stop
filter(var_tb, EFFECT == "stop_lost" | EFFECT == "stop_gained")

# Filtering option 2 to select for effect on stop
filter(var_tb, EFFECT %in% c("stop_lost", "stop_gained"))

# Filtering on effect and selected columns
filter(var_tb, EFFECT %in% c("stop_lost", "stop_gained")) %>% select(SAMPLE, CHROM, GENE, EFFECT)

# Read depth per position
# Define your variable
p_DP_POS <- ggplot(data = var_tb, aes(x=POS, y=DP, fill= SAMPLE)) + scale_fill_brewer(palette="RdBu") + labs(title="DP_per_Position") + theme(legend.position="bottom")

# Plot
p_DP_POS + geom_point(shape = 21, size = 5)

# Plot with transparency options
p_DP_POS + geom_point(shape = 21, size = 5, alpha = 0.7)

# Using what you learned on ggplot2 syntax, scatter plots, and the var_tb data, plot the depth (DP) over alternative base depth (ALT_DP). No further formatting is required at this stage.
p_DP_ALTdp <- ggplot(data=var_tb, aes(x=ALT_DP,y=DP))+labs(title="DP over ALT_DP")+theme(legend.position="bottom")
p_DP_ALTdp + geom_point()

# Using the var_tb data, plot the depth (DP) over alternative base depth (ALT_DP), while showing the samples in different colors, with the following points characteristics: (shape = 21, size = 5). Please use the palette “RdBu”.
p_DP_ALTdp <- ggplot(data=var_tb, aes(x=ALT_DP,y=DP, fill=SAMPLE))+scale_fill_brewer(palette="RdBu")
p_DP_ALTdp + geom_point(shape=21, size=5)

# Create a new variable called “test”, and place in it the code you just created. Then add to this code a title (title=”DP_per_ALT_DP”) and place the legend at the bottom (legend.position=”bottom”)
test <- ggplot(data=var_tb, aes(x=ALT_DP,y=DP, fill=SAMPLE))+scale_fill_brewer(palette="RdBu")+labs(title="DP over ALT_DP")+theme(legend.position="bottom")
test + geom_point(shape=21, size=5)

