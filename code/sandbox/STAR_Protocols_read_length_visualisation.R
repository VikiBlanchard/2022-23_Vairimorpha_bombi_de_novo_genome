#!/usr/bin/env Rscript
# Please specify your working directory using setwd
setwd("/path/to/Input_CSV_file")
library(ggplot2)
library(dplyr)
library(cowplot)
# Import the read-length distribution table
read_length_df <- read.csv("length.csv")
# Organize the imported read-length table
# You can replace the level arguments for your platform, species, or strains
read_length_df$platform <- as.factor(read_length_df$platform)
read_length_df$platform <- factor(read_length_df$platform,level = c("PacBio_CLR",
"PacBio_HiFi","ONT"))
# Calculate the average read-lengths for each platform
summary_df <- ddply(read_length_df, "platform", summarise, grp.mean=mean(length))
# Draw a read-length distribution plot for all reads
total.length.plot <- ggplot(read_length_df, aes(x=length, fill=platform, color=plat
form)) +
geom_histogram(binwidth=100, alpha=0.5, position="dodge") +
geom_vline(data=summary_df, aes(xintercept=grp.mean, color=platform), linetype="-
dashed", size =0.2) +
scale_x_continuous(labels = comma) +
scale_y_continuous(labels = comma) +
labs(x = "Read length (bp)", y = "Count") +
theme_bw()
# Draw a read-length distribution plot for reads % 20 kb in length
20 kb.length.plot <- ggplot(read_length_df, aes(x=length, fill=platform, color=platform)) +
geom_histogram(binwidth=50, alpha=0.5, position="dodge") +
geom_vline(data=summary_df, aes(xintercept=grp.mean, color=platform), linetype=
"dashed", size=0.2) +
scale_x_continuous(labels = comma, limit = c(0,20000)) +
scale_y_continuous(labels = comma) +
labs(x = "Read length (bp)", y = "Count") +
theme_bw()
# Merge both the read-length distribution plots
plot <- plot_grid(total.length.plot, 20 kb.length.plot, ncol = 1)
# Save the figure using the file name, ‘‘read.length.pdf’’
pdf("read.length.pdf",width=6,height=8,paper=’special’)
print(plot)
dev.off()