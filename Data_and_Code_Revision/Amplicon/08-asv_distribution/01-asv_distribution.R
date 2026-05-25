### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set Seeds
set.seed(1994)

# Import Packages
library(ggplot2)
library(dplyr)


# Create Directory
dir_name <- "01-asv_distribution"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}
# ------------------------------------------------------------------------------


### Define Variables ------------------------------------------------------------
amplicon <- c("16S", "ITS", "Protist")
# ------------------------------------------------------------------------------


### Import Data ----------------------------------------------------------------
core_tax_all <- read.csv("../01-sort_data/02-taxonomy/All_core_ASV_taxonomy.csv", row.names = 1)
# ------------------------------------------------------------------------------


### Get Results ----------------------------------------------------------------
all_statistic_df <- all_top_df <- core_top_df <- data.frame()

# Step1: 16S
amp <- "16S"
kin <- ifelse(amp == "16S", "Bacteria", ifelse(amp == "ITS", "Fungi", "Protist"))

absolute_df <- read.csv(paste0("../01-sort_data/06-feature_abundance/", amp, "_feature_absolute_abundance.csv"), row.names = 1)
tax_df <- read.csv(paste0("../01-sort_data/02-taxonomy/", amp, "_ASV_taxonomy.csv"))

core_tax_df <- subset(core_tax_all, Clade == kin)
core_absolute_df <- absolute_df[core_tax_df$FeatureID,]

# Step1.1 Pie Plot
all_number <- nrow(absolute_df)
core_number <- nrow(core_tax_df)
all_abundance <- sum(absolute_df)
core_abundance <- sum(core_absolute_df)

statistic_df <- data.frame(
    Feature = c("Core ASVs", "Other ASVs", "Core ASVs", "Other ASVs"),
    Type = c("ASV numbers", "ASV numbers", "ASV abundance", "ASV abundance"),
    Value = c(core_number, all_number - core_number, all_abundance, all_abundance - core_abundance)
)
statistic_df <- statistic_df %>%
    group_by(Type) %>%
    mutate(
        Percent = Value / sum(Value),
    )
statistic_df$Type <- factor(statistic_df$Type, levels = c("ASV numbers", "ASV abundance"))
statistic_df$Feature <- factor(statistic_df$Feature, levels = c("Core ASVs", "Other ASVs"))

p_pie_16s <- ggplot(statistic_df, aes(x = "", y = Percent, fill = Feature)) +
    geom_bar(stat = "identity") +
    labs(
        title = paste0(kin, " (", all_number, " ASVs)")
    ) + 
    geom_text(
        aes(label = scales::percent(Percent, 0.1)),
        color = "white", position = position_stack(vjust = 0.5),
        size = 7 / 2.835
    ) +
    coord_polar(theta = "y") +
    facet_wrap(~ Type, strip.position = "bottom") +
    theme_void() +
    scale_fill_manual(values = c('#c01f25', '#57657c')) +
    theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        panel.border = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks = element_blank(), 
        plot.title = element_text(size = 7, hjust = 0.5),
        panel.spacing = unit(-1, "lines"),
        strip.text = element_text(size = 7),
        legend.title = element_text(size = 7, color = 'black'), 
        legend.text = element_text(size = 6,  color = 'black'), 
        legend.key.size = unit(0.25, 'cm'),
        legend.position = "none"
    ) 
statistic_df$Kingdom <- amp
all_statistic_df <- rbind(all_statistic_df, statistic_df)

# Step1.2: Top All ASV
all_tax_df <- tax_df["Order"]

raw_df <- data.frame(table(all_tax_df$Order))
names(raw_df) <- c("Order", "Value")
raw_df$Percent <- raw_df$Value / sum(raw_df$Value)
raw_df <- raw_df[order(-raw_df$Value),]

res_df <- raw_df[!(raw_df$Order %in% c("", "unidentified")),]

top_df <- res_df[1:10,]
top_df$Order <- as.character(top_df$Order)
top_df[top_df$Order == "0319-6G20", "Order"] <- "Oligoflexia_0319-6G20"
top_df$Kingdom <- kin

all_top_df <- rbind(all_top_df, top_df)

# Step1.3: Top Core ASV
core_tax_df <- core_tax_df["Order"]

raw_df <- data.frame(table(core_tax_df$Order))
names(raw_df) <- c("Order", "Value")
raw_df$Percent <- raw_df$Value / sum(raw_df$Value)
raw_df <- raw_df[order(-raw_df$Value),]

res_df <- raw_df[raw_df$Order != "",]

top_df <- res_df[1:10,]
top_df$Order <- as.character(top_df$Order)
top_df$Kingdom <- kin

core_top_df <- rbind(core_top_df, top_df)

# Step2: ITS
amp <- "ITS"
kin <- ifelse(amp == "16S", "Bacteria", ifelse(amp == "ITS", "Fungi", "Protist"))

absolute_df <- read.csv(paste0("../01-sort_data/06-feature_abundance/", amp, "_feature_absolute_abundance.csv"), row.names = 1)
tax_df <- read.csv(paste0("../01-sort_data/02-taxonomy/", amp, "_ASV_taxonomy.csv"))

core_tax_df <- subset(core_tax_all, Clade == kin)
core_absolute_df <- absolute_df[core_tax_df$FeatureID,]

# Step2.1 Pie Plot
all_number <- nrow(absolute_df)
core_number <- nrow(core_tax_df)
all_abundance <- sum(absolute_df)
core_abundance <- sum(core_absolute_df)

statistic_df <- data.frame(
    Feature = c("Core ASVs", "Other ASVs", "Core ASVs", "Other ASVs"),
    Type = c("ASV numbers", "ASV numbers", "ASV abundance", "ASV abundance"),
    Value = c(core_number, all_number - core_number, all_abundance, all_abundance - core_abundance)
)
statistic_df <- statistic_df %>%
    group_by(Type) %>%
    mutate(
        Percent = Value / sum(Value),
    )
statistic_df$Type <- factor(statistic_df$Type, levels = c("ASV numbers", "ASV abundance"))
statistic_df$Feature <- factor(statistic_df$Feature, levels = c("Core ASVs", "Other ASVs"))

p_pie_its <- ggplot(statistic_df, aes(x = "", y = Percent, fill = Feature)) +
    geom_bar(stat = "identity") +
    labs(
        title = paste0(kin, " (", all_number, " ASVs)")
    ) + 
    geom_text(
        aes(label = scales::percent(Percent, 0.1)),
        color = "white", position = position_stack(vjust = 0.5),
        size = 7 / 2.835
    ) +
    coord_polar(theta = "y") +
    facet_wrap(~ Type, strip.position = "bottom") +
    theme_void() +
    scale_fill_manual(values = c('#c01f25', '#57657c')) +
    theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        panel.border = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks = element_blank(), 
        strip.text = element_text(size = 7),
        plot.title = element_text(size = 7, hjust = 0.5),
        panel.spacing = unit(-1, "lines"),
        legend.title = element_text(size = 7, color = 'black'), 
        legend.text = element_text(size = 6,  color = 'black'), 
        legend.key.size = unit(0.25, 'cm'),
        legend.position = "none"
    ) 

statistic_df$Kingdom <- amp
all_statistic_df <- rbind(all_statistic_df, statistic_df)

# Step2.2: Top All ASV
all_tax_df <- tax_df["Order"]

raw_df <- data.frame(table(all_tax_df$Order))
names(raw_df) <- c("Order", "Value")
raw_df$Percent <- raw_df$Value / sum(raw_df$Value)
raw_df <- raw_df[order(-raw_df$Value),]

res_df <- raw_df[!(raw_df$Order %in% c("", "unidentified")),]

top_df <- res_df[1:10,]
top_df$Order <- as.character(top_df$Order)
top_df$Kingdom <- kin

all_top_df <- rbind(all_top_df, top_df)

# Step2.3: Top Core ASV
core_tax_df <- core_tax_df["Order"]

raw_df <- data.frame(table(core_tax_df$Order))
names(raw_df) <- c("Order", "Value")
raw_df$Percent <- raw_df$Value / sum(raw_df$Value)
raw_df <- raw_df[order(-raw_df$Value),]

res_df <- raw_df[raw_df$Order != "",]

top_df <- res_df[1:10,]
top_df$Order <- as.character(top_df$Order)
top_df$Kingdom <- kin

core_top_df <- rbind(core_top_df, top_df)


# Step3: Protist
amp <- "Protist"
kin <- ifelse(amp == "16S", "Bacteria", ifelse(amp == "ITS", "Fungi", "Protist"))

absolute_df <- read.csv(paste0("../01-sort_data/06-feature_abundance/", amp, "_feature_absolute_abundance.csv"), row.names = 1)
tax_df <- read.csv(paste0("../01-sort_data/02-taxonomy/", amp, "_ASV_taxonomy.csv"))

core_tax_df <- subset(core_tax_all, Clade == kin)
core_absolute_df <- absolute_df[core_tax_df$FeatureID,]

# Step3.1 Pie Plot
all_number <- nrow(absolute_df)
core_number <- nrow(core_tax_df)
all_abundance <- sum(absolute_df)
core_abundance <- sum(core_absolute_df)

statistic_df <- data.frame(
    Feature = c("Core ASVs", "Other ASVs", "Core ASVs", "Other ASVs"),
    Type = c("ASV numbers", "ASV numbers", "ASV abundance", "ASV abundance"),
    Value = c(core_number, all_number - core_number, all_abundance, all_abundance - core_abundance)
)
statistic_df <- statistic_df %>%
    group_by(Type) %>%
    mutate(
        Percent = Value / sum(Value),
    )
statistic_df$Type <- factor(statistic_df$Type, levels = c("ASV numbers", "ASV abundance"))
statistic_df$Feature <- factor(statistic_df$Feature, levels = c("Core ASVs", "Other ASVs"))

p_pie_pro <- ggplot(statistic_df, aes(x = "", y = Percent, fill = Feature)) +
    geom_bar(stat = "identity") +
    labs(
        title = paste0(kin, " (", all_number, " ASVs)")
    ) + 
    geom_text(
        aes(label = scales::percent(Percent, 0.1)),
        color = "white", position = position_stack(vjust = 0.5),
        size = 7 / 2.835
    ) +
    coord_polar(theta = "y") +
    facet_wrap(~ Type, strip.position = "bottom") +
    theme_void() +
    scale_fill_manual(values = c('#c01f25', '#57657c')) +
    theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        panel.border = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks = element_blank(), 
        plot.title = element_text(size = 7, hjust = 0.5),
        strip.text = element_text(size = 7),
        panel.spacing = unit(-1, "lines"),
        legend.title = element_text(size = 7, color = 'black'), 
        legend.text = element_text(size = 6,  color = 'black'), 
        legend.key.size = unit(0.25, 'cm'),
        legend.position = "none"
    ) 

statistic_df$Kingdom <- amp
all_statistic_df <- rbind(all_statistic_df, statistic_df)

name <- paste0(dir_name, "/Richness_abundance")
write.csv(all_statistic_df, paste0(name, ".csv"), quote = F, row.names = F)

p_legend <- ggplot(statistic_df, aes(x = "", y = Percent, fill = Feature)) +
    geom_bar(stat = "identity") +
    labs(
        title = paste0(kin, " (", all_number, " ASVs)")
    ) + 
    geom_text(
        aes(label = scales::percent(Percent, 0.1)),
        color = "white", position = position_stack(vjust = 0.5),
        size = 7 / 2.835
    ) +
    coord_polar(theta = "y") +
    facet_wrap(~ Type, strip.position = "bottom") +
    theme_void() +
    scale_fill_manual(values = c('#c01f25', '#57657c'), name = 'Feature') +
    theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        panel.border = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks = element_blank(), 
        plot.title = element_text(size = 7, hjust = 0.5),
        strip.text = element_text(size = 7),
        panel.spacing = unit(-1, "lines"),
        legend.title = element_text(size = 7, color = 'black'), 
        legend.text = element_text(size = 6,  color = 'black'), 
        legend.key.size = unit(0.25, 'cm')
    ) 

legend_p <- cowplot::get_legend(p_legend)


p1 <- cowplot::plot_grid(
    legend_p,
    p_pie_16s,
    p_pie_its,
    p_pie_pro,
    align = "hv", axis = "tblr", hjust = 0, vjust = 0,
    ncol = 4, nrow = 1, rel_widths = c(1, 4, 4, 4) 
)


width <- 17
height <- 5

ggsave(paste0(name, ".pdf"), p1, width = width, height = height, units = "cm")

# Step3.2: Top All ASV
all_tax_df <- tax_df["Order"]

raw_df <- data.frame(table(all_tax_df$Order))
names(raw_df) <- c("Order", "Value")
raw_df$Percent <- raw_df$Value / sum(raw_df$Value)
raw_df <- raw_df[order(-raw_df$Value),]

res_df <- raw_df[!(raw_df$Order %in% c("", "unidentified")),]

top_df <- res_df[1:10,]
top_df$Order <- as.character(top_df$Order)
top_df$Kingdom <- kin

all_top_df <- rbind(all_top_df, top_df)

# Step3.3: Top Core ASV
core_tax_df <- core_tax_df["Order"]

raw_df <- data.frame(table(core_tax_df$Order))
names(raw_df) <- c("Order", "Value")
raw_df$Percent <- raw_df$Value / sum(raw_df$Value)
raw_df <- raw_df[order(-raw_df$Value),]

res_df <- raw_df[raw_df$Order != "",]

top_df <- res_df[1:10,]
top_df$Order <- as.character(top_df$Order)
top_df[top_df$Order == "ATCC50593-Flamella-WIM80-lineage", "Order"] <- "ATCC50593-Flamella-WIM80"
top_df$Kingdom <- kin

core_top_df <- rbind(core_top_df, top_df)


# Step4: All Top
data_df <- all_top_df
data_df$Order <- factor(data_df$Order, levels = data_df$Order)
data_df[data_df$Kingdom == "Protist", "Kingdom"] <- "Protists"
data_df$Kingdom <- factor(data_df$Kingdom, levels = c("Bacteria", "Fungi", "Protists"))

p2 <- ggplot(data_df, aes(x = Order, y = Percent)) +
    geom_bar(
        stat = "identity",
        position = "stack",
        width = 0.68,
        fill = '#57657c'
    ) +
    labs(
        x = '',
        y = 'Proportion of all ASVs'
    ) + 
    scale_y_continuous(labels = scales::percent) +
    
    facet_grid(
        ~ Kingdom,
        scales = "free_x",
        switch = NULL,
        space = "fixed",
    ) +
    theme_bw() + theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 7, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
        panel.spacing = unit(0.1, "cm"),
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right"
    )

name <- paste0(dir_name, "/Relative_ASV_richness")
width <- 17
height <- 7
ggsave(paste0(name, ".pdf"), p2, width = width, height = height, units = "cm")
write.csv(data_df, paste0(name, ".csv"), row.names = F)

# Step5: Core Top
data_df <- core_top_df
data_df$Order <- factor(data_df$Order, levels = data_df$Order)
data_df[data_df$Kingdom == "Protist", "Kingdom"] <- "Protists"
data_df$Kingdom <- factor(data_df$Kingdom, levels = c("Bacteria", "Fungi", "Protists"))

p3 <- ggplot(data_df, aes(x = Order, y = Percent)) +
    geom_bar(
        stat = "identity",
        position = "stack",
        width = 0.68,
        fill = '#c01f25'
    ) +
    labs(
        x = '',
        y = 'Proportion of core ASVs'
    ) + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    
    facet_grid(
        ~ Kingdom,
        scales = "free_x",
        switch = NULL,
        space = "fixed",
    ) +
    
    theme_bw() + theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 7, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
        panel.spacing = unit(0.1, "cm"),
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right"
    )

name <- paste0(dir_name, "/Taxonomic_composition_of_core_ASV")
width <- 17
height <- 7.4
ggsave(paste0(name, ".pdf"), p3, width = width, height = height, units = "cm")
write.csv(data_df, paste0(name, ".csv"), row.names = F)
# ------------------------------------------------------------------------------
