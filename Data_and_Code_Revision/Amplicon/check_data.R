### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)


### Define variable -----------------------------------------------------------
kingdom <- c("16S", "ITS", "Protist")
# ------------------------------------------------------------------------------


### Sort data ------------------------------------------------------------------
# 1. Reads Numbers
for (i in kingdom) {
    df <- read.csv(paste0("00-rawdata/data_stats/", i, "_stats_spike.csv"), header = T)
    print(paste0("=== ", i, " (", round(sum(df$Input) / 1000000), ") ==="))
}

# 2. Feature Numbers
for (i in kingdom) {
    df <- read.csv(paste0("01-sort_data/02-taxonomy/", i, "_ASV_taxonomy.csv"), row.names = 1)
    print(paste0("=== ", i, " (", nrow(df), ") ==="))
}

# 3. Relative Abundance
# 3.1 Class
check_list <- NULL
check_list[["16S_Class"]] <- c("Gammaproteobacteria", "Alphaproteobacteria", "Actinobacteria")
check_list[["ITS_Class"]] <- c("Sordariomycetes", "Dothideomycetes", "Eurotiomycetes")
check_list[["Protist_Class"]] <- c("Gregarinomorphea", "Filosa-Sarcomonadea", "Colpodea")

for (i in kingdom) {
    df <- read.csv(paste0("01-sort_data/04-barplot_data/", i, "_barplot_data_relative_sorted.csv"))
    for (j in check_list[[paste0(i, "_Class")]]) {
        value <- df[df$Level == "Class" & df$Taxonomy == j, "Mean"]
        print(paste0("=== ", j, " (", round(value * 100, 2), ") ==="))
    }
}
# 3.1 Order
check_list <- NULL
check_list[["16S_Order"]] <- "Rhizobiales"
check_list[["ITS_Order"]] <- "Pleosporales"
check_list[["Protist_Order"]] <- "Eugregarinorida"
for (i in kingdom) {
    df <- read.csv(paste0("01-sort_data/04-barplot_data/", i, "_barplot_data_relative_sorted.csv"))
    j <- check_list[[paste0(i, "_Order")]]
    value <- df[df$Level == "Order" & df$Taxonomy == j, "Mean"]
    print(paste0("=== ", j, " (", round(value * 100, 2), ") ==="))
}

# 4. Unidentified Genus ASV
for (i in kingdom) {
    df <- read.csv(paste0("01-sort_data/02-taxonomy/", i, "_ASV_taxonomy.csv"))
    value <- nrow(df[df$Genus %in% c("", "Unknown_Family", "uncultured", "unidentified"),]) / nrow(df)
    print(paste0("=== ", i, " (", round(value * 100, 2), ") ==="))
    
    if (i == "Protist") {
        value <- nrow(df[df$Supergroup == "",]) / nrow(df)
        print(paste0("=== ", i, " SuperGroup (", round(value * 100, 2), ") ==="))
    }
}

# 5. ASV Numbers and Relative Abundance
check_list <- NULL
check_list[["16S_Order"]] <- c("Rhizobiales", "Pseudomonadales", "Chlamydiales")
i <- "16S"

# 5.1 Relative Abundance
df <- read.csv(paste0("01-sort_data/04-barplot_data/", i, "_barplot_data_relative_sorted.csv"))
for (j in check_list[[paste0(i, "_Order")]]) {
    value <- df[df$Level == "Order" & df$Taxonomy == j, "Mean"]
    print(paste0("=== ", j, " (", round(value * 100, 2), ") ==="))
}

# 5.2 ASV Numbers
df <- read.csv(paste0("01-sort_data/02-taxonomy/", i, "_ASV_taxonomy.csv"))
for (j in check_list[[paste0(i, "_Order")]]) {
    value <- nrow(df[df$Order == j,]) / nrow(df)
    print(paste0("=== ", j, " (", round(value * 100, 2), ") ==="))
}

# 6. Absolute Abundance
df <- read.csv("09-absolute_abundance/01-absolute_abundance/Sample_absolute_abundance.csv", header = T)
for (i in c("16S", "ITS", "18S")) {
    sub_df <- subset(df, Group == i)
    value <- mean(sub_df$Value)
    print(paste0("=== ", i, " (", format(value, scientific = TRUE, digits = 3), ") ==="))
}
# ------------------------------------------------------------------------------

