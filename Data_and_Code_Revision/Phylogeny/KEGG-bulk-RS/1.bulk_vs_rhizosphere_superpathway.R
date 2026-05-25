pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)
library(tidyverse)
library(patchwork) 

kingdoms <- c("Bacteria", "Fungi", "Protist")

metadata <- read.csv("./data/metadata.csv", check.names = F)
Superpathway_color <- read.csv("./data/Superpathway_color.csv", check.names = F)
Superpathway_color <- Superpathway_color[Superpathway_color$Superpathway != "Unassign", ]
color_map <- setNames(Superpathway_color$Color, Superpathway_color$Superpathway)

# 存储所有kingdom的绘图列表
all_plot_list <- list()

for (kingdom in kingdoms) {
  ## Import Data
  bs_df_path <- paste0("./data/",kingdom,"_bulk_tpm_superpathway_abundance.csv")
  rs_df_path <- paste0("./data/",kingdom,"_RS_tpm_superpathway_abundance.csv")
  
  bs_df <- read.csv(bs_df_path, check.names = F, row.names = 1)
  rs_df <- read.csv(rs_df_path, check.names = F, row.names = 1)
  
  # 去掉Unassign
  bs_df <- bs_df[row.names(bs_df) != "Unassign", ]
  rs_df <- rs_df[row.names(rs_df) != "Unassign", ]
  
  ## Data Processing
  bs_core <- gsub("^S-", "", colnames(bs_df))
  names(bs_core) <- colnames(bs_df)
  rs_core <- gsub("^RS-", "", colnames(rs_df))
  names(rs_core) <- colnames(rs_df)
  common_ids <- intersect(bs_core, rs_core)
  
  bs_filtered <- bs_df[, names(bs_core)[bs_core %in% common_ids], drop = FALSE]
  rs_filtered <- rs_df[, names(rs_core)[rs_core %in% common_ids], drop = FALSE]
  
  all_rows <- union(rownames(bs_filtered), rownames(rs_filtered))
  bs_complete <- matrix(0, nrow = length(all_rows), ncol = ncol(bs_filtered),
                        dimnames = list(all_rows, colnames(bs_filtered)))
  bs_complete[rownames(bs_filtered), ] <- as.matrix(bs_filtered)
  bs_complete <- as.data.frame(bs_complete)
  
  rs_complete <- matrix(0, nrow = length(all_rows), ncol = ncol(rs_filtered),
                        dimnames = list(all_rows, colnames(rs_filtered)))
  rs_complete[rownames(rs_filtered), ] <- as.matrix(rs_filtered)
  rs_complete <- as.data.frame(rs_complete)
  
  combined_df <- cbind(bs_complete, rs_complete)
  
  long_df <- combined_df %>%
    rownames_to_column("Superpathway") %>%
    pivot_longer(-Superpathway, names_to = "SampleID", values_to = "Abundance") %>%
    mutate(CoreID = gsub("^S-|^RS-", "", SampleID)) %>%
    mutate(SampleType = ifelse(grepl("^S-", SampleID), "Bulk", "Rhizo"))
  
  # 按CoreID排序，确保对应编号位置一致
  core_order <- sort(unique(long_df$CoreID))
  long_df$CoreID <- factor(long_df$CoreID, levels = core_order)
  long_df$SampleID <- factor(long_df$SampleID, 
                             levels = c(paste0("S-", core_order), paste0("RS-", core_order)))
  
  long_df$Superpathway <- factor(
    long_df$Superpathway,
    levels = Superpathway_color$Superpathway
  )
  
  ## Plot
  p <- ggplot(long_df, aes(x = SampleID, y = Abundance, fill = Superpathway)) +
    geom_col(position = "fill", width = 0.8) + 
    facet_wrap(~SampleType, scales = "free_x", nrow = 1) +
    theme_bw() +
    scale_fill_manual(values = color_map) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.title.x = element_text(size = 8),  
      axis.title.y = element_text(size = 8), 
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 8, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
      panel.spacing = unit(0.5, "lines")
    ) +
    labs(
      x = "",
      y = "Relative abundance",
      title = kingdom
    ) +
    scale_x_discrete(expand = c(0.02, 0.02))
  
  all_plot_list[[kingdom]] <- p
}

combined_plot <- wrap_plots(all_plot_list, ncol = 1, heights = c(1, 1, 1)) +
  plot_layout(guides = 'collect') &
  guides(fill = guide_legend(ncol = 3)) &
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.2, "cm"),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
  )


# 保存图片
ggsave("combined_superpathway_abundance.png", combined_plot, 
       width = 14, height = 27, dpi = 600, units = "cm")
ggsave("combined_superpathway_abundance.pdf", combined_plot, 
       width = 14, height = 27, units = "cm")

# 显示图形
print(combined_plot)