# README

This directory contains the R scripts and real example data used to generate the figures and tables reported in the manuscript. It is a code/scripts submission, not a standalone software package. The Nature Research “software/code” checklist is addressed below.

---

## 1. System requirements

### 1.1 Tested environment

| Item | Version / configuration |
|---|---|
| Operating system | Microsoft Windows 11 Home (Chinese edition), version 10.0.26200, Build 26200, 64-bit |
| R | 4.5.2 |
| RStudio | Recommended; scripts use `rstudioapi::getSourceEditorContext()` to detect their own path. Tested with `rstudioapi` 0.18.0. |
| CPU used for testing | Intel Core i7-12700, 12 cores / 20 threads |
| Memory used for testing | 64 GB RAM |
| Non-standard hardware | None required. A standard 64-bit desktop is sufficient. |

Only the configuration above was tested. Compatibility with other operating systems or R versions has not been evaluated.

### 1.2 Software dependencies

All analyses were performed in R. The following package versions were installed and tested under R 4.5.2.

| R package | Tested version |
|---|---:|
| `ComplexHeatmap` | 2.26.1 |
| `GUniFrac` | 1.9 |
| `RColorBrewer` | 1.1.3 |
| `ReporterScore` | 0.2.5 |
| `Vennerable` | 3.1.0.9000 |
| `agricolae` | 1.3.7 |
| `ape` | 5.8.1 |
| `circlize` | 0.4.17 |
| `clusterProfiler` | 4.18.4 |
| `cowplot` | 1.2.0 |
| `devtools` | 2.5.0 |
| `dplyr` | 1.2.0 |
| `geosphere` | 1.6.5 |
| `ggnewscale` | 0.5.2 |
| `ggplot2` | 4.0.2 |
| `ggplotify` | 0.1.3 |
| `ggpubr` | 0.6.3 |
| `ggraph` | 2.2.2 |
| `ggtree` | 4.0.5 |
| `ggtreeExtra` | 1.20.1 |
| `grid` | 4.5.2 |
| `gridExtra` | 2.3 |
| `igraph` | 2.2.2 |
| `magrittr` | 2.0.4 |
| `microeco` | 2.0.0 |
| `pROC` | 1.19.0.1 |
| `pacman` | 0.5.1 |
| `patchwork` | 1.3.2 |
| `phytools` | 2.5.2 |
| `picante` | 1.8.2 |
| `purrr` | 1.2.1 |
| `randomForest` | 4.7.1.2 |
| `rbiom` | 2.2.1 |
| `reshape2` | 1.4.5 |
| `rstatix` | 0.7.3 |
| `rstudioapi` | 0.18.0 |
| `scales` | 1.4.0 |
| `stats` | 4.5.2 |
| `stringr` | 1.6.0 |
| `tibble` | 3.3.1 |
| `tidyr` | 1.3.2 |
| `tidyverse` | 2.0.0 |
| `tools` | 4.5.2 |
| `treeio` | 1.34.0 |
| `varSelRF` | 0.7.9 |
| `vegan` | 2.7.3 |

---

## 2. Installation guide

1. Install R 4.5.2 from CRAN.
2. Install a recent RStudio release (recommended, not strictly required).
3. Install the R packages listed in Section 1.2. CRAN packages can be installed with `install.packages()`; Bioconductor packages such as `ComplexHeatmap`, `clusterProfiler`, and `circlize` should be installed with `BiocManager::install()`. For any package not available from CRAN or Bioconductor, install it from its official source at the version listed above.
4. Unzip this directory while preserving the folder structure. All scripts use relative paths and expect to be run from their own module folder.

Typical installation time on a normal desktop with a broadband connection is approximately 15–30 min: about 3–5 min for R, 2–3 min for RStudio, and 10–20 min for the R packages, depending on network speed.

---

## 3. Demo

### 3.1 Demo dataset

A separate simulated dataset is not provided because this is a code/scripts submission. The real example data included in this directory are the demo dataset.

### 3.2 Instructions to run on the demo data

The quickest demo is the tree/barplot script:

1. Open RStudio.
2. Open `tree_barplot/1.tree_microbiome_barplot.R`.
3. Click **Source**, or run the whole script.

The script detects its own location through `rstudioapi` and sets the working directory automatically. If running outside RStudio, first set the working directory to the script folder and replace the `rstudioapi::getSourceEditorContext()` line with `pwd <- getwd()`.

### 3.3 Expected output

The demo script regenerates the same PDF/PNG figures and CSV result tables that are already present in the corresponding module folder. The complete set of regenerated files is provided in `../check/`, which mirrors the structure of this directory exactly and can be used for side-by-side comparison.

### 3.4 Expected run time

On the tested desktop, the demo script `tree_barplot/1.tree_microbiome_barplot.R` takes about 11 s. Other quick scripts take approximately 1–30 s, medium-sized scripts take approximately 1–10 min, and permutation-based, random-forest, or phylogenetic-signal scripts can take several hours each. A single complete pass over all 57 scripts took approximately 7 h; allowing the time-consuming scripts to run longer took approximately 14 h in total.

### 3.5 Large raw inputs intentionally not included

The following large raw input files are not included. Their published downstream result tables and figures are retained, so the manuscript results are fully represented. Re-running the two upstream scripts that require these raw files is possible only after the original raw files are supplied:

- `kegg_abundance/absolute/reads/Bacteria_KO_number.tsv`
- `kegg_abundance/absolute/reads/Fungi_KO_number.tsv`
- `kegg_abundance/absolute/reads/Protist_KO_number.tsv`
- `funguild/rhizosphere_dada2_sequences_blast_with_Fungi_no_spike_taxonomy.tsv`
- `funguild/rhizosphere_dada2_table_blast_with_Fungi_no_spike.tsv`

---

## 4. Instructions for use

### 4.1 Running the code on your own data

Replace the input files in `metadata/`, `feature_table/`, and the module-specific `data/` folders with files that have the same column names and sample identifiers. Keep the original file names where possible; otherwise update the corresponding `read.csv()` or `read.delim()` paths. Then run the scripts module by module, following the order in the mapping table below. All outputs are written back into the corresponding module folders.

For `distance_decay/`, run `1.distance_decay-ASV.R` and `1.distance_decay_for_all_orders-ASV.R` in the same R session and in that order, because the second script uses distance matrices created by the first.

### 4.2 Reproduction instructions

The `../check/` folder is a complete mirror of this directory containing the regenerated outputs. Compare each file in `../check/` with its counterpart here to verify that the published figures and tables are reproduced. The module-to-figure mapping below identifies which output belongs to which manuscript figure or table.

---

## 5. Code modules ↔ paper figures/tables

| Code module (under `02_Lineage_specific_Analysis_Code_Data/`) | Corresponding paper figure/table | Key scripts |
|---|---|---|
| `plant_phylo_micro_composition/` | Fig. 1a (plant ITS tree with class-level relative abundances of bacteria, fungi and protists); ITS outlier detection and removal of eight inconsistent samples used for plant-tree plotting | `Figure1A_revised.R`, `ITS_outlier_check.R` |
| `faithpd_treedistance_geodistances_across_plants/` | Fig. 1c; Suppl. Fig. 14c–e | `1.FaithPD_across_PlantOrders.R`, `3.TreePhyloDistance_across_PlantOrders-revised.R`, `3.GeoDistance_across_PlantOrders_revised.R`, `5.Root_Diameter_all_class_top10orders.R` |
| `mantel_test/` | Fig. 2a,b; Suppl. Fig. 13, 14a,b; Table S7a,b | `1.ASV_beta_mantel_test.R`, `2.KO_beta_mantel_test.R` |
| `tree_barplot/` | Fig. 2c–f (Asparagales, Arecales, Lamiales, Sapindales) and Suppl. Fig. 13a–f (Fabales, Rosales, Gentianales, Malpighiales, Malvales, Myrtales): plant phylogenetic trees (left) + absolute abundances of bacterial, fungal and protistan classes (right) | `1.tree_microbiome_barplot.R` |
| `random_forest/` | Fig. 3; indicator taxa in Table S8 | `amplicons_randomforest.R`, `kegg_randomforest.R`, `tree_clades_amplicon_randomforest.R`, `tree_clades_kegg_randomforest.R` |
| `coevolution_analysis_based_on_rf_and_lambda/tree_heatmap_dot/` | Fig. 3a–c (circular-tree dot plots) | `Bacterial/Fungal/Protistan_coevolved_ASV_tree_dotplot.R` |
| `coevolution_analysis_based_on_rf_and_lambda/plant_tree_lambda_rf/asv_qa/` | Fig. 3a–c; Suppl. Fig. 16; Table S8a–c | `Bacterial/Fungal/Protistan_coevolved_ASVs_stats.R`, `*_barplot.R`, `Plant_tree_sig_asv_and_barplot.R` |
| `phylogenetic_signal_lambda/` | Fig. 3a–f; Table S8a–f (Pagel's λ) | `amplicon/01-core_ASV_phylogenetic_signal.R`, `metagenomic/01-core_KO_phylogenetic_signal.R` |
| `coevolution_analysis_based_on_rf_and_lambda/bipartite_network/` | Fig. 3d–f (bipartite networks linking phylogenetically conserved, host-specialized bacterial (d), fungal (e) and protistan (f) KOs to plant orders) | `bipartite_network_revised.R` |
| `coevolution_analysis_based_on_rf_and_lambda/plant_tree_lambda_rf/ko_qa/` | Fig. 3d–h; ED Table 1; Suppl. Fig. 17, 18; Table S8d–f | `Enrichment_coevolved_KOs.R`, `Plant_coevolved_*_KO_enrichment_1KPM_ref.R`, `Plant_tree_*_sig_KO_barplot-revised.R` |
| `metastats/metastat-ASV/` | ED Fig. 4; Table S6b,c | `1.amplicon_metastats_analysis.R`, `2.metastats_pairwise_heatmap_plot.R`, `4.all_pairwise_enrichment_plot.R` |
| `metastats/metastat-KO/` | ED Fig. 5; Table S6d,e | `01.ko_metastat_analysis.R`, `02.metastat_pairwise_heatmap_plot.R`, `04.all_pairwise_enrichment_plot.R` |
| `ReporterScore/` | Suppl. Fig. 11; Table S6f | `get_the_reporter_score.R`, `ReporterScore_plot.R` |
| `distance_decay/` | Suppl. Fig. 12 (taxonomic distance–decay for four plant orders, panels a–d) | `1.distance_decay-ASV.R`, `1.distance_decay_for_all_orders-ASV.R` |
| `funguild/` | Suppl. Fig. 14f–h | `1.FUNGuildDatabase_anno.R`, `2.FUNGuild_plots.R`, `add_root_diameter_to_metadata.R` |
| `coevolution_analysis_based_on_rf_and_lambda/lambda_boxplot/` | Suppl. Fig. 15 | `lambda_boxplot.R` |
| `permanova/` | Table S6a (PERMANOVA statistics) | `1.PERMANOVA_stats.R` |
| `metadata/`, `feature_table/`, `asv_qa_alpha_beta/`, `ko_qa_alpha_beta/`, `kegg_abundance/` | Shared upstream inputs/outputs used by most modules; no direct paper figure/table generated in these modules | `01-get_absolute_alpha_diversity.R`, `02-get_absolute_beta_diversity.R`, `diversity/diversity.R`, `relative/pathway_tpm.R`, `absolute/quantative_profile.R` |

---

If you have any questions, please contact: **xlwang@scau.edu.cn**
