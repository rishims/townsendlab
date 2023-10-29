# Melanoma Stage-Specific Analysis
# Rishi Shah
# August 10, 2023

# Load necessary libraries
library(cancereffectsizeR)
library(BiocIO)
library(data.table)
library(ggplot2)
library(TCGAbiolinks)
library(car)
library(MASS)
library(ggpubr)
library(tidyverse)
library(dplyr)

# Import TCGA genomic and clinical data
tcga_skcm = 'TCGA-SKCM_with_met_exclude_multisample.maf.gz'
tcga_clinical <- GDCquery_clinic(project = "TCGA-SKCM", type = "clinical")
if (! file.exists(tcga_skcm)) {
  tmp_maf = paste0(tempfile(), '.maf')
  get_TCGA_project_MAF(project = 'TCGA-SKCM', filename = tcga_skcm, exclude_TCGA_nonprimary = FALSE)
  
  # 2 patients have 2 samples each. For simplicity, we'll exclude all.
  maf_to_edit = fread(tcga_skcm)
  stopifnot(maf_to_edit[multisample_patient == T, 
                        uniqueN(Tumor_Sample_Barcode) == 4 && uniqueN(Unique_Patient_Identifier) == 2])
  maf_to_edit = maf_to_edit[multisample_patient == FALSE]
  fwrite(maf_to_edit, tcga_skcm, sep = "\t")
  unlink(tcga_skcm)
}

tcga_maf <- preload_maf(maf = maf_to_edit, refset = "ces.refset.hg19", 
                        chain_file = "hg38ToHg19.over.chain", detect_hidden_mnv = FALSE)

tcga_maf <-  tcga_maf[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]

# Merge genomic and clinical data
tcga_clinical$Unique_Patient_Identifier <- tcga_clinical$bcr_patient_barcode
tcga_maf <- merge(tcga_maf, tcga_clinical, by = 'Unique_Patient_Identifier')

# Add staging and risk information to TCGA data
tcga_maf$Stage <- tcga_maf$ajcc_pathologic_stage

tcga_maf$Stage <- gsub("Stage ", "", tcga_maf$Stage)
tcga_maf$Stage <- gsub("A|B|C", "", tcga_maf$Stage)
tcga_maf$Stage <- gsub("Not Reported", NA, tcga_maf$Stage)

tcga_maf <- subset(tcga_maf, !is.na(Stage))

recode <- car::recode
tcga_maf$SAMPLE_TYPE <- recode(tcga_maf$Stage, "c(0, 'I', 'II', 'III') = 'Primary'; 'IV' = 'Metastasis'")

# Create CESA and load maf file
merged_cesa <- cancereffectsizeR::CESAnalysis(refset = "ces.refset.hg19") # using ref human genome 19

merged_cesa <- load_maf(cesa = merged_cesa, maf = tcga_maf, maf_name = "tcga_SKCM", sample_data_cols = c("SAMPLE_TYPE"))

# Create variable to hold table of primary/metastatic samples
tcga_breakdown <- table(merged_cesa$samples$Risk)

# Import GENIE genomic and clinical data
genie_clin <- read_tsv("data_clinical_sample.txt", comment="#")
genie_coverage <- read_tsv("genomic_information.txt")
genie_maf_all <- read_tsv("data_mutations_extended.txt")

genie_seq_IDs <- genie_clin %>% 
  filter(ONCOTREE_CODE == "SKCM")

# Find patients with more than one row in the clinical data,
# signifying they have more than one sample and maybe more than one assay
patients_w_more_than_one_row <- names(table(genie_seq_IDs$PATIENT_ID)[table(genie_seq_IDs$PATIENT_ID) > 1])
for(tumor_ind in 1:length(patients_w_more_than_one_row)){
  # find all the sequencing assays used in this patient
  these_panels <- genie_seq_IDs %>%
    filter(PATIENT_ID == patients_w_more_than_one_row[tumor_ind]) %>%
    pull(SEQ_ASSAY_ID)
  if(length(unique(these_panels)) > 1){
    # find panel with the most gene names
    panel_w_most_genes <- genie_coverage %>%
      filter(SEQ_ASSAY_ID %in% these_panels) %>%
      distinct(SEQ_ASSAY_ID,Hugo_Symbol) %>%
      count(SEQ_ASSAY_ID) %>%
      arrange(desc(n)) %>%
      dplyr::slice(1) %>%
      pull(SEQ_ASSAY_ID)
    # remove low coverage panels from clinical data from this patient, so we only
    # eventually load in tumor data from the largest coverage
    genie_seq_IDs <- genie_seq_IDs %>%
      filter(!(PATIENT_ID == patients_w_more_than_one_row[tumor_ind] &
                 SEQ_ASSAY_ID != panel_w_most_genes))
  }
}

patient_sample_data <- genie_seq_IDs[, c(1:2, 5)]

patient_sample_data <- patient_sample_data %>% 
  filter(!SAMPLE_TYPE %in% c("Unspecified", "Not Applicable or Heme", "Not Collected"))

genie_maf <- genie_maf_all %>%
  filter(Tumor_Sample_Barcode %in% patient_sample_data$SAMPLE_ID) %>%
  left_join(patient_sample_data, by = c("Tumor_Sample_Barcode" = "SAMPLE_ID"))

genie_maf_preload <- cancereffectsizeR::preload_maf(
  maf = genie_maf,
  refset = "ces.refset.hg19",
  sample_col = "PATIENT_ID",
  keep_extra_columns = c("SAMPLE_TYPE", "Tumor_Sample_Barcode"),
  detect_hidden_mnv = FALSE
)

genie_maf_preload <- genie_maf_preload %>%
  filter(germline_variant_site == FALSE) %>%
  filter(repetitive_region == FALSE | cosmic_site_tier %in% 1:3)

# Find unique assays we will query
these_seq_IDs <- unique(genie_seq_IDs$SEQ_ASSAY_ID)

# For every sequence assay ID we have corresponding to our tumors...
for(seq_ind in 1:length(these_seq_IDs)) {
  # ... Name this assay
  this_seq_assay <- these_seq_IDs[seq_ind]
  # Find all the samples using this assay
  samples_w_this_assay <- genie_seq_IDs %>%
    filter(SEQ_ASSAY_ID == this_seq_assay) %>%
    pull(SAMPLE_ID)
  # Find all the patients using this assay
  patients_w_this_assay <- genie_seq_IDs %>%
    filter(SEQ_ASSAY_ID == this_seq_assay) %>%
    pull(PATIENT_ID)
  # Make a dataset of just these data with this coverage
  this_maf <- genie_maf_preload %>%
    filter(Unique_Patient_Identifier %in% patients_w_this_assay)
  this_maf <- this_maf %>%
    group_by(Unique_Patient_Identifier) %>%
    filter(n_distinct(SAMPLE_TYPE) == 1) %>%
    ungroup()
  
  missing_tumors <- NULL
  
  # If there is data with this coverage
  if (nrow(this_maf) > 0) {
    # Pull out the coverage data
    this_coverage <- genie_coverage %>%
      filter(SEQ_ASSAY_ID == this_seq_assay) %>%
      filter(!str_detect(string = Chromosome, pattern = "Un_"))
    # Convert to GRanges
    this_coverage_granges <- GRanges(
      seqnames = this_coverage$Chromosome,
      ranges = IRanges(
        start = this_coverage$Start_Position,
        end = this_coverage$End_Position
      )
    )
    merged_cesa <- cancereffectsizeR::load_maf(
      cesa = merged_cesa,
      maf = this_maf,
      coverage = "targeted",
      covered_regions = this_coverage_granges,
      covered_regions_name = this_seq_assay,
      covered_regions_padding = 20,
      sample_data_cols = c("SAMPLE_TYPE")
    )
  } else {
    missing_tumors <- c(missing_tumors, this_seq_assay)
    message(paste(this_seq_assay, "is missing from the maf data"))
  }
  print(paste0(
    seq_ind,
    " out of ",
    length(these_seq_IDs),
    " GENIE MAFs loaded in..."
  ))
}

# Determine mutational signatures to exclude
signature_exclusions_skcm <- suggest_cosmic_signature_exclusions(cancer_type = "SKCM", treatment_naive = TRUE)

# Calculate tri-nucleotide mutation rates
merged_cesa <- trinuc_mutation_rates(merged_cesa,
                                     signature_set = ces.refset.hg19$signatures$COSMIC_v3.2,
                                     signature_exclusions = signature_exclusions_skcm)


# Calculate gene mutation rates
primary <- merged_cesa$samples[merged_cesa$samples$SAMPLE_TYPE == "Primary"]
metastatic <- merged_cesa$samples[merged_cesa$samples$SAMPLE_TYPE == "Metastasis"]

merged_cesa <- gene_mutation_rates(merged_cesa, covariates = ces.refset.hg19$covariates$SKCM, 
                                   samples = primary, save_all_dndscv_output = T)
merged_cesa <- gene_mutation_rates(merged_cesa, covariates = ces.refset.hg19$covariates$SKCM, 
                                   samples = metastatic, save_all_dndscv_output = T)


# Incorporate new sequential selection model
RefCDS = ces.refset.hg19$RefCDS
dndscv_gene_names <- merged_cesa$gene_rates$gene
nsyn_sites = sapply(RefCDS[dndscv_gene_names], function(x) colSums(x[["L"]])[1])

# Selecting mutation rate data for samples in the primary group
primary_samples <- length(unique(merged_cesa$dNdScv_results$rate_grp_1$annotmuts$sampleID))

# Selecting mutation rate data for samples in the metastatic group
metastatic_samples <- length(unique(merged_cesa$dNdScv_results$rate_grp_2$annotmuts$sampleID))

# Creating a data frame with mutation rate data for primary and metastatic groups
mut_rate_df <- tibble(gene = merged_cesa$dNdScv_results$rate_grp_1$genemuts$gene_name,
                      exp_primary_mu = merged_cesa$dNdScv_results$rate_grp_1$genemuts$exp_syn_cv,
                      exp_metastasis_mu = merged_cesa$dNdScv_results$rate_grp_2$genemuts$exp_syn_cv)

mut_rate_df$n_syn_sites = nsyn_sites[mut_rate_df$gene]

mut_rate_df %>% 
  mutate(primary_mu = (exp_primary_mu / n_syn_sites) / primary_samples) %>%
  mutate(metastasis_mu = (exp_metastasis_mu / n_syn_sites) / metastatic_samples) %>%
  mutate(cancer_greater = metastasis_mu > primary_mu) -> 
  mut_rate_df

# Defining rate 1 and rate 2 as mutation rates for primary and metastatic groups, respectively
rate_1 <- mut_rate_df|>
  select(gene, primary_mu)
rate_2 <- mut_rate_df|>
  select(gene, metastasis_mu)

# Change in mutation rate across stages
mut_rate_df <- mut_rate_df %>% 
  select(gene, primary_mu, metastasis_mu) %>% 
  mutate(p_1 = primary_mu / metastasis_mu) %>% 
  mutate(p_2 = 1 - p_1)

# Saving "last" gene mutation rates into separate data frame, "last" rates meaning from last stage metastasis_mu
set_cancer_rates <- mut_rate_df %>%
  select(gene, metastasis_mu) %>%
  data.table::setDT()

# Clear the gene rates in the cesa object 
merged_cesa <- clear_gene_rates(cesa = merged_cesa)

# Setting gene rates to highest rates from metastasis_mu
merged_cesa <- set_gene_rates(cesa = merged_cesa, rates = set_cancer_rates, missing_genes_take_nearest = T)

# Calculate selection intensities, use sequential model
source("new_sequential_lik.R")

for(ind in 1:length(merged_cesa$variants)){
  
  this_var <- merged_cesa$variants[ind, ]
  
  this_gene <- this_var$gene
  these_props <- mut_rate_df[mut_rate_df$gene == this_gene, c("p_1","p_2")]
  these_props <- c(these_props$p_1, these_props$p_2)
  
  merged_cesa <- ces_variant(cesa = merged_cesa, 
                             model = sequential_lik_dev, ordering_col = 'SAMPLE_TYPE', 
                             ordering = c('Primary', 'Metastasis'),
                             lik_args = list(sequential_mut_prop = these_props),
                             run_name = "seq_mod")
  
}

# Optional Plotting

# Plotting single-nucleotide somatic variant profile
snv_counts <- merged_cesa$mutational_signatures$snv_counts

summed_snv_by_group <- data.table()
receptor_groups <- unique(na.omit(merged_cesa$samples$SAMPLE_TYPE))
samples_with_snvs <- merged_cesa$samples[colnames(snv_counts), on = "Unique_Patient_Identifier"]
for (grp in receptor_groups) {
  curr_samples <- samples_with_snvs[grp, Unique_Patient_Identifier, on = "SAMPLE_TYPE"]
  curr_snv_sum <- rowSums(snv_counts[, curr_samples])
  summed_snv_by_group[, (grp) := curr_snv_sum]
}

summed_snv_by_group <- as.matrix(summed_snv_by_group)
rownames(summed_snv_by_group) <- rownames(snv_counts)
# pdf(file = "tcga_skcm_tri_nucleotide_mutations_plot.pdf", height = 10, width = 10)
MutationalPatterns::plot_96_profile(summed_snv_by_group, ymax = 0.25)
# dev.off()

# Plotting gene mutation rates
primary_plot <- ggplot(merged_cesa$gene_rates, aes(x = rate_grp_1, y = (..count../max(..count..)))) + xlim(0, 1.2e-04) +
  geom_density(color = "red", linewidth = 1.1) + labs(x = "Mutation rate in primary tumors", y = "Density") + 
  theme_pubr() + geom_vline(aes(xintercept = median(rate_grp_1)), color = "blue", linetype = "dashed", linewidth = 1)

mets_plot <- ggplot(merged_cesa$gene_rates, aes(x = rate_grp_2, y = (..count../max(..count..)))) + xlim(0, 1.2e-04) +
  geom_density(color = "red", linewidth = 1.1) + labs(x = "Mutation rate in metastatic tumors", y = "Density") + 
  theme_pubr() + geom_vline(aes(xintercept = median(rate_grp_1)), color = "blue", linetype = "dashed", linewidth = 1)

primary_vs_mets_risk_plot <- ggplot(merged_cesa$gene_rates, aes(x = rate_grp_1, y = rate_grp_2)) + ylim(0, 1.2e-04) +
  xlim(0, 1.2e-04) + geom_point(alpha = 0.5, color = "red") +
  labs(x = "Mutation rate in primary tumors", y = "Mutation rate in metastatic tumors") + 
  geom_abline(slope = 1, linetype = "dashed", color = "maroon") + theme_pubr()

# Plotting selection results
top <- merged_cesa$selection$skcm
top <- top[order(-si_Primary)][1:20] # take top 20 by SI
top <- top[order(si_Primary)] # will plot lowest to highest (left to right)
top[, display_name := gsub("_", "\n", variant_name)]
top[, display_levels := factor(display_name, levels = display_name, ordered = T)]

plot_title <- "Top Cancer Effects in Skin Cutaneous Melanoma"
n.dodge <- 2

breaks <- unique(as.numeric(round(quantile(top$si_Primary, probs = c(0, .5, .75, 1)))))

top$si_Primary <- as.numeric(as.character(top$si_Primary))

selection_plots <- ggplot(top, aes(x = display_levels, y = si_Primary)) +
  geom_errorbar(aes(ymin = ci_low_95_si_Primary, ymax = ci_high_95_si_Primary), width = .2, color = "darkgrey") +
  geom_point(aes(color = display_levels), size = 3) +
  scale_x_discrete(guide = guide_axis(n.dodge = n.dodge)) + scale_y_log10() + 
  scale_color_viridis_c(
    name = "variant prevalence", guide = "colorbar", trans = "log10",
    option = "plasma", breaks = 0:300
  ) +
  skcm_plot_m + xlab(element_blank()) +
  ylab(expression("cancer effect" ~ scriptstyle(~ ~ (log[10])))) +
  ggtitle(plot_title) +
  guides(color = guide_colourbar(ticks = FALSE)) +
  theme_minimal() +
  theme(
    text = element_text(family = "Verdana"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 8),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )