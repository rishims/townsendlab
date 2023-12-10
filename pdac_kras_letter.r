# Calculate cancer effect sizes in PDAC data with epistasis analysis
# October 29, 2023

# load libraries -----

library(cancereffectsizeR)
library(tidyverse)
library(ces.refset.hg19)
library(TCGAbiolinks)
library(ggplot2)
library(data.table)
library(dplyr)
library(patchwork)
library(cowplot)

# the code below can be run to replicate our analysis or alternatively, the complete CESA can be loaded in as is done directly below this comment
pdac_cesa_complete <- load_cesa('pdac_kras_letter_cesa.rds')

# define CESA object -----
pdac_cesa <- cancereffectsizeR::CESAnalysis(refset = ces.refset.hg19)

# import whole genome / exome ----- 
#' Files were manually downloaded from links provided

# generate TCGA data -----
# Generate TCGA maf file, match records and find matching IDs
tcga_paad <- "TCGA-PAAD.maf.gz"
tcga_clinical <- GDCquery_clinic(project = "TCGA-PAAD", type = "clinical")

if (!file.exists(tcga_paad)) {
  get_TCGA_project_MAF(project = "PAAD", filename = tcga_paad, 
                       exclude_TCGA_nonprimary = FALSE)
}

maf_tcga_PAAD <- preload_maf(maf = tcga_paad, refset = "ces.refset.hg19", chain_file = "hg38ToHg19.over.chain", detect_hidden_mnv = FALSE)
maf_tcga_PAAD <- maf_tcga_PAAD[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]


# + UTSW -----
# Witkiewicz, A. K., McMillan, E. A., Balaji, U., Baek, G., Lin, W. C., Mansour, J., et al.
  # Whole-exome sequencing of pancreatic cancer defines genetic diversity and therapeutic targets. Nature communications, 6, 6744.
# https://cbioportal-datahub.s3.amazonaws.com/paad_utsw_2015.tar.gz 
maf_utsw_PAAD <- preload_maf(maf = "paad_utsw_2015_data_mutations_extended\\paad_utsw_2015\\data_mutations_extended.txt",
                             refset = ces.refset.hg19)
maf_utsw_PAAD = maf_utsw_PAAD[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]

# + ICGC ---
# Biankin, A. V., Waddell, N., Kassahn, K. S., Gingras, M. C., Muthuswamy, L. B., Johns, A. L., et al.
  # Pancreatic cancer genomes reveal aberrations in axon guidance pathway genes. Nature, 491(7424), 399–405.
# https://cbioportal-datahub.s3.amazonaws.com/paad_icgc.tar.gz 
maf_icgc_PAAD <- preload_maf(maf = "paad_icgc_data_mutations_extended\\paad_icgc\\data_mutations.txt",
                             refset = ces.refset.hg19)
maf_icgc_PAAD = maf_icgc_PAAD[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]

# + QCMG ---- 
# Bailey, P., Chang, D. K., Nones, K., Johns, A. L., Patch, A. M., Gingras, M. C., et al.
  # Genomic analyses identify molecular subtypes of pancreatic cancer. Nature, 531(7592), 47–52.
# https://cbioportal-datahub.s3.amazonaws.com/paad_qcmg_uq_2016.tar.gz 
maf_qcmg_PAAD <- preload_maf(maf = "paad_qcmg_uq_data_mutations_extended\\paad_qcmg_uq_2016\\data_mutations_extended.txt",
                             refset = ces.refset.hg19)
maf_qcmg_PAAD = maf_qcmg_PAAD[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]

# + CPTAC 
# Cao, L., Huang, C., Cui Zhou, D., Hu, Y., Lih, T. M., Savage, S. R., et al.
  # Proteogenomic characterization of pancreatic ductal adenocarcinoma. Cell, 184(19), 5031–5052.e26.
# https://cbioportal-datahub.s3.amazonaws.com/paad_cptac_2021.tar.gz
maf_paad_cptac_2021 <- preload_maf("cptac\\paad_cptac_2021\\data_mutations.txt",refset = ces.refset.hg19)
maf_paad_cptac_2021 = maf_paad_cptac_2021[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]


# + Yale-Gilead collaboration
# https://pubmed.ncbi.nlm.nih.gov/30365005/
maf_yg <- preload_maf(maf = "mutationsTN_26_Pancreatic_Cancer.maf",chr_col = "Chrom",sample_col = "Patient_ID",refset = ces.refset.hg19)
maf_yg = maf_yg[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]


# import targeted sequencing -----

# + GENIE-----
  # The AACR Project GENIE Consortium. AACR Project GENIE: Powering Precision Medicine Through An International Consortium, Cancer Discov. 2017 Aug;7(8):818-831
# downloaded using synapse package within "/data_download/synapse_downloads.R"

genie_clin <- read_tsv("data_clinical_sample.txt", comment="#")
genie_coverage <- read_tsv("genomic_information.txt")
genie_maf <- read_tsv("data_mutations_extended.txt")

genie_clin_paad <- genie_clin %>% 
  filter(CANCER_TYPE == "Pancreatic Cancer" &
           (CANCER_TYPE_DETAILED == "Pancreatic Adenocarcinoma" |
              CANCER_TYPE_DETAILED == "Adenosquamous Carcinoma of the Pancreas" |
              CANCER_TYPE_DETAILED == "Acinar Cell Carcinoma of the Pancreas"))

genie_maf <- genie_maf %>% 
  filter(Tumor_Sample_Barcode %in% genie_clin_paad$SAMPLE_ID)

genie_maf <- genie_maf %>%
  filter(Chromosome != "X")

genie_maf %>% 
  mutate(Chromosome = case_when(
    Chromosome == "23" ~ "X" , 
    TRUE ~ Chromosome)) -> genie_maf

# adding in the patient
genie_maf <- left_join(x = genie_maf, y=genie_clin_paad[,c("PATIENT_ID","SAMPLE_ID")],
                       by = c("Tumor_Sample_Barcode" = "SAMPLE_ID"))

# find the samples with the most substitutions
genie_patients_to_keep <- genie_maf %>% 
  count(PATIENT_ID, Tumor_Sample_Barcode) %>% 
  group_by(PATIENT_ID) %>% 
  slice_max(order_by = n, n=1,with_ties = F) %>% 
  pull(Tumor_Sample_Barcode)

# filter for samples with the most substitutions 
genie_maf %>% 
  filter(Tumor_Sample_Barcode %in% genie_patients_to_keep) -> 
  genie_maf

genie_maf <- preload_maf(maf = genie_maf,refset = ces.refset.hg19,keep_extra_columns = c("PATIENT_ID"))

# filter for tumor subtypes ------ 

## QCMG

QCMG_clin <- read_tsv(file = "paad_qcmg_uq_data_mutations_extended\\paad_qcmg_uq_2016\\paad_qcmg_uq_2016_clinical_data.tsv")

QCMG_clin %>% 
  filter(`Tumor Other Histologic Subtype` %in% "Pancreatic Ductal Adenocarcinoma") %>% 
  mutate(same_patient_sample = `Patient ID` == `Sample ID`) %>%
  pull(`Patient ID`) -> 
  qcmg_patients_to_keep

maf_qcmg_PAAD %>% 
  filter(Unique_Patient_Identifier %in% qcmg_patients_to_keep) -> 
  maf_qcmg_PAAD

## TCGA

tcga_biospecimen <- read.csv("TCGA clinical\\gdac.broadinstitute.org_PAAD.Merge_Clinical.Level_1.2016012800.0.0\\PAAD.clin.merged.txt",sep = "\t",header = F)

tcga_biospecimen <- as.data.frame(t(tcga_biospecimen))
colnames(tcga_biospecimen) <- tcga_biospecimen[1,]

tcga_biospecimen <- tcga_biospecimen[-1,]


tcga_biospecimen %>% 
  filter(patient.histological_type == "pancreas-adenocarcinoma ductal type") %>% 
  pull(patient.bcr_patient_barcode) %>% 
  toupper() -> 
  tcga_samples_needed 

maf_tcga_PAAD %>% 
  mutate(bcr_patient_barcode = str_sub(Unique_Patient_Identifier,1,12)) %>% 
  filter(bcr_patient_barcode %in% tcga_samples_needed) -> 
  maf_tcga_PAAD

# Check for duplications -----

unique(maf_icgc_PAAD$Unique_Patient_Identifier)
unique(maf_qcmg_PAAD$Unique_Patient_Identifier)

maf_icgc_PAAD <- maf_icgc_PAAD %>% 
  mutate(icgc_identifier = str_replace(string = Unique_Patient_Identifier,pattern = "_TD",replacement = ""))

length(which(unique(maf_qcmg_PAAD$Unique_Patient_Identifier) %in% unique(maf_icgc_PAAD$icgc_identifier)))

# filter out identical names 
maf_icgc_PAAD <- maf_icgc_PAAD %>%
  filter(!icgc_identifier %in% maf_qcmg_PAAD$Unique_Patient_Identifier)





all_mafs <- list(maf1 = maf_icgc_PAAD,
                 maf2 =  maf_qcmg_PAAD,
                 maf3 = maf_utsw_PAAD, 
                 maf4 = maf_tcga_PAAD, 
                 maf5  = maf_paad_cptac_2021)


combined_maf <- data.table::rbindlist(all_mafs, idcol="source",fill=T)



combined_maf <- combined_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

possible_dups <- check_sample_overlap(combined_maf)

possible_dups %>%
  filter(variants_shared > 5)

to_remove <- c("ICGC_0229",
               "PCSI0022_T",
               "ICGC_0150",
               "TCGA-US-A774-01A-21D-A32N-08",
               "PCSI0007_T",
               "ICGC_0097",
               "PCSI0024_T",
               "ICGC_0164")

# load in substitution data ---- 


# + Panel data ------ 

## + GENIE-----

genie_coverage <- read_tsv(file = "genomic_information.txt",  col_types = list( Chromosome = col_character()))

genie_coverage <- genie_coverage %>% 
  filter(Chromosome!="Un_gl000228")

# loop to load in genie data with correct panel 
for(seq_assay_ind in 1:length(unique(genie_clin_paad$SEQ_ASSAY_ID))){
  
  # get sequencing assay data
  this_seq_assay <- unique(genie_clin_paad$SEQ_ASSAY_ID)[seq_assay_ind]
  
  # get tumor IDs that use this assay 
  these_tumor_names <- genie_clin_paad %>% 
    filter(SEQ_ASSAY_ID == this_seq_assay) %>% 
    pull(SAMPLE_ID) %>% 
    unique()
  
  # get mutation data of these tumors
  these_tumors_maf <- genie_maf %>% 
    filter(Unique_Patient_Identifier %in% these_tumor_names)
  
  # load in each set up variants with the correct panel 
  if(nrow(these_tumors_maf) > 0 & !("GENIE-YALE-TPL478-1" %in% these_tumors_maf$Unique_Patient_Identifier)){
    
    this_coverage <- genie_coverage %>% 
      filter(SEQ_ASSAY_ID == this_seq_assay)
    
    this_coverage_granges <- GRanges(seqnames = this_coverage$Chromosome, ranges = IRanges(start = this_coverage$Start_Position,end = this_coverage$End_Position))
    
    # load into cesa
    
    pdac_cesa <- load_maf(cesa = pdac_cesa, maf = these_tumors_maf,
                          coverage = "targeted",
                          covered_regions = this_coverage_granges,
                          covered_regions_name = this_seq_assay)
  }
}

# + Whole exome data ---- 

pdac_cesa <- load_maf(cesa = pdac_cesa, 
                      maf = maf_icgc_PAAD %>% 
                        filter(!Unique_Patient_Identifier %in% to_remove))




pdac_cesa <- load_maf(cesa = pdac_cesa, maf = maf_qcmg_PAAD %>% 
                        filter(! Unique_Patient_Identifier %in% to_remove))





pdac_cesa <- load_maf(cesa = pdac_cesa, maf = maf_utsw_PAAD %>% 
                        filter(! Unique_Patient_Identifier %in% to_remove))




pdac_cesa <- load_maf(cesa = pdac_cesa, maf = maf_paad_cptac_2021 %>% 
                        filter(! Unique_Patient_Identifier %in% to_remove))





pdac_cesa <- load_maf(cesa = pdac_cesa, maf = maf_tcga_PAAD %>% 
                        filter(!Unique_Patient_Identifier %in% to_remove))




# Mutation rate estimates---- 

## trinucleotide rates ----
pdac_cesa <- cancereffectsizeR::trinuc_mutation_rates(cesa = pdac_cesa,signature_set = "COSMIC_v3.2",
                                                      signatures_to_remove = cancereffectsizeR::suggest_cosmic_signatures_to_remove(cancer_type = "PAAD"),
                                                      cores = 2)

## gene by gene mutation rates  ---- 
pdac_cesa <- cancereffectsizeR::gene_mutation_rates(cesa = pdac_cesa, 
                                                    covariates = "pancreas")



cancereffectsizeR::save_cesa(cesa = pdac_cesa, file = "pdac_full_analysis_2023.rds")

# generate graph of top 20 selection intensities -----
pdac_cesa <- ces_variant(cesa = pdac_cesa, run_name = "kras")

# extract selection results from CESAnalysis and take top variants for visualization
top <- pdac_cesa$selection$kras
top <- top[order(-selection_intensity)][1:20] # take top 20 by SI
top <- top[order(selection_intensity)] # will plot lowest to highest (left to right)

top[, display_name := gsub("_", "\n", variant_name)]
top[, display_levels := factor(display_name, levels = display_name, ordered = T)]

plot_title <- ""
n.dodge <- 2 

breaks <- unique(as.numeric(round(quantile(top$included_with_variant, probs = c(0, .5, .75, 1)))))

# create and save plot
pdac_si_plot <- ggplot(top, aes(x = display_levels, y = log10(selection_intensity))) +
  geom_errorbar(aes(ymin = log10(ci_low_95), ymax = log10(ci_high_95)), width = .2, color = "darkgrey") +
  geom_point(aes(color = included_with_variant), size = 3) +
  scale_x_discrete(guide = guide_axis(n.dodge = n.dodge)) +
  scale_y_log10() +
  scale_color_viridis_c(
    name = "variant prevalence", guide = "colorbar", trans = "log10",
    option = "plasma", breaks = breaks
  ) +
  xlab("variant name") +
  ylab(expression("cancer effect" ~ scriptstyle(~ ~ (log[10])))) +
  ggtitle(plot_title) +
  guides(color = guide_colourbar(ticks = FALSE)) +
  theme_minimal() +
  theme(
    text = element_text(family = "Verdana"),
    axis.title.x = element_text(size = 11),
    axis.text.x = element_text(size = 8),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave("pdac_si.jpeg", pdac_plot, width = 8, height = 5)

# perform pairwise epistasis analysis

# create list of all other genes mutated in combined PDAC dataset besides KRAS
other_genes = pdac_cesa$gene_rates$gene[pdac_cesa$gene_rates$gene != 'KRAS']
gene_pairs = lapply(other_genes, c, 'KRAS')

# calculate pairwise epistatic selection intensitites for each gene paired with KRAS
pdac_cesa <- ces_gene_epistasis(pdac_cesa, genes = gene_pairs, variants = "recurrent", run_name = "kras_epi")

# subset results for know driver genes and IDH1, which was found to have particularly significant suppressive effects on the selection of KRAS mutations
gene_ep_results <- subset(pdac_cesa$epistasis$kras_epi, 
                          variant_A %in% c('BRAF', 'NRAS', 'IDH1', 'SMAD4', 'TP53', 'CTNNB1', 'CDKN2A.p16INK4a'))

gene_ep_results[2]$variant_A <- 'CDKN2A'

# define function to plot pairwise epistatic effects
plot_epistasis_results = function(ep_results, current_gene, upper_limit1, upper_limit2) {
    data_for_compare <- ep_results %>%
    mutate(variant1 = gsub("\\.1.*", "", variant_A)) %>%
    mutate(variant2 = gsub("\\.1.*", "", variant_B)) %>%
    filter(variant1 == current_gene | variant2 == current_gene) #filter for GOI results
  
    data_for_compare <- data_for_compare %>%
    mutate(gene_of_interest = current_gene) %>%
    mutate(other_gene = case_when(
      variant1 == current_gene ~ variant2,
      variant2 == current_gene ~ variant1)) %>%
    mutate(ces_GOI = case_when(
      variant1 == current_gene ~ ces_A0,
      variant2 == current_gene ~ ces_B0)) %>%
    mutate(ces_OTHER = case_when(
      variant2 == current_gene ~ ces_A0,
      variant1 == current_gene ~ ces_B0)) %>%
    mutate(ces_GOI_after_OTHER = case_when(
      variant1 == current_gene ~ ces_A_on_B,
      variant2 == current_gene ~ ces_B_on_A)) %>%
    mutate(ces_OTHER_after_GOI = case_when(
      variant2 == current_gene ~ ces_A_on_B,
      variant1 == current_gene ~ ces_B_on_A)) %>%
    mutate(joint_cov_samples_just_GOI = case_when(
      variant2 == current_gene ~ nB0,
      variant1 == current_gene ~ nA0)) %>%
    mutate(joint_cov_samples_just_OTHER = case_when(
      variant1 == current_gene ~ nB0,
      variant2 == current_gene ~ nA0)) %>%
    mutate(p_epistasis_scientific = formatC(p_epistasis, format = "e", digits = 2)) %>%
    mutate(significance = case_when(p_epistasis <= 0.001 ~ "***", 
                                    p_epistasis <= 0.01 ~ "**",
                                    p_epistasis <= 0.05 ~ "*", 
                                    p_epistasis > 0.05 ~ "ns")) %>%
    mutate(change_GOI = abs(ces_GOI - ces_GOI_after_OTHER)) %>%
    mutate(change_other = abs(ces_OTHER - ces_OTHER_after_GOI))
  
    data_for_compare$other_gene <- factor(data_for_compare$other_gene, levels = unique(data_for_compare$other_gene))
    data_for_compare$other_gene <- reorder(data_for_compare$other_gene, -data_for_compare$change_other)
  
    text_size <- 18
    geom_text_size <- text_size * (5/14)
  
    GOI_plot <- ggplot(data = data_for_compare) +
    geom_segment(aes(x=as.numeric(other_gene)-0.1,xend=as.numeric(other_gene)+0.1, y=ces_GOI, yend=ces_GOI_after_OTHER), 
                 arrow = arrow(length = unit(0.1, "inches"), type ="closed")) +
    geom_point(aes(x=as.numeric(other_gene)-0.1, y=ces_GOI), size=2) +
    scale_x_discrete(limits = levels(data_for_compare$other_gene)) +
    theme_classic() +
    geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), lwd = 0.5, color = "lightgrey") +
    geom_hline(yintercept = 0, lwd = 0.5, color = "lightgrey", linetype = "dotted") +
    labs(x = paste0("Gene (paired with ", current_gene, ")"), y = paste0("\nScaled selection \ncoefficients of \n", current_gene, "\nwhen paired \ngene is \nwildtype (\u25cf) \nand mutated (\u25b8)")) +
    theme(text = element_text(family = "Verdana"),
          axis.text = element_text(size = text_size), 
          axis.title = element_text(size = text_size), 
          axis.title.x = element_text(vjust = -5),
          axis.title.y = element_text(angle = 0, margin = margin(r = 30)),
          plot.margin = margin(0, 1, 2, 1, "cm")) +
    annotate("text", x = as.numeric(data_for_compare$other_gene), y = -(upper_limit1/2.8), label = data_for_compare$significance, size = geom_text_size) +
    coord_cartesian(ylim = c(0, upper_limit1), clip = "off")

    OTHER_plot <- ggplot(data = data_for_compare) +
    geom_segment(aes(x=as.numeric(other_gene)-0.1, xend=as.numeric(other_gene)+0.1, y=ces_OTHER, yend=ces_OTHER_after_GOI), 
                 arrow = arrow(length = unit(0.1, "inches"), type ="closed")) +
    geom_point(aes(x=as.numeric(other_gene)-0.1, y=ces_OTHER), size=2) +
    scale_x_discrete(limits = levels(data_for_compare$other_gene)) +
    theme_classic() +
    theme(text = element_text(family = "Verdana"),
          axis.text = element_text(size = text_size), 
          axis.title = element_text(size = text_size), 
          axis.title.x = element_text(vjust = -5),
          axis.title.y = element_text(angle=0, margin = margin(r = 30)),
          plot.margin = margin(1, 1, 2, 1, "cm")) +
    geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), lwd = 0.5, color = "lightgrey") +
    geom_hline(yintercept = 0, lwd = 0.5, color = "lightgrey", linetype = "dotted") +
    labs(x = paste0("Gene (paired with ",current_gene, ")\n"), y = paste0("\nScaled selection \ncoefficients of \npaired gene when\n", current_gene, "\nis wildtype (\u25cf) \nand mutated (\u25b8)"))+
    scale_color_discrete(name="p-value", breaks=c("yes", "no"), labels=c("p < 0.05", "p ≥ 0.05")) +
    annotate("text", x = as.numeric(data_for_compare$other_gene), y = -(upper_limit2/2.8), label = data_for_compare$significance, size = geom_text_size) +
    coord_cartesian(ylim = c(0, upper_limit2), clip = "off")
  
    both_plots <- (OTHER_plot/GOI_plot)
  
    return(both_plots)
}

# create and save plot
kras_epi_plot <- plot_epistasis_results_with_ci(gene_ep_results, "KRAS", 525000, 30000)
ggsave("kras_epi_fig.jpeg", kras_ep_plot, width = 18, height = 14)
