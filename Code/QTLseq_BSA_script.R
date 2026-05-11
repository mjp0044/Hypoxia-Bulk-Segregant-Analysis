##QTLseq analysis script for BSA project
#Matthew Powers 
#powerm3@oregonstate.edu
#Combines count files 
#Last edited 3-16-2026
#Built originally on R version 4.4.2 
#Last tested on 4.5.1

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.21") # Use 3.21 with R 4.5.1 or newer

# #library(remotes) #To install packages from github using 'install_github' function as below
# devtools::install_github("bmansfeld/QTLseqr")
# 
# #Use to install biomanager packages
# BiocManager::install("QTLseqr")

#Load required and preferred packages
library(MESS) #data wrangling and modeling
library(stats) #data wrangling and modeling
library(ggpubr) #Also loads basic ggplot2
library(ggridges) # For ridgeline plots
library(gggenes) # to make gene tracks
library(patchwork) # Stitch plots
library(cowplot) #Pretty ggplots!
library(dplyr) #Data wrangling
library(tidyverse) #Data wrangling
library(stringr) #data wrangling for splitting column characters
library(MetBrewer) #Pretty colors!
library(MexBrewer) # Colores bonito!
library(scales) # Viewing pretty colors!
library(QTLseqr) #QTL analysis
library(data.table) # For manipulating large sets
#library(fuzzyjoin) # Use to map SNP positions in results to gff information with gene ID
library(GenomicRanges) # Use to map SNP positions in results to gff information with gene ID
library(rtracklayer) #For querying with genomicranges package
library(writexl) # Write results to xlsx format
library(reactable) #Generate html files for interactable tables
library(reactablefmtr) #Expansion to reactable functionality
library(vcfR) # Read in and manipulate vcf files
library(seqinr) # For functions to help calculate dn/ds
library(VariantAnnotation) #To parse vcf from snpeff
library(flextable) # #Export data frames as tables
library(officer) #For use with flextable objects being inserted into word
library(openxlsx) # Exporting data frames as xlsx files

#Go term testing
library(topGO) #GO analysis
library(rrvgo) #Vizualize results from GO enrichment analysis
library(Rgraphviz) #Visualize results from topGO

#install.packages("./org.Tcalifornicus.eg.db", type = "source", repos=NULL)
library(org.Tcalifornicus.eg.db) #Load custom OrgDb package for T. californicus
Tcalif_orgdb_object <- org.Tcalifornicus.eg.db #Assign package to an object for use in rrvgo

#Set theme globally
theme_set(theme_cowplot())


#Read in count tables
# Will remove these later after we filter them down to interesting snps from the analysis
  low <- read.table("Bottom_S1_SNP_FINAL_COUNTS.txt")
  head(low)
  
  high <- read.table("Top_S2_SNP_FINAL_COUNTS.txt")
  head(high)
  

#Merge the two count tables by chromosome and SNP (just keeps rows shared between two data frames)
  comb <- merge(low, high, by=c("V1","V2"))
 #Order the data frame in ascending order 
  comb <- comb[order(comb$V1, comb$V2),]
  head(comb)

#Sanity check to make sure reference allele and alternative allele match
  identical(comb$V3.x, comb$V3.y)
  identical(comb$V4.x, comb$V4.y)
  comb[which(comb$V4.x != comb$V4.y), ] #Print the line(s) with the mismatch
# One row had a different alternative allele in V4.y so removed that manually.
  comb <- comb[comb$V4.x == comb$V4.y, ]

#remove columns 8 and 9, which are repetitive of 3 and 4.
  comb=comb[c(1:6,10:11)]
  
  head(comb)
  
  
#Rename the columns to match the required nomenclature for QTLseqr
  colnames(comb) <- c("CHROM", "POS", "REF", "ALT", 
                      "AD_REF.BOTTOM", "AD_ALT.BOTTOM", 
                      "AD_REF.TOP", "AD_ALT.TOP")
  
  head(comb)
  
#Rename chromosome names before exporting
  comb %>%
    mutate(CHROM = fct_recode(CHROM,
                          "Chr_1" = "SD_Chromosome_1",
                          "Chr_2" = "SD_Chromosome_2",
                          "Chr_3" = "SD_Chromosome_3",
                          "Chr_4" = "SD_Chromosome_4",
                          "Chr_5" = "SD_Chromosome_5",
                          "Chr_6" = "SD_Chromosome_6",
                          "Chr_7" = "SD_Chromosome_7",
                          "Chr_8" = "SD_Chromosome_8",
                          "Chr_9" = "SD_Chromosome_9",
                          "Chr_10" = "SD_Chromosome_10",
                          "Chr_11" = "SD_Chromosome_11",
                          "Chr_12" = "SD_Chromosome_12")) -> comb
  
  
  write.csv(comb, "Combined_counts_forQTLseqr.txt", row.names=FALSE, quote=FALSE)
  
  
# remove comb data frame to save R space
  rm(comb, high, low)
  
  #Create vector of chromosomes to be used. The "Chr_" needs to match exactly what is in the data files.
  Chroms <- paste0(rep("Chr_", 1), 1:12)
  Chroms
  
  
#Import the data file, pre-formatted according to the manual.
  df <- importFromTable(file = "Combined_counts_forQTLseqr.txt",
                        highBulk = "TOP",
                        lowBulk = "BOTTOM",
                        chromList = Chroms)
  head(df)
  
  
#Before filtering, check total read depth, total reference allele frequency, and per bulk SNP index in our data set
  # Total read depth
    ggplot(data = df) +
      geom_histogram(aes(x = DP.HIGH + DP.LOW)) +
      scale_y_continuous(labels = label_comma())+
      xlim(0,1000)
  
  # Total reference allele frequency
    ggplot(data = df) +
      geom_histogram(aes(x = REF_FRQ))+
      scale_y_continuous(labels = label_comma())
  
 # Per bulk SNP-index
  #High
    ggplot(data = df) +
      geom_histogram(aes(x = SNPindex.HIGH))+
      scale_y_continuous(labels = label_comma())
  #Low
    ggplot(data = df) +
      geom_histogram(aes(x = SNPindex.LOW))+
      scale_y_continuous(labels = label_comma())
    
#Filtering
  df_filt <-
    filterSNPs(
      SNPset = df,
      refAlleleFreq = 0.1,
      minTotalDepth = 80,
      maxTotalDepth = 400,
      depthDifference = 100,
      minSampleDepth = 40,
      verbose = TRUE
    )

  
# Remove df dataframe to save R space
  rm(df)
  
#G' analysis. This is the main analysis.

## NOTE about windowSize: Depending on genome size and data complexity, only a certain
## range of window sizes can be used. Going outside of that range will produce 
## memory errors and the command does not complete. Simply try a few values to find out the
## working window size. For our copepod data, I've found that the lowest allowable size is 200,000, 
## and the best range balancing resolution and accuracy is 500,000 to 1,000,000.

### Also about windowSize: large differences in window size will affect statistical thresholds
### for significance of peaks. That is partly because smaller windows will produce more point estimates,
### and hence suffer more impact of multiple testing. 
  df_filt2 <- runGprimeAnalysis(df_filt,
                                windowSize = 750000,
                                outlierFilter = "deltaSNP",
                                filterThreshold = 0.1)
  
  # Export bed graph of our Gprime results
    # Remap Chr_1 -> Chromosome_1
    df_filt2$CHROM_igv <- gsub("Chr_", "Chromosome_", df_filt2$CHROM)
    
    # Export bedGraph
    gprime_bedgraph <- data.frame(
      chrom      = df_filt2$CHROM_igv,
      chromStart = df_filt2$POS - 1,  # 0-based
      chromEnd   = df_filt2$POS,
      score      = df_filt2$Gprime
    )
    
    # Convert to integer
    gprime_bedgraph$chromStart <- as.integer(gprime_bedgraph$chromStart)
    gprime_bedgraph$chromEnd   <- as.integer(gprime_bedgraph$chromEnd)
    
    # round Gprime
    gprime_bedgraph$score <- round(df_filt2$Gprime, 4)
    
    # Export
    write.table(gprime_bedgraph, "gprime.bedgraph",
                sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
  
  
#PLOTS:

## Number of SNPs: Plot of number of SNPs along the genome (optional; not terribly useful, but
## can be used to check whether some regions are very poorly sequenced)
  nsnps <- plotQTLStats(SNPset = df_filt2, var = "nSNPs", scaleChroms = TRUE, col="darkblue") + 
    theme(strip.text = element_text(size = 9.5))
  nsnps
  #Export graph
  jpeg(filename = "nSNPs.jpg", units= "in", width = 12, height = 4, res = 300)
  nsnps
  dev.off()
  
  
## G'. Main plot. Change 'q' as needed for the statistical threshold line.
  gprime.plot <- plotQTLStats(SNPset = df_filt2, var = "Gprime", plotThreshold = TRUE, q = 0.01, col="darkblue", scaleChroms = TRUE)+
    theme(strip.text = element_text(size = 9.5))
  gprime.plot
  
  #Export graph
  jpeg(filename = "Gprime.jpg", units= "in", width = 12, height = 4, res = 300)
  gprime.plot
  dev.off()
  
## The plots can also be subsetted to show only one or few chromosomes.
 chroms_with_peaks <- plotQTLStats(SNPset = df_filt2, var = "Gprime", plotThreshold = TRUE, q = 0.01,
               subset = c("Chr_3", "Chr_4", "Chr_8",
                          "Chr_9", "Chr_10", "Chr_11",
                          "Chr_12"), col = "darkblue")+
    theme(strip.text = element_text(size = 9.5))
 chroms_with_peaks
  #Export graph
  jpeg(filename = "Gprime_chroms_with_peaks.jpg", units= "in", width = 8, height = 3, res = 300)
  chroms_with_peaks
  dev.off()
  
  
  
  
# QTLseq analysis and plot: 
#This is an alternative analysis, but really only useful to be able to plot
# the nice deltaSNP plot with confidence intervals.

# ## NOTE 1: Make sure to use the original object (df_filt) and not the one from the G' analysis (df_filt2)
# ## NOTE 2: Pay attention to popStruc (which hybrid generation?), and bulkSize (number of individuals pooled for sequencing?)
#   df_filt3 <- runQTLseqAnalysis(SNPset = df_filt, windowSize = 750000, popStruc = "F30", 
#                                 bulkSize = 80, replications = 10000, intervals = c(95, 99) )
#   
#   # Generate plot with deltaSNP values
#   deltaSNP <- plotQTLStats(SNPset = df_filt3, var = "deltaSNP", plotIntervals = TRUE)+
#     theme(strip.text = element_text(size = 9.5))
#  deltaSNP
#  #Export graph
#  jpeg(filename = "deltaSNP.jpg", units= "in", width = 12, height = 4, res = 300)
#  deltaSNP
#  dev.off()
  
  
  
#Export significant regions

## If desired, the main G' prime output (df_filt2) can be exported as a table.
## This contains the stats for every single SNP, so it can be used for other types of
## analyses and for custom plots.

# For getting coordinates of the main peaks, use the function below and export as table.
 # Auto exports result as csv
  results <- getQTLTable(SNPset = df_filt2, method = "Gprime",alpha = 0.01, export = TRUE)
  results
  
  
# For getting the coordinates of all the SNPs
  all.snps <- getSigRegions(SNPset = df_filt2, method = "Gprime", alpha = 0.01)
  all.snps.comb <- bind_rows(all.snps) #Combine all the list elements into one data frame
  rm(all.snps) # Get rid of the intermediate list 
  
#For looking at G' prime distribution during QC
  plotGprimeDist(SNPset =df_filt2, outlierFilter = "deltaSNP", filterThreshold = 0.1)
  

########## Exporting vcf files for SNPeff  ##########
  
  # All SNPs  
  ### Format and export combined info to a vcf for SNPeff in bash
  vcf_df <- all.snps.comb %>%
    transmute(
      CHROM  = CHROM,
      POS    = POS,
      ID     = paste(CHROM, POS, REF, ALT, sep = "_"),
      REF    = toupper(REF),
      ALT    = toupper(ALT),
      QUAL   = ".",
      FILTER = "PASS",
      INFO   = "."
    )
  
  # Write the vcf file
  vcf_file <- "all_snps.vcf"
  # Create the header for the file 
  writeLines(c("##fileformat=VCFv4.2",
               "##source=all.snps.comb",
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), con = vcf_file)
  
  
  write.table(vcf_df, file = vcf_file, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  
  
  ## Subset top SNPs per QTL based on Gprime, regardless if exonic
  # For each QTL, select the SNP(s) with the maximum Gprime
  top_snps <- all.snps.comb %>%
    group_by(qtl) %>%
    filter(Gprime == max(Gprime, na.rm = TRUE)) %>%
    ungroup()
  
  # Prepare VCF from top exonic SNPs
  vcf_top_snps <- top_snps %>%
    transmute(CHROM = CHROM,
              POS = POS,
              ID = ".",
              REF = REF,
              ALT = ALT,
              QUAL = ".",
              FILTER = ".",
              INFO = ".")
  
  # Write the VCF file
  vcf_top_file <- "top_snps.vcf"
  
  # VCF header
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##source=all.snps.comb",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  ), con = vcf_top_file)
  
  # Append SNP data
  write.table(vcf_top_snps,file = vcf_top_file,append = TRUE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
  
  
  
  
  
  
   

########## Make summary of results from SNPeff output ###############
   
  # Read the SnpEff-annotated VCF
  # This object stores:
  #  Genomic coordinates (rowRanges)
  # Reference/alternate alleles (ref(), alt())
  # INFO fields (including ANN) in a structured way
  vcf <- readVcf("all_snps_annotated/all_snps.ann.vcf")
  
  # Inspect the raw ANN field for one variant
  info(vcf)$ANN[1]
  
  vcf_core <- as.data.frame(rowRanges(vcf))[, c("seqnames", "start")] #Extract core variant coordinates
  vcf_core$REF <- as.character(ref(vcf)) #Extracts the reference allele for each variant
  vcf_core$ALT <- as.character(unlist(alt(vcf))) #Extracts alternate alleles
  
  colnames(vcf_core)[1:2] <- c("CHROM", "POS") #Renames columns to conventional VCF-style naming
  
  #Extract the ANN field as a list
  ann_list <- info(vcf)$ANN
  
 # Convert ANN into long format (one row per annotation)
  ann_long <- data.frame(
    row_id = rep(seq_along(ann_list), lengths(ann_list)), # get how many entries and repeat SNP ID that many times
    ANN = unlist(ann_list, use.names = FALSE), # flatten annotation into a single vector
    stringsAsFactors = FALSE
  )
  
  # Reattach variant-level data to each annotation
  SNPeff_ann_long <- cbind(
    vcf_core[ann_long$row_id, ],
    ANN = ann_long$ANN
  )
  
  #Define expected ANN fields (SnpEff spec)
  ann_fields <- c(
    "Allele",
    "Effect",
    "Impact",
    "Gene_Name",
    "Gene_ID",
    "Feature_Type",
    "Feature_ID",
    "Transcript_BioType",
    "Rank",
    "HGVS_c",
    "HGVS_p",
    "cDNA_pos",
    "CDS_pos",
    "Protein_pos",
    "Distance",
    "Errors"
  )
  
  # Split ANN strings into columns
  # Splits each annotation string on |
  ann_matrix <- do.call(
    rbind,
    strsplit(SNPeff_ann_long$ANN, "\\|", fixed = FALSE)
  )
  
  # Normalize to exactly 16 columns
  # Some ANN entries are shorter (missing trailing fields)
  # This line:
  # Truncates extra fields if present
  # Pads missing fields with NA
  # Ensures consistent column alignment
  ann_matrix <- ann_matrix[, seq_len(16), drop = FALSE]
  colnames(ann_matrix) <- ann_fields
  
  # Build final tidy annotation table
  SNPeff_ann_long <- cbind(
    SNPeff_ann_long[, c("CHROM", "POS", "REF", "ALT")], #Applies meaningful column names
    as.data.frame(ann_matrix, stringsAsFactors = FALSE)
  )
  
  # Validation checks
  # multiple genes per SNP? I.e, we don't leave out any cases of multiple genes affected by a single snp
  any(duplicated(SNPeff_ann_long[, c("CHROM", "POS")]))
  
  # gene IDs correct format?
  head(unique(SNPeff_ann_long$Gene_ID))
  
  # Clean up intermediate files
  rm(vcf_core, ann_long, ann_list, ann_fields)
  
  # Filter SNPs before making gene level summary
  SNPeff_filtered <- SNPeff_ann_long %>%
    # Keep coding, intron, UTR, and upstream within 1000 bp
    filter(
      !(Effect == "downstream_gene_variant") &
        !(Effect == "upstream_gene_variant" & as.numeric(Distance) > 1000) &
        !(Effect == "intergenic_region")  # optionally also remove intergenic
    )
  
  
  library(dplyr)
  
  # Join SNPeff annotations to the SNP statistics
  # Use CHROM + POS as the key
  SNPeff_stats <- SNPeff_ann_long %>%
    left_join(
      all.snps.comb %>%
        dplyr::select(CHROM, POS, qtl, Gprime, qvalue),
      by = c("CHROM", "POS")
    )
  
  # Filter out SNPs we don't care about:
  # So the ones downstream or upstream > 1000 bp
  SNPeff_stats_filtered <- SNPeff_stats %>%
    filter(
      !(Effect == "downstream_gene_variant") &
        !(Effect == "upstream_gene_variant" & as.numeric(Distance) > 1000) &
        !(Effect == "intergenic_region")
    )
  
  # Build gene-level summary
  # Create list of effects labels corresponding to exonic snps
  exonic_effects <- c(
    "missense_variant",
    "synonymous_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "stop_retained_variant",
    "start_retained_variant",
    "missense_variant&splice_region_variant",
    "stop_gained&splice_region_variant",
    "stop_lost&splice_region_variant",
    "splice_region_variant&synonymous_variant",
    "splice_region_variant&stop_retained_variant",
    "3_prime_UTR_variant",
    "5_prime_UTR_variant",
    "5_prime_UTR_premature_start_codon_gain_variant"
  )
  # Make data frame summary
  gene_summary <- SNPeff_stats_filtered %>%
    group_by(Gene_ID) %>%
    summarise(
      n_SNPs = n(),
      n_in_exons = sum(Effect %in% exonic_effects, na.rm = TRUE),
      n_introns = sum(Effect %in% c("intron_variant", "splice_region_variant", "splice_acceptor_variant&intron_variant",
                                    "splice_region_variant&intron_variant", "splice_donor_variant&intron_variant"), na.rm = TRUE),
      n_upstream = sum(Effect == "upstream_gene_variant", na.rm = TRUE),
      n_modifer = sum(Impact == "MODIFIER", na.rm = TRUE),
      n_low_impact = sum(Impact == "LOW", na.rm = TRUE),
      n_mod_impact = sum(Impact == "MODERATE", na.rm = TRUE),
      n_high_impact = sum(Impact == "HIGH", na.rm = TRUE),
      effects = paste(unique(Effect), collapse = ";"),
      max_Gprime = max(Gprime, na.rm = TRUE),
      min_qvalue = min(qvalue, na.rm = TRUE),
      max_upstream_dist = ifelse(any(Effect == "upstream_gene_variant"),
                                 max(as.numeric(Distance[Effect == "upstream_gene_variant"]), na.rm = TRUE),
                                 NA_real_),
      min_upstream_dist = ifelse(any(Effect == "upstream_gene_variant"),
                                 min(as.numeric(Distance[Effect == "upstream_gene_variant"]), na.rm = TRUE),
                                 NA_real_),
      CHROM = unique(CHROM),
      QTL = unique(qtl)
    ) %>%
    ungroup()
  # After the fact, realized we should add n_in_gene column for modeling questions
  gene_summary$n_in_gene <- gene_summary$n_in_exons + gene_summary$n_introns
  
  
 # Sanity check to make n_snps matches n_in_exons + n_introns + n_upstream
  gene_summary_check <- gene_summary %>%
    mutate(
      classified_sum = n_in_exons + n_introns + n_upstream,
      diff = n_SNPs - classified_sum
    )
  
  table(gene_summary_check$diff) # Should indicate 0 for all 721 genes
  rm(gene_summary_check) # I don't want to save this 
  
# Read in gene effect summary from SNPeff in case we want to compare
  gene_summary_from_SNPeff <- read.delim("all_snps_annotated/all_snps.snpeff.genes.txt", 
                                         comment.char = "#", 
                                         stringsAsFactors = FALSE,
                                         header = TRUE)
  
  
# Check overlap between grabbing all genes in qtl with our snpeff summaries
  # Convert QTLseqr results to GRanges
    qtl_gr <- GRanges(
      seqnames = results$CHROM,
      ranges   = IRanges(start = results$start,
                         end   = results$end),
      QTL_ID   = results$qtl
    )
    
  # Import the GFF 
    gff_gr <- import("SDv2.2_Genes_NCBI.gff")
    
    # If they don't change one to fix it. 
    # I chose the gff to change
    seqlevels(gff_gr) <- sub("^Chromosome_", "Chr_", seqlevels(gff_gr))
    
   # Keep only gene features
    genes_gr <- gff_gr[gff_gr$type == "gene"]
    colnames(mcols(genes_gr))
    
    # Check to see if chromosome names match
    unique(seqnames(qtl_gr))
    unique(seqnames(genes_gr))
    
    # Find overlaps between QTLs and genes
    hits <- findOverlaps(qtl_gr, genes_gr)
    
    qtl_gene_table <- data.frame(
      QTL    = mcols(qtl_gr)$QTL_ID[queryHits(hits)],
      GeneID = mcols(genes_gr)$ID[subjectHits(hits)],
      stringsAsFactors = FALSE
    )
  
  # Check to see which genes don't match between genes found in snpeff and those found by genomic ranges in the QTLs
    setdiff(qtl_gene_table$GeneID, gene_summary_from_SNPeff$GeneId)

    
    
# While we are working with the gff, grab the length of the significant genes from our genomic ranges summary table
    # Add gene length to GRanges metadata as a new column 
    mcols(genes_gr)$gene_length <- GenomicRanges::width(genes_gr)
    # Pull the lengths out as a data frame
    genes_lengths <- as.data.frame(genes_gr) %>%
      dplyr::select(ID, gene_length)
    # Merge with gene_summary frame
    gene_summary <- gene_summary %>%
      left_join(genes_lengths, by = c("Gene_ID" = "ID"))
    # get rid of intermediate
    rm(genes_lengths)
    # sanity check, should be 0 NAs
    sum(is.na(gene_summary$gene_length))
    
    # Duplicate column in units of kb
    gene_summary$gene_length_kb <- gene_summary$gene_length / 1000
    # Make column with number of snps per kb of gene length
    gene_summary$n_in_gene_per_kb <- gene_summary$n_in_gene / gene_summary$gene_length_kb
    
    # Calculate ratio of snps upstream to snps in gene
    gene_summary$up_to_in_ratio <- gene_summary$n_upstream / gene_summary$n_in_gene
    gene_summary$up_to_in_ratio[is.infinite(gene_summary$up_to_in_ratio)] <- NA # Get rid of Inf cases where we divide by zero
    gene_summary$up_to_in_ratio[is.nan(gene_summary$up_to_in_ratio)] <- NA # Get rid of Inf cases where both are zero
    
    # Calculate ratio of snps upstream to snps in exons per kb
    gene_summary$up_to_in_per_kb_ratio <- gene_summary$n_upstream / (gene_summary$n_in_exons/ gene_summary$gene_length_kb)
    gene_summary$up_to_in_per_kb_ratio[is.infinite(gene_summary$up_to_in_per_kb_ratio)] <- NA # Get rid of Inf cases where we divide by zero
    gene_summary$up_to_in_per_kb_ratio[is.nan(gene_summary$up_to_in_per_kb_ratio)] <- NA # Get rid of Inf cases where both are zero

## Combine with gene annotations
# Read in latest blast annotation
  annot <- read.csv("SDv2.2_BlastTable.csv", header = TRUE)
  
  # Adjust sequence name to match my naming scheme (remove the -PA and the IF from TCALIF)
  annot$SeqName <- str_split_i(annot$SeqName, "-", 1)
  
  # Rename seqname to have the same column header as combined_sigs frame
  names(annot)[1] <- "Gene_ID"   
  
  # Merge the annotation with our combined_sigs data frame, keeping only the annotations that match to our entries
  gene_summary <- merge(gene_summary, annot, by = "Gene_ID", all.x = TRUE)
  
  
  # Read in significant gene list from RNA seq study
  time_series_sigs <- read.csv("Time Series Master Significant Genes list with Annotation.csv", header = TRUE)
  names(time_series_sigs)[1] <- "Gene_ID" # Rename first column
  
  # If you want to filter to be just  genes for some reason
  #time_series_sigs <- time_series_sigs[!is.na(time_series_sigs$clusters_),]
  
  # Add a column to our significant gene list indicating whether the gene was identified in our time series experiment
  gene_summary$In_time_series <- gene_summary$Gene_ID %in% time_series_sigs$Gene_ID
  
  # Reorder the columns for easier viewing
  gene_summary <- gene_summary[, c(15,16,1,23,24,2,17,19,20,21,22,3,4,5,6,7,8,9,10,11,12,13,14)]
  
  
  # Sanity check to make sure all genes in our summary from the SNP annotation are actually in the gene summary from SNPeff
  genes_missing_from_snpeff <- setdiff(unique(gene_summary$Gene_ID),unique(gene_summary_from_SNPeff$GeneId))
  # If 0, means all are present
  length(genes_missing_from_snpeff)
  
  
  # Check overlap with pcrit only genes from time series
    pcrit_only_genes <- read.csv("Pcrit_only_from_time_series.csv", header = T)
  # Add to frame
    gene_summary$In_pcrit_only <- gene_summary$Gene_ID %in% pcrit_only_genes$Gene
  
  
# Export significant genes and a list of the same number of unsignificant genes at text files
  # Get Id's in a list
  candidate_genes <- unique(gene_summary$Gene_ID)
  length(candidate_genes)  # should be 721
  # Add the -PA notation to match genome files
  candidate_genes <- paste0(candidate_genes, "-PA")
  
  
  write.table(candidate_genes,file = "candidate_genes_721.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  
  # Randomly sample 721 genes not associated with snps
  # Get list of genes not associated with our snps
  background_gene_pool <- setdiff(
    unique(annot$Gene_ID),
    unique(gene_summary_from_SNPeff$GeneId)
  )
  
  length(background_gene_pool) # Should be 14826
  
  # Randomly sample genes
  set.seed(1974)  # reproducibility (Go Sounders!)
  random_background_genes <- sample(
    background_gene_pool,
    size = 721,
    replace = FALSE
  )
  
  length(random_background_genes)  # should be 721
  
  # Add the -PA notation to match genome files
  random_background_genes <- paste0(random_background_genes, "-PA")
  
  write.table(
    random_background_genes,
    file = "random_background_genes_721.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  # Check for no overlap with *any* SNP-associated genes
  # Should be character(0)
  intersect(random_background_genes, unique(gene_summary_from_SNPeff$GeneId))
       
########## Functional analysis of genes in QTLs #############
       
       
    ## Go analysis 
      #Read in gene to GO mappings
       geneID2GO <- readMappings(file = "Genes2GOterms_SDv2.2_noPA.txt")
       
       geneUniverse <- names(geneID2GO) #assign gene names to a list. This will be our global list of genes. 
      
       
      #Get first column from results w/ gene names and pvalues. Send to list. 
       GeneInt.all <- as.character(gene_summary$Gene_ID) 
       
      #Make named vector with the gene universe where genes of interest are coded with a 1 so the 'new' function knows to focus on them. 
       geneList.all <- factor(as.integer(geneUniverse %in% GeneInt.all)) 
       names(geneList.all) <- geneUniverse
       
     # We now have all data necessary to build an object of type topGOdata. This object will contain all gene
     # identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information
     # needed to perform the desired enrichment analysis.
       all_GO <- new("topGOdata", description = "GO analysis of all sig genes", ontology = "BP", nodeSize = 10,
                          allGenes = geneList.all, annot = annFUN.gene2GO, gene2GO = geneID2GO)
       
     # Once we have an object of class topGOdata we can start with the enrichment analysis.
     # Run the classic fisher exact test to find enriched go terms
       resultFisher.all <- runTest(all_GO, algorithm = "weight01", statistic = "fisher")
       
       resultFisher.all #View results summary
       
     # GenTable is an easy to use function for analysing the most significant GO terms and the corresponding p
     # values. We list the top 150 significant GO terms 
     #Adjusted to look just at the classic fisher results since we have just counts
       allRes.all <- GenTable(all_GO, fisher = resultFisher.all, ranksOf = "fisher", topNodes = 150, numChar = 500)
       colnames(allRes.all)[6] <- "p-value"
       allRes.all$`p-value` <- as.numeric(allRes.all$`p-value`)
       
     # Filter nodes with pvalue greater than 0.05 
       allRes.all <- allRes.all[allRes.all$`p-value` < 0.05,]
       
       
     #Extract names of significant genes in GenTable result. Add to column in allRes.all
       allRes.all$genes <- sapply(allRes.all$GO.ID, function(x)
       {
         genes<-genesInTerm(all_GO, x)
         genes[[1]][genes[[1]] %in% GeneInt.all] # myGenes is the queried gene list
       })
       
     # Create and export enrichment result as a formatted and interactable table
       all.allsigs.table <- reactable(allRes.all[,1:6], defaultPageSize = 150, theme = pff(centered = FALSE, font_color = "black"), 
                                           wrap = FALSE, bordered = TRUE, compact = TRUE, striped = TRUE, highlight = TRUE, fullWidth = TRUE,
                                           columns = list(Term = colDef(minWidth = 280),
                                                          GO.ID = colDef(cell = pill_buttons(colors = "darkgreen"), minWidth = 120)
                                           )
       ) %>%
         add_title("All genes - all significant genes", font_size = 20) %>% 
         add_subtitle("443 of 728 significant genes", font_size = 16, font_style = "italic") 
       
       
       all.allsigs.table #View table
       
     #Save table as an html file
       save_reactable_test(all.allsigs.table, "All sig genes table.html")
       
     # investigate how the significant GO terms are distributed over the GO graph
       jpeg(filename = "TopGO subgraph all genes.jpg", width = 10, height = 10, units = "in", res = 900)
       par(cex = 0.9)
       showSigOfNodes(all_GO, score(resultFisher.all), firstSigNodes = 25, useInfo = 'all')
       dev.off()   
       
       
  # Make sure you load the T.californicus OrgDb package before doing this. 
  # Cluster 1 from 
      #Create sim matrix using method Wang since it is most recent and keytype GID since we dont have ENTREZIDs
       simMatrix.all <- calculateSimMatrix(allRes.all$GO.ID,
                                                orgdb = Tcalif_orgdb_object,
                                                ont="BP",
                                                method="Wang", 
                                                keytype = "GID")
       
      #Create groupings of reduced terms for easier visualization
       reducedTerms.all <- reduceSimMatrix(simMatrix.all,
                                                threshold=0.8,
                                                orgdb=Tcalif_orgdb_object,
                                                keytype = "GID")
       
      # Make scatter plot depicting groups and distance between terms  
       jpeg(filename = "rrvgo scatter all.jpg", width = 16, height = 11, units = "in", res = 300)
       scatterPlot(simMatrix.all, reducedTerms.all)
       dev.off()
       
      # Make treemap of GO terms clustered under their parent terms
       jpeg(filename = "rrvgo treemap all.jpg", width = 10, height = 6, units = "in", res = 300)
       treemapPlot(reducedTerms.all)
       dev.off()
       
       
   # Explore genes associated with different GO terms
    # Send genes of interest to a list (right now just overwriting it each time before 
    # going to command below to pull them out of sig_genes_by_qtl
       
       # GO:0006112 energy reserve metabolic process
       int.list <- allRes.all$genes$`GO:0006112`
       
       # GO:0043462 regulation of ATP-dependent activity
       int.list <- allRes.all$genes$`GO:0043462`
       
       # GO:0010632 regulation of epithelial cell migration
       int.list <- allRes.all$genes$`GO:0010632`
       
       # GO:0043406 positive regulation of MAP kinase activity
       int.list <- allRes.all$genes$`GO:0043406`
       
       # GO:0070585 protein localization to mitochondrion
       int.list <- allRes.all$genes$`GO:0070585`
       
       # GO:0071333 cellular response to glucose stimulus
       int.list <- allRes.all$genes$`GO:0071333`
       
       # GO:0080164 regulation of nitric oxide metabolic process
       int.list <- allRes.all$genes$`GO:0080164`
       
       # GO:1901135 carbohydrate derivative metabolic process
       int.list <- allRes.all$genes$`GO:1901135`
       
       # GO:0009651 response to salt stress
       int.list <- allRes.all$genes$`GO:0009651`
       
       # GO:0009719 response to endogenous stimulus
       int.list <- allRes.all$genes$`GO:0009719`
       
       # GO:0007041 lysosomal transport
       int.list <- allRes.all$genes$`GO:0007041`
       
      # Pull out parts of the combined_sigs master results data frame that correspond to go terms of interest saved from above
       gene_summary[gene_summary$Gene_ID %in% int.list, c(1,7)]
       
       
###### Exploring overlap of genes from QTL with functional groups ########
    # Read in cuticle genes
      # Read in data and reformat gene names to remove PA notation if present
       cuticles <- read.csv(file = "Functional groups lists/Chitin_and_cuticle_processes_gene_list_SDv2.2_noPA.csv", header = TRUE)
       names(cuticles)[1] <- "Gene_ID"
       cuticles$tip <- paste0(cuticles$Gene_ID,": ", cuticles$Description)
       
      # Add coloring instructions for lines based on whether gene is in significant list
       cuticles$sig <- NA
       for(i in 1:nrow(cuticles)) {
         if (cuticles[i,1] %in% gene_summary$Gene_ID) {
           cuticles[i,7] <- "sig"
         } else {
           cuticles[i,7] <- "nonsig"
         }
       }
       
       # Enrichment test
       # Minus the proteins of unknown function (178 of these in sig genes and 3901 of these in blast annotation total)
       # So there are 543 sig genes to work with and 11676 genes in annotated universe
       # We create a matrix(c(a, b, c, d)), 
       # where a = significant genes in our set that overlap with genes in cuticle group
       # b = significant genes in our set that don't overlap
       # c = cuticle genes not overlapping with our set
       # d = all other non-sig genes from either study (the rest of the gene universe)
       # Run fisher test to get Odds ratio
       cuticles_overlap <- matrix(c(11, 532, 212, 10921), nrow =2, byrow = TRUE)
       fisher.test(cuticles_overlap)
       
       
      # Read in data 
       glycolysis <- read.csv(file = "Functional groups lists/Glycolysis_and_related_processes_gene_list_SDv2.2_noPA.csv", header = TRUE)
       names(glycolysis)[1] <- "Gene_ID"
       glycolysis$tip <- paste0(glycolysis$Gene_ID,": ", glycolysis$Description)
      
      # Add coloring instructions for lines based on whether gene is in significant list
       glycolysis$sig <- NA
       for(i in 1:nrow(glycolysis)) {
         if (glycolysis[i,1] %in% gene_summary$Gene_ID) {
           glycolysis[i,8] <- "sig"
         } else {
           glycolysis[i,8] <- "nonsig"
         }
       }
       
       # Enrichment test
       # Minus the proteins of unknown function (178 of these in sig genes and 3901 of these in blast annotation total)
       # So there are 543 sig genes to work with and 11676 genes in annotated universe
       # We create a matrix(c(a, b, c, d)), 
       # where a = significant genes in our set that overlap with genes in glycolysis group
       # b = significant genes in our set that don't overlap
       # c = glycolysis genes not overlapping with our set (out of 196 unique genes in list)
       # d = all other non-sig genes from either study (the rest of the gene universe)
       # Run fisher test to get Odds ratio
       glycolysis_overlap <- matrix(c(5, 538, 191, 10942), nrow =2, byrow = TRUE)
       fisher.test(glycolysis_overlap)
       
    # Antioxidant related genes
       
      # Read in data
       antioxidants <- read.csv(file = "Functional groups lists/Antioxidant_genes_noPA.csv", header = TRUE)
       names(antioxidants)[1] <- "Gene_ID" # redundant for this file I believe but its fine
       antioxidants$tip <- paste0(antioxidants$Gene_ID,": ", antioxidants$Description)
       
      # Add coloring instructions for lines based on whether gene is in significant list
       antioxidants$sig <- NA
       for(i in 1:nrow(antioxidants)) {
         if (antioxidants[i,1] %in% gene_summary$Gene_ID) {
           antioxidants[i,7] <- "sig"
         } else {
           antioxidants[i,7] <- "nonsig"
         }
       }
       
       
    # Pigment processing genes
       
      # Read in data
       carotenoids <- read.csv(file = "Functional groups lists/Pigment_processes_gene_list_SDv2.2_noPA_final.csv", header = TRUE)
       names(carotenoids)[1] <- "Gene_ID"
       carotenoids$tip <- paste0(carotenoids$Gene_ID,": ", carotenoids$Description)
       
      # Add coloring instructions for lines based on whether gene is in significant list
       carotenoids$sig <- NA
       for(i in 1:nrow(carotenoids)) {
         if (carotenoids[i,1] %in% gene_summary$Gene_ID) {
           carotenoids[i,8] <- "sig"
         } else {
           carotenoids[i,8] <- "nonsig"
         }
       }
       
       
    # Mito targeted genes
      
      # Read in data 
       mitos <- read.csv(file = "Functional groups lists/Mito-target_proteins_SDv2.2_June2025.csv", header = TRUE)
       names(mitos)[c(1,2,4)] <- c("Gene_ID", "Description", "Process")
       mitos$tip <- paste0(mitos$Gene,": ", mitos$Description)
       
     # Add coloring instructions for lines based on whether gene is in significant list
       mitos$sig <- NA
       for(i in 1:nrow(mitos)) {
         if (mitos[i,1] %in% gene_summary$Gene_ID) {
           mitos[i,6] <- "sig"
         } else {
           mitos[i,6] <- "nonsig"
         }
       }
       
       # Enrichment test
       # Minus the proteins of unknown function (178 of these in sig genes and 3901 of these in blast annotation total)
       # So there are 543 sig genes to work with and 11676 genes in annotated universe
       # We create a matrix(c(a, b, c, d)), 
       # where a = significant genes in our set that overlap with genes in mitos group
       # b = significant genes in our set that don't overlap
       # c = mito genes not overlapping with our set
       # d = all other non-sig genes from either study (the rest of the gene universe)
       # Run fisher test to get Odds ratio
       mito_overlap <- matrix(c(24, 518, 588, 10547), nrow =2, byrow = TRUE)
       fisher.test(mito_overlap)
       
       

   
       
       
###### dn/ds analysis results #####
  # Read in results from PAML
   omegas <- read.xlsx("Genes for dnds analysis/paml_tableOUT.xlsx", colNames = TRUE)
   
  # Remove any instance of dn/ds being 99 or 0
   omegas_filt <- omegas[omegas$omega != 99 & omegas$omega != 0,]
# check group sizes
   table(omegas_filt$group)
 # check distribution for fun
   hist(omegas_filt$omega, breaks = 100)   # Yep, right skewed. 
   
  # Run glm with gamma distribution 
   dnds.mod <- glm(omega ~ group, family = Gamma(link = "log"), data = omegas_filt)
   summary(dnds.mod)
   exp(-0.09353)
   
 # Generate boxplot
   dnds.plot <- ggplot(omegas_filt, aes(x = group, y = omega, fill = group))+
    geom_boxplot()+
    scale_y_log10()+
     scale_fill_manual(values = c("lightsteelblue", "limegreen"))+
     scale_x_discrete(labels = c("Background genes", "Genes in QTLs"))+
     ylab("dN/dS (log scaled)")+xlab("")+
     theme(legend.position = "none")+
     labs(subtitle = expression(OR == 0.91 ~ "," ~ italic(p) == 0.343))
  dnds.plot   
  
  
  
  # Export
  jpeg("dn_ds_boxplot.jpg", width = 5, height = 5, units = 'in', res = 300)
  dnds.plot
  dev.off()
  
 # Merge dn/ds values into gene summary table
  # Adjust sequence name to match my naming scheme (remove the -PA and the IF from TCALIF)
  omegas$gene <- str_split_i(omegas$gene, "-", 1)
  
  omegas <- omegas %>% dplyr::rename(Gene_ID = gene)
  
  gene_summary <- merge(gene_summary, omegas[,c(1,3,4,5)], by = "Gene_ID", )
  
 
###### Comparison to time series ######
  
  # First, basic Enrichmet test for gene overlap (OR - fishers)
  # We create a matrix(c(a, b, c, d)), 
  # where a = significant genes in our set that overlap with genes in time series study
  # b = significant genes in our set that don't overlap
  # c = time series study genes not overlapping with our set
  # d = all other non-sig genes from either study (the rest of the gene universe)
  # Run fisher test to get Odds ratio
  timeseries_overlap <- matrix(c(101, 620, 1797, 13059), nrow =2, byrow = TRUE)
  fisher.test(timeseries_overlap)
  
  
  # Next, we make paired distributions of data for time series and non-time series genes
    # Will also run models to see if distributions are different on average
    
      # Define color palette
        pal_1 <- mex.brewer("Revolucion")
        show_col(pal_1)
        my_pal_ts <- pal_1[c(3,8)]
        show_col(my_pal_ts)
        
      # Make summary data frame for overlaying mean's and CIs
        # Make sqrt transformed variable for n_in_gene_per_kb
        gene_summary$n_in_gene_per_kb_tf <- sqrt(gene_summary$n_in_gene_per_kb)
        
        # Define our groups
        ridge_vars <- c(
          "n_upstream",
          "n_in_gene",
          "n_introns",
          "n_in_exons",
          "gene_length_kb",
          "n_in_gene_per_kb_tf",
          "n_modifer",
          "n_low_impact",
          "n_mod_impact",
          "n_high_impact"
        )
        
        # Pivot long to make summary
        ridge_long <- gene_summary %>%
          dplyr::select(In_time_series, all_of(ridge_vars)) %>%
          pivot_longer(
            cols = all_of(ridge_vars),
            names_to = "metric",
            values_to = "value"
          ) %>%
          filter(!is.infinite(value) & !is.na(value))
        
        # Generate summary
        ridge_summary <- ridge_long %>%
          group_by(metric, In_time_series) %>%
          summarise(
            n = sum(!is.na(value)),
            mean_val = mean(value, na.rm = TRUE),
            sd_val = sd(value, na.rm = TRUE),
            se = sd_val / sqrt(n),
            ci_low = mean_val - qt(0.975, df = n - 1) * se,
            ci_high = mean_val + qt(0.975, df = n - 1) * se,
            .groups = "drop"
          )
        
        # Define function to auto call the correct elements when plotting
        get_ridge_ci <- function(metric_name) {
          ridge_summary %>%
            filter(metric == metric_name)
        }
        
      
        
      # Upstream genes
        # Run glm with poisson distribution 
        upstream.mod <- glm(n_upstream ~ In_time_series + n_in_gene_per_kb, family = quasipoisson, data = gene_summary)
        summary(upstream.mod)
        exp(0.053168)
        # Make plot
        ridge_upstream <- ggplot(gene_summary, aes(x = n_upstream, y = In_time_series, fill = In_time_series)) +
          geom_density_ridges(alpha = 0.6, scale = 2, color = "white", bandwidth = 2, rel_min_height =0.01,
                              jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0), 
                              point_shape = "|", point_size = 4, point_color = "grey40") +
          geom_errorbar(data = get_ridge_ci("n_upstream"), aes(y = In_time_series, xmin = ci_low,xmax = ci_high), 
            orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black")+
          geom_point(data = get_ridge_ci("n_upstream"), aes(x = mean_val, y = In_time_series),
                     inherit.aes = FALSE,size = 3,shape = 21,fill = "red",color = "black")+
          labs(x = "Number of SNPs 1000 bp upstream", y = "", title = "Distribution of Upstream SNP Counts") +
          scale_fill_manual(values = my_pal_ts)+
          scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE"  = "DE"))+
          coord_cartesian(xlim = c(0, max(gene_summary$n_upstream))) +
          annotate("text", x = 28, y = 2.5, size = 5,
                   label = "OR==1.06 ~ ';' ~ italic(p)==0.632",
                   parse = TRUE)+
          # annotate("text", x = 28, y = 3.1, size = 3,
          #          label = "# SNPs upstream ~ Time Series overlap (Yes/No) + # SNPs per kb in gene",
          #          parse = FALSE)+
          theme(legend.position = "none")
        ridge_upstream
        
      # In gene
        # Run glm with poisson distribution 
        ingene.mod <- glm(n_in_gene ~ In_time_series + offset(log(gene_length_kb)), family = quasipoisson, data = gene_summary)
        summary(ingene.mod)
        exp(0.18906)
        # Make plot
          ridge_in_gene <- ggplot(gene_summary, aes(x = n_in_gene, y = In_time_series, fill = In_time_series)) +
            geom_density_ridges(alpha = 0.6, scale = 2, color = "white", bandwidth = 5, rel_min_height =0.01,
                                jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0), 
                                point_shape = "|", point_size = 4, point_color = "grey40") +
            geom_errorbar(data = get_ridge_ci("n_in_gene"), aes(y = In_time_series, xmin = ci_low,xmax = ci_high),
                          orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black")+
            geom_point(data = get_ridge_ci("n_in_gene"), aes(x = mean_val, y = In_time_series),
                       inherit.aes = FALSE,size = 3,shape = 21,fill = "red",color = "black")+
            labs(x = "Number of SNPs in Genes", y = "", title = "Distribution of SNP counts in genes") +
            scale_fill_manual(values = my_pal_ts)+
            scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE"  = "DE"))+
            coord_cartesian(xlim = c(0, max(gene_summary$n_in_gene))) +
            annotate("text", x = 250, y = 2.5, size = 5,
                     label = "OR==1.21 ~ ';' ~ italic(p)<0.001",
                     parse = TRUE)+
            # annotate("text", x = 250, y = 3, size = 3,
            #          label = "# SNPs in gene ~ Time Series overlap (Yes/No) + offset by gene length",
            #          parse = FALSE)+
            theme(legend.position = "none")
          ridge_in_gene
          
          
        # In gene per kb
          ingeneperkb.mod <- lm(n_in_gene_per_kb_tf ~ In_time_series, data = gene_summary)
          summary(ingeneperkb.mod)
          plot(ingeneperkb.mod)
          # Make plot
          ridge_in_gene_per_kb <- ggplot(gene_summary, aes(x = n_in_gene_per_kb_tf, y = In_time_series, fill = In_time_series)) +
            geom_density_ridges(alpha = 0.6, scale = 2, color = "white", bandwidth = 0.2, rel_min_height =0.01,
                                jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0), 
                                point_shape = "|", point_size = 4, point_color = "grey40") +
            geom_errorbar(data = get_ridge_ci("n_in_gene_per_kb_tf"), aes(y = In_time_series, xmin = ci_low,xmax = ci_high), 
                          orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black")+
            geom_point(data = get_ridge_ci("n_in_gene_per_kb_tf"), aes(x = mean_val, y = In_time_series),
                       inherit.aes = FALSE,size = 3,shape = 21,fill = "red",color = "black")+
            labs(x = "Number of SNPs in genes per kb", y = "", title = "Distribution of SNP counts in genes per kb") +
            scale_fill_manual(values = my_pal_ts)+
            scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE"  = "DE"))+
            coord_cartesian(xlim = c(0, max(gene_summary$n_in_gene_per_kb_tf))) +
            annotate("text", x = 4.4, y = 2.9, size = 5,
                     label = "italic(beta)==0.29 ~ ';' ~ italic(p)==0.003",
                     parse = TRUE)+
            # annotate("text", x = 4.4, y = 3.4, size = 3,
            #          label = "# SNPs per kb in gene ~ Time Series overlap (Yes/No)",
            #          parse = FALSE)+
            theme(legend.position = "none")
          ridge_in_gene_per_kb
          
      # Gene length
          # Run glm with poisson distribution 
          length.mod <- glm(gene_length_kb ~ In_time_series, family = Gamma(link = "log"), data = gene_summary)
          summary(length.mod)
          exp(0.17360)
          # Make plot
          ridge_gene_length <- ggplot(gene_summary, aes(x = gene_length_kb, y = In_time_series, fill = In_time_series)) +
            geom_density_ridges(alpha = 0.6, scale = 2, color = "white", bandwidth = 5, rel_min_height =0.01,
                                jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0), 
                                point_shape = "|", point_size = 4, point_color = "grey40") +
            geom_errorbar(data = get_ridge_ci("gene_length_kb"), aes(y = In_time_series, xmin = ci_low,xmax = ci_high), 
                          orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black")+
            geom_point(data = get_ridge_ci("gene_length_kb"), aes(x = mean_val, y = In_time_series),
                       inherit.aes = FALSE,size = 3,shape = 21,fill = "red",color = "black")+
            labs(x = "Gene Length in kb", y = "", title = "Distribution of Gene Lengths") +
            scale_fill_manual(values = my_pal_ts)+
            scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE"  = "DE"))+
            coord_cartesian(xlim = c(0, max(gene_summary$gene_length_kb))) +
            annotate("text", x = 75, y = 2.5, size = 5,
                     label = "OR==1.19 ~ ';' ~ italic(p)==0.187",
                     parse = TRUE)+
            # annotate("text", x = 75, y = 3, size = 3,
            #          label = "Gene Length ~ Time Series overlap (Yes/No)",parse = FALSE)+
            theme(legend.position = "none")
          ridge_gene_length
      
      # In exons
          # Run glm with poisson distribution 
          inexon.mod <- glm(n_in_exons ~ In_time_series + offset(log(gene_length_kb)), family = quasipoisson, data = gene_summary)
          summary(inexon.mod)
          exp(-0.02967)
          # Make plot
          ridge_in_exon <- ggplot(gene_summary, aes(x = n_in_exons, y = In_time_series, fill = In_time_series)) +
            geom_density_ridges(alpha = 0.6, scale = 2, color = "white", bandwidth = 2, rel_min_height =0.01,
                                jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0), 
                                point_shape = "|", point_size = 4, point_color = "grey40") +
            geom_errorbar(data = get_ridge_ci("n_in_exons"), aes(y = In_time_series, xmin = ci_low,xmax = ci_high), 
                          orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black")+
            geom_point(data = get_ridge_ci("n_in_exons"), aes(x = mean_val, y = In_time_series),
                       inherit.aes = FALSE,size = 3,shape = 21,fill = "red",color = "black")+
            labs(x = "# SNPs in exons", y = "", title = "Distribution of SNP counts in Exons") +
            scale_fill_manual(values = my_pal_ts)+
            scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE"  = "DE"))+
            coord_cartesian(xlim = c(0, max(gene_summary$n_in_exons))) +
            annotate("text", x = 70, y = 2.5, size = 5,
                     label = "OR==0.97 ~ ';' ~ italic(p)==0.821",
                     parse = TRUE)+
            # annotate("text", x = 70, y = 3, size = 3,
            #          label = "# SNPs in exons ~ Time Series overlap (Yes/No) + offset by gene length",
            #          parse = FALSE)+
            theme(legend.position = "none")
          ridge_in_exon
          
          
        # In introns
          # Run glm with poisson distribution 
          inintron.mod <- glm(n_introns ~ In_time_series + offset(log(gene_length_kb)), family = quasipoisson, data = gene_summary)
          summary(inintron.mod)
          exp(0.27580)
          # Make plot
          ridge_in_intron <- ggplot(gene_summary, aes(x = n_introns, y = In_time_series, fill = In_time_series)) +
            geom_density_ridges(alpha = 0.6, scale = 2, color = "white", bandwidth = 4, rel_min_height =0.01,
                                jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0), 
                                point_shape = "|", point_size = 4, point_color = "grey40") +
            geom_errorbar(data = get_ridge_ci("n_introns"), aes(y = In_time_series, xmin = ci_low,xmax = ci_high), 
                          orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black")+
            geom_point(data = get_ridge_ci("n_introns"), aes(x = mean_val, y = In_time_series),
                       inherit.aes = FALSE,size = 3,shape = 21,fill = "red",color = "black")+
            labs(x = "# SNPs in introns", y = "", title = "Distribution of SNP counts in Introns") +
            scale_fill_manual(values = my_pal_ts)+
            scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE"  = "DE"))+
            coord_cartesian(xlim = c(0, max(gene_summary$n_introns))) +
            annotate("text", x = 220, y = 2.5, size = 5,
                     label = "OR==1.32 ~ ';' ~ italic(p)<0.001",
                     parse = TRUE)+
            # annotate("text", x = 218, y = 3, size = 3,
            #          label = "# of SNPs in introns ~ Time Series overlap (Yes/No) + offset by gene length",
            #          parse = FALSE)+
            theme(legend.position = "none")
          ridge_in_intron
      
     ridges_patch <- (ridge_in_gene_per_kb / ridge_in_exon / ridge_in_intron / ridge_upstream / ridge_gene_length) 
      
     ggsave(filename = "ridgeline_column_plot.png", plot = ridges_patch, width = 7, height = 12, dpi = 300)
     
     
     
  # G' comparison between time series groupings
     gprime.mod <- glm(max_Gprime ~ In_time_series, family = Gamma(link = "log"), data = gene_summary)
     plot(gprime.mod)
     summary(gprime.mod)
     exp(-0.002365)
     
     # Generate boxplot
     gprime.plot <- ggplot(gene_summary, aes(x = In_time_series, y = max_Gprime, fill = In_time_series))+
       geom_boxplot()+
       scale_y_log10()+
       scale_fill_manual(values = c("#d9792e", "#368990"))+
       scale_x_discrete(labels = c("Not DE", "DE"))+
       ylab("Maximum G' per gene")+xlab("")+
       theme(legend.position = "none", axis.text.y = element_text(size = 17), axis.text.x = element_text(size = 16),
             axis.title.y = element_text(size=19),
             plot.subtitle = element_text(size =17))+
       labs(subtitle = expression(OR == 1.00 ~ "," ~ italic(p) == 0.849))
     gprime.plot   
     
     # Export
     ggsave("gprime_boxplot_time_series.jpg", gprime.plot, width = 4, height = 4, dpi = 300)
     
     
     
  ## Do genes in time series tend to have higher ratio of SNPs upstream to in exons?
     mod.ratio.1 <- lm(log1p(up_to_in_per_kb_ratio) ~ In_time_series, data = gene_summary)
     summary(mod.ratio.1)
    
  # Do genes with a greater ratio of upstream to coding SNPs have higher or lower max G primes?
     mod.ratio.2 <- lm(log1p(up_to_in_per_kb_ratio) ~ max_Gprime, data = gene_summary)
     summary(mod.ratio.2) 
     
   
  # Do genes that overlap with the time series tend to have more of any of the types of impact categories for the SNPs?
     low.mod <- glm(n_low_impact ~ In_time_series + offset(log(gene_length_kb)), family = quasipoisson, data = gene_summary)
     summary(low.mod)
     exp(-0.02440)
     
     moderat.mod <- glm(n_mod_impact ~ In_time_series + offset(log(gene_length_kb)), family = quasipoisson, data = gene_summary)
     summary(moderat.mod)
     exp(-0.15444)
     
     high.mod <- glm(n_high_impact ~ In_time_series + offset(log(gene_length_kb)), family = quasipoisson, data = gene_summary)
     summary(high.mod)
     exp(0.3417)
     
     modifier.mod <- glm(n_modifer ~ In_time_series + offset(log(gene_length_kb)), family = quasipoisson, data = gene_summary)
     summary(modifier.mod)
     exp(0.19342)
     
     # Make plot for the mpdifier impact
     ridge_modifier_impact <- ggplot(gene_summary, aes(x = n_modifer, y = In_time_series, fill = In_time_series)) +
       geom_density_ridges(alpha = 0.6, scale = 2, color = "white", bandwidth = 5, rel_min_height =0.01,
                           jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0), 
                           point_shape = "|", point_size = 4, point_color = "grey40") +
       geom_errorbar(data = get_ridge_ci("n_modifer"), aes(y = In_time_series, xmin = ci_low,xmax = ci_high), 
                     orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black")+
       geom_point(data = get_ridge_ci("n_modifer"), aes(x = mean_val, y = In_time_series),
                  inherit.aes = FALSE,size = 3,shape = 21,fill = "red",color = "black")+
       labs(x = "# SNPs with 'modifier' impact", y = "", title = "Distribution of SNP counts with 'modifier' impact") +
       scale_fill_manual(values = my_pal_ts)+
       scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE"  = "DE"))+
       coord_cartesian(xlim = c(0, max(gene_summary$n_modifer))) +
       annotate("text", x = 220, y = 2.5, size = 5,
                label = "OR==1.21 ~ ';' ~ italic(p)==0.002",
                parse = TRUE)+
       # annotate("text", x = 218, y = 3, size = 3,
       #          label = "# of SNPs with modifier impact ~ Time Series overlap (Yes/No) + offset by gene length",
       #          parse = FALSE)+
       theme(legend.position = "none")
     ridge_modifier_impact
     
     # Make plot for the low impact
     ridge_low_impact <- ggplot(gene_summary, aes(x = n_low_impact, y = In_time_series, fill = In_time_series)) +
       geom_density_ridges(alpha = 0.6, scale = 2, color = "white", bandwidth = 2, rel_min_height =0.01,
                           jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0), 
                           point_shape = "|", point_size = 4, point_color = "grey40") +
       geom_errorbar(data = get_ridge_ci("n_low_impact"), aes(y = In_time_series, xmin = ci_low,xmax = ci_high), 
                     orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black")+
       geom_point(data = get_ridge_ci("n_low_impact"), aes(x = mean_val, y = In_time_series),
                  inherit.aes = FALSE,size = 3,shape = 21,fill = "red",color = "black")+
       labs(x = "# SNPs with low impact", y = "", title = "Distribution of SNP counts with low impact") +
       scale_fill_manual(values = my_pal_ts)+
       scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE"  = "DE"))+
       coord_cartesian(xlim = c(0, max(gene_summary$n_low_impact))) +
       annotate("text", x = 55, y = 2.5, size = 5,
                label = "OR==0.98 ~ ';' ~ italic(p)==0.866",
                parse = TRUE)+
       # annotate("text", x = 55, y = 3, size = 3,
       #          label = "# of SNPs with low impact ~ Time Series overlap (Yes/No) + offset by gene length",
       #          parse = FALSE)+
       theme(legend.position = "none")
     ridge_low_impact
     
     # Make plot for the moderate impact
     ridge_moderate_impact <- ggplot(gene_summary, aes(x = n_mod_impact, y = In_time_series, fill = In_time_series)) +
       geom_density_ridges(alpha = 0.6, scale = 2, color = "white", bandwidth = 2, rel_min_height =0.01,
                           jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0), 
                           point_shape = "|", point_size = 4, point_color = "grey40") +
       geom_errorbar(data = get_ridge_ci("n_mod_impact"), aes(y = In_time_series, xmin = ci_low,xmax = ci_high), 
                     orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black")+
       geom_point(data = get_ridge_ci("n_mod_impact"), aes(x = mean_val, y = In_time_series),
                  inherit.aes = FALSE,size = 3,shape = 21,fill = "red",color = "black")+
       labs(x = "# SNPs with moderate impact", y = "", title = "Distribution of SNP counts with moderate impact") +
       scale_fill_manual(values = my_pal_ts)+
       scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE"  = "DE"))+
       coord_cartesian(xlim = c(0, max(gene_summary$n_mod_impact))) +
       annotate("text", x = 24, y = 2.5, size = 5,
                label = "OR==0.86 ~ ';' ~ italic(p)==0.471",
                parse = TRUE)+
       # annotate("text", x = 24, y = 3, size = 3,
       #          label = "# of SNPs with moderate impact ~ Time Series overlap (Yes/No) + offset by gene length",
       #          parse = FALSE)+
       theme(legend.position = "none")
     ridge_moderate_impact
     
     # Make plot for the high impact
     ridge_high_impact <- ggplot(gene_summary, aes(x = n_high_impact, y = In_time_series, fill = In_time_series)) +
       geom_density_ridges(alpha = 0.6, scale = 2, color = "white", bandwidth = 0.05, rel_min_height =0.01,
                           jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0), 
                           point_shape = "|", point_size = 4, point_color = "grey40") +
       geom_errorbar(data = get_ridge_ci("n_high_impact"), aes(y = In_time_series, xmin = ci_low,xmax = ci_high), 
                     orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black")+
       geom_point(data = get_ridge_ci("n_high_impact"), aes(x = mean_val, y = In_time_series),
                  inherit.aes = FALSE,size = 3,shape = 21,fill = "red",color = "black")+
       labs(x = "# SNPs with high impact", y = "", title = "Distribution of SNP counts with high impact") +
       scale_fill_manual(values = my_pal_ts)+
       scale_y_discrete(labels = c("FALSE" = "Not DE", "TRUE"  = "DE"))+
       coord_cartesian(xlim = c(0, max(gene_summary$n_high_impact))) +
       annotate("text", x = 1.4, y = 2.5, size = 5,
                label = "OR==1.41 ~ ';' ~ italic(p)==0.637",
                parse = TRUE)+
       # annotate("text", x = 1.4, y = 3, size = 3,
       #          label = "# of SNPs with high impact ~ Time Series overlap (Yes/No) + offset by gene length",
       #          parse = FALSE)+
       theme(legend.position = "none")
     ridge_high_impact
     
     ridges_patch_impact <- (ridge_modifier_impact / ridge_low_impact / ridge_moderate_impact / ridge_high_impact)
     
     ggsave(filename = "ridgeline_column_plot_impacts.png", plot = ridges_patch_impact, width = 7, height = 10, dpi = 300)
     
     
 
  # Different GO processes? 
     
     # In time series 
     #Get first column from  results w/ gene names and pvalues. Send to list. 
     GeneInt.ts <- gene_summary %>%
       filter(In_time_series == "TRUE") %>%
       pull(Gene_ID) 
     
     
     #Make named vector with the gene universe where genes of interest are coded with a 1 so the 'new' function knows to focus on them. 
     geneList.ts <- factor(as.integer(geneUniverse %in% GeneInt.ts)) 
     names(geneList.ts) <- geneUniverse
     
     # We now have all data necessary to build an object of type topGOdata. This object will contain all gene
     # identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information
     # needed to perform the desired enrichment analysis.
     ts_GO <- new("topGOdata", description = "GO analysis of all sig genes", ontology = "BP", nodeSize = 10,
                  allGenes = geneList.ts, annot = annFUN.gene2GO, gene2GO = geneID2GO)
     
     # Once we have an object of class topGOdata we can start with the enrichment analysis.
     # Run the classic fisher exact test to find enriched go terms
     resultFisher.ts <- runTest(ts_GO, algorithm = "weight01", statistic = "fisher")
     
     resultFisher.ts #View results summary
     
     # GenTable is an easy to use function for analysing the most significant GO terms and the corresponding p
     # values. We list the top 150 significant GO terms.
     #Adjusted to look just at the classic fisher results since we have just counts
     allRes.ts <- GenTable(ts_GO, fisher = resultFisher.ts, ranksOf = "fisher", topNodes = 150, numChar = 500)
     colnames(allRes.ts)[6] <- "p-value"
     allRes.ts$`p-value` <- as.numeric(allRes.ts$`p-value`)
     
     # Filter nodes with pvalue greater than 0.05 
     allRes.ts <- allRes.ts[allRes.ts$`p-value` < 0.05,]
     
     
     #Extract names of significant genes in GenTable result. Add to column in allRes.ts
     allRes.ts$genes <- sapply(allRes.ts$GO.ID, function(x)
     {
       genes<-genesInTerm(ts_GO, x)
       genes[[1]][genes[[1]] %in% GeneInt.ts] # myGenes is the queried gene list
     })
     
     # Create and export enrichment result as a formatted and interactable table
     ts.sigs.table <- reactable(allRes.ts[,1:6], defaultPageSize = 150, theme = pff(centered = FALSE, font_color = "black"), 
                                wrap = FALSE, bordered = TRUE, compact = TRUE, striped = TRUE, highlight = TRUE, fullWidth = TRUE,
                                columns = list(Term = colDef(minWidth = 280),
                                               GO.ID = colDef(cell = pill_buttons(colors = "darkgreen"), minWidth = 120)
                                )
     ) %>%
       add_title("QTL genes in DE dataset", font_size = 20) %>% 
       add_subtitle("55 of 9701 genes post filtering", font_size = 16, font_style = "italic") 
     
     
     ts.sigs.table #View table
     
     #Save table as an html file
     save_reactable_test(ts.sigs.table, "QTL genes in time series.html")
     
     # investigate how the significant GO terms are distributed over the GO graph
     jpeg(filename = "TopGO subgraph ts genes.jpg", width = 10, height = 10, units = "in", res = 900)
     par(cex = 0.9)
     showSigOfNodes(ts_GO, score(resultFisher.ts), firstSigNodes = 25, useInfo = 'all')
     dev.off()   
     
     
     # Make sure you load the T.californicus OrgDb package before doing this. 
     # QTL genes in time series
     #Create sim matrix using method Wang since it is most recent and keytype GID since we dont have ENTREZIDs
     simMatrix.ts <- calculateSimMatrix(allRes.ts$GO.ID,
                                        orgdb = Tcalif_orgdb_object,
                                        ont="BP",
                                        method="Wang", 
                                        keytype = "GID")
     
     #Create groupings of reduced terms for easier visualization
     reducedTerms.ts <- reduceSimMatrix(simMatrix.ts,
                                        threshold=0.8,
                                        orgdb=Tcalif_orgdb_object,
                                        keytype = "GID")
     
     # Make scatter plot depicting groups and distance between terms  
     jpeg(filename = "rrvgo scatter ts.jpg", width = 16, height = 11, units = "in", res = 300)
     scatterPlot(simMatrix.ts, reducedTerms.ts)
     dev.off()
     
     # Make treemap of GO terms clustered under their parent terms
     jpeg(filename = "rrvgo treemap ts.jpg", width = 9, height = 6, units = "in", res = 300)
     treemapPlot(reducedTerms.ts, fontsize.labels = c(18, 12),  # c(group label, item label)
                 fontsize.title = 20, title = "QTL genes in DE dataset")
     dev.off()
     
     # Query GO terms for gene IDs in that term
     # GO:0007041 lysosomal transport
     int.list <- allRes.ts$genes$`GO:0120035`
     # Pull out parts of the combined_sigs master results data frame that correspond to go terms of interest saved from above
     gene_summary[gene_summary$Gene_ID %in% int.list, c(2,3)]
     
     
    # QTL genes not in time series
     #Get first column from  results w/ gene names and pvalues. Send to list. 
     GeneInt.nits <- gene_summary %>%
       filter(In_time_series == "FALSE") %>%
       pull(Gene_ID) 
     
     #Make named vector with the gene universe where genes of interest are coded with a 1 so the 'new' function knows to focus on them. 
     geneList.nits <- factor(as.integer(geneUniverse %in% GeneInt.nits)) 
     names(geneList.nits) <- geneUniverse
     
     # We now have all data necessary to build an object of type topGOdata. This object will contain all gene
     # identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information
     # needed to perform the desired enrichment analysis.
     nits_GO <- new("topGOdata", description = "GO analysis of genes not in time series", ontology = "BP", nodeSize = 10,
                  allGenes = geneList.nits, annot = annFUN.gene2GO, gene2GO = geneID2GO)
     
     # Once we have an object of class topGOdata we can start with the enrichment analysis.
     # Run the classic fisher exact test to find enriched go terms
     resultFisher.nits <- runTest(nits_GO, algorithm = "weight01", statistic = "fisher")
     
     resultFisher.nits #View results summary
     
     # GenTable is an easy to use function for analysing the most significant GO terms and the corresponding p
     # values. We list the top 10 significant GO terms
     #Adjusted to look just at the classic fisher results since we have just counts
     allRes.nits <- GenTable(nits_GO, fisher = resultFisher.nits, ranksOf = "fisher", topNodes = 150, numChar = 500)
     colnames(allRes.nits)[6] <- "p-value"
     allRes.nits$`p-value` <- as.numeric(allRes.nits$`p-value`)
     
     # Filter nodes with pvalue greater than 0.05 
     allRes.nits <- allRes.nits[allRes.nits$`p-value` < 0.05,]
     
     
     #Extract names of significant genes in GenTable result. Add to column in allRes.nits
     allRes.nits$genes <- sapply(allRes.nits$GO.ID, function(x)
     {
       genes<-genesInTerm(nits_GO, x)
       genes[[1]][genes[[1]] %in% GeneInt.nits] # myGenes is the queried gene list
     })
     
     # Create and export enrichment result as a formatted and interactable table
     nits.sigs.table <- reactable(allRes.nits[,1:6], defaultPageSize = 150, theme = pff(centered = FALSE, font_color = "black"), 
                                wrap = FALSE, bordered = TRUE, compact = TRUE, striped = TRUE, highlight = TRUE, fullWidth = TRUE,
                                columns = list(Term = colDef(minWidth = 280),
                                               GO.ID = colDef(cell = pill_buttons(colors = "darkgreen"), minWidth = 120)
                                )
     ) %>%
       add_title("QTL genes not in DE dataset", font_size = 20) %>% 
       add_subtitle("388 of 9701 genes post filtering", font_size = 16, font_style = "italic") 
     
     
     nits.sigs.table #View table
     
     #Save table as an html file
     save_reactable_test(nits.sigs.table, "QTL genes not in time series.html")
     
     # investigate how the significant GO terms are distributed over the GO graph
     jpeg(filename = "TopGO subgraph nits genes.jpg", width = 10, height = 10, units = "in", res = 900)
     par(cex = 0.9)
     showSigOfNodes(nits_GO, score(resultFisher.nits), firstSigNodes = 25, useInfo = 'all')
     dev.off()   
     
     # Make sure you load the T.californicus OrgDb package before doing this. 
     # QTL genes in time series
     #Create sim matrix using method Wang since it is most recent and keytype GID since we dont have ENTREZIDs
     simMatrix.nits <- calculateSimMatrix(allRes.nits$GO.ID,
                                        orgdb = Tcalif_orgdb_object,
                                        ont="BP",
                                        method="Wang", 
                                        keytype = "GID")
     
     #Create groupings of reduced terms for easier visualization
     reducedTerms.nits <- reduceSimMatrix(simMatrix.nits,
                                        threshold=0.8,
                                        orgdb=Tcalif_orgdb_object,
                                        keytype = "GID")
     
     # Make scatter plot depicting groups and distance between terms  
     jpeg(filename = "rrvgo scatter nits.jpg", width = 16, height = 11, units = "in", res = 300)
     scatterPlot(simMatrix.nits, reducedTerms.nits)
     dev.off()
     
     # Make treemap of GO terms clustered under their parent terms
     jpeg(filename = "rrvgo treemap nits.jpg", width = 9, height = 6, units = "in", res = 300)
     treemapPlot(reducedTerms.nits, fontsize.labels = c(18, 14),  # c(group label, item label)
                 fontsize.title = 20, title = "QTL genes not in DE dataset")
     dev.off()
    
     
     # Query GO terms for gene IDs in that term
     # GO:0007041 lysosomal transport
     int.list <- allRes.nits$genes$`GO:0015980`
     # Pull out parts of the combined_sigs master results data frame that correspond to go terms of interest saved from above
     gene_summary[gene_summary$Gene_ID %in% int.list, c(2,3)]
     
     
  # dnds difference between two groups? 
     # Run glm with gamma distribution 
     dnds.ts.mod <- glm(omega ~ In_time_series, family = Gamma(link = "log"), data = subset(gene_summary, omega != 99 & omega != 0))
     summary(dnds.ts.mod)
     exp(-0.06498)
     
     # Generate boxplot
     dnds.ts.plot <- ggplot(subset(gene_summary, omega != 99 & omega != 0), aes(x = In_time_series, y = omega, fill = In_time_series))+
       geom_boxplot()+
       scale_y_log10()+
       scale_fill_manual(values = my_pal_ts)+
       scale_x_discrete(labels = c("FALSE" = "Not DE", "TRUE"  = "DE"))+
       ylab("dN/dS (log scaled)")+xlab("")+
       theme(legend.position = "none", axis.text.y = element_text(size = 17), axis.text.x = element_text(size = 16),
             axis.title.y = element_text(size=19),
             plot.subtitle = element_text(size =17))+
       labs(subtitle = expression(OR == 0.94 ~ "," ~ italic(p) == 0.630))
     dnds.ts.plot   
     
    
     # Export
     ggsave("dnds_by_expression_overlap.jpg", dnds.ts.plot, height = 4, width = 4, dpi = 300)
     
     
     
###### Which allele (SD or SH) increased in frequency in low-Pcrit group? ######
   # Start with using our all.snps.comb data frame because it has the high and low info
     qtl_snps_subset <- all.snps.comb[,c(2:8,10:12)]
     
     colnames(qtl_snps_subset)[5:10] <- c("sd_low", "sh_low", "total_low", "sd_high", "sh_high", "total_high")
     
     library(tidyverse)
     
    # Pivot subset to long format for plotting
     
   qtl_snps_subset_long <- qtl_snps_subset %>%
     dplyr::select(-REF, -ALT) %>% # First select all columns except the ref and alt alleles
     # pivot all cols expect CHROM and POS, separate col names to new summary columns 
     pivot_longer(cols = -c(CHROM, POS), names_to  = c("population", "treatment"), names_sep = "_") %>% 
     # Separate back out the data from individual pops so we can calculate allele frequency
     pivot_wider(names_from  = population, values_from = value) %>%
     # Now that we calculated the frequency, go longer again for plotting by putting the pops back in their own column 
     # but NOT the total, so it gets duplicated so each pop has the total associated with it
     pivot_longer(cols = c(sd, sh), names_to  = "population", values_to = "count") %>%
     mutate(frequency = count / total, population = toupper(population))  # Use toupper to make sd, sh to SD, SH
   
  # Realized after the fact I could use the QTL number for plotting purposes
   qtl_snps_subset_long <- qtl_snps_subset_long %>%
     left_join(all.snps.comb %>% dplyr::select(CHROM, POS, qtl), by = c("CHROM", "POS"))
     
   br_pal <- met.brewer("Signac")
   br_pal2 <- met.brewer("Renoir")
   my_pal2 <- c(br_pal2[8], br_pal[3], br_pal2[c(2,12)])
   show_col(my_pal2)
  
   
  # Make plot where we group by population 
   AF_by_QTL_and_bulk_plot <- qtl_snps_subset_long %>%
     dplyr::mutate(POS_Mb = POS / 1e6) %>%
     ggplot(aes(x = POS_Mb, y = frequency, color = population)) +
     geom_point(alpha = 0.2, size = 0.4) +
     geom_smooth(method = "loess", se = FALSE, aes(group = interaction(population, qtl)), color = "black", linewidth = 1.5) +
     geom_smooth(method = "loess", se = TRUE, linewidth = 1, aes(group = interaction(population, qtl))) +
     scale_color_manual(values = c("gold", "deepskyblue"))+
     facet_grid(treatment ~ CHROM, scales = "free_x", space = "free_x") +
     labs(x = "Position (Mb)", y = "Allele Frequency", color = "Population") +
     theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text  = element_text(face = "bold"))+
     ggtitle(label = "Allele frequencies by population and bulk in chromosomes subsetted by QTL")
   AF_by_QTL_and_bulk_plot
   
   ggsave("AF_by_QTL_and_bulk_plot.jpg", plot = AF_by_QTL_and_bulk_plot, width = 13, height = 6, dpi = 300)
   
  # Make plot of just SD allele frequency
   AF_by_QTL_and_bulk_plot_SD_only <- qtl_snps_subset_long %>%
     filter(population == "SD") %>%
     dplyr::mutate(POS_Mb = POS / 1e6) %>%
     ggplot(aes(x = POS_Mb, y = frequency, color = treatment)) +
     geom_point(alpha = 0.2, size = 0.4) +
     geom_smooth(method = "loess", se = FALSE, aes(group = interaction(treatment, qtl)), color = "black", linewidth = 1.5) +
     geom_smooth(method = "loess", se = TRUE, aes(group = interaction(treatment, qtl)), linewidth = 1) +
     scale_color_manual(values = c("gold", "forestgreen"))+
     facet_grid(~CHROM, scales = "free_x", space = "free_x") +
     labs(x = "Position (Mb)", y = "Allele Frequency", color = "Bulk") +
     theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text  = element_text(face = "bold"))+
     ggtitle(label = "SD allele frequencies by bulk in chromosomes subsetted by QTL")
   AF_by_QTL_and_bulk_plot_SD_only
   
   ggsave("AF_of_SD_by_bulk_plot.jpg", plot = AF_by_QTL_and_bulk_plot_SD_only, width = 13, height = 4, dpi = 300)
   
   
   
   ## Check on genome wide allele frequency in the low to compare to QTLs
   # Genome-wide background SD allele frequency in the low bulk
   mean_freq <- sum(df_filt2$AD_REF.HIGH, na.rm = TRUE) / sum(df_filt2$DP.HIGH, na.rm = TRUE)
   cat("Genome-wide SD allele frequency (low bulk):", round(mean_freq, 4), "\n")
   
   # Per-chromosome breakdown across all genome
   df_filt2 %>%
     group_by(CHROM) %>%
     summarise(
       SD_freq = sum(AD_REF.LOW, na.rm = TRUE) / sum(DP.LOW, na.rm = TRUE),
       n_sites = n()
     ) %>%
     arrange(CHROM)
   
   # Per-chromosome breakdown of SD allele freq in tolerant bulk across just the QTLs
   all.snps.comb %>%
     group_by(CHROM) %>%
     summarise(
       SD_freq = sum(AD_REF.LOW, na.rm = TRUE) / sum(DP.LOW, na.rm = TRUE),
       n_sites = n()
     ) %>%
     arrange(CHROM)
   
   # Calculate the SD allele frequency in the chromosome regions not in the QTLs
   qlt_chroms <- c("Chr_3", "Chr_8", "Chr_9", "Chr_11", "Chr_12")
   
   df_filt2 %>%
     filter(CHROM %in% qlt_chroms) %>%
     anti_join(all.snps.comb, by = c("CHROM", "POS")) %>%
     group_by(CHROM) %>%
     summarise(
       SD_freq_background = sum(AD_REF.LOW, na.rm = TRUE) / sum(DP.LOW, na.rm = TRUE),
       n_sites = n()
     ) %>%
     arrange(CHROM)
   
  
   # Do chi-square test to see if departures are significant based on the counts
   # Get QTL counts
   qtl_counts <- all.snps.comb %>%
     group_by(CHROM) %>%
     summarise(
       SD_qtl = sum(AD_REF.LOW, na.rm = TRUE),
       total_qtl = sum(DP.LOW, na.rm = TRUE)
     )
   
   # Get non-QTL counts
   background_counts <- df_filt2 %>%
     filter(CHROM %in% qlt_chroms) %>%
     anti_join(all.snps.comb, by = c("CHROM", "POS")) %>%
     group_by(CHROM) %>%
     summarise(
       SD_bg = sum(AD_REF.LOW, na.rm = TRUE),
       total_bg = sum(DP.LOW, na.rm = TRUE)
     )
   
   # Join and run chi-square per chromosome
   qtl_counts %>%
     inner_join(background_counts, by = "CHROM") %>%
     rowwise() %>%
     mutate(
       p_value = chisq.test(matrix(
         c(SD_qtl, total_qtl - SD_qtl,
           SD_bg,  total_bg  - SD_bg),
         nrow = 2
       ))$p.value
     ) %>%
     arrange(CHROM)
   
   # Check is whether the raw deltaSNP distributions differ between QTL sets:
   all.snps.comb %>%
     mutate(QTL_type = ifelse(CHROM %in% c("Chr_8", "Chr_9"), "SH_biased", "SD_biased")) %>%
     group_by(QTL_type) %>%
     summarise(
       mean_deltaSNP = mean(deltaSNP, na.rm = TRUE),
       mean_Gprime = mean(Gprime, na.rm = TRUE),
       n_sites = n()
     )
   
   
###### Now that we see which QTLs are SD or SH biased in the low-pcrit group, lets pool genes for those and ask questions ######
  # We know that Chrom 8 and 9 are the SH QTLs and Chroms 3, 11, and 12 have the SD QTLs
  # So add a new quantifier column
   gene_summary <- gene_summary %>%
     dplyr::mutate(bias = case_when(
       CHROM %in% c("Chr_8", "Chr_9")          ~ "SH_bias",
       CHROM %in% c("Chr_3", "Chr_11", "Chr_12") ~ "SD_bias",
       TRUE                                      ~ NA_character_
     ))

  # dnds difference between two groups? 
   # Run glm with gamma distribution 
   dnds.popbias.mod <- glm(omega ~ bias, family = Gamma(link = "log"), data = subset(gene_summary, omega != 99 & omega != 0))
   summary(dnds.popbias.mod)
   exp(-0.20865)
   
   # Generate boxplot
   dnds.popbias.plot <- ggplot(subset(gene_summary, omega != 99 & omega != 0), aes(x = bias, y = omega, fill = bias))+
     geom_boxplot()+
     scale_y_log10()+
     scale_fill_manual(values = c("sienna2", "seagreen2"))+
     scale_x_discrete(labels = c("SD-bias", "SH-bias"))+
     ylab("dN/dS")+xlab("")+
     theme(legend.position = "none", axis.text.y = element_text(size = 17), axis.text.x = element_text(size = 16),
           axis.title.y = element_text(size=19),
           plot.subtitle = element_text(size =17))+
     labs(subtitle = expression(OR == 0.81 ~ "," ~ italic(p) == 0.048))
   dnds.popbias.plot   
   
   # Export
   ggsave("dnds_by_pop_bias.jpg", dnds.popbias.plot, height = 4, width = 4, dpi = 300)
   
  # More represented in time series? 
   # Run glm with gamma distribution 
   ts.popbias.mod <- glm(In_time_series ~ bias, family = binomial, data = gene_summary)
   summary(ts.popbias.mod)
   
   # Get predicted probabilities with confidence intervals
   pred_data <- data.frame(bias = c("SD_bias", "SH_bias"))
   pred <- predict(ts.popbias.mod, newdata = pred_data, type = "response", se.fit = TRUE)
   pred_data$fit <- pred$fit
   # Get confidence limits
   pred_data$lower <- pred$fit - 1.96 * pred$se.fit
   pred_data$upper <- pred$fit + 1.96 * pred$se.fit
   
  # Make plot
  ts.popbias.mod<- ggplot(pred_data, aes(x = bias, y = fit, color = bias)) +
     geom_point(size = 4) +
     geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 1) +
    ylim(c(0,1))+
     scale_color_manual(values = c("sienna2", "seagreen2"))+
     labs(x = "", y = "Probability of Being a DE gene", color = "Bias") +
     theme(strip.text = element_text(face = "bold")) +
     labs(subtitle = expression(beta == -0.07 ~ "," ~ italic(p) == 0.774))
  ts.popbias.mod
  
   # Export
   ggsave("time_series_by_pop_bias.jpg", ts.popbias.mod, height = 5, width = 5, dpi = 300)
   
 #Repeat summary table code from time series grouping comparison but with pop bias groupings
  # Pivot long to make summary
   ridge_long_popbias <- gene_summary %>%
     dplyr::select(bias, all_of(ridge_vars)) %>%
     pivot_longer(
       cols = all_of(ridge_vars),
       names_to = "metric",
       values_to = "value"
     ) %>%
     filter(!is.infinite(value) & !is.na(value))
   
   # Generate summary
   ridge_summary_popbias <- ridge_long_popbias %>%
     group_by(metric, bias) %>%
     summarise(
       n = sum(!is.na(value)),
       mean_val = mean(value, na.rm = TRUE),
       sd_val = sd(value, na.rm = TRUE),
       se = sd_val / sqrt(n),
       ci_low = mean_val - qt(0.975, df = n - 1) * se,
       ci_high = mean_val + qt(0.975, df = n - 1) * se,
       .groups = "drop"
     )
   
   # Define function to auto call the correct elements when plotting
   get_ridge_ci_pb <- function(metric_name) {
     ridge_summary_popbias %>%
       filter(metric == metric_name)
   }
   
   
  # Upstream genes
   # Run glm with poisson distribution 
   upstream.mod.pb <- glm(n_upstream ~ bias + n_in_gene_per_kb, family = quasipoisson, data = gene_summary)
   summary(upstream.mod.pb)
   exp(0.054018)
   # Make plot
   ridge_upstream.pb <- ggplot(gene_summary, aes(x = n_upstream, y = bias, fill = bias)) +
     geom_density_ridges(alpha = 0.6, scale = 2, color = "white", bandwidth = 2, rel_min_height =0.01,
                         jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0), 
                         point_shape = "|", point_size = 4, point_color = "grey40") +
     geom_errorbar(data = get_ridge_ci_pb("n_upstream"), aes(y = bias, xmin = ci_low,xmax = ci_high), 
                   orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black")+
     geom_point(data = get_ridge_ci_pb("n_upstream"), aes(x = mean_val, y = bias),
                inherit.aes = FALSE,size = 3,shape = 21,fill = "red",color = "black")+
     labs(x = "Number of SNPs 1000 bp upstream", y = "", title = "Distribution of Upstream SNP Counts") +
     scale_fill_manual(values = c("sienna2", "seagreen2"))+
     scale_y_discrete(labels = c("SD_bias" = "SD-biased QTLs", "SH_bias"  = "SH-biased QTLs"))+
     coord_cartesian(xlim = c(0, max(gene_summary$n_upstream))) +
     annotate("text", x = 28, y = 2.5, size = 5,
              label = "OR==1.06 ~ ';' ~ italic(p)==0.528",
              parse = TRUE)+
     # annotate("text", x = 28, y = 3.1, size = 3,
     #          label = "# SNPs upstream ~ AF bias (SD/SH) + # SNPs per kb in gene",
     #          parse = FALSE)+
     theme(legend.position = "none")
   ridge_upstream.pb
   
  # In gene
   # Run glm with poisson distribution 
   ingene.mod.pb <- glm(n_in_gene ~ bias + offset(log(gene_length_kb)), family = quasipoisson, data = gene_summary)
   summary(ingene.mod.pb)
   exp(-0.11770)
   # Make plot
   ridge_in_gene.pb <- ggplot(gene_summary, aes(x = n_in_gene, y = bias, fill = bias)) +
     geom_density_ridges(alpha = 0.6, scale = 2, color = "white", bandwidth = 5, rel_min_height =0.01,
                         jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0), 
                         point_shape = "|", point_size = 4, point_color = "grey40") +
     geom_errorbar(data = get_ridge_ci_pb("n_in_gene"), aes(y = bias, xmin = ci_low,xmax = ci_high),
                   orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black")+
     geom_point(data = get_ridge_ci_pb("n_in_gene"), aes(x = mean_val, y = bias),
                inherit.aes = FALSE,size = 3,shape = 21,fill = "red",color = "black")+
     labs(x = "Number of SNPs in Genes", y = "", title = "Distribution of SNP counts in genes") +
     scale_fill_manual(values = c("sienna2", "seagreen2"))+
     scale_y_discrete(labels = c("SD_bias" = "SD-biased QTLs", "SH_bias"  = "SH-biased QTLs"))+
     coord_cartesian(xlim = c(0, max(gene_summary$n_in_gene))) +
     annotate("text", x = 250, y = 2.5, size = 5,
              label = "OR==0.89 ~ ';' ~ italic(p)==0.017",
              parse = TRUE)+
     # annotate("text", x = 250, y = 3, size = 3,
     #          label = "# SNPs in gene ~ AF bias (SD/SH) + offset by gene length",
     #          parse = FALSE)+
     theme(legend.position = "none")
   ridge_in_gene.pb
   
   
   # In gene per kb
   ingeneperkb.mod.pb <- lm(n_in_gene_per_kb_tf ~ bias, data = gene_summary)
   summary(ingeneperkb.mod.pb)
   # Make plot
   ridge_in_gene_per_kb.pb <- ggplot(gene_summary, aes(x = n_in_gene_per_kb_tf, y = bias, fill = bias)) +
     geom_density_ridges(alpha = 0.6, scale = 2, color = "white", bandwidth = 0.2, rel_min_height =0.01,
                         jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0), 
                         point_shape = "|", point_size = 4, point_color = "grey40") +
     geom_errorbar(data = get_ridge_ci_pb("n_in_gene_per_kb_tf"), aes(y = bias, xmin = ci_low,xmax = ci_high), 
                   orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black")+
     geom_point(data = get_ridge_ci_pb("n_in_gene_per_kb_tf"), aes(x = mean_val, y = bias),
                inherit.aes = FALSE,size = 3,shape = 21,fill = "red",color = "black")+
     labs(x = "Number of SNPs in genes per kb", y = "", title = "Distribution of SNP counts in genes per kb") +
     scale_fill_manual(values = c("sienna2", "seagreen2"))+
     scale_y_discrete(labels = c("SD_bias" = "SD-biased QTLs", "SH_bias"  = "SH-biased QTLs"))+
     coord_cartesian(xlim = c(0, max(gene_summary$n_in_gene_per_kb_tf))) +
     annotate("text", x = 4.4, y = 2.9, size = 5,
              label = "italic(beta)==-0.22 ~ ';' ~ italic(p)==0.003",
              parse = TRUE)+
     # annotate("text", x = 4.4, y = 3.4, size = 3,
     #          label = "# SNPs per kb in gene ~ AF bias (SD/SH)",
     #          parse = FALSE)+
     theme(legend.position = "none")
   ridge_in_gene_per_kb.pb
   
   # Gene length
   # Run glm with poisson distribution 
   length.mod.pb <- glm(gene_length_kb ~ bias, family = Gamma(link = "log"), data = gene_summary)
   summary(length.mod.pb)
   exp(-0.09853)
   # Make plot
   ridge_gene_length.pb <- ggplot(gene_summary, aes(x = gene_length_kb, y = bias, fill = bias)) +
     geom_density_ridges(alpha = 0.6, scale = 2, color = "white", bandwidth = 5, rel_min_height =0.01,
                         jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0), 
                         point_shape = "|", point_size = 4, point_color = "grey40") +
     geom_errorbar(data = get_ridge_ci_pb("gene_length_kb"), aes(y = bias, xmin = ci_low,xmax = ci_high), 
                   orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black")+
     geom_point(data = get_ridge_ci_pb("gene_length_kb"), aes(x = mean_val, y = bias),
                inherit.aes = FALSE,size = 3,shape = 21,fill = "red",color = "black")+
     labs(x = "Gene Length in kb", y = "", title = "Distribution of Gene Lengths") +
     scale_fill_manual(values = c("sienna2", "seagreen2"))+
     scale_y_discrete(labels = c("SD_bias" = "SD-biased QTLs", "SH_bias"  = "SH-biased QTLs"))+
     coord_cartesian(xlim = c(0, max(gene_summary$gene_length_kb))) +
     annotate("text", x = 75, y = 2.5, size = 5,
              label = "OR==0.91 ~ ';' ~ italic(p)==0.311",
              parse = TRUE)+
     # annotate("text", x = 75, y = 3, size = 3,
     #          label = "Gene Length ~ AF bias (SD/SH)",parse = FALSE)+
     theme(legend.position = "none")
   ridge_gene_length.pb
   
   # In exons
   # Run glm with poisson distribution 
   inexon.mod.pb <- glm(n_in_exons ~ bias + offset(log(gene_length_kb)), family = quasipoisson, data = gene_summary)
   summary(inexon.mod.pb)
   exp(-0.004567)
   # Make plot
   ridge_in_exon.pb <- ggplot(gene_summary, aes(x = n_in_exons, y = bias, fill = bias)) +
     geom_density_ridges(alpha = 0.6, scale = 2, color = "white", bandwidth = 2, rel_min_height =0.01,
                         jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0), 
                         point_shape = "|", point_size = 4, point_color = "grey40") +
     geom_errorbar(data = get_ridge_ci_pb("n_in_exons"), aes(y = bias, xmin = ci_low,xmax = ci_high), 
                   orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black")+
     geom_point(data = get_ridge_ci_pb("n_in_exons"), aes(x = mean_val, y = bias),
                inherit.aes = FALSE,size = 3,shape = 21,fill = "red",color = "black")+
     labs(x = "# SNPs in exons", y = "", title = "Distribution of SNP counts in Exons") +
     scale_fill_manual(values = c("sienna2", "seagreen2"))+
     scale_y_discrete(labels = c("SD_bias" = "SD-biased QTLs", "SH_bias"  = "SH-biased QTLs"))+
     coord_cartesian(xlim = c(0, max(gene_summary$n_in_exons))) +
     annotate("text", x = 70, y = 2.5, size = 5,
              label = "OR==1.00 ~ ';' ~ italic(p)==0.966",
              parse = TRUE)+
     # annotate("text", x = 70, y = 3, size = 3,
     #          label = "# SNPs in exons ~ AF bias (SD/SH) + offset by gene length",
     #          parse = FALSE)+
     theme(legend.position = "none")
   ridge_in_exon.pb
   
   
   # In introns
   # Run glm with poisson distribution 
   inintron.mod.pb <- glm(n_introns ~ bias + offset(log(gene_length_kb)), family = quasipoisson, data = gene_summary)
   summary(inintron.mod.pb)
   exp(-0.16900)
   # Make plot
   ridge_in_intron.pb <- ggplot(gene_summary, aes(x = n_introns, y = bias, fill = bias)) +
     geom_density_ridges(alpha = 0.6, scale = 2, color = "white", bandwidth = 4, rel_min_height =0.01,
                         jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0), 
                         point_shape = "|", point_size = 4, point_color = "grey40") +
     geom_errorbar(data = get_ridge_ci_pb("n_introns"), aes(y = bias, xmin = ci_low,xmax = ci_high), 
                   orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black")+
     geom_point(data = get_ridge_ci_pb("n_introns"), aes(x = mean_val, y = bias),
                inherit.aes = FALSE,size = 3,shape = 21,fill = "red",color = "black")+
     labs(x = "# SNPs in introns", y = "", title = "Distribution of SNP counts in Introns") +
     scale_fill_manual(values = c("sienna2", "seagreen2"))+
     scale_y_discrete(labels = c("SD_bias" = "SD-biased QTLs", "SH_bias"  = "SH-biased QTLs"))+
     coord_cartesian(xlim = c(0, max(gene_summary$n_introns))) +
     annotate("text", x = 220, y = 2.5, size = 5,
              label = "OR==0.85 ~ ';' ~ italic(p)==0.006",
              parse = TRUE)+
     # annotate("text", x = 218, y = 3, size = 3,
     #          label = "# of SNPs in introns ~ AF bias (SD/SH) + offset by gene length",
     #          parse = FALSE)+
     theme(legend.position = "none", axis.text.y = element_text(size = 17), axis.text.x = element_text(size = 16),
           axis.title.x = element_text(size=17),
           plot.title = element_text(size =17))
   ridge_in_intron.pb
   
   ridges_patch.pb <- (ridge_in_gene_per_kb.pb / ridge_in_exon.pb / ridge_in_intron.pb / ridge_upstream.pb / ridge_gene_length.pb)
   
   ggsave(filename = "ridgeline_column_plot_pb.png", plot = ridges_patch.pb, width = 8, height = 12, dpi = 300)
   
   
  ## Do genes in one pop bias tend to have higher ratio of SNPs upstream to in exons?
   mod.ratio.1.pb <- lm(log1p(up_to_in_per_kb_ratio) ~ bias, data = gene_summary)
   summary(mod.ratio.1.pb)
   
   
   # Do genes that overlap with one pop bias tend to have more of any of the types of impact categories for the SNPs?
   low.mod.pb <- glm(n_low_impact ~ bias + offset(log(gene_length_kb)), family = quasipoisson, data = gene_summary)
   summary(low.mod.pb)
   exp(-0.01487)
   
   moderat.mod.pb <- glm(n_mod_impact ~ bias + offset(log(gene_length_kb)), family = quasipoisson, data = gene_summary)
   summary(moderat.mod.pb)
   exp(-0.15919)
   
   high.mod.pb <- glm(n_high_impact ~ bias + offset(log(gene_length_kb)), family = quasipoisson, data = gene_summary)
   summary(high.mod.pb)
   exp(-0.6241)
   
   modifier.mod.pb <- glm(n_modifer ~ bias + offset(log(gene_length_kb)), family = quasipoisson, data = gene_summary)
   summary(modifier.mod.pb)
   exp(-0.07998)
   
   # Make plot for the mpdifier impact
   ridge_modifier_impact.pb <- ggplot(gene_summary, aes(x = n_modifer, y = bias, fill = bias)) +
     geom_density_ridges(alpha = 0.6, scale = 2, color = "white", bandwidth = 5, rel_min_height =0.01,
                         jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0), 
                         point_shape = "|", point_size = 4, point_color = "grey40") +
     geom_errorbar(data = get_ridge_ci_pb("n_modifer"), aes(y = bias, xmin = ci_low,xmax = ci_high), 
                   orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black")+
     geom_point(data = get_ridge_ci_pb("n_modifer"), aes(x = mean_val, y = bias),
                inherit.aes = FALSE,size = 3,shape = 21,fill = "red",color = "black")+
     labs(x = "# SNPs with 'modifier' impact", y = "", title = "Distribution of SNP counts with 'modifier' impact") +
     scale_fill_manual(values = c("sienna2", "seagreen2"))+
     scale_y_discrete(labels = c("SD_bias" = "SD-biased QTLs", "SH_bias"  = "SH-biased QTLs"))+
     coord_cartesian(xlim = c(0, max(gene_summary$n_modifer))) +
     annotate("text", x = 220, y = 2.5, size = 5,
              label = "OR==0.92 ~ ';' ~ italic(p)==0.152",
              parse = TRUE)+
     # annotate("text", x = 218, y = 3, size = 3,
     #          label = "# of SNPs with modifier impact ~ AF bias (SD/SH) + offset by gene length",
     #          parse = FALSE)+
     theme(legend.position = "none")
   ridge_modifier_impact.pb
   
   # Make plot for the low impact
   ridge_low_impact.pb <- ggplot(gene_summary, aes(x = n_low_impact, y = bias, fill = bias)) +
     geom_density_ridges(alpha = 0.6, scale = 2, color = "white", bandwidth = 2, rel_min_height =0.01,
                         jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0), 
                         point_shape = "|", point_size = 4, point_color = "grey40") +
     geom_errorbar(data = get_ridge_ci_pb("n_low_impact"), aes(y = bias, xmin = ci_low,xmax = ci_high), 
                   orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black")+
     geom_point(data = get_ridge_ci_pb("n_low_impact"), aes(x = mean_val, y = bias),
                inherit.aes = FALSE,size = 3,shape = 21,fill = "red",color = "black")+
     labs(x = "# SNPs with low impact", y = "", title = "Distribution of SNP counts with low impact") +
     scale_fill_manual(values = c("sienna2", "seagreen2"))+
     scale_y_discrete(labels = c("SD_bias" = "SD-biased QTLs", "SH_bias"  = "SH-biased QTLs"))+
     coord_cartesian(xlim = c(0, max(gene_summary$n_low_impact))) +
     annotate("text", x = 55, y = 2.5, size = 5,
              label = "OR==0.99 ~ ';' ~ italic(p)==0.899",
              parse = TRUE)+
     # annotate("text", x = 55, y = 3, size = 3,
     #          label = "# of SNPs with low impact ~ AF bias (SD/SH) + offset by gene length",
     #          parse = FALSE)+
     theme(legend.position = "none")
   ridge_low_impact.pb
   
   # Make plot for the moderate impact
   ridge_moderate_impact.pb <- ggplot(gene_summary, aes(x = n_mod_impact, y = bias, fill = bias)) +
     geom_density_ridges(alpha = 0.6, scale = 2, color = "white", bandwidth = 2, rel_min_height =0.01,
                         jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0), 
                         point_shape = "|", point_size = 4, point_color = "grey40") +
     geom_errorbar(data = get_ridge_ci_pb("n_mod_impact"), aes(y = bias, xmin = ci_low,xmax = ci_high), 
                   orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black")+
     geom_point(data = get_ridge_ci_pb("n_mod_impact"), aes(x = mean_val, y = bias),
                inherit.aes = FALSE,size = 3,shape = 21,fill = "red",color = "black")+
     labs(x = "# SNPs with moderate impact", y = "", title = "Distribution of SNP counts with moderate impact") +
     scale_fill_manual(values = c("sienna2", "seagreen2"))+
     scale_y_discrete(labels = c("SD_bias" = "SD-biased QTLs", "SH_bias"  = "SH-biased QTLs"))+
     coord_cartesian(xlim = c(0, max(gene_summary$n_mod_impact))) +
     annotate("text", x = 24, y = 2.5, size = 5,
              label = "OR==0.85 ~ ';' ~ italic(p)==0.355",
              parse = TRUE)+
     # annotate("text", x = 24, y = 3, size = 3,
     #          label = "# of SNPs with moderate impact ~ AF bias (SD/SH) + offset by gene length",
     #          parse = FALSE)+
     theme(legend.position = "none")
   ridge_moderate_impact.pb
   
   # Make plot for the high impact
   ridge_high_impact.pb <- ggplot(gene_summary, aes(x = n_high_impact, y = bias, fill = bias)) +
     geom_density_ridges(alpha = 0.6, scale = 2, color = "white", bandwidth = 0.05, rel_min_height =0.01,
                         jittered_points = TRUE, position = position_points_jitter(width = 0.1, height = 0), 
                         point_shape = "|", point_size = 4, point_color = "grey40") +
     geom_errorbar(data = get_ridge_ci_pb("n_high_impact"), aes(y = bias, xmin = ci_low,xmax = ci_high), 
                   orientation = "y", inherit.aes = FALSE, height = 0.15, linewidth = 1.1, color = "black")+
     geom_point(data = get_ridge_ci_pb("n_high_impact"), aes(x = mean_val, y = bias),
                inherit.aes = FALSE,size = 3,shape = 21,fill = "red",color = "black")+
     labs(x = "# SNPs with high impact", y = "", title = "Distribution of SNP counts with high impact") +
     scale_fill_manual(values = c("sienna2", "seagreen2"))+
     scale_y_discrete(labels = c("SD_bias" = "SD-biased QTLs", "SH_bias"  = "SH-biased QTLs"))+
     coord_cartesian(xlim = c(0, max(gene_summary$n_high_impact))) +
     annotate("text", x = 1.4, y = 2.5, size = 5,
              label = "OR==0.54 ~ ';' ~ italic(p)==0.41",
              parse = TRUE)+
     # annotate("text", x = 1.4, y = 3, size = 3,
     #          label = "# of SNPs with high impact ~ AF bias (SD/SH) + offset by gene length",
     #          parse = FALSE)+
     theme(legend.position = "none")
   ridge_high_impact.pb
   
   ridges_patch_impact.pb <- (ridge_modifier_impact.pb / ridge_low_impact.pb / ridge_moderate_impact.pb / ridge_high_impact.pb)
   
   ggsave(filename = "ridgeline_column_plot_impacts.pb.png", plot = ridges_patch_impact.pb, width = 7, height = 10, dpi = 300)
   
   
   
   
   
   
   # G' comparison between time series groupings
   gprime.mod.pb <- glm(max_Gprime ~ bias, family = Gamma(link = "log"), data = gene_summary)
   summary(gprime.mod.pb)
   exp(-0.088774)
   
   # Generate boxplot
   gprime.plot.pb <- ggplot(gene_summary, aes(x = bias, y = max_Gprime, fill = bias))+
     geom_boxplot()+
     scale_y_log10()+
     scale_fill_manual(values = c("sienna2", "seagreen2"))+
     scale_x_discrete(labels = c("SD-bias", "SH-bias"))+
     ylab("Maximum G' per gene")+xlab("")+
     theme(legend.position = "none", axis.text.y = element_text(size = 17), axis.text.x = element_text(size = 16),
           axis.title.y = element_text(size=19),
           plot.subtitle = element_text(size =17))+
     labs(subtitle = expression(OR == 0.92 ~ "," ~ italic(p) < 0.001))
   gprime.plot.pb   
   
   ggsave("Gprime_by_popbias.jpg", gprime.plot.pb, height = 4, width = 4, dpi = 300)
   
   
   
# Different GO processes? 
   
   # Make reduced GeneUniverse that is just our significant genes
   geneUniverse_just_sigs <- geneUniverse[geneUniverse %in% gene_summary$Gene_ID]
   
  # SD-biased QTL genes 
   #Get first column from  results w/ gene names and pvalues. Send to list. 
   GeneInt.SD <- gene_summary %>%
     filter(bias == "SD_bias") %>%
     pull(Gene_ID) 
   
   
   #Make named vector with the gene universe where genes of interest are coded with a 1 so the 'new' function knows to focus on them. 
   geneList.SD <- factor(as.integer(geneUniverse_just_sigs %in% GeneInt.SD)) 
   names(geneList.SD) <- geneUniverse_just_sigs
   
   # We now have all data necessary to build an object of type topGOdata. This object will contain all gene
   # identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information
   # needed to perform the desired enrichment analysis.
   SD_GO <- new("topGOdata", description = "GO analysis of all sig genes", ontology = "BP", nodeSize = 10,
                 allGenes = geneList.SD, annot = annFUN.gene2GO, gene2GO = geneID2GO)
   
   # Once we have an object of class topGOdata we can start with the enrichment analysis.
   # Run the classic fisher exact test to find enriched go terms
   resultFisher.SD <- runTest(SD_GO, algorithm = "weight01", statistic = "fisher")
   
   resultFisher.SD #View results summary
   
   # GenTable is an easy to use function for analysing the most significant GO terms and the corresponding p
   # values. In the following example, we list the top 10 significant GO terms identified by the elim method. At
   # the same time we also compare the ranks and the p-values of these GO terms with the ones obtained by the
   # classic method   
   #Adjusted to look just at the classic fisher results since we have just counts
   allRes.SD <- GenTable(SD_GO, fisher = resultFisher.SD, ranksOf = "fisher", topNodes = 150, numChar = 500)
   colnames(allRes.SD)[6] <- "p-value"
   allRes.SD$`p-value` <- as.numeric(allRes.SD$`p-value`)
   
   # Filter nodes with pvalue greater than 0.05 
   allRes.SD <- allRes.SD[allRes.SD$`p-value` < 0.05,]
   
   
   #Extract names of significant genes in GenTable result. Add to column in allRes.SD
   allRes.SD$genes <- sapply(allRes.SD$GO.ID, function(x)
   {
     genes<-genesInTerm(SD_GO, x)
     genes[[1]][genes[[1]] %in% GeneInt.SD] # myGenes is the queried gene list
   })
   
   # Create and export enrichment result as a formatted and interactable table
   SD.sigs.table <- reactable(allRes.SD[,1:6], defaultPageSize = 150, theme = pff(centered = FALSE, font_color = "black"), 
                                  wrap = FALSE, bordered = TRUE, compact = TRUE, striped = TRUE, highlight = TRUE, fullWidth = TRUE,
                                  columns = list(Term = colDef(minWidth = 280),
                                                 GO.ID = colDef(cell = pill_buttons(colors = "darkgreen"), minWidth = 120)
                                  )
   ) %>%
     add_title("SD-biased QTL genes", font_size = 20) %>% 
     add_subtitle("290 of 443 genes post filtering", font_size = 16, font_style = "italic") 
   
   
   SD.sigs.table #View table
   
   #Save table as an html file
   save_reactable_test(SD.sigs.table, "SD-biased QTL genes.html")
   
   # investigate how the significant GO terms are distributed over the GO graph
   jpeg(filename = "TopGO subgraph SD genes.jpg", width = 10, height = 10, units = "in", res = 900)
   par(cex = 0.9)
   showSigOfNodes(SD_GO, score(resultFisher.SD), firstSigNodes = 25, useInfo = 'all')
   dev.off()   
   
   
   # Make sure you load the T.californicus OrgDb package before doing this. 
   # SD-biased QTL genes
   #Create sim matrix using method Wang since it is most recent and keytype GID since we dont have ENTREZIDs
   simMatrix.SD <- calculateSimMatrix(allRes.SD$GO.ID,
                                       orgdb = Tcalif_orgdb_object,
                                       ont="BP",
                                       method="Wang", 
                                       keytype = "GID")
   
   #Create groupings of reduced terms for easier visualization
   reducedTerms.SD <- reduceSimMatrix(simMatrix.SD,
                                       threshold=0.8,
                                       orgdb=Tcalif_orgdb_object,
                                       keytype = "GID")
   
   # Make scatter plot depicting groups and distance between terms  
   jpeg(filename = "rrvgo scatter SD.jpg", width = 16, height = 11, units = "in", res = 300)
   scatterPlot(simMatrix.SD, reducedTerms.SD)
   dev.off()
   
   # Make treemap of GO terms clustered under their parent terms
   jpeg(filename = "rrvgo treemap SD.jpg", width = 8, height = 5, units = "in", res = 300)
   treemapPlot(reducedTerms.SD, fontsize.labels = c(22, 22),  # c(group label, item label)
               fontsize.title = 20, title = "SD-biased QTL genes")
   dev.off()
   
   # Query GO terms for gene IDs in that term
   # GO:0007041 lysosomal transport
   int.list <- allRes.SD$genes$`GO:0120035`
   # Pull out parts of the combined_sigs master results data frame that correspond to go terms of interest saved from above
   gene_summary[gene_summary$Gene_ID %in% int.list, c(2,3)]
   
   
# SH-biased QTL genes 
   #Get first column from  results w/ gene names and pvalues. Send to list. 
   GeneInt.SH <- gene_summary %>%
     filter(bias == "SH_bias") %>%
     pull(Gene_ID) 
   
   #Make named vector with the gene universe where genes of interest are coded with a 1 so the 'new' function knows to focus on them. 
   geneList.SH <- factor(as.integer(geneUniverse_just_sigs %in% GeneInt.SH)) 
   names(geneList.SH) <- geneUniverse_just_sigs
   
   # We now have all data necessary to build an object of type topGOdata. This object will contain all gene
   # identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information
   # needed to perform the desired enrichment analysis.
   SH_GO <- new("topGOdata", description = "GO analysis of all sig genes", ontology = "BP", nodeSize = 10,
                allGenes = geneList.SH, annot = annFUN.gene2GO, gene2GO = geneID2GO)
   
   # Once we have an object of class topGOdata we can start with the enrichment analysis.
   # Run the classic fisher exact test to find enriched go terms
   resultFisher.SH <- runTest(SH_GO, algorithm = "weight01", statistic = "fisher")
   
   resultFisher.SH #View results summary
   
   # GenTable is an easy to use function for analysing the most significant GO terms and the corresponding p
   # values. In the following example, we list the top 10 significant GO terms identified by the elim method. At
   # the same time we also compare the ranks and the p-values of these GO terms with the ones obtained by the
   # classic method   
   #Adjusted to look just at the classic fisher results since we have just counts
   allRes.SH <- GenTable(SH_GO, fisher = resultFisher.SH, ranksOf = "fisher", topNodes = 150, numChar = 500)
   colnames(allRes.SH)[6] <- "p-value"
   allRes.SH$`p-value` <- as.numeric(allRes.SH$`p-value`)
   
   # Filter nodes with pvalue greater than 0.05 
   allRes.SH <- allRes.SH[allRes.SH$`p-value` < 0.05,]
   
   
   #Extract names of significant genes in GenTable result. Add to column in allRes.SH
   allRes.SH$genes <- sapply(allRes.SH$GO.ID, function(x)
   {
     genes<-genesInTerm(SH_GO, x)
     genes[[1]][genes[[1]] %in% GeneInt.SH] # myGenes is the queried gene list
   })
   
   # Create and export enrichment result as a formatted and interactable table
   SH.sigs.table <- reactable(allRes.SH[,1:6], defaultPageSize = 150, theme = pff(centered = FALSE, font_color = "black"), 
                              wrap = FALSE, bordered = TRUE, compact = TRUE, striped = TRUE, highlight = TRUE, fullWidth = TRUE,
                              columns = list(Term = colDef(minWidth = 280),
                                             GO.ID = colDef(cell = pill_buttons(colors = "darkgreen"), minWidth = 120)
                              )
   ) %>%
     add_title("SH-biased QTL genes", font_size = 20) %>% 
     add_subtitle("153 of 443 genes post filtering", font_size = 16, font_style = "italic") 
   
   
   SH.sigs.table #View table
   
   #Save table as an html file
   save_reactable_test(SH.sigs.table, "SH-biased QTL genes.html")
   
   # investigate how the significant GO terms are distributed over the GO graph
   jpeg(filename = "TopGO subgraph SH genes.jpg", width = 10, height = 10, units = "in", res = 900)
   par(cex = 0.9)
   showSigOfNodes(SH_GO, score(resultFisher.SH), firstSigNodes = 25, useInfo = 'all')
   dev.off()   
   
   
   # Make sure you load the T.californicus OrgDb package before doing this. 
   # SH-biased QTL genes
   #Create sim matrix using method Wang since it is most recent and keytype GID since we dont have ENTREZIDs
   simMatrix.SH <- calculateSimMatrix(allRes.SH$GO.ID,
                                      orgdb = Tcalif_orgdb_object,
                                      ont="BP",
                                      method="Wang", 
                                      keytype = "GID")
   
   #Create groupings of reduced terms for easier visualization
   reducedTerms.SH <- reduceSimMatrix(simMatrix.SH,
                                      threshold=0.8,
                                      orgdb=Tcalif_orgdb_object,
                                      keytype = "GID")
   
   # Make scatter plot depicting groups and distance between terms  
   jpeg(filename = "rrvgo scatter SH.jpg", width = 16, height = 11, units = "in", res = 300)
   scatterPlot(simMatrix.SH, reducedTerms.SH)
   dev.off()
   
   # Make treemap of GO terms clustered under their parent terms
   jpeg(filename = "rrvgo treemap SH.jpg", width = 8, height = 5, units = "in", res = 300)
   treemapPlot(reducedTerms.SH, fontsize.labels = c(22, 22),  # c(group label, item label)
               fontsize.title = 20, title = "SH-biased QTL genes")
   dev.off()
   
   
   # Query GO terms for gene IDs in that term
   # GO:0007041 lysosomal transport
   int.list <- allRes.SH$genes$`GO:0015980`
   # Pull out parts of the combined_sigs master results data frame that correspond to go terms of interest saved from above
   gene_summary[gene_summary$Gene_ID %in% int.list, c(2,3)]
   
   
   # Save gene_sumamry to Excel file
   write_xlsx(gene_summary, path = "QTL figures and results/Table S4 Significant genes.xlsx")
   
   
   
   # Now that we have the gene summary done, add its data to our functional group data frames and export those too
   mitos <- merge(mitos, gene_summary[,-c(3)], by = "Gene_ID", all.x = TRUE)
   glycolysis <- merge(glycolysis, gene_summary[,-c(3)], by = "Gene_ID", all.x = TRUE)
   carotenoids <- merge(carotenoids, gene_summary[,-c(3)], by = "Gene_ID", all.x = TRUE)
   cuticles <- merge(cuticles, gene_summary[,-c(3)], by = "Gene_ID", all.x = TRUE)
   antioxidants <- merge(antioxidants, gene_summary[,-c(3)], by = "Gene_ID", all.x = TRUE)
   
   # Export significant gene lists as xlsx file with tabs along with overlapping genes from previous study
   functional.groups.wb <- createWorkbook()
   
   addWorksheet(functional.groups.wb, "Cuticle genes")
   writeData(functional.groups.wb, "Cuticle genes", cuticles, rowNames = F)
   
   addWorksheet(functional.groups.wb, "Glycolysis genes")
   writeData(functional.groups.wb, "Glycolysis genes", glycolysis, rowNames = F)
   
   addWorksheet(functional.groups.wb, "Antioxidant genes")
   writeData(functional.groups.wb, "Antioxidant genes", antioxidants, rowNames = F)
   
   addWorksheet(functional.groups.wb, "Carotenoid genes")
   writeData(functional.groups.wb, "Carotenoid genes", carotenoids, rowNames = F)
   
   addWorksheet(functional.groups.wb, "Mitochondrial genes")
   writeData(functional.groups.wb, "Mitochondrial genes", mitos, rowNames = F)
   
   saveWorkbook(functional.groups.wb, "Functional groups results.xlsx", overwrite = TRUE)
   
   
   
   
###### Now lets make a figure with genes of interest ######
   
  # Start with LCC2 for fun
    # Change ID for any gene of interest
   target_gene <- gff_gr[grep("TCAL_08505", mcols(gff_gr)$ID)]
   target_gene
   
   # Define a window around the target gene (adjustable)
    # Use position values from target_gene object 
   window_start <- 13366328 - 1000
   window_end   <- 13367932 + 1000
   
   # Get all features in that window on the chromosome of interest and save as GR object
   region_gr <- gff_gr[
     seqnames(gff_gr) == "Chr_11" &
       start(gff_gr) >= window_start &
       end(gff_gr) <= window_end
   ]
   
   # How many genes are in this window?
   region_genes <- region_gr[mcols(region_gr)$type == "gene"]
   mcols(region_genes)$ID
   
   # Get the G' and SNP data for this window
   window_snps <- df_filt2[
     df_filt2$CHROM == "Chr_11" &
       df_filt2$POS >= window_start &
       df_filt2$POS <= window_end,
   ]
   
   # Export as example data
   write.csv(window_snps, "window_snps.csv", row.names = FALSE)
   
   # Check data frame components
   nrow(window_snps)
   head(window_snps)
   
   
   ## --- SNP track ---
   p_snps <- ggplot(window_snps, aes(x = POS, y = deltaSNP, color = deltaSNP, alpha = ifelse(qvalue < 0.05, 1, 0.3))) +
     geom_point(size = 1) +
     scale_color_gradient2(low = "blue", mid = "grey80", high = "red", midpoint = 0, name = "ΔSNPindex") +
     scale_alpha_identity() +
     scale_x_continuous(limits = c(window_start, window_end), labels = scales::comma) +
     geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
     labs(y = "ΔSNPindex", x = NULL) +
     theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
   
   ## --- 3. Gene model track ---
   # Extract CDS features for the 5 genes in the window
   region_ids <- mcols(region_genes)$ID
   # Get mRNA IDs (Parent of CDS is mRNA, which ends in -PA)
   region_mrna_ids <- paste0(region_ids, "-PA")
   
   # Subset region_gr to keep only features that are CDS and match gene ID(s) in window and save as GR object
   cds_gr <- region_gr[mcols(region_gr)$type == "CDS" & as.character(mcols(region_gr)$Parent) %in% region_mrna_ids]
   
   # Make data frame for plotting
   cds_df <- data.frame(
     molecule  = as.character(mcols(cds_gr)$Parent),
     gene      = sub("-PA", "", as.character(mcols(cds_gr)$Parent)),
     start     = start(cds_gr),
     end       = end(cds_gr),
     forward   = as.character(strand(cds_gr)) == "+",
     highlight = sub("-PA", "", as.character(mcols(cds_gr)$Parent)) == "TCAL_08505"
   )
   
  # Export as example data
   write.csv(cds_df, "cds_df.csv", row.names = FALSE)
   
   # Make gggenes plot
   p_genes <- ggplot(cds_df, aes(xmin = start, xmax = end, y = gene, fill = highlight, forward = forward)) +
     geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width  = unit(2, "mm")) +
     scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "grey70"), guide = "none") +
     scale_x_continuous(limits = c(window_start, window_end), labels = scales::comma) +
     labs(x = "Position (Chr_11)", y = NULL) +
     theme_genes() +
     theme(axis.text.y = element_text(size = 8))
   
   ## --- Assemble with patchwork ---
  gene_of_interest<- p_snps / p_genes +
     plot_layout(heights = c(1.5, 1)) +
     plot_annotation(title = "TCAL_08505 (cyp6a14) — Chr_11",
                     subtitle = "+/- 1 kb",
                     theme = theme(plot.title = element_text(size = 12, face = "bold"), 
                                   plot.subtitle = element_text(size = 10)
                                   ))
  
  ggsave("Gene Tracks/cyp6a14 only gene track with delta snp.jpg", gene_of_interest, width = 7, height = 3, dpi = 300)
   
###### GO enrichment of top 10 significant annotated genes per QTL #####
  
# Get the subset of top 10 most significant genes per qtl from gene_summary
  top10_genes_per_QTL <- gene_summary %>%
    filter(!grepl("Protein of unknown function", Description, ignore.case = TRUE)) %>%
    group_by(QTL) %>%
    arrange(desc(max_Gprime), .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup()

# Top 10 QTL genes 
  #Get first column from results w/ gene names and pvalues. Send to list. 
  GeneInt.top10 <- top10_genes_per_QTL %>%
    pull(Gene_ID) 
  
  
#Make named vector with the gene universe where genes of interest are coded with a 1 so the 'new' function knows to focus on them. 
  geneList.top10 <- factor(as.integer(geneUniverse_just_sigs %in% GeneInt.top10)) 
  names(geneList.top10) <- geneUniverse_just_sigs
  
# We now have all data necessary to build an object of type topGOdata. This object will contain all gene
# identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information
# needed to perform the desired enrichment analysis.
  top10_GO <- new("topGOdata", description = "GO analysis of top 10 genes", ontology = "BP", nodeSize = 10,
               allGenes = geneList.top10, annot = annFUN.gene2GO, gene2GO = geneID2GO)
  
  # Once we have an object of class topGOdata we can start with the enrichment analysis.
  # Run the classic fisher exact test to find enriched go terms
  resultFisher.top10 <- runTest(top10_GO, algorithm = "weight01", statistic = "fisher")
  
  resultFisher.top10 #View results summary
  
  # GenTable is an easy to use function for analysing the most significant GO terms and the corresponding p
  # values. In the following example, we list the top 10 significant GO terms identified by the elim method. At
  # the same time we also compare the ranks and the p-values of these GO terms with the ones obtained by the
  # classic method   
  #Adjusted to look just at the classic fisher results since we have just counts
  allRes.top10 <- GenTable(top10_GO, fisher = resultFisher.top10, ranksOf = "fisher", topNodes = 150, numChar = 500)
  colnames(allRes.top10)[6] <- "p-value"
  allRes.top10$`p-value` <- as.numeric(allRes.top10$`p-value`)
  
  # Filter nodes with pvalue greater than 0.05 
  allRes.top10 <- allRes.top10[allRes.top10$`p-value` < 0.05,]
  
  
  #Extract names of significant genes in GenTable result. Add to column in allRes.top10
  allRes.top10$genes <- sapply(allRes.top10$GO.ID, function(x)
  {
    genes<-genesInTerm(top10_GO, x)
    genes[[1]][genes[[1]] %in% GeneInt.top10] # myGenes is the queried gene list
  })
  
  # Create and export enrichment result as a formatted and interactable table
  top10.sigs.table <- reactable(allRes.top10[,1:6], defaultPageSize = 150, theme = pff(centered = FALSE, font_color = "black"), 
                             wrap = FALSE, bordered = TRUE, compact = TRUE, striped = TRUE, highlight = TRUE, fullWidth = TRUE,
                             columns = list(Term = colDef(minWidth = 280),
                                            GO.ID = colDef(cell = pill_buttons(colors = "darkgreen"), minWidth = 120)
                             )
  ) %>%
    add_title("top10-biased QTL genes", font_size = 20) %>% 
    add_subtitle("60 of 443 genes post filtering", font_size = 16, font_style = "italic") 
  
  
  top10.sigs.table #View table
  
  #Save table as an html file
  save_reactable_test(top10.sigs.table, "top10-biased QTL genes.html")
  
  # investigate how the significant GO terms are distributed over the GO graph
  jpeg(filename = "TopGO subgraph top10 genes.jpg", width = 10, height = 10, units = "in", res = 900)
  par(cex = 0.9)
  showSigOfNodes(top10_GO, score(resultFisher.top10), firstSigNodes = 25, useInfo = 'all')
  dev.off()   
  
  
  # Make sure you load the T.californicus OrgDb package before doing this. 
  # top10 QTL genes
  #Create sim matrix using method Wang since it is most recent and keytype GID since we dont have ENTREZIDs
  simMatrix.top10 <- calculateSimMatrix(allRes.top10$GO.ID,
                                     orgdb = Tcalif_orgdb_object,
                                     ont="BP",
                                     method="Wang", 
                                     keytype = "GID")
  
  #Create groupings of reduced terms for easier visualization
  reducedTerms.top10 <- reduceSimMatrix(simMatrix.top10,
                                     threshold=0.8,
                                     orgdb=Tcalif_orgdb_object,
                                     keytype = "GID")
  
  # Make scatter plot depicting groups and distance between terms  
  jpeg(filename = "rrvgo scatter top10.jpg", width = 16, height = 11, units = "in", res = 300)
  scatterPlot(simMatrix.top10, reducedTerms.top10)
  dev.off()
  
  # Make treemap of GO terms clustered under their parent terms
  jpeg(filename = "rrvgo treemap top10.jpg", width = 9, height = 6, units = "in", res = 300)
  treemapPlot(reducedTerms.top10, fontsize.labels = c(22, 22),  # c(group label, item label)
              fontsize.title = 20, title = "Top 10 QTL genes by G'")
  dev.off()
  
  # Query GO terms for gene IDs in that term
  # GO:0007041 lysosomal transport
  int.list <- allRes.top10$genes$`GO:0006950`
  # Pull out parts of the combined_sigs master results data frame that correspond to go terms of interest saved from above
  gene_summary[gene_summary$Gene_ID %in% int.list, c(1,2,3,4)]
    
  