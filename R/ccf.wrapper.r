# This wrapper is designed to work with any single sample data
# The two main requirements are:
# a mutation table (in the format described below) which was originally obtained from https://confluence.broadinstitute.org/display/GDAC/DCC+MAFs
# but has been modified to fit description below
# ASCAT copy number data (a segmented matrix)
# source("~/Work/GD.functions.R")
# source("~/Work/STM_final_scripts/main.functions.r")

suppressPackageStartupMessages(library(gdata))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(sequenza))

#' clonality.estimation
#' clonality.estimation integrates the read depth and copy number information at single nucleotide variant locations to determine the timing and clonality of each mutation.
#' @usage clonality.estimation(mutation.table.loc, seg.mat.loc, data.type, TCGA.barcode, ANALYSIS.DIR, sub.clonal.cut.off = 1, min.var.prop = 0.05, min.alt.reads = 5, min.depth = 30, plotting = TRUE)
#' @param mutation.table.loc The location of a mutation table. 
#'    The table should contain the following columns:
#'                      c('Patient'
#'                      ,'Chr'
#'                      ,'Start_position'
#'                      ,'End_position'
#'                      ,'Reference'
#'                      ,'Alternate'
#'                      ,'Variant_freq'
#'                      ,'Ref_freq'
#'                      ,'Hugo_Symbol'
#'                      ,'Variant_Classification')
#' @param seg.mat.loc The location of a segmented matrix (in RData format). This should be in the format of a list. 
#' @param data.type Data type, e.g. TCGA_BRCA
#' @param TCGA.barcode Sample identifier, e.g. 'TCGA-AN-A0FF'
#' @param ANALYSIS.DIR Where do you want the folders to be saved?
#' @param sub.clonal.cut.off What cut off should be used for the subclonal fraction, defaults to 1 (i.e. 95\% CI overlaps 1)
#' @param min.var.prop Minimum variant proportion (defaults to 0.05)
#' @param min.alt.reads Minimum number of alternate reads (defaults to 5)
#' @param min.depth Minimum depth at mutated base (defaults to 30)
#' @param plotting Do you want plotting? (defaults to TRUE)
#' 
#' @details
#' 
#' Nothing is returned directly. However, three files (if plotting TRUE) should be returned within a folder in the specified directory. Patient.earlylate.tsv is a tab-separated table that contains a number of columns. Description of the columns is as follows:
#' [1] "patient"  #Patient identifier;
#' [2] "TCGA.purity"  # ASCAT tumor purity estimate;
#' [3] "mutation_id" # mutation id;
#' [4] "Reference_Base" # reference base; 
#' [5] "Alternate_Base" # alternate base;
#' [6] "ref_counts" # number of reads supporting reference base;
#' [7] "var_counts" # number of reads supporting variant base;
#' [8] "obs.VAF" # observed variant allele frequency;
#' [9] "normal_cn" # local normal copy number at mutation site;
#' [10] "minor_cn" # local minor allele copy number at mutation site;
#'[11] "major_cn" # local major allele copy number at mutation site;
#' [12] "mut.multi" # mutation multiplicity (number of copies of mutant allele);
#'[13] "absolute.ccf" # estimated cancer cell fraction;
#'[14] "absolute.ccf.0.05" # CI of ccf;
#'[15] "absolute.ccf.0.95" # CI of ccf;
#'[16] "prob.clonal";      
#'[17] "prob.subclonal";    
#'[18] "GD.status";        
#'[19] "timing" # timing based on copy number;
#'[20] "comb.timing" # timing based on clonality and copy number;
#' @examples 
#' 
#' clonality.estimation(mutation.table.loc="BRCA.mutation.table.txt"
#'                     ,seg.mat.loc="tcga.brca.seg.hg19.rdata"
#'                     ,data.type='TCGA_BRCA'
#'                     ,TCGA.barcode="TCGA-AR-A2LR"
#'                     ,ANALYSIS.DIR="example/")
#' 
#' @export


clonality.estimation <- function(mutation.table.loc,
                                 seg.mat.loc,
                                 data.type,
                                 TCGA.barcode,
                                 ANALYSIS.DIR,
                                 sub.clonal.cut.off = 1,
                                 min.var.prop = 0.05,
                                 min.alt.reads = 5,
                                 min.depth = 30,
                                 plotting = TRUE) {
  # The following function is the main wrapper for clonality estimation
  # The input required is as follows
  # mutation.table.loc # mutation tables as in ../input
  # seg.mat.loc # segmented matrices in list format as in ../input
  # data type # e.g. TCGA_BRCA
  # ANALYSIS.DIR #what directory do you want to create as a base directory, e.g. example/


  # create the root directory
  data.folder <- paste(ANALYSIS.DIR, data.type, sep = "")

  if (!file.exists(data.folder)) {
    if (!dir.create(data.folder, recursive = TRUE)) {
      stop("Unable to create root directory.\n")
    } 
  }

  # create a folder for the specific patient
  patient.folder <- paste(data.folder, "/", TCGA.barcode, sep = "")

  if (!file.exists(patient.folder)) {
    if (!dir.create(patient.folder, recursive = TRUE)) {
      stop("Unable to create tmp directory.\n")
    }
  }


  # load somatic mutations ####
  if (!file.exists(mutation.table.loc)) {
    stop(paste("\nUnable to locate:\n", mutation.table.loc, " .\n"))
  }

  mutation.table <- read.table(mutation.table.loc,
    header = TRUE,
    stringsAsFactors = FALSE,
    sep = "\t",
    fill = TRUE,
    quote = '"'
  )

  req.mut.colnames <- c(
    "Patient",
    "Chr",
    "Start_position",
    "End_position",
    "Reference",
    "Alternate",
    "Variant_freq",
    "Ref_freq",
    "Hugo_Symbol",
    "Variant_Classification"
  )

  if (!all(req.mut.colnames %in% colnames(mutation.table))) {
    stop(paste("Mutation table not in correct format\nThe following columns are required:\n",
      PasteVector(req.mut.colnames, sep = "\n"),
      sep = ""
    ))
  }

  # load somatic copy numbers ####
  if (!file.exists(seg.mat.loc)) {
    stop(paste("\nUnable to locate:\n", seg.mat.loc, " .\n"))
  }

  seg.mat.copy <- load(seg.mat.loc)
  seg.mat.copy.list <- get(seg.mat.copy)
  seg.mat.copy <- seg.mat.copy.list$segments

  # combine copy number and mutation data #####
  # first of all, select only samples where both available
  barcodes <- intersect(unique(seg.mat.copy$SampleID), unique(mutation.table[, 1]))
  cat(paste("\n\nThere are ", length(barcodes), " patients with both copy number and mutation data. \nOnly these will be used "))

  if (!TCGA.barcode %in% barcodes) {
    stop("TCGA barcode not found in either mutation table or seg.mat.copy")
  }


  # remove sex chromosomes
  seg.mat.copy <- seg.mat.copy[seg.mat.copy$SampleID %in% barcodes, ]
  seg.mat.copy <- seg.mat.copy[seg.mat.copy$Chr %in% c(1:22), ]
  mutation.table <- mutation.table[mutation.table[, 1] %in% barcodes, ]
  mutation.table <- mutation.table[mutation.table$Chr %in% c(1:22), ]
  mutation.table$Chr <- as.numeric(mutation.table$Chr)

  # select patient specific data
  sub.mat.copy <- seg.mat.copy[seg.mat.copy$SampleID == TCGA.barcode, , drop = FALSE]
  sub.mat.mut <- mutation.table[mutation.table[, 1] == TCGA.barcode, , drop = FALSE]

  # check whether there is only one copy number profile
  if (length(unique(as.character(sub.mat.copy$ChipNames))) > 1) {
    stop("You have more than one copy number profile")
  }

  if (length(unique(as.character(sub.mat.copy$NormalFile))) > 1) {
    stop("You have more than one normal copy number profile")
  }


  sub.mat.copy <- unique(sub.mat.copy)

  # combine mutation and copy number
  mut.table <- data.frame(t(sapply(1:nrow(sub.mat.mut), identify.mut.copy.number.ascat, sub.mat.mut, sub.mat.copy)),
    stringsAsFactors = FALSE
  )
  mut.table <- mut.table[!is.na(mut.table$minor_cn), ]
  mut.table <- mut.table[!is.na(mut.table$ref_counts), ]
  mut.table <- mut.table[!duplicated(mut.table$mutation_id), ]

  ### APPLY filters ####

  # min.alt.reads.cut.off
  mut.table <- mut.table[as.numeric(mut.table$var_counts) >= min.alt.reads, ]

  # min.cov.at.site
  mut.table <- mut.table[as.numeric(mut.table$var_counts) + as.numeric(mut.table$ref_counts) >= min.depth, ]

  # min variant proportion
  mut.table <- mut.table[(as.numeric(mut.table$var_counts) / (as.numeric(mut.table$var_counts) + as.numeric(mut.table$ref_counts))) >= min.var.prop, ]


  # calculating genome doubling estimates
  sub.mat.copy.arm <- get.seg.mat.arm(sub.mat.copy)
  # seg.mat.minor    - segmented matrix with minor allele (e.g. seg.mat.LOH)
  #                  - columns should be as follows: c("Sample","Chrom","Start","End","Num.probes","val")
  # seg.mat.copy     - segmented total copy number matrix (columns as above)
  sub.mat.minor <- cbind(
    sub.mat.copy.arm$SampleID,
    sub.mat.copy.arm$Chr,
    sub.mat.copy.arm$Start,
    sub.mat.copy.arm$End,
    sub.mat.copy.arm$nProbes,
    apply(cbind(sub.mat.copy.arm$nA, sub.mat.copy.arm$nB), 1, min)
  )
  colnames(sub.mat.minor) <- c("Sample", "Chrom", "Start", "End", "Num.probes", "val")

  sub.mat.cn <- cbind(
    sub.mat.copy.arm$SampleID,
    sub.mat.copy.arm$Chr,
    sub.mat.copy.arm$Start,
    sub.mat.copy.arm$End,
    sub.mat.copy.arm$nProbes,
    sub.mat.copy.arm$cn
  )

  colnames(sub.mat.cn) <- c("Sample", "Chrom", "Start", "End", "Num.probes", "val")

  GD <- genome.doub.sig(
    sample = TCGA.barcode, seg.mat.minor = sub.mat.minor,
    seg.mat.copy = sub.mat.cn, number.of.sim = 10000
  )
  
  GD.status <- fun.GD.status(
    GD.pval = GD$pval,
    ploidy.val = round(sub.mat.copy$Ploidy[1])
  )
  
  GD.pval = GD$pval
  
  #write.GD
  GD.tsv <- paste(patient.folder, "/", TCGA.barcode, ".GD.tsv", sep = "")
  GD.out <- data.frame(
    patient = TCGA.barcode,
    prop.major.even.obs = round(GD$prop.major.even.obs,4),
    ploidy = round(sub.mat.copy$Ploidy[1], 3),
    GD.pval = GD$pval,
    GD.status = GD.status
    
  )
  
  write.table(GD.out,
              sep = "\t",
              quote = FALSE,
              col.names = TRUE,
              row.names = FALSE,
              file = GD.tsv
  )
  
  # Timing the GD.
  GD.time <- genome.doub.timing(sample = TCGA.barcode, seg.mat.minor = sub.mat.minor,
                     seg.mat.copy = sub.mat.cn, loss.level = 3)
  
  #write chromosomal arm levels status.
  write.table(GD.time$sub.mat.arms,
              sep = "\t", quote = FALSE,col.names = TRUE, row.names = FALSE,
              file = paste(patient.folder, "/", TCGA.barcode, ".arms.tsv", sep = "")
  )
  
  #write chromosomal arm levels status.
  write.table(GD.time$GD.timing,
              sep = "\t", quote = FALSE,col.names = TRUE, row.names = FALSE,
              file = paste(patient.folder, "/", TCGA.barcode, ".GD.Timing.tsv", sep = "")
  )
  
  
  TCGA.purity <- as.character(unique(sub.mat.copy[, grep("Aberrant", colnames(sub.mat.copy))]))

  if (TCGA.purity > 1) {
    stop("\nYour tumour content is greater than 1! \nYou probably don't have an ASCAT estimate. \nThere is an alernative script for this scenario.")
  }


  #### EarlyorLate implementation #####

  TCGA.earlyLate <- earlyORlate(patient = TCGA.barcode, complete.mutation.table = mut.table, purity = TCGA.purity)
  TCGA.earlyLate <- cbind(TCGA.earlyLate, GD.pval, GD.status, TCGA.purity)
  
  # removing NA's data.
  if( any(is.na(TCGA.earlyLate$ccf)) ){
    message("Removing NA's: ", sum(is.na(TCGA.earlyLate$ccf)) )
    TCGA.earlyLate = TCGA.earlyLate[ !is.na(TCGA.earlyLate$ccf) ,]
  }



  # Let's choose the columns of interest
  TCGA.earlyLate.out <- TCGA.earlyLate
  TCGA.earlyLate.out$comb.timing <- NA
  TCGA.earlyLate.out[TCGA.earlyLate.out$absolute.ccf.0.95 >= sub.clonal.cut.off & !TCGA.earlyLate.out$timing %in% c("late"), "comb.timing"] <- "Early"
  TCGA.earlyLate.out[TCGA.earlyLate.out$absolute.ccf.0.95 < sub.clonal.cut.off | TCGA.earlyLate.out$timing %in% c("late"), "comb.timing"] <- "Late"
  col.names <- c("patient", "TCGA.purity", "mutation_id", "Reference_Base", "Alternate_Base", "ref_counts", "var_counts", "obs.VAF", "normal_cn", "minor_cn", "major_cn", "mut.multi", "Exp.Cpn", "Exp.Cpn.Likelihood", "absolute.ccf", "absolute.ccf.0.05", "absolute.ccf.0.95", "prob.clonal", "prob.subclonal", "GD.status", "timing", "comb.timing")

  TCGA.earlyLate.out <- TCGA.earlyLate.out[, col.names]

  # let's write the file to a location

  earlylate.tsv <- paste(patient.folder, "/", TCGA.barcode, ".earlylate.tsv", sep = "")
  earlylate.out <- apply(TCGA.earlyLate.out, 2, as.character)
  write.table(earlylate.out,
    sep = "\t",
    quote = FALSE,
    col.names = TRUE,
    row.names = FALSE,
    file = earlylate.tsv
  )

  # plot what you've just found for the given patient

  early.late.pdf <- paste(patient.folder, "/", TCGA.barcode, ".earlylate.pdf", sep = "")

  TCGA.earlyLate <- TCGA.earlyLate[!is.na(TCGA.earlyLate$absolute.ccf.0.95), ]

  # add extra column to

  sub.mat.copy <- cbind(sub.mat.copy, sub.mat.copy$nA + sub.mat.copy$nB)
  colnames(sub.mat.copy)[ncol(sub.mat.copy)] <- "Copy.Number"
  sub.mat.copy <- cbind(sub.mat.copy, apply(cbind(sub.mat.copy$nA, sub.mat.copy$nB), 1, min))
  colnames(sub.mat.copy)[ncol(sub.mat.copy)] <- "min.allele"



  colnames(sub.mat.copy)[2] <- "Chromosome"
  colnames(sub.mat.copy)[3] <- "StartPosition"
  colnames(sub.mat.copy)[4] <- "EndPosition"
  colnames(sub.mat.copy)[5] <- "nr.probes"

  if (plotting) {
    pdf(early.late.pdf)
    par(mar = c(13, 4, 13, 4))
    par(lend = 1)
    plot.EarlyOrLate(
      seg.mat.patient = sub.mat.copy,
      TCGA.earlyLate = TCGA.earlyLate,
      TCGA.purity = TCGA.purity,
      TCGA.barcode = TCGA.barcode,
      sub.clonal = sub.clonal.cut.off
    )

    dev.off()

    ccf.plot <- paste(patient.folder, "/", TCGA.barcode, ".ccfplot.pdf", sep = "")

    pdf(ccf.plot)
    plot.TCGA.sample.ccfs(TCGA.earlyLate = TCGA.earlyLate, clonal.cut.off = 1)
    dev.off()
  }
}
