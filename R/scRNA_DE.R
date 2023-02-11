if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("DESeq2")
#devtools::install_github('kevinblighe/EnhancedVolcano')

require(data.table)
library(tidyr)
library(tidyverse)
library(dplyr)
library(DESeq2)
library(EnhancedVolcano)
library(optparse)
library(future)

options(future.globals.maxSize=1677721600000)
option_list <- list(
  make_option(c("-l", "--snvlist"), action = "store", type ="character", 
              default = NULL, help = "file needs to be header with SNV, and
              row with the SNV's genomic coordinates."),
  make_option(c("-m", "--genematrix"), action = "store", type ="character", 
              default = NULL, help = "file needs to be GE matrix."),
  make_option(c("-r", "--screadcounts"), action = "store", type ="character", 
              default = NULL, help = "file needs to be an scReadCounts type output."),
  make_option(c("-s", "--samplename"), action = "store", type = "character",
              default = NULL, help = "file needs to simply be the sample corresponding to the 
              GE matrix/scReadcounts file, EX: SRR100000")
  #make_option(c("-l", "--legend"), action = "store", type ="character", 
  #            default = "TRUE",
  #            help= "Gererate the graph with/or without legend, Default = TRUE",
              )#)

opt <- parse_args(OptionParser(option_list=option_list, description = "-s, -m, -r, -l are necessary!", 
                               usage = "usage: Rscript DE_scReadCount_071322.R -l <snv list> -m <gene matrix> -r <screadcounts> -s <sample name>"))
if (is.null(opt$samplename)) {
  stop("Please make sure you have provided the file containing the sample name.")
}
if (is.null(opt$genematrix)) {
  stop("Please make sure you have provided the GE matrix.")
}
if (is.null(opt$screadcounts)) {
  stop("Please make sure you have provided the scReadCounts file.")
}
if (is.null(opt$snvlist)) {
  stop("Please make sure you have provided the file containing the SNV coordinates.")
}


sample_list <- fread(opt$samplename, header = F)
sample_id <- sample_list$V1
snv_input <- fread(opt$snvlist, header = T)
screadcounts <- na.omit(fread(opt$screadcounts, header = T))
express_matrix <- na.omit(fread(opt$genematrix, header = T))

# changed express_matrix so col and row numbers match
express_matrix <- data.frame(express_matrix, row.names="gene_id")

snv_input = snv_input %>% select(1:4)
snv_input = snv_input %>% mutate(SNV = paste0(CHROM, ":", POS, "_", REF, ">", ALT)) %>% select(SNV) %>% distinct()

snv_list <- snv_input$SNV

#list of barcodes to use for filtering
barcodelist <- colnames(express_matrix) 

#separate out/expand the Cells column
filt_screadcounts <- screadcounts %>% filter(ReadGroup %in% barcodelist)

#keep only cells found in the GE barcode list and remove screadcounts
vaf_cells <- filt_screadcounts[filt_screadcounts$VAF == 0, ]
rm(screadcounts)

# make directory for output files: not in for loop
folder <- "outputs"

# don't recreate folder if it already exists
if (file.exists(folder)) {
  cat("The folder already exists")
} else {
  dir.create(folder)
}

# set wd to outputs
setwd(folder)

####


#DESEQ2 GE using cells with specific SNVs
#set two for loops to progress through the selected SNVs from input list
#second loop gets the SNV into the correct format for the first loop to match the SNV to the screadcounts data
plan("multiprocess", workers = 4)
future.seed = TRUE
for (input in snv_list) {
  for (string in input) {
    match <- as.data.frame(str_split(input, ":|_|>", simplify = TRUE))
  }

  snv <- filt_screadcounts[filt_screadcounts$CHROM %in% match$V1, ]
  snv <- snv[snv$POS %in% match$V2, ]
  snv <- snv[snv$REF %in% match$V3, ]
  snv <- snv[snv$ALT %in% match$V4, ]
  snv_vaf <- snv[snv$VAF > 0, ]
  total_vaf <- sum(snv_vaf$VAF)
  
  #if selected SNV coordinates correspond to VAF of 0 then the loop will stop, and start over on the next input
  if(total_vaf ==0){
    print("VAF = 0 for selected SNV. Check your input and SNV coordinates. 
          Skipping and moving to next input.")
    next
  }
  
  # A is cells that have the selected SNV with a VAF > 0
  snv_cells_A <- snv_vaf$ReadGroup
  #B is all others/control
  snv_cells_B <- vaf_cells$ReadGroup
  
  #filtering out only the selected snvs and setting the condition for deseq as A
  snv_cells_A <- as.data.frame(snv_cells_A)
  snv_cells_B <- as.data.frame(snv_cells_B)
  snv_cells_A$condition <- "A"
  snv_cells_B$condition <- "B"
  names(snv_cells_A)[1] <- "barcode"
  names(snv_cells_B)[1] <- "barcode"
  
  #binding the two conditional dataframes together 
  dseq_meta <- rbind(snv_cells_A, snv_cells_B)
  dseq_meta$condition <- as.factor(dseq_meta$condition)
  
  #create the GE matrix barcode dataframe for reordering 
  #match the GE barcode to the metadata barcode and get an index for reordering
  coldata <- data.frame(barcode = barcodelist, row.names = "barcode")
  dseq_meta_index <- match(rownames(coldata), dseq_meta$barcode)
  dseq_meta_reorder <- as.data.frame(dseq_meta[dseq_meta_index, ])

  # changed to remove NAs from first row: Kai
  dseq_meta_reorder <- data.frame(dseq_meta_reorder[complete.cases(dseq_meta_reorder),], row.names = "barcode")

  #reordering and dataframe manipulation to get into deseq2 format
  
  ###metadata creation
  
  #create the deseq object
  dds <- DESeqDataSetFromMatrix(countData = express_matrix, colData = dseq_meta_reorder, 
                                design = ~ condition)
  
  #set reference condition
  #here we set B reference/control
  dds$condition <- relevel(dds$condition, ref = "B")
  
  #initialize the differential expression
  dds <- DESeq(dds)

  #initialize the results, with no filtering
  res <- results(dds, independentFiltering = FALSE)
  
  #results by order of adjusted p-value
  res <- res[order(res$padj),]
  #view summary of results that match alpha <0.05
  summary(results(dds, alpha=0.05))
  #normalize counts variable
  normalized_counts <- counts(dds, normalized=TRUE)
  
  #output genes ordered by upregulated genes, set decreasing=FALSE for downregulated
  output_gene <- as.data.frame(res[order(res$log2FoldChange, decreasing=TRUE),])

  # create directory for each SNV data, makes it less cluttered
  SNV_folder <- paste0(gsub(":|>", "_",input))
  
  # don't recreate folder if it already exists
  if (file.exists(SNV_folder)) {
    cat("The folder already exists")
  } else {
    dir.create(SNV_folder)
  }
  
  #  input already has chromosome number, unnecessary to add it again here
  file_name_genes <- paste0("DE_genes_", input, ".tsv")
  file_name_genes <- gsub(":|>", "_", file_name_genes)
  
  # change directory to SNV folder
  cwd <- getwd()
  setwd(SNV_folder)

  #file name format for all plots and gene table, removing ":" and ">"
  file_name_plotMA_shrink <- paste0("DE_MA_Shrink_", input, ".pdf")
  file_name_plotMA_shrink <- gsub(":|>", "_",file_name_plotMA_shrink)
  
  file_name_plotMA <- paste0("DE_MA_", input, ".pdf")
  file_name_plotMA <- gsub(":|>", "_", file_name_plotMA)
  
  file_name_volcano_shrink <- paste0("DE_Volcano_", input, ".pdf")
  file_name_volcano_shrink <- gsub(":|>", "_", file_name_volcano_shrink)

  # write gene names to file
  file.create(file_name_genes, showWarnings = TRUE)
  write.table(x = output_gene, file = file_name_genes, quote = FALSE, sep = "\t", col.names = NA)

  # apply lfcShrink to results
  resLFC<- lfcShrink(dds, coef="condition_A_vs_B", type = "normal", lfcThreshold = 1)
  
  # changed file creation 
  file.create(file_name_volcano_shrink, showWarnings = TRUE)
  pdf(file = file_name_volcano_shrink)
  p <- EnhancedVolcano(resLFC, lab=rownames(resLFC), x="log2FoldChange", y= "pvalue", 
                       pCutoff = 1, title = "Normal vs. SNV" ,xlab = bquote(~Log[2]~ "fold change"), 
                       colAlpha = .5, pointSize = 1.0, labSize = 3.0, col = c("black", "blue", "green", "red"), 
                       legendPosition = "right",
                       legendLabSize = 10.0,
                       legendIconSize = 3.0)
  print(p)
  dev.off()
  
  # added MA plots
  file.create(file_name_plotMA_shrink, showWarnings = TRUE)
  pdf(file = file_name_plotMA_shrink)
  plotMA(resLFC, ylim=c(-2,2), cex=.4)
  dev.off()
  
  print(paste0("SNV ", input, " finished."))
  
  # change back to outputs directory
  setwd(cwd)
  
}

####
print("DONE")

