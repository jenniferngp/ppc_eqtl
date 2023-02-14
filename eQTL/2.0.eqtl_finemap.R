setwd("/projects/PPC/analysis/ppc_eqtls")

source("scripts/packages.R"  )
source("scripts/functions.R" )
source("scripts/input_data.R")

library(coloc)

option_list = list(make_option("--tiss", type = "character", help = "tissue"),
                   make_option("--taskid", type = "integer", help = "tadk id"),
                   make_option("--type", type = "character", help = "eqtl or sqtl"))

opt_parser = OptionParser(option_list = option_list)
opt        = parse_args(opt_parser)
tiss       = opt$tiss
taskid     = opt$taskid
type       = opt$type

if (type == "eqtl")
{
    if (tiss %like% "Islets")
    {
        genes = readLines("input/eqtl/egene_lists/Islets.egenes.txt")
        gene_id = genes[taskid]

        qtl_file = paste("input/eqtl/eqtls_by_gene", paste(tiss, gene_id, "txt", sep = "."), sep = "/")
    } else 
    {
        genes = readLines(paste("input/eqtl/egene_lists", paste(tiss, "v8.egenes.txt", sep = "."), sep = "/"))
        gene_id = genes[taskid]

        number_samples = fread("/frazer01/home/jennifer/references/gtex/GTEx_v8_eQTL_Number_Samples.txt", data.table = F)   

        qtl_file = paste("input/eqtl/eqtls_by_gene", paste(tiss, "v8_eQTL", gene_id, "txt", sep = "."), sep = "/")
    }
    
} else 
{
    if (tiss %like% "Islets")
    {
        genes = readLines("input/eqtl/egene_lists/Islets.exon_egenes.txt")
        gene_id = genes[taskid]

        qtl_file = paste("input/eqtl/eqtls_by_gene", paste(tiss, gene_id, "txt", sep = "."), sep = "/")
    } else 
    {
        genes = readLines(paste("input/eqtl/egene_lists", paste(tiss, "v8.sgenes.txt", sep = "."), sep = "/"))
        gene_id = genes[taskid]

        number_samples = fread("/frazer01/home/jennifer/references/gtex/GTEx_v8_sQTL_Number_Samples.txt", data.table = F) 
        qtl_file = paste("input/eqtl/eqtls_by_gene", paste(tiss, "v8_eQTL", gene_id, sep = "."), sep = "/")
    }
    
}

if (tiss %like% "Islets")
{
   out_file = paste("pipeline/3.4.coloc_adult/eqtl_finemap", paste(tiss, gene_id, "txt", sep = "."), sep = "/") 
} else 
{
    out_file = paste("pipeline/3.4.coloc_adult/eqtl_finemap", paste(tiss, "v8_eQTL", gene_id, "txt", sep = "."), sep = "/")
}


#if (file.exists(out_file) == T) { stop(paste(out_file, "done")) }

# Get N
if (!tiss %like% "Islets")
{
    N = number_samples[number_samples$tiss == tiss,]$nsamples 
} else 
{
    N = 420
}

message(paste("running", tiss, gene_id))
message(paste("N samples:", N))
message(paste("opening", qtl_file))

# Get QTL data
if (!tiss %like% "Islets")
{
    qtl_data = fread(qtl_file, data.table = F)
    
    if (nrow(qtl_data) == 0) { stop(paste("no data", gene_id)) }
    
    qtl_data$variant_id = unlist(lapply(qtl_data$variant_id, function(x) { paste(unlist(strsplit(x, "_"))[1:4], collapse = "_") }))
    colnames(qtl_data)[which(colnames(qtl_data) == "af")] = "maf"
} else 
{
    qtl_data = fread(qtl_file, data.table = F)
    
    if (nrow(qtl_data) == 0) { stop("no QTL data") }
    
    if (tiss %like% "Exons")
    {
        colnames(qtl_data) = c("GeneID", "ExonsID", "GeneChr", "TSS", "DistanceToGene", "rsID", "SNPchr", "SNPposition", "REF", "ALT", "FreqREF", "FreqALT",
                           "NA", "Pvalue", "Slope", "Lead", "T-stat", "SE")
        
        qtl_data$FreqALT = as.numeric(qtl_data$FreqALT)
        qtl_data$FreqREF = as.numeric(qtl_data$FreqREF)
        qtl_data$maf = ifelse(qtl_data$FreqALT < qtl_data$FreqREF, qtl_data$FreqALT, qtl_data$FreqREF)
        qtl_data$variant_id = paste(qtl_data$SNPchr, qtl_data$SNPposition, qtl_data$REF, qtl_data$ALT, sep = "_")
        qtl_data$pval_nominal = qtl_data$Pvalue
    }
}

qtl_data$maf = as.numeric(qtl_data$maf)

# Make sure MAF is between 0 and 1
qtl_data = qtl_data[qtl_data$maf > 0 & qtl_data$maf < 1 & !is.na(qtl_data$maf),] %>% arrange(pval_nominal)
qtl_data = qtl_data[duplicated(qtl_data$variant_id) == F,]
qtl_data = qtl_data %>% filter(!variant_id %like% "NA")

# Finemap
out = finemap.abf(list(snp = qtl_data$variant_id,
                   pvalues = qtl_data$pval_nominal, 
                   N = N, 
                   MAF = qtl_data$maf, 
                   type = "quant"),
                  p1 = 1 / nrow(qtl_data))


# Label snps in 99% credible set
cred = suppressWarnings(out %>% arrange(desc(SNP.PP)) %>% dplyr::slice(seq_len(which.max(cumsum(SNP.PP) >= 0.99))))
out$in_credible = ifelse(out$snp %in% cred$snp, T, F)
out$changed_prior = F

message(paste("topsnp pp:", max(out$SNP.PP)))
message(paste(capture.output(table(out$in_credible)), collapse = "\n"))

# Write
fwrite(out, out_file, row.names = F, sep = "\t")
message(paste(out_file, "written!"))
