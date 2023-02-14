# qsub -N coloc -o logs -e logs -V -cwd -pe smp 1 -l short -tc 300 -t 1-4065:1 scripts/3.4.coloc_adult/allpair_coloc.sh

setwd("/frazer01/projects/PPC/analysis/ppc_eqtls")

source("scripts/packages.R"  )
source("scripts/input_data.R")
source("scripts/functions.R" )
source("scripts/coloc_functions.R")
source("scripts/3.4.coloc_adult/functions.R")

suppressMessages(library(coloc        ))

option_list = list(make_option("--taskid", type = "integer", help = "SGE task id"),
                   make_option("--type"  , type = "character", help = "gene or isoform"))

opt_parser = OptionParser(option_list = option_list)
opt        = parse_args(opt_parser)
taskid     = opt$taskid
type       = opt$type

message("opening maf")
ppc_maf = fread("input/genotype/ppc_formated.frq")
colnames(ppc_maf)[which(colnames(ppc_maf) == "maf")] = "af"

suppressWarnings(dir.create("pipeline/3.4.coloc_adult/allpairs"))
suppressWarnings(dir.create("pipeline/3.4.coloc_adult/allpairs/meta"))
suppressWarnings(dir.create("pipeline/3.4.coloc_adult/allpairs/data"))
suppressWarnings(dir.create("pipeline/3.4.coloc_adult/allpairs/summary"))

# Making gene info
geneinfo = fread("/reference/private/Gencode.v34lift37/gene_info.txt") %>% 
    mutate(transcript_id = gene_id, 
           transcript_id = simplify_id(transcript_id), 
           gene_id = simplify_id(gene_id)) %>% 
    dplyr::select(transcript_id, gene_id, gene_name, chrom, start, end)
isofinfo = fread("/reference/private/Gencode.v34lift37/isoform_info.txt") %>% 
    mutate(transcript_id = simplify_id(transcript_id), 
           gene_id = simplify_id(gene_id)) %>% 
    dplyr::select(transcript_id, gene_id, gene_name, chrom, start, end)

geneinfo = rbind(geneinfo, isofinfo)
colnames(geneinfo) = c("transcriptid", "geneid", "gene_name", "chrom", "start", "end")

# Reading Inputs
if (type == "gene")
{
    outfile = "pipeline/3.4.coloc_adult/allpairs/meta/allpairs.genes_to_coloc.txt" # See 3.4.coloc_adult/08.PPCcoloc
    number_samples = fread("/frazer01/home/jennifer/references/gtex/GTEx_v8_eQTL_Number_Samples.txt", data.table = F)

} else 
{
    outfile = "pipeline/3.4.coloc_adult/allpairs/meta/allpairs.isofs_to_coloc.txt" # See 3.4.coloc_adult/08.PPCcoloc
    number_samples = fread("/frazer01/home/jennifer/references/gtex/GTEx_v8_sQTL_Number_Samples.txt", data.table = F)
}

# Getting Gene ID
genes   = readLines(outfile)
transcript_id = genes[taskid]
gene_id = geneinfo[geneinfo$transcriptid == transcript_id,]$geneid

message(paste("Transcript_ID:", transcript_id))
message(paste("Gene_ID:", gene_id))

# Get All Genes to Test
Chr           = geneinfo[geneinfo$transcriptid == gene_id,]$chrom
Start         = max(0, geneinfo[geneinfo$transcriptid == gene_id,]$start - 2e6)
End           = geneinfo[geneinfo$transcriptid == gene_id,]$start + 2e6
genes_to_test = geneinfo[geneinfo$start >= Start & geneinfo$end <= End & geneinfo$chrom == Chr & geneinfo$transcriptid %like% "ENSG",]
genes_to_test = simplify_id(unique(genes_to_test$transcriptid))

message(paste("running", length(unique(genes_to_test)), "genes"))

ppc_qtl = get_qtl_data("ppc", transcript_id, type) %>% select(chrom, pos, ref, alt, af, pval, beta, type, id) %>% mutate(id = paste("VAR", chrom, pos, sep = "_"))

qtl_files = list.files("/projects/PPC/analysis/ppc_eqtls/input/eqtl/eqtls_by_gene", full.names = T)

gene_list = list()
for (gene2 in genes_to_test)
{
    message(paste(Sys.time(), gene2, "========="))
    
    if (type == "gene") 
    { 
        tissues = c("Islets", "Pancreas.v8_eQTL.hg19")
    } else 
    {
        tissues = c("Islets_Exons", "Pancreas.v8_sQTL.hg19")
    }
    
    coloc_list = list()
    coloc_list[["ppc"]] = ppc_qtl
    
	for (tiss in tissues)
    {
        if (type == "gene")
        {
            files = paste("input/eqtl/eqtls_by_gene", paste(tiss, gene2, "txt", sep = "."), sep = "/")
        } else 
        {
            files = qtl_files[which(qtl_files %like% tiss & qtl_files %like% gene2)]
        }
        
        for (file in files)
        {
            if (file.exists(file) == T)
            {
                message(file)
                out       = fread(file, data.table = F)
                colnames(out)[which(colnames(out) == "maf")] = "af"
                out$chrom = gsub("chr", "", unlist(lapply(out$variant_id, function(x) { unlist(strsplit(x, "_"))[1] })))
                out$pos   = unlist(lapply(out$variant_id, function(x) { unlist(strsplit(x, "_"))[2] }))
                out$ref   = unlist(lapply(out$variant_id, function(x) { unlist(strsplit(x, "_"))[3] }))
                out$alt   = unlist(lapply(out$variant_id, function(x) { unlist(strsplit(x, "_"))[4] }))
                out$id    = paste("VAR", out$chrom, out$pos, out$ref, out$alt, sep = "_")
                out$pval  = out$pval_nominal
                out$beta  = out$slope
                out$se    = out$slope_se

                out = out %>% select(chrom, pos, ref, alt, af, pval, beta, id)  
                out$id = paste("VAR", out$chrom, out$pos, sep = "_")
                
                id = gsub("/projects/PPC/analysis/ppc_eqtls/input/eqtl/eqtls_by_gene/", "", file)
                id = gsub("input/eqtl/eqtls_by_gene/", "", file)
                id = gsub(".txt", "", id)
                id = gsub(paste0(tiss, ".v8_eQTL."), "", id)
                id = gsub(paste0(tiss, "."), "", id)
                
                coloc_list[[paste(tiss, id)]] = out
            }  
        }
    }

    coloc_out_list = list()

    for (adult_id in names(coloc_list)[which(names(coloc_list) != "ppc")])
    {
        adult = unlist(strsplit(adult_id, " "))[1]
        gene = unlist(strsplit(adult_id, " "))[2]
        
        if (!adult %like% "Islets")
        {
            N_adult = number_samples[number_samples$tiss == "Pancreas",]$nsamples
        } else
        {
            N_adult = 420
        }

        message(paste("Colocalizing with", adult_id, N_adult))

        totest = merge(coloc_list[["ppc"]], coloc_list[[adult_id]], by = "id", suffixes = 1:2)

        if (nrow(totest) > 0)
        {
            totest = check_beta(totest)
            totest$snp_id = paste("VAR", totest$chrom1, totest$pos1, totest$ref1, totest$alt1, sep = "_")

            types = as.vector(unique(totest$type))

            type_list = list()
            for (qtl_type in types)
            {
                message(paste(adult_id, qtl_type))
                tocoloc = totest[totest$type == qtl_type,] %>% mutate(af1 = as.numeric(af1), af2 = as.numeric(af2))
                tocoloc = tocoloc[tocoloc$af1 > 0 & 
                                  tocoloc$af1 < 1 & 
                                  tocoloc$af2 > 0 & 
                                  tocoloc$af2 < 1 &
                                  !is.na(tocoloc$af1) & 
                                  !is.na(tocoloc$af2) &
                                  duplicated(tocoloc$snp_id) == F,]

                coloc = suppressWarnings(coloc.abf(dataset1 = list(snp = tocoloc$snp_id, pvalues = tocoloc$pval1, N = 107    , type = "quant", MAF = tocoloc$af1),
                                                   dataset2 = list(snp = tocoloc$snp_id, pvalues = tocoloc$pval2, N = N_adult, type = "quant", MAF = tocoloc$af2)))
                
                type_list[[paste("type", qtl_type)]] = coloc
            }
            coloc_out_list[[adult_id]] = type_list
        }
    }
    gene_list[[gene2]] = coloc_out_list
    
}

outfile = paste("pipeline/3.4.coloc_adult/allpairs/data", paste("ppc", transcript_id, "robj", sep = "."), sep = "/")
save(gene_list, file = outfile)

message(paste(outfile, "saved!"))