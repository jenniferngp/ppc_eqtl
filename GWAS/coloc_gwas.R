setwd("/frazer01/projects/PPC/analysis/ppc_eqtls")

# Load libraries
source("scripts/packages.R"  )
source("scripts/input_data.R")
source("scripts/functions.R" )
source("scripts/3.3.gwas/functions.R")

suppressPackageStartupMessages(library(coloc))

# Load coloc functions
source("scripts/3.3.gwas/run_ukbb.R")
source("scripts/3.3.gwas/run_magic.R")
source("scripts/3.3.gwas/run_t1d.R")
source("scripts/3.3.gwas/run_mahajan.R")

# Read arguments
option_list   = list(make_option("--taskid"          , type="integer"  , default=0      , help="SGE task ID"       , metavar="character"),
					 make_option("--analysis"        , type="character", default="gene" , help="gene or isoform"   , metavar="character")
					) 

opt_parser        = OptionParser(option_list=option_list)
opt               = parse_args(opt_parser)
taskid            = opt$taskid
analysis          = opt$analysis

# Get Gene Info
if (analysis == "gene")
{
    geneinfo   = fread("pipeline/1.2.expression/gene_info.txt", data.table = F)
    genes      = readLines("pipeline/3.3.gwas/meta/ppc.peer20.egenes.txt")
    gene_id    = genes[taskid]
    geneinfo   = geneinfo[geneinfo$transcript_id %like% gene_id,]

} else {
    geneinfo   = fread("pipeline/1.2.expression/isoform_info.txt", data.table = F)
    genes      = readLines("pipeline/3.3.gwas/meta/ppc.peer20.eisofs.txt")
    gene_id    = genes[taskid]
    geneinfo   = geneinfo[geneinfo$transcript_id %like% gene_id,]
}

# Other input
manifest = fread("pipeline/3.3.gwas/meta/manifest_to_use_April2022.txt", data.table = F)
ppc_maf  = fread("input/genotype/ppc_formated.frq", data.table = F)

# Run colocatlization
outlist = list()
for (trait in manifest$trait_id)
{
    gene_id = unlist(strsplit(gene_id, "[.]"))[1] 
    outlist[[trait]] = run_coloc_ukbb(gene_id, trait, ppc_maf, geneinfo, manifest, analysis)
}

outfile = paste("pipeline/3.3.gwas/coloc_height", paste("ppc", gene_id, "robj", sep = "."), sep = "/")
save(outlist, file = outfile)

message(paste(outfile, "written!"))

