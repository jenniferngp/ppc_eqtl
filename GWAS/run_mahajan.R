# gene_id = gene identifier
# gwas = trait id
# ppc_maf = MAF for each variant
# geneinfo = gene information
# manifest = contains file path for gwas summary statistic
# analysis = gene or isoform

run_coloc_mahajan = function(gene_id, gwas, ppc_maf, geneinfo, manifest, analysis, qtl_data, N)
{
    gwas_file = manifest[manifest$trait_id == gwas,]$filename

    coord    = paste(geneinfo$chrom, ":", max(c(0, geneinfo$start - 500000)), "-", geneinfo$end + 500000, sep = "")

    tmp_file = paste("pipeline/3.3.gwas/scratch", paste("TMP", gene_id, gwas, "txt", sep = "."), sep = "/")

    command = paste("tabix", gwas_file, coord, ">", tmp_file)
    message(command)
    system(command)

    gwas_data           = fread(tmp_file, sep = "\t", header = FALSE , data.table = FALSE)
    colnames(gwas_data) = c("Chr", "Pos", "SNP", "EA", "NEA", "EAF", "Beta", "SE", "Pvalue", "Neff")
    gwas_data$Chr       = unlist(lapply(gwas_data$Chr, function(x) { gsub("chr", "", x) }))

    tryCatch(
    {
        #N = 107

        #qtl_data = get_qtl_data(gene_id, analysis)

        gwas_data$id = paste("VAR", gwas_data$Chr, gwas_data$Pos, sep = "_")
        qtl_data$id  = paste("VAR", qtl_data$chrom, qtl_data$pos, sep = "_")

        totest       = merge(qtl_data, gwas_data, by = "id", suffixes = 1:2)

        message(paste("intersecting snps:", length(intersect(gwas_data$id, qtl_data$id))))

        if (nrow(totest) > 0)
        {
            # Check if ref and alt are the same
            totest1 = totest[totest$ref == totest$EA & totest$alt == totest$NEA,]
            totest2 = totest[totest$ref == totest$NEA & totest$alt == totest$EA,]

            totest2$Beta =   - totest2$Beta
            totest = rbind(totest1, totest2)

            # SNP ID
            totest$snp_id = paste("VAR", totest$chrom, totest$pos, totest$ref, totest$alt, sep = "_")

            # Get sample size
            n = 74124 + 824006
            cases_fr = 74124 / n

            type_list = list()
            for (qtl_type in unique(totest$type))
            {
                tocoloc = totest[totest$type == qtl_type,]
                tocoloc = tocoloc[tocoloc$af > 0 & tocoloc$af < 1 & !is.na(tocoloc$af) & duplicated(tocoloc$snp_id) == F & 
                                  tocoloc$EAF > 0 & tocoloc$EAF < 1 & !is.na(tocoloc$EAF),]

                coloc_mapped = suppressWarnings(coloc.abf(dataset1 = list(snp = tocoloc$snp_id, 
                                                                          pvalues = tocoloc$pval, 
                                                                          N = N, 
                                                                          type = "quant", 
                                                                          MAF = tocoloc$af),
                                                          dataset2 = list(snp = tocoloc$snp_id, 
                                                                          pvalues = tocoloc$Pvalue, 
                                                                          N = n, 
                                                                          type = "cc", 
                                                                          s = cases_fr, 
                                                                          MAF = tocoloc$EAF)))


                type_list[[paste("type", qtl_type)]] = coloc_mapped
            }
            
            return(type_list)
        
        } else 
        {
            message("merging gwas and eqtl resulted in no overlap")
        }
    }, error = function(cond)
    {
        message(paste(gene_id, cond))
    }, warning = function(cond)
    {
        message(cond)
    })

}