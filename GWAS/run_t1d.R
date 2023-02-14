run_coloc_gaulton_t1d = function(gene_id, gwas, ppc_maf, geneinfo, manifest, analysis, qtl_data, N)
{
    chr = geneinfo$chrom
    
    coord    = paste(sub("chr", "", geneinfo$chrom), ":", max(c(0, geneinfo$start - 1e6)), "-", geneinfo$end + 1e6, sep = "")

    gwas_file = paste("/reference/public/GCST90014023/hg19", paste0(chr, "_t1d_sorted.txt.gz"), sep = "/")
    
    tmp_file = paste("pipeline/3.3.gwas/scratch", paste("TMP", gene_id, gwas, "txt", sep = "."), sep = "/")
    
    command = paste("tabix", gwas_file, coord, ">", tmp_file)
    message(command)
    system(command)
    
    gwas_data        = fread(tmp_file, sep = "\t", header = FALSE , data.table = FALSE)
    myhead           = unlist(strsplit(system(paste("zcat", gwas_file, "|", "head", "-n", 1), intern = TRUE), split = "\t"))
    colnames(gwas_data) = myhead

    gwas_data$base_pair_location = gwas_data$hg19_position

    tryCatch(
    {
        #N = 107

        #qtl_data = get_qtl_data(gene_id, analysis)

        gwas_data$id = paste("VAR", gwas_data$chromosome, gwas_data$base_pair_location, sep = "_")
        qtl_data$id  = paste("VAR", qtl_data$chrom, qtl_data$pos, sep = "_")

        totest       = merge(qtl_data, gwas_data, by = "id", suffixes = 1:2)

        # Check if ref and alt are the same
        totest1 = totest[totest$ref == totest$effect_allele & totest$alt == totest$other_allele,]

        totest2 = totest[totest$ref == totest$other_allele & totest$alt == totest$effect_allele,]
        totest2$beta2 =   - totest2$beta2
        totest2$effect_allele_frequency = 1 - totest2$effect_allele_frequency
        totest = rbind(totest1, totest2)

        # SNP ID
        totest$snp_id = paste("VAR", totest$chrom, totest$pos, totest$ref, totest$alt, sep = "_")

        type_list = list()
        for (qtl_type in unique(totest$type))
        {
            tocoloc = totest[totest$type == qtl_type,]
            tocoloc = tocoloc[tocoloc$effect_allele_frequency > 0 & tocoloc$effect_allele_frequency < 1 & tocoloc$af > 0 & tocoloc$af < 1 & 
                              !is.na(tocoloc$effect_allele_frequency) & !is.na(tocoloc$af) & duplicated(tocoloc$snp_id) == F,]

            coloc_mapped = suppressWarnings(coloc.abf(dataset1 = list(snp = tocoloc$snp_id, 
                                                                      pvalues = tocoloc$pval, 
                                                                      N = N, 
                                                                      MAF = tocoloc$af, 
                                                                      type = "quant"),
                                                      dataset2 = list(snp = tocoloc$snp_id, 
                                                                      pvalues = tocoloc$p_value, 
                                                                      N = tocoloc$sample_size, 
                                                                      MAF = tocoloc$effect_allele_frequency, 
                                                                      type = "quant")))

            type_list[[paste("type", qtl_type)]] = coloc_mapped
        }

        return(type_list)

    }, error = function(cond)
    {
        message(paste(gene_id, cond))
    }, warning = function(cond)
    {
        message("")
    })
    
    system(paste("rm", tmp_file))
}