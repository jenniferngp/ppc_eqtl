run_coloc_ukbb = function(gene_id, gwas, ppc_maf, geneinfo, manifest, analysis, qtl_data, N)
{   
    tmp_file  = paste("pipeline/3.3.gwas/scratch", paste("TMP", gene_id, gwas, "txt", sep = "."), sep = "/")
    gwas_file = manifest[manifest$trait_id == gwas,]$filename
    coord     = paste(sub("chr", "", geneinfo$chrom), ":", max(c(0, geneinfo$start - 500000)), "-", geneinfo$end + 500000, sep = "")
    command   = paste("tabix", gwas_file, coord, ">", tmp_file)
    message(command)
    system(command)

    trait_type = manifest[manifest$trait_id == gwas, "trait_type"]

    if (file.size(tmp_file) > 0)
    {
        gwas_data           = fread(tmp_file, sep = "\t", header = FALSE , data.table = FALSE)
        myhead              = unlist(strsplit(system(paste("zcat", gwas_file, "|", "head", "-n", 1), intern = TRUE), split = "\t"))
        colnames(gwas_data) = myhead
    } else {
        message("no snps in gwas")
        stop()
    }

    #N = 107

    tryCatch(
    {
        #qtl_data = get_qtl_data(gene_id, analysis)

        gwas_data$id   = paste("VAR", gwas_data$chr, gwas_data$pos, sep = "_")
        qtl_data$id    = paste("VAR", qtl_data$chrom, qtl_data$pos, sep = "_")

        totest         = merge(qtl_data, gwas_data, by = "id", suffixes = 1:2)
        
        if (nrow(totest) > 0)
        {
            n = rowSums(manifest[manifest$trait_id == gwas, c("n_cases_full_cohort_both_sexes", "n_controls_total")])[[1]] 

            # Filter Heterogeneity
            totest = totest[ totest$pval_heterogeneity > 1e-6,]

            if(!trait_type %in% c("biomarkers", "continuous"))
            {
                totest$af_meta = totest$af_controls_meta
                type2coloc     = "cc"
                cases_fr       = manifest[manifest$trait_id == gwas, "n_cases_full_cohort_both_sexes"] / n
            }

            # Check if ref and alt are the same
            totest1 = totest[totest$ref1 == totest$ref2 & totest$alt1 == totest$alt2,]
            totest2 = totest[totest$ref1 == totest$alt2 & totest$alt1 == totest$ref2,]

            totest2[,paste("beta", "meta", sep = "_")] =   - totest2[,paste("beta", "meta", sep = "_")]
            totest2[,paste("af"  , "meta", sep = "_")] = 1 - totest2[,paste("af"  , "meta", sep = "_")]

            totest = rbind(totest1, totest2)

            # SNP ID
            totest$snp_id = paste("VAR", totest$chr, totest$pos1, totest$ref1, totest$alt1, sep = "_")

            type_list = list()
            for (qtl_type in unique(totest$type)[which(unique(totest$type) %in% c(0,1,2,3,4,5))])
            {
                tocoloc = totest[totest$type == qtl_type,]

                message(paste("QTL type:", qtl_type, nrow(tocoloc), "=========="))

                tocoloc = tocoloc[tocoloc$af_meta > 0 & tocoloc$af_meta < 1 & tocoloc$af > 0 & tocoloc$af < 1 & 
                              !is.na(tocoloc$af_meta) & !is.na(tocoloc$af) & duplicated(tocoloc$snp_id) == F,]

                if(trait_type != "phecode")
                {
                    coloc_mapped = suppressWarnings(coloc.abf(dataset1 = list(snp = tocoloc$snp_id, 
                                                                              pvalues = tocoloc$pval, 
                                                                              N = N, 
                                                                              MAF = tocoloc$af, 
                                                                              type = "quant"),
                                                              dataset2 = list(snp = tocoloc$snp_id, 
                                                                              pvalues = tocoloc$pval_meta, 
                                                                              N = n, 
                                                                              MAF = tocoloc$af_meta, 
                                                                              type = "quant")))
                }else
                {
                    tocoloc$af_meta = tocoloc$af_controls_meta
                    coloc_mapped = suppressWarnings(coloc.abf(dataset1 = list(snp = tocoloc$snp_id, 
                                                                              pvalues = tocoloc$pval, 
                                                                              N = N, 
                                                                              MAF = tocoloc$af, 
                                                                              type = "quant"),
                                                              dataset2 = list(snp = tocoloc$snp_id, 
                                                                              pvalues = tocoloc$pval_meta, 
                                                                              N = n, 
                                                                              MAF = tocoloc$af_meta, 
                                                                              type = "cc"   , 
                                                                              s = cases_fr)))
                }
                type_list[[paste("type", qtl_type)]] = coloc_mapped

            }

            return(type_list)
        } else
        {
            message("no snps between gwas and qtl")
        }
                                                
    }, error = function(cond)
    {
        message(paste(gene_id, cond))
    })
}