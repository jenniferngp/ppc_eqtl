### Version history:
### V02: 6/10/2020
### V01: 4/17/2020

liftover  = function(bed)
{
    suppressMessages(library(rtracklayer  ))
    suppressMessages(library(GenomicRanges))
    
    path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
    ch   = import.chain(path)
    
    gr38 = makeGRangesFromDataFrame(bed, T)
    gr19 = liftOver(gr38, ch)
    return(gr19) 
}

psize = function(w, h)
{
    options(repr.plot.height = h, repr.plot.width = w)
}

# General functions

plot_signal_ld = function(data, col, ylab, title)
{
    p = ggplot(data, aes(x = pos / 1e6, y = -log10(pval), color = ld_group)) + geom_point(size = 5) + 
        jn_theme + 
        scale_color_manual(values = c("red", "orange", "green", "skyblue", "navy")) + 
#         scale_color_gradientn(colors = c("lightgrey", col)) + 
#         geom_point(data = data %>% filter(ld_group == "r2 < 0.2"), aes(x = pos / 1e6, y = -log10(pval), color = ld_group), size = 4) + 
        xlab("Position (Mb)") + ylab(ylab) + 
        geom_point(data = data, aes(x = pos / 1e6, y = -log10(pval)), shape = 1, color = "black", stroke = 0.5, size = 5) + 
        geom_hline(yintercept = -log10(0.1 / nrow(data)), linetype = "dashed")  +
        geom_vline(xintercept = as.numeric(unlist(strsplit(as.character(data[data$lead == T,]$id), "_"))[3]) / 1e6 , linetype = "dashed")
    
    if (title != "NULL") { p = p + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))}
    
    return(p)
}


get_ldgroups = function(data)
{
    data$ld_group = ifelse(data$LD >= 0.8, "r2 > 0.8", NA)
    data$ld_group = ifelse(data$LD >= 0.6 & data$LD < 0.8, "0.6 > r2 > 0.8", data$ld_group)
    data$ld_group = ifelse(data$LD >= 0.4 & data$LD < 0.6, "0.8 > r2 > 0.6", data$ld_group)
    data$ld_group = ifelse(data$LD >= 0.2 & data$LD < 0.4, "0.4 > r2 > 0.2", data$ld_group)
    data$ld_group = ifelse(data$LD < 0.2, "r2 < 0.2", data$ld_group)
    
    data$ld_group = factor(data$ld_group, levels = rev(c("r2 < 0.2", "0.4 > r2 > 0.2", "0.8 > r2 > 0.6", "0.6 > r2 > 0.8", "r2 > 0.8")))
    return(data)
}

get_qtl = function(transcript_id, qtl_type)
{
    if (transcript_id %like% "ENST")
    {
        egene_dir   = "pipeline/3.2.eqtls/eqtls_by_gene/ppc_eqtls_peer20.isoform"
    } else 
    {
        egene_dir   = "pipeline/3.2.eqtls/eqtls_by_gene/ppc_eqtls_peer20.gene"
    }
    
    egene_files = list.files(egene_dir)
    qtl = fread(paste(egene_dir, egene_files[which(egene_files %like% transcript_id & egene_files %like% "qtl")], sep = "/"), data.table = F) %>% 
        dplyr::filter(type == qtl_type) 
    lead = qtl %>% dplyr::arrange(pval) %>% head(1) %>% dplyr::select(id) %>% unlist()
    qtl$lead = ifelse(qtl$id == lead, T, F)
    
    qtl = merge(get_ld(transcript_id, lead), qtl, by = "id")
    
    return(qtl)
    
}

get_ld = function(gene, ref)
{
    analysis = ifelse(gene %like% "ENST", "use_isoform", "tpm_gene")
    
    dir = paste("pipeline/1.3.genotype", analysis, sep = "/")
    files = list.files(dir)
    file = paste(dir, files[which(files %like% gene & files %like% "gt_data")], sep = "/")
    
    message(paste("Opening", file))
    gt_data = add_rownames(fread(file, data.table = F))
    other = rownames(gt_data)[which(!rownames(gt_data) == ref)]
    
    message("Computing LD")
    out = melt(cor(t(gt_data))) %>% dplyr::mutate(value = ifelse(is.na(value), 0, abs(value)) )
    colnames(out) = c("snp1", "id", "LD")
    out = out %>% dplyr::filter(snp1 == ref) %>% dplyr::select(id, LD)
    
    message("Done.")
    return(out)
    
}

plot_isoform_structure = function(genes, start, end)
{
    exoninfo = fread("/reference/private/Gencode.v34lift37/exon_info.txt", data.table = F) %>% 
        mutate(start = start / 1e6, 
               end = end / 1e6, 
               transcript_id = simplify_id(transcript_id)) 

    exoninfo = exoninfo %>% filter(transcript_id %in% genes)

    exoninfo = exoninfo[exoninfo$start >= start& exoninfo$end <= end,]

    g = ggplot() + jn_theme

    # Add exons
    for (row in c(1:nrow(exoninfo)))
    {
        this = exoninfo[row,]
        g = g + annotate("segment", x = this$start, xend = this$end, y = this$transcript_id, yend = this$transcript_id, lwd = 15)
    }

    # Connect exons
    for (transcript in unique(exoninfo$transcript_id))
    {
        this = exoninfo[exoninfo$transcript_id == transcript,]
        #min = max(start, min(this$start))
        #max = min(end, max(this$end))
        min = min(this$start)
        max = max(this$end)

        g = g + annotate("segment", x = min , xend = max, y = transcript, yend = transcript, size = 1)
    }

    options(repr.plot.height = 8, repr.plot.width = 12)

    g = g  + xlim(start, end) + theme(panel.border = element_blank(), axis.ticks = element_blank(), axis.line.x = element_line(size = 0.8)) + 
        xlab("Position (Mb)") + ylab("")

    return(g)
}
    
    

simplify_id = function(vector)
{
    unlist(lapply(vector, function(x) { unlist(strsplit(x, "[.]"))[1] }))
}

plot_signal = function(data, col, ylab, title)
{
    p = ggplot(data, aes(x = pos / 1e6, y = -log10(pval), color = LD)) + geom_point(size = 4) + 
        jn_theme + 
        scale_color_gradientn(colors = c("lightgrey", col)) + 
        geom_point(data = data %>% filter(LD >= 0.7), aes(x = pos / 1e6, y = -log10(pval), color = LD), size = 4) + 
        geom_point(data = data, aes(x = pos / 1e6, y = -log10(pval)), shape = 1, color = "black", stroke = 1.2, size = 4) + 
        xlab("Position (Mb)") + ylab(ylab) + geom_hline(yintercept = -log10(0.05 / nrow(data)), linetype = "dashed") + 
        geom_vline(xintercept = as.numeric(unlist(strsplit(as.character(data[data$lead == T,]$id), "_"))[3]) / 1e6 , linetype = "dashed")
    
    if (title != "NULL") { p = p + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))}
    
    return(p)
}

compare_pval = function(data1, data2, xlab, ylab, coloc_lead)
{
    data = merge(data1, data2, by = "pos")
    p = ggplot(data, aes(x = -log10(pval.x), y = -log10(pval.y)) ) + geom_point(size = 3) +
        jn_theme  + xlab(xlab) + ylab(ylab)
    
    if (coloc_lead == "NULL")
    {
        p = p + geom_point(data = data %>% filter(lead.x == T) %>% arrange(desc(pval.y)) %>% head(1), aes(x = -log10(pval.x), y = -log10(pval.y)), color = "red", size = 5) + 
                geom_point(data = data %>% filter(lead.y == T) %>% arrange(desc(pval.x)) %>% head(1), aes(x = -log10(pval.x), y = -log10(pval.y)), color = "red", size = 5) 
    } else 
    {
        p = p + geom_point(data = data %>% filter(pos == coloc_lead) %>% arrange(pval.x, pval.y) %>% head(1),aes(x = -log10(pval.x), y = -log10(pval.y)), color = "red", size = 5)
    }
    
    return(p)
}

plot_coloc = function(hyp, data, isof_gene)
{
    # Plot
    col_pal = list("#ffbf00", "#04646a", "#e06666", "#6c4681", "#567b40")
    names(col_pal) = paste0("PP.H", c(0:4), ".abf")
    hyps = paste0("PP.H", c(0:4), ".abf")
    hyp_order = hyps
    
    data = data %>% filter(likely_hyp %like% hyp)
    
    this          = unique(melt(data %>% select(paste0("PP.H", c(0:4), ".abf"), pair)))
    this          = this %>% arrange(desc(value))
    
    gene_order    = this %>% filter(variable %like% hyp) %>% arrange(desc(value)) %>% select(pair) %>% unlist()

    this$variable = factor(this$variable, levels = rev(c(hyp, hyps[which(hyps != hyp)])) )
    this$pair     = factor(this$pair, levels = gene_order)

    p = ggplot(this, aes(x = pair, y = value, fill = variable)) + geom_bar(stat = "identity", width = 1.01) + 
                jn_theme + 
                theme(axis.text.x = element_blank(), 
                      axis.ticks.x = element_blank(),
                      legend.text = element_text(size = 25),
                      panel.border = element_blank(), 
                      axis.line.x = element_line(size = 1.2),
                      plot.title = element_text(hjust = 0.5, size = 29), 
                      axis.title = element_text(size = 25, vjust = -1)
                     ) + 
                scale_fill_manual(values = col_pal[hyp_order], name = "") + ylab("PPA") + 
                guides(fill = guide_legend(override.aes = list(size = 10))) + 
                scale_y_continuous(expand = c(0,0)) + 
                xlab(paste(length(unique(this$pair)), "pairs")) + 
                geom_hline(yintercept = 0.8, linetype = "dashed", size = 1.2) + 
                ylim(0,1.01)
    
    if (hyp != "PP.H0.abf") { p = p + theme(legend.position = "none", plot.margin = unit(c(0.3, 0.2, 0.5, 0.2), "cm")) }
    if (hyp != "PP.H4.abf") { p = p + theme(plot.margin = unit(c(0.3, 0.2, 0.5, -1.5), "cm"), axis.ticks.y = element_blank(), axis.text.y = element_blank() ) + ylab("") }
    if (hyp == "PP.H4.abf") { p = p + theme(axis.line.y = element_line(size = 1.2), plot.margin = unit(c(0.3, 0.2, 0.5, 0.5), "cm")) }
    
    if (isof_gene == T) { p = p + xlab(paste(length(unique(data$transcript_id.2)), "eIsoforms\n", length(unique(data$transcript_id.1)), "eGenes"))  }

    return(p)    
}


jn_theme = theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 20)) + theme(plot.margin = unit(c(1,1,1,3), "lines"),  title = element_text(size = 20), panel.background = element_blank(), axis.ticks.length=unit(0.50, "cm"), axis.text = element_text(color = "black"), axis.title.x = element_text(vjust = -1), axis.title.y = element_text(vjust = 3), axis.ticks = element_line(size = 1.25), panel.border = element_rect(color = "black", size = 2, fill = NA))

add_rownames = function(x) # add rownames to fread
{
	rownames(x) = x[,1]
	x[,1]       = NULL
	return(x)
}

### Gene expression normalization
my.invnorm = function(x)
{
    res = rank(x)
    res = qnorm(res/(length(res)+0.5))
    return(res)
}

transform_standard_normal = function(df)
{
    df                           = as.matrix(df)
    data_valid_expressed_full_qn = normalize.quantiles(df, copy=FALSE)
    input_mat                    = as.data.frame(t(apply(t(data_valid_expressed_full_qn), 2, my.invnorm)))
    
    return(input_mat)
}

normalize_tpm = function(outfolder, tpm, geneinfo, normalize = TRUE, gene_ids = NULL, filter_exp = FALSE, min_tpm = 2, min_samples = 0.1, use = FALSE)
{
    if(is.null(gene_ids) == FALSE){tpm = tpm[gene_ids,]}
    message("Normalizing expression...")
    message(paste("Total genes/peaks", nrow(tpm), sep = " = "))
    message(paste("Total samples"    , ncol(tpm), sep = " = "))
    
    expressed = as.matrix(tpm)
    
    if(filter_exp == TRUE)
    {
        expressed[as.matrix(tpm) <  min_tpm] = 0
        expressed[as.matrix(tpm) >= min_tpm] = 1

        tpm_f = tpm[rowSums(expressed) >= (min_samples * ncol(tpm)),]
        
    }else
    {
        tpm_f = tpm
    }
    
    if(use == TRUE)
    {
        geneinfo   = geneinfo[geneinfo$transcript_id %in% rownames(tpm_f),]
		gene_table = table(geneinfo$gene_id)
		geneinfo   = geneinfo[geneinfo$gene_id %in% names(gene_table[gene_table > 1]),]
		
        tpm_f = as.data.frame(rbindlist(lapply(sort(unique(geneinfo$gene_id)), function(gene_id)
        {
            this_transcript_ids = geneinfo[geneinfo$gene_id == gene_id, "transcript_id"]
            this                = tpm_f[this_transcript_ids, ]
            mysums              = unlist(lapply(colSums(this), function(x){max(c(100, x))}))
            out                 = as.data.frame(100 * t(t(as.matrix(this)) / mysums))
            out$transcript_id   = this_transcript_ids
            
            return(out)
        })), stringsAsFactors = FALSE)
        
        rownames(tpm_f)     = tpm_f$transcript_id
        tpm_f$transcript_id = NULL
    }
    
    message(paste("Total expressed genes"    , nrow(tpm_f), sep = " = "))
    
    if(normalize == TRUE)
    {
        tpm_f_std_norm = transform_standard_normal(tpm_f)
    }else
    {
        tpm_f_std_norm = tpm_f
    }
    
    fwrite    (tpm_f          , file = paste(outfolder, "expressed" , "txt", sep = "."), sep = "\t", row.names = TRUE , col.names = TRUE)
    fwrite    (tpm_f_std_norm , file = paste(outfolder, "normalized", "txt", sep = "."), sep = "\t", row.names = TRUE , col.names = TRUE)
    writeLines(rownames(tpm_f), con  = paste(outfolder, "gene_ids"  , "txt", sep = "."))
    
    return(list(tpm_f = tpm_f, tpm_f_std_norm = tpm_f_std_norm, gene_ids = rownames(tpm_f), sample_ids = colnames(tpm_f)))
}

calculate_peer_factors = function(outfolder, expdata, peer_factor_n = 10, n_genes = 0)
{
    tpm_f_std_norm    = expdata$tpm_f_std_norm
    tpm_f             = expdata$tpm_f
    
    if (n_genes > 0)
    {
        totest         = data.frame(gene_id = rownames(tpm_f), sd = unlist(apply(tpm_f, 1, sd)))
        totest         = totest[order(totest$sd, decreasing = TRUE), "gene_id"]
        tpm_f_std_norm = tpm_f_std_norm[totest[1:n_genes],]
    }
   
    model = PEER()
    PEER_setPhenoMean (model, t(as.matrix(tpm_f_std_norm)))
    PEER_setNk        (model, peer_factor_n               )
    PEER_update       (model                              )

    factors   = PEER_getX        (model)
    weights   = PEER_getW        (model)
    precision = PEER_getAlpha    (model)
    residuals = PEER_getResiduals(model)
    precision = data.frame(peer = paste("peer", 1:peer_factor_n, sep = "") , precision = precision, stringsAsFactors = FALSE)
    residuals = t(residuals)

    rownames(factors  ) = colnames(tpm_f_std_norm)
    rownames(weights  ) = rownames(tpm_f_std_norm)
    rownames(residuals) = rownames(tpm_f_std_norm)
    colnames(residuals) = colnames(tpm_f_std_norm)
    colnames(factors  ) = paste("peer", 1:peer_factor_n, sep = "")
    colnames(weights  ) = paste("peer", 1:peer_factor_n, sep = "")

    factors          = as.data.frame(factors  )
    weights          = as.data.frame(weights  )
    residuals        = as.data.frame(residuals)
    factors$assay_id = rownames     (factors  )
    
    fwrite(tpm_f_std_norm, file = paste(outfolder, "peer", "genes", "txt", sep = "."), sep = "\t", row.names = TRUE, col.names = TRUE)
    fwrite(factors       , file = paste(outfolder, "peer", "factors"  , "txt", sep = "."), sep = "\t", row.names = TRUE , col.names = TRUE)
    fwrite(weights       , file = paste(outfolder, "peer", "weights"  , "txt", sep = "."), sep = "\t", row.names = TRUE , col.names = TRUE)
    fwrite(precision     , file = paste(outfolder, "peer", "precision", "txt", sep = "."), sep = "\t", row.names = TRUE , col.names = TRUE)
    fwrite(residuals     , file = paste(outfolder, "peer", "residuals", "txt", sep = "."), sep = "\t", row.names = TRUE , col.names = TRUE)
    
    return(factors)
}

divide_phenotypes_by_gene = function(expdata, outfolder)
{
    message("Dividing phenotype data by gene/peak...")
    
    dir.create(outfolder, showWarnings = FALSE)

    gene_ids   = expdata$gene_ids
    sample_ids = expdata$sample_ids
    rawexp     = expdata$tpm_f
    normexp    = expdata$tpm_f_std_norm
    
    invisible(lapply(gene_ids, function(gene)
    {
        outdata = data.frame(sample_id = sample_ids, raw = as.numeric(rawexp[gene,sample_ids]), norm = as.numeric(normexp[gene,sample_ids]))
        fwrite(outdata, paste(outfolder, paste(gene, "txt", sep = "."), sep = "/"), sep = "\t", row.names = FALSE, col.names = TRUE)
    }))
}

divide_phenotypes_tissue = function(outfolder, tpm, geneinfo, run_peer = FALSE, normalize = FALSE, gene_ids = NULL, filter_exp = FALSE, phenotype_min_value = 1, phenotype_min_samples = 0.1, peer_factor_n = 10, n_genes = 1000, use = FALSE)
{
    message("---------------------------")
    
    expdata = normalize_tpm(outfolder, tpm, geneinfo, normalize, gene_ids, filter_exp, min_tpm = phenotype_min_value, min_samples = phenotype_min_samples, use)
    
    if(run_peer == TRUE){message("Computing PEER..."); peerdata = calculate_peer_factors(outfolder, expdata, peer_factor_n, n_genes)}
    
    divide_phenotypes_by_gene(expdata, outfolder)
}

# eQTL analysis
write_h5_file = function(h5_file, analysis, gtinfo, gtdata, gene_id, expdata)
{
    invisible(suppressWarnings(file.remove(h5_file)))

    invisible(h5createFile (h5_file))
    invisible(h5createGroup(h5_file, "genotype"))
    invisible(h5createGroup(h5_file, "genotype/col_header"))
    invisible(h5createGroup(h5_file, "genotype/row_header"))
    invisible(h5createGroup(h5_file, "phenotype"))
    invisible(h5createGroup(h5_file, "phenotype/col_header"))
    invisible(h5createGroup(h5_file, "phenotype/row_header"))
    
    genedata = geneinfo[geneinfo$transcript_id == gene_id,]

    h5write(gtinfo$chrom           , file = h5_file, name="genotype/col_header/chrom"  )
    h5write(gtinfo$pos             , file = h5_file, name="genotype/col_header/pos"    )
    h5write(gtinfo$pos             , file = h5_file, name="genotype/col_header/pos_cum")
    h5write(as.matrix(gtdata)      , file = h5_file, name="genotype/matrix"            )
    h5write(colnames (gtdata)      , file = h5_file, name="genotype/row_header/sample_ID")
    h5write(c(gene_id)             , file = h5_file, name="phenotype/col_header/gene_ID"     )
    h5write(genedata[, "chrom"    ], file = h5_file, name="phenotype/col_header/gene_chrom"  )
    h5write(genedata[, "end"      ], file = h5_file, name="phenotype/col_header/gene_end"    )
    h5write(genedata[, "start"    ], file = h5_file, name="phenotype/col_header/gene_start"  )
    h5write(genedata[, "strand"   ], file = h5_file, name="phenotype/col_header/gene_strand" )
    h5write(genedata[, "gene_name"], file = h5_file, name="phenotype/col_header/phenotype_ID")
    h5write(as.matrix(expdata)     , file = h5_file, name="phenotype/matrix"                 )
    h5write(colnames( expdata)     , file = h5_file, name="phenotype/row_header/sample_ID")
}

run_eigenmt = function(gt_file, geneloc_file, snploc_file, qtl_file, fdr_file, chrom)
{
    command = paste("python", paste(getwd(),"scripts/eigenMT.py", sep = "/"), 
                    "--CHROM"   , chrom,
                     "--QTL"    , qtl_file,
                     "--GEN"    , gt_file,
                     "--GENPOS" , snploc_file,
                     "--PHEPOS" , geneloc_file,
                     "--OUT"    , fdr_file
                   )
    
	message(command)			   
    system(command)
}


run_eqtl = function(gene_id, analysis, tmp_folder, expdata, gtinfo, gtdata, h5_file, cov_file, kin_file, qtl_file_tmp, geneloc_file, snploc_file, gt_file, chrom)
{
    suppressWarnings(write_h5_file(h5_file, analysis, gtinfo, gtdata, gene_id, expdata))
	
    command_py = paste("python", paste(getwd(), "scripts", "run_limix.py", sep = "/"), tmp_folder, gene_id, "scan", "normal")
    
    system(command_py)
    
    indata                   = add_rownames(fread(qtl_file_tmp, sep = ",", header = TRUE, data.table = FALSE))
	colnames(indata)         = c("beta", "se", "pval")
    indata                   = cbind(gtinfo, indata)
    indata$gene_id           = gene_id
    indata$bonferroni        = p.adjust(indata$pval, method = "bonferroni")
	
	if(nrow(indata[indata$se > 100, ]) > 0){indata[indata$se > 100, "se"] = 100}
	
    if(min(indata$bonferroni) == 1)
    {
        fdrdata       = indata[which.min(indata$pval), ]
        fdrdata$fdr   = 1
        fdrdata$tests = nrow(indata)
    }else
    {
        eigenmt_input            = indata[,c("id", "gene_id", "beta", "se", "pval", "bonferroni")]
        colnames(eigenmt_input)  = c("SNP", "gene", "beta", "t-stat", "p-value", "FDR")
        qtl_file_emt_in          = sub("qtl.csv", "eigenmt_input.txt" , qtl_file_tmp)
        qtl_file_emt_out         = sub("qtl.csv", "eigenmt_output.txt", qtl_file_tmp)
	
        fwrite(eigenmt_input, qtl_file_emt_in , sep = "\t", row.names = FALSE, col.names = TRUE)
	
        run_eigenmt(gt_file, geneloc_file, snploc_file, qtl_file_emt_in, qtl_file_emt_out, chrom)
	
        fdrdata           = fread(qtl_file_emt_out, sep = "\t", header = TRUE, data.table = FALSE)
        colnames(fdrdata) = c("id", "gene_id", "beta", "se", "pval", "bonferroni", "fdr", "tests")
    }
    
    fdrdata           = merge(gtinfo, fdrdata)
    fdrdata           = fdrdata[,c(colnames(indata), "fdr", "tests")]
    return(list(lead = fdrdata, qtl = indata))
	return(indata)
}



# Find interactions
find_interactions = function(totest, var1, gene_id, analysis, tmp_folder, expdata, gtinfo, gtdata, h5_file, cov_file, kin_file, qtl_file_tmp, geneloc_file, snploc_file, gt_file, chrom)
{
    suppressWarnings(write_h5_file(h5_file, analysis, gtinfo, gtdata, gene_id, expdata))
	
    command_py = paste("python", paste(getwd(), "script", "run_limix.py", sep = "/"), tmp_folder, gene_id, var1, "normal")
    
    system(command_py)
    
    indata             = add_rownames(fread(qtl_file_tmp, sep = ",", header = TRUE, data.table = FALSE))
	colnames(indata)   = c("beta", "se", "pval")
	indata             = cbind(totest, indata)
	indata$interaction = var1
	
	return(indata)
}

# Merge eQTLs

merge_qtls = function(infolder, geneinfo)
{
    infiles = list.files(paste("pipeline/3.2.eqtls/eqtls_by_gene", infolder, sep = "/"), pattern = "^fdr", full.names = TRUE)
    
    indata     = suppressWarnings(as.data.frame(rbindlist(lapply(infiles, function(x){fread(x, sep = "\t", header = TRUE, data.table = FALSE)})), stringsAsFactors = FALSE))
    outdata    = list()
    gene_ids   = unique(indata$gene_id)
    
    for(type in sort(unique(indata$type)))
    {
        if(length(gene_ids) > 0)
        {
            this = indata[indata$type == type & indata$gene_id %in% gene_ids,]

            #this$qval  = qvalue(this$fdr)$qvalues
            this$qval  = p.adjust(this$fdr, method = "BH")
            this$egene = FALSE

            this[this$qval <= 0.05, "egene"] = TRUE

            gene_ids            = this[this$egene == TRUE, "gene_id"]
            outdata[[type + 1]] = this
        }else
        {
            break
        }
    }
    
    out          = as.data.frame(rbindlist(outdata), stringsAsFactors = FALSE)
    out          = merge(geneinfo[,c("transcript_id", "gene_id", "gene_name", "gene_type", "start", "end", "strand")], out, by.x = "transcript_id", by.y = "gene_id")
    out$distance = out$pos - out$start

    out[out$strand == "-", "distance"] = out[out$strand == "-", "end"] - out[out$strand == "-", "pos"]
    
    fwrite(out, paste("pipeline/3.2.eqtls/eqtls", paste(infolder, "egenes", "txt", sep = "."), sep = "/"), sep = "\t", col.names = TRUE, row.names = FALSE)
    
    return(out)
}

colPanel = c(
		"grey", "#E31A1C", "#FFD700", "#771122", "#777711", "#1F78B4", "#68228B", "#AAAA44",
		"#60CC52", "#771155", "#DDDD77", "#774411", "#AA7744", "#AA4455", "#117744", 
		"#000080", "#44AA77", "#AA4488", "#DDAA77", "#D9D9D9", "#BC80BD", "#FFED6F",
	    "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17",
	    "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
	    "#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
	    "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928", "#FBB4AE", "#B3CDE3",
	    "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2",
	    "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC",
	    "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FFFF33", "#A65628",
	    "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
	    "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
	    "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5"
	   )


