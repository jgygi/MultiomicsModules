#' Get multi-omics module scores for a list of analytes.
#' @param name What version of `scores_df` to use (see `ModObj$scores$...`)? Added via `$add_analyte_scores(...)`. Defaults to `"table1"`.
#' @param method What method to derive scores? Defaults to `"fisher"`.
#' @param p.adjust Method to adjust p.values? Defaults to `"BH"`. Can also be `"none"`.
#' @param min.analytes Minimum number of analytes needed to calculate a score (per assay). Defaults to `0`.
#' @examples
#' get_module_scores(...)
#' 
#'@export
get_module_scores = function(name = NULL, method = "fisher", p.adjust = "BH", min.analytes = 0){
  # check name:
  if(is.null(name)){
    name = names(self$scores)
  }
  
  # Get scores for all names:
  all.p.vals <- list()
  mod.df <- self$modules
  for(cur.name in name){
    if(!self$params$quiet){cat("Calculating scores for name =", cur.name, "...\n")}
    tmp.assay.res <- list()
    for(assay in c("Transcriptomics", "Proteomics", "Metabolomics", "Total")){
      tmp.p.val.df <- data.frame()
      for(mod.id in unique(mod.df$ModuleID)){
        cur.mod.df <- dplyr::filter(mod.df, ModuleID == mod.id)
        if(assay != "Total"){
          cur.mod.df <- dplyr::filter(cur.mod.df, Assay == assay)
        }
        rownames(cur.mod.df) <- cur.mod.df$AnalyteID
        if(nrow(cur.mod.df) == 0){
          # No analytes for this module...
          next
        }
        cur.mod.df$P.value <- NA
        # Get relative scores (for module):
        p.val.scores <- sapply(1:nrow(self$scores[[cur.name]]), function(i){
          assay = self$scores[[cur.name]]$Assay[i]
          analyte.id = self$scores[[cur.name]]$ConvertedID[i]
          direction = self$scores[[cur.name]]$Direction[i]
          p.value = self$scores[[cur.name]]$P.value[i]
          if(!is.na(analyte.id)){
            if(analyte.id %in% cur.mod.df$AnalyteID){
              return(ifelse(direction == "Pos", p.value, -p.value))
            } else {
              return(NA)
            }
          } else {
            return(NA)
          }
        })
        names(p.val.scores) <- self$scores[[cur.name]]$ConvertedID
        p.val.scores <- p.val.scores[!is.na(p.val.scores)]
        cur.mod.df[names(p.val.scores),"P.value"] <- p.val.scores
        found.cur.mod.df <- dplyr::filter(cur.mod.df, !is.na(P.value))
        found.cur.mod.df$FoundDir <- ifelse(found.cur.mod.df$P.value >= 0, "Pos", "Neg")
        found.cur.mod.df$P.value <- abs(found.cur.mod.df$P.value)
        
        # NOTE: change here for different types of combinations:
        # also: add check for min num analytes
        tmp <- dplyr::filter(found.cur.mod.df, Direction == FoundDir)
        reg.p.val <- poolr::fisher(tmp$P.value)$p
        tmp <- dplyr::filter(found.cur.mod.df, Direction != FoundDir)
        inv.p.val <- poolr::fisher(tmp$P.value)$p
        tmp <- found.cur.mod.df
        comb.p.val <- poolr::fisher(tmp$P.value)$p
        
        # add to results
        tmp.p.val.df <- rbind(tmp.p.val.df, data.frame(
          Pos_P.val = reg.p.val,
          Neg_P.val = inv.p.val,
          All_P.val = comb.p.val, row.names = mod.id
        )) 
      }
      res = as.matrix(tmp.p.val.df)
      if(p.adjust != "none"){
        tmp.res <- matrix(p.adjust(res, method = "BH"), ncol = 3)
        rownames(tmp.res) <- rownames(res)
        colnames(tmp.res) <- colnames(res)
        res <- tmp.res
      }
      tmp.assay.res[[assay]] <- res
    }
    self$reports$scores[[cur.name]] <- tmp.assay.res
    all.p.vals[[cur.name]] <- tmp.assay.res
  }
  return(all.p.vals)
}



#' Plot the representative analytes of each assay (Proteomics, Metabolomics, Transcriptomics) vs scores
#' @param name What version of `scores_df` to use (see `ModObj$scores$...`)? Added via `$add_analyte_scores(...)`. Defaults to all (`NULL`).
#' @param assay What assay to show? Can be `"total"` (default) or `"proteomics"`, `"metabolomics"`, or `"transcriptomics"`.
#' @examples
#' # run first:
#' ModObj$get_module_scores(...)
#' ModObj$plot_module_scores()
#' 
#'@export
plot_module_scores = function(name = NULL, assay = "Total"){
  if(is.null(name)){
    name = names(self$scores)
  }
  if(length(name) > 1){
    # Multi plot:
    all.df <- NULL
    for(cur.name in name){
      df <- self$reports$scores[[cur.name]][[assay]]
      if(is.null(all.df)){
        all.df <- df[,3,drop=F]
      } else {
        all.df <- cbind(all.df, df[,3])
      }
    }
    colnames(all.df) <- name
    df.melt <- reshape2::melt(as.matrix(all.df))
    
    p <- ggplot(df.melt) +
      geom_point(aes(x = Var2, y = Var1, size = -log10(value), color = Var2)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = .5))
    
    return(p)
  } else {
    # Single plot:
    if(is.null(self$reports$scores[[name]])){
      stop("ERROR: No scores found for name: ", name, ". Please run $get_module_scores(...) first...")
    }
    df <- self$reports$scores[[name]][[assay]]
    df.melt <- reshape2::melt(df)
    
    ggplot(df.melt) +
      geom_point(aes(x = Var2, y = Var1, size = -log10(value), color = Var2)) +
      scale_color_manual(values = c(
        Pos_P.val = "red",
        Neg_P.val = "blue",
        All_P.val = "grey"
      ))
  }
}



#' Add analyte scores to MOModules Object.
#' @param scores_df Data.frame of scores including the following columns: AnalyteID (name of analyte), AnalyteIDtype (type of analyte for processing), Assay (one of `Proteomics`, `Metabolomics` or `Transcriptomics`), Score (a numeric), Direction (either `Pos` for Positive or `Neg` for Negative), and Significant (boolean, `TRUE` or `FALSE`)
#' @param name What to store the `scores_df` under (see `ModObj$scores$...`)? Defaults to `table1` for first scores, `table2`, and so on.
#' @param map Try to map the data to the multi-omics modules? Defaults to `TRUE`
#' @examples
#' ModObj <- make.multiomics.modules(...)
#' 
#' ModObj$get_module_scores(...)
#' 
#'@export
add_analyte_scores = function(scores_df = NULL, name = NULL, map = TRUE){
  # Add name (if null):
  if(is.null(name)){
    names.with.tab <- grep("table", names(self$scores))
    if(length(names.with.tab) > 0){
      tab.num = max(as.integer(gsub("table", "", names(self$scores)[grepl("table", names(self$scores))]))) + 1
    } else {
      tab.num = 1
    }
    name = paste0("table", tab.num)
  }
  
  if(is.null(scores_df)){stop("ERROR: Please supply the 'scores_df' parameter appropriately (see documentation).")}
  # check if columns are correct
  # 7 columns:
  if(class(scores_df) != "data.frame"){
    scores_df <- as.data.frame(scores_df)
    warning("WARNING: 'scores_df' supplied is not of class 'data.frame'. Coercing to data.frame...")
  }
  if(ncol(scores_df) != 7){stop("ERROR: `scores_df` does not have the required 7 columns (AnalyteID, AnalyteIDtype, Assay, Score, Direction, P.value, and Significant). All are required.")}
  # quickly check four columns:
  # column 2: must be one of:
  if(any(!scores_df[,2] %in% unique(unlist(sapply(self$ID_tables, colnames))))){
    stop(paste0("ERROR: One or more 'AnalyteIDtype' are not allowed. Must be exactly one of the following: '", paste(unique(unlist(sapply(self$ID_tables, colnames))), collapse = "', '"), "'"))
  }
  # column 3: must be one of Proteomics, Transcriptomics, or Metabolomics
  if(any(!scores_df[,3] %in% c("Proteomics", "Transcriptomics", "Metabolomics"))){
    stop("ERROR: One or more 'Assay' are not allowed. Must all be one of: 'Transcriptomics', 'Proteomics', 'Metabolomics'.")
  }
  # column 5: must be one of Pos, Neg
  if(any(!scores_df[,5] %in% c("Pos", "Neg"))){
    stop("ERROR: One or more 'Direction' are not allowed. Must all be one of: 'Pos', 'Neg'.")
  }
  # column 7: must be one of TRUE, FALSE
  if(any(!scores_df[,7] %in% c(TRUE, FALSE))){
    stop("ERROR: One or more 'Significant' are not allowed. Must all be one of: 'TRUE', 'FALSE'.")
  }
  # fix colnames:
  colnames(scores_df) <- c("AnalyteID", "AnalyteIDtype", "Assay", "Score", "Direction", "P.value", "Significant")
  
  # Map to data?
  if(map){
    mapped.ids <- sapply(1:nrow(scores_df), function(i){
      idx = which(self$ID_tables[[scores_df$Assay[i]]][, scores_df$AnalyteIDtype[i]] == scores_df$AnalyteID[i])
      if(length(idx) == 0){
        return(NA)
      } else {
        return(idx)
      }
    })
    scores_df$MappedIDX = mapped.ids
    converted.ids <- sapply(1:nrow(scores_df), function(i){
      if(is.na(scores_df$MappedIDX[i])){return(NA)}
      if(scores_df$Assay[i] == "Transcriptomics"){
        return(self$ID_tables$Transcriptomics$ENSEMBL[scores_df$MappedIDX[i]])
      } else if(scores_df$Assay[i] == "Metabolomics"){
        return(self$ID_tables$Metabolomics$IMPACC_ID[scores_df$MappedIDX[i]])
      } else if(scores_df$Assay[i] == "Proteomics"){
        return(self$ID_tables$Proteomics$GENE_SYMBOL[scores_df$MappedIDX[i]])
      }
    })
    scores_df$ConvertedID = converted.ids
    # Summarize:
    total.genes = sum(scores_df$Assay == "Transcriptomics")
    if(total.genes > 0){
      mapped.genes = sum(!is.na(scores_df[scores_df$Assay == "Transcriptomics",]$MappedIDX))
    } else {mapped.genes = 0}
    total.mets = sum(scores_df$Assay == "Metabolomics")
    if(total.mets > 0){
      mapped.mets = sum(!is.na(scores_df[scores_df$Assay == "Metabolomics",]$MappedIDX))
    } else {mapped.mets = 0}
    total.prots = sum(scores_df$Assay == "Proteomics")
    if(total.prots > 0){
      mapped.prots = sum(!is.na(scores_df[scores_df$Assay == "Proteomics",]$MappedIDX))
    } else {mapped.prots = 0}
    
    map.df <- data.frame(
      MappedGenes = mapped.genes,
      TotalGenes = total.genes,
      PropGenes = mapped.genes/total.genes,
      MappedProts = mapped.prots,
      TotalProts = total.prots,
      PropProts = mapped.prots/total.prots,
      MappedMets = mapped.mets,
      TotalMets = total.mets,
      PropMets = mapped.mets/total.mets
    )
    
    self$reports$map[[name]] <- map.df
    
  } else {
    scores_df$MappedIDX = NA
    self$reports$map[[name]] <- data.frame(
      MappedGenes = NA,
      TotalGenes = NA,
      PropGenes = NA,
      MappedProts = NA,
      TotalProts = NA,
      PropProts = NA,
      MappedMets = NA,
      TotalMets = NA,
      PropMets = NA
    )
  }
  
  # add scores:
  self$scores[[name]] <- scores_df
  if(!self$params$quiet){
    cat(paste0("Scores successfully added under $scores$", name, "\n"), sep =)
    if(map){
      cat("Analytes:    \tMapped:\tPercentage:\n")
      cat("Genes       |\t", mapped.genes, "/", total.genes, "\t(", round(100*mapped.genes/total.genes, 2), "%)\n", sep = "")
      cat("Proteins    |\t", mapped.prots, "/", total.prots, "\t(", round(100*mapped.prots/total.prots, 2), "%)\n", sep = "")
      cat("Metabolites |\t", mapped.mets, "/", total.mets, "\t(", round(100*mapped.mets/total.mets, 2), "%)\n", sep = "")
    }
  }
  
  return(invisible(self))
}


#' Plot the overlap per analyte of each assay (Proteomics, Metabolomics, Transcriptomics)
#' @param name What version of `scores_df` to use (see `ModObj$scores$...`)? Added via `$add_analyte_scores(...)`. Defaults to `"table1"`.
#' @examples
#' ModObj$plot_mapping()
#' 
#' ModObj$plot_mapping(name = "my_scores")
#' 
#'@export
plot_mapping = function(name = "table1"){
  tmp <- self$modules
  tmp$Mapped <- "Missing"
  tmp$ScoreDirection <- NA
  for(i in 1:nrow(self$scores[[name]])){
    if(!is.na(self$scores[[name]]$MappedIDX[i])){
      assay = self$scores[[name]]$Assay[i]
      analyte.id = self$scores[[name]]$AnalyteID[i]
      analyte.id.type = self$scores[[name]]$AnalyteIDtype[i]
      tmp$Mapped[which(tmp$Assay == assay & tmp$AnalyteID == self$scores[[name]]$ConvertedID[i])] <- assay
      tmp$ScoreDirection[which(tmp$Assay == assay & tmp$AnalyteID == self$scores[[name]]$ConvertedID[i])] <- self$scores[[name]]$Direction[i]
    }
  }
  tmp$MatchesDirection = tmp$Direction == tmp$ScoreDirection
  
  p <- ggplot2::ggplot(tmp) +
    ggplot2::geom_bar(ggplot2::aes(x = factor(ModuleID, levels = paste0("Mod", 1:100)), 
                                   fill = factor(Mapped, levels = c("Missing", "Metabolomics", "Proteomics", "Transcriptomics"))), 
                      stat = "count", position = "fill") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = .5),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "#EDEDED")) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("% Mapped") +
    ggplot2::scale_fill_manual(values = c(Missing = "#D5E1E8", self$params$color_scheme)) +
    ggplot2::guides(fill = "none") +
    ggplot2::facet_wrap(ggplot2::vars(factor(Assay, levels = c("Transcriptomics", "Metabolomics", "Proteomics"))), ncol = 1) +
    ggplot2::scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0%", "25%", "50%", "75%", "100%"))
  
  return(p)
}



#' Convert an analyte id to be used within the IMPACC scores:
#' @param analyte.id ID of analyte to convert.
#' @param analyte.id.type Type of ID (e.g. ENSEMBL, HMBD, KEGG, GENE_SYMBOL, ...)
#' @param assay Which assay? (Proteomics, Transcriptomics, or Metabolomics)
#' @examples
#' ModObj$convert_id("IL6", "GENE_SYMBOL", "Proteomics")
#'  
#'@export
convert_id <- function(analyte.id, analyte.id.type, assay){
  if(assay == "Proteomics"){
    return(self$ID_tables[[assay]]$GENE_SYMBOL[which(self$ID_tables[[assay]][[analyte.id.type]] == analyte.id)])
  } else if(assay == "Transcriptomics"){
    return(self$ID_tables[[assay]]$ENSEMBL[which(self$ID_tables[[assay]][[analyte.id.type]] == analyte.id)])
  } else if(assay == "Metabolomics"){
    return(self$ID_tables[[assay]]$IMPACC_ID[which(self$ID_tables[[assay]][[analyte.id.type]] == analyte.id)])
  }
}



#' Perform enrichment analysis using the determined method on the chosen scores data.frame.
#' @param name What version of `scores_df` to use (see `ModObj$scores$...`)? Added via `$add_analyte_scores(...)`. Defaults to `"table1"`.
#' @param method What method of enrichment to use? Defaults to `"match"`
#' @param return Return the table stored under `$reports$...`? Defaults to `TRUE`
#' @examples
#' ModObj$do_enrichment()
#' 
#'@export
do_enrichment = function(name = "table1", method = "match", return = TRUE){
  all.res.df <- data.frame()
  cur.sig.df <- dplyr::filter(self$scores[[name]], !is.na(MappedIDX), Significant)
  cur.nosig.df <- dplyr::filter(self$scores[[name]], !is.na(MappedIDX), !Significant)
  gene.background <- dplyr::filter(self$scores[[name]], Assay == "Transcriptomics", !is.na(MappedIDX))
  protein.background <- dplyr::filter(self$scores[[name]], Assay == "Proteomics", !is.na(MappedIDX))
  metabolite.background <- dplyr::filter(self$scores[[name]], Assay == "Metabolomics", !is.na(MappedIDX))
  
  if(method == "match"){
    for(mod.id in unique(self$modules$ModuleID)){
      cur.mod.df <- dplyr::filter(self$modules, ModuleID == mod.id)
      print(mod.id)
      # Genes
      if(nrow(gene.background) > 0){
        cur.mod.df.assay <- dplyr::filter(cur.mod.df, Assay == "Transcriptomics")
        mod.analytes <- cur.mod.df.assay$AnalyteID
        mod.dirs <- cur.mod.df.assay$Direction
        cur.sig.df.assay <- dplyr::filter(cur.sig.df, Assay == "Transcriptomics")
        sig.analytes <- cur.sig.df.assay$ConvertedID
        sig.dirs <- cur.sig.df.assay$Direction
        cur.nosig.df.assay <- dplyr::filter(cur.nosig.df, Assay == "Transcriptomics")
        nosig.analytes <- cur.nosig.df.assay$ConvertedID
        nosig.dirs <- cur.nosig.df.assay$Direction
        
        #all.res.df <- rbind(all.res.df, data.frame(
        #  Module = mod.id,
        #  Assay = "Transcriptomics",
        #  sig.in.pathway = length(intersect(sig.analytes, mod.analytes)),
        #  nosig.in.pathway = length(intersect(nosig.analytes, mod.analytes))
        #))
      }
      # Prots
      if(nrow(protein.background) > 0){
        cur.mod.df.assay <- dplyr::filter(cur.mod.df, Assay == "Proteomics")
        mod.analytes <- cur.mod.df.assay$AnalyteID
        cur.sig.df.assay <- dplyr::filter(cur.sig.df, Assay == "Proteomics")
        sig.analytes <- cur.sig.df.assay$ConvertedID
        if(length(intersect(mod.analytes, sig.analytes)) > 0){
          mod.dirs <- cur.mod.df.assay$Direction[which(mod.analytes %in% sig.analytes)]
          sig.dirs <- cur.sig.df.assay$Direction[which(sig.analytes %in% mod.analytes)]
          mod.match.ids <- mod.analytes[which(mod.analytes %in% sig.analytes)]
          sig.match.ids <- sig.analytes[which(sig.analytes %in% mod.analytes)]
          sig.to.mod <- sapply(sig.match.ids, function(sig.id){return(which(mod.match.ids %in% sig.id))})
          sig.match.ids <- sig.match.ids[sig.to.mod]
          sig.dirs <- sig.dirs[sig.to.mod]
          cur.nosig.df.assay <- dplyr::filter(cur.nosig.df, Assay == "Proteomics")
          nosig.analytes <- cur.nosig.df.assay$ConvertedID
          nosig.dirs <- cur.nosig.df.assay$Direction
          
          # Calculate:
          sig.in.pathway = length(intersect(sig.analytes, mod.analytes))
          nosig.in.pathway = length(intersect(nosig.analytes, mod.analytes))
          total.in.pathway = sig.in.pathway + nosig.in.pathway
          sig.ratio = sig.in.pathway/total.in.pathway
          match.dir = sum(sig.dirs == mod.dirs)
          dir.ratio = match.dir/sig.in.pathway
        } else {
          sig.in.pathway = NA
          nosig.in.pathway = NA
          sig.ratio = NA
          match.dir = NA
          dir.ratio = NA
        }
        
        
        all.res.df <- rbind(all.res.df, data.frame(
          Module = mod.id,
          Assay = "Proteomics",
          sig.in.pathway = sig.in.pathway,
          nosig.in.pathway = nosig.in.pathway,
          sig.ratio = sig.ratio,
          match.dir = match.dir,
          dir.ratio = dir.ratio
        ))
      }
      # Mets
      if(nrow(metabolite.background) > 0){
        cur.mod.df.assay <- dplyr::filter(cur.mod.df, Assay == "Metabolomics")
        mod.analytes <- cur.mod.df.assay$AnalyteID
        mod.dirs <- cur.mod.df.assay$Direction
        cur.sig.df.assay <- dplyr::filter(cur.sig.df, Assay == "Metabolomics")
        sig.analytes <- cur.sig.df.assay$ConvertedID
        sig.dirs <- cur.sig.df.assay$Direction
        cur.nosig.df.assay <- dplyr::filter(cur.nosig.df, Assay == "Metabolomics")
        nosig.analytes <- cur.nosig.df.assay$ConvertedID
        nosig.dirs <- cur.nosig.df.assay$Direction
        
        #all.res.df <- rbind(all.res.df, data.frame(
        #  Module = mod.id,
        #  Assay = "Metabolomics",
        #  sig.in.pathway = length(intersect(sig.analytes, mod.analytes)),
        #  nosig.in.pathway = length(intersect(nosig.analytes, mod.analytes))
        #))
      }
    }
  }
  if(method == "fisher"){
    for(mod.id in unique(self$modules$ModuleID)){
      cur.mod.df <- dplyr::filter(self$modules, ModuleID == mod.id)
      # Transcriptomics
      if(nrow(gene.background) > 0){
        cur.mod.df.assay <- dplyr::filter(cur.mod.df, Assay == "Transcriptomics")
        mod.analytes <- cur.mod.df.assay$AnalyteID
        cur.sig.df.assay <- dplyr::filter(cur.sig.df, Assay == "Transcriptomics")
        sig.analytes <- cur.sig.df.assay$ConvertedID
        cur.nosig.df.assay <- dplyr::filter(cur.nosig.df, Assay == "Transcriptomics")
        nosig.analytes <- cur.nosig.df.assay$ConvertedID
        f.sig.in.pathway = length(intersect(sig.analytes, mod.analytes))
        f.nosig.in.pathway = length(intersect(nosig.analytes, mod.analytes))
        f.sig.no.pathway = length(intersect(sig.analytes, gene.background$ConvertedID)) - f.sig.in.pathway
        f.nosig.no.pathway = length(intersect(nosig.analytes, gene.background$ConvertedID)) - f.nosig.in.pathway
        ct <- matrix(c(
          sp=f.sig.in.pathway,
          sn=f.sig.no.pathway,
          np=f.nosig.in.pathway,
          ns=f.nosig.no.pathway),2,2,dimnames=list(c("MOD","NOMOD"),c("SIG","NO SIG")))
        p.val.gene = fisher.test(ct ,alt="less")$p.value
        #cat(paste0("Testing Genes: SigPath = ", f.sig.in.pathway, ", SigNopath = ", f.sig.no.pathway, ", NosigPath = ", f.nosig.in.pathway, ", NosigNopath = ", f.nosig.no.pathway, " | p.val = ", p.val.gene, "\n"))
      } else {
        p.val.gene <- NA
      }
      # Proteomics
      if(nrow(protein.background) > 0){
        cur.mod.df.assay <- dplyr::filter(cur.mod.df, Assay == "Proteomics")
        mod.analytes <- cur.mod.df.assay$AnalyteID
        cur.sig.df.assay <- dplyr::filter(cur.sig.df, Assay == "Proteomics")
        sig.analytes <- cur.sig.df.assay$ConvertedID
        cur.nosig.df.assay <- dplyr::filter(cur.nosig.df, Assay == "Proteomics")
        nosig.analytes <- cur.nosig.df.assay$ConvertedID
        f.sig.in.pathway = length(intersect(sig.analytes, mod.analytes))
        f.nosig.in.pathway = length(intersect(nosig.analytes, mod.analytes))
        f.sig.no.pathway = length(intersect(sig.analytes, protein.background$ConvertedID)) - f.sig.in.pathway
        f.nosig.no.pathway = length(intersect(nosig.analytes, protein.background$ConvertedID)) - f.nosig.in.pathway
        ct <- matrix(c(
          sp=f.sig.in.pathway,
          sn=f.sig.no.pathway,
          np=f.nosig.in.pathway,
          ns=f.nosig.no.pathway),2,2,dimnames=list(c("MOD","NOMOD"),c("SIG","NO SIG")))
        p.val.prot = fisher.test(ct ,alt="less")$p.value
        #cat(paste0("Testing Prots: SigPath = ", f.sig.in.pathway, ", SigNopath = ", f.sig.no.pathway, ", NosigPath = ", f.nosig.in.pathway, ", NosigNopath = ", f.nosig.no.pathway, " | p.val = ", p.val.prot, "\n"))
      } else {
        p.val.prot = NA
      }
      # Metabolomics
      if(nrow(metabolite.background) > 0){
        cur.mod.df.assay <- dplyr::filter(cur.mod.df, Assay == "Metabolomics")
        mod.analytes <- cur.mod.df.assay$AnalyteID
        cur.sig.df.assay <- dplyr::filter(cur.sig.df, Assay == "Metabolomics")
        sig.analytes <- cur.sig.df.assay$ConvertedID
        cur.nosig.df.assay <- dplyr::filter(cur.nosig.df, Assay == "Metabolomics")
        nosig.analytes <- cur.nosig.df.assay$ConvertedID
        f.sig.in.pathway = length(intersect(sig.analytes, mod.analytes))
        f.nosig.in.pathway = length(intersect(nosig.analytes, mod.analytes))
        f.sig.no.pathway = length(intersect(sig.analytes, metabolite.background$ConvertedID)) - f.sig.in.pathway
        f.nosig.no.pathway = length(intersect(nosig.analytes, metabolite.background$ConvertedID)) - f.nosig.in.pathway
        ct <- matrix(c(
          sp=f.sig.in.pathway,
          sn=f.sig.no.pathway,
          np=f.nosig.in.pathway,
          ns=f.nosig.no.pathway),2,2,dimnames=list(c("MOD","NOMOD"),c("SIG","NO SIG")))
        p.val.met = fisher.test(ct ,alt="less")$p.value
        #cat(paste0("Testing Mets: SigPath = ", f.sig.in.pathway, ", SigNopath = ", f.sig.no.pathway, ", NosigPath = ", f.nosig.in.pathway, ", NosigNopath = ", f.nosig.no.pathway, " | p.val = ", p.val.met, "\n"))
      } else {
        p.val.met = NA
      }
      all.res.df <- rbind(all.res.df, data.frame(
        ModuleID = mod.id,
        TRANS.pval = p.val.gene,
        PROT.pval = p.val.prot,
        MET.pval = p.val.met
      ))
    }
  }
  self$reports$enrichment[[name]][[method]] <- all.res.df
  return(all.res.df)
}



#' Plot the representative analytes of each assay (Proteomics, Metabolomics, Transcriptomics) vs scores
#' @param name What version of `scores_df` to use (see `ModObj$scores$...`)? Added via `$add_analyte_scores(...)`. Defaults to all (`NULL`).
#' @param module Which module to plot (as integer)? Defaults to `1`.
#' @param use Which variable to plot for Data? Can be `"score"` (score provided) or `"dir"` (direction)
#' @param invert Should the module directions be inverted? Defaults to `FALSE`
#' @param only.show.sig Only show significant analytes. Defaults to `TRUE`
#' @param all.assays Show assays even if no analytes are present? Defaults to `FALSE`
#' @param show.analyte.names Include x-axis names (analytes)? Defaults to `FALSE`
#' @examples
#' ModObj$plot_mapping()
#' 
#' ModObj$plot_mapping(name = "my_scores")
#' 
#'@export
plot_fingerprint = function(name = NULL, module = 1, use = "dir", invert = FALSE, only.show.sig = TRUE, all.assays = FALSE, show.analyte.names = FALSE){
  all.cur.mod.df <- data.frame()
  if(is.null(name)){
    name = names(self$scores)
  }
  for(cur.name in name){
    mod.id = paste0("Mod", module)
    mod.df <- self$modules
    cur.mod.df <- dplyr::filter(mod.df, ModuleID == mod.id)
    # invert?
    if(invert){
      cur.mod.df$Direction <- ifelse(cur.mod.df$Direction == "Pos", "Neg", "Pos")
    }
    cur.mod.df$Expr <- NA
    rownames(cur.mod.df) <- cur.mod.df$AnalyteID
    rel.scores <- sapply(1:nrow(self$scores[[cur.name]]), function(i){
      assay = self$scores[[cur.name]]$Assay[i]
      analyte.id = self$scores[[cur.name]]$ConvertedID[i]
      direction = ifelse(use == "dir", self$scores[[cur.name]]$Direction[i], self$scores[[cur.name]]$Score[i])
      significant = self$scores[[cur.name]]$Significant[i]
      if(only.show.sig & !significant){
        return(NA)
      }
      if(!is.na(analyte.id)){
        if(analyte.id %in% cur.mod.df$AnalyteID){
          return(direction)
        } else {
          return(NA)
        }
      } else {
        return(NA)
      }
    })
    names(rel.scores) <- self$scores[[cur.name]]$ConvertedID
    rel.scores <- rel.scores[!is.na(rel.scores)]
    cur.mod.df[names(rel.scores),"Expr"] <- rel.scores
    
    # arrange:
    cur.mod.df <- dplyr::arrange(cur.mod.df, desc(Assay), Direction)
    cur.mod.df$AnalyteID <- factor(cur.mod.df$AnalyteID, levels = cur.mod.df$AnalyteID)
    cur.mod.df$Name <- cur.name
    # Remove irrelevant assays:
    if(!all.assays){
      cur.mod.df <- dplyr::filter(cur.mod.df, Assay %in% unique(self$scores[[cur.name]]$Assay))
    }
    all.cur.mod.df <- rbind(all.cur.mod.df, cur.mod.df)
  }
  # Check if not important:
  if(nrow(cur.mod.df) == 0){
    stop("No significant analytes found mapped to this module. Change only.show.sig = FALSE and all.assays = TRUE to force show plot.")
  }
  
  # Add Mod:
  cur.mod.df$Name <- ifelse(invert, "Mod (Inv)", "Module")
  # Switch to numeric if using 'score':
  if(use == "score"){
    biggest.val = max(c(abs(max(all.cur.mod.df$Expr, na.rm = TRUE)), min(all.cur.mod.df$Expr, na.rm = TRUE)))
    cur.mod.df$Direction <- ifelse(cur.mod.df$Direction == "Pos", biggest.val, -(biggest.val))
  }
  cur.mod.df$Expr <- cur.mod.df$Direction
  all.cur.mod.df <- rbind(all.cur.mod.df, cur.mod.df)

  p = ggplot(all.cur.mod.df) +
        geom_tile(aes(x = AnalyteID, y = Name, fill = Expr), color = "black") +
        theme_bw() +
        theme(panel.grid = element_blank(), panel.border = element_blank()) +
        facet_wrap(vars(factor(Assay, levels = c("Transcriptomics", "Proteomics", "Metabolomics"))), scales = "free_x") +
        #scale_y_discrete(labels = c(name, )) +
        guides(fill = "none") +
        ylab(NULL) +
        xlab(NULL) +
        labs(color = NULL)
  
  if(!show.analyte.names){
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    p <- p + theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5))
  }
  
  if(use == "dir"){
    p <- p + scale_fill_manual(values = c(Pos = "red", Neg = "blue"), na.value = "grey80")
  } else if(use == "score"){
    p <- p + scale_fill_gradient2(high = "red", low = "blue", na.value = "grey80")
  }
  
  return(p)
}