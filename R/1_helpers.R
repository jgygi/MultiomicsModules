#' Get multi-omics module scores for a list of analytes.
#' @param name What version of `scores_df` to use (see `ModObj$scores$...`)? Added via `$add_analyte_scores(...)`. Defaults to `"table1"`.
#' @param method What method to derive scores? Defaults to `"fischer"`.
#' @param correction How to correct p.values (if method == `"fischer"`)? Defaults to `"BH"`. Use `"none"` to prevent correction. 
#' @param min.analytes Minimum number of analytes needed to calculate a score (per assay). Defaults to `0`.
#' @examples
#' get_module_scores(...)
#' 
#'@export
get_module_scores = function(name = NULL, method = "fischer", correction = "BH", min.analytes = 0){
  # check name:
  if(is.null(name)){
    name = names(self$scores)
  }
  if(method == "percentage"){
    for(cur.name in name){
      if(!self$params$quiet){cat("Calculating scores for name =", cur.name, "...\n")}
      all.cur.mod.df = data.frame()
      for(mod.id in unique(self$modules$ModuleID)){
        mod.df <- self$modules
        cur.mod.df <- dplyr::filter(mod.df, ModuleID == mod.id)
        cur.mod.df$Expr <- NA
        rownames(cur.mod.df) <- cur.mod.df$AnalyteID
        rel.scores <- sapply(1:nrow(self$scores[[cur.name]]), function(i){
          assay = self$scores[[cur.name]]$Assay[i]
          analyte.id = self$scores[[cur.name]]$ConvertedID[i]
          direction = self$scores[[cur.name]]$Direction[i]
          significant = self$scores[[cur.name]]$Significant[i]
          if(!significant){
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
        cur.mod.df$ModuleID <- mod.id

        # Get scores:
        assay.scores <- data.frame()
        for(assay in c("Transcriptomics", "Proteomics", "Metabolomics", "Total")){
          if(assay != "Total"){
            tmp <- dplyr::filter(cur.mod.df, Assay == assay)
          } else {
            tmp <- cur.mod.df
          }
          if(nrow(tmp) > 0 & sum(!is.na(tmp$Expr)) != 0 & sum(!is.na(tmp$Expr)) >= min.analytes){
            total.mapped <- self$scores[[cur.name]]
            if(assay != "Total"){total.mapped <- dplyr::filter(total.mapped, Assay == assay)}
            
            
            # TODO: Change jaccard total.mapped to just be the number of overlapping??? Very low...
            
            
            # Upregulated:
            jaccard_similarity.up <- sum(tmp$Expr == tmp$Direction, na.rm = TRUE) / 
              (nrow(tmp) + nrow(total.mapped) - sum(tmp$Expr == tmp$Direction, na.rm = TRUE))
            per.up <- sum(tmp$Expr == tmp$Direction, na.rm = TRUE)/sum(!is.na(tmp$Expr))
            # Downregulated:
            jaccard_similarity.down <- sum(tmp$Expr != tmp$Direction, na.rm = TRUE) / 
              (nrow(tmp) + nrow(total.mapped) - sum(tmp$Expr != tmp$Direction, na.rm = TRUE))
            per.down <- sum(tmp$Expr == ifelse(tmp$Direction == "Pos", "Neg", "Pos"), na.rm = TRUE)/sum(!is.na(tmp$Expr))
              
            assay.scores <- rbind(assay.scores, data.frame(
              ModuleID = mod.id,
              Assay = assay,
              PercentageUp = per.up,
              PercentageDown = per.down,
              JaccardUp = jaccard_similarity.up,
              JaccardDown = jaccard_similarity.down,
              CorrectedUp = per.up * jaccard_similarity.up,
              CorrectedDown = per.down * jaccard_similarity.down
            ))
          } else {
            assay.scores <- rbind(assay.scores, data.frame(
              ModuleID = mod.id,
              Assay = assay,
              PercentageUp = NA,
              PercentageDown = NA,
              JaccardUp = NA,
              JaccardDown = NA,
              CorrectedUp = NA,
              CorrectedDown = NA
            ))
          }
        }
        all.cur.mod.df <- rbind(all.cur.mod.df, assay.scores)
      }
      self$reports$percentage[[cur.name]] <- all.cur.mod.df
    }
    cat("Done! Results saved under $reports$scores$...\n")
    return(invisible(self))
  } else if(method == "fischer"){
    # Get scores for all names:
    module.scores <- list()
    mod.df <- self$modules
    for(cur.name in name){
      if(!self$params$quiet){cat("Calculating scores for name =", cur.name, "...\n")}
      all.p.val.df <- data.frame()
      for(assay in c("Transcriptomics", "Proteomics", "Metabolomics", "Combined")){
        tmp.p.val.df <- data.frame()
        for(mod.id in unique(mod.df$ModuleID)){
          if(assay != "Combined"){
            if(!assay %in% unique(self$scores[[cur.name]]$Assay)){
              # No analytes for this assay provided...
              next
            }
            cur.mod.df <- dplyr::filter(mod.df, ModuleID == mod.id, Assay == assay)
            total.num <- nrow(cur.mod.df)
          } else {
            cur.mod.df <- dplyr::filter(mod.df, ModuleID == mod.id)
            total.num <- nrow(cur.mod.df)
          }
          rownames(cur.mod.df) <- cur.mod.df$AnalyteID
          if(nrow(cur.mod.df) == 0){
            # No analytes for this module...
            next
          }
          cur.mod.df$P.value <- NA
          # Get p-values from all scores data.frames:
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
          # NOTE:: What to do if missing?? Just treat as 1's?
          # cur.mod.df$P.value[is.na(cur.mod.df$P.value)] <- 1
          
          found.cur.mod.df <- dplyr::filter(cur.mod.df, !is.na(P.value))
          found.cur.mod.df$FoundDir <- ifelse(found.cur.mod.df$P.value >= 0, "Pos", "Neg")
          found.cur.mod.df$P.value <- abs(found.cur.mod.df$P.value)

          # NOTE: change here for different types of combinations:
          # also: add check for min num analytes
          tmp <- dplyr::filter(found.cur.mod.df, Direction == FoundDir)
          if(nrow(tmp) > 0){
            reg.p.val <- poolr::fisher(tmp$P.value)$p
            reg.num <- nrow(tmp)
          } else {
            reg.p.val <- NA
            reg.num <- 0
          }
          tmp <- dplyr::filter(found.cur.mod.df, Direction != FoundDir)
          if(nrow(tmp) > 0){
            inv.p.val <- poolr::fisher(tmp$P.value)$p
            inv.num <- nrow(tmp)
          } else {
            inv.p.val <- NA
            inv.num <- 0
          }
          tmp <- found.cur.mod.df
          if(nrow(tmp) > 0){
            comb.p.val <- poolr::fisher(tmp$P.value)$p
            comb.num <- nrow(tmp)
          } else {
            comb.p.val <- NA
            comb.num <- 0
          }
          
          # add to results
          tmp.p.val.df <- rbind(tmp.p.val.df, data.frame(
            ModuleID = mod.id,
            Assay = assay,
            Pos_P.val = reg.p.val,
            Neg_P.val = inv.p.val,
            Comb_P.val = comb.p.val,
            TotalNum = total.num,
            PosNum = reg.num,
            InvNum = inv.num,
            CombNum = comb.num
          )) 
        }
        all.p.val.df <- rbind(all.p.val.df, tmp.p.val.df)
      }
      # correct p.values (if desired):
      if(correction != "none"){
        # Question: correct all together? Or just correct positive, negative, etc.? Maybe better to choose "best" (reg or inv) then correct?
        M <- as.matrix(dplyr::select(all.p.val.df, Pos_P.val, Neg_P.val, Comb_P.val))
        M[] <- p.adjust(M, method = correction)
        all.p.val.df$Pos_P.adj = M[,1]
        all.p.val.df$Neg_P.adj = M[,2]
        all.p.val.df$Comb_P.adj = M[,3]
      }
      
      self$reports$fischer[[cur.name]] <- all.p.val.df
    }
    cat("Done! Results saved under $reports$scores$...\n")
    return(invisible(self))
  }
}



#' Plot the representative analytes of each assay (Proteomics, Metabolomics, Transcriptomics) vs scores
#' @param name What version of `scores_df` to use (see `ModObj$scores$...`)? Added via `$add_analyte_scores(...)`. Defaults to all (`NULL`).
#' @param module What modules to plot? Defaults to all (`NULL`).
#' @param assay What assays to plot? Defaults to all (`NULL`).
#' @param type Type of report (generated from `get_module_scores`), defaults to `"fischer"`.
#' @param cap Should p-value magnitudes be capped? Useful when a few very significant p-values drive the visualization. Defaults to `1e-5`
#' @param map.cutoff Required proportion of analytes to be mapped? Defaults to `0.15` (15%). Use `0` to ignore cutoff.
#' @param show.adjust Show adjusted p.values? Defaults to `TRUE`.
#' @param show.all Show all p.values? Defaults to `FALSE`.
#' @examples
#' # run first:
#' ModObj$get_module_scores(...)
#' ModObj$plot_module_scores()
#' 
#'@export
plot_report = function(name = NULL, module = NULL, assay = NULL, type = "fischer", cap = 1e-20, map.cutoff = .15, show.adjust = TRUE,  show.all = FALSE){
  if(is.null(name)){
    name = names(self$scores)
  }
  if(is.null(module)){
    module = paste0("Mod", 1:100)
  } else {
    module = paste0("Mod", module)
  }
  if(is.null(assay)){
    assay = c("Transcriptomics", "Proteomics", "Metabolomics", "Combined")
  }
  
  ### Fischer plot:
  if(type == "fischer"){
    all.df <- data.frame()
    
    for(cur.name in name){
      tmp.df <- self$reports$fischer[[cur.name]]
      tmp.df <- dplyr::filter(tmp.df, ModuleID %in% module)
      tmp.df <- dplyr::filter(tmp.df, Assay %in% assay)
      if(nrow(tmp.df) > 0){
        tmp.df$Name <- cur.name
        all.df <- rbind(all.df, tmp.df)
      }
    }
    
    if(map.cutoff > 0){
      pos.props <- all.df$PosNum/all.df$TotalNum
      neg.props <- all.df$InvNum/all.df$TotalNum
      comb.props <- all.df$CombNum/all.df$TotalNum
      pos.keep <- pos.props >= map.cutoff
      neg.keep <- neg.props >= map.cutoff
      comb.keep <- comb.props >= map.cutoff
      if(show.adjust){
        all.df$Pos_P.adj[!pos.keep] <- NA
        all.df$Neg_P.adj[!neg.keep] <- NA
        all.df$Comb_P.adj[!comb.keep] <- NA
      } else {
        all.df$Pos_P.val[!pos.keep] <- NA
        all.df$Neg_P.val[!neg.keep] <- NA
        all.df$Comb_P.val[!comb.keep] <- NA
      }
    }
    
    if(show.adjust){
      all.df <- dplyr::select(all.df, Name, ModuleID, Assay, Pos_P.adj, Neg_P.adj, Comb_P.adj)
    } else {
      all.df <- dplyr::select(all.df, Name, ModuleID, Assay, Pos_P.val, Neg_P.val, Comb_P.val)
    }
    
    if(!show.all){
      all.df$P.adj <- sapply(1:nrow(all.df), function(i){
        p.pos = ifelse(show.adjust, all.df$Pos_P.adj[i], all.df$Pos_P.val[i])
        p.neg = ifelse(show.adjust, -all.df$Neg_P.adj[i], -all.df$Neg_P.val[i])
        # check for NA p-values:
        if(is.na(p.pos)){
          if(is.na(p.neg)){
            return(NA)
          } else {
            return(p.neg)
          }
        } else if(is.na(p.neg)){
          return(p.pos)
        }
        # If both aren't NA, take the best one:
        if(p.pos < abs(p.neg)){
          return(p.pos)
        } else {
          return(p.neg)
        }
      })
      all.df <- dplyr::select(all.df, Name, ModuleID, Assay, P.adj)
    } else {
      if(show.adjust){
        all.df$Neg_P.adj <- -all.df$Neg_P.adj
      } else {
        all.df$Neg_P.val <- -all.df$Neg_P.val
      }
    }
    
    tmp.melt <- reshape2::melt(all.df, id.vars = c("Name", "ModuleID", "Assay"))
    
    # cap magnitude?
    if(!is.null(cap)){
      tmp.melt$value[tmp.melt$value > 0 & tmp.melt$value < cap] <- cap
      tmp.melt$value[tmp.melt$value < 0 & tmp.melt$value > -cap] <- -cap
    }
    
    tmp.melt$SignedLogValue <- ifelse(tmp.melt$value > 0, -log10(tmp.melt$value), log10(-tmp.melt$value))
    
    
    if(show.all){
      # Plot to show all p.values (pos and inv)
      
      p <- ggplot2::ggplot(tmp.melt) +
        ggplot2::geom_tile(ggplot2::aes(x = variable, 
                                        y = factor(ModuleID, levels = paste0("Mod", 1:100)),
                                        fill = SignedLogValue), color = "black") +
        ggplot2::scale_fill_gradient2(low = "#2959E3", high = "#CE2525", mid = "white", na.value = "grey80") +
        ggplot2::ylab("Module") +
        ggplot2::xlab(NULL) +
        ggplot2::facet_grid(cols = ggplot2::vars(factor(Assay, levels = c("Transcriptomics", "Proteomics", "Metabolomics", "Combined"))),
                            rows = ggplot2::vars(Name)) +
        ggplot2::scale_x_discrete(labels = c("Pos", "Inv", "All")) 
    } else {
      # Plot to only show values per name...
      
      p <- ggplot2::ggplot(tmp.melt) +
        ggplot2::geom_tile(ggplot2::aes(x = Name, 
                                        y = factor(ModuleID, levels = paste0("Mod", 1:100)),
                                        fill = SignedLogValue), color = "black") +
        ggplot2::scale_fill_gradient2(low = "#2959E3", high = "#CE2525", mid = "white", na.value = "grey80") +
        ggplot2::ylab("Module") +
        ggplot2::xlab(NULL) +
        ggplot2::facet_grid(cols = ggplot2::vars(factor(Assay, levels = c("Transcriptomics", "Proteomics", "Metabolomics", "Combined")))) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = .5),
                       panel.background = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_blank())

    }

      
      return(p)
  } else {
    stop("ERROR: type argument not recognized or supported...")
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
      cat("Metabolites |\t", mapped.mets, "/", total.mets, "\t(", round(100*mapped.mets/total.mets, 2), "%)\n\n", sep = "")
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



#' Plot the representative analytes of each assay (Proteomics, Metabolomics, Transcriptomics) vs scores
#' @param name What version of `scores_df` to use (see `ModObj$scores$...`)? Added via `$add_analyte_scores(...)`. Defaults to all (`NULL`).
#' @param module Which module to plot (as integer)? Defaults to `1`.
#' @param use Which variable to plot for Data? Can be `"p.value"` (p.value, default), `"score"` (score provided) or `"dir"` (direction)
#' @param cap Provide a cap for the plot? Defaults to `NULL` (no cap). Alternatively, provide a numeric.
#' @param invert Should the module directions be inverted? Defaults to `FALSE`
#' @param only.show.sig Only show significant analytes. Defaults to `TRUE`
#' @param all.assays Show assays even if no analytes are present? Defaults to `FALSE`
#' @param show.analyte.names Include x-axis names (analytes)? Defaults to `TRUE`
#' @examples
#' ModObj$plot_mapping()
#' 
#' ModObj$plot_mapping(name = "my_scores")
#' 
#'@export
plot_fingerprint = function(name = NULL, module = 1, use = "p.value", cap = NULL, invert = FALSE, only.show.sig = TRUE, all.assays = FALSE, show.analyte.names = TRUE){
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
      if(use == "dir"){
        direction = self$scores[[cur.name]]$Direction[i]
      } else if(use == "score"){
        direction = self$scores[[cur.name]]$Score[i]
      } else if(use == "p.value"){
        direction = ifelse(self$scores[[cur.name]]$Direction[i] == "Pos", self$scores[[cur.name]]$P.value[i], -self$scores[[cur.name]]$P.value[i])
      }
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
    cur.mod.df$Expr <- ifelse(cur.mod.df$Direction == "Pos", biggest.val, -(biggest.val))
    # Add cap (if provided)
    if(!is.null(cap)){
      all.cur.mod.df$Expr[all.cur.mod.df$Expr > cap] <- cap
      all.cur.mod.df$Expr[all.cur.mod.df$Expr < -cap] <- -cap
    }
  } else if(use == "p.value"){
    biggest.val = max(-log10(abs(all.cur.mod.df$Expr)), na.rm = TRUE)
    all.cur.mod.df$Expr <- ifelse(all.cur.mod.df$Expr > 0, -log10(abs(all.cur.mod.df$Expr)), log10(abs(all.cur.mod.df$Expr)))
    # Add cap (if provided)
    if(!is.null(cap)){
      all.cur.mod.df$Expr[all.cur.mod.df$Expr > cap] <- cap
      all.cur.mod.df$Expr[all.cur.mod.df$Expr < -cap] <- -cap
    }
    cur.mod.df$Expr <- ifelse(cur.mod.df$Direction == "Pos", biggest.val, -(biggest.val))
  } else if(use == "dir"){
    cur.mod.df$Expr <- cur.mod.df$Direction
  }
  all.cur.mod.df <- rbind(all.cur.mod.df, cur.mod.df)
  
  p = ggplot2::ggplot(all.cur.mod.df) +
    ggplot2::geom_tile(ggplot2::aes(x = AnalyteID, y = Name, fill = Expr), color = "black") +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_blank()) +
        ggplot2::facet_wrap(ggplot2::vars(factor(Assay, levels = c("Transcriptomics", "Proteomics", "Metabolomics"))), scales = "free_x") +
        #ggplot2::scale_y_discrete(labels = c(name, )) +
        ggplot2::ylab(NULL) +
        ggplot2::xlab(NULL) +
        ggplot2::labs(color = NULL)
  
  if(!show.analyte.names){
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  } else {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0.5))
  }
  
  if(use == "dir"){
    p <- p + ggplot2::scale_fill_manual(values = c(Pos = "red", Neg = "blue"), na.value = "grey80")
  } else if(use == "score" | use == "p.value"){
    label = ifelse(use == "score", "Score", "Signed -log10(p.val)")
    p <- p + ggplot2::scale_fill_gradient2(high = "red", low = "blue", na.value = "grey80") +
      ggplot2::labs(fill = label) +
      ggplot2::theme(legend.position = "bottom")
  }
  
  return(p)
} 





