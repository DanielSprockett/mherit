#' @title Agglomerate and merge a phyloseq object
#'
#' @description Agglomerates by \code{rank_names()} and then merges them into a phyloseq object.
#'
#' @param ps A \code{phyloseq} object that contains a \code{\link[phyloseq]{tax_table}}.
#'
#' @param ranks A vector of taxonomic rank names found in \code{\link[phyloseq]{rank_names}}. Defaults to all available ranks (\code{Phylum} through \code{Genus}).
#'
#' @param replace_ASV_names Logical. If TRUE, ASVs will be named with a unique and the lowest taxonomic rank available.
#'
#' @param verbose Logical. If FALSE, messages will be supressed.
#'
#' @return This function returns the input \code{phyloseq} object with additional taxa that are the
#' agglomerated taxa at the taxonomic levels provided.
#'
#' @examples
#' ps <- add_ranks(ps, ranks = c("Genus", "Family"))

add_ranks <- function(ps = ps, ranks = c(rank_names(ps)), replace_ASV_names = TRUE, verbose = TRUE){

  # Check Inputs ----------

  if (is.null(ps) | class(ps)[1] != "phyloseq") {
    message("Phyloseq object not found or is not an object of class 'phyloseq'.")
  }

  if (is.null(tax_table(ps))) {
    message("Phyloseq object does not contain tax_table().")
  }

  # Replace Names ----------

  if(replace_ASV_names == TRUE){
    rank_list <- rev(rank_names(ps))
    taxa_list <- as.character(unname(tax_table(ps)[,rank_list[1]]))
    if(sum(is.na(taxa_list)) > 0){
      for(i in seq(2, c(length(rank_list) - 1))){
        taxa_list[is.na(taxa_list)] <- as.character(unname(tax_table(ps)[is.na(taxa_list),rank_list[i + 1]]))
      }
    }
    taxa_names(ps) <- paste0("ASV_", seq(1,length(taxa_names(ps))), "_", taxa_list)
    message("ASV names have been replaced.")
  }

  ranks <- ranks[which(!ranks %in% c("lowest_rank", "Species", "species", "Kingdom", "kingdom"))]

  # Check Ranks ----------

  null_nanks <- ranks[!ranks %in% rank_names(ps)]
  if(length(null_nanks) > 0){
    ranks <- ranks[which(!ranks %in% null_nanks)]
    message(paste0("WARNING: ", null_nanks, " not found in rank_names() and will be ignored."))
  }

  # Remove phy_tree() and refseqs() ----------
  if(!is.null(phy_tree(ps, errorIfNULL = FALSE)) | !is.null(refseq(ps, errorIfNULL = FALSE))){
    # remove phy_tree and
    ps <- phyloseq(otu_table(ps), sample_data(ps), tax_table(ps))
    message(paste0("WARNING: The phy_tree() and/or refseq() in your phyloseq object will be removed."))
  }

  ps_list <- vector("list", length = length(ranks) + 1)
  names(ps_list) <- c("ASV", rev(ranks))
  ps_list[["ASV"]] <- ps

  # Agglomerate By Ranks ----------

  for(i in rev(ranks)){
    i_ps <-  tax_glom(ps, i)
    taxa_names(i_ps) <- make.unique(paste0(i, "_", tax_table(i_ps)[,i]))
    ps_list[[i]] <- i_ps
    if(verbose == TRUE){message(paste0("Agglomerated by ", i, " (n = ", ntaxa(i_ps), ")."))}
  }

  # merge phyloseq objects----------

  merged_ps <- do.call(merge_phyloseq, ps_list)

  # Correct Taxa Names (if needed) ----------

  if(sum(!grepl("^[a-zA-Z0-9_.]*$", taxa_names(merged_ps))) > 0 ){
    if(verbose == TRUE){message("Replacing invalid charaters in taxa_names() with .")}
    taxa_names(merged_ps) <- make.names(taxa_names(merged_ps), unique = TRUE)
  }

  return(merged_ps)
}



#' @title Normalize count data in a phyloseq object
#'
#' @description Normalizes count data in a phyloseq object, first by adding a pseudocount (+1), then transforming to relative abundances, then by normalizes via the Box-Cox transformation in the \code{car} package.
#'
#' @param ps A \code{phyloseq} object that contains \code{\link[phyloseq]{sample_data()}}.
#'
#' @param covariates A vector of covariates found in \code{\link[phyloseq]{sample_data}}. Defaults to `x ~ 1`.
#'
#' @param normalization_plots Logical. If TRUE, creates a directory and saves plots of ASV distributions pre- and post-normalization. Defaults to FALSE.
#'
#' @param verbose Logical. If FALSE, messages will be supressed.
#'
#' @return This function returns the input \code{phyloseq} object with the \code{otu_table()} counts normalized.
#'
#' @examples
#' ps <- normalize_for_heritability(ps, covariates = c("Plate", "Age", "Sex"), normalization_plots = TRUE)
normalize_for_heritability <- function(ps, covariates = c(), normalization_plots = FALSE, verbose = TRUE){

  # Check Inputs ----------

  if (is.null(ps) | class(ps)[1] != "phyloseq") {
    message("phyloseq object not found or is not an object of class 'phyloseq'")
  }

  # Edit Taxa Names to avoid issues later
  if(sum(!grepl("^[a-zA-Z0-9_.]*$", taxa_names(ps))) > 0 ){
    if(verbose == TRUE){message("Replacing invalid charaters in taxa_names() with .")}
    taxa_names(ps) <- make.names(taxa_names(ps), unique = TRUE)
  }

  # Add a Pseduocount (+ 1) ----------

  if(verbose == TRUE){message("Adding pseduocount (+1) to ASV table.")}
  otu_table(ps) <- otu_table(ps) +1

  # Transform to Relative Abundance ----------

  if(verbose == TRUE){message("Transforming ASV table to relative abundances.")}
  ps <- transform_sample_counts(ps, function(x) x/sum(x))

  # Identify Covariates ----------

  if(length(covariates) == 0){message("No Covariates found. Defaulting to `x ~ 1`.")}

  if(length(covariates) > 0){

    message("Covariates found:")
    cov_list <- c()
    cov_df <- data.frame(matrix(ncol = length(covariates), nrow = nsamples(ps)))
    colnames(cov_df) <- covariates
    rownames(cov_df) <- sample_names(ps)

    for(cov in covariates){
      if(cov %in% sample_variables(ps)){
        var <- get_variable(ps, cov)
        cl <- class(var)
        if(is.numeric(var)){
          cov_df[,cov] = scale(var, center = TRUE)
          message(paste0("     ", cov, " (class: ", cl,", centered and scaled)"))
        } else {
          cov_df[,cov] = as.factor(var)
          message(paste0("     ", cov, " (class: ", cl,")"))
        }

        cov_list <- c(cov_list, cov)

      } else {
        cov_df[,cov] <- NULL
        message(paste0("     ", cov, " not found in sample_variables(), and will be ignored."))
      }
    }
  }
    cov_cc <- cov_df[complete.cases(cov_df),]
    if(length(rownames(cov_cc)) == 0){message("None of the supplied covariates are present for every sample.")}
    removed_samples <- rownames(cov_df)[!rownames(cov_df) %in% rownames(cov_cc)]
    #print(removed_samples)

    if(length(removed_samples) > 0) {
      cc_ps <- subset_samples(ps, !sample_names(ps) %in% removed_samples)
      ASVs <- data.frame(otu_table(cc_ps)@.Data)
      tASVs <- ASVs
      message(paste0("WARNING: ", length(removed_samples), " samples were removed due to missing covariate data."))
      if(length(removed_samples) < 5 ){
        message("Those samples are ", paste(removed_samples, collapse=", "))
      } else {
        message("The first 5 of those removed samples are ", paste(removed_samples[1:5], collapse=", "))
      }
    } else {
      ASVs <- data.frame(otu_table(ps)@.Data)
      tASVs <- ASVs
    }


  # Print Normalization Plots ----------

  if(normalization_plots == TRUE){
    if(verbose == TRUE){message("Printing normalization plots.")}
    dir.create(file.path(getwd(), "Normalization_Plots"), showWarnings = FALSE)
  }


  # Model Taxa  ----------

  for(i in taxa_names(ps)){
    if(verbose == TRUE){message(paste0("Modeling Feature: ", i))}
    asv <- tASVs[,i]
    # Get lambda
    lambda <- powerTransform(asv)$roundlam
    # Box-Cox Transformation
    bcp <- bcPower(asv,lambda)
    df <- cbind(bcp, cov_cc)
    if(length(covariates) == 0){cov_list <- c(1)}
    formula <- as.formula(paste("bcp", paste(cov_list, collapse=" + "), sep=" ~ "))
    model <- lm(formula, data = df)
    # Store Residuals in data.frame
    tASVs[,i] <- summary(model)$residuals

    if(normalization_plots == TRUE){
      df <- data.frame("Pre_Normalization" = ASVs[,i])
      p <- ggplot(df)
      p1 <- p + geom_density(aes(x = Pre_Normalization), fill = "#2B4162", alpha = .5)
      p1qq <- p + geom_qq(aes(sample = Pre_Normalization), color = "#2B4162")  + geom_qq_line(aes(sample = Pre_Normalization), color = "#2B4162")

      df <- data.frame("Post_Normalization" = tASVs[,i])
      p <- ggplot(df)
      p2 <- p + geom_density(aes(x = Post_Normalization), fill = "#0B6E4F", alpha = .5)
      p2qq <- p + geom_qq(aes(sample = Post_Normalization), color = "#0B6E4F") + geom_qq_line(aes(sample = Post_Normalization), color = "#0B6E4F")

      file_name <- paste0(i, "_Normalization_Plot")
      pdf(paste0(getwd(), "/Normalization_Plots/", file_name, ".pdf"), width=6,height=5)
      grid.arrange(p1, p2, p1qq, p2qq, ncol=2, nrow=2, top = i)
      invisible(dev.off())
    }
  }

  otu_table(ps) <- otu_table(tASVs, taxa_are_rows = FALSE)

  message("The otu_table() of the phyloseq output contains Box-Cox transformed residuals following linear regression. ")
  return(ps)

}


.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to the mherit R package. Have a wonderful day!")
}




