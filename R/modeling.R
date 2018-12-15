#' @title Estimate the A.C.E. model from a phyloseq object
#'
#' @description Estimates the A.C.E. model of heritabily (with confidence intervals) from a phyloseq object using the \code{\link[mets]{twinlm}} function.
#'
#' @param ps A \code{phyloseq} object that contains \code{\link[phyloseq]{sample_data}}.
#'
#' @param zygosity The name of the column in sample_data() that contains the \code{zygosity} variable (either MZ or DZ).
#'
#' @param dz The level in the \code{zygosity} variable corresponding to dyzogitic twins.
#'
#' @param family_id The name of the column in sample_data() that contains the \code{family_id} variable. In other words, both members of a twin pair should have the same value in this column.
#'
#' @param verbose Logical. If FALSE, messages will be supressed.
#'
#' @return This function returns a \code{data.frame} with A.C.E. values (and confidence intervals) for each ASV.
#'
#' @examples
#' df <- calculate_ACE(ps, zygosity = "Zygosity", dz = "DZ", family_id = "Family_ID", verbose = TRUE)
calculate_ACE <- function(ps, zygosity, dz, family_id, verbose = TRUE){

  # Check Inputs ----------

  if(is.null(ps) | class(ps)[1] != "phyloseq"){
    message("Phyloseq object not found or is not an object of class 'phyloseq'.")
  }

  if(is.null(zygosity)){
    message("The `zygosity` variable is missing with no default. This is name of the column in sample_data() that contains the `zygosity` variable.")
  }

  if(!zygosity %in% sample_variables(ps)){
    message("The `zygosity` variable ", zygosity," is not found in sample_data().")
  }

  if(is.null(dz)){
    message("The `dz` variable is missing with no default. This is the level in the `zygosity` variable corresponding to the dyzogitic twins.")
  }

  if(!dz %in% levels(get_variable(ps, zygosity))){
    message("The `dz` variable is not found in levels(get_variable(", deparse(substitute(ps)) , ", ", zygosity, ").
            The `dz` variable should be the level in the `zygosity` variable corresponding to the dyzogitic twins.")
  }

  if(is.null(family_id)){
    message("The `family_id` variable is missing with no default. This is name of the column in sample_data() that contains the `family_id` variable.")
  }

  if(!family_id %in% sample_variables(ps)){
    message("The `family_id` variable ", family_id," is not found in sample_data().")
  }

  # Create data data.frames ----------

  ASVs <- data.frame(otu_table(ps)@.Data)
  sd <- data.frame(sample_data(ps))
  sd <- sd[,c(zygosity, family_id)]

  if(sum(ASVs %% 1 == 0) > 0){
    message("Warning, this looks like count data.
            If you would like normalize first, use the normalize_for_heritability() function.")
  }

  # Create holder data.frame ----------

  ACE_df <- data.frame(matrix(nrow = ntaxa(ps), ncol = 9))
  rownames(ACE_df) <- taxa_names(ps)
  colnames(ACE_df) <- c("A_Estimate", "C_Estimate", "E_Estimate", "A_2.5%", "C_2.5%", "E_2.5%", "A_97.5%", "C_97.5%", "E_97.5%")

  # Calculate Heritibility ----------

  for(i in taxa_names(ps)){
    if(verbose == TRUE){message(paste0("Modeling Feature: ", i))}
    asv <- ASVs[,i]
    model <- twinlm(asv ~ 1, data = sd, DZ = dz, zyg= zygosity , id = family_id)
    # Store coefficient in data.frame
    ACE_df[i,] <- as.numeric(as.character(summary(model)$coef))
  }

  # Add rank Column ----------

  ACE_df$rank <- "Other"
  ranks <- rev(c("ASV", "Genus", "Family", "Order", "Class", "Phylum"))
  for(i in ranks){
    ACE_df$rank[grep(i, rownames(ACE_df))] <- i
  }
  ACE_df$rank <- factor(ACE_df$rank, levels = c("ASV", "Genus", "Family", "Order", "Class", "Phylum"))

  tt <- tax_table(ps)@.Data
  ACE_tt <- cbind(tt, ACE_df)

  if(verbose == TRUE){message(paste0("The data.frame output contains A.C.E. Estimates and Confidence Intervals for ", ntaxa(ps), " taxonomic features."))}
  return(ACE_tt)
}


