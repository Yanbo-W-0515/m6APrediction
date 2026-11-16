#' Encode DNA sequences by nucleotide position
#'
#' Splits DNA sequence strings into individual nucleotides and stores
#' each position as a factor with levels A, T, C, and G.
#'
#' @param dna_strings A character vector of DNA sequences. Sequences should
#'   have the same length.
#' @return A data frame with one row per sequence and one factor column per
#'   nucleotide position.
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}
#' Predict m6A sites for multiple sequences
#'
#' Uses a fitted random forest model to predict the probability and status
#' of m6A sites from a feature table.
#'
#' @param ml_fit A fitted random forest model.
#' @param feature_df A data frame containing the required feature columns:
#'   gc_content, RNA_type, RNA_region, exon_length, distance_to_junction,
#'   evolutionary_conservation, and DNA_5mer.
#' @param positive_threshold Numeric probability threshold for a Positive
#'   prediction. Scores must be strictly greater than this value; scores equal
#'   to it are classified as Negative.
#' @return The input data frame with two added columns: predicted_m6A_prob
#'   and predicted_m6A_status.
#' @import randomForest
#' @importFrom stats predict
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                               package = "m6APrediction", mustWork = TRUE))
#' example_df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                                    package = "m6APrediction", mustWork = TRUE))
#' prediction_multiple(ml_fit, head(example_df, 3))
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df)))
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df)
}
#' Predict m6A status for a single sequence
#'
#' Uses a fitted random forest model to predict the probability and status
#' of m6A for one sequence and its associated features.
#'
#' @param ml_fit A fitted random forest model.
#' @param gc_content Numeric GC content feature.
#' @param RNA_type RNA type: mRNA, lincRNA, lncRNA, or pseudogene.
#' @param RNA_region RNA region: CDS, intron, 3'UTR, or 5'UTR.
#' @param exon_length Numeric exon length feature, using the scale expected
#'   by the model.
#' @param distance_to_junction Numeric distance to the nearest junction,
#'   using the scale expected by the model.
#' @param evolutionary_conservation Numeric evolutionary conservation
#'   feature.
#' @param DNA_5mer A five-nucleotide DNA sequence with A at the central
#'   position.
#' @param positive_threshold Numeric probability threshold. Scores must be
#'   strictly greater than this value to be classified as Positive.
#' @return A named character vector with `predicted_m6A_prob` and
#'   `predicted_m6A_status`.
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                               package = "m6APrediction", mustWork = TRUE))
#' prediction_single(ml_fit, gc_content = 0.6, RNA_type = "mRNA",
#'                   RNA_region = "CDS", exon_length = 12,
#'                   distance_to_junction = 5,
#'                   evolutionary_conservation = 0.8,
#'                   DNA_5mer = "GGACA", positive_threshold = 0.5)
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector)
}

