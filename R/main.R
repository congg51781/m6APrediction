#' Encode DNA Sequences into Position-Specific Factors
#'
#' This function converts DNA sequences into a data frame
#' where each column represents a specific nucleotide position.
#' Each position is stored as a factor with predefined
#' nucleotide levels. This uses factor encoding
#' rather than one-hot encoding.
#'
#' @param dna_strings A character vector containing DNA
#' sequences, typically 5-nucleotide DNA_5mer sequences.
#' All sequences should have the same length
#' and contain uppercase nucleotide letters (A, T, C, G).
#'
#' @return A data frame in which each column represents
#' one nucleotide position. Column names are generated
#' as nt_pos1, nt_pos2, and so on. Each column is a factor
#' with levels A, T, C, and G in that order.
#'
dna_encoding <- function(dna_strings) {
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(
    unlist(strsplit(dna_strings, "")),
    ncol = nn,
    byrow = TRUE
  )
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(
    seq_df,
    factor,
    levels = c("A", "T", "C", "G")
  )
  return(seq_df)
}






#' Predict m6A Status for Multiple Samples
#'
#' This function uses a trained Random Forest classification
#' model to predict m6A probabilities and classification
#' results for multiple samples. It processes RNA-related
#' features, sets predefined factor levels for categorical
#' variables, and encodes DNA sequences into nucleotide-
#' position factors. Samples are classified as Positive or
#' Negative based on the specified probability threshold.
#'
#' @param ml_fit A trained randomForest classification model
#' used to predict m6A probabilities.
#'
#' @param feature_df A data frame containing the input features:
#' gc_content, RNA_type, RNA_region, exon_length,
#' distance_to_junction, evolutionary_conservation,
#' and DNA_5mer. RNA_type, RNA_region, and DNA_5mer values should
#' match the categorical levels used during model training.
#'
#' @param positive_threshold A numeric value between 0 and 1
#' used as the classification threshold. Samples with
#' predicted probabilities greater than this value are
#' classified as Positive. Default is 0.5.
#'
#' @return A data frame containing the original input features
#' together with two additional columns:
#' predicted_m6A_prob, representing the predicted m6A probability,
#' and predicted_m6A_status, representing the predicted
#' classification (Positive or Negative).
#'
#' @examples
#' ml_fit <- readRDS(
#'   system.file("extdata", "rf_fit.rds",
#'               package = "m6APrediction")
#' )
#'
#' example_df <- read.csv(
#'   system.file("extdata", "m6A_input_example.csv",
#'               package = "m6APrediction")
#' )
#'
#' prediction_multiple(
#'   ml_fit,
#'   example_df,
#'   positive_threshold = 0.6
#' )
#'
#' @import randomForest
#' @importFrom stats predict
#' @export
prediction_multiple <- function(
    ml_fit, feature_df, positive_threshold = 0.5) {

  stopifnot(all(c(
    "gc_content",
    "RNA_type",
    "RNA_region",
    "exon_length",
    "distance_to_junction",
    "evolutionary_conservation",
    "DNA_5mer"
  ) %in% colnames(feature_df)))

  new_data <- feature_df

  new_data$RNA_type <- factor(
    new_data$RNA_type,
    levels = c(
      "mRNA", "lincRNA", "lncRNA", "pseudogene"
    )
  )

  new_data$RNA_region <- factor(
    new_data$RNA_region,
    levels = c(
      "CDS", "intron", "3'UTR", "5'UTR"
    )
  )

  new_data <- cbind(
    new_data,
    dna_encoding(new_data$DNA_5mer)
  )

  predicted_m6A_prob <- predict(
    ml_fit,
    newdata = new_data,
    type = "prob"
  )[, "Positive"]

  predicted_m6A_status <- ifelse(
    predicted_m6A_prob > positive_threshold,
    "Positive",
    "Negative"
  )

  feature_df$predicted_m6A_prob <- predicted_m6A_prob

  feature_df$predicted_m6A_status <- predicted_m6A_status

  return(feature_df)
}








#' Predict m6A Status for a Single Sample
#'
#' This function predicts the m6A probability and classification
#' status for a single sample. It converts the input features
#' into a data frame and uses prediction_multiple() to
#' perform the prediction.
#'
#' @param ml_fit A trained randomForest classification model
#' used for prediction.
#'
#' @param gc_content A numeric value representing the GC content
#' of the DNA sequence.
#'
#' @param RNA_type A character string specifying the RNA type.
#' Accepted values are mRNA, lincRNA, lncRNA, and pseudogene.
#'
#' @param RNA_region A character string specifying the RNA region.
#' Accepted values are CDS, intron, 3'UTR, and 5'UTR.
#'
#' @param exon_length A numeric value representing the exon length.
#'
#' @param distance_to_junction A numeric value representing the
#' distance to the nearest exon junction.
#'
#' @param evolutionary_conservation A numeric value representing
#' the evolutionary conservation score.
#'
#' @param DNA_5mer A character string representing the DNA 5-mer
#' sequence.
#'
#' @param positive_threshold A numeric value between 0 and 1
#' used as the classification threshold. Samples with predicted
#' probabilities greater than this value are classified as
#' Positive; otherwise, they are classified as Negative.
#' Default is 0.5.
#'
#' @return A named vector containing two elements:
#' predicted_m6A_prob, representing the predicted probability,
#' and predicted_m6A_status, representing the predicted
#' classification (Positive or Negative). The vector elements
#' are named predicted_m6A_prob and predicted_m6A_status.
#'
#' @examples
#' ml_fit <- readRDS(
#'   system.file("extdata", "rf_fit.rds",
#'               package = "m6APrediction")
#' )
#'
#' prediction_single(
#'   ml_fit = ml_fit,
#'   gc_content = 0.6,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 12,
#'   distance_to_junction = 5,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "GGACA",
#'   positive_threshold = 0.5
#' )
#'
#'
#' @export
prediction_single <- function(
    ml_fit,
    gc_content,
    RNA_type,
    RNA_region,
    exon_length,
    distance_to_junction,
    evolutionary_conservation,
    DNA_5mer,
    positive_threshold = 0.5) {

  single_df <- data.frame(
    gc_content = gc_content,

    RNA_type = factor(
      RNA_type,
      levels = c(
        "mRNA", "lincRNA", "lncRNA", "pseudogene"
      )
    ),

    RNA_region = factor(
      RNA_region,
      levels = c(
        "CDS", "intron", "3'UTR", "5'UTR"
      )
    ),
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer
  )
  prediction_result <- prediction_multiple(
    ml_fit = ml_fit,
    feature_df = single_df,
    positive_threshold = positive_threshold
  )
  returned_vector <- c(
    prediction_result$predicted_m6A_prob[1],
    prediction_result$predicted_m6A_status[1]
  )
  names(returned_vector) <- c(
    "predicted_m6A_prob",
    "predicted_m6A_status"
  )
  return(returned_vector)
}

