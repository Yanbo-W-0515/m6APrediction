pkg_dir <- normalizePath(".")
model <- readRDS(file.path(pkg_dir, "inst", "extdata", "rf_fit.rds"))
oob <- data.frame(
  truth = as.character(model$y),
  prob_positive = model$votes[, "Positive"],
  stringsAsFactors = FALSE
)
saveRDS(oob, file.path(pkg_dir, "data-raw", "oob_predictions.rds"))
