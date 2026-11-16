dat <- read.csv("data-raw/m6A_imbalanced.csv", stringsAsFactors = FALSE)
dat$m6A_status <- factor(dat$m6A_status, levels = c("Negative", "Positive"))
model <- readRDS("inst/extdata/rf_fit.rds")
oob <- readRDS("data-raw/oob_predictions.rds")
truth <- oob$truth == "Positive"
p <- oob$prob_positive
pos_n <- sum(truth)
neg_n <- sum(!truth)
fig_dir <- "man/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

png(file.path(fig_dir, "fig01_dataset_overview.png"), width = 1800, height = 1200, res = 180)
par(mfrow = c(2, 2), mar = c(5, 4.5, 3, 1), oma = c(0, 0, 1, 0))
barplot(table(dat$m6A_status), col = c("#4C78A8", "#E45756"),
        main = "Class balance", ylab = "Candidate sites", las = 1)
barplot(sort(table(dat$RNA_type), decreasing = TRUE), col = "#72B7B2",
        main = "Transcript biotype", ylab = "Candidate sites", las = 2, cex.names = 0.8)
barplot(sort(table(dat$RNA_region), decreasing = TRUE), col = "#F2CF5B",
        main = "Transcript region", ylab = "Candidate sites", las = 2, cex.names = 0.8)
boxplot(m6A_level ~ m6A_status, data = dat, col = c("#4C78A8", "#E45756"),
        main = "Modification level", xlab = "Class", ylab = "m6A level")
dev.off()

png(file.path(fig_dir, "fig02_feature_distributions.png"), width = 1800, height = 1200, res = 180)
par(mfrow = c(2, 2), mar = c(5, 4.5, 3, 1))
for (v in c("gc_content", "exon_length", "distance_to_junction", "evolutionary_conservation")) {
  boxplot(dat[[v]] ~ dat$m6A_status, col = c("#4C78A8", "#E45756"),
          main = v, xlab = "Class", ylab = v, outline = FALSE)
}
dev.off()

imp <- randomForest::importance(model, type = 2)[, 1]
imp <- sort(imp, decreasing = TRUE)
png(file.path(fig_dir, "fig03_feature_importance.png"), width = 1600, height = 1100, res = 180)
par(mar = c(5, 12, 3, 1))
barplot(rev(imp), horiz = TRUE, las = 1, col = "#4C78A8",
        main = "Random forest feature importance", xlab = "Mean decrease in Gini")
dev.off()

ord <- order(p, decreasing = TRUE)
tp <- cumsum(truth[ord])
fp <- cumsum(!truth[ord])
tpr <- tp / pos_n
fpr <- fp / neg_n
recall <- tpr
precision <- tp / (tp + fp)
auc <- sum(diff(c(0, fpr)) * (c(0, tpr)[-1] + c(0, tpr)[-length(c(0, tpr))]) / 2)
ap <- sum(diff(c(0, recall)) * precision)
prevalence <- mean(truth)

png(file.path(fig_dir, "fig04_roc_curve.png"), width = 1500, height = 1100, res = 180)
par(mar = c(5, 5, 3, 1))
plot(c(0, fpr), c(0, tpr), type = "l", lwd = 3, col = "#2463A6",
     xlim = c(0, 1), ylim = c(0, 1), xlab = "False positive rate",
     ylab = "True positive rate", main = "Out-of-bag ROC curve")
abline(0, 1, lty = 2, col = "#888888")
legend("bottomright", sprintf("OOB AUC = %.3f", auc), bty = "n")
dev.off()

png(file.path(fig_dir, "fig05_prc_curve.png"), width = 1500, height = 1100, res = 180)
par(mar = c(5, 5, 3, 1))
plot(c(0, recall), c(1, precision), type = "l", lwd = 3, col = "#B34242",
     xlim = c(0, 1), ylim = c(0, 1), xlab = "Recall", ylab = "Precision",
     main = "Out-of-bag precision-recall curve")
abline(h = prevalence, lty = 2, col = "#888888")
legend("topright", legend = c(sprintf("Average precision = %.3f", ap),
       sprintf("Prevalence = %.3f", prevalence)), col = c("#B34242", "#888888"),
       lty = c(1, 2), lwd = c(3, 1), bty = "n")
dev.off()

thresholds <- seq(0.05, 0.95, by = 0.05)
sweep <- do.call(rbind, lapply(thresholds, function(th) {
  pred <- p > th
  TP <- sum(pred & truth); FP <- sum(pred & !truth)
  TN <- sum(!pred & !truth); FN <- sum(!pred & truth)
  data.frame(threshold = th, sensitivity = TP / (TP + FN),
             specificity = TN / (TN + FP),
             precision = ifelse(TP + FP == 0, NA, TP / (TP + FP)),
             f1 = ifelse(2 * TP + FP + FN == 0, NA, 2 * TP / (2 * TP + FP + FN)))
}))
png(file.path(fig_dir, "fig06_threshold_tradeoff.png"), width = 1700, height = 1000, res = 180)
matplot(sweep$threshold, sweep[, c("sensitivity", "specificity", "precision", "f1")],
        type = "l", lwd = 2.5, lty = 1, col = c("#E45756", "#4C78A8", "#72B7B2", "#F2CF5B"),
        ylim = c(0, 1), xlab = "Positive-class probability threshold", ylab = "Metric",
        main = "Threshold trade-off")
legend("bottomleft", legend = c("Sensitivity", "Specificity", "Precision", "F1"),
       col = c("#E45756", "#4C78A8", "#72B7B2", "#F2CF5B"), lty = 1, lwd = 2.5, bty = "n")
dev.off()

bin <- ceiling(rank(p, ties.method = "first") / length(p) * 12)
cal <- aggregate(cbind(predicted = p, observed = as.integer(truth)),
                 list(bin = bin), mean)
png(file.path(fig_dir, "fig07_calibration.png"), width = 1500, height = 1100, res = 180)
par(mar = c(5, 5, 3, 1))
plot(cal$predicted, cal$observed, type = "b", pch = 19, lwd = 2,
     col = "#2463A6", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Mean predicted probability", ylab = "Observed positive rate",
     main = "Out-of-bag calibration (12 bins)")
abline(0, 1, lty = 2, col = "#888888")
dev.off()

pred <- factor(ifelse(p > 0.5, "Positive", "Negative"), levels = c("Negative", "Positive"))
truth_f <- factor(ifelse(truth, "Positive", "Negative"), levels = c("Negative", "Positive"))
cm <- table(Truth = truth_f, Predicted = pred)
png(file.path(fig_dir, "fig08_confusion_matrix.png"), width = 1300, height = 1100, res = 180)
par(mar = c(5, 5, 4, 2))
image(x = 1:2, y = 1:2, z = t(cm), col = c("#EAF2F8", "#4C78A8"),
      axes = FALSE, xlab = "Predicted class", ylab = "True class",
      main = "Confusion matrix at threshold 0.5")
axis(1, at = 1:2, labels = colnames(cm))
axis(2, at = 1:2, labels = rev(rownames(cm)), las = 1)
text(rep(1:2, each = 2), rep(2:1, 2), labels = as.vector(cm), cex = 1.5)
dev.off()

cat(sprintf("OOB observations: %d\nOOB positives: %d\nOOB negatives: %d\nAUC: %.4f\nAverage precision: %.4f\n",
            length(p), pos_n, neg_n, auc, ap))
