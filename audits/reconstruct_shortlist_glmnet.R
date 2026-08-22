#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(data.table); library(glmnet); library(pROC)})

args <- commandArgs(trailingOnly = TRUE)
data_path <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260803)
d <- fread(data_path)

fixed_endpoint <- function(time, event, h) {
  ifelse(event == 1 & time <= h, 1, ifelse(time >= h, 0, NA_real_))
}
group_folds <- function(groups, k, seed) {
  set.seed(seed); ids <- sample(unique(as.character(groups)))
  map <- setNames(rep(seq_len(k), length.out = length(ids)), ids)
  unname(map[as.character(groups)])
}
auc_safe <- function(y, s) as.numeric(pROC::auc(y, s, quiet = TRUE, direction = "<"))
chi_rank <- function(x, y) {
  vapply(seq_len(ncol(x)), function(j) {
    z <- x[, j] > 0
    tab <- table(factor(z, levels = c(FALSE, TRUE)), factor(y, levels = 0:1))
    suppressWarnings(chisq.test(tab, correct = FALSE)$statistic[[1]])
  }, numeric(1))
}

any_cols <- grep("^bp_[0-9]+[pq]_bin[0-9]{2}$", names(d), value = TRUE)
typed_cols <- grep("^bp_(loss|gain)_[0-9]+[pq]_bin[0-9]{2}$", names(d), value = TRUE)
restriction <- with(d, primary == 1 & gtr == 1 & prior_rt == 0 & adj_rt == 0 & rt_induced == 0 & angiomatous == 0)
rows <- list(); selected_rows <- list(); z <- 0L; zs <- 0L

for (family in c("all")) {
  cols0 <- if (family == "any") any_cols else c(any_cols, typed_cols)
  x0 <- as.matrix(d[, ..cols0]); storage.mode(x0) <- "numeric"; x0[!is.finite(x0)] <- 0
  for (support_basis in c("patient")) {
    support <- if (support_basis == "tumor") colSums(x0 > 0) else {
      ux <- rowsum(1 * (x0 > 0), group = as.character(d$patient_id), reorder = FALSE)
      colSums(ux > 0)
    }
    keep <- support >= 20
    x <- 1 * (x0[, keep, drop = FALSE] > 0); cols <- cols0[keep]
    for (k_top in c(8L)) for (lambda_rule in c("lambda.1se", "lambda.min")) {
      for (outer_scheme in c("cyclic", "stored")) {
        reps <- if (outer_scheme == "stored") 1L else 5L
        scores <- list()
        fold_scores <- list()
        for (h in c(3L, 5L, 7L)) {
          ep <- fixed_endpoint(d$rfs_time_years, d$rfs_event, h); known <- is.finite(ep)
          xx <- x[known, , drop = FALSE]; y <- ep[known]; groups <- d$patient_id[known]
          rest <- restriction[known]; stored <- as.integer(d$cv_fold[known])
          pred_mat <- matrix(NA_real_, nrow(xx), reps)
          for (r in seq_len(reps)) {
            outer <- if (outer_scheme == "stored") stored else group_folds(groups, 5L, 20260803L + r)
            for (fold in 1:5) {
              tr <- outer != fold; te <- outer == fold
              stat <- chi_rank(xx[tr, , drop = FALSE], y[tr]); chosen <- head(order(stat, decreasing = TRUE, na.last = NA), k_top)
              inner <- group_folds(groups[tr], 5L, 20260803L + h * 1000L + r * 100L + fold)
              fit <- cv.glmnet(xx[tr, chosen, drop = FALSE], y[tr], family = "binomial", alpha = .25,
                               foldid = inner, standardize = TRUE, type.measure = "deviance")
              pred_mat[te, r] <- drop(predict(fit, xx[te, chosen, drop = FALSE], s = lambda_rule, type = "link"))
              for (j in chosen) {zs <- zs + 1L; selected_rows[[zs]] <- data.table(family, support_basis, k_top, lambda_rule, outer_scheme, h, r, fold, feature = cols[j])}
            }
          }
          scores[[paste0(h, "b")]] <- mean(apply(pred_mat, 2, function(s) auc_safe(y, s)))
          scores[[paste0(h, "r")]] <- mean(apply(pred_mat, 2, function(s) auc_safe(y[rest], s[rest])))
          scores[[paste0(h, "ba")]] <- auc_safe(y, rowMeans(pred_mat))
          scores[[paste0(h, "ra")]] <- auc_safe(y[rest], rowMeans(pred_mat)[rest])
        }
        for (aggregation in c("mean_repeat", "averaged_prediction")) {
          suff <- if (aggregation == "mean_repeat") c("b", "r") else c("ba", "ra")
          vals <- c(scores[[paste0(3,suff[1])]],scores[[paste0(5,suff[1])]],scores[[paste0(7,suff[1])]],
                    scores[[paste0(3,suff[2])]],scores[[paste0(5,suff[2])]],scores[[paste0(7,suff[2])]])
          target <- c(.702,.735,.752,.715,.753,.757)
          z <- z + 1L
          rows[[z]] <- data.table(family,support_basis,n_candidates=ncol(x),k_top,lambda_rule,outer_scheme,aggregation,
                                  rmse=sqrt(mean((vals-target)^2)),max_abs_error=max(abs(vals-target)),
                                  auc_3y_broad=vals[1],auc_5y_broad=vals[2],auc_7y_broad=vals[3],
                                  auc_3y_restricted=vals[4],auc_5y_restricted=vals[5],auc_7y_restricted=vals[6])
        }
      }
    }
  }
}
result <- rbindlist(rows)[order(rmse,max_abs_error)]
fwrite(result, file.path(out_dir,"glmnet_configuration_results.csv"))
fwrite(rbindlist(selected_rows), file.path(out_dir,"glmnet_selected_features.csv"))
print(head(result,20))
