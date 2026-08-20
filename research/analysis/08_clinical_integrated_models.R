#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(glmnet)
})

root <- Sys.getenv("MENIBREAK_PROJECT_ROOT")
if (!nzchar(root)) stop("Set MENIBREAK_PROJECT_ROOT to the analysis project root; see config/path_config.example.sh")
out_dir <- file.path(root, "analysis_folders", "Clinical_Integrated_Model_Final_Figures_2026")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

feature_path <- file.path(root, "analysis_folders", "ml_breakpoint_survival_stratification_2026",
                          "ml_breakpoint_model_dataset_2026.csv")
outcome_path <- file.path(root, "software", "mbpr", "analysis", "bakeoff_2x2", "outcomes.csv")

raw <- read_csv(feature_path, show_col_types = FALSE)
outcomes <- read_csv(outcome_path, show_col_types = FALSE)

scna_features <- grep("^cnv_|^scna_(loss_count|gain_count|altered_arm_count)$",
                      names(raw), value = TRUE)
breakpoint_candidates <- c(
  "total_breakpoint_boundaries", "breakpoint_arms_touched",
  "breakpoint_route_entropy", "max_arm_breakpoint_fraction",
  "telomeric_breakpoint_fraction", "central_breakpoint_fraction",
  "centromeric_breakpoint_fraction", "mean_rel_tel_to_cent",
  "sd_rel_tel_to_cent", "cna_segment_count", "loss_segment_count",
  "gain_segment_count", "median_segment_length", "small_segment_count_lt5mb",
  "dynamic_hotspot_breakpoints_all", "dynamic_hotspot_breakpoints_non1p",
  "dynamic_hotspot_arms_used_all", "dynamic_hotspot_arms_used_non1p"
)
breakpoint_features <- intersect(breakpoint_candidates, names(raw))
clinical_features <- c("age", "log1p_mitoses", "gtr", "prior_rt", "adj_rt")

dat <- raw %>%
  transmute(
    sample_id = cg,
    age = as.numeric(age),
    log1p_mitoses = log1p(as.numeric(mitoses)),
    gtr = as.numeric(gtr),
    prior_rt = as.numeric(prior_rt),
    adj_rt = as.numeric(adj_rt),
    who_grade = as.integer(who_grade),
    integrated_grade = as.integer(integrated_grade),
    across(all_of(c(scna_features, breakpoint_features)))
  ) %>%
  inner_join(outcomes, by = "sample_id")

fixed_horizon_endpoint <- function(time, event, horizon) {
  ifelse(event == 1 & time <= horizon, 1,
         ifelse(time >= horizon, 0, NA_real_))
}

auc_safe <- function(y, score) {
  ok <- is.finite(y) & is.finite(score)
  y <- y[ok]; score <- score[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(score, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

group_folds <- function(groups, k, seed) {
  set.seed(seed)
  ids <- sample(unique(as.character(groups)))
  map <- setNames(rep(seq_len(k), length.out = length(ids)), ids)
  unname(map[as.character(groups)])
}

fit_enet_predict <- function(train, test, feature_names, seed) {
  xtr <- as.matrix(train[, feature_names, drop = FALSE])
  xte <- as.matrix(test[, feature_names, drop = FALSE])
  storage.mode(xtr) <- storage.mode(xte) <- "double"
  med <- apply(xtr, 2, function(z) median(z[is.finite(z)], na.rm = TRUE))
  med[!is.finite(med)] <- 0
  for (j in seq_along(feature_names)) {
    xtr[!is.finite(xtr[, j]), j] <- med[[j]]
    xte[!is.finite(xte[, j]), j] <- med[[j]]
  }
  foldid <- group_folds(train$patient_id, min(5L, n_distinct(train$patient_id)), seed)
  fit <- cv.glmnet(xtr, train$endpoint, family = "binomial", alpha = 0.25,
                   foldid = foldid, type.measure = "deviance", standardize = TRUE)
  drop(predict(fit, xte, s = "lambda.1se", type = "link"))
}

clinical_scna_features <- c(scna_features, clinical_features)
clinical_full_features <- c(scna_features, breakpoint_features, clinical_features)

prediction_cache <- file.path(out_dir, "clinical_integrated_oof_predictions.csv")
reuse_cache <- identical(Sys.getenv("MBPR_REUSE_OOF", "0"), "1")
if (reuse_cache && file.exists(prediction_cache)) {
  pred <- read_csv(prediction_cache, show_col_types = FALSE)
} else {
  predictions <- list()
  set.seed(20260802)
  for (h in c(3, 5, 7)) {
    d <- dat %>% mutate(endpoint = fixed_horizon_endpoint(recurrence_time, recurrence_event, h)) %>%
      filter(!is.na(endpoint))
    for (r in seq_len(10L)) {
      outer <- group_folds(d$patient_id, 5L, 20260703L + h * 10000L + r)
      p_scna <- p_full <- rep(NA_real_, nrow(d))
      for (fold in seq_len(5L)) {
        tr <- d[outer != fold, , drop = FALSE]
        te <- d[outer == fold, , drop = FALSE]
        p_scna[outer == fold] <- fit_enet_predict(
          tr, te, clinical_scna_features,
          20260703L + h * 100000L + r * 1000L + fold
        )
        p_full[outer == fold] <- fit_enet_predict(
          tr, te, clinical_full_features,
          20260703L + h * 200000L + r * 1000L + fold
        )
      }
      predictions[[length(predictions) + 1L]] <- tibble(
        sample_id = d$sample_id, patient_id = d$patient_id,
        horizon_years = h, repetition = r, endpoint = d$endpoint,
        who_grade = d$who_grade, integrated_grade = d$integrated_grade,
        raw_clinical_scna = p_scna, raw_clinical_full = p_full
      )
    }
  }
  pred <- bind_rows(predictions)
  write_csv(pred, prediction_cache)
}

bootstrap_auc <- function(data, B = 2000L, seed = 20260802L) {
  set.seed(seed)
  patients <- sort(unique(data$patient_id))
  reps <- sort(unique(data$repetition))
  bind_rows(lapply(seq_len(B), function(b) {
    sampled <- sample(patients, length(patients), replace = TRUE)
    wt <- table(sampled)
    per_rep <- bind_rows(lapply(reps, function(r) {
      dr <- data[data$repetition == r, , drop = FALSE]
      times <- as.integer(wt[match(dr$patient_id, names(wt))]); times[is.na(times)] <- 0L
      idx <- rep(seq_len(nrow(dr)), times = times)
      tibble(
        auc_scna = auc_safe(dr$endpoint[idx], dr$raw_clinical_scna[idx]),
        auc_full = auc_safe(dr$endpoint[idx], dr$raw_clinical_full[idx])
      )
    }))
    tibble(draw = b, auc_scna = mean(per_rep$auc_scna), auc_full = mean(per_rep$auc_full),
           delta = mean(per_rep$auc_full - per_rep$auc_scna))
  }))
}

summarize_horizons <- function(data) {
  bind_rows(lapply(c(3, 5, 7), function(h) {
    d <- filter(data, horizon_years == h)
    point <- d %>% group_by(repetition) %>% summarise(
      n = n(), events = sum(endpoint),
      auc_scna = auc_safe(endpoint, raw_clinical_scna),
      auc_full = auc_safe(endpoint, raw_clinical_full), .groups = "drop")
    boot <- bootstrap_auc(d, seed = 20260802L + h)
    tibble(
      horizon_years = h, n = first(point$n), events = first(point$events),
      auc_scna = mean(point$auc_scna), auc_full = mean(point$auc_full),
      delta = mean(point$auc_full - point$auc_scna),
      scna_low = quantile(boot$auc_scna, .025), scna_high = quantile(boot$auc_scna, .975),
      full_low = quantile(boot$auc_full, .025), full_high = quantile(boot$auc_full, .975),
      delta_low = quantile(boot$delta, .025), delta_high = quantile(boot$delta, .975),
      delta_p = 2 * min(mean(boot$delta <= 0), mean(boot$delta >= 0))
    )
  }))
}

summary_main <- summarize_horizons(pred)
write_csv(summary_main, file.path(out_dir, "clinical_integrated_auc_summary.csv"))

model_cols <- c(
  "Recorded WHO grade" = "#7F7F7F",
  "Prespecified integrated grade" = "#7876B1",
  "Clinical + SCNA" = "#0072B5",
  "Clinical + SCNA + breakpoint architecture" = "#BC3C29"
)
model_shapes <- c("Recorded WHO grade"=16, "Prespecified integrated grade"=17,
                  "Clinical + SCNA"=15, "Clinical + SCNA + breakpoint architecture"=18)
model_linetypes <- c("Recorded WHO grade"="dotdash", "Prespecified integrated grade"="dashed",
                     "Clinical + SCNA"="solid", "Clinical + SCNA + breakpoint architecture"="solid")
theme_pub <- theme_classic(base_family = "sans", base_size = 12.5) +
  theme(plot.title=element_text(face="bold",size=13),
        plot.subtitle=element_text(size=10.5,color="grey35"),
        legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=9.5),
        legend.key.width=unit(1.05,"cm"),panel.grid.major.y=element_line(color="grey90",linewidth=.35),
        plot.margin=margin(10,12,10,10))

make_figure4 <- function(s) {
  perf <- bind_rows(
    transmute(s,horizon_years,model="Clinical + SCNA",auc=auc_scna,low=scna_low,high=scna_high),
    transmute(s,horizon_years,model="Clinical + SCNA + breakpoint architecture",auc=auc_full,low=full_low,high=full_high)
  )
  p1 <- ggplot(perf,aes(horizon_years,auc,color=model,shape=model,linetype=model,group=model)) +
    geom_line(linewidth=.9)+geom_point(size=2.6)+geom_errorbar(aes(ymin=low,ymax=high),width=.1,linewidth=.55)+
    scale_color_manual(values=model_cols)+scale_shape_manual(values=model_shapes)+
    scale_linetype_manual(values=model_linetypes)+scale_x_continuous(breaks=c(3,5,7))+
    scale_y_continuous(limits=c(.55,.90),breaks=seq(.55,.90,.05))+
    labs(title="A  Recurrence discrimination",x="Recurrence horizon (years)",y="Cross-validated AUC")+theme_pub
  p2 <- ggplot(s,aes(horizon_years,delta))+
    geom_hline(yintercept=0,linetype=2,color="grey55")+
    geom_errorbar(aes(ymin=delta_low,ymax=delta_high),width=.1,color=model_cols[[4]],linewidth=.7)+
    geom_point(color=model_cols[[4]],size=2.6)+
    geom_label(aes(y=delta_high+.004,label=sprintf("%+.3f",delta)),family="sans",fontface="bold",size=3.3,
               color=model_cols[[4]],fill="white",label.size=0,label.padding=unit(.08,"lines"))+
    scale_x_continuous(breaks=c(3,5,7))+
    coord_cartesian(ylim=range(c(0,s$delta_low,s$delta_high))+c(-.012,.020))+
    labs(title="B  Incremental breakpoint value",x="Recurrence horizon (years)",
         y="AUC difference vs clinical + SCNA")+theme_pub+theme(legend.position="none")
  denom <- s %>% mutate(nonevents=n-events) %>% select(horizon_years,events,nonevents) %>%
    pivot_longer(c(events,nonevents),names_to="status",values_to="tumors") %>%
    mutate(status=factor(status,levels=c("events","nonevents"),labels=c("Recurrence","No recurrence by horizon")))
  p3 <- ggplot(denom,aes(factor(horizon_years),tumors,fill=status))+geom_col(width=.64)+
    geom_text(data=s,aes(x=factor(horizon_years),y=n,label=paste0("n = ",n)),inherit.aes=FALSE,
              vjust=-.4,family="sans",size=3.2,fontface="bold")+
    scale_fill_manual(values=c("Recurrence"="#E18727","No recurrence by horizon"="#D9DEE3"))+
    scale_y_continuous(expand=expansion(mult=c(0,.1)))+
    labs(title="C  Fixed-horizon cohorts",x="Recurrence horizon (years)",y="Tumors")+theme_pub
  fig <- (p1|p2)/p3+plot_layout(heights=c(1,.72),guides="collect")+
    plot_annotation(title="Clinical-integrated SCNA-plus-breakpoint model performance",
      subtitle="Age, log mitotic count, prior/adjuvant radiotherapy, and extent of resection included in both models",
      theme=theme(plot.title=element_text(family="sans",face="bold",size=16),
                  plot.subtitle=element_text(family="sans",size=10.5,color="grey35"),
                  plot.margin=margin(12,14,10,14))) & theme(legend.position="bottom")
  ggsave(file.path(out_dir,"figure4_clinical_integrated_recurrence_discrimination.png"),fig,width=9.4,height=7.4,dpi=400,bg="white")
  ggsave(file.path(out_dir,"figure4_clinical_integrated_recurrence_discrimination.pdf"),fig,width=9.4,height=7.4)
}

bootstrap_common <- function(data, B=2000L, seed=20260822L) {
  models <- names(model_cols)
  cols <- c("who_grade","integrated_grade","raw_clinical_scna","raw_clinical_full")
  point <- bind_rows(lapply(seq_along(models),function(i) data %>% group_by(horizon_years,repetition) %>%
    summarise(n=n(),events=sum(endpoint),auc=auc_safe(endpoint,.data[[cols[[i]]]]),.groups="drop") %>%
    group_by(horizon_years) %>% summarise(model=models[[i]],n=first(n),events=first(events),auc=mean(auc),.groups="drop")))
  set.seed(seed); draws <- list()
  for(h in c(3,5,7)) {
    d <- filter(data,horizon_years==h); patients <- unique(d$patient_id)
    for(b in seq_len(B)) {
      wt <- table(sample(patients,length(patients),replace=TRUE))
      per_rep <- lapply(sort(unique(d$repetition)),function(r) {
        dr <- d[d$repetition==r,,drop=FALSE]
        times <- as.integer(wt[match(dr$patient_id,names(wt))]);times[is.na(times)]<-0L
        idx <- rep(seq_len(nrow(dr)),times=times)
        vapply(cols,function(v) auc_safe(dr$endpoint[idx],dr[[v]][idx]),numeric(1))
      })
      aucs <- rowMeans(do.call(cbind,per_rep),na.rm=TRUE)
      draws[[length(draws)+1L]] <- tibble(horizon_years=h,draw=b,model=models,auc=aucs,delta=aucs[[4]]-aucs)
    }
  }
  boot <- bind_rows(draws)
  point %>% left_join(boot %>% group_by(horizon_years,model) %>% summarise(
    low=quantile(auc,.025),high=quantile(auc,.975),delta_vs_full=median(delta),
    delta_low=quantile(delta,.025),delta_high=quantile(delta,.975),.groups="drop"),by=c("horizon_years","model"))
}

common <- pred %>% filter(who_grade %in% 1:3, integrated_grade %in% 1:3)
common_summary <- bootstrap_common(common)
write_csv(common_summary,file.path(out_dir,"clinical_integrated_common_comparator_auc_summary.csv"))

make_figure5 <- function(s) {
  s <- mutate(s,model=factor(model,levels=names(model_cols)))
  delta <- s %>% filter(model!="Clinical + SCNA + breakpoint architecture") %>%
    mutate(x_plot=horizon_years+case_when(model=="Recorded WHO grade"~-.16,
                                          model=="Prespecified integrated grade"~0,TRUE~.16))
  p1 <- ggplot(s,aes(horizon_years,auc,color=model,shape=model,linetype=model,group=model))+
    geom_line(linewidth=.85)+geom_point(size=2.35)+geom_errorbar(aes(ymin=low,ymax=high),width=.09,linewidth=.5)+
    scale_color_manual(values=model_cols,drop=FALSE)+scale_shape_manual(values=model_shapes,drop=FALSE)+
    scale_linetype_manual(values=model_linetypes,drop=FALSE)+scale_x_continuous(breaks=c(3,5,7))+
    scale_y_continuous(limits=c(.48,.92),breaks=seq(.5,.9,.1))+
    labs(title="A  Head-to-head discrimination",x="Recurrence horizon (years)",y="AUC")+theme_pub+
    guides(color=guide_legend(nrow=2,byrow=TRUE),shape=guide_legend(nrow=2,byrow=TRUE),linetype=guide_legend(nrow=2,byrow=TRUE))
  p2 <- ggplot(delta,aes(x_plot,delta_vs_full,color=model,shape=model,linetype=model,group=model))+
    geom_hline(yintercept=0,linetype=2,color="grey60")+geom_line(linewidth=.75)+geom_point(size=2.2)+
    geom_errorbar(aes(ymin=delta_low,ymax=delta_high),width=.08,linewidth=.45)+
    scale_color_manual(values=model_cols,drop=FALSE)+scale_shape_manual(values=model_shapes,drop=FALSE)+
    scale_linetype_manual(values=model_linetypes,drop=FALSE)+scale_x_continuous(breaks=c(3,5,7),limits=c(2.6,7.4))+
    labs(title="B  Incremental AUC over each comparator",x="Recurrence horizon (years)",y="AUC difference")+
    theme_pub+guides(color="none",shape="none",linetype="none")
  denom <- s %>% filter(model=="Recorded WHO grade") %>% mutate(nonevents=n-events) %>%
    select(horizon_years,events,nonevents) %>% pivot_longer(c(events,nonevents),names_to="status",values_to="tumors") %>%
    mutate(status=factor(status,levels=c("events","nonevents"),labels=c("Recurrence","No recurrence by horizon")))
  p3 <- ggplot(denom,aes(factor(horizon_years),tumors,fill=status))+geom_col(width=.64)+
    geom_text(data=s %>% filter(model=="Recorded WHO grade"),aes(x=factor(horizon_years),y=n,label=paste0("n = ",n)),
              inherit.aes=FALSE,vjust=-.4,family="sans",size=3.2,fontface="bold")+
    scale_fill_manual(values=c("Recurrence"="#E18727","No recurrence by horizon"="#D9DEE3"),guide="none")+
    scale_y_continuous(expand=expansion(mult=c(0,.1)))+
    labs(title="C  Common-comparator cohorts",x="Recurrence horizon (years)",y="Tumors")+theme_pub
  fig <- (p1|p2)/p3+plot_layout(heights=c(1,.7),guides="collect")+
    plot_annotation(title="Clinical-integrated model versus WHO and integrated grade",
      subtitle="Identical tumors; clinical variables are added to both genomic models",
      theme=theme(plot.title=element_text(family="sans",face="bold",size=16),
                  plot.subtitle=element_text(family="sans",size=10.5,color="grey35"),
                  plot.margin=margin(12,14,10,14))) & theme(legend.position="bottom")
  ggsave(file.path(out_dir,"figure5_clinical_integrated_common_comparator_auc.png"),fig,width=9.4,height=7.4,dpi=400,bg="white")
  ggsave(file.path(out_dir,"figure5_clinical_integrated_common_comparator_auc.pdf"),fig,width=9.4,height=7.4)
}

make_figure4(summary_main)
make_figure5(common_summary)

write_lines(c(
  "Clinical covariates added identically to both genomic models:",
  "- age (continuous)",
  "- log(1 + mitotic count)",
  "- gross-total resection indicator",
  "- prior radiotherapy indicator",
  "- postoperative/adjuvant radiotherapy indicator",
  "Missing values were median-imputed from the training portion of each inner/outer fit only.",
  "Patient grouping, fixed-horizon eligibility, 10 repeats, 5 outer folds, elastic-net alpha, and lambda.1se selection match the locked analysis."
), file.path(out_dir,"analysis_specification.txt"))

print(summary_main)
print(common_summary)
