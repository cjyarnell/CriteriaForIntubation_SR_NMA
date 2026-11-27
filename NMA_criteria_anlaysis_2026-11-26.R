# Network meta-analysis and network meta-regression of intubation criteria in
# randomized trials of non-invasive respiratory support for acute hypoxemic
# respiratory failure 
#
# Input files, in the working directory:
#   Corrected_Data_Extraction 2026-11-18.xlsx     general and outcome data
#   intubation_criteria_by_study.csv   study-level intubation criteria and
#                                      work-of-breathing criteria flags
#   
#

########## Packages and working drive ##########

# Packages
library(tidyverse)
library(readxl)
library(gemtc)
library(igraph)
library(meta)
library(flextable)
library(officer)

set.seed(2026)

# setwd("path/to/input/files")

data_file     <- list.files(pattern = "^Corrected.Data.Extraction.*\\.xlsx$")[1]
criteria_file <- list.files(pattern = "^intubation.criteria.by.study.*\\.csv$")[1]
if (is.na(data_file))     stop("Extraction workbook not found in ", getwd())
if (is.na(criteria_file)) stop("Criteria CSV not found in ", getwd())

data_file     <- normalizePath(data_file)
criteria_file <- normalizePath(criteria_file)

dir.create("N01_priors", showWarnings = FALSE)
setwd("N01_priors")
dir.create("Figures", showWarnings = FALSE)

# Colours

c_ice   <- "#AAC0CE"
c_light <- "#DCBCBC"
c_mid   <- "#B97C7C"
c_dark  <- "#8F2727"


###################### Functions ###############################

# gemtc puts vague N(0, (15 * om.scale)^2) priors on its basic parameters (the
# edges of a spanning tree of the network, not necessarily comparisons with
# SOT) and Cauchy(0, om.scale) priors on regression coefficients.
# set_n01_priors() rewrites the generated JAGS code so that the log odds ratio
# of every treatment vs SOT has a N(0, 1) prior, and the criteria coefficient
# has a N(0, 1) prior. 

set_n01_priors <- function(model, reference = "SOT") {
  code <- model$code
  d_prior <- "d\\.[A-Za-z0-9_]+\\.[A-Za-z0-9_]+ ~ dnorm\\(0, prior\\.prec\\)"
  params <- sub(" ~.*", "", regmatches(code, gregexpr(d_prior, code))[[1]])
  if (length(params) == 0) stop("Relative-effect priors not found in the gemtc model code")

  # in node-splitting models the direct-evidence parameter keeps its own prior
  split_param <- if (!is.null(model$split)) paste("d", model$split[1], model$split[2], sep = ".") else character(0)
  tree <- setdiff(params, split_param)
  ends <- do.call(rbind, strsplit(sub("^d\\.", "", tree), ".", fixed = TRUE))
  trts <- setdiff(unique(c(ends)), reference)
  if (!(reference %in% ends) || length(tree) != length(trts)) {
    stop("Basic parameters do not form a spanning tree that includes ", reference)
  }

  # basic parameters become differences of effects vs the reference
  eff <- function(t) ifelse(t == reference, "0", paste0("e.", t))
  defs <- paste0(tree, " <- ", eff(ends[, 2]), " - ", eff(ends[, 1]))
  new_priors <- paste0("e.", trts, " ~ dnorm(0, 1)", collapse = "\n\t")
  for (i in seq_along(tree)) {
    line <- paste0(tree[i], " ~ dnorm(0, prior.prec)")
    code <- sub(line, if (i == 1) paste0(new_priors, "\n\t", defs[i]) else defs[i], code, fixed = TRUE)
  }
  for (p in split_param) {
    code <- sub(paste0(p, " ~ dnorm(0, prior.prec)"), paste0(p, " ~ dnorm(0, 1)"), code, fixed = TRUE)
  }

  code <- gsub("B ~ dt(0, reg.prior.prec, 1)", "B ~ dnorm(0, 1)", code, fixed = TRUE)
  code <- gsub("beta[k] ~ dt(0, reg.prior.prec, 1)", "beta[k] ~ dnorm(0, 1)", code, fixed = TRUE)
  if (grepl(d_prior, code) || grepl("dt(0, reg.prior.prec", code, fixed = TRUE)) {
    stop("Some treatment or covariate priors were not replaced")
  }
  model$code <- code

  # initial values: convert the basic-parameter inits to effects vs the reference
  A <- matrix(0, length(tree), length(trts), dimnames = list(tree, trts))
  for (i in seq_along(tree)) {
    if (ends[i, 2] != reference) A[i, ends[i, 2]] <- 1
    if (ends[i, 1] != reference) A[i, ends[i, 1]] <- -1
  }
  model$inits <- lapply(model$inits, function(init) {
    e <- solve(A, unlist(init[tree]))
    init[tree] <- NULL
    for (t in trts) init[[paste0("e.", t)]] <- unname(e[t])
    init
  })
  model
}

mtc_model_n01 <- function(...) set_n01_priors(gemtc::mtc.model(...))

# mtc.nodesplit() with the same priors (mirrors gemtc::mtc.nodesplit)
mtc_nodesplit_n01 <- function(network, comparisons = mtc.nodesplit.comparisons(network), ...) {
  results <- lapply(seq_len(nrow(comparisons)), function(i) {
    mtc.run(mtc_model_n01(network, type = "nodesplit",
                          t1 = as.character(comparisons$t1[i]),
                          t2 = as.character(comparisons$t2[i]), ...))
  })
  names(results) <- paste("d", comparisons$t1, comparisons$t2, sep = ".")
  results$consistency <- mtc.run(mtc_model_n01(network, type = "consistency", ...))
  class(results) <- "mtc.nodesplit"
  results
}


# Posterior draws from relative.effect(), all chains stacked, labelled by group
post_draws <- function(rel_effect, label) {
  as.data.frame(as.matrix(rel_effect$samples)) %>%
    mutate(Criteria = label)
}

# Node-splitting estimates on the odds ratio scale
extract_nodesplit_ors <- function(nodesplit) {

  bind_rows(
    mutate(summary(nodesplit)$dir.effect, Model = "Direct"),
    mutate(summary(nodesplit)$ind.effect, Model = "Indirect"),
    mutate(summary(nodesplit)$cons.effect, Model = "Network")
  ) %>%
    mutate(flip = t1,
           t1 = t2,
           t2 = flip,
           pe = -pe,
           flip2 = ci.u,
           ci.u = -ci.l,
           ci.l = -flip2) %>%
    select(-flip, -flip2) %>%
    mutate(
      pe = exp(pe),
      ci.l = exp(ci.l),
      ci.u = exp(ci.u),
      label = paste0(Model, ": ", round(pe, 2), " (", round(ci.l, 2), ", ", round(ci.u, 2), ")")
    ) %>%
    mutate(
      Comparison = paste(t1,t2, sep = ".")
    )
}

# Forest plot of node-splitting results
nodesplit_forestplot <- function(nodesplit, outcome, ylim_upper = 3) {

  nodesplit_results <- extract_nodesplit_ors(nodesplit)

  comparison_labels <- c(
    "HFNC.Bilevel_NIPPV" = "Bilevel NIPPV vs HFNC",
    "SOT.Bilevel_NIPPV" = "Bilevel NIPPV vs SOT",
    "HFNC.CPAP" = "CPAP vs HFNC",
    "SOT.CPAP" = "CPAP vs SOT",
    "SOT.HFNC" = "HFNC vs SOT"
  )
  nodesplit_results$Comparison <- factor(nodesplit_results$Comparison,
                                         levels = names(comparison_labels),
                                         labels = comparison_labels)

  nodesplit_results$Model <- factor(nodesplit_results$Model,
                                    levels = c("Network",
                                               "Indirect",
                                               "Direct"),
                                    ordered = T)

  ggplot(nodesplit_results, aes(x = Comparison, y = pe, ymin = ci.l, ymax = ci.u, group = Model)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
    geom_pointrange(position = position_dodge(width = 0.5), size = 0.8) +
    geom_point(aes(color = Model),
               position = position_dodge(width = 0.5), size = 3) +
    geom_text(aes(label = label),
              y = log(0.02), # fixed positions are on the log scale
              position = position_dodge(width = 0.5),
              size = 3,
              hjust = 0) +
    scale_y_continuous(
      name = "Odds Ratio (95% CrI)",
      trans = "log",
      breaks = c(0.5, 1, 2, 5),
      limits = c(0.02, ylim_upper)
    ) +
    scale_x_discrete(name = "") +
    scale_color_manual(
      values = c("Network" = "black", "Indirect" = "grey25", "Direct" = "grey90"),
      guide = "none"
    ) +
    labs(
      title = paste0("Node-Splitting Analysis for ", outcome),
      color = "Estimate Type"
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    coord_flip()
}

# Rank probabilities from posterior draws of log odds ratios vs SOT
# (rank 1 = lowest odds of the outcome)
rank_probs <- function(ranks){

  ranksamples <-
    as.data.frame(as.matrix(ranks)) %>%
    select(d.SOT.Bilevel_NIPPV,
           d.SOT.CPAP,
           d.SOT.HFNC) %>%
    mutate(SOT = 0)

  as.data.frame(t(apply(ranksamples, 1, rank))) %>%
    pivot_longer(cols = everything()) %>%
    group_by(name, value) %>%
    summarise(p = n()/nrow(ranksamples), .groups = "drop") %>%
    arrange(value)
}

rankplot <- function(rank_df1, rank_df0, title = ""){
  rank_df1 <- mutate(rank_df1,
                     Criteria = "Criteria")
  rank_df0 <- mutate(rank_df0,
                     Criteria = "No criteria")

  bind_rows(rank_df1,
            rank_df0) %>%
    mutate(Treatment = factor(name,
                              levels = c("SOT",
                                         "d.SOT.Bilevel_NIPPV",
                                         "d.SOT.HFNC",
                                         "d.SOT.CPAP"),
                              labels = c("SOT",
                                         "Bilevel NIPPV",
                                         "HFNC",
                                         "CPAP"))) %>%
    rename(Rank = value) %>%
    ggplot(aes(x = Rank,
               y = p,
               fill = Treatment)) +
    geom_col() +
    facet_grid(.~Criteria) +
    scale_fill_manual(values = c("grey90", c_light, c_mid, c_dark)) +
    theme_minimal()+
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank()) +
    labs(y = "Probability", title = title)

}

# Forest plots of odds ratios vs SOT, by covariate value
gemtc_forestplot <- function(result_criteria,
                             outcome,
                             ylim_upper = 2){

  rel_effect_1 <-  relative.effect(result_criteria, t1 = "SOT", covariate = 1)
  rel_effect_0 <-  relative.effect(result_criteria, t1 = "SOT", covariate = 0)

  rel_effect_df <-
    bind_rows(post_draws(rel_effect_1, "Yes"),
              post_draws(rel_effect_0, "No"))

  summary_data <- rel_effect_df %>%
    pivot_longer(
      cols = starts_with("d."),
      names_to = "Comparison",
      values_to = "LogOdds"
    ) %>%
    group_by(Comparison, Criteria) %>%
    summarize(
      mean = exp(mean(LogOdds)),
      lower = exp(quantile(LogOdds, 0.025)),
      upper = exp(quantile(LogOdds, 0.975)),
      .groups = "drop"
    )

  comparison_labels <- c(
    "d.SOT.HFNC" = "HFNC",
    "d.SOT.CPAP" = "CPAP",
    "d.SOT.Bilevel_NIPPV" = "Bilevel NIPPV"
  )
  summary_data$Comparison <- factor(summary_data$Comparison, levels = names(comparison_labels), labels = comparison_labels)

  summary_data$Criteria <- factor(summary_data$Criteria, labels = c("No Criteria", "Criteria"))

  summary_data <- summary_data %>%
    mutate(label = sprintf("%.2f (%.2f, %.2f)", mean, lower, upper))

  ggplot(summary_data, aes(x = Comparison, y = mean, ymin = lower, ymax = upper, group = Criteria)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
    geom_pointrange(position = position_dodge(width = 0.5), size = 0.8) +
    geom_point(aes(color = Criteria),
               position = position_dodge(width = 0.5),
               size = 3) +
    geom_text(
      aes(label = label),
      y = log(0.1), # fixed positions are on the log scale
      position = position_dodge(width = 0.5),
      size = 3,
      hjust = 0
    ) +
    geom_text(label = "Odds ratio (95% CrI)",
              x = 3.4,
              y = log(0.1),
              hjust = 0,
              size = 3) +
    geom_vline(xintercept = 0) +
    scale_y_continuous(
      name = "Odds ratio (95% credible interval [CrI])",
      trans = "log",
      breaks = c(0.6, 1,1.5),
      limits = c(0.1, ylim_upper)
    ) +
    scale_x_discrete(name = "") +
    scale_color_manual(
      values = c("black", "grey75"),
      labels = c("No Criteria", "Criteria")
    ) +
    labs(title = paste0("Odds ratios for ", outcome, " compared with standard oxygen therapy"),
         color = "") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    guides(color = guide_legend(reverse=TRUE)) +
    coord_flip() +

    annotate("text", x = 0.5, y = 0.9, label = "Favors Treatment", hjust = 1, size = 4, fontface = "bold") +
    annotate("text", x = 0.5, y = 1/0.9, label = "Favors SOT", hjust = 0, size = 4, fontface = "bold")
}

gemtc_forestplot2 <- function(result_criteria,
                              result_nocovar,
                              outcome,
                              ylim_upper = 2){

  rel_effect <- relative.effect(result_nocovar, t1 = "SOT")
  rel_effect_1 <-  relative.effect(result_criteria, t1 = "SOT", covariate = 1)
  rel_effect_0 <-  relative.effect(result_criteria, t1 = "SOT", covariate = 0)

  rel_effect_df <-
    bind_rows(post_draws(rel_effect, "Overall"),
              post_draws(rel_effect_1, "Criteria"),
              post_draws(rel_effect_0, "No criteria")) %>%
    mutate(Criteria = factor(Criteria,
                             levels = rev(c("Overall",
                                            "Criteria",
                                            "No criteria")),
                             ordered = T))

  summary_data <- rel_effect_df %>%
    pivot_longer(
      cols = starts_with("d."),
      names_to = "Comparison",
      values_to = "LogOdds"
    ) %>%
    group_by(Comparison, Criteria) %>%
    summarize(
      mean = exp(mean(LogOdds)),
      lower = exp(quantile(LogOdds, 0.025)),
      upper = exp(quantile(LogOdds, 0.975)),
      .groups = "drop"
    )

  comparison_labels <- c(
    "d.SOT.HFNC" = "HFNC",
    "d.SOT.CPAP" = "CPAP",
    "d.SOT.Bilevel_NIPPV" = "Bilevel NIPPV"
  )
  summary_data$Comparison <- factor(summary_data$Comparison, levels = names(comparison_labels), labels = comparison_labels)

  summary_data <- summary_data %>%
    mutate(label = sprintf("%.2f (%.2f, %.2f)", mean, lower, upper))

  ggplot(summary_data, aes(x = Comparison, y = mean, ymin = lower, ymax = upper, group = Criteria)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
    geom_pointrange(position = position_dodge(width = 0.5), size = 0.8) +
    geom_point(aes(color = Criteria),
               position = position_dodge(width = 0.5),
               size = 3) +
    geom_text(
      aes(label = label),
      y = log(0.1), # fixed positions are on the log scale
      position = position_dodge(width = 0.5),
      size = 3,
      hjust = 0
    ) +
    geom_text(label = "Odds ratio (95% CrI)",
              x = 3.4,
              y = log(0.1),
              hjust = 0,
              size = 3) +
    geom_vline(xintercept = 0) +
    scale_y_continuous(
      name = "Odds ratio (95% credible interval [CrI])",
      trans = "log",
      breaks = c(0.6, 1,1.5),
      limits = c(0.1, ylim_upper)
    ) +
    scale_x_discrete(name = "") +
    scale_color_manual(
      values = c(c_light, c_ice, "grey90")
    ) +
    labs(title = paste0("Odds ratios for ", outcome, " compared with standard oxygen therapy"),
         color = "") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    guides(color = guide_legend(reverse=TRUE)) +
    coord_flip() +

    annotate("text", x = 0.5, y = 0.9, label = "Favors Treatment", hjust = 1, size = 4, fontface = "bold") +
    annotate("text", x = 0.5, y = 1/0.9, label = "Favors SOT", hjust = 0, size = 4, fontface = "bold")
}

gemtc_forestplot3 <- function(result_criteria,
                              outcome,
                              ylim_upper = 2,
                              label1,
                              label0){

  rel_effect_1 <-  relative.effect(result_criteria, t1 = "SOT", covariate = 1)
  rel_effect_0 <-  relative.effect(result_criteria, t1 = "SOT", covariate = 0)

  rel_effect_df <-
    bind_rows(post_draws(rel_effect_1, "Yes"),
              post_draws(rel_effect_0, "No"))

  summary_data <- rel_effect_df %>%
    pivot_longer(
      cols = starts_with("d."),
      names_to = "Comparison",
      values_to = "LogOdds"
    ) %>%
    group_by(Comparison, Criteria) %>%
    summarize(
      mean = exp(mean(LogOdds)),
      lower = exp(quantile(LogOdds, 0.025)),
      upper = exp(quantile(LogOdds, 0.975)),
      .groups = "drop"
    )

  comparison_labels <- c(
    "d.SOT.HFNC" = "HFNC",
    "d.SOT.CPAP" = "CPAP",
    "d.SOT.Bilevel_NIPPV" = "Bilevel NIPPV"
  )
  summary_data$Comparison <- factor(summary_data$Comparison, levels = names(comparison_labels), labels = comparison_labels)

  summary_data$Criteria <- factor(summary_data$Criteria, labels = c("No Criteria", "Criteria"))

  summary_data <- summary_data %>%
    mutate(label = sprintf("%.2f (%.2f, %.2f)", mean, lower, upper))

  ggplot(summary_data, aes(x = Comparison, y = mean, ymin = lower, ymax = upper, group = Criteria)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
    geom_pointrange(position = position_dodge(width = 0.5), size = 0.8) +
    geom_point(aes(color = Criteria),
               position = position_dodge(width = 0.5),
               size = 3) +
    geom_text(
      aes(label = label),
      y = log(0.1), # fixed positions are on the log scale
      position = position_dodge(width = 0.5),
      size = 3,
      hjust = 0
    ) +
    geom_text(label = "Odds ratio (95% CrI)",
              x = 3.4,
              y = log(0.1),
              hjust = 0,
              size = 3) +
    geom_vline(xintercept = 0) +
    scale_y_continuous(
      name = "Odds ratio (95% credible interval [CrI])",
      trans = "log",
      breaks = c(0.6, 1,1.5),
      limits = c(0.1, ylim_upper)
    ) +
    scale_x_discrete(name = "") +
    scale_color_manual(
      values = c("black", "grey75"),
      labels = c(label0, label1)
    ) +
    labs(title = paste0("Odds ratios for ", outcome, " compared with standard oxygen therapy"),
         color = "") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    guides(color = guide_legend(reverse=TRUE)) +
    coord_flip() +

    annotate("text", x = 0.5, y = 0.9, label = "Favors Treatment", hjust = 1, size = 4, fontface = "bold") +
    annotate("text", x = 0.5, y = 1/0.9, label = "Favors SOT", hjust = 0, size = 4, fontface = "bold")
}

gemtc_forestplot4 <- function(result_criteria,
                              result_nocovar,
                              outcome,
                              ylim_upper = 2,
                              network){


  # events/N per comparison with SOT, overall and by criteria, for the labels
  arms <- network$data.ab %>%
    as.data.frame() %>%
    mutate(
      study = as.character(study),
      treatment = as.character(treatment)
    ) %>%
    group_by(study, treatment) %>%
    summarise(
      responders = sum(responders, na.rm = TRUE),
      sampleSize = sum(sampleSize, na.rm = TRUE),
      .groups = "drop"
    )

  study_pairs <- arms %>%
    group_by(study) %>%
    nest() %>%
    mutate(
      rows = map2(study, data, function(st, df) {
        trts <- sort(unique(df$treatment))
        if (length(trts) < 2) return(tibble())

        pairs <- combn(trts, 2, simplify = FALSE)

        map_dfr(pairs, function(p) {
          t1 <- p[1]; t2 <- p[2]
          a1 <- df %>% filter(treatment == t1) %>% slice(1)
          a2 <- df %>% filter(treatment == t2) %>% slice(1)

          tibble(
            Treatment_1 = t1,
            Treatment_2 = t2,
            Treatment_1_responders = a1$responders,
            Treatment_1_sampleSize = a1$sampleSize,
            Treatment_2_responders = a2$responders,
            Treatment_2_sampleSize = a2$sampleSize
          )
        })
      })
    ) %>%
    select(rows) %>%
    unnest(rows) %>%
    ungroup()

  comparison_table_summed <- study_pairs %>%
    group_by(Treatment_1, Treatment_2) %>%
    summarise(
      Treatment_1_responders = sum(Treatment_1_responders, na.rm = TRUE),
      Treatment_1_sampleSize = sum(Treatment_1_sampleSize, na.rm = TRUE),
      Treatment_2_responders = sum(Treatment_2_responders, na.rm = TRUE),
      Treatment_2_sampleSize = sum(Treatment_2_sampleSize, na.rm = TRUE),
      n_studies = dplyr::n(),
      .groups = "drop"
    ) %>%
    mutate(comparison = paste(Treatment_1, "vs", Treatment_2)) %>%
    relocate(comparison, .before = Treatment_1)

  Ntable_overall <-
    comparison_table_summed %>%
    filter(Treatment_2 == "SOT") %>%
    mutate(TreatmentLabel = paste0(Treatment_1_responders, "/", Treatment_1_sampleSize),
           SOTLabel = paste0(Treatment_2_responders, "/", Treatment_2_sampleSize),
           Criteria = "Overall"
    )%>%
    select(Treatment_1, TreatmentLabel, SOTLabel, Criteria)

  comparison_table_summed_criteria <-
    network$data.ab %>%
    as.data.frame() %>%
    mutate(
      study = as.character(study),
      treatment = as.character(treatment)
    ) %>%
    left_join(
      network$studies %>%
        mutate(study = as.character(study)),
      by = "study"
    ) %>%
    mutate(
      Criteria_label = ifelse(Criteria == 1, "Criteria", "No criteria")
    ) %>%
    group_by(study, Criteria_label, treatment) %>%
    summarise(
      responders = sum(responders, na.rm = TRUE),
      sampleSize = sum(sampleSize, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(study, Criteria_label) %>%
    nest() %>%
    mutate(
      rows = map(data, function(df) {
        trts <- sort(unique(df$treatment))
        if (length(trts) < 2) return(tibble())

        pairs <- combn(trts, 2, simplify = FALSE)

        map_dfr(pairs, function(p) {
          t1 <- p[1]; t2 <- p[2]
          a1 <- df %>% filter(treatment == t1) %>% slice(1)
          a2 <- df %>% filter(treatment == t2) %>% slice(1)

          tibble(
            Treatment_1 = t1,
            Treatment_2 = t2,
            Treatment_1_responders = a1$responders,
            Treatment_1_sampleSize = a1$sampleSize,
            Treatment_2_responders = a2$responders,
            Treatment_2_sampleSize = a2$sampleSize
          )
        })
      })
    ) %>%
    select(Criteria_label, rows) %>%
    unnest(rows) %>%
    group_by(Criteria_label, Treatment_1, Treatment_2) %>%
    summarise(
      Treatment_1_responders = sum(Treatment_1_responders, na.rm = TRUE),
      Treatment_1_sampleSize = sum(Treatment_1_sampleSize, na.rm = TRUE),
      Treatment_2_responders = sum(Treatment_2_responders, na.rm = TRUE),
      Treatment_2_sampleSize = sum(Treatment_2_sampleSize, na.rm = TRUE),
      n_studies = n(),
      .groups = "drop"
    ) %>%
    mutate(
      comparison = paste(Treatment_1, "vs", Treatment_2),
      Treatment_1_resp_n = paste0(Treatment_1_responders, "/", Treatment_1_sampleSize),
      Treatment_2_resp_n = paste0(Treatment_2_responders, "/", Treatment_2_sampleSize)
    ) %>%
    relocate(comparison, .after = Criteria_label)

  Ntable <-
    comparison_table_summed_criteria %>%
    filter(Treatment_2 == "SOT") %>%
    mutate(TreatmentLabel = paste0(Treatment_1_responders, "/", Treatment_1_sampleSize),
           SOTLabel = paste0(Treatment_2_responders, "/", Treatment_2_sampleSize),
           Criteria = Criteria_label
    )%>%
    select(Treatment_1, TreatmentLabel, SOTLabel, Criteria) %>%
    bind_rows(Ntable_overall) %>%
    mutate(Comparison = case_when(Treatment_1 == "Bilevel_NIPPV" ~ "Bilevel NIPPV",
                                  T ~ Treatment_1)) %>%
    select(-Treatment_1) %>%
    mutate(Criteria = factor(Criteria,
                             levels = rev(c("Overall",
                                            "Criteria",
                                            "No criteria")),
                             ordered = T)) %>%
    arrange(Comparison, Criteria)

  rel_effect <- relative.effect(result_nocovar, t1 = "SOT")
  rel_effect_1 <-  relative.effect(result_criteria, t1 = "SOT", covariate = 1)
  rel_effect_0 <-  relative.effect(result_criteria, t1 = "SOT", covariate = 0)

  rel_effect_df <-
    bind_rows(post_draws(rel_effect, "Overall"),
              post_draws(rel_effect_1, "Criteria"),
              post_draws(rel_effect_0, "No criteria")) %>%
    mutate(Criteria = factor(Criteria,
                             levels = rev(c("Overall",
                                            "Criteria",
                                            "No criteria")),
                             ordered = T))

  summary_data <- rel_effect_df %>%
    pivot_longer(
      cols = starts_with("d."),
      names_to = "Comparison",
      values_to = "LogOdds"
    ) %>%
    group_by(Comparison, Criteria) %>%
    summarize(
      mean = exp(mean(LogOdds)),
      lower = exp(quantile(LogOdds, 0.025)),
      upper = exp(quantile(LogOdds, 0.975)),
      .groups = "drop"
    )

  comparison_labels <- c(
    "d.SOT.HFNC" = "HFNC",
    "d.SOT.CPAP" = "CPAP",
    "d.SOT.Bilevel_NIPPV" = "Bilevel NIPPV"
  )
  summary_data$Comparison <- factor(summary_data$Comparison, levels = names(comparison_labels), labels = comparison_labels)

  # rows are matched by position, so check they line up
  stopifnot(identical(as.character(summary_data$Comparison), Ntable$Comparison),
            identical(as.character(summary_data$Criteria), as.character(Ntable$Criteria)))

  summary_data <- summary_data %>%
    mutate(label = sprintf("%.2f (%.2f, %.2f)", mean, lower, upper)) %>%
    bind_cols(select(Ntable, -Criteria, -Comparison))

  ggplot(summary_data, aes(x = Comparison, y = mean, ymin = lower, ymax = upper, group = Criteria)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
    geom_pointrange(position = position_dodge(width = 0.5), size = 0.8) +
    geom_point(aes(color = Criteria),
               position = position_dodge(width = 0.5),
               size = 3) +
    geom_text(
      aes(label = label),
      y = log(0.1), # fixed positions are on the log scale
      position = position_dodge(width = 0.5),
      size = 3,
      hjust = 0
    ) +
    geom_text(
      aes(label = TreatmentLabel),
      y = log(0.04), # fixed positions are on the log scale
      position = position_dodge(width = 0.5),
      size = 3,
      hjust = 1
    ) +
    geom_text(
      aes(label = SOTLabel),
      y = log(0.08), # fixed positions are on the log scale
      position = position_dodge(width = 0.5),
      size = 3,
      hjust = 1
    ) +
    geom_text(label = "Odds ratio (95% CrI)",
              x = 3.4,
              y = log(0.1),
              hjust = 0,
              size = 3) +
    geom_text(label = "Treatment",
              x = 3.4,
              y = log(0.04),
              hjust = 1,
              size = 3) +
    geom_text(label = "SOT",
              x = 3.4,
              y = log(0.08),
              hjust = 1,
              size = 3) +
    geom_vline(xintercept = 0) +
    scale_y_continuous(
      name = "Odds ratio (95% credible interval [CrI])",
      trans = "log",
      breaks = c(0.2, 0.6, 1,1.5),
      limits = c(0.025, ylim_upper)
    ) +
    scale_x_discrete(name = "") +
    scale_color_manual(
      values = c(c_light, c_ice, "grey90")
    ) +
    labs(title = paste0("Odds ratios for ", outcome, " compared with standard oxygen therapy"),
         color = "") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    guides(color = guide_legend(reverse=TRUE)) +
    coord_flip() +

    annotate("text", x = 0.5, y = 0.9, label = "Favors Treatment", hjust = 1, size = 4, fontface = "bold") +
    annotate("text", x = 0.5, y = 1/0.9, label = "Favors SOT", hjust = 0, size = 4, fontface = "bold")
}


# League table from posterior samples, summarised as exp(mean log OR) with
# 95% CrI (2.5th and 97.5th percentiles), to match the forest plots.
# gemtc's relative.effect.table() reports medians, so it is not used.
# Cell [row, column] = odds ratio of the COLUMN treatment vs the ROW treatment.
league_table_mean <- function(result) {
  re <- relative.effect(result, t1 = "SOT", preserve.extra = FALSE)
  s  <- as.matrix(re$samples)
  d  <- cbind(SOT = 0, s[, grep("^d\\.SOT\\.", colnames(s)), drop = FALSE])
  colnames(d) <- sub("^d\\.SOT\\.", "", colnames(d))
  trts <- colnames(d)
  out <- matrix("", length(trts), length(trts), dimnames = list(trts, trts))
  for (i in trts) for (j in trts) {
    if (i == j) { out[i, j] <- i; next }
    lor <- d[, j] - d[, i]
    out[i, j] <- sprintf("%.2f (%.2f, %.2f)",
                         exp(mean(lor)),
                         exp(quantile(lor, 0.025)),
                         exp(quantile(lor, 0.975)))
  }
  as.data.frame(out) %>% rownames_to_column("Row treatment")
}

# Save a network figure; report and skip it if plotting fails
safe_svg <- function(filename, plot_fun, width = 5, height = 6) {
  opened <- FALSE
  tryCatch({
    svg(filename = filename, width = width, height = height)
    opened <- TRUE
    plot_fun()
  },
  error = function(e) {
    message("Skipping figure ", filename, ": ", conditionMessage(e))
  },
  finally = {
    if (opened) dev.off()
  })
  invisible(NULL)
}

########### Functions to read the local extraction workbook #################

clean_names <- function(x) {
  vctrs::vec_as_names(str_squish(x), repair = "unique", quiet = TRUE)
}

# Read an outcome sheet as text and convert the count columns to numbers;
# non-numeric entries (e.g. "N/a") become NA.
read_outcomes <- function(sheet) {
  read_excel(data_file,
             sheet = sheet,
             col_types = "text",
             .name_repair = clean_names) %>%
    mutate(
      across(c(Study, Intervention, `Intervention (simplified)`), str_squish),
      across(c(`Number of participants`,
               `Mortality (n)`,
               `Number of participants - intubation outcome`,
               `Intubation (n)`),
             ~ suppressWarnings(as.numeric(.x)))
    ) %>%
    filter(!is.na(Study)) %>%
    filter(!(Study == "Added studies")) %>%
    # not in the networks: control arm was usual care (SOT, NIPPV or HFNC)
    filter(!(Study == "Arabi, 2022")) %>%
    # not in the networks: NIPPV interface and mode not specified
    filter(!(Study == "Elkheshen, 2024"))
}

read_general <- function() {
  read_excel(data_file,
             sheet = "General Data Extraction - Aggre",
             guess_max = 2000,
             .name_repair = clean_names) %>%
    mutate(Study = str_squish(Study))
}

########## Network Meta-analysis with Study-level Covariates ################

########## Load outcome data and join criteria data #########################

general_data <- read_general()

outcome_data <- read_outcomes("Outcomes - Aggregated")

# order interventions so all comparisons have the same direction
interventions <- factor(outcome_data$Intervention)
unique(interventions)
order_interventions <- rev(c("H Bilevel",
                             "FM Bilevel",
                             "H CPAP",
                             "FM CPAP",
                             "HFNC",
                             "SOT"))

interventions_simplified <- factor(outcome_data$`Intervention (simplified)`)
unique(interventions_simplified)
order_interventions_simplified <- rev(c("Bilevel NIPPV",
                                        "CPAP",
                                        "HFNC",
                                        "SOT"))

outcome_data$Intervention <- factor(outcome_data$Intervention,
                                    levels = order_interventions,
                                    labels = order_interventions,
                                    ordered = T)

outcome_data$`Intervention (simplified)` <- factor(outcome_data$`Intervention (simplified)`,
                                                   levels = order_interventions_simplified,
                                                   labels = order_interventions_simplified,
                                                   ordered = T)
outcome_data <-
  outcome_data %>%
  arrange(Study, Intervention)

# collapse arms into the four nodes (SOT, HFNC, CPAP, bilevel NIPPV)

outcome_data <-
  outcome_data %>%
  group_by(Study, `Intervention (simplified)`) %>%
  summarise(`Number of participants` = sum(`Number of participants`),
            `Mortality (n)` = sum(`Mortality (n)`),
            `Number of participants - intubation outcome` =
              sum(`Number of participants - intubation outcome`),
            `Intubation (n)` = sum(`Intubation (n)`))


# Study-level flags: Criteria (any intubation criteria) and
# Dyspnea/work of breathing (any dyspnea / work-of-breathing criterion).
# The CSV uses Table 1 study names; map them to the outcome-sheet names:
#   Lemiale (JAMA) = Lemiale, 2015a; Lemiale (Crit Care) = Lemiale, 2015b
#   Perkins is split into 2022a (CPAP vs SOT) and 2022b (HFNC vs SOT)

criteria_raw <- read_csv(criteria_file, show_col_types = FALSE)

criteria_by_study <-
  bind_rows(
    criteria_raw %>% filter(study != "Perkins, 2022"),
    criteria_raw %>% filter(study == "Perkins, 2022") %>% mutate(study = "Perkins, 2022a"),
    criteria_raw %>% filter(study == "Perkins, 2022") %>% mutate(study = "Perkins, 2022b")
  ) %>%
  mutate(study = recode(study,
                        "Lemiale, 2015 (JAMA)" = "Lemiale, 2015a",
                        "Lemiale, 2015 (Crit Care)" = "Lemiale, 2015b")) %>%
  transmute(Study = study,
            Criteria = intubation_criteria,
            `Dyspnea/work of breathing` = wob_criteria)

outcome_data_criteria <-
  outcome_data %>%
  ungroup() %>%
  left_join(criteria_by_study,
            by = "Study")

# stop if any study has no criteria coding
missing_criteria <-
  outcome_data_criteria %>%
  filter(!is.na(`Intervention (simplified)`), is.na(Criteria)) %>%
  distinct(Study)

if (nrow(missing_criteria) > 0) {
  stop("No criteria match for: ", paste(missing_criteria$Study, collapse = "; "))
}

# quick check of study-level coding
outcome_data_criteria %>%
  distinct(Study, Criteria, `Dyspnea/work of breathing`) %>%
  count(Criteria, `Dyspnea/work of breathing`)


################ Network plots ###########################

# one row per within-study comparison (for the network plots)
data_mortality <-
  outcome_data_criteria %>%
  filter(!is.na(`Mortality (n)`)) %>%
  select(Study:`Mortality (n)`,
         Criteria, `Dyspnea/work of breathing`) %>%
  rename(events = `Mortality (n)`,
         total = `Number of participants`) %>%
  group_by(Study) %>%
  # studies with only one node after classification
  filter(!(Study %in% c("Patel, 2016",
                        "Saxena, 2022"))) %>%
  reframe(
    comparison = combn(`Intervention (simplified)`, 2, simplify = FALSE),
    events1 = combn(events, 2, FUN = function(x) x[1], simplify = TRUE),
    events2 = combn(events, 2, FUN = function(x) x[2], simplify = TRUE),
    total1 = combn(total, 2, FUN = function(x) x[1], simplify = TRUE),
    total2 = combn(total, 2, FUN = function(x) x[2], simplify = TRUE),
    Criteria = max(Criteria)
  ) %>%
  mutate(
    intervention1 = map_int(comparison, ~ as.integer(.x[1])),
    intervention2 = map_int(comparison, ~ as.integer(.x[2]))
  ) %>%
  ungroup() %>%
  select(-comparison) %>%
  filter(intervention1 != intervention2) %>%
  relocate(Study, intervention1, intervention2) %>%
  mutate(intervention1 = factor(intervention1,
                                levels = 1:4,
                                labels = order_interventions_simplified),
         intervention2 = factor(intervention2,
                                levels = 1:4,
                                labels = order_interventions_simplified))


data_intubation <-
  outcome_data_criteria %>%
  filter(!is.na(`Intubation (n)`)) %>%
  select(Study, `Intervention (simplified)`, `Number of participants - intubation outcome`,
         `Intubation (n)`, Criteria, `Dyspnea/work of breathing`) %>%
  rename(events = `Intubation (n)`,
         total = `Number of participants - intubation outcome`) %>%
  group_by(Study) %>%
  # studies with only one node after classification
  filter(!(Study %in% c("Patel, 2016",
                        "Saxena, 2022"))) %>%
  reframe(
    comparison = combn(`Intervention (simplified)`, 2, simplify = FALSE),
    events1 = combn(events, 2, FUN = function(x) x[1], simplify = TRUE),
    events2 = combn(events, 2, FUN = function(x) x[2], simplify = TRUE),
    total1 = combn(total, 2, FUN = function(x) x[1], simplify = TRUE),
    total2 = combn(total, 2, FUN = function(x) x[2], simplify = TRUE),
    Criteria = max(Criteria)
  ) %>%
  mutate(
    intervention1 = map_int(comparison, ~ as.integer(.x[1])),
    intervention2 = map_int(comparison, ~ as.integer(.x[2]))
  ) %>%
  ungroup() %>%
  select(-comparison) %>%
  filter(intervention1 != intervention2) %>%
  relocate(Study, intervention1, intervention2) %>%
  mutate(intervention1 = factor(intervention1,
                                levels = 1:4,
                                labels = order_interventions_simplified),
         intervention2 = factor(intervention2,
                                levels = 1:4,
                                labels = order_interventions_simplified))


plot_network <- function(df_pairs, title){

  edge_counts <- df_pairs %>%
    count(intervention1, intervention2, Criteria)

  g <- graph_from_data_frame(edge_counts, directed = FALSE)

  E(g)$weight <- edge_counts$n

  E(g)$color <- ifelse(edge_counts$Criteria == 1,
                       c_ice,
                       c_light)

  E(g)$label.x <- c(0.5,
                    0.3,
                    -0.5,
                    -0.5,
                    0.3,
                    0.5,
                    -0.5,
                    -0.1,
                    0.1)
  E(g)$label.y <- c(0.7,0.5,
                    0.2, -0.2,
                    -0.5, -0.7,
                    0.5, -0.5,
                    -0.5)

  plot(
    g,
    layout = layout_in_circle,
    vertex.color = "grey90",
    vertex.size = 40,
    vertex.label.color = "black",
    vertex.label.cex = 1.2,
    edge.width = 5+E(g)$weight,
    edge.label = E(g)$weight,
    edge.alpha = 0.5,
    edge.label.color = "black",
    edge.label.cex = 1,
    main = title
  )
}

safe_svg("Figures/Network_Mortality_Criteria.svg",
         function() plot_network(data_mortality, "Mortality"))

safe_svg("Figures/Network_Intubation_Criteria.svg",
         function() plot_network(data_intubation, "Intubation"))

# Networks for work-of-breathing criteria

# one row per within-study comparison (for the network plots)
data_mortality <-
  outcome_data_criteria %>%
  filter(!is.na(`Mortality (n)`)) %>%
  select(Study:`Mortality (n)`,
         Criteria, `Dyspnea/work of breathing`) %>%
  rename(events = `Mortality (n)`,
         total = `Number of participants`) %>%
  group_by(Study) %>%
  filter(!(Study %in% c("Patel, 2016",
                        "Saxena, 2022"))) %>%
  reframe(
    comparison = combn(`Intervention (simplified)`, 2, simplify = FALSE),
    events1 = combn(events, 2, FUN = function(x) x[1], simplify = TRUE),
    events2 = combn(events, 2, FUN = function(x) x[2], simplify = TRUE),
    total1 = combn(total, 2, FUN = function(x) x[1], simplify = TRUE),
    total2 = combn(total, 2, FUN = function(x) x[2], simplify = TRUE),
    Criteria = max(`Dyspnea/work of breathing`)
  ) %>%
  mutate(
    intervention1 = map_int(comparison, ~ as.integer(.x[1])),
    intervention2 = map_int(comparison, ~ as.integer(.x[2]))
  ) %>%
  ungroup() %>%
  select(-comparison) %>%
  filter(intervention1 != intervention2) %>%
  relocate(Study, intervention1, intervention2) %>%
  mutate(intervention1 = factor(intervention1,
                                levels = 1:4,
                                labels = order_interventions_simplified),
         intervention2 = factor(intervention2,
                                levels = 1:4,
                                labels = order_interventions_simplified))


data_intubation <-
  outcome_data_criteria %>%
  filter(!is.na(`Intubation (n)`)) %>%
  select(Study, `Intervention (simplified)`, `Number of participants - intubation outcome`,
         `Intubation (n)`, Criteria, `Dyspnea/work of breathing`) %>%
  rename(events = `Intubation (n)`,
         total = `Number of participants - intubation outcome`) %>%
  group_by(Study) %>%
  filter(!(Study %in% c("Patel, 2016",
                        "Saxena, 2022"))) %>%
  reframe(
    comparison = combn(`Intervention (simplified)`, 2, simplify = FALSE),
    events1 = combn(events, 2, FUN = function(x) x[1], simplify = TRUE),
    events2 = combn(events, 2, FUN = function(x) x[2], simplify = TRUE),
    total1 = combn(total, 2, FUN = function(x) x[1], simplify = TRUE),
    total2 = combn(total, 2, FUN = function(x) x[2], simplify = TRUE),
    Criteria = max(`Dyspnea/work of breathing`)
  ) %>%
  mutate(
    intervention1 = map_int(comparison, ~ as.integer(.x[1])),
    intervention2 = map_int(comparison, ~ as.integer(.x[2]))
  ) %>%
  ungroup() %>%
  select(-comparison) %>%
  filter(intervention1 != intervention2) %>%
  relocate(Study, intervention1, intervention2) %>%
  mutate(intervention1 = factor(intervention1,
                                levels = 1:4,
                                labels = order_interventions_simplified),
         intervention2 = factor(intervention2,
                                levels = 1:4,
                                labels = order_interventions_simplified))

plot_network_wob <- function(df_pairs, title){

  edge_counts <- df_pairs %>%
    count(intervention1, intervention2, Criteria)

  g <- graph_from_data_frame(edge_counts, directed = FALSE)

  E(g)$weight <- edge_counts$n

  E(g)$color <- ifelse(edge_counts$Criteria == 1,
                       c_ice,
                       c_light)

  E(g)$label.x <- c(0.5,
                    0.3,
                    -0.5,
                    -0.5,
                    0.3,
                    0.5,
                    -0.5,
                    -0.3,
                    -0.1,
                    0.1)
  E(g)$label.y <- c(0.7,0.5,
                    0.2, -0.2,
                    -0.5, -0.7,
                    0.7,
                    0.5, -0.5,
                    -0.5)

  plot(
    g,
    layout = layout_in_circle,
    vertex.color = "grey90",
    vertex.size = 40,
    vertex.label.color = "black",
    vertex.label.cex = 1.2,
    edge.width = 5+E(g)$weight,
    edge.label = E(g)$weight,
    edge.alpha = 0.5,
    edge.label.color = "black",
    edge.label.cex = 1,
    main = title
  )
}

safe_svg("Figures/Network_Mortality_WOBCriteria.svg",
         function() plot_network_wob(data_mortality, "Mortality, work-of-breathing criteria"))

safe_svg("Figures/Network_Intubation_WOBCriteria.svg",
         function() plot_network_wob(data_intubation, "Intubation, work-of-breathing criteria"))


#################### Meta-Analysis, Pairwise and Network ######################

# Heterogeneity priors: empirical log-normal priors from Turner et al. 2012
# (doi: 10.1093/ije/dys041), non-pharmacological interventions.
# JAGS dlnorm takes a precision, so the second parameter is 1/sd^2.
# gemtc's default priors are used for all other parameters.

# semi-objective outcome: log-normal(mean -2.89, sd 1.92)
hy.prior.ett <- mtc.hy.prior(type="std.dev", distr="dlnorm", -2.89, 0.27)
# all-cause mortality: log-normal(mean -3.93, sd 1.51)
hy.prior.mort <- mtc.hy.prior(type="std.dev", distr="dlnorm", -3.93, 0.44)

############# Data - Intubation, any criteria #######################


df_ett <-
  outcome_data_criteria %>%
  filter(!is.na(`Intubation (n)`)) %>%
  select(Study, `Intubation (n)`,
         `Number of participants - intubation outcome`,
         `Intervention (simplified)`,
         Criteria, `Dyspnea/work of breathing`) %>%
  rename(events = `Intubation (n)`,
         treatment = `Intervention (simplified)`,
         sample_size = `Number of participants - intubation outcome`) %>%
  select(Study, events, sample_size, treatment, Criteria)

df_ett %>% ungroup() %>% summarise(sum(sample_size))

df_ett %>% ungroup() %>% distinct(Study) %>% print(n = 100)

df_ett %>% print(n=100)

# combine arms that fall in the same node after simplification
df_mtc <- df_ett %>%
  dplyr::rename(
    study = Study,
    responders = events,
    sampleSize = sample_size,
    treatment = treatment) %>%
  mutate(Criteria = Criteria == 1,
         treatment = dplyr::recode(treatment, "Bilevel NIPPV" = "Bilevel_NIPPV")) %>%
  dplyr::group_by(study, treatment) %>%
  dplyr::summarise(
    responders = sum(responders),
    sampleSize = sum(sampleSize),
    .groups = "drop"
  )

# make study-level covariate table
studies <-
  df_ett %>%
  group_by(Study) %>%
  summarise(Criteria = max(Criteria)) %>%
  rename(study = Study) %>%
  filter(!(study %in% c(
    "Patel, 2016",
    "Saxena, 2022"
  )))


regressor <- list(coefficient = 'shared',
                  variable = 'Criteria',
                  control = 'SOT')


network <- mtc.network(data.ab = df_mtc,
                       studies = studies)

# how many patients?
sum(network$data.ab$sampleSize)


############## Pairwise meta-analysis, Intubation outcome ###############


data_by_comparison <- list()

i <- 1

for (treat1 in network$treatments$id) {
  for (treat2 in network$treatments$id) {
    if (treat1 < treat2) {
      df_tmp <-
        network$data %>%
        dplyr::filter(treatment %in% c(treat1, treat2)) %>%
        dplyr::mutate(
          treatment = ifelse(treatment == treat1, "treat1", "treat2")
        ) %>%
        tidyr::pivot_wider(
          names_from  = treatment,
          values_from = responders:sampleSize
        ) %>%
        dplyr::filter(
          !is.na(responders_treat1) & !is.na(responders_treat2)
        ) %>%
        dplyr::rename(
          event1 = responders_treat1,
          n1     = sampleSize_treat1,
          event2 = responders_treat2,
          n2     = sampleSize_treat2
        ) %>%
        dplyr::left_join(network$studies, by = "study") %>%
        dplyr::mutate(
          Criteria = ifelse(Criteria == 1, "Yes", "No")
        )

      if (nrow(df_tmp) > 0) {
        data_by_comparison[[i]] <- list(
          comparison = c(treat1 = treat1, treat2 = treat2),
          data       = df_tmp
        )
        i <- i + 1
      }
    }
  }
}

results      <- list()
loo_results  <- list()
loo_summaries <- list()

for (i in seq_along(data_by_comparison)) {

  comp_name <- paste(
    data_by_comparison[[i]]$comparison[1],
    data_by_comparison[[i]]$comparison[2],
    sep = "_vs_"
  )

  dat_i <- data_by_comparison[[i]]$data

  res <-
    metabin(
      event1, n1, event2, n2,
      data    = dat_i,
      studlab = study,
      sm      = "OR",
      method  = "Inverse",
      random  = TRUE,
      subgroup = dat_i$Criteria
    )

  results[[comp_name]] <- res

  meta::forest(
    res,
    label.e    = data_by_comparison[[i]]$comparison[1],
    label.c    = data_by_comparison[[i]]$comparison[2],
    label.left = paste0("Favours ",
                        data_by_comparison[[i]]$comparison[1]),
    label.right = paste0("Favours ",
                         data_by_comparison[[i]]$comparison[2]),
    file  = paste0(
      "Figures/pairwise_forest_intubation_",
      data_by_comparison[[i]]$comparison[1],
      "_vs_",
      data_by_comparison[[i]]$comparison[2],
      ".svg"
    ),
    width = 10,
    text.addline1 = "Outcome = Intubation"
  )

    # leave-one-out (random effects), when there are at least 2 studies
  if (res$k > 1) {

    loo <- metainf(res, pooled = "random")
    loo_results[[comp_name]] <- loo

    n_loo <- length(loo$studlab)

    loo_df <- data.frame(
      comparison     = rep(comp_name, n_loo),
      study_left_out = loo$studlab,
      TE             = loo$TE,
      seTE           = loo$seTE,
      lower          = loo$lower,
      upper          = loo$upper
    )

    loo_df$OR    <- exp(loo_df$TE)
    loo_df$OR_l  <- exp(loo_df$lower)
    loo_df$OR_u  <- exp(loo_df$upper)

    loo_summaries[[comp_name]] <- loo_df

    meta::forest(
      loo,
      sortvar = loo$TE,
      leftcols = "studlab",
      leftlabs = "Study omitted",
      lab.e = data_by_comparison[[i]]$comparison[1],
      lab.c = data_by_comparison[[i]]$comparison[2],
      file = paste0(
        "Figures/LOO_intubation_",
        data_by_comparison[[i]]$comparison[1],
        "_vs_",
        data_by_comparison[[i]]$comparison[2],
        ".svg"
      ),
      width = 10,
      text.addline1 = "Leave-one-out random-effects OR"
    )
  }
}

loo_all <- dplyr::bind_rows(loo_summaries)
print(loo_all)

write_csv(loo_all, "Intubation_LOO.csv")

# Pairwise league table: cell [row, column] = OR of the row treatment vs the
# column treatment (the opposite orientation to league_table_mean()).

all_treatments <- sort(unique(unlist(lapply(data_by_comparison, function(x) x$comparison))))

league_mat <- matrix(NA, nrow = length(all_treatments), ncol = length(all_treatments),
                     dimnames = list(all_treatments, all_treatments))

for (name in names(results)) {
  split_name <- strsplit(name, "_vs_")[[1]]
  t1 <- split_name[1]
  t2 <- split_name[2]
  res <- results[[name]]

  est <- res$TE.random
  lb  <- res$lower.random
  ub  <- res$upper.random

  or   <- exp(est)
  lb_or <- exp(lb)
  ub_or <- exp(ub)

  cell_value <- sprintf("%.2f [%.2f, %.2f]", or, lb_or, ub_or)

  league_mat[t1, t2] <- cell_value
  league_mat[t2, t1] <- sprintf("%.2f [%.2f, %.2f]", 1/or, 1/ub_or, 1/lb_or)
}

league_table_df <- as.data.frame(league_mat) %>%
  rownames_to_column("Treatment")

ft <- flextable(league_table_df)
ft <- autofit(ft)
ft <- set_caption(ft, caption = "League Table of Pairwise Odds Ratios for Intubation with 95% CIs")

doc <- read_docx()
doc <- body_add_flextable(doc, ft)
print(doc, target = "intubation_league_table.docx")

############## Bayesian network meta-analysis - Intubation #############

# model without covariate
model_nocovar <- mtc_model_n01(
  network,
  type="consistency",
  hy.prior = hy.prior.ett,
  likelihood = "binom",
  link = "logit")

result_nocovar <- mtc.run(model_nocovar)

saveRDS(result_nocovar,
        file = "intubation_nocovar_nma.rds")

result_nocovar <- readRDS("intubation_nocovar_nma.rds")


model <- mtc_model_n01(
  network,
  type="regression",
  regressor=regressor,
  hy.prior = hy.prior.ett,
  likelihood = "binom",
  link = "logit")

result_criteria <- mtc.run(model)

saveRDS(result_criteria,
        file = "intubation_criteria_nma.rds")

result_criteria <- readRDS("intubation_criteria_nma.rds")


summary(result_criteria)

rank_df1 <- rank_probs(relative.effect(result_criteria, t1 = "SOT",
                                       covariate = 1)$samples)

rank_df0 <- rank_probs(relative.effect(result_criteria, t1 = "SOT",
                                       covariate = 0)$samples)


rankplot(rank_df1, rank_df0, title = "Treatment rankings: intubation")
ggsave("Figures/intubation_treatment_rankings.svg",
       width = 7, height = 5)

rank_df1
rank_df0


gemtc_forestplot2(result_criteria,
                  result_nocovar,
                  "intubation")
ggsave("Figures/Forest_Intubation_Criteria.svg",
       width = 7, height = 4, units = "in")

gemtc_forestplot4(result_criteria,
                  result_nocovar,
                  "intubation",
                  network = network)

ggsave("Figures/Forest_Intubation_Criteria4.svg",
       width = 7, height = 4, units = "in")


# league table (no-covariate model)
league_table_mean(result_nocovar) %>%
  write_csv("gemtc_intubation_overall_league_table.csv")


# Check node-splitting
nodesplit <- mtc_nodesplit_n01(network,
                           hy.prior = hy.prior.ett)

saveRDS(nodesplit,
        file = "nodesplit_criteria_intubation.rds")

summary(readRDS("nodesplit_criteria_intubation.rds"))

nodesplit_ett <- readRDS("nodesplit_criteria_intubation.rds")

nodesplit_forestplot(nodesplit_ett, "Intubation")
ggsave("Figures/nodesplitting_forest_intubation.svg",
       width = 6, height = 6)

######### Sensitivity - Intubation: coefficient = exchangeable ###########

# The main model uses a shared criteria coefficient: the same effect of criteria
# on every comparison with SOT. Here the coefficients are exchangeable (drawn
# from a common normal distribution); below, they are unrelated (estimated
# separately for each comparison).


regressor <- list(coefficient = 'exchangeable',
                  variable = 'Criteria',
                  control = 'SOT')

model <- mtc_model_n01(
  network,
  type="regression",
  regressor=regressor,
  hy.prior = hy.prior.ett,
  likelihood = "binom",
  link = "logit")

result_criteria_exch <- mtc.run(model)

saveRDS(result_criteria_exch,
        file = "intubation_criteria_nma_exchangeable.rds")

result_criteria_exch <- readRDS("intubation_criteria_nma_exchangeable.rds")

rank_df1_exch <- rank_probs(relative.effect(result_criteria_exch, t1 = "SOT",
                                            covariate = 1)$samples)

rank_df0_exch <- rank_probs(relative.effect(result_criteria_exch, t1 = "SOT",
                                            covariate = 0)$samples)


rankplot(rank_df1_exch, rank_df0_exch, title = "Treatment rankings: intubation, exchangeable coefficient")
ggsave("Figures/intubation_treatment_rankings_exch.svg",
       width = 7, height = 5)


gemtc_forestplot(result_criteria_exch,
                 "intubation")
ggsave("Figures/Forest_Intubation_Criteria_exchangeable.svg",
       width = 7, height = 4, units = "in")


######### Sensitivity - Intubation: coefficient = unrelated ###########


regressor <- list(coefficient = 'unrelated',
                  variable = 'Criteria',
                  control = 'SOT')

model <- mtc_model_n01(
  network,
  type="regression",
  regressor=regressor,
  hy.prior = hy.prior.ett,
  likelihood = "binom",
  link = "logit")

result_criteria_unr <- mtc.run(model)

saveRDS(result_criteria_unr,
        file = "intubation_criteria_nma_unrelated.rds")

result_criteria_unr <- readRDS("intubation_criteria_nma_unrelated.rds")

rank_df1_unr <- rank_probs(relative.effect(result_criteria_unr, t1 = "SOT",
                                           covariate = 1)$samples)

rank_df0_unr <- rank_probs(relative.effect(result_criteria_unr, t1 = "SOT",
                                           covariate = 0)$samples)
rankplot(rank_df1_unr, rank_df0_unr, title = "Treatment rankings: intubation, unrelated coefficient")
ggsave("Figures/intubation_treatment_rankings_unr.svg",
       width = 7, height = 5)

gemtc_forestplot(result_criteria_unr,
                 "intubation",
                 ylim_upper = 2.55)
ggsave("Figures/Forest_Intubation_Criteria_unrelated.svg",
       width = 7, height = 4, units = "in")


############# Data - Mortality, any criteria #######################

df_mort <-
  outcome_data_criteria %>%
  filter(!is.na(`Mortality (n)`)) %>%
  select(Study, `Mortality (n)`,
         `Number of participants`,
         `Intervention (simplified)`,
         Criteria, `Dyspnea/work of breathing`) %>%
  rename(events = `Mortality (n)`,
         treatment = `Intervention (simplified)`,
         sample_size = `Number of participants`) %>%
  select(Study, events, sample_size, treatment, Criteria)

# combine arms that fall in the same node after simplification
df_mtc <- df_mort %>%
  dplyr::rename(
    study = Study,
    responders = events,
    sampleSize = sample_size,
    treatment = treatment) %>%
  mutate(Criteria = Criteria == 1,
         treatment = dplyr::recode(treatment, "Bilevel NIPPV" = "Bilevel_NIPPV")) %>%
  dplyr::group_by(study, treatment) %>%
  dplyr::summarise(
    responders = sum(responders),
    sampleSize = sum(sampleSize),
    .groups = "drop"
  )

# make study-level covariate table
studies <-
  df_mort %>%
  group_by(Study) %>%
  summarise(Criteria = max(Criteria)) %>%
  rename(study = Study) %>%
  filter(!(study %in% c(
    "Patel, 2016",
    "Saxena, 2022"
  )))


network <- mtc.network(data.ab = df_mtc,
                       studies = studies)

sum(network$data.ab$sampleSize)


regressor <- list(coefficient = 'shared',
                  variable = 'Criteria',
                  control = 'SOT')


############## Pairwise meta-analysis, Mortality outcome ###############


data_by_comparison <- list()

i = 1

for (treat1 in network$treatments$id){
  for (treat2 in network$treatments$id){
    if (treat1 < treat2){
      df_tmp <-
        network$data %>%
        filter(treatment %in% c(treat1, treat2)) %>%
        mutate(treatment = ifelse(treatment == treat1, "treat1",
                                  "treat2")) %>%
        pivot_wider(names_from = treatment,
                    values_from = responders:sampleSize) %>%
        filter(!is.na(responders_treat1) & !is.na(responders_treat2)) %>%
        rename(event1 = responders_treat1,
               n1 = sampleSize_treat1,
               event2 = responders_treat2,
               n2 = sampleSize_treat2) %>%
        left_join(network$studies,
                  by = "study") %>%
        mutate(Criteria = ifelse(Criteria == 1,
                                 "Yes",
                                 "No"))

      if(nrow(df_tmp) > 0){
        data_by_comparison[[i]] <- list(
          comparison = c(treat1 = treat1, treat2 = treat2),
          data = df_tmp
        )

        i = i+1
      }
    }
  }
}

results      <- list()
loo_results  <- list()
loo_summaries <- list()

for (i in seq_along(data_by_comparison)) {

  comp_name <- paste(
    data_by_comparison[[i]]$comparison[1],
    data_by_comparison[[i]]$comparison[2],
    sep = "_vs_"
  )

  dat_i <- data_by_comparison[[i]]$data

  res <-
    metabin(
      event1, n1, event2, n2,
      data    = dat_i,
      studlab = study,
      sm      = "OR",
      method  = "Inverse",
      random  = TRUE,
      subgroup = dat_i$Criteria
    )

  results[[comp_name]] <- res

  meta::forest(
    res,
    label.e    = data_by_comparison[[i]]$comparison[1],
    label.c    = data_by_comparison[[i]]$comparison[2],
    label.left = paste0("Favours ",
                        data_by_comparison[[i]]$comparison[1]),
    label.right = paste0("Favours ",
                         data_by_comparison[[i]]$comparison[2]),
    file  = paste0(
      "Figures/pairwise_forest_mortality_",
      data_by_comparison[[i]]$comparison[1],
      "_vs_",
      data_by_comparison[[i]]$comparison[2],
      ".svg"
    ),
    width = 10,
    text.addline1 = "Outcome = Mortality"
  )

    # leave-one-out (random effects), when there are at least 2 studies
  if (res$k > 1) {

    loo <- metainf(res, pooled = "random")
    loo_results[[comp_name]] <- loo

    n_loo <- length(loo$studlab)

    loo_df <- data.frame(
      comparison     = rep(comp_name, n_loo),
      study_left_out = loo$studlab,
      TE             = loo$TE,
      seTE           = loo$seTE,
      lower          = loo$lower,
      upper          = loo$upper
    )

    loo_df$OR    <- exp(loo_df$TE)
    loo_df$OR_l  <- exp(loo_df$lower)
    loo_df$OR_u  <- exp(loo_df$upper)

    loo_summaries[[comp_name]] <- loo_df

    meta::forest(
      loo,
      sortvar = loo$TE,
      leftcols = "studlab",
      leftlabs = "Study omitted",
      lab.e = data_by_comparison[[i]]$comparison[1],
      lab.c = data_by_comparison[[i]]$comparison[2],
      file = paste0(
        "Figures/LOO_mortality_",
        data_by_comparison[[i]]$comparison[1],
        "_vs_",
        data_by_comparison[[i]]$comparison[2],
        ".svg"
      ),
      width = 10,
      text.addline1 = "Leave-one-out random-effects OR"
    )
  }
}

loo_all <- dplyr::bind_rows(loo_summaries)
print(loo_all)

write_csv(loo_all, "Mortality_LOO.csv")

# Pairwise league table: cell [row, column] = OR of the row treatment vs the
# column treatment (the opposite orientation to league_table_mean()).

all_treatments <- sort(unique(unlist(lapply(data_by_comparison, function(x) x$comparison))))

league_mat <- matrix(NA, nrow = length(all_treatments), ncol = length(all_treatments),
                     dimnames = list(all_treatments, all_treatments))

for (name in names(results)) {
  split_name <- strsplit(name, "_vs_")[[1]]
  t1 <- split_name[1]
  t2 <- split_name[2]
  res <- results[[name]]

  est <- res$TE.random
  lb  <- res$lower.random
  ub  <- res$upper.random

  or   <- exp(est)
  lb_or <- exp(lb)
  ub_or <- exp(ub)

  cell_value <- sprintf("%.2f [%.2f, %.2f]", or, lb_or, ub_or)

  league_mat[t1, t2] <- cell_value
  league_mat[t2, t1] <- sprintf("%.2f [%.2f, %.2f]", 1/or, 1/ub_or, 1/lb_or)
}

league_table_df <- as.data.frame(league_mat) %>%
  rownames_to_column("Treatment")

ft <- flextable(league_table_df)
ft <- autofit(ft)
ft <- set_caption(ft, caption = "League Table of Pairwise Odds Ratios for Mortality with 95% CIs")

doc <- read_docx()
doc <- body_add_flextable(doc, ft)
print(doc, target = "mortality_league_table.docx")


############## Bayesian network meta-analysis - Mortality #############


# model without covariate
model_nocovar <- mtc_model_n01(
  network,
  type="consistency",
  hy.prior = hy.prior.mort,
  likelihood = "binom",
  link = "logit")

result_nocovar <- mtc.run(model_nocovar)

saveRDS(result_nocovar,
        file = "mortality_nocovar_nma.rds")

result_nocovar <- readRDS("mortality_nocovar_nma.rds")


model <- mtc_model_n01(
  network,
  type="regression",
  hy.prior = hy.prior.mort,
  regressor=regressor,
  likelihood = "binom",
  link = "logit")

result_criteria <- mtc.run(model)

saveRDS(result_criteria,
        file = "mortality_criteria_nma.rds")

result_criteria <- readRDS("mortality_criteria_nma.rds")

rank_df1 <- rank_probs(relative.effect(result_criteria, t1 = "SOT",
                                       covariate = 1)$samples)

rank_df0 <- rank_probs(relative.effect(result_criteria, t1 = "SOT",
                                       covariate = 0)$samples)
rankplot(rank_df1, rank_df0, title = "Treatment rankings: mortality")
ggsave("Figures/mortality_treatment_rankings.svg",
       width = 7, height = 5)


summary(result_criteria)

gemtc_forestplot2(result_criteria,
                  result_nocovar,
                  "mortality")
ggsave("Figures/Forest_Mortality_Criteria.svg",
       width = 7, height = 4, units = "in")

gemtc_forestplot4(result_criteria,
                  result_nocovar,
                  "mortality",
                  network = network)

ggsave("Figures/Forest_Mortality_Criteria4.svg",
       width = 7, height = 4, units = "in")


# league table (no-covariate model)
league_table_mean(result_nocovar) %>%
  write_csv("gemtc_mortality_overall_league_table.csv")


# Check node-splitting
nodesplit <- mtc_nodesplit_n01(network,
                           hy.prior = hy.prior.mort)

saveRDS(nodesplit,
        file = "nodesplit_criteria_mortality.rds")

summary(readRDS("nodesplit_criteria_mortality.rds"))

nodesplit_mort <- readRDS("nodesplit_criteria_mortality.rds")

nodesplit_forestplot(nodesplit_mort, "Mortality")
ggsave("Figures/nodesplitting_forest_mortality.svg",
       width = 6, height = 6)

######### Sensitivity - Mortality: coefficient = exchangeable ###########

# The main model uses a shared criteria coefficient: the same effect of criteria
# on every comparison with SOT. Here the coefficients are exchangeable (drawn
# from a common normal distribution); below, they are unrelated (estimated
# separately for each comparison).


regressor <- list(coefficient = 'exchangeable',
                  variable = 'Criteria',
                  control = 'SOT')

model <- mtc_model_n01(
  network,
  type="regression",
  regressor=regressor,
  hy.prior = hy.prior.mort,
  likelihood = "binom",
  link = "logit")

result_criteria_exch <- mtc.run(model)


saveRDS(result_criteria_exch,
        file = "mortality_criteria_nma_exchangeable.rds")

result_criteria_exch <- readRDS("mortality_criteria_nma_exchangeable.rds")

rank_df1_exch <- rank_probs(relative.effect(result_criteria_exch, t1 = "SOT",
                                            covariate = 1)$samples)

rank_df0_exch <- rank_probs(relative.effect(result_criteria_exch, t1 = "SOT",
                                            covariate = 0)$samples)
rankplot(rank_df1_exch, rank_df0_exch, title = "Treatment rankings: mortality, exchangeable coefficient")
ggsave("Figures/mortality_treatment_rankings_exchangeable.svg",
       width = 7, height = 5)


gemtc_forestplot(result_criteria_exch,
                 "mortality",
                 ylim_upper = 2.7)
ggsave("Figures/Forest_Mortality_Criteria_exchangeable.svg",
       width = 7, height = 4, units = "in")


######### Sensitivity - Mortality: coefficient = unrelated ###########


regressor <- list(coefficient = 'unrelated',
                  variable = 'Criteria',
                  control = 'SOT')

model <- mtc_model_n01(
  network,
  type="regression",
  regressor=regressor,
  hy.prior = hy.prior.mort,
  likelihood = "binom",
  link = "logit")

result_criteria_unr <- mtc.run(model)

saveRDS(result_criteria_unr,
        file = "mortality_criteria_nma_unrelated.rds")

result_criteria_unr <- readRDS("mortality_criteria_nma_unrelated.rds")

rank_df1_unr <- rank_probs(relative.effect(result_criteria_unr, t1 = "SOT",
                                           covariate = 1)$samples)

rank_df0_unr <- rank_probs(relative.effect(result_criteria_unr, t1 = "SOT",
                                           covariate = 0)$samples)
rankplot(rank_df1_unr, rank_df0_unr, title = "Treatment rankings: mortality, unrelated coefficient")
ggsave("Figures/mortality_treatment_rankings_unrelated.svg",
       width = 7, height = 5)


gemtc_forestplot(result_criteria_unr,
                 "mortality",
                 ylim_upper = 3.3)
ggsave("Figures/Forest_Mortality_Criteria_unrelated.svg",
       width = 7, height = 4, units = "in")


############# Model - Intubation, work-of-breathing criteria #######################

df_ett <-
  outcome_data_criteria %>%
  filter(!is.na(`Intubation (n)`)) %>%
  select(Study, `Intubation (n)`,
         `Number of participants - intubation outcome`,
         `Intervention (simplified)`,
         Criteria, `Dyspnea/work of breathing`) %>%
  rename(events = `Intubation (n)`,
         treatment = `Intervention (simplified)`,
         sample_size = `Number of participants - intubation outcome`) %>%
  select(Study, events, sample_size, treatment, `Dyspnea/work of breathing`) %>%
  rename(Criteria = `Dyspnea/work of breathing`)

# combine arms that fall in the same node after simplification
df_mtc <- df_ett %>%
  dplyr::rename(
    study = Study,
    responders = events,
    sampleSize = sample_size,
    treatment = treatment) %>%
  mutate(Criteria = Criteria == 1,
         treatment = dplyr::recode(treatment, "Bilevel NIPPV" = "Bilevel_NIPPV")) %>%
  dplyr::group_by(study, treatment) %>%
  dplyr::summarise(
    responders = sum(responders),
    sampleSize = sum(sampleSize),
    .groups = "drop"
  )

# make study-level covariate table
studies <-
  df_ett %>%
  group_by(Study) %>%
  summarise(Criteria = max(Criteria)) %>%
  rename(study = Study) %>%
  filter(!(study %in% c(
    "Patel, 2016",
    "Saxena, 2022"
  )))


network <- mtc.network(data.ab = df_mtc,
                       studies = studies)


regressor <- list(coefficient = 'shared',
                  variable = 'Criteria',
                  control = 'SOT')


model <- mtc_model_n01(
  network,
  type="regression",
  hy.prior = hy.prior.ett,
  regressor=regressor,
  likelihood = "binom",
  link = "logit")

result_criteria <- mtc.run(model)

saveRDS(result_criteria,
        file = "intubation_wobcriteria_nma.rds")

result_criteria <- readRDS("intubation_wobcriteria_nma.rds")

summary(result_criteria)


gemtc_forestplot(result_criteria,
                 "intubation") +
  labs(title = "Odds ratios for intubation compared with SOT, work-of-breathing criteria")
ggsave("Figures/Forest_Intubation_WOBCriteria.svg",
       width = 7, height = 4, units = "in")


############# Model - Mortality, work-of-breathing criteria #######################

df_mort <-
  outcome_data_criteria %>%
  filter(!is.na(`Mortality (n)`)) %>%
  select(Study, `Mortality (n)`,
         `Number of participants`,
         `Intervention (simplified)`,
         Criteria, `Dyspnea/work of breathing`) %>%
  rename(events = `Mortality (n)`,
         treatment = `Intervention (simplified)`,
         sample_size = `Number of participants`) %>%
  select(Study, events, sample_size, treatment, `Dyspnea/work of breathing`) %>%
  rename(Criteria = `Dyspnea/work of breathing`)

# combine arms that fall in the same node after simplification
df_mtc <- df_mort %>%
  dplyr::rename(
    study = Study,
    responders = events,
    sampleSize = sample_size,
    treatment = treatment) %>%
  mutate(Criteria = Criteria == 1,
         treatment = dplyr::recode(treatment, "Bilevel NIPPV" = "Bilevel_NIPPV")) %>%
  dplyr::group_by(study, treatment) %>%
  dplyr::summarise(
    responders = sum(responders),
    sampleSize = sum(sampleSize),
    .groups = "drop"
  )

# make study-level covariate table
studies <-
  df_mort %>%
  group_by(Study) %>%
  summarise(Criteria = max(Criteria)) %>%
  rename(study = Study) %>%
  filter(!(study %in% c(
    "Patel, 2016",
    "Saxena, 2022"
  )))


network <- mtc.network(data.ab = df_mtc,
                       studies = studies)


regressor <- list(coefficient = 'shared',
                  variable = 'Criteria',
                  control = 'SOT')

model <- mtc_model_n01(
  network,
  type="regression",
  hy.prior = hy.prior.mort,
  regressor=regressor,
  likelihood = "binom",
  link = "logit")

result_criteria <- mtc.run(model)

saveRDS(result_criteria,
        file = "mortality_wobcriteria_nma.rds")

summary(result_criteria)

result_criteria <- readRDS("mortality_wobcriteria_nma.rds")

gemtc_forestplot(result_criteria,
                 "mortality") +
  labs(title = "Odds ratios for mortality compared with SOT, work-of-breathing criteria")
ggsave("Figures/Forest_Mortality_WOBCriteria.svg",
       width = 7, height = 4, units = "in")


############# Model - Mortality, only criteria studies #######################

# post hoc, because of incoherence in the mortality network

df_mort <-
  outcome_data_criteria %>%
  filter(!is.na(`Mortality (n)`)) %>%
  select(Study, `Mortality (n)`,
         `Number of participants`,
         `Intervention (simplified)`,
         Criteria, `Dyspnea/work of breathing`) %>%
  rename(events = `Mortality (n)`,
         treatment = `Intervention (simplified)`,
         sample_size = `Number of participants`) %>%
  select(Study, events, sample_size, treatment, Criteria) %>%
  filter(Criteria == 1) %>%
  select(-Criteria)

# combine arms that fall in the same node after simplification
df_mtc <- df_mort %>%
  dplyr::rename(
    study = Study,
    responders = events,
    sampleSize = sample_size,
    treatment = treatment) %>%
  mutate(treatment = dplyr::recode(treatment, "Bilevel NIPPV" = "Bilevel_NIPPV")) %>%
  dplyr::group_by(study, treatment) %>%
  dplyr::summarise(
    responders = sum(responders),
    sampleSize = sum(sampleSize),
    .groups = "drop"
  )


network <- mtc.network(data.ab = df_mtc)

nodesplit_sens <- mtc_nodesplit_n01(network,
                                hy.prior = hy.prior.mort)

saveRDS(nodesplit_sens,
        file = "nodesplit_mortality_onlycriteriastudies.rds")

nodesplit_sens <- readRDS("nodesplit_mortality_onlycriteriastudies.rds")

summary(readRDS("nodesplit_mortality_onlycriteriastudies.rds"))

nodesplit_forestplot(nodesplit_sens, "Mortality, Criteria Studies")
ggsave("Figures/nodesplitting_forest_mortality_onlycriteriastudies.svg",
       width = 6, height = 6)


############# Data - Intubation, COVID/Immunocompromised/ward #######################

general_data <-
  general_data %>%
  mutate(Study =
           case_when(Study == "Lemiale, 2015 (JAMA)" ~ "Lemiale, 2015a",
                     Study == "Lemiale, 2015 (Crit. Care)" ~ "Lemiale, 2015b",
                     TRUE ~ Study))

outcome_data <- read_outcomes("Outcomes - Aggregated, with Mai")

# order interventions so all comparisons have the same direction
interventions <- factor(outcome_data$Intervention)
unique(interventions)
order_interventions <- rev(c("H Bilevel",
                             "FM Bilevel",
                             "H CPAP",
                             "FM CPAP",
                             "HFNC",
                             "SOT"))

interventions_simplified <- factor(outcome_data$`Intervention (simplified)`)
unique(interventions_simplified)
order_interventions_simplified <- rev(c("Bilevel NIPPV",
                                        "CPAP",
                                        "HFNC",
                                        "SOT"))

outcome_data$Intervention <- factor(outcome_data$Intervention,
                                    levels = order_interventions,
                                    labels = order_interventions,
                                    ordered = T)

outcome_data$`Intervention (simplified)` <- factor(outcome_data$`Intervention (simplified)`,
                                                   levels = order_interventions_simplified,
                                                   labels = order_interventions_simplified,
                                                   ordered = T)
outcome_data <-
  outcome_data %>%
  arrange(Study, Intervention)

# collapse arms into the four nodes (SOT, HFNC, CPAP, bilevel NIPPV)

outcome_data <-
  outcome_data %>%
  group_by(Study, `Intervention (simplified)`) %>%
  summarise(`Number of participants` = sum(`Number of participants`),
            `Mortality (n)` = sum(`Mortality (n)`),
            `Number of participants - intubation outcome` =
              sum(`Number of participants - intubation outcome`),
            `Intubation (n)` = sum(`Intubation (n)`))


df_ett <-
  outcome_data %>%
  filter(!is.na(`Intervention (simplified)`)) %>%
  left_join(
    general_data %>%
      distinct(Study,
               `COVID (1 = 100% COVID, 0 = everything else)`,
               `Immunocompromised (1 = 100% immunocompromised, 0 = everything else)`,
               `Setting (1 = 100% ward, 0 = everything else)`) %>%
      filter(!is.na(Study)),
    by = "Study"
  ) %>%
  rename(COVID19 = `COVID (1 = 100% COVID, 0 = everything else)`,
         Immunocompromised = `Immunocompromised (1 = 100% immunocompromised, 0 = everything else)`,
         Ward = `Setting (1 = 100% ward, 0 = everything else)`) %>%
  filter(!is.na(`Intubation (n)`)) %>%
  select(Study, `Intubation (n)`,
         `Number of participants - intubation outcome`,
         `Intervention (simplified)`,
         COVID19, Immunocompromised, Ward) %>%
  rename(events = `Intubation (n)`,
         treatment = `Intervention (simplified)`,
         sample_size = `Number of participants - intubation outcome`) %>%
  select(Study, events, sample_size, treatment, COVID19, Immunocompromised, Ward)

missing_covariates <-
  df_ett %>%
  filter(is.na(COVID19) | is.na(Immunocompromised) | is.na(Ward)) %>%
  distinct(Study)

if (nrow(missing_covariates) > 0) {
  stop("No COVID / immunocompromised / ward coding for: ",
       paste(missing_covariates$Study, collapse = "; "))
}

df_ett %>% ungroup() %>% summarise(sum(sample_size))

df_ett %>% ungroup() %>% distinct(Study) %>% print(n = 100)

df_ett %>% print(n=100)

# combine arms that fall in the same node after simplification
df_mtc <- df_ett %>%
  dplyr::rename(
    study = Study,
    responders = events,
    sampleSize = sample_size,
    treatment = treatment) %>%
  mutate(treatment = dplyr::recode(treatment, "Bilevel NIPPV" = "Bilevel_NIPPV")) %>%
  dplyr::group_by(study, treatment) %>%
  dplyr::summarise(
    responders = sum(responders),
    sampleSize = sum(sampleSize),
    .groups = "drop"
  )

# make study-level covariate table
studies <-
  df_ett %>%
  distinct(Study, COVID19, Immunocompromised, Ward) %>%
  mutate(COVID19 = COVID19 == 1,
         Immunocompromised = Immunocompromised == 1,
         Ward = Ward==1) %>%
  rename(study = Study) %>%
  filter(!(study %in% c(
    "Patel, 2016",
    "Saxena, 2022"
  )))


regressorCOVID <- list(coefficient = 'shared',
                       variable = 'COVID19',
                       control = 'SOT')
regressorImmuno <- list(coefficient = 'shared',
                        variable = 'Immunocompromised',
                        control = 'SOT')
regressorWard <- list(coefficient = 'shared',
                      variable = 'Ward',
                      control = 'SOT')


network <- mtc.network(data.ab = df_mtc,
                       studies = studies)


#### Bayesian network meta-regression - COVID; Intubation #############

plot(network)


modelCOVID <- mtc_model_n01(
  network,
  type="regression",
  regressor=regressorCOVID,
  hy.prior = hy.prior.ett,
  likelihood = "binom",
  link = "logit")

resultCOVID <- mtc.run(modelCOVID)

saveRDS(resultCOVID,
        file = "intubation_covid_nma.rds")

modelImmuno <- mtc_model_n01(
  network,
  type="regression",
  regressor=regressorImmuno,
  hy.prior = hy.prior.ett,
  likelihood = "binom",
  link = "logit")

resultImmuno <- mtc.run(modelImmuno)

saveRDS(resultImmuno,
        file = "intubation_immuno_nma.rds")

modelWard <- mtc_model_n01(
  network,
  type="regression",
  regressor=regressorWard,
  hy.prior = hy.prior.ett,
  likelihood = "binom",
  link = "logit")

resultWard <- mtc.run(modelWard)

saveRDS(resultWard,
        file = "intubation_ward_nma.rds")


resultCOVID <- readRDS("intubation_covid_nma.rds")
resultImmuno <- readRDS("intubation_immuno_nma.rds")
resultWard <- readRDS("intubation_ward_nma.rds")

summary(resultCOVID)

gemtc_forestplot3(resultCOVID,
                  "intubation",
                  label1 = "COVID-19",
                  label0 = "Not COVID-19")
ggsave("Figures/Forest_Intubation_COVID.svg",
       width = 7, height = 4, units = "in")

summary(resultImmuno)

gemtc_forestplot3(resultImmuno,
                  "intubation",
                  label1 = "Immunocompromised",
                  label0 = "Not immunocompromised")
ggsave("Figures/Forest_Intubation_Immuno.svg",
       width = 7, height = 4, units = "in")

summary(resultWard)

gemtc_forestplot3(resultWard,
                  "intubation",
                  label1 = "Ward-based",
                  label0 = "Not ward-based")
ggsave("Figures/Forest_Intubation_Ward.svg",
       width = 7, height = 4, units = "in")

############# Data - Mortality, COVID/Immunocompromised/ward #######################

general_data <-
  general_data %>%
  mutate(Study =
           case_when(Study == "Lemiale, 2015 (JAMA)" ~ "Lemiale, 2015a",
                     Study == "Lemiale, 2015 (Crit. Care)" ~ "Lemiale, 2015b",
                     TRUE ~ Study))

outcome_data <- read_outcomes("Outcomes - Aggregated, with Mai")

# order interventions so all comparisons have the same direction
interventions <- factor(outcome_data$Intervention)
unique(interventions)
order_interventions <- rev(c("H Bilevel",
                             "FM Bilevel",
                             "H CPAP",
                             "FM CPAP",
                             "HFNC",
                             "SOT"))

interventions_simplified <- factor(outcome_data$`Intervention (simplified)`)
unique(interventions_simplified)
order_interventions_simplified <- rev(c("Bilevel NIPPV",
                                        "CPAP",
                                        "HFNC",
                                        "SOT"))

outcome_data$Intervention <- factor(outcome_data$Intervention,
                                    levels = order_interventions,
                                    labels = order_interventions,
                                    ordered = T)

outcome_data$`Intervention (simplified)` <- factor(outcome_data$`Intervention (simplified)`,
                                                   levels = order_interventions_simplified,
                                                   labels = order_interventions_simplified,
                                                   ordered = T)
outcome_data <-
  outcome_data %>%
  arrange(Study, Intervention)

# collapse arms into the four nodes (SOT, HFNC, CPAP, bilevel NIPPV)

outcome_data <-
  outcome_data %>%
  group_by(Study, `Intervention (simplified)`) %>%
  summarise(`Number of participants` = sum(`Number of participants`),
            `Mortality (n)` = sum(`Mortality (n)`),
            `Number of participants - intubation outcome` =
              sum(`Number of participants - intubation outcome`),
            `Intubation (n)` = sum(`Intubation (n)`))


df_mort <-
  outcome_data %>%
  filter(!is.na(`Mortality (n)`)) %>%
  select(Study, `Mortality (n)`,
         `Number of participants`,
         `Intervention (simplified)`) %>%
  left_join(
    general_data %>%
      distinct(Study,
               `COVID (1 = 100% COVID, 0 = everything else)`,
               `Immunocompromised (1 = 100% immunocompromised, 0 = everything else)`,
               `Setting (1 = 100% ward, 0 = everything else)`) %>%
      filter(!is.na(Study)),
    by = "Study"
  ) %>%
  rename(COVID19 = `COVID (1 = 100% COVID, 0 = everything else)`,
         Immunocompromised = `Immunocompromised (1 = 100% immunocompromised, 0 = everything else)`,
         Ward = `Setting (1 = 100% ward, 0 = everything else)`) %>%
  filter(!is.na(`Mortality (n)`)) %>%
  select(Study, `Mortality (n)`,
         `Number of participants`,
         `Intervention (simplified)`,
         COVID19, Immunocompromised, Ward) %>%
  rename(events = `Mortality (n)`,
         treatment = `Intervention (simplified)`,
         sample_size = `Number of participants`) %>%
  select(Study, events, sample_size, treatment, COVID19, Immunocompromised, Ward)

missing_covariates <-
  df_mort %>%
  filter(is.na(COVID19) | is.na(Immunocompromised) | is.na(Ward)) %>%
  distinct(Study)

if (nrow(missing_covariates) > 0) {
  stop("No COVID / immunocompromised / ward coding for: ",
       paste(missing_covariates$Study, collapse = "; "))
}

# combine arms that fall in the same node after simplification
df_mtc <- df_mort %>%
  dplyr::rename(
    study = Study,
    responders = events,
    sampleSize = sample_size,
    treatment = treatment) %>%
  mutate(treatment = dplyr::recode(treatment, "Bilevel NIPPV" = "Bilevel_NIPPV")) %>%
  dplyr::group_by(study, treatment) %>%
  dplyr::summarise(
    responders = sum(responders),
    sampleSize = sum(sampleSize),
    .groups = "drop"
  )

# make study-level covariate table
studies <-
  df_mort %>%
  distinct(Study, COVID19, Immunocompromised, Ward) %>%
  mutate(COVID19 = COVID19 == 1,
         Immunocompromised = Immunocompromised == 1,
         Ward = Ward==1) %>%
  rename(study = Study) %>%
  filter(!(study %in% c(
    "Patel, 2016",
    "Saxena, 2022"
  )))


regressorCOVID <- list(coefficient = 'shared',
                       variable = 'COVID19',
                       control = 'SOT')
regressorImmuno <- list(coefficient = 'shared',
                        variable = 'Immunocompromised',
                        control = 'SOT')
regressorWard <- list(coefficient = 'shared',
                      variable = 'Ward',
                      control = 'SOT')


network <- mtc.network(data.ab = df_mtc,
                       studies = studies)


#### Bayesian network meta-regression - COVID/Immuno/Ward; Mortality #############

plot(network)


modelCOVID <- mtc_model_n01(
  network,
  type="regression",
  regressor=regressorCOVID,
  hy.prior = hy.prior.mort,
  likelihood = "binom",
  link = "logit")

resultCOVID <- mtc.run(modelCOVID)

saveRDS(resultCOVID,
        file = "mortality_covid_nma.rds")

modelImmuno <- mtc_model_n01(
  network,
  type="regression",
  regressor=regressorImmuno,
  hy.prior = hy.prior.ett,
  likelihood = "binom",
  link = "logit")

resultImmuno <- mtc.run(modelImmuno)

saveRDS(resultImmuno,
        file = "mortality_immuno_nma.rds")

modelWard <- mtc_model_n01(
  network,
  type="regression",
  regressor=regressorWard,
  hy.prior = hy.prior.ett,
  likelihood = "binom",
  link = "logit")

resultWard <- mtc.run(modelWard)

saveRDS(resultWard,
        file = "mortality_ward_nma.rds")


resultCOVID <- readRDS("mortality_covid_nma.rds")
resultImmuno <- readRDS("mortality_immuno_nma.rds")
resultWard <- readRDS("mortality_ward_nma.rds")

summary(resultCOVID)

gemtc_forestplot3(resultCOVID,
                  "mortality",
                  label1 = "COVID-19",
                  label0 = "Not COVID-19")
ggsave("Figures/Forest_Mortality_COVID.svg",
       width = 7, height = 4, units = "in")

summary(resultImmuno)

gemtc_forestplot3(resultImmuno,
                  "mortality",
                  label1 = "Immunocompromised",
                  label0 = "Not immunocompromised")
ggsave("Figures/Forest_Mortality_Immuno.svg",
       width = 7, height = 4, units = "in")

summary(resultWard)

gemtc_forestplot3(resultWard,
                  "mortality",
                  label1 = "Ward-based",
                  label0 = "Not ward-based")
ggsave("Figures/Forest_Mortality_Ward.svg",
       width = 7, height = 4, units = "in")


##### COVID and Immuno Networks ############


# one row per within-study comparison (for the network plots)
data_mortality <-
  outcome_data %>%
  filter(!is.na(`Mortality (n)`)) %>%
  select(Study, `Mortality (n)`,
         `Number of participants`,
         `Intervention (simplified)`) %>%
  left_join(
    general_data %>%
      distinct(Study,
               `COVID (1 = 100% COVID, 0 = everything else)`,
               `Immunocompromised (1 = 100% immunocompromised, 0 = everything else)`,
               `Setting (1 = 100% ward, 0 = everything else)`) %>%
      filter(!is.na(Study)),
    by = "Study"
  ) %>%
  rename(COVID19 = `COVID (1 = 100% COVID, 0 = everything else)`,
         Immunocompromised = `Immunocompromised (1 = 100% immunocompromised, 0 = everything else)`,
         Ward = `Setting (1 = 100% ward, 0 = everything else)`) %>%
  filter(!is.na(`Mortality (n)`)) %>%
  rename(events = `Mortality (n)`,
         total = `Number of participants`) %>%
  group_by(Study) %>%
  # studies with only one node after classification
  filter(!(Study %in% c("Patel, 2016",
                        "Saxena, 2022"))) %>%
  reframe(
    comparison = combn(`Intervention (simplified)`, 2, simplify = FALSE),
    events1 = combn(events, 2, FUN = function(x) x[1], simplify = TRUE),
    events2 = combn(events, 2, FUN = function(x) x[2], simplify = TRUE),
    total1 = combn(total, 2, FUN = function(x) x[1], simplify = TRUE),
    total2 = combn(total, 2, FUN = function(x) x[2], simplify = TRUE),
    COVID19 = max(COVID19),
    Immunocompromised = max(Immunocompromised),
    Ward = max(Ward)
  ) %>%
  mutate(
    intervention1 = map_int(comparison, ~ as.integer(.x[1])),
    intervention2 = map_int(comparison, ~ as.integer(.x[2]))
  ) %>%
  ungroup() %>%
  select(-comparison) %>%
  filter(intervention1 != intervention2) %>%
  relocate(Study, intervention1, intervention2) %>%
  mutate(intervention1 = factor(intervention1,
                                levels = 1:4,
                                labels = order_interventions_simplified),
         intervention2 = factor(intervention2,
                                levels = 1:4,
                                labels = order_interventions_simplified))


data_intubation <-
  outcome_data %>%
  filter(!is.na(`Intervention (simplified)`)) %>%
  left_join(
    general_data %>%
      distinct(Study,
               `COVID (1 = 100% COVID, 0 = everything else)`,
               `Immunocompromised (1 = 100% immunocompromised, 0 = everything else)`,
               `Setting (1 = 100% ward, 0 = everything else)`) %>%
      filter(!is.na(Study)),
    by = "Study"
  ) %>%
  rename(COVID19 = `COVID (1 = 100% COVID, 0 = everything else)`,
         Immunocompromised = `Immunocompromised (1 = 100% immunocompromised, 0 = everything else)`,
         Ward = `Setting (1 = 100% ward, 0 = everything else)`) %>%
  filter(!is.na(`Intubation (n)`)) %>%
  select(Study, `Intubation (n)`,
         `Number of participants - intubation outcome`,
         `Intervention (simplified)`,
         COVID19, Immunocompromised, Ward) %>%
  rename(events = `Intubation (n)`,
         total = `Number of participants - intubation outcome`) %>%
  group_by(Study) %>%
  # studies with only one node after classification
  filter(!(Study %in% c("Patel, 2016",
                        "Saxena, 2022"))) %>%
  reframe(
    comparison = combn(`Intervention (simplified)`, 2, simplify = FALSE),
    events1 = combn(events, 2, FUN = function(x) x[1], simplify = TRUE),
    events2 = combn(events, 2, FUN = function(x) x[2], simplify = TRUE),
    total1 = combn(total, 2, FUN = function(x) x[1], simplify = TRUE),
    total2 = combn(total, 2, FUN = function(x) x[2], simplify = TRUE),
    COVID19 = max(COVID19),
    Immunocompromised = max(Immunocompromised),
    Ward = max(Ward)
  ) %>%
  mutate(
    intervention1 = map_int(comparison, ~ as.integer(.x[1])),
    intervention2 = map_int(comparison, ~ as.integer(.x[2]))
  ) %>%
  ungroup() %>%
  select(-comparison) %>%
  filter(intervention1 != intervention2) %>%
  relocate(Study, intervention1, intervention2) %>%
  mutate(intervention1 = factor(intervention1,
                                levels = 1:4,
                                labels = order_interventions_simplified),
         intervention2 = factor(intervention2,
                                levels = 1:4,
                                labels = order_interventions_simplified))


plot_network_covid = function(df_pairs, title){
  edge_counts <- df_pairs %>%
    count(intervention1, intervention2, Criteria)

  g <- graph_from_data_frame(edge_counts, directed = FALSE)

  E(g)$weight <- edge_counts$n

  E(g)$color <- ifelse(edge_counts$Criteria == 1,
                       c_ice,
                       c_light)

  E(g)$label.x <- c(0.6,
                    0.4,
                    0,
                    0,
                    0.5,
                    -0.6,
                    -0.4,
                    -0.2,
                    0.2)
  E(g)$label.y <- c(0.6,
                    0.4,
                    0.2,
                    -0.2,
                    -0.5,
                    0.6,
                    0.4,
                    0,
                    0)

  plot(
    g,
    layout = layout_in_circle,
    vertex.color = "grey90",
    vertex.size = 40,
    vertex.label.color = "black",
    vertex.label.cex = 1.2,
    edge.width = 5+E(g)$weight,
    edge.label = E(g)$weight,
    edge.alpha = 0.5,
    edge.label.color = "black",
    edge.label.cex = 1,
    main = title
  )
}

safe_svg("Figures/Network_Mortality_COVID.svg",
         function() plot_network_covid(data_mortality %>% rename(Criteria = COVID19), "Mortality - COVID"))

safe_svg("Figures/Network_Intubation_COVID19.svg",
         function() plot_network_covid(data_intubation %>% rename(Criteria = COVID19), "Intubation - COVID"))

plot_network_immuno = function(df_pairs, title){
  edge_counts <- df_pairs %>%
    count(intervention1, intervention2, Criteria)

  g <- graph_from_data_frame(edge_counts, directed = FALSE)

  E(g)$weight <- edge_counts$n

  E(g)$color <- ifelse(edge_counts$Criteria == 1,
                       c_ice,
                       c_light)

  E(g)$label.x <- c(0.6,
                    0.4,
                    0,
                    0,
                    0.4,
                    0.6,
                    -0.5,
                    -0.2,
                    0.2)
  E(g)$label.y <- c(0.6,
                    0.4,
                    0.2,
                    -0.2,
                    -0.4,
                    -0.6,
                    0.5,
                    0,
                    0)

  plot(
    g,
    layout = layout_in_circle,
    vertex.color = "grey90",
    vertex.size = 40,
    vertex.label.color = "black",
    vertex.label.cex = 1.2,
    edge.width = 5+E(g)$weight,
    edge.label = E(g)$weight,
    edge.alpha = 0.5,
    edge.label.color = "black",
    edge.label.cex = 1,
    main = title
  )
}

safe_svg("Figures/Network_Mortality_Immuno.svg",
         function() plot_network_immuno(data_mortality %>% rename(Criteria = Immunocompromised), "Mortality - Immunocompromised"))

safe_svg("Figures/Network_Intubation_Immuno.svg",
         function() plot_network_immuno(data_intubation %>% rename(Criteria = Immunocompromised), "Intubation - Immunocompromised"))


################################################################################
############ POST-HOC: Interface network (SOT / HFNC / Facemask NIV / Helmet NIV)
################################################################################
# Splits NIPPV by interface rather than mode (CPAP and bilevel combined within
# each interface). Run last: it redefines outcome_data and
# order_interventions_simplified.

gemtc_forestplot2_interface <- function(result_criteria,
                                        result_nocovar,
                                        outcome,
                                        ylim_upper = 2,
                                        ylim_lower = 0.1){

  rel_effect <- relative.effect(result_nocovar, t1 = "SOT")
  rel_effect_1 <-  relative.effect(result_criteria, t1 = "SOT", covariate = 1)
  rel_effect_0 <-  relative.effect(result_criteria, t1 = "SOT", covariate = 0)

  rel_effect_df <-
    bind_rows(post_draws(rel_effect, "Overall"),
              post_draws(rel_effect_1, "Criteria"),
              post_draws(rel_effect_0, "No criteria")) %>%
    mutate(Criteria = factor(Criteria,
                             levels = rev(c("Overall",
                                            "Criteria",
                                            "No criteria")),
                             ordered = T))

  summary_data <- rel_effect_df %>%
    pivot_longer(
      cols = starts_with("d."),
      names_to = "Comparison",
      values_to = "LogOdds"
    ) %>%
    group_by(Comparison, Criteria) %>%
    summarize(
      mean = exp(mean(LogOdds)),
      lower = exp(quantile(LogOdds, 0.025)),
      upper = exp(quantile(LogOdds, 0.975)),
      .groups = "drop"
    )

  comparison_labels <- c(
    "d.SOT.HFNC" = "HFNC",
    "d.SOT.Facemask_NIV" = "Facemask NIV",
    "d.SOT.Helmet_NIV" = "Helmet NIV"
  )
  summary_data$Comparison <- factor(summary_data$Comparison, levels = names(comparison_labels), labels = comparison_labels)

  summary_data <- summary_data %>%
    mutate(label = sprintf("%.2f (%.2f, %.2f)", mean, lower, upper))

  ggplot(summary_data, aes(x = Comparison, y = mean, ymin = lower, ymax = upper, group = Criteria)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
    geom_pointrange(position = position_dodge(width = 0.5), size = 0.8) +
    geom_point(aes(color = Criteria),
               position = position_dodge(width = 0.5),
               size = 3) +
    geom_text(
      aes(label = label),
      y = log(ylim_lower), # fixed positions are on the log scale
      position = position_dodge(width = 0.5),
      size = 3,
      hjust = 0
    ) +
    geom_text(label = "Odds ratio (95% CrI)",
              x = 3.4,
              y = log(ylim_lower),
              hjust = 0,
              size = 3) +
    geom_vline(xintercept = 0) +
    scale_y_continuous(
      name = "Odds ratio (95% credible interval [CrI])",
      trans = "log",
      breaks = c(0.1,0.3, 1),
      limits = c(ylim_lower, ylim_upper)
    ) +
    scale_x_discrete(name = "") +
    scale_color_manual(
      values = c(c_light, c_ice, "grey90")
    ) +
    labs(title = paste0("Odds ratios for ", outcome, " compared with standard oxygen therapy"),
         color = "") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    guides(color = guide_legend(reverse=TRUE)) +
    coord_flip() +

    annotate("text", x = 0.5, y = 0.9, label = "Favors Treatment", hjust = 1, size = 4, fontface = "bold") +
    annotate("text", x = 0.5, y = 1/0.9, label = "Favors SOT", hjust = 0, size = 4, fontface = "bold")
}

########## Load outcome data (interface network) and join criteria data ##########

outcome_data <- read_outcomes("Outcomes - Aggregated")

order_interventions <- rev(c("H Bilevel",
                             "FM Bilevel",
                             "H CPAP",
                             "FM CPAP",
                             "HFNC",
                             "SOT"))

# classify by interface; anything else (e.g. Nair's mixed "Helmet/ FM Bilevel")
# becomes NA and drops out
outcome_data <-
  outcome_data %>%
  mutate(`Intervention (simplified)` =
           case_when(Intervention %in% c("H Bilevel","H CPAP") ~ "Helmet NIV",
                     Intervention %in% c("FM Bilevel", "FM CPAP")~ "Facemask NIV",
                     Intervention == "HFNC" ~ "HFNC",
                     Intervention == "SOT" ~ "SOT"))

interventions_simplified <- factor(outcome_data$`Intervention (simplified)`)
unique(interventions_simplified)
order_interventions_simplified <- rev(c("Helmet NIV",
                                        "Facemask NIV",
                                        "HFNC",
                                        "SOT"))

outcome_data$Intervention <- factor(outcome_data$Intervention,
                                    levels = order_interventions,
                                    labels = order_interventions,
                                    ordered = T)

outcome_data$`Intervention (simplified)` <- factor(outcome_data$`Intervention (simplified)`,
                                                   levels = order_interventions_simplified,
                                                   labels = order_interventions_simplified,
                                                   ordered = T)
outcome_data <-
  outcome_data %>%
  arrange(Study, Intervention)

outcome_data <-
  outcome_data %>%
  group_by(Study, `Intervention (simplified)`) %>%
  summarise(`Number of participants` = sum(`Number of participants`),
            `Mortality (n)` = sum(`Mortality (n)`),
            `Number of participants - intubation outcome` =
              sum(`Number of participants - intubation outcome`),
            `Intubation (n)` = sum(`Intubation (n)`))

# same study-level criteria flags as the main analysis (criteria_by_study from above)
outcome_data_criteria <-
  outcome_data %>%
  ungroup() %>%
  left_join(criteria_by_study,
            by = "Study")

missing_criteria <-
  outcome_data_criteria %>%
  filter(!is.na(`Intervention (simplified)`), is.na(Criteria)) %>%
  distinct(Study)

if (nrow(missing_criteria) > 0) {
  stop("No criteria match for: ", paste(missing_criteria$Study, collapse = "; "))
}

################ Network plots (interface network) ###########################

data_mortality <-
  outcome_data_criteria %>%
  filter(!is.na(`Mortality (n)`)) %>%
  select(Study:`Mortality (n)`,
         Criteria, `Dyspnea/work of breathing`) %>%
  rename(events = `Mortality (n)`,
         total = `Number of participants`) %>%
  group_by(Study) %>%
  # studies with only one node after classification
  filter(!(Study %in% c("Nair, 2021"))) %>%
  reframe(
    comparison = combn(`Intervention (simplified)`, 2, simplify = FALSE),
    events1 = combn(events, 2, FUN = function(x) x[1], simplify = TRUE),
    events2 = combn(events, 2, FUN = function(x) x[2], simplify = TRUE),
    total1 = combn(total, 2, FUN = function(x) x[1], simplify = TRUE),
    total2 = combn(total, 2, FUN = function(x) x[2], simplify = TRUE),
    Criteria = max(Criteria)
  ) %>%
  mutate(
    intervention1 = map_int(comparison, ~ as.integer(.x[1])),
    intervention2 = map_int(comparison, ~ as.integer(.x[2]))
  ) %>%
  ungroup() %>%
  select(-comparison) %>%
  filter(intervention1 != intervention2) %>%
  relocate(Study, intervention1, intervention2) %>%
  mutate(intervention1 = factor(intervention1,
                                levels = 1:4,
                                labels = order_interventions_simplified),
         intervention2 = factor(intervention2,
                                levels = 1:4,
                                labels = order_interventions_simplified))

data_intubation <-
  outcome_data_criteria %>%
  filter(!is.na(`Intubation (n)`)) %>%
  select(Study, `Intervention (simplified)`, `Number of participants - intubation outcome`,
         `Intubation (n)`, Criteria, `Dyspnea/work of breathing`) %>%
  rename(events = `Intubation (n)`,
         total = `Number of participants - intubation outcome`) %>%
  group_by(Study) %>%
  filter(!(Study %in% c("Nair, 2021"))) %>%
  reframe(
    comparison = combn(`Intervention (simplified)`, 2, simplify = FALSE),
    events1 = combn(events, 2, FUN = function(x) x[1], simplify = TRUE),
    events2 = combn(events, 2, FUN = function(x) x[2], simplify = TRUE),
    total1 = combn(total, 2, FUN = function(x) x[1], simplify = TRUE),
    total2 = combn(total, 2, FUN = function(x) x[2], simplify = TRUE),
    Criteria = max(Criteria)
  ) %>%
  mutate(
    intervention1 = map_int(comparison, ~ as.integer(.x[1])),
    intervention2 = map_int(comparison, ~ as.integer(.x[2]))
  ) %>%
  ungroup() %>%
  select(-comparison) %>%
  filter(intervention1 != intervention2) %>%
  relocate(Study, intervention1, intervention2) %>%
  mutate(intervention1 = factor(intervention1,
                                levels = 1:4,
                                labels = order_interventions_simplified),
         intervention2 = factor(intervention2,
                                levels = 1:4,
                                labels = order_interventions_simplified))

plot_network_interface <- function(df_pairs, title){

  edge_counts <- df_pairs %>%
    count(intervention1, intervention2, Criteria)

  g <- graph_from_data_frame(edge_counts, directed = FALSE)

  E(g)$weight <- edge_counts$n

  E(g)$color <- ifelse(edge_counts$Criteria == 1,
                       c_ice,
                       c_light)

  E(g)$label.y <- c(0.6,
                    0.4,
                    0.2,
                    -0.2,
                    -0.5,
                    0.6,
                    0.4,
                    0,
                    -0.5)
  E(g)$label.x <- c(0.6,
                    0.4,
                    0.2,
                    0.2,
                    0.5,
                    -0.6,
                    -0.4,
                    0,
                    -0.5)

  plot(
    g,
    layout = layout_in_circle,
    vertex.color = "grey90",
    vertex.size = 40,
    vertex.label.color = "black",
    vertex.label.cex = 1.2,
    edge.width = 5+E(g)$weight,
    edge.label = E(g)$weight,
    edge.alpha = 0.5,
    edge.label.color = "black",
    edge.label.cex = 1,
    main = title
  )
}

safe_svg("Figures/Network2_Mortality_Criteria.svg",
         function() plot_network_interface(data_mortality, "Mortality"))

safe_svg("Figures/Network2_Intubation_Criteria.svg",
         function() plot_network_interface(data_intubation, "Intubation"))

############# Data - Intubation, interface network #######################

df_ett <-
  outcome_data_criteria %>%
  filter(!is.na(`Intubation (n)`)) %>%
  select(Study, `Intubation (n)`,
         `Number of participants - intubation outcome`,
         `Intervention (simplified)`,
         Criteria) %>%
  rename(events = `Intubation (n)`,
         treatment = `Intervention (simplified)`,
         sample_size = `Number of participants - intubation outcome`) %>%
  select(Study, events, sample_size, treatment, Criteria)

df_ett %>% ungroup() %>% summarise(sum(sample_size))

df_ett %>% ungroup() %>% distinct(Study) %>% print(n = 100)

df_mtc <- df_ett %>%
  dplyr::rename(
    study = Study,
    responders = events,
    sampleSize = sample_size,
    treatment = treatment) %>%
  mutate(Criteria = Criteria == 1,
         treatment = dplyr::recode(treatment, "Helmet NIV" = "Helmet_NIV")) %>%
  mutate(treatment = dplyr::recode(treatment, "Facemask NIV" = "Facemask_NIV")) %>%
  dplyr::group_by(study, treatment) %>%
  dplyr::summarise(
    responders = sum(responders),
    sampleSize = sum(sampleSize),
    .groups = "drop"
  ) %>%
  filter(!(study == "Nair, 2021"))

studies <-
  df_ett %>%
  group_by(Study) %>%
  summarise(Criteria = max(Criteria)) %>%
  rename(study = Study) %>%
  filter(!(study %in% c(
    "Nair, 2021")))

regressor <- list(coefficient = 'shared',
                  variable = 'Criteria',
                  control = 'SOT')

network <- mtc.network(data.ab = df_mtc,
                       studies = studies)

sum(network$data.ab$sampleSize)

############## Bayesian network meta-analysis - Intubation, interface network #############

model_nocovar <- mtc_model_n01(
  network,
  type="consistency",
  hy.prior = hy.prior.ett,
  likelihood = "binom",
  link = "logit")

result_nocovar <- mtc.run(model_nocovar)

saveRDS(result_nocovar,
        file = "network2_intubation_nocovar_nma.rds")

result_nocovar <- readRDS("network2_intubation_nocovar_nma.rds")

model <- mtc_model_n01(
  network,
  type="regression",
  regressor=regressor,
  hy.prior = hy.prior.ett,
  likelihood = "binom",
  link = "logit")

result_criteria <- mtc.run(model)

saveRDS(result_criteria,
        file = "network2_intubation_criteria_nma.rds")

result_criteria <- readRDS("network2_intubation_criteria_nma.rds")

summary(result_criteria)

gemtc_forestplot2_interface(result_criteria,
                            result_nocovar,
                            "intubation",
                            ylim_upper = 2,
                            ylim_lower = 0.01)
ggsave("Figures/network2_Forest_Intubation_Criteria.svg",
       width = 7, height = 4, units = "in")

# league table (no-covariate model)
league_table_mean(result_nocovar) %>%
  write_csv("network2_gemtc_intubation_overall_league_table.csv")


############# Data - Mortality, interface network #######################

df_mort <-
  outcome_data_criteria %>%
  filter(!is.na(`Mortality (n)`)) %>%
  select(Study, `Mortality (n)`,
         `Number of participants`,
         `Intervention (simplified)`,
         Criteria) %>%
  rename(events = `Mortality (n)`,
         treatment = `Intervention (simplified)`,
         sample_size = `Number of participants`) %>%
  select(Study, events, sample_size, treatment, Criteria)

df_mtc <-
  df_mort %>%
  dplyr::rename(
    study = Study,
    responders = events,
    sampleSize = sample_size,
    treatment = treatment) %>%
  mutate(Criteria = Criteria == 1,
         treatment = dplyr::recode(treatment, "Helmet NIV" = "Helmet_NIV")) %>%
  mutate(treatment = dplyr::recode(treatment, "Facemask NIV" = "Facemask_NIV")) %>%
  dplyr::group_by(study, treatment) %>%
  dplyr::summarise(
    responders = sum(responders),
    sampleSize = sum(sampleSize),
    .groups = "drop"
  ) %>%
  filter(!(study == "Nair, 2021"))

studies <-
  df_mort %>%
  group_by(Study) %>%
  summarise(Criteria = max(Criteria)) %>%
  rename(study = Study) %>%
  filter(!(study %in% c(
    "Nair, 2021")))

network <- mtc.network(data.ab = df_mtc,
                       studies = studies)

sum(network$data.ab$sampleSize)

regressor <- list(coefficient = 'shared',
                  variable = 'Criteria',
                  control = 'SOT')

############## Bayesian network meta-analysis - Mortality, interface network #############

model_nocovar <- mtc_model_n01(
  network,
  type="consistency",
  hy.prior = hy.prior.mort,
  likelihood = "binom",
  link = "logit")

result_nocovar <- mtc.run(model_nocovar)

saveRDS(result_nocovar,
        file = "network2_mortality_nocovar_nma.rds")

result_nocovar <- readRDS("network2_mortality_nocovar_nma.rds")

model <- mtc_model_n01(
  network,
  type="regression",
  hy.prior = hy.prior.mort,
  regressor=regressor,
  likelihood = "binom",
  link = "logit")

result_criteria <- mtc.run(model)

saveRDS(result_criteria,
        file = "network2_mortality_criteria_nma.rds")

result_criteria <- readRDS("network2_mortality_criteria_nma.rds")

summary(result_criteria)

gemtc_forestplot2_interface(result_criteria,
                            result_nocovar,
                            "mortality",
                            ylim_lower = 0.05)
ggsave("Figures/Network2_Forest_Mortality_Criteria.svg",
       width = 7, height = 4, units = "in")

# league table (no-covariate model)
league_table_mean(result_nocovar) %>%
  write_csv("network2_gemtc_mortality_overall_league_table.csv")
