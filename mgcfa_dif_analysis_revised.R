# =============================================================================
# DIF TESTING PIPELINE FOR A PLANNED-MISSINGNESS (MULTI-STUDY) DESIGN
# Multigroup ordinal CFA (WLSMV) with anchor purification (Woods, 2009)
#
# -----------------------------------------------------------------------------


library(lavaan)
library(dplyr)
library(tidyr)
library(ggplot2)

LV <- "distress"

# =============================================================================
# HELPERS
# =============================================================================

# Prepare the data for one pair: reference = study1, focal = study2.
prepare_pair_data <- function(data, items, study1, study2, study_var = "Study") {
  pd <- data[data[[study_var]] %in% c(study1, study2), , drop = FALSE]
  pd$group <- factor(pd[[study_var]], levels = c(study1, study2))
  for (it in items) {
    x <- pd[[it]]
    pd[[it]] <- if (is.factor(x)) droplevels(as.ordered(x)) else ordered(x)
  }
  pd
}

# Number of thresholds per item = observed categories in the pair - 1
count_thresholds <- function(pd, items) {
  setNames(vapply(items, function(it) nlevels(droplevels(pd[[it]])) - 1L, integer(1)),
           items)
}

# Categories that are observed in one group but empty in the other
# (lavaan cannot estimate shared thresholds for those).
empty_category_items <- function(pd, items) {
  items[vapply(items, function(it) {
    tab <- table(pd[[it]], pd$group)
    any(tab == 0)
  }, logical(1))]
}

# -----------------------------------------------------------------------------
# Model builder with explicit identification
# -----------------------------------------------------------------------------
build_mgcfa_model <- function(items, constrained_items, n_thresholds,
                              free_focal = TRUE, lv = LV) {

  stopifnot(all(constrained_items %in% items))
  if (free_focal && length(constrained_items) == 0) {
    stop("At least one constrained item is required to estimate the focal ",
         "group's latent mean and variance.")
  }

  load_terms <- vapply(items, function(it) {
    if (it %in% constrained_items) {
      sprintf("c(lambda_%s, lambda_%s)*%s", it, it, it)
    } else {
      sprintf("c(NA, NA)*%s", it)
    }
  }, character(1))

  lines <- paste0(lv, " =~ ", paste(load_terms, collapse = " + "))

  # Latent distribution: reference fixed (0, 1); focal free
  if (free_focal) {
    lines <- c(lines,
               sprintf("%s ~~ c(1, NA)*%s", lv, lv),
               sprintf("%s ~ c(0, NA)*1", lv))
  } else {
    lines <- c(lines,
               sprintf("%s ~~ c(1, 1)*%s", lv, lv),
               sprintf("%s ~ c(0, 0)*1", lv))
  }

  for (it in items) {
    k <- n_thresholds[[it]]
    if (it %in% constrained_items) {
      thr <- paste(sprintf("c(tau%d_%s, tau%d_%s)*t%d", 1:k, it, 1:k, it, 1:k),
                   collapse = " + ")
      # METHODOLOGICAL JUDGEMENT: the focal scale factor is freed only for
      # constrained items with >= 2 thresholds. For a binary item, freeing it
      # alongside a free latent mean/variance is not separately identified
      # (Millsap & Yun-Tein, 2004), so it is held at 1. Distress items here
      # are polytomous, so this branch should rarely matter.
      sf <- if (free_focal && k >= 2) "c(1, NA)" else "c(1, 1)"
    } else {
      thr <- paste(sprintf("c(NA, NA)*t%d", 1:k), collapse = " + ")
      sf  <- "c(1, 1)"
    }
    lines <- c(lines,
               sprintf("%s | %s", it, thr),
               sprintf("%s ~*~ %s*%s", it, sf, it))
  }

  paste(lines, collapse = "\n")
}

# -----------------------------------------------------------------------------
# Fit wrapper: never throws; records convergence, admissibility and warnings
# -----------------------------------------------------------------------------
fit_pair_model <- function(model, pd, items, missing = "listwise") {
  warns <- character()
  err   <- NULL
  fit <- withCallingHandlers(
    tryCatch(
      cfa(model, data = pd, group = "group", ordered = items,
          estimator = "WLSMV", parameterization = "delta",
          std.lv = FALSE, auto.fix.first = FALSE, missing = missing),
      error = function(e) { err <<- conditionMessage(e); NULL }
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  converged <- !is.null(fit) && isTRUE(lavInspect(fit, "converged"))
  admissible <- FALSE
  if (converged) {
    pe <- parTable(fit)
    focal_var <- pe$est[pe$lhs == LV & pe$op == "~~" & pe$rhs == LV & pe$group == 2]
    admissible <- length(focal_var) == 1 && focal_var > 0
  }

  list(fit = fit, ok = converged && admissible, converged = converged,
       admissible = admissible, error = err,
       warnings = if (length(warns)) paste(unique(warns), collapse = " | ") else NA_character_)
}

# Focal-group latent mean and variance (reference fixed at 0 / 1)
extract_latent_distribution <- function(fit) {
  pe <- parameterEstimates(fit)
  m  <- pe[pe$lhs == LV & pe$op == "~1" & pe$group == 2, ]
  v  <- pe[pe$lhs == LV & pe$op == "~~" & pe$rhs == LV & pe$group == 2, ]
  list(mean_diff = m$est, mean_diff_se = m$se,
       focal_var = v$est, focal_var_se = v$se)
}

# Scaled difference test, consistently ordered (more constrained model first)
scaled_diff_test <- function(fit_constrained, fit_free) {
  lr <- lavTestLRT(fit_constrained, fit_free, method = "satorra.2000")
  list(chi_sq = lr$`Chisq diff`[2], df = lr$`Df diff`[2], p = lr$`Pr(>Chisq)`[2])
}

# =============================================================================
# STEP 1: ANCHOR COVERAGE — now records ALL pairs, including untestable ones
# =============================================================================

identify_anchor_coverage <- function(data, anchor_items, study_var = "Study",
                                     min_responses = 10, min_shared = 3) {

  cat("STEP 1: ANALYSING ANCHOR ITEM COVERAGE\n")
  cat("======================================\n\n")

  coverage_matrix <- data %>%
    select(all_of(c(study_var, anchor_items))) %>%
    group_by(!!sym(study_var)) %>%
    summarise(across(all_of(anchor_items), ~ sum(!is.na(.x)) >= min_responses),
              .groups = "drop")

  studies <- as.character(unique(data[[study_var]]))
  cov_mat <- as.matrix(coverage_matrix[, anchor_items])
  rownames(cov_mat) <- as.character(coverage_matrix[[study_var]])

  all_pairs <- list()
  for (i in 1:(length(studies) - 1)) {
    for (j in (i + 1):length(studies)) {
      s1 <- studies[i]; s2 <- studies[j]
      shared <- anchor_items[cov_mat[s1, ] & cov_mat[s2, ]]
      all_pairs[[paste(s1, s2, sep = "_vs_")]] <- list(
        study1 = s1, study2 = s2, shared_anchors = shared,
        n_anchors = length(shared), testable = length(shared) >= min_shared
      )
    }
  }

  pair_table <- data.frame(
    comparison = names(all_pairs),
    study1     = vapply(all_pairs, `[[`, "", "study1"),
    study2     = vapply(all_pairs, `[[`, "", "study2"),
    n_shared   = vapply(all_pairs, `[[`, 0L, "n_anchors"),
    testable   = vapply(all_pairs, `[[`, TRUE, "testable"),
    shared_items = vapply(all_pairs, function(p) paste(p$shared_anchors, collapse = ";"), ""),
    stringsAsFactors = FALSE, row.names = NULL
  )

  cat("Total anchor items:", length(anchor_items), "\n")
  cat("Total study pairs:", nrow(pair_table), "\n")
  cat("Testable pairs (>=", min_shared, "shared items):", sum(pair_table$testable), "\n")
  cat("Untestable by design (<", min_shared, "shared items):", sum(!pair_table$testable), "\n\n")

  list(coverage_matrix = coverage_matrix,
       pairwise_tests  = all_pairs[pair_table$testable],
       all_pairs       = all_pairs,
       pair_table      = pair_table,
       all_studies     = studies)
}

# =============================================================================
# STEP 2: PURIFICATION — each item freed in turn against all others as anchors
# =============================================================================

purify_anchors <- function(data, anchor_items, study1, study2,
                           study_var = "Study", min_anchors = 3,
                           cfi_threshold = 0.01, tli_threshold = 0.01,
                           purification_criterion = c("fit_indices", "lrt", "rank"),
                           alpha = 0.05, missing = "listwise") {

  # purification_criterion — how items are cleared as anchors:
  #   "fit_indices": not flagged if |dCFI| <= cfi_threshold and |dTLI| <=
  #                  tli_threshold (the ORIGINAL rule; kept for a like-for-like
  #                  rerun). Simulation shows it is insensitive to single-item
  #                  DIF once the latent mean is correctly identified.
  #   "lrt":         not flagged if BH-adjusted scaled chi-square p >= alpha.
  #                  With large samples, a DIF item left among the anchors
  #                  biases the others, so clean items get flagged too.
  #   "rank":        Woods (2009): the min_anchors items with the smallest
  #                  chi-square/df are the anchors; all others are studied.
  #                  Never triggers the fallback.
  # METHODOLOGICAL JUDGEMENT: which rule to report is a substantive choice.
  purification_criterion <- match.arg(purification_criterion)
  cat(sprintf("\nSTEP 2: ANCHOR PURIFICATION FOR %s vs %s\n", study1, study2))

  out <- list(status = NA_character_, purified_anchors = character(),
              forced_anchors = character(), fallback_used = FALSE,
              valid_anchors = character(), excluded_items = character(),
              fit_differences = NULL, invariant_model = NULL)

  pd <- prepare_pair_data(data, anchor_items, study1, study2, study_var)

  # Items with insufficient data/variation in the pair
  ok_item <- vapply(anchor_items, function(it) {
    x <- pd[[it]]
    sum(!is.na(x)) >= 20 && nlevels(droplevels(x)) >= 2
  }, logical(1))
  valid <- anchor_items[ok_item]

  # Items with a category observed in one study only cannot share thresholds
  empty_cat <- empty_category_items(pd, valid)
  if (length(empty_cat) > 0) {
    cat("  Excluding items with an empty category in one study:",
        paste(empty_cat, collapse = ", "), "\n")
    # METHODOLOGICAL JUDGEMENT: exclusion (rather than collapsing categories)
    # keeps the item set honest but reduces the shared item count. Collapsing
    # sparse top categories within the pair is an alternative.
  }
  valid <- setdiff(valid, empty_cat)
  out$valid_anchors  <- valid
  out$excluded_items <- setdiff(anchor_items, valid)

  if (length(valid) < min_anchors) {
    out$status <- "insufficient_valid_items"
    cat("  Insufficient valid items for testing\n")
    return(out)
  }

  n_thr <- count_thresholds(pd, valid)

  # Reference model: all items constrained, focal mean/variance free.
  # This is also the "full loading+threshold invariance" model, from which
  # the invariance-based latent mean difference is reported.
  ref <- fit_pair_model(build_mgcfa_model(valid, valid, n_thr), pd, valid, missing)
  if (!ref$ok) {
    out$status <- "nonconverged_purification_reference"
    cat("  Purification reference model failed:",
        if (!is.null(ref$error)) ref$error else "not converged / inadmissible", "\n")
    return(out)
  }
  out$invariant_model <- ref
  ref_idx <- unclass(fitMeasures(ref$fit, c("cfi", "tli")))

  res <- data.frame(item = valid, chi_sq = NA_real_, df = NA_real_,
                    p_value = NA_real_, delta_cfi = NA_real_, delta_tli = NA_real_,
                    mean_diff_item_freed = NA_real_, converged = FALSE,
                    warnings = NA_character_, stringsAsFactors = FALSE)

  for (i in seq_along(valid)) {
    it  <- valid[i]
    cmp <- fit_pair_model(build_mgcfa_model(valid, valid[-i], n_thr), pd, valid, missing)
    res$warnings[i] <- cmp$warnings
    if (!cmp$ok) next
    tst <- scaled_diff_test(ref$fit, cmp$fit)
    ci  <- unclass(fitMeasures(cmp$fit, c("cfi", "tli")))
    res$chi_sq[i]    <- tst$chi_sq
    res$df[i]        <- tst$df
    res$p_value[i]   <- tst$p
    # Positive = freeing the item improved fit
    res$delta_cfi[i] <- ci["cfi"] - ref_idx["cfi"]
    res$delta_tli[i] <- ci["tli"] - ref_idx["tli"]
    res$mean_diff_item_freed[i] <- extract_latent_distribution(cmp$fit)$mean_diff
    res$converged[i] <- TRUE
  }

  conv <- res[res$converged, ]
  if (nrow(conv) == 0) {
    out$status <- "nonconverged_purification_items"
    out$fit_differences <- res
    return(out)
  }

  conv <- conv %>%
    mutate(p_bh          = p.adjust(p_value, method = "BH"),
           composite_dif = (abs(delta_cfi) + abs(delta_tli)) / 2,
           chisq_per_df  = chi_sq / df,
           flagged = switch(purification_criterion,
             fit_indices = abs(delta_cfi) > cfi_threshold | abs(delta_tli) > tli_threshold,
             lrt         = p_bh < alpha,
             rank        = rank(chisq_per_df, ties.method = "first") > min_anchors)) %>%
    arrange(if (purification_criterion == "rank") chisq_per_df else composite_dif)

  not_flagged <- conv$item[!conv$flagged]

  if (length(not_flagged) >= min_anchors) {
    out$purified_anchors <- not_flagged
  } else {
    # Forced-anchor fallback (retained from the original procedure, now
    # recorded). The items with the smallest composite fit change are used as
    # anchors EVEN THOUGH some were flagged. Results for such pairs rest on an
    # anchor set that the procedure itself did not clear.
    n_select <- min(min_anchors, nrow(conv))
    chosen   <- head(conv$item, n_select)
    out$purified_anchors <- chosen
    out$forced_anchors   <- intersect(chosen, conv$item[conv$flagged])
    out$fallback_used    <- TRUE
    cat("  Forced-anchor fallback used; flagged items kept as anchors:",
        paste(out$forced_anchors, collapse = ", "), "\n")
  }

  # Items that did not converge at purification cannot be cleared as anchors
  out$fit_differences <- res %>%
    left_join(conv %>% select(item, p_bh, composite_dif, flagged), by = "item")
  out$n_thresholds <- n_thr
  out$status <- "purified"

  cat(sprintf("  Flagged at purification: %d of %d converged items; anchors: %s\n",
              sum(conv$flagged), nrow(conv), paste(out$purified_anchors, collapse = ", ")))
  out
}

# =============================================================================
# STEP 3: FORMAL DIF TESTING WITH PURIFIED ANCHORS
# =============================================================================

formal_dif_testing <- function(data, purified_result, study1, study2,
                               study_var = "Study", alpha = 0.05,
                               cfi_threshold = 0.01, tli_threshold = 0.01,
                               missing = "listwise") {

  cat(sprintf("STEP 3: FORMAL DIF TESTING FOR %s vs %s\n", study1, study2))

  anchors <- purified_result$purified_anchors
  # Items that failed to converge at purification are neither anchors nor
  # studied items; they are reported as not converged.
  conv_items <- purified_result$fit_differences$item[purified_result$fit_differences$converged]
  studied <- setdiff(conv_items, anchors)
  items   <- c(anchors, studied)
  n_thr   <- purified_result$n_thresholds[items]

  out <- list(status = NA_character_, dif_results = NULL, purified_anchors = anchors,
              studied_items = studied, dif_items = character(), final_model = NULL)

  pd <- prepare_pair_data(data, items, study1, study2, study_var)

  if (length(studied) == 0) {
    # Nothing flagged at purification: the pair has NO DIF (not a failure).
    out$status <- "no_dif"
    out$final_model <- fit_pair_model(build_mgcfa_model(items, items, n_thr),
                                      pd, items, missing)
    cat("  No items flagged at purification -> no DIF\n")
    return(out)
  }

  ref <- fit_pair_model(build_mgcfa_model(items, anchors, n_thr), pd, items, missing)
  if (!ref$ok) {
    out$status <- "nonconverged_formal_reference"
    return(out)
  }
  ref_idx <- unclass(fitMeasures(ref$fit, c("cfi", "tli")))

  res <- data.frame(item = studied, chi_sq = NA_real_, df = NA_real_,
                    p_value = NA_real_, delta_cfi = NA_real_, delta_tli = NA_real_,
                    converged = FALSE, warnings = NA_character_,
                    stringsAsFactors = FALSE)

  for (i in seq_along(studied)) {
    cmp <- fit_pair_model(build_mgcfa_model(items, c(anchors, studied[i]), n_thr),
                          pd, items, missing)
    res$warnings[i] <- cmp$warnings
    if (!cmp$ok) next
    tst <- scaled_diff_test(cmp$fit, ref$fit)
    ci  <- unclass(fitMeasures(cmp$fit, c("cfi", "tli")))
    res$chi_sq[i]    <- tst$chi_sq
    res$df[i]        <- tst$df
    res$p_value[i]   <- tst$p
    # Positive = constraining the item worsened fit
    res$delta_cfi[i] <- ref_idx["cfi"] - ci["cfi"]
    res$delta_tli[i] <- ref_idx["tli"] - ci["tli"]
    res$converged[i] <- TRUE
  }

  res <- res %>%
    mutate(p_bh = ifelse(converged, p.adjust(ifelse(converged, p_value, NA), method = "BH"), NA),
           significant_BH = converged & !is.na(p_bh) & p_bh < alpha,
           dif_by_fit_indices = converged & (abs(delta_cfi) > cfi_threshold |
                                               abs(delta_tli) > tli_threshold))

  out$dif_results <- res
  out$dif_items   <- res$item[res$significant_BH]

  # Final partial-invariance model: anchors + non-DIF studied items constrained,
  # DIF items free. Its focal mean is the reported (DIF-adjusted) mean difference.
  # METHODOLOGICAL JUDGEMENT: studied items that did not reach BH significance
  # are constrained here (treated as invariant); items that failed to converge
  # in the formal stage are left free so that no untested constraint is imposed.
  free_final <- union(out$dif_items, res$item[!res$converged])
  out$final_model <- fit_pair_model(build_mgcfa_model(items, setdiff(items, free_final), n_thr),
                                    pd, items, missing)
  out$status <- if (all(res$converged)) "dif_tested" else "dif_tested_partial_convergence"

  cat(sprintf("  Studied: %d; significant DIF (BH): %d\n", length(studied), length(out$dif_items)))
  out
}

# =============================================================================
# STEP 4: IMPACT — compare scores ignoring vs accounting for DIF
# =============================================================================

assess_dif_impact <- function(data, purified_result, formal_results, study1, study2,
                              study_var = "Study", missing = "listwise") {

  dif_items <- formal_results$dif_items
  if (length(dif_items) == 0) return(NULL)

  model1 <- purified_result$invariant_model           # all constrained
  model2 <- formal_results$final_model                # DIF items free
  if (is.null(model1) || is.null(model2) || !model1$ok || !model2$ok) return(NULL)

  # Model 1 was fitted on the purification item set; if any item dropped out
  # through non-convergence the item sets differ, so refit Model 1 on the
  # formal item set to compare like with like.
  items2 <- lavNames(model2$fit, "ov")
  items1 <- lavNames(model1$fit, "ov")
  if (!setequal(items1, items2)) {
    pd <- prepare_pair_data(data, items2, study1, study2, study_var)
    model1 <- fit_pair_model(build_mgcfa_model(items2, items2,
                                               purified_result$n_thresholds[items2]),
                             pd, items2, missing)
    if (!model1$ok) return(NULL)
  }

  get_scores <- function(fit, name) {
    sc  <- lavPredict(fit)
    idx <- lavInspect(fit, "case.idx")
    grp <- lavInspect(fit, "group.label")
    bind_rows(lapply(seq_along(grp), function(g)
      setNames(data.frame(idx[[g]], as.numeric(sc[[g]]), grp[g]),
               c("case_idx", name, "group"))))
  }

  s <- inner_join(get_scores(model1$fit, "score1"),
                  get_scores(model2$fit, "score2"), by = c("case_idx", "group"))
  s$group <- factor(s$group, levels = c(study1, study2))

  r    <- cor(s$score1, s$score2)
  bias <- s$score2 - s$score1
  mad  <- mean(abs(s$score1 - s$score2))
  rmse <- sqrt(mean((s$score1 - s$score2)^2))
  # Lin's concordance correlation (replaces the original "ICC", which omitted
  # the squared mean difference and so could not detect a shift).
  ccc  <- 2 * cov(s$score1, s$score2) /
    (var(s$score1) + var(s$score2) + (mean(s$score1) - mean(s$score2))^2)

  ld1 <- extract_latent_distribution(model1$fit)
  ld2 <- extract_latent_distribution(model2$fit)

  # Per-group mean shift of individual scores between models
  shift_by_group <- s %>% group_by(group) %>%
    summarise(mean_shift = mean(score2 - score1), .groups = "drop")

  impact <- if (r > 0.95 && mad < 0.2) "Minimal" else
            if (r > 0.90 && mad < 0.5) "Moderate" else "Substantial"

  p1 <- ggplot(s, aes(score1, score2, colour = group)) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = sprintf("Factor scores: %s vs %s", study1, study2),
         subtitle = sprintf("r = %.3f, MAD = %.3f, CCC = %.3f; latent mean diff: invariant %.3f, partial %.3f",
                            r, mad, ccc, ld1$mean_diff, ld2$mean_diff),
         x = "Model 1: all items constrained", y = "Model 2: DIF items free",
         colour = NULL) +
    theme_minimal() + theme(legend.position = "bottom")

  # Difference plot: makes a uniform per-group offset (bands) vs a
  # trait-dependent pattern directly visible.
  p2 <- ggplot(s, aes((score1 + score2) / 2, score2 - score1, colour = group)) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_smooth(method = "loess", se = FALSE, formula = y ~ x) +
    labs(title = sprintf("Score difference by trait level: %s vs %s", study1, study2),
         x = "Mean of the two scores", y = "Model 2 - Model 1", colour = NULL) +
    theme_minimal() + theme(legend.position = "bottom")

  list(correlation = r, bias = bias, mad = mad, rmse = rmse, ccc = ccc, impact = impact,
       mean_diff_invariant = ld1$mean_diff, mean_diff_partial = ld2$mean_diff,
       shift_by_group = shift_by_group, n_compared = nrow(s),
       comparison_plot = p1, difference_plot = p2, score_data = s)
}

# =============================================================================
# ITEM-LEVEL TABLE — every shared item in a pair, with a final status
# =============================================================================

build_item_table <- function(pair_name, pair_info, purified, formal) {
  shared <- pair_info$shared_anchors
  it <- data.frame(comparison = pair_name, study1 = pair_info$study1,
                   study2 = pair_info$study2, item = shared,
                   stringsAsFactors = FALSE)

  pf <- purified$fit_differences
  if (!is.null(pf)) {
    it <- it %>% left_join(
      pf %>% transmute(item, pur_chi_sq = chi_sq, pur_df = df, pur_p = p_value,
                       pur_p_bh = if ("p_bh" %in% names(pf)) p_bh else NA_real_,
                       pur_delta_cfi = delta_cfi, pur_delta_tli = delta_tli,
                       pur_converged = converged,
                       pur_flagged = if ("flagged" %in% names(pf)) flagged else NA),
      by = "item")
  }
  if (!is.null(formal) && !is.null(formal$dif_results)) {
    it <- it %>% left_join(
      formal$dif_results %>% transmute(item, formal_chi_sq = chi_sq, formal_df = df,
                                       formal_p = p_value, formal_p_bh = p_bh,
                                       formal_delta_cfi = delta_cfi,
                                       formal_delta_tli = delta_tli,
                                       formal_converged = converged,
                                       dif_significant_BH = significant_BH),
      by = "item")
  }
  for (col in c("pur_converged", "pur_flagged", "formal_converged", "dif_significant_BH"))
    if (!col %in% names(it)) it[[col]] <- NA

  it$role <- case_when(
    it$item %in% purified$excluded_items          ~ "excluded",
    it$item %in% purified$forced_anchors          ~ "anchor_forced",
    it$item %in% purified$purified_anchors        ~ "anchor",
    !is.na(it$pur_converged) & !it$pur_converged  ~ "nonconverged_purification",
    TRUE                                          ~ "studied"
  )

  it$final_status <- case_when(
    it$role == "excluded"                                     ~ "excluded (insufficient data/empty category)",
    it$role == "anchor"                                       ~ "no DIF (anchor, not flagged)",
    it$role == "anchor_forced"                                ~ "assumed invariant (forced anchor, flagged)",
    it$role == "nonconverged_purification"                    ~ "not converged",
    it$role == "studied" & isTRUE(formal$status == "no_dif")  ~ "no DIF (anchor, not flagged)",
    it$role == "studied" & !is.na(it$formal_converged) & !it$formal_converged ~ "not converged",
    it$role == "studied" & it$dif_significant_BH %in% TRUE    ~ "DIF",
    it$role == "studied" & it$dif_significant_BH %in% FALSE   ~ "no DIF (flagged, not significant)",
    TRUE                                                      ~ "not tested (pair failed)"
  )
  it
}

# =============================================================================
# MASTER PIPELINE
# =============================================================================

run_comprehensive_dif_pipeline <- function(data, anchor_items,
                                           study_var = "Study",
                                           min_responses = 10,
                                           min_anchors = 3,
                                           alpha = 0.05,
                                           cfi_threshold = 0.01,
                                           tli_threshold = 0.01,
                                           purification_criterion = "fit_indices",
                                           missing = "listwise",
                                           output_dir = "dif_results",
                                           save_plots = TRUE) {

  cat("DIF PIPELINE (revised: focal latent mean/variance freely estimated)\n\n")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  coverage <- identify_anchor_coverage(data, anchor_items, study_var,
                                       min_responses, min_shared = min_anchors)

  pair_rows <- list(); item_rows <- list(); all_results <- list()

  for (pair_name in coverage$pair_table$comparison) {
    pinfo <- coverage$all_pairs[[pair_name]]
    row <- list(comparison = pair_name, study1 = pinfo$study1, study2 = pinfo$study2,
                reference_group = pinfo$study1, focal_group = pinfo$study2,
                n_shared = pinfo$n_anchors)

    if (!pinfo$testable) {
      row$status <- "untestable_fewer_than_min_shared"
      pair_rows[[pair_name]] <- row
      next
    }

    cat("\n----", pair_name, "----\n")
    purified <- purify_anchors(data, pinfo$shared_anchors, pinfo$study1, pinfo$study2,
                               study_var, min_anchors, cfi_threshold, tli_threshold,
                               purification_criterion, alpha, missing)
    formal <- NULL; impact <- NULL

    row$n_valid        <- length(purified$valid_anchors)
    row$fallback_used  <- purified$fallback_used
    row$forced_anchors <- paste(purified$forced_anchors, collapse = ";")
    row$anchors        <- paste(purified$purified_anchors, collapse = ";")

    if (!is.null(purified$invariant_model) && purified$invariant_model$ok) {
      f  <- purified$invariant_model$fit
      ld <- extract_latent_distribution(f)
      nn <- lavInspect(f, "nobs")
      row$n_used_ref <- nn[1]; row$n_used_focal <- nn[2]
      row$mean_diff_invariant    <- ld$mean_diff
      row$mean_diff_invariant_se <- ld$mean_diff_se
      row$focal_var_invariant    <- ld$focal_var
      row$cfi_invariant   <- as.numeric(fitMeasures(f, "cfi"))
      row$rmsea_invariant <- as.numeric(fitMeasures(f, "rmsea"))
    }

    if (purified$status != "purified") {
      row$status <- purified$status
    } else {
      formal <- formal_dif_testing(data, purified, pinfo$study1, pinfo$study2,
                                   study_var, alpha, cfi_threshold, tli_threshold, missing)
      row$status <- formal$status
      row$n_flagged_purification <- sum(purified$fit_differences$flagged %in% TRUE)
      row$n_tested_formal <- length(formal$studied_items)
      row$n_dif <- length(formal$dif_items)
      row$dif_items <- paste(formal$dif_items, collapse = ";")

      fm <- formal$final_model
      if (!is.null(fm) && fm$ok) {
        ld <- extract_latent_distribution(fm$fit)
        row$mean_diff_partial    <- ld$mean_diff
        row$mean_diff_partial_se <- ld$mean_diff_se
        row$focal_var_partial    <- ld$focal_var
      } else if (!is.null(fm)) {
        row$status <- paste0(row$status, "; final_model_nonconverged")
      }

      if (length(formal$dif_items) > 0) {
        impact <- assess_dif_impact(data, purified, formal, pinfo$study1, pinfo$study2,
                                    study_var, missing)
        if (!is.null(impact)) {
          row[c("correlation", "bias", "mad", "rmse", "ccc", "s", "impact")] <-
            impact[c("correlation", "bias", "mad", "rmse", "ccc", "s", "impact")]
          if (save_plots) {
            ggsave(file.path(output_dir, paste0(pair_name, "_score_comparison.png")),
                   impact$comparison_plot, width = 8, height = 6)
            ggsave(file.path(output_dir, paste0(pair_name, "_score_difference.png")),
                   impact$difference_plot, width = 8, height = 6)
          }
        }
      }
    }

    items_tbl <- build_item_table(pair_name, pinfo, purified, formal)
    item_rows[[pair_name]] <- items_tbl
    row$n_dif_items_all_shared <- sum(items_tbl$final_status == "DIF")
    row$prop_dif_of_shared     <- row$n_dif_items_all_shared / pinfo$n_anchors

    pair_rows[[pair_name]] <- row
    all_results[[pair_name]] <- list(purified = purified, formal = formal, impact = impact)
  }

  pair_summary <- bind_rows(lapply(pair_rows, as.data.frame, stringsAsFactors = FALSE))
  item_level   <- bind_rows(item_rows)

  if (!"impact" %in% names(pair_summary)) pair_summary$impact <- NA_character_
  no_dif_pair <- pair_summary$status %in% "no_dif" |
    (grepl("^dif_tested", pair_summary$status) & pair_summary$n_dif %in% 0)
  pair_summary$impact[is.na(pair_summary$impact) & no_dif_pair] <- "Not applicable (no DIF)"

  write.csv(pair_summary, file.path(output_dir, "dif_pair_summary.csv"), row.names = FALSE)
  write.csv(item_level,   file.path(output_dir, "dif_item_level.csv"),   row.names = FALSE)
  write.csv(pair_summary[pair_summary$status == "untestable_fewer_than_min_shared",
                         c("comparison", "study1", "study2", "n_shared")],
            file.path(output_dir, "untestable_pairs.csv"), row.names = FALSE)

  cat("\n=========== SUMMARY ===========\n")
  print(table(pair_summary$status))
  cat("\nPairs using forced-anchor fallback:", sum(pair_summary$fallback_used %in% TRUE), "\n")
  if (nrow(item_level)) {
    cat("\nItem-level final status (all shared items in testable pairs):\n")
    print(table(item_level$final_status))
  }

  invisible(list(coverage = coverage, pair_summary = pair_summary,
                 item_level = item_level, pairwise_results = all_results))
}

# =============================================================================
# EXAMPLE USAGE (wrapped so that sourcing this file only defines functions)
# =============================================================================
if (FALSE) {
  data_dif     <- results[["data_preparation"]][["reference_wave_data"]]
  anchor_items <- results$anchor_analysis$all_anchors

  dif_results <- run_comprehensive_dif_pipeline(
    data = data_dif, anchor_items = anchor_items, study_var = "Study",
    min_responses = 10, min_anchors = 3, alpha = 0.05,
    cfi_threshold = 0.001, tli_threshold = 0.001,
    purification_criterion = "fit_indices",   # as in the original analysis
    missing = "listwise",                     # see note on planned missingness
    output_dir = "dif_results_revised2"
  )

  # Before running all pairs, confirm the identification on one pair:
  # pt <- parTable(dif_results$pairwise_results[[1]]$purified$invariant_model$fit)
  # pt[pt$lhs == "distress" & pt$op %in% c("~1", "~~"), ]
  # -> group 1: fixed 0 / 1; group 2: free (free > 0)
}
