# =============================================================================
# COMPREHENSIVE DIF TESTING PIPELINE FOR PLANNED MISSINGNESS DESIGN
# Using Multigroup CFA Approach with Purification (Woods, 2009)
# =============================================================================

library(lavaan)
library(dplyr)
library(tidyr)
library(ggplot2)

# =============================================================================
# STEP 1: IDENTIFY STUDY-SPECIFIC ANCHOR AVAILABILITY
# =============================================================================

identify_anchor_coverage <- function(data, anchor_items, study_var = "Study", 
                                     min_responses = 10) {
  
  cat("STEP 1: ANALYZING ANCHOR ITEM COVERAGE\n")
  cat("======================================\n\n")
  
  # Calculate coverage matrix
  coverage_matrix <- data %>%
    select(all_of(c(study_var, anchor_items))) %>%
    group_by(!!sym(study_var)) %>%
    summarise(across(all_of(anchor_items), 
                     ~ sum(!is.na(.x)) >= min_responses,
                     .names = "{.col}"),
              .groups = "drop")
  
  # Identify pairwise testable study combinations
  studies <- unique(data[[study_var]])
  pairwise_tests <- list()
  
  for (i in 1:(length(studies)-1)) {
    for (j in (i+1):length(studies)) {
      study1 <- studies[i]
      study2 <- studies[j]
      
      # Find shared anchors with sufficient data
      shared_anchors <- anchor_items[
        coverage_matrix[coverage_matrix[[study_var]] == study1, anchor_items] &
          coverage_matrix[coverage_matrix[[study_var]] == study2, anchor_items]
      ]
      
      if (length(shared_anchors) >= 3) {  # Minimum 3 anchors
        pairwise_tests[[paste(study1, study2, sep = "_vs_")]] <- list(
          study1 = study1,
          study2 = study2,
          shared_anchors = shared_anchors,
          n_anchors = length(shared_anchors)
        )
      }
    }
  }
  
  cat("Total anchor items:", length(anchor_items), "\n")
  cat("Testable study pairs:", length(pairwise_tests), "\n\n")
  
  # Print coverage summary
  for (pair_name in names(pairwise_tests)) {
    pair <- pairwise_tests[[pair_name]]
    cat(sprintf("%s: %d shared anchors\n", 
                pair_name, pair$n_anchors))
  }
  cat("\n")
  
  return(list(
    coverage_matrix = coverage_matrix,
    pairwise_tests = pairwise_tests,
    all_studies = studies
  ))
}

# =============================================================================
# STEP 2: PURIFICATION STAGE - IDENTIFY DIF-FREE ANCHORS
# =============================================================================

purify_anchors <- function(data, anchor_items, study1, study2,
                           study_var = "Study", min_anchors = 3,
                           cfi_threshold = 0.01, tli_threshold = 0.01) {

  cat(sprintf("\nSTEP 2: ANCHOR PURIFICATION FOR %s vs %s\n", study1, study2))
  cat("=================================================\n\n")

  # Filter to relevant studies
  pair_data <- data %>%
    filter(!!sym(study_var) %in% c(study1, study2)) %>%
    mutate(group = as.factor(!!sym(study_var)))

  # Remove items with insufficient variation
  valid_anchors <- anchor_items[sapply(anchor_items, function(item) {
    if (!item %in% names(pair_data)) return(FALSE)
    n_valid <- sum(!is.na(pair_data[[item]]))
    n_unique <- length(unique(pair_data[[item]][!is.na(pair_data[[item]])]))
    return(n_valid >= 20 && n_unique >= 2)
  })]

  if (length(valid_anchors) < min_anchors) {
    cat("Insufficient valid anchors for testing\n")
    return(NULL)
  }

  # --- Estimate shared reference model ONCE (all valid anchors constrained) ---
  # Each item is then tested by comparing a model with that item freed against
  # this single reference, so all LRTs share the same baseline.
  # Sign convention: ΔCFI = comparison CFI - reference CFI.
  # Positive ΔCFI means freeing the item improved fit → DIF present.
  ref_model_purification <- build_mgcfa_model(
    items        = valid_anchors,
    anchor_items = valid_anchors,
    ordinal      = TRUE,
    data         = pair_data
  )

  cat("Fitting purification reference model (all anchors constrained)...\n")
  ref_fit_purification <- tryCatch(
    cfa(ref_model_purification, data = pair_data, group = "group",
        ordered = TRUE, estimator = "WLSMV", std.lv = TRUE),
    error = function(e) {
      cat(sprintf("Reference model error: %s\n", e$message))
      NULL
    }
  )

  if (is.null(ref_fit_purification) || !lavInspect(ref_fit_purification, "converged")) {
    cat("Purification reference model did not converge. Cannot proceed.\n")
    return(NULL)
  }

  ref_indices_pur <- fitMeasures(ref_fit_purification, c("cfi", "tli"))
  cat(sprintf("Reference model: CFI = %.4f, TLI = %.4f\n\n",
              ref_indices_pur["cfi"], ref_indices_pur["tli"]))

  # Results table — one row per anchor item
  fit_differences <- data.frame(
    item           = valid_anchors,
    chi_sq         = NA_real_,
    df             = NA_real_,
    p_value        = NA_real_,
    delta_cfi      = NA_real_,
    delta_tli      = NA_real_,
    cfi_reference  = NA_real_,
    cfi_comparison = NA_real_,
    tli_reference  = NA_real_,
    tli_comparison = NA_real_,
    converged      = FALSE,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(valid_anchors)) {
    tested_item  <- valid_anchors[i]
    temp_anchors <- valid_anchors[-i]

    cat(sprintf("Testing item %d/%d: %s\n", i, length(valid_anchors), tested_item))

    # Comparison model: full item set, tested item freed across groups
    comparison_model <- build_mgcfa_model(
      items        = valid_anchors,
      anchor_items = temp_anchors,
      ordinal      = TRUE,
      data         = pair_data
    )

    tryCatch({
      fit_comparison <- cfa(comparison_model, data = pair_data,
                            group = "group", ordered = TRUE,
                            estimator = "WLSMV", std.lv = TRUE)

      if (lavInspect(fit_comparison, "converged")) {

        # LRT: reference (more constrained) vs comparison (item freed)
        lr_test <- lavTestLRT(ref_fit_purification, fit_comparison,
                              method = "satorra.2000")

        fit_differences$chi_sq[i]  <- lr_test$`Chisq diff`[2]
        fit_differences$df[i]      <- lr_test$`Df diff`[2]
        fit_differences$p_value[i] <- lr_test$`Pr(>Chisq)`[2]

        comp_indices <- fitMeasures(fit_comparison, c("cfi", "tli"))

        fit_differences$cfi_reference[i]  <- ref_indices_pur["cfi"]
        fit_differences$cfi_comparison[i] <- comp_indices["cfi"]
        fit_differences$tli_reference[i]  <- ref_indices_pur["tli"]
        fit_differences$tli_comparison[i] <- comp_indices["tli"]

        # ΔCFI = comparison - reference (positive = freeing item improved fit = DIF)
        fit_differences$delta_cfi[i] <- comp_indices["cfi"] - ref_indices_pur["cfi"]
        fit_differences$delta_tli[i] <- comp_indices["tli"] - ref_indices_pur["tli"]
        fit_differences$converged[i] <- TRUE

        cat(sprintf("  χ²(%.0f) = %.2f, p = %.4f, ΔCFI = %.4f, ΔTLI = %.4f\n",
                    fit_differences$df[i], fit_differences$chi_sq[i],
                    fit_differences$p_value[i],
                    fit_differences$delta_cfi[i],
                    fit_differences$delta_tli[i]))
      } else {
        cat("  Comparison model did not converge\n")
      }
    }, error = function(e) {
      cat(sprintf("  Error: %s\n", e$message))
    })
  }

  # Filter to converged items
  converged_results <- fit_differences %>%
    filter(converged)

  if (nrow(converged_results) == 0) {
    cat("\nNo items converged successfully\n")
    return(NULL)
  }

  # Calculate composite DIF indicator (average of absolute differences)
  converged_results <- converged_results %>%
    mutate(
      abs_delta_cfi = abs(delta_cfi),
      abs_delta_tli = abs(delta_tli),
      composite_dif = (abs_delta_cfi + abs_delta_tli) / 2
    ) %>%
    arrange(composite_dif)

  # Identify items with potential DIF
  # Items show DIF if ΔCFI or ΔTLI exceed thresholds
  converged_results <- converged_results %>%
    mutate(
      shows_dif = abs_delta_cfi > cfi_threshold | abs_delta_tli > tli_threshold
    )

  n_no_dif   <- sum(!converged_results$shows_dif)
  n_with_dif <- sum(converged_results$shows_dif)

  cat("\nPURIFICATION RESULTS:\n")
  cat(sprintf("Items tested: %d\n", nrow(converged_results)))
  cat(sprintf("Items with no DIF (ΔCFI ≤ %.3f, ΔTLI ≤ %.3f): %d\n",
              cfi_threshold, tli_threshold, n_no_dif))
  cat(sprintf("Items with potential DIF: %d\n\n", n_with_dif))

  # Select purified anchors
  if (n_no_dif >= min_anchors) {
    purified_anchors <- converged_results %>%
      filter(!shows_dif) %>%
      pull(item)
    cat(sprintf("Using %d items with no DIF as purified anchors:\n", length(purified_anchors)))
  } else {
    n_select <- min(min_anchors, nrow(converged_results))
    purified_anchors <- converged_results %>%
      slice_head(n = n_select) %>%
      pull(item)
    cat(sprintf("All items show potential DIF. Selecting %d items with smallest differences:\n",
                n_select))
  }

  purified_info <- converged_results %>%
    filter(item %in% purified_anchors) %>%
    select(item, delta_cfi, delta_tli, composite_dif)

  for (i in 1:nrow(purified_info)) {
    cat(sprintf("  %s: ΔCFI = %.4f, ΔTLI = %.4f, Composite = %.4f\n",
                purified_info$item[i],
                purified_info$delta_cfi[i],
                purified_info$delta_tli[i],
                purified_info$composite_dif[i]))
  }
  cat("\n")

  return(list(
    purified_anchors = purified_anchors,
    fit_differences  = converged_results,
    valid_anchors    = valid_anchors,
    n_no_dif         = n_no_dif,
    n_with_dif       = n_with_dif
  ))
}

# =============================================================================
# STEP 3: FORMAL DIF TESTING WITH PURIFIED ANCHORS
# =============================================================================

formal_dif_testing <- function(data, purified_result, study1, study2,
                               study_var = "Study", alpha = 0.05,
                               cfi_threshold = 0.01, tli_threshold = 0.01) {

  cat(sprintf("\nSTEP 3: FORMAL DIF TESTING FOR %s vs %s\n", study1, study2))
  cat("===========================================\n\n")

  if (is.null(purified_result)) {
    cat("No purified anchors available\n")
    return(NULL)
  }

  purified_anchors <- purified_result$purified_anchors
  all_items        <- purified_result$valid_anchors
  studied_items    <- setdiff(all_items, purified_anchors)

  if (length(studied_items) == 0) {
    cat("No items left to test after purification\n")
    return(NULL)
  }

  # Filter to relevant studies
  pair_data <- data %>%
    filter(!!sym(study_var) %in% c(study1, study2)) %>%
    mutate(group = as.factor(!!sym(study_var)))

  full_item_set <- c(studied_items, purified_anchors)

  # --- Estimate shared reference model ONCE ---
  # Reference: anchors constrained, ALL studied items free.
  # Each item is then tested by additionally constraining it (comparison model)
  # and comparing that more-constrained model back to this reference.
  # Sign convention differs from purification: here comparison is MORE constrained,
  # so ΔCFI = reference CFI - comparison CFI. Positive = constraining item worsened
  # fit = DIF present. Both stages yield positive ΔCFI when DIF is detected.
  ref_model_formal <- build_mgcfa_model(
    items        = full_item_set,
    anchor_items = purified_anchors,
    ordinal      = TRUE,
    data         = pair_data
  )

  cat("Fitting formal testing reference model (anchors constrained, studied items free)...\n")
  ref_fit_formal <- tryCatch(
    cfa(ref_model_formal, data = pair_data, group = "group",
        ordered = TRUE, estimator = "WLSMV", std.lv = TRUE),
    error = function(e) {
      cat(sprintf("Reference model failed for %s vs %s: %s\n", study1, study2, e$message))
      NULL
    }
  )

  if (is.null(ref_fit_formal) || !lavInspect(ref_fit_formal, "converged")) {
    cat(sprintf("Formal testing reference model did not converge for %s vs %s.\n",
                study1, study2))
    return(NULL)
  }

  ref_indices <- fitMeasures(ref_fit_formal, c("cfi", "tli"))
  cat(sprintf("Reference model: CFI = %.4f, TLI = %.4f\n\n",
              ref_indices["cfi"], ref_indices["tli"]))

  # Results table — one row per studied item
  dif_results <- data.frame(
    item           = studied_items,
    chi_sq         = NA_real_,
    df             = NA_real_,
    p_value        = NA_real_,
    delta_cfi      = NA_real_,
    delta_tli      = NA_real_,
    cfi_reference  = NA_real_,
    cfi_comparison = NA_real_,
    tli_reference  = NA_real_,
    tli_comparison = NA_real_,
    converged      = FALSE,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(studied_items)) {
    tested_item <- studied_items[i]

    cat(sprintf("Testing item %d/%d: %s\n", i, length(studied_items), tested_item))

    # Comparison model: same full item set, tested item additionally constrained
    comparison_model <- build_mgcfa_model(
      items        = full_item_set,
      anchor_items = c(tested_item, purified_anchors),
      ordinal      = TRUE,
      data         = pair_data
    )

    tryCatch({
      fit_comparison <- cfa(comparison_model, data = pair_data,
                            group = "group", ordered = TRUE,
                            estimator = "WLSMV", std.lv = TRUE)

      if (lavInspect(fit_comparison, "converged")) {

        # LRT: comparison (more constrained) vs reference (less constrained)
        lr_test <- lavTestLRT(fit_comparison, ref_fit_formal, method = "satorra.2000")

        dif_results$chi_sq[i]  <- lr_test$`Chisq diff`[2]
        dif_results$df[i]      <- lr_test$`Df diff`[2]
        dif_results$p_value[i] <- lr_test$`Pr(>Chisq)`[2]

        comp_indices <- fitMeasures(fit_comparison, c("cfi", "tli"))

        dif_results$cfi_reference[i]  <- ref_indices["cfi"]
        dif_results$cfi_comparison[i] <- comp_indices["cfi"]
        dif_results$tli_reference[i]  <- ref_indices["tli"]
        dif_results$tli_comparison[i] <- comp_indices["tli"]

        # ΔCFI = reference - comparison (positive = constraining item worsened fit = DIF)
        dif_results$delta_cfi[i] <- ref_indices["cfi"] - comp_indices["cfi"]
        dif_results$delta_tli[i] <- ref_indices["tli"] - comp_indices["tli"]
        dif_results$converged[i] <- TRUE

        cat(sprintf("  χ²(%.0f) = %.2f, p = %.4f, ΔCFI = %.4f, ΔTLI = %.4f\n",
                    dif_results$df[i], dif_results$chi_sq[i], dif_results$p_value[i],
                    dif_results$delta_cfi[i], dif_results$delta_tli[i]))
      } else {
        cat("  Comparison model did not converge\n")
      }
    }, error = function(e) {
      cat(sprintf("  Error: %s\n", e$message))
    })
  }

  # Filter to converged items
  dif_results <- dif_results %>%
    filter(converged, !is.na(p_value))

  if (nrow(dif_results) == 0) {
    cat("\nNo items converged successfully\n")
    return(NULL)
  }

  # Apply Benjamini-Hochberg correction
  dif_results <- dif_results %>%
    arrange(p_value) %>%
    mutate(
      rank = row_number(),
      bh_threshold = (rank / n()) * alpha,
      significant_BH = p_value < bh_threshold,
      # Also flag items with substantial fit differences
      substantial_cfi_diff = abs(delta_cfi) > cfi_threshold,
      substantial_tli_diff = abs(delta_tli) > tli_threshold,
      dif_by_fit_indices   = substantial_cfi_diff | substantial_tli_diff
    )

  n_sig_lr  <- sum(dif_results$significant_BH, na.rm = TRUE)
  n_sig_fit <- sum(dif_results$dif_by_fit_indices, na.rm = TRUE)

  cat(sprintf("\nItems with significant DIF (BH-adjusted α = %.2f): %d\n",
              alpha, n_sig_lr))
  cat(sprintf("Items with substantial fit differences (ΔCFI > %.3f or ΔTLI > %.3f): %d\n",
              cfi_threshold, tli_threshold, n_sig_fit))

  if (n_sig_lr > 0) {
    cat("\nSignificant DIF items (LR test):\n")
    sig_items <- dif_results %>% filter(significant_BH)
    for (i in 1:nrow(sig_items)) {
      cat(sprintf("  %s: p = %.4f (threshold = %.4f), ΔCFI = %.4f, ΔTLI = %.4f\n",
                  sig_items$item[i], sig_items$p_value[i],
                  sig_items$bh_threshold[i],
                  sig_items$delta_cfi[i], sig_items$delta_tli[i]))
    }
  }

  if (n_sig_fit > 0 && n_sig_fit != n_sig_lr) {
    cat("\nItems with substantial fit differences:\n")
    fit_items <- dif_results %>% filter(dif_by_fit_indices, !significant_BH)
    for (i in 1:nrow(fit_items)) {
      cat(sprintf("  %s: ΔCFI = %.4f, ΔTLI = %.4f (p = %.4f)\n",
                  fit_items$item[i],
                  fit_items$delta_cfi[i], fit_items$delta_tli[i],
                  fit_items$p_value[i]))
    }
  }
  cat("\n")

  return(list(
    dif_results      = dif_results,
    purified_anchors = purified_anchors,
    studied_items    = studied_items
  ))
}

# =============================================================================
# STEP 4: ASSESS DIF IMPACT WITH FACTOR SCORES
# =============================================================================


assess_dif_impact <- function(data, formal_results, study1, study2,
                              study_var = "Study") {
  
  cat(sprintf("\nSTEP 4: ASSESSING DIF IMPACT FOR %s vs %s\n", study1, study2))
  cat("=============================================\n\n")
  
  if (is.null(formal_results) || nrow(formal_results$dif_results) == 0) {
    cat("No DIF results available\n")
    return(NULL)
  }
  
  # Filter to relevant studies
  pair_data <- data %>%
    filter(!!sym(study_var) %in% c(study1, study2)) %>%
    mutate(group = as.factor(!!sym(study_var)))
  
  purified_anchors <- formal_results$purified_anchors
  all_items <- c(purified_anchors, formal_results$studied_items)
  
  # Identify items with significant DIF
  dif_items <- formal_results$dif_results %>%
    filter(significant_BH) %>%
    pull(item)
  
  if (length(dif_items) == 0) {
    cat("No significant DIF items found - impact assessment not needed\n")
    return(NULL)
  }
  
  cat("Comparing two models:\n")
  cat("  Model 1: All items constrained (ignore DIF)\n")
  cat("  Model 2: DIF items free across groups\n\n")
  
  # Model 1: All constrained
  model1 <- build_mgcfa_model(
    items = all_items,
    anchor_items = all_items,
    ordinal = TRUE,
    data = pair_data
  )

  fit1 <- cfa(model1, data = pair_data, group = "group",
              ordered = TRUE, estimator = "WLSMV", std.lv = TRUE)

  # Model 2: DIF items free
  model2 <- build_mgcfa_model(
    items = all_items,
    anchor_items = setdiff(all_items, dif_items),
    ordinal = TRUE,
    data = pair_data
  )
  
  fit2 <- cfa(model2, data = pair_data, group = "group",
              ordered = TRUE, estimator = "WLSMV", std.lv = TRUE)
  
  # Check convergence
  if (!lavInspect(fit1, "converged") || !lavInspect(fit2, "converged")) {
    cat("One or both models did not converge\n")
    return(NULL)
  }
  
  # Extract factor scores - lavPredict returns a list for multi-group models
  scores1_raw <- lavPredict(fit1)
  scores2_raw <- lavPredict(fit2)
  
  # Get case indices for each group
  case_idx1 <- lavInspect(fit1, "case.idx")
  case_idx2 <- lavInspect(fit2, "case.idx")
  
  # Process multi-group results
  if (is.list(scores1_raw)) {
    # Multi-group: combine across groups
    group_names <- names(scores1_raw)
    
    score_list1 <- list()
    score_list2 <- list()
    
    for (g in group_names) {
      # Model 1 scores for this group
      if (!is.null(scores1_raw[[g]]) && length(scores1_raw[[g]]) > 0) {
        group_scores1 <- data.frame(
          case_idx = case_idx1[[g]],
          score1 = as.numeric(scores1_raw[[g]]),
          group = g,
          stringsAsFactors = FALSE
        )
        score_list1[[g]] <- group_scores1
      }
      
      # Model 2 scores for this group
      if (!is.null(scores2_raw[[g]]) && length(scores2_raw[[g]]) > 0) {
        group_scores2 <- data.frame(
          case_idx = case_idx2[[g]],
          score2 = as.numeric(scores2_raw[[g]]),
          group = g,
          stringsAsFactors = FALSE
        )
        score_list2[[g]] <- group_scores2
      }
    }
    
    # Combine groups
    all_scores1 <- bind_rows(score_list1)
    all_scores2 <- bind_rows(score_list2)
    
  } else {
    # Single group (shouldn't happen with group= specified, but handle it)
    all_scores1 <- data.frame(
      case_idx = case_idx1,
      score1 = as.numeric(scores1_raw),
      stringsAsFactors = FALSE
    )
    
    all_scores2 <- data.frame(
      case_idx = case_idx2,
      score2 = as.numeric(scores2_raw),
      stringsAsFactors = FALSE
    )
  }
  
  # Merge scores from both models
  score_data <- all_scores1 %>%
    inner_join(all_scores2, by = c("case_idx", "group"))
  
  cat(sprintf("Valid comparisons: %d participants\n", nrow(score_data)))
  cat(sprintf("  Model 1 scored: %d participants\n", nrow(all_scores1)))
  cat(sprintf("  Model 2 scored: %d participants\n", nrow(all_scores2)))
  cat(sprintf("  Both models: %d participants\n\n", nrow(score_data)))
  
  if (nrow(score_data) < 100) {
    cat("WARNING: Very few participants with scores from both models\n")
    cat("Results may be unreliable\n\n")
  }
  
  # Calculate agreement metrics
  correlation <- cor(score_data$score1, score_data$score2, use = "complete.obs")
  
  # Mean absolute difference
  mad <- mean(abs(score_data$score1 - score_data$score2), na.rm = TRUE)
  
  # RMSE
  rmse <- sqrt(mean((score_data$score1 - score_data$score2)^2, na.rm = TRUE))
  
  # Intraclass correlation
  icc_data <- score_data %>% 
    filter(!is.na(score1) & !is.na(score2))
  
  if (nrow(icc_data) > 0 && var(icc_data$score1) > 0 && var(icc_data$score2) > 0) {
    icc <- 2 * cov(icc_data$score1, icc_data$score2) /
      (var(icc_data$score1) + var(icc_data$score2))
  } else {
    icc <- NA
  }
  
  cat("IMPACT ASSESSMENT RESULTS:\n")
  cat("==========================\n")
  cat(sprintf("Correlation: %.4f\n", correlation))
  cat(sprintf("Mean Absolute Difference: %.4f\n", mad))
  cat(sprintf("RMSE: %.4f\n", rmse))
  if (!is.na(icc)) {
    cat(sprintf("ICC: %.4f\n", icc))
  } else {
    cat("ICC: Could not be calculated\n")
  }
  cat("\n")
  
  # Interpretation
  if (correlation > 0.95 && mad < 0.2) {
    cat("Interpretation: DIF has MINIMAL impact on scores\n")
  } else if (correlation > 0.90 && mad < 0.5) {
    cat("Interpretation: DIF has MODERATE impact on scores\n")
  } else {
    cat("Interpretation: DIF has SUBSTANTIAL impact on scores\n")
  }
  cat("\n")
  
  # Create comparison plot
  p <- ggplot(score_data, aes(x = score1, y = score2, color = group)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = sprintf("Factor Score Comparison: %s vs %s", study1, study2),
         subtitle = sprintf("r = %.3f, MAD = %.3f, RMSE = %.3f (n = %d)", 
                            correlation, mad, rmse, nrow(score_data)),
         x = "Model 1: All Items Constrained",
         y = "Model 2: DIF Items Free") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p)
  
  return(list(
    score_data = score_data,
    correlation = correlation,
    mad = mad,
    rmse = rmse,
    icc = icc,
    comparison_plot = p,
    dif_items = dif_items,
    n_compared = nrow(score_data),
    n_model1 = nrow(all_scores1),
    n_model2 = nrow(all_scores2)
  ))
}

# =============================================================================
# HELPER FUNCTION: BUILD MGCFA MODEL SYNTAX
# =============================================================================

build_mgcfa_model <- function(items, anchor_items, ordinal = TRUE, data = NULL) {

  # Factor loadings
  model <- paste0("distress =~ ", paste(items, collapse = " + "), "\n")

  # Group invariance constraints for anchors
  if (length(anchor_items) > 0) {
    for (item in items) {
      if (item %in% anchor_items) {
        # Constrain loadings
        model <- paste0(model, "distress =~ c(lambda_", item, ", lambda_", item, ") * ", item, "\n")

        # If ordinal, constrain thresholds — number derived from observed categories
        if (ordinal) {
          if (!is.null(data) && item %in% names(data)) {
            item_values <- data[[item]]
            if ("group" %in% names(data)) {
              for (g in unique(data$group)) {
                g_vals <- item_values[data$group == g]
                if (length(unique(g_vals[!is.na(g_vals)])) < 2) {
                  warning(sprintf("Item '%s' has fewer than 2 observed categories in group '%s'", item, g))
                }
              }
            }
            n_cats      <- length(unique(item_values[!is.na(item_values)]))
            n_thresholds <- max(1L, n_cats - 1L)
          } else {
            n_thresholds <- 4L
          }
          for (t in seq_len(n_thresholds)) {
            model <- paste0(model, item, " | c(tau", t, "_", item, ", tau", t, "_", item, ") * t", t, "\n")
          }
        }
      }
    }
  }

  return(model)
}

# =============================================================================
# SUMMARY REPORT GENERATION
# =============================================================================

generate_summary_report <- function(all_results, coverage_analysis, output_dir) {
  
  cat("\nGenerating summary report...\n")
  
  # Create summary dataframe
  summary_df <- data.frame(
    comparison = character(),
    study1 = character(),
    study2 = character(),
    n_shared_anchors = integer(),
    n_purified_anchors = integer(),
    n_tested_items = integer(),
    n_sig_dif_lr = integer(),
    n_sig_dif_fit = integer(),
    prop_dif = numeric(),
    correlation = numeric(),
    mad = numeric(),
    rmse = numeric(),
    icc = numeric(),
    impact = character(),
    stringsAsFactors = FALSE
  )
  
  for (pair_name in names(all_results)) {
    result <- all_results[[pair_name]]
    
    n_sig_lr <- sum(result$dif_results$significant_BH, na.rm = TRUE)
    n_sig_fit <- sum(result$dif_results$dif_by_fit_indices, na.rm = TRUE)
    n_tested <- nrow(result$dif_results)
    
    # Determine impact level
    if (!is.null(result$impact_assessment)) {
      corr <- result$impact_assessment$correlation
      mad_val <- result$impact_assessment$mad
      
      if (corr > 0.95 && mad_val < 0.2) {
        impact <- "Minimal"
      } else if (corr > 0.90 && mad_val < 0.5) {
        impact <- "Moderate"
      } else {
        impact <- "Substantial"
      }
    } else {
      impact <- "Not assessed"
    }
    
    summary_df <- rbind(summary_df, data.frame(
      comparison = pair_name,
      study1 = result$study1,
      study2 = result$study2,
      n_shared_anchors = length(result$shared_anchors),
      n_purified_anchors = length(result$purified_anchors),
      n_tested_items = n_tested,
      n_sig_dif_lr = n_sig_lr,
      n_sig_dif_fit = n_sig_fit,
      prop_dif = ifelse(n_tested > 0, n_sig_lr / n_tested, NA),
      correlation = ifelse(!is.null(result$impact_assessment), 
                           result$impact_assessment$correlation, NA),
      mad = ifelse(!is.null(result$impact_assessment), 
                   result$impact_assessment$mad, NA),
      rmse = ifelse(!is.null(result$impact_assessment), 
                    result$impact_assessment$rmse, NA),
      icc = ifelse(!is.null(result$impact_assessment), 
                   result$impact_assessment$icc, NA),
      impact = impact,
      stringsAsFactors = FALSE
    ))
  }
  
  # Save summary
  write.csv(summary_df, 
            file.path(output_dir, "dif_summary_all_comparisons.csv"),
            row.names = FALSE)
  
  # Print summary statistics
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("SUMMARY ACROSS ALL PAIRWISE COMPARISONS\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat(sprintf("Total comparisons: %d\n", nrow(summary_df)))
  cat(sprintf("Mean items tested per comparison: %.1f\n", mean(summary_df$n_tested_items, na.rm = TRUE)))
  cat(sprintf("Mean significant DIF items: %.1f\n", mean(summary_df$n_sig_dif_lr, na.rm = TRUE)))
  cat(sprintf("Mean DIF proportion: %.2f%%\n", mean(summary_df$prop_dif, na.rm = TRUE) * 100))
  
  if (any(!is.na(summary_df$correlation))) {
    cat(sprintf("\nScore agreement (n=%d comparisons with impact data):\n", 
                sum(!is.na(summary_df$correlation))))
    cat(sprintf("  Mean correlation: %.3f\n", mean(summary_df$correlation, na.rm = TRUE)))
    cat(sprintf("  Mean MAD: %.3f\n", mean(summary_df$mad, na.rm = TRUE)))
    cat(sprintf("  Mean RMSE: %.3f\n", mean(summary_df$rmse, na.rm = TRUE)))
    
    # Impact breakdown
    impact_table <- table(summary_df$impact[summary_df$impact != "Not assessed"])
    cat("\nDIF Impact Distribution:\n")
    for (level in names(impact_table)) {
      cat(sprintf("  %s: %d comparisons\n", level, impact_table[level]))
    }
  }
  
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  cat("Summary report saved to:", file.path(output_dir, "dif_summary_all_comparisons.csv"), "\n")
  
  return(summary_df)
}

# =============================================================================
# MASTER PIPELINE FUNCTION
# =============================================================================

run_comprehensive_dif_pipeline <- function(data, anchor_items, 
                                           study_var = "Study",
                                           min_responses = 10,
                                           min_anchors = 3,
                                           alpha = 0.05,
                                           cfi_threshold = 0.01,
                                           tli_threshold = 0.01,
                                           output_dir = "dif_results") {
  
  cat("╔════════════════════════════════════════════════════════════╗\n")
  cat("║  COMPREHENSIVE DIF TESTING PIPELINE                        ║\n")
  cat("║  CFI/TLI-Based Approach for WLSMV Estimation               ║\n")
  cat("╚════════════════════════════════════════════════════════════╝\n\n")
  
  cat(sprintf("Settings:\n"))
  cat(sprintf("  CFI threshold: %.3f\n", cfi_threshold))
  cat(sprintf("  TLI threshold: %.3f\n", tli_threshold))
  cat(sprintf("  Minimum anchors: %d\n", min_anchors))
  cat(sprintf("  Alpha (BH correction): %.2f\n\n", alpha))
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Step 1: Identify anchor coverage
  coverage_analysis <- identify_anchor_coverage(data, anchor_items, 
                                                study_var, min_responses)
  
  # Initialize results storage
  all_results <- list()
  
  # Step 2-4: Process each pairwise comparison
  for (pair_name in names(coverage_analysis$pairwise_tests)) {
    
    cat("\n")
    cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
    cat(sprintf("PROCESSING: %s\n", pair_name))
    cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
    
    pair_info <- coverage_analysis$pairwise_tests[[pair_name]]
    study1 <- pair_info$study1
    study2 <- pair_info$study2
    shared_anchors <- pair_info$shared_anchors
    
    # Purification with CFI/TLI
    purified <- purify_anchors(data, shared_anchors, study1, study2, 
                               study_var, min_anchors, 
                               cfi_threshold, tli_threshold)
    
    if (is.null(purified)) {
      cat("Skipping formal testing due to purification failure\n")
      next
    }
    
    # Formal DIF testing
    formal <- formal_dif_testing(data, purified, study1, study2, 
                                 study_var, alpha,
                                 cfi_threshold, tli_threshold)
    
    if (is.null(formal)) {
      cat("Skipping impact assessment due to testing failure\n")
      next
    }
    
    # Impact assessment
    impact <- assess_dif_impact(data, formal, study1, study2, study_var)
    
    # Store results
    all_results[[pair_name]] <- list(
      study1 = study1,
      study2 = study2,
      shared_anchors = shared_anchors,
      purified_anchors = purified$purified_anchors,
      purification_details = purified$fit_differences,
      dif_results = formal$dif_results,
      impact_assessment = impact
    )
    
    # Save pairwise results
    write.csv(purified$fit_differences, 
              file.path(output_dir, paste0(pair_name, "_purification.csv")),
              row.names = FALSE)
    
    write.csv(formal$dif_results, 
              file.path(output_dir, paste0(pair_name, "_dif_results.csv")),
              row.names = FALSE)
    
    if (!is.null(impact)) {
      ggsave(file.path(output_dir, paste0(pair_name, "_score_comparison.png")),
             plot = impact$comparison_plot, width = 8, height = 6)
    }
  }
  
  # Generate summary report
  generate_summary_report(all_results, coverage_analysis, output_dir)
  
  cat("\n")
  cat("╔════════════════════════════════════════════════════════════╗\n")
  cat("║  PIPELINE COMPLETE                                         ║\n")
  cat(sprintf("║  Results saved to: %-40s║\n", output_dir))
  cat("╚════════════════════════════════════════════════════════════╝\n")
  
  return(list(
    coverage_analysis = coverage_analysis,
    pairwise_results = all_results
  ))
}

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# # Load your data
# load("df with LSAC duplicate SIDs removed.RData")

data_dif<- results[["data_preparation"]][["reference_wave_data"]]

# # Load anchor items from harmonization results
# # results <- harmonize_studies_reference_wave_concurrent(...)
# anchor_items <- results$anchor_analysis$all_anchors
# 
# # Run DIF pipeline with CFI/TLI approach
 dif_results <- run_comprehensive_dif_pipeline(
   data = data_dif,
   anchor_items = anchor_items,
   study_var = "Study",
   min_responses = 10,
   min_anchors = 3,
   alpha = 0.05,
   cfi_threshold = 0.01,  # Conventional threshold
   tli_threshold = 0.01,
   output_dir = "dif_cfi_tli_results"
 )


 

 