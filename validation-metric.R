# =============================================================================
# VALIDATION: HARMONIZED SCORES VS ORIGINAL INSTRUMENT SUM SCORES
# =============================================================================
# This code validates the IRT harmonized common metric by comparing it to
# the original raw sum scores from each distress measure within each study.
# High correlations provide evidence of construct validity.
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(corrplot)
library(gridExtra)

# =============================================================================
# STEP 1: IDENTIFY INSTRUMENTS USED IN EACH STUDY
# =============================================================================

identify_study_instruments <- function(data, study_var = "Study") {
  
  cat("STEP 1: IDENTIFYING INSTRUMENTS BY STUDY\n")
  cat("=========================================\n\n")
  
  # Define instrument patterns
  instrument_patterns <- list(
    PHQ9 = "^phq9i",
    GAD7 = "^gad7i",
    CESD = "^cesdi",
    K10 = "^k10i",
    SF36 = "^sf36i",
    GHQ12 = "^ghq12i",
    SDQ = "^sdqi",
    DASS = "^dassi",
    SMFQ = "^smfqi",
    CHQ = "^chqi",
    YSR = "^ysri",
    HBSC = "^hbsci",
    GADS = "^gadsi",
    SCAS = "^scasi",
    PedsQL = "^pedsqlgwi"
  )
  
  # Get all column names
  all_cols <- names(data)
  
  # Initialize results
  study_instruments <- list()
  
  studies <- unique(data[[study_var]])
  
  for (study in studies) {
    
    study_data <- data %>% filter(!!sym(study_var) == study)
    
    study_instruments[[study]] <- list()
    
    for (instrument_name in names(instrument_patterns)) {
      
      pattern <- instrument_patterns[[instrument_name]]
      
      # Find items matching this instrument
      instrument_items <- all_cols[grepl(pattern, all_cols, ignore.case = TRUE)]
      
      if (length(instrument_items) > 0) {
        
        # Check how many participants have data
        n_with_data <- study_data %>%
          select(all_of(instrument_items)) %>%
          filter(if_any(everything(), ~ !is.na(.x))) %>%
          nrow()
        
        # Check how many items have data
        items_with_data <- study_data %>%
          select(all_of(instrument_items)) %>%
          summarise(across(everything(), ~ sum(!is.na(.x)))) %>%
          select(where(~ .x >= 10)) %>%
          names()
        
        if (length(items_with_data) >= 3 && n_with_data >= 10) {
          study_instruments[[study]][[instrument_name]] <- list(
            items = items_with_data,
            n_items = length(items_with_data),
            n_participants = n_with_data
          )
        }
      }
    }
  }
  
  # Print summary
  cat("Instruments identified by study:\n")
  cat("================================\n\n")
  
  for (study in names(study_instruments)) {
    instruments <- names(study_instruments[[study]])
    if (length(instruments) > 0) {
      cat(sprintf("%s: %s\n", study, paste(instruments, collapse = ", ")))
      for (inst in instruments) {
        cat(sprintf("  - %s: %d items, %d participants\n", 
                    inst, 
                    study_instruments[[study]][[inst]]$n_items,
                    study_instruments[[study]][[inst]]$n_participants))
      }
    } else {
      cat(sprintf("%s: No instruments identified\n", study))
    }
    cat("\n")
  }
  
  return(study_instruments)
}

# =============================================================================
# STEP 2: CALCULATE SUM SCORES FOR EACH INSTRUMENT IN EACH STUDY
# =============================================================================

calculate_instrument_sum_scores <- function(data, study_instruments, 
                                            study_var = "Study",
                                            min_items_for_score = 0.5) {
  
  cat("STEP 2: CALCULATING INSTRUMENT SUM SCORES\n")
  cat("==========================================\n\n")
  
  # Add ID columns if they exist
  id_cols <- c("Study", "Wave", "GID", "SID", "Study_Wave")
  id_cols <- id_cols[id_cols %in% names(data)]
  
  sum_scores_list <- list()
  
  for (study in names(study_instruments)) {
    
    study_data <- data %>% filter(!!sym(study_var) == study)
    
    # Start with ID columns
    study_sum_scores <- study_data %>%
      select(all_of(id_cols))
    
    instruments <- study_instruments[[study]]
    
    if (length(instruments) == 0) {
      cat(sprintf("%s: No instruments to score\n", study))
      next
    }
    
    cat(sprintf("%s:\n", study))
    
    for (instrument_name in names(instruments)) {
      
      items <- instruments[[instrument_name]]$items
      n_items <- length(items)
      min_required <- ceiling(n_items * min_items_for_score)
      
      cat(sprintf("  %s: %d items (min %d required)\n", 
                  instrument_name, n_items, min_required))
      
      # Calculate sum score
      item_data <- study_data %>%
        select(all_of(items)) %>%
        mutate(across(everything(), as.numeric))
      
      # Count valid responses
      n_valid <- rowSums(!is.na(item_data))
      
      # Calculate sum (only if minimum items present)
      sum_score <- rowSums(item_data, na.rm = TRUE)
      sum_score[n_valid < min_required] <- NA
      
      # Calculate mean score (average of available items)
      mean_score <- rowMeans(item_data, na.rm = TRUE)
      mean_score[n_valid < min_required] <- NA
      
      # Calculate proportion maximum possible (POMP)
      # Assumes items are scored from min to max observed in data
      min_vals <- sapply(item_data, min, na.rm = TRUE)
      max_vals <- sapply(item_data, max, na.rm = TRUE)
      
      # Calculate actual range for each item
      ranges <- max_vals - min_vals
      
      # POMP score
      pomp_scores <- sapply(1:nrow(item_data), function(i) {
        row_data <- as.numeric(item_data[i, ])
        if (sum(!is.na(row_data)) < min_required) return(NA)
        
        valid_items <- which(!is.na(row_data))
        item_scores <- row_data[valid_items]
        item_mins <- min_vals[valid_items]
        item_ranges <- ranges[valid_items]
        
        # Normalize each item to 0-1 scale
        normalized <- (item_scores - item_mins) / item_ranges
        
        # Average normalized scores
        mean(normalized, na.rm = TRUE) * 100
      })
      
      # Add to dataframe
      col_prefix <- paste0(instrument_name, "_")
      study_sum_scores[[paste0(col_prefix, "sum")]] <- sum_score
      study_sum_scores[[paste0(col_prefix, "mean")]] <- mean_score
      study_sum_scores[[paste0(col_prefix, "pomp")]] <- pomp_scores
      study_sum_scores[[paste0(col_prefix, "n_items")]] <- n_valid
      
      # Summary
      n_scored <- sum(!is.na(sum_score))
      cat(sprintf("    Scored: %d participants (%.1f%%)\n", 
                  n_scored, n_scored / nrow(study_data) * 100))
      if (n_scored > 0) {
        cat(sprintf("    Sum range: %.1f - %.1f (M=%.1f, SD=%.1f)\n",
                    min(sum_score, na.rm = TRUE),
                    max(sum_score, na.rm = TRUE),
                    mean(sum_score, na.rm = TRUE),
                    sd(sum_score, na.rm = TRUE)))
        cat(sprintf("    POMP range: %.1f - %.1f (M=%.1f, SD=%.1f)\n",
                    min(pomp_scores, na.rm = TRUE),
                    max(pomp_scores, na.rm = TRUE),
                    mean(pomp_scores, na.rm = TRUE),
                    sd(pomp_scores, na.rm = TRUE)))
      }
    }
    
    sum_scores_list[[study]] <- study_sum_scores
    cat("\n")
  }
  
  # Combine all studies
  all_sum_scores <- bind_rows(sum_scores_list)
  
  cat(sprintf("Total sum scores generated: %d participants\n\n", nrow(all_sum_scores)))
  
  return(all_sum_scores)
}

# =============================================================================
# STEP 3: MERGE WITH HARMONIZED SCORES
# =============================================================================

merge_harmonized_and_sum_scores <- function(harmonized_scores, sum_scores,
                                            merge_keys = c("Study", "Wave", "GID", "SID")) {
  
  cat("STEP 3: MERGING HARMONIZED AND SUM SCORES\n")
  cat("==========================================\n\n")
  
  # Identify available merge keys
  available_keys <- merge_keys[merge_keys %in% names(harmonized_scores) & 
                                 merge_keys %in% names(sum_scores)]
  
  if (length(available_keys) == 0) {
    stop("No common merge keys found between harmonized_scores and sum_scores")
  }
  
  cat("Merge keys:", paste(available_keys, collapse = ", "), "\n")
  
  # Merge
  merged <- harmonized_scores %>%
    inner_join(sum_scores, by = available_keys, suffix = c("", "_sum"))
  
  cat(sprintf("Harmonized scores: %d participants\n", nrow(harmonized_scores)))
  cat(sprintf("Sum scores: %d participants\n", nrow(sum_scores)))
  cat(sprintf("Successfully merged: %d participants\n", nrow(merged)))
  cat(sprintf("Match rate: %.1f%%\n\n", nrow(merged) / nrow(harmonized_scores) * 100))
  
  return(merged)
}

# =============================================================================
# STEP 4: EXAMINE CORRELATIONS WITHIN EACH STUDY
# =============================================================================

examine_correlations_by_study <- function(merged_data, study_instruments,
                                          study_var = "Study",
                                          theta_var = "Theta_Standardized") {
  
  cat("STEP 4: EXAMINING CORRELATIONS BY STUDY\n")
  cat("========================================\n\n")
  
  correlation_results <- list()
  
  for (study in names(study_instruments)) {
    
    study_data <- merged_data %>% filter(!!sym(study_var) == study)
    
    if (nrow(study_data) == 0) {
      cat(sprintf("%s: No merged data\n\n", study))
      next
    }
    
    instruments <- names(study_instruments[[study]])
    
    if (length(instruments) == 0) {
      cat(sprintf("%s: No instruments\n\n", study))
      next
    }
    
    cat(sprintf("%s (n=%d):\n", study, nrow(study_data)))
    cat(paste(rep("=", 60), collapse = ""), "\n")
    
    study_correlations <- data.frame(
      study = character(),
      instrument = character(),
      score_type = character(),
      n = integer(),
      correlation = numeric(),
      r_squared = numeric(),
      rmse = numeric(),
      mae = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (instrument in instruments) {
      
      cat(sprintf("\n%s:\n", instrument))
      cat(paste(rep("-", 40), collapse = ""), "\n")
      
      # Get score columns for this instrument
      sum_col <- paste0(instrument, "_sum")
      mean_col <- paste0(instrument, "_mean")
      pomp_col <- paste0(instrument, "_pomp")
      
      theta_values <- study_data[[theta_var]]
      
      # Sum score correlation
      if (sum_col %in% names(study_data)) {
        sum_values <- study_data[[sum_col]]
        valid <- !is.na(theta_values) & !is.na(sum_values)
        
        if (sum(valid) >= 10) {
          cor_sum <- cor(theta_values[valid], sum_values[valid])
          r2_sum <- cor_sum^2
          
          # RMSE and MAE (with standardized theta)
          rmse_sum <- sqrt(mean((scale(theta_values[valid])[,1] - 
                                   scale(sum_values[valid])[,1])^2))
          mae_sum <- mean(abs(scale(theta_values[valid])[,1] - 
                                scale(sum_values[valid])[,1]))
          
          cat(sprintf("Sum Score:  r = %.3f (R² = %.3f), RMSE = %.3f, MAE = %.3f (n=%d)\n",
                      cor_sum, r2_sum, rmse_sum, mae_sum, sum(valid)))
          
          study_correlations <- rbind(study_correlations, data.frame(
            study = study,
            instrument = instrument,
            score_type = "sum",
            n = sum(valid),
            correlation = cor_sum,
            r_squared = r2_sum,
            rmse = rmse_sum,
            mae = mae_sum,
            stringsAsFactors = FALSE
          ))
        }
      }
      
      # Mean score correlation
      if (mean_col %in% names(study_data)) {
        mean_values <- study_data[[mean_col]]
        valid <- !is.na(theta_values) & !is.na(mean_values)
        
        if (sum(valid) >= 10) {
          cor_mean <- cor(theta_values[valid], mean_values[valid])
          r2_mean <- cor_mean^2
          
          rmse_mean <- sqrt(mean((scale(theta_values[valid])[,1] - 
                                    scale(mean_values[valid])[,1])^2))
          mae_mean <- mean(abs(scale(theta_values[valid])[,1] - 
                                 scale(mean_values[valid])[,1]))
          
          cat(sprintf("Mean Score: r = %.3f (R² = %.3f), RMSE = %.3f, MAE = %.3f (n=%d)\n",
                      cor_mean, r2_mean, rmse_mean, mae_mean, sum(valid)))
          
          study_correlations <- rbind(study_correlations, data.frame(
            study = study,
            instrument = instrument,
            score_type = "mean",
            n = sum(valid),
            correlation = cor_mean,
            r_squared = r2_mean,
            rmse = rmse_mean,
            mae = mae_mean,
            stringsAsFactors = FALSE
          ))
        }
      }
      
      # POMP score correlation
      if (pomp_col %in% names(study_data)) {
        pomp_values <- study_data[[pomp_col]]
        valid <- !is.na(theta_values) & !is.na(pomp_values)
        
        if (sum(valid) >= 10) {
          cor_pomp <- cor(theta_values[valid], pomp_values[valid])
          r2_pomp <- cor_pomp^2
          
          rmse_pomp <- sqrt(mean((scale(theta_values[valid])[,1] - 
                                    scale(pomp_values[valid])[,1])^2))
          mae_pomp <- mean(abs(scale(theta_values[valid])[,1] - 
                                 scale(pomp_values[valid])[,1]))
          
          cat(sprintf("POMP Score: r = %.3f (R² = %.3f), RMSE = %.3f, MAE = %.3f (n=%d)\n",
                      cor_pomp, r2_pomp, rmse_pomp, mae_pomp, sum(valid)))
          
          study_correlations <- rbind(study_correlations, data.frame(
            study = study,
            instrument = instrument,
            score_type = "pomp",
            n = sum(valid),
            correlation = cor_pomp,
            r_squared = r2_pomp,
            rmse = rmse_pomp,
            mae = mae_pomp,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    correlation_results[[study]] <- study_correlations
    cat("\n")
  }
  
  # Combine all results
  all_correlations <- bind_rows(correlation_results)
  
  # Overall summary
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("OVERALL SUMMARY\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  if (nrow(all_correlations) > 0) {
    
    # By score type
    cat("By Score Type:\n")
    score_summary <- all_correlations %>%
      group_by(score_type) %>%
      summarise(
        n_comparisons = n(),
        mean_r = mean(correlation, na.rm = TRUE),
        median_r = median(correlation, na.rm = TRUE),
        min_r = min(correlation, na.rm = TRUE),
        max_r = max(correlation, na.rm = TRUE),
        mean_r2 = mean(r_squared, na.rm = TRUE),
        .groups = "drop"
      )
    print(score_summary)
    cat("\n")
    
    # By instrument
    cat("By Instrument:\n")
    inst_summary <- all_correlations %>%
      filter(score_type == "sum") %>%  # Use sum scores for this summary
      group_by(instrument) %>%
      summarise(
        n_studies = n(),
        mean_r = mean(correlation, na.rm = TRUE),
        median_r = median(correlation, na.rm = TRUE),
        min_r = min(correlation, na.rm = TRUE),
        max_r = max(correlation, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(mean_r))
    print(inst_summary)
    cat("\n")
    
    # Quality assessment
    cat("Quality Assessment:\n")
    cat(sprintf("  Excellent (r > 0.90): %d (%.1f%%)\n",
                sum(all_correlations$correlation > 0.90),
                mean(all_correlations$correlation > 0.90) * 100))
    cat(sprintf("  Good (r > 0.80): %d (%.1f%%)\n",
                sum(all_correlations$correlation > 0.80),
                mean(all_correlations$correlation > 0.80) * 100))
    cat(sprintf("  Acceptable (r > 0.70): %d (%.1f%%)\n",
                sum(all_correlations$correlation > 0.70),
                mean(all_correlations$correlation > 0.70) * 100))
    cat(sprintf("  Moderate (r > 0.60): %d (%.1f%%)\n",
                sum(all_correlations$correlation > 0.60),
                mean(all_correlations$correlation > 0.60) * 100))
    cat(sprintf("  Weak (r ≤ 0.60): %d (%.1f%%)\n",
                sum(all_correlations$correlation <= 0.60),
                mean(all_correlations$correlation <= 0.60) * 100))
    
  } else {
    cat("No correlations computed\n")
  }
  
  return(all_correlations)
}

# =============================================================================
# STEP 4b: PAIRWISE INTER-INSTRUMENT CORRELATIONS WITHIN STUDIES
# =============================================================================

examine_inter_instrument_correlations <- function(sum_scores,
                                                  study_instruments,
                                                  study_var = "Study",
                                                  min_n = 30) {

  cat("\nSTEP 4b: INTER-INSTRUMENT CORRELATIONS WITHIN STUDIES\n")
  cat("======================================================\n\n")

  inter_instrument_results <- data.frame(
    study        = character(),
    instrument1  = character(),
    instrument2  = character(),
    n            = integer(),
    correlation  = numeric(),
    r_squared    = numeric(),
    stringsAsFactors = FALSE
  )

  for (study in names(study_instruments)) {

    study_rows <- sum_scores[[study_var]] == study
    study_ss   <- sum_scores[study_rows, , drop = FALSE]

    # Identify instruments with a valid _pomp column for this study
    available_instruments <- character(0)
    for (inst in names(study_instruments[[study]])) {
      pomp_col <- paste0(inst, "_pomp")
      if (pomp_col %in% names(study_ss)) {
        n_valid <- sum(!is.na(study_ss[[pomp_col]]))
        if (n_valid >= min_n) {
          available_instruments <- c(available_instruments, inst)
        }
      }
    }

    if (length(available_instruments) < 2) {
      cat(sprintf("%s: only %d instrument(s) available, skipping\n",
                  study, length(available_instruments)))
      next
    }

    # Build POMP matrix — rows with at least one non-NA value
    pomp_cols   <- paste0(available_instruments, "_pomp")
    pomp_matrix <- as.matrix(study_ss[, pomp_cols, drop = FALSE])
    colnames(pomp_matrix) <- available_instruments
    pomp_matrix <- pomp_matrix[rowSums(!is.na(pomp_matrix)) > 0, , drop = FALSE]

    # Pairwise correlations (pairwise complete observations)
    cor_matrix <- cor(pomp_matrix, use = "pairwise.complete.obs")

    # Sample sizes per pair
    n_inst   <- ncol(pomp_matrix)
    n_matrix <- matrix(NA_integer_, n_inst, n_inst,
                       dimnames = list(colnames(pomp_matrix), colnames(pomp_matrix)))
    for (i in seq_len(n_inst)) {
      for (j in seq_len(n_inst)) {
        n_matrix[i, j] <- sum(complete.cases(pomp_matrix[, c(i, j)]))
      }
    }

    # Extract upper triangle into long format, instrument1 < instrument2 alphabetically
    inst_names  <- colnames(pomp_matrix)
    study_pairs <- data.frame(
      study       = character(),
      instrument1 = character(),
      instrument2 = character(),
      n           = integer(),
      correlation = numeric(),
      r_squared   = numeric(),
      stringsAsFactors = FALSE
    )

    for (i in seq_len(n_inst - 1)) {
      for (j in seq(i + 1, n_inst)) {
        pair        <- sort(c(inst_names[i], inst_names[j]))
        study_pairs <- rbind(study_pairs, data.frame(
          study       = study,
          instrument1 = pair[1],
          instrument2 = pair[2],
          n           = n_matrix[i, j],
          correlation = cor_matrix[i, j],
          r_squared   = cor_matrix[i, j]^2,
          stringsAsFactors = FALSE
        ))
      }
    }

    inter_instrument_results <- rbind(inter_instrument_results, study_pairs)

    # Print study results
    cat(sprintf("\n%s (instruments: %s)\n", study,
                paste(available_instruments, collapse = ", ")))
    cat(paste(rep("-", 50), collapse = ""), "\n")
    for (k in seq_len(nrow(study_pairs))) {
      cat(sprintf("  %s vs %s: r = %.3f (n = %d)\n",
                  study_pairs$instrument1[k], study_pairs$instrument2[k],
                  study_pairs$correlation[k], study_pairs$n[k]))
    }
    cat(sprintf("  Mean inter-instrument r: %.3f\n",
                mean(study_pairs$correlation, na.rm = TRUE)))
  }

  # Overall summary across all studies
  cat("\nINTER-INSTRUMENT CORRELATION SUMMARY (across all studies)\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")

  if (nrow(inter_instrument_results) > 0) {
    pair_summary <- inter_instrument_results %>%
      group_by(instrument1, instrument2) %>%
      summarise(
        n_studies = n(),
        mean_r    = mean(correlation, na.rm = TRUE),
        median_r  = median(correlation, na.rm = TRUE),
        min_r     = min(correlation, na.rm = TRUE),
        max_r     = max(correlation, na.rm = TRUE),
        .groups   = "drop"
      ) %>%
      arrange(desc(mean_r))
    print(pair_summary)
  } else {
    cat("No inter-instrument pairs found\n")
    pair_summary <- data.frame(
      instrument1 = character(), instrument2 = character(),
      n_studies   = integer(),   mean_r      = numeric(),
      median_r    = numeric(),   min_r       = numeric(),
      max_r       = numeric(),
      stringsAsFactors = FALSE
    )
  }

  return(list(
    pair_correlations = inter_instrument_results,
    pair_summary      = pair_summary
  ))
}

# =============================================================================
# STEP 4c: BACK-TRANSLATION METRICS (OBSERVED vs IRT-EXPECTED SUM SCORES)
# =============================================================================

compute_backtranslation_metrics <- function(
  data,
  raw_data,
  model,
  model_items,
  study_instruments,
  study_var          = "Study",
  theta_var          = "Theta_Score",
  target_instruments = c("K10", "SF36", "GHQ12", "SDQ"),
  min_items_in_model = 3,
  min_n_instrument   = 30
) {

  cat("\nBACK-TRANSLATION: Observed vs IRT-Expected Sum Scores\n")
  cat(paste(rep("-", 60), collapse = ""), "\n\n")

  # data (merged_data from Step 3) contains theta_var but not raw item columns.
  # raw_data is the original item-level dataset which has item columns but may
  # not have theta_var. We join theta from data onto raw_data below.
  if (!theta_var %in% names(data)) {
    available_numeric <- names(data)[sapply(data, is.numeric)]
    cat(sprintf(
      "  ERROR: theta_var '%s' not found in data.\n  Available numeric columns: %s\n",
      theta_var,
      paste(head(available_numeric, 10), collapse = ", ")
    ))
    return(NULL)
  }

  if (!is.data.frame(raw_data) || nrow(raw_data) == 0) {
    cat("  ERROR: raw_data is missing or empty.\n")
    return(NULL)
  }

  if (!requireNamespace("mirt", quietly = TRUE)) {
    stop("mirt package required for back-translation metrics")
  }

  # Check whether raw_data already contains theta_var (common when the
  # caller passes a pre-merged object that has both item columns and theta).
  # If so, use raw_data directly — no join needed and no collision risk.
  # Only join from data if theta_var is genuinely absent from raw_data.

  if (theta_var %in% names(raw_data)) {

    cat(sprintf(
      "  theta_var '%s' already present in raw_data — using directly, no join required.\n",
      theta_var))

    working_data <- raw_data %>%
      filter(!is.na(.data[[theta_var]]))

    cat(sprintf(
      "  Working data: %d rows with valid theta\n",
      nrow(working_data)))

  } else {

    cat(sprintf(
      "  theta_var '%s' not in raw_data — joining from data.\n",
      theta_var))

    possible_keys  <- c("Study", "Wave", "GID", "SID")
    available_keys <- possible_keys[
      possible_keys %in% names(raw_data) &
      possible_keys %in% names(data)
    ]

    if (length(available_keys) == 0) {
      cat("  ERROR: No common join keys between raw_data and data.\n")
      return(NULL)
    }

    # Select ONLY keys + theta from data to prevent any column collision
    theta_to_join <- data %>%
      select(all_of(c(available_keys, theta_var))) %>%
      distinct(across(all_of(available_keys)), .keep_all = TRUE)

    working_data <- raw_data %>%
      left_join(theta_to_join,
                by           = available_keys,
                relationship = "many-to-one") %>%
      filter(!is.na(.data[[theta_var]]))

    cat(sprintf(
      "  Working data: %d rows with valid theta (joined from data)\n",
      nrow(working_data)))

  }

  if (nrow(working_data) == 0) {
    cat(sprintf(
      "  ERROR: No rows with valid theta. theta_var = '%s'\n",
      theta_var))
    cat(sprintf(
      "  Columns in working_data: %s\n",
      paste(head(names(working_data), 15), collapse = ", ")))
    return(NULL)
  }

  # Spot-check that item columns survived
  item_pattern_present <- any(grepl(
    "^k10i|^sf36i|^ghq12i|^sdqi", names(working_data)))
  if (!item_pattern_present) {
    cat("  WARNING: No target item columns found in working_data.\n")
    cat(sprintf("  First 20 columns: %s\n",
                paste(head(names(working_data), 20), collapse = ", ")))
  }

  # Map item names to 1-based model column indices for mirt::expected.test()
  item_index_lookup <- setNames(seq_along(model_items), model_items)

  person_level <- data.frame(
    study         = character(),
    instrument    = character(),
    theta         = numeric(),
    observed_sum  = numeric(),
    expected_sum  = numeric(),
    n_model_items = integer(),
    stringsAsFactors = FALSE
  )

  for (study in names(study_instruments)) {
    for (instrument in target_instruments) {

      if (!instrument %in% names(study_instruments[[study]])) {
        cat(sprintf("  %s / %s: not administered, skipping\n", study, instrument))
        next
      }

      all_items        <- study_instruments[[study]][[instrument]]$items
      model_inst_items <- intersect(all_items, model_items)

      if (length(model_inst_items) < min_items_in_model) {
        cat(sprintf("  %s / %s: only %d item(s) in model, skipping\n",
                    study, instrument, length(model_inst_items)))
        next
      }

      item_indices <- as.integer(item_index_lookup[model_inst_items])

      sd <- working_data %>% filter(!!sym(study_var) == study)

      # Complete cases on model item subset only — ensures observed and
      # expected scores span exactly the same item set (no partial sums)
      response_mat  <- sd[, model_inst_items, drop = FALSE]
      complete_rows <- complete.cases(response_mat)
      sd_complete   <- sd[complete_rows, , drop = FALSE]

      if (nrow(sd_complete) < min_n_instrument) {
        cat(sprintf("  %s / %s: only %d complete cases, skipping\n",
                    study, instrument, nrow(sd_complete)))
        next
      }

      # Vectorised expected scores: build lookup on unique rounded thetas to
      # avoid calling expected.test() once per person (slow for large studies)
      theta_vals    <- sd_complete[[theta_var]]
      unique_thetas <- unique(round(theta_vals, 4))

      expected_lookup <- tryCatch(
        setNames(
          sapply(unique_thetas, function(th)
            mirt::expected.test(
              model,
              Theta       = matrix(th, nrow = 1),
              which.items = item_indices
            )
          ),
          as.character(unique_thetas)
        ),
        error = function(e) {
          cat(sprintf("  %s / %s: expected.test() failed: %s\n",
                      study, instrument, e$message))
          NULL
        }
      )
      if (is.null(expected_lookup)) next

      expected_vals <- expected_lookup[as.character(round(theta_vals, 4))]
      observed_vals <- rowSums(
        as.matrix(sd_complete[, model_inst_items]),
        na.rm = FALSE  # complete cases — no NAs present
      )

      cat(sprintf(
        "  %s / %s: %d complete cases, %d model items of %d instrument items\n",
        study, instrument, nrow(sd_complete),
        length(model_inst_items), length(all_items)))

      person_level <- rbind(person_level, data.frame(
        study         = study,
        instrument    = instrument,
        theta         = theta_vals,
        observed_sum  = observed_vals,
        expected_sum  = as.numeric(expected_vals),
        n_model_items = length(model_inst_items),
        stringsAsFactors = FALSE
      ))
    }
  }

  if (nrow(person_level) == 0) {
    cat("  No back-translation data generated.\n")
    return(NULL)
  }

  # Pooled metrics per instrument (primary output)
  pooled_metrics <- person_level %>%
    group_by(instrument) %>%
    summarise(
      n        = n(),
      n_studies = n_distinct(study),
      r        = cor(observed_sum, expected_sum),
      rmse     = sqrt(mean((observed_sum - expected_sum)^2)),
      mad      = mean(abs(observed_sum - expected_sum)),
      bias     = mean(observed_sum - expected_sum),
      obs_mean = mean(observed_sum),
      obs_sd   = sd(observed_sum),
      exp_mean = mean(expected_sum),
      exp_sd   = sd(expected_sum),
      .groups  = "drop"
    )

  # Per-study-instrument metrics (secondary, for heterogeneity inspection)
  study_metrics <- person_level %>%
    group_by(study, instrument) %>%
    summarise(
      n        = n(),
      r        = cor(observed_sum, expected_sum),
      rmse     = sqrt(mean((observed_sum - expected_sum)^2)),
      mad      = mean(abs(observed_sum - expected_sum)),
      bias     = mean(observed_sum - expected_sum),
      obs_mean = mean(observed_sum),
      obs_sd   = sd(observed_sum),
      exp_mean = mean(expected_sum),
      exp_sd   = sd(expected_sum),
      .groups  = "drop"
    )

  # Print pooled summary
  cat("\nBACK-TRANSLATION METRICS: OBSERVED vs EXPECTED SUM SCORES\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("Pooled across studies, restricted to model item subset\n\n")

  for (i in seq_len(nrow(pooled_metrics))) {
    inst <- pooled_metrics$instrument[i]
    n_mi <- max(person_level$n_model_items[person_level$instrument == inst])
    cat(sprintf("%s (n=%d across %d studies, %d model items):\n",
                inst, pooled_metrics$n[i], pooled_metrics$n_studies[i], n_mi))
    cat(sprintf("  r (observed vs expected):  %.3f\n", pooled_metrics$r[i]))
    cat(sprintf("  RMSE (original scale):     %.3f\n", pooled_metrics$rmse[i]))
    cat(sprintf("  MAD  (original scale):     %.3f\n", pooled_metrics$mad[i]))
    cat(sprintf("  Bias (obs - exp):          %.3f\n", pooled_metrics$bias[i]))
    cat(sprintf("  Observed  M (SD): %.2f (%.2f)\n",
                pooled_metrics$obs_mean[i], pooled_metrics$obs_sd[i]))
    cat(sprintf("  Expected  M (SD): %.2f (%.2f)\n",
                pooled_metrics$exp_mean[i], pooled_metrics$exp_sd[i]))
    cat("\n")
  }

  return(list(
    person_level   = person_level,
    pooled_metrics = pooled_metrics,
    study_metrics  = study_metrics
  ))
}

# =============================================================================
# STEP 5: CREATE VISUALIZATIONS
# =============================================================================

create_validation_plots <- function(merged_data, correlation_results,
                                    study_instruments, study_var = "Study",
                                    theta_var = "Theta_Standardized",
                                    inter_instrument_results = NULL,
                                    backtranslation_results = NULL) {
  
  cat("\nSTEP 5: CREATING VISUALIZATIONS\n")
  cat("================================\n\n")
  
  plots <- list()
  
  # 1. Overall correlation heatmap
  if (nrow(correlation_results) > 0) {
    
    # Reshape for heatmap
    heatmap_data <- correlation_results %>%
      filter(score_type == "sum") %>%
      select(study, instrument, correlation) %>%
      pivot_wider(names_from = instrument, values_from = correlation, values_fill = NA)
    
    if (ncol(heatmap_data) > 1) {
      
      # Convert to matrix
      heatmap_matrix <- as.matrix(heatmap_data[, -1])
      rownames(heatmap_matrix) <- heatmap_data$study
      
      plots$heatmap <- corrplot(heatmap_matrix, 
                                method = "color",
                                type = "full",
                                is.corr = FALSE,
                                cl.lim = c(0, 1),
                                col = colorRampPalette(c("white", "lightblue", "blue"))(100),
                                addCoef.col = "black",
                                number.cex = 0.7,
                                tl.col = "black",
                                tl.srt = 45,
                                title = "Correlations: Harmonized Theta vs Original Sum Scores",
                                mar = c(0, 0, 2, 0))
      
      cat("Created correlation heatmap\n")
    }
  }
  
  # 2. Distribution of correlations
  if (nrow(correlation_results) > 0) {
    
    p_dist <- ggplot(correlation_results %>% filter(score_type == "sum"), 
                     aes(x = correlation)) +
      geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black", alpha = 0.7) +
      geom_vline(xintercept = c(0.70, 0.80, 0.90), linetype = "dashed", color = "red") +
      labs(title = "Distribution of Correlations (Sum Scores)",
           subtitle = "Red lines at r = 0.70, 0.80, 0.90",
           x = "Correlation (r)",
           y = "Count") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
    
    plots$distribution <- p_dist
    cat("Created correlation distribution plot\n")
  }
  
  # 3. Correlations by instrument
  if (nrow(correlation_results) > 0) {
    
    inst_plot_data <- correlation_results %>%
      filter(score_type == "sum") %>%
      group_by(instrument) %>%
      summarise(
        mean_r = mean(correlation, na.rm = TRUE),
        se_r = sd(correlation, na.rm = TRUE) / sqrt(n()),
        n_studies = n(),
        .groups = "drop"
      ) %>%
      arrange(desc(mean_r))
    
    p_inst <- ggplot(inst_plot_data, aes(x = reorder(instrument, mean_r), y = mean_r)) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
      geom_errorbar(aes(ymin = mean_r - se_r, ymax = mean_r + se_r), 
                    width = 0.2) +
      geom_hline(yintercept = c(0.70, 0.80, 0.90), linetype = "dashed", 
                 color = "red", alpha = 0.5) +
      geom_text(aes(label = paste0("n=", n_studies)), 
                vjust = -0.5, size = 3) +
      coord_flip() +
      labs(title = "Mean Correlation by Instrument",
           subtitle = "Error bars show SE; dashed lines at r = 0.70, 0.80, 0.90",
           x = "Instrument",
           y = "Mean Correlation (r)") +
      ylim(0, 1) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
    
    plots$by_instrument <- p_inst
    cat("Created by-instrument plot\n")
  }
  
  # 4. Scatterplots for each study-instrument combination
  scatterplots <- list()
  
  for (study in unique(correlation_results$study)) {
    
    study_data <- merged_data %>% filter(!!sym(study_var) == study)
    instruments <- names(study_instruments[[study]])
    
    for (instrument in instruments) {
      
      sum_col <- paste0(instrument, "_sum")
      
      if (sum_col %in% names(study_data)) {
        
        plot_data <- study_data %>%
          filter(!is.na(!!sym(theta_var)) & !is.na(!!sym(sum_col))) %>%
          select(theta = !!sym(theta_var), sum_score = !!sym(sum_col))
        
        if (nrow(plot_data) >= 10) {
          
          cor_val <- cor(plot_data$theta, plot_data$sum_score)
          
          p_scatter <- ggplot(plot_data, aes(x = sum_score, y = theta)) +
            geom_point(alpha = 0.4, color = "steelblue") +
            geom_smooth(method = "lm", color = "red", se = TRUE) +
            labs(title = paste0(study, " - ", instrument),
                 subtitle = sprintf("r = %.3f, n = %d", cor_val, nrow(plot_data)),
                 x = paste(instrument, "Sum Score"),
                 y = "Harmonized Theta (Standardized)") +
            theme_minimal() +
            theme(plot.title = element_text(face = "bold"))
          
          plot_name <- paste0(study, "_", instrument)
          scatterplots[[plot_name]] <- p_scatter
        }
      }
    }
  }
  
  plots$scatterplots <- scatterplots
  cat(sprintf("Created %d scatterplots\n", length(scatterplots)))

  # Per-study inter-instrument correlation heatmaps (>= 3 unique instruments only)
  if (!is.null(inter_instrument_results) &&
      nrow(inter_instrument_results$pair_correlations) > 0) {

    plots$inter_instrument_heatmaps <- list()

    pair_df <- inter_instrument_results$pair_correlations
    studies_with_pairs <- unique(pair_df$study)

    for (current_study in studies_with_pairs) {
      study_pairs <- pair_df %>% filter(study == current_study)

      # Only plot when >= 3 unique instruments (>= 3 pairs)
      if (nrow(study_pairs) < 3) next

      instruments_in_study <- unique(c(study_pairs$instrument1, study_pairs$instrument2))
      n_inst <- length(instruments_in_study)
      mat <- matrix(1, n_inst, n_inst,
                    dimnames = list(instruments_in_study, instruments_in_study))
      for (k in seq_len(nrow(study_pairs))) {
        i1 <- study_pairs$instrument1[k]
        i2 <- study_pairs$instrument2[k]
        mat[i1, i2] <- study_pairs$correlation[k]
        mat[i2, i1] <- study_pairs$correlation[k]
      }

      tryCatch({
        corrplot(mat,
                 method       = "color",
                 type         = "upper",
                 is.corr      = TRUE,
                 addCoef.col  = "black",
                 number.cex   = 0.8,
                 tl.col       = "black",
                 tl.srt       = 45,
                 title        = paste0("Inter-instrument correlations: ", current_study),
                 mar          = c(0, 0, 2, 0))
        p <- recordPlot()
        plots$inter_instrument_heatmaps[[current_study]] <- p
      }, error = function(e) {
        cat(sprintf("  Could not create heatmap for %s: %s\n", current_study, e$message))
      })
    }
    cat(sprintf("Created %d inter-instrument heatmaps\n",
                length(plots$inter_instrument_heatmaps)))

    # Summary dot plot across all studies
    if (nrow(inter_instrument_results$pair_summary) > 0) {
      p_pair_summary <- ggplot(
        inter_instrument_results$pair_summary,
        aes(x = mean_r,
            y = reorder(paste(instrument1, "vs", instrument2), mean_r),
            size = n_studies)) +
        geom_point(color = "steelblue", alpha = 0.8) +
        geom_errorbarh(aes(xmin = min_r, xmax = max_r), height = 0.2, alpha = 0.5) +
        geom_vline(xintercept = c(0.50, 0.70), linetype = "dashed",
                   color = "red", alpha = 0.5) +
        labs(title    = "Inter-Instrument Correlations Across Studies",
             subtitle = "Points = mean r; bars = range across studies; dashed lines at r = 0.50 and 0.70",
             x        = "Pearson r (POMP scores)",
             y        = "Instrument Pair",
             size     = "N studies") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"))

      plots$inter_instrument_summary <- p_pair_summary
      cat("Created inter-instrument correlation summary plot\n")
    }
  }

  # Back-translation scatter plots and pooled summary bar chart
  if (!is.null(backtranslation_results) &&
      !is.null(backtranslation_results$person_level) &&
      nrow(backtranslation_results$person_level) > 0) {

    plots$backtranslation <- list()

    for (inst in unique(backtranslation_results$person_level$instrument)) {
      inst_data    <- backtranslation_results$person_level %>%
        filter(instrument == inst)
      inst_metrics <- backtranslation_results$pooled_metrics %>%
        filter(instrument == inst)

      p_bt <- ggplot(inst_data,
                     aes(x = expected_sum, y = observed_sum, colour = study)) +
        geom_point(alpha = 0.3, size = 1) +
        geom_abline(intercept = 0, slope = 1,
                    linetype = "dashed", colour = "black", linewidth = 0.8) +
        geom_smooth(aes(group = 1), method = "lm",
                    colour = "red", se = TRUE, linewidth = 0.8) +
        labs(
          title    = sprintf("Back-Translation: %s", inst),
          subtitle = sprintf(
            "r = %.3f  |  RMSE = %.3f  |  MAD = %.3f  |  Bias = %.3f  |  n = %d",
            inst_metrics$r, inst_metrics$rmse,
            inst_metrics$mad, inst_metrics$bias, inst_metrics$n),
          x      = "Expected Sum Score (from IRT theta)",
          y      = "Observed Sum Score",
          colour = "Study"
        ) +
        theme_minimal() +
        theme(plot.title      = element_text(face = "bold"),
              plot.subtitle   = element_text(size = 9),
              legend.position = "right")

      plots$backtranslation[[inst]] <- p_bt
    }

    # Summary bar chart: r, RMSE, MAD side by side per instrument
    summary_long <- backtranslation_results$pooled_metrics %>%
      select(instrument, r, rmse, mad) %>%
      pivot_longer(cols      = c(r, rmse, mad),
                   names_to  = "metric",
                   values_to = "value") %>%
      mutate(metric = factor(metric,
                             levels = c("r", "rmse", "mad"),
                             labels = c("Correlation (r)", "RMSE", "MAD")))

    p_bt_summary <- ggplot(summary_long,
                            aes(x = instrument, y = value, fill = metric)) +
      geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
      geom_hline(data = data.frame(
                   metric     = factor("Correlation (r)",
                                       levels = levels(summary_long$metric)),
                   yintercept = 0.90),
                 aes(yintercept = yintercept),
                 linetype = "dashed", colour = "red") +
      facet_wrap(~ metric, scales = "free_y") +
      labs(
        title    = "Back-Translation Metrics by Instrument (Pooled)",
        subtitle = "Dashed red line on correlation panel = r 0.90 threshold",
        x        = "Instrument",
        y        = "Value",
        fill     = "Metric"
      ) +
      theme_minimal() +
      theme(plot.title      = element_text(face = "bold"),
            legend.position = "none")

    plots$backtranslation_summary <- p_bt_summary
    cat(sprintf("Created %d back-translation scatter plots\n",
                length(plots$backtranslation)))
    cat("Created back-translation summary bar chart\n")
  }

  return(plots)
}

# =============================================================================
# STEP 6: EXPORT RESULTS
# =============================================================================

export_validation_results <- function(sum_scores, correlation_results, plots,
                                      output_dir = "validation_results",
                                      inter_instrument_results = NULL,
                                      backtranslation_results = NULL) {
  
  cat("\nSTEP 6: EXPORTING RESULTS\n")
  cat("=========================\n\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Export sum scores
  write.csv(sum_scores, 
            file.path(output_dir, "instrument_sum_scores.csv"),
            row.names = FALSE)
  cat("✓ Exported instrument_sum_scores.csv\n")
  
  # Export correlations
  write.csv(correlation_results,
            file.path(output_dir, "validation_correlations.csv"),
            row.names = FALSE)
  cat("✓ Exported validation_correlations.csv\n")
  
  # Export plots
  if ("heatmap" %in% names(plots)) {
    png(file.path(output_dir, "correlation_heatmap.png"), 
        width = 10, height = 8, units = "in", res = 300)
    plots$heatmap
    dev.off()
    cat("✓ Exported correlation_heatmap.png\n")
  }
  
  if ("distribution" %in% names(plots)) {
    ggsave(file.path(output_dir, "correlation_distribution.png"),
           plot = plots$distribution, width = 8, height = 6)
    cat("✓ Exported correlation_distribution.png\n")
  }
  
  if ("by_instrument" %in% names(plots)) {
    ggsave(file.path(output_dir, "correlations_by_instrument.png"),
           plot = plots$by_instrument, width = 8, height = 6)
    cat("✓ Exported correlations_by_instrument.png\n")
  }
  
  # Export scatterplots
  if ("scatterplots" %in% names(plots) && length(plots$scatterplots) > 0) {
    
    scatterplot_dir <- file.path(output_dir, "scatterplots")
    if (!dir.exists(scatterplot_dir)) {
      dir.create(scatterplot_dir)
    }
    
    for (plot_name in names(plots$scatterplots)) {
      ggsave(file.path(scatterplot_dir, paste0(plot_name, ".png")),
             plot = plots$scatterplots[[plot_name]], 
             width = 6, height = 5)
    }
    
    cat(sprintf("✓ Exported %d scatterplots to scatterplots/\n", 
                length(plots$scatterplots)))
  }
  
  # Export inter-instrument correlation results
  if (!is.null(inter_instrument_results)) {

    write.csv(inter_instrument_results$pair_correlations,
              file.path(output_dir, "inter_instrument_pair_correlations.csv"),
              row.names = FALSE)
    cat("✓ Exported inter_instrument_pair_correlations.csv\n")

    write.csv(inter_instrument_results$pair_summary,
              file.path(output_dir, "inter_instrument_pair_summary.csv"),
              row.names = FALSE)
    cat("✓ Exported inter_instrument_pair_summary.csv\n")

    if ("inter_instrument_summary" %in% names(plots) &&
        !is.null(plots$inter_instrument_summary)) {
      ggsave(file.path(output_dir, "inter_instrument_summary.png"),
             plot = plots$inter_instrument_summary, width = 10, height = 8)
      cat("✓ Exported inter_instrument_summary.png\n")
    }

    if ("inter_instrument_heatmaps" %in% names(plots) &&
        length(plots$inter_instrument_heatmaps) > 0) {
      heatmap_dir <- file.path(output_dir, "inter_instrument_heatmaps")
      if (!dir.exists(heatmap_dir)) dir.create(heatmap_dir)
      for (study_name in names(plots$inter_instrument_heatmaps)) {
        tryCatch({
          png(file.path(heatmap_dir, paste0(study_name, "_inter_instrument.png")),
              width = 800, height = 700)
          replayPlot(plots$inter_instrument_heatmaps[[study_name]])
          dev.off()
        }, error = function(e) {
          cat(sprintf("  Could not save heatmap for %s: %s\n",
                      study_name, e$message))
        })
      }
      cat(sprintf("✓ Exported %d inter-instrument heatmaps\n",
                  length(plots$inter_instrument_heatmaps)))
    }
  }

  # Export back-translation results
  if (!is.null(backtranslation_results)) {

    write.csv(backtranslation_results$pooled_metrics,
              file.path(output_dir, "backtranslation_pooled_metrics.csv"),
              row.names = FALSE)
    cat("✓ Exported backtranslation_pooled_metrics.csv\n")

    write.csv(backtranslation_results$study_metrics,
              file.path(output_dir, "backtranslation_study_metrics.csv"),
              row.names = FALSE)
    cat("✓ Exported backtranslation_study_metrics.csv\n")

    write.csv(backtranslation_results$person_level,
              file.path(output_dir, "backtranslation_person_level.csv"),
              row.names = FALSE)
    cat("✓ Exported backtranslation_person_level.csv\n")

    if ("backtranslation_summary" %in% names(plots) &&
        !is.null(plots$backtranslation_summary)) {
      ggsave(file.path(output_dir, "backtranslation_summary.png"),
             plot  = plots$backtranslation_summary,
             width = 10, height = 6)
      cat("✓ Exported backtranslation_summary.png\n")
    }

    if ("backtranslation" %in% names(plots) &&
        length(plots$backtranslation) > 0) {
      bt_dir <- file.path(output_dir, "backtranslation_plots")
      if (!dir.exists(bt_dir)) dir.create(bt_dir)
      for (inst_name in names(plots$backtranslation)) {
        tryCatch(
          ggsave(file.path(bt_dir, paste0(inst_name, "_backtranslation.png")),
                 plot  = plots$backtranslation[[inst_name]],
                 width = 8, height = 6),
          error = function(e)
            cat(sprintf("  Could not save plot for %s: %s\n",
                        inst_name, e$message))
        )
      }
      cat(sprintf("✓ Exported %d back-translation scatter plots\n",
                  length(plots$backtranslation)))
    }
  }

  cat(sprintf("\nAll results saved to: %s\n", output_dir))
}

# =============================================================================
# MASTER VALIDATION FUNCTION
# =============================================================================

validate_harmonized_scores <- function(data, harmonized_scores,
                                       study_var = "Study",
                                       theta_var = "Theta_Standardized",
                                       min_items_for_score = 0.5,
                                       output_dir = "validation_results",
                                       model = NULL,
                                       model_items = NULL) {
  
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║  HARMONIZED SCORE VALIDATION: COMPARISON TO ORIGINAL SCALES  ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n\n")
  
  # Step 1: Identify instruments
  study_instruments <- identify_study_instruments(data, study_var)
  
  # Step 2: Calculate sum scores
  sum_scores <- calculate_instrument_sum_scores(data, study_instruments, 
                                                study_var, min_items_for_score)
  
  # Step 3: Merge with harmonized scores
  merged_data <- merge_harmonized_and_sum_scores(harmonized_scores, sum_scores)
  
  # Step 4: Examine correlations
  correlation_results <- examine_correlations_by_study(merged_data, study_instruments,
                                                       study_var, theta_var)

  # Step 4b: Inter-instrument correlations within studies
  inter_instrument_results <- examine_inter_instrument_correlations(
    sum_scores        = sum_scores,
    study_instruments = study_instruments,
    study_var         = study_var,
    min_n             = 30
  )

  # Step 4c: Back-translation metrics (observed vs IRT-expected sum scores)
  backtranslation_results <- NULL
  if (!is.null(model) && !is.null(model_items)) {
    cat("\nSTEP 4c: BACK-TRANSLATION VALIDATION\n")
    cat("=====================================\n")
    backtranslation_results <- compute_backtranslation_metrics(
      data               = merged_data,
      raw_data           = data,
      model              = model,
      model_items        = model_items,
      study_instruments  = study_instruments,
      study_var          = study_var,
      theta_var          = theta_var,
      target_instruments = c("K10", "SF36", "GHQ12", "SDQ"),
      min_items_in_model = 3,
      min_n_instrument   = 30
    )
  } else {
    cat("\nSkipping back-translation validation (model not supplied)\n")
  }

  # Step 5: Create visualizations
  plots <- create_validation_plots(merged_data, correlation_results,
                                   study_instruments, study_var, theta_var,
                                   inter_instrument_results,
                                   backtranslation_results)

  # Step 6: Export results
  export_validation_results(sum_scores, correlation_results, plots,
                            output_dir, inter_instrument_results,
                            backtranslation_results)

  cat("\n╔══════════════════════════════════════════════════════════════╗\n")
  cat("║  VALIDATION COMPLETE                                         ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")

  return(list(
    study_instruments       = study_instruments,
    sum_scores              = sum_scores,
    merged_data             = merged_data,
    correlations            = correlation_results,
    inter_instrument        = inter_instrument_results,
    backtranslation         = backtranslation_results,
    plots                   = plots
  ))
}

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# After running your harmonization:
# results <- harmonize_studies_reference_wave_concurrent(...)
# 
# # Get harmonized scores
# harmonized_scores <- results$final_scores

har_score_ref <- harmonized_scores %>%
  filter(Wave==1)

# # Run validation
 validation_results <- validate_harmonized_scores(
   data                = merged_data,  # Your original data
   harmonized_scores   = har_score_ref,
   study_var           = "Study",
   theta_var           = "Theta_Score",
   min_items_for_score = 0.5,  # Require 50% of items for valid sum score
   output_dir          = "validation_results_ref",
   model               = results$concurrent_calibration$model,
   model_items         = results$concurrent_calibration$final_items
 )

# --- Wave-invariance validation ---
# Subsequent wave scores come from fixed-parameter scoring off the Wave 1 model,
# so their validity cannot be assumed from Wave 1 results alone.
har_score_subseq <- harmonized_scores %>%
  filter(Wave > 1)

if (nrow(har_score_subseq) > 0) {

  validation_results_subseq <- validate_harmonized_scores(
    data = merged_data,
    harmonized_scores = har_score_subseq,
    study_var = "Study",
    theta_var = "Theta_Score",
    min_items_for_score = 0.5,
    output_dir = "validation_results_subsequent"
  )

  # Extract mean r by instrument (sum scores) from each validation run
  inst_summary_ref <- validation_results$correlations %>%
    filter(score_type == "sum") %>%
    group_by(instrument) %>%
    summarise(mean_r_wave1 = mean(correlation, na.rm = TRUE), .groups = "drop")

  inst_summary_subseq <- validation_results_subseq$correlations %>%
    filter(score_type == "sum") %>%
    group_by(instrument) %>%
    summarise(mean_r_subsequent = mean(correlation, na.rm = TRUE), .groups = "drop")

  wave_invariance_comparison <- inst_summary_ref %>%
    inner_join(inst_summary_subseq, by = "instrument") %>%
    mutate(r_diff = mean_r_subsequent - mean_r_wave1)

  # Side-by-side comparison: Wave 1 vs subsequent waves mean r by instrument
  cat("\n--- WAVE-INVARIANCE VALIDATION: Wave 1 vs Subsequent Waves ---\n")
  print(wave_invariance_comparison)
  cat("\n")

  # Flag instruments where subsequent wave mean r drops more than 0.05 below Wave 1
  flagged <- wave_invariance_comparison %>% filter(r_diff < -0.05)
  if (nrow(flagged) > 0) {
    for (i in seq_len(nrow(flagged))) {
      warning(sprintf(
        "Wave-invariance concern: %s subsequent mean r (%.3f) is >0.05 below Wave 1 mean r (%.3f)",
        flagged$instrument[i], flagged$mean_r_subsequent[i], flagged$mean_r_wave1[i]
      ))
    }
  }

} else {
  cat("No subsequent wave scores found; skipping subsequent wave validation.\n")
}

#
# # Examine results
# summary(validation_results$correlations)
# 
# # View specific scatterplot
# print(validation_results$plots$scatterplots$LSAC_SDQ)

# =============================================================================
# INTERPRETATION GUIDE
# =============================================================================

# CORRELATION STRENGTH INTERPRETATION:
# - r > 0.90: Excellent - very strong agreement with original measure
# - r > 0.80: Good - strong agreement, minor differences
# - r > 0.70: Acceptable - adequate agreement, some divergence
# - r > 0.60: Moderate - noticeable differences, investigate
# - r < 0.60: Weak - substantial differences, major concern
#
# WHAT TO LOOK FOR:
# 1. Are correlations generally high (>0.80) across studies?
# 2. Are there specific instruments with consistently lower correlations?
# 3. Are there specific studies with problematic correlations?
# 4. Do POMP scores (normalized) correlate better than raw sum scores?
#
# POTENTIAL ISSUES:
# - Low correlations may indicate:
#   a) IRT model not capturing the construct well
#   b) DIF problems affecting harmonization
#   c) Original sum scores unreliable (many missing items)
#   d) Different constructs being measured
#
# RECOMMENDATIONS:
# - Overall mean r > 0.80: Harmonization successful
# - Overall mean r 0.70-0.80: Generally acceptable, note limitations
# - Overall mean r < 0.70: Investigate causes, consider refinements