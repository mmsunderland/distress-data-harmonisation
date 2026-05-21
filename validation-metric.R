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
# STEP 5: CREATE VISUALIZATIONS
# =============================================================================

create_validation_plots <- function(merged_data, correlation_results,
                                    study_instruments, study_var = "Study",
                                    theta_var = "Theta_Standardized") {
  
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
  
  return(plots)
}

# =============================================================================
# STEP 6: EXPORT RESULTS
# =============================================================================

export_validation_results <- function(sum_scores, correlation_results, plots,
                                      output_dir = "validation_results") {
  
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
  
  cat(sprintf("\nAll results saved to: %s\n", output_dir))
}

# =============================================================================
# MASTER VALIDATION FUNCTION
# =============================================================================

validate_harmonized_scores <- function(data, harmonized_scores,
                                       study_var = "Study",
                                       theta_var = "Theta_Standardized",
                                       min_items_for_score = 0.5,
                                       output_dir = "validation_results") {
  
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
  
  # Step 5: Create visualizations
  plots <- create_validation_plots(merged_data, correlation_results, 
                                   study_instruments, study_var, theta_var)
  
  # Step 6: Export results
  export_validation_results(sum_scores, correlation_results, plots, output_dir)
  
  cat("\n╔══════════════════════════════════════════════════════════════╗\n")
  cat("║  VALIDATION COMPLETE                                         ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  
  return(list(
    study_instruments = study_instruments,
    sum_scores = sum_scores,
    merged_data = merged_data,
    correlations = correlation_results,
    plots = plots
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
   data = data_dif,  # Your original data
   harmonized_scores = har_score_ref,
   study_var = "Study",
   theta_var = "Theta_Score",
   min_items_for_score = 0.5,  # Require 50% of items for valid sum score
   output_dir = "validation_results_ref"
 )
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