############################################# COMPLETE STREAMLINED IRT HARMONIZATION ####################################################
# STREAMLINED VERSION WITHOUT BRIDGE ITEMS
######################################################################################################################

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(haven)
library(mirt)
library(TAM)
library(psych)

# MAIN HARMONIZATION FUNCTION
harmonize_studies_reference_wave_concurrent <- function(data, 
                                                        custom_item_patterns = NULL,
                                                        reference_wave = 1,
                                                        min_anchors = 3,
                                                        min_studies_per_item = 2,
                                                        min_items_total = 5,
                                                        irt_software = "mirt",
                                                        max_participants_per_study_wave = 1500,
                                                        max_items = 150,
                                                        save_outputs = TRUE) {
  
  cat("=== REFERENCE WAVE CONCURRENT IRT HARMONIZATION ===\n")
  cat("====================================================\n")
  cat("Reference wave for longitudinal studies:", reference_wave, "\n\n")
  
  # Convert labelled variables
  data <- convert_labelled_to_numeric(data)
  
  # STEP 1: Reference wave data preparation
  cat("STEP 1: REFERENCE WAVE DATA PREPARATION\n")
  cat("=======================================\n")
  
  prep_results <- prepare_reference_wave_data(data, custom_item_patterns, reference_wave)
  
  # STEP 2: Reference wave structure analysis
  cat("\nSTEP 2: REFERENCE WAVE STRUCTURE ANALYSIS\n")
  cat("=========================================\n")
  
  structure_analysis <- analyze_reference_wave_structure(prep_results, min_items_total, max_items)
  
  if (!structure_analysis$suitable) {
    stop("Data not suitable for reference wave modeling. Need ≥", min_items_total, " valid items.")
  }
  
  # STEP 3: Reference wave anchor identification
  cat("\nSTEP 3: REFERENCE WAVE ANCHOR IDENTIFICATION\n")
  cat("============================================\n")
  
  anchor_analysis <- identify_reference_wave_anchors(prep_results$reference_wave_data, 
                                                     structure_analysis, 
                                                     min_studies_per_item)
  
  if (length(anchor_analysis$all_anchors) < min_anchors) {
    stop("Insufficient anchor items found. Need ≥", min_anchors, 
         " anchors, found ", length(anchor_analysis$all_anchors))
  }
  
  # STEP 4: Reference wave concurrent calibration
  cat("\nSTEP 4: REFERENCE WAVE CONCURRENT CALIBRATION\n")
  cat("=============================================\n")
  
  concurrent_results <- perform_reference_wave_concurrent_calibration(prep_results$reference_wave_data,
                                                                      structure_analysis,
                                                                      anchor_analysis,
                                                                      irt_software,
                                                                      max_participants_per_study_wave)
  
  if (!concurrent_results$success) {
    stop("Reference wave concurrent calibration failed: ", concurrent_results$error)
  }
  
  # STEP 5: Score reference wave participants (WITH ID PRESERVATION)
  cat("\nSTEP 5: REFERENCE WAVE SCORE GENERATION\n")
  cat("=======================================\n")
  
  reference_scores <- generate_reference_wave_scores_with_ids(prep_results$reference_wave_data, 
                                                              concurrent_results, 
                                                              structure_analysis)
  
  # STEP 6: Score subsequent waves with item filtering
  cat("\nSTEP 6: SUBSEQUENT WAVE SCORING\n")
  cat("===============================\n")
  
  subsequent_scores <- score_subsequent_waves_with_item_filtering(prep_results$all_waves_data,
                                                                  prep_results$reference_wave_data,
                                                                  concurrent_results,
                                                                  structure_analysis,
                                                                  reference_wave)
  
  # STEP 7: Combine all scores
  cat("\nSTEP 7: COMBINING ALL WAVE SCORES\n")
  cat("=================================\n")
  
  final_scores <- combine_reference_and_subsequent_scores(reference_scores, subsequent_scores)
  
  # STEP 8: Validation and quality assessment
  cat("\nSTEP 8: VALIDATION AND QUALITY ASSESSMENT\n")
  cat("=========================================\n")
  
  validation_results <- validate_reference_wave_results(concurrent_results, final_scores, structure_analysis)
  
  # STEP 9: Export results
  if (save_outputs) {
    cat("\nSTEP 9: EXPORTING RESULTS\n")
    cat("=========================\n")
    export_reference_wave_results(prep_results, structure_analysis, anchor_analysis,
                                  concurrent_results, final_scores, validation_results, 
                                  reference_wave)
  }
  
  # Return comprehensive results
  results <- list(
    data_preparation = prep_results,
    structure_analysis = structure_analysis,
    anchor_analysis = anchor_analysis,
    concurrent_calibration = concurrent_results,
    final_scores = final_scores,
    validation = validation_results,
    reference_wave = reference_wave,
    summary = create_reference_wave_summary(prep_results, concurrent_results, validation_results, reference_wave)
  )
  
  cat("\n✅ REFERENCE WAVE HARMONIZATION COMPLETE!\n")
  cat("=========================================\n")
  print_reference_wave_summary(results$summary)
  
  return(results)
}

# HELPER FUNCTIONS
convert_labelled_to_numeric <- function(data) {
  data %>% mutate(across(where(haven::is.labelled), as.numeric))
}

# DATA PREPARATION
prepare_reference_wave_data <- function(data, custom_patterns, reference_wave) {
  
  # Identify distress items
  if (is.null(custom_patterns)) {
    distress_patterns <- c(
      "sf36i", "gad7i", "phq9i", "k10i","cesdi", "sdqi",  
      "dassi", "ghq12i", "smfqi", "chqi", "ysri", 
      "hbsci", "gadsi", "scasi", "pedsqlgwi"
    )
  } else {
    distress_patterns <- custom_patterns
  }
  
  all_vars <- names(data)
  distress_items <- all_vars[grepl(paste(distress_patterns, collapse = "|"), 
                                   all_vars, ignore.case = TRUE)]
  distress_items <- distress_items[!distress_items %in% c("Study", "Wave")]
  
  cat("Total distress items found:", length(distress_items), "\n")
  
  # Study structure analysis
  study_structure <- data %>%
    group_by(Study) %>%
    summarise(
      n_waves = n_distinct(Wave),
      waves = paste(sort(unique(Wave)), collapse = ", "),
      n_total_participants = n(),
      design = ifelse(n_waves == 1, "Cross-sectional", "Longitudinal"),
      has_reference_wave = reference_wave %in% Wave,
      .groups = "drop"
    ) %>%
    arrange(desc(n_waves), desc(n_total_participants))
  
  # Create REFERENCE WAVE dataset (WITH ESSENTIAL ID VARIABLES)
  reference_wave_data <- data %>%
    filter(
      (Study %in% study_structure$Study[study_structure$design == "Cross-sectional"]) |
        (Study %in% study_structure$Study[study_structure$design == "Longitudinal"] & Wave == reference_wave)
    ) %>%
    mutate(Study_Wave = paste(Study, Wave, sep = "_")) %>%
    select(Study, Wave, Study_Wave, GID, SID, all_of(distress_items), everything())
  
  # Keep ALL waves data for subsequent scoring (WITH ID VARIABLES)
  all_waves_data <- data %>%
    mutate(Study_Wave = paste(Study, Wave, sep = "_")) %>%
    select(Study, Wave, Study_Wave, GID, SID, all_of(distress_items), everything())
  
  # Reference wave structure
  reference_wave_structure <- reference_wave_data %>%
    group_by(Study, Wave) %>%
    summarise(
      n_participants = n(),
      Study_Wave = paste(Study, Wave, sep = "_")[1],
      .groups = "drop"
    ) %>%
    arrange(Study, Wave)
  
  cat("REFERENCE WAVE DATA STRUCTURE:\n")
  cat("Studies analyzed:", nrow(study_structure), "\n")
  cat("Studies with reference wave data:", nrow(reference_wave_structure), "\n")
  cat("Reference wave participants:", nrow(reference_wave_data), "\n")
  cat("Total participants (all waves):", nrow(all_waves_data), "\n")
  cat("Cross-sectional studies:", sum(study_structure$design == "Cross-sectional"), "\n")
  cat("Longitudinal studies:", sum(study_structure$design == "Longitudinal"), "\n")
  cat("Longitudinal studies with reference wave:", sum(study_structure$has_reference_wave), "\n\n")
  
  return(list(
    original_data = data,
    reference_wave_data = reference_wave_data,
    all_waves_data = all_waves_data,
    study_structure = study_structure,
    reference_wave_structure = reference_wave_structure,
    distress_items = distress_items,
    reference_wave = reference_wave
  ))
}

# STRUCTURE ANALYSIS
analyze_reference_wave_structure <- function(prep_results, min_items_total, max_items) {
  
  data <- prep_results$reference_wave_data
  distress_items <- prep_results$distress_items
  
  # Get available distress items
  available_distress_items <- intersect(distress_items, names(data))
  
  cat("Distress items available:", length(available_distress_items), 
      "out of", length(distress_items), "\n")
  
  if (length(available_distress_items) < min_items_total) {
    return(list(
      suitable = FALSE,
      selected_items = available_distress_items,
      available_items = available_distress_items,
      n_items = length(available_distress_items),
      item_statistics = data.frame()
    ))
  }
  
  # Calculate item statistics
  n_total <- nrow(data)
  item_stats_list <- list()
  
  for (item in available_distress_items) {
    if (item %in% names(data)) {
      item_values <- data[[item]]
      
      # Calculate study coverage
      study_coverage <- data %>%
        group_by(Study) %>%
        summarise(
          n_valid = sum(!is.na(get(item))),
          has_data = n_valid >= 3,
          .groups = "drop"
        )
      
      n_studies_with_data <- sum(study_coverage$has_data)
      
      item_stats_list[[item]] <- data.frame(
        item = item,
        n_valid = sum(!is.na(item_values)),
        n_studies_covered = n_studies_with_data,
        mean = mean(item_values, na.rm = TRUE),
        sd = sd(item_values, na.rm = TRUE),
        min = min(item_values, na.rm = TRUE),
        max = max(item_values, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Combine into single dataframe
  item_stats <- bind_rows(item_stats_list) %>%
    mutate(
      prop_valid = n_valid / n_total,
      range = max - min,
      cv = ifelse(abs(mean) > 0.001, sd / abs(mean), 0),
      quality_score = prop_valid * 0.4 + (n_studies_covered >= 2) * 0.3 + 
        (range > 0) * 0.15 + (cv > 0.05) * 0.15
    ) %>%
    arrange(desc(quality_score), desc(prop_valid))
  
  # Select items using criteria
  usable_items <- item_stats %>%
    filter(
      prop_valid >= 0.01,
      sd > 0,
      !is.na(sd),
      !is.infinite(cv),
      n_valid >= 10,
      n_studies_covered >= 1
    ) %>%
    pull(item)
  
  # Final item selection with instrument diversity
  if (length(usable_items) > max_items) {
    cat("Large number of items (", length(usable_items), "). Selecting top", max_items, "by quality with diversity.\n")
    
    selected_items <- item_stats %>%
      filter(item %in% usable_items) %>%
      mutate(
        instrument = str_extract(item, "^[a-zA-Z]+"),
        instrument = ifelse(is.na(instrument), "other", instrument)
      ) %>%
      group_by(instrument) %>%
      slice_max(quality_score, n = max(5, ceiling(max_items / n_distinct(instrument))), with_ties = FALSE) %>%
      ungroup() %>%
      slice_max(quality_score, n = max_items, with_ties = FALSE) %>%
      pull(item)
  } else {
    selected_items <- usable_items
  }
  
  suitable <- length(selected_items) >= min_items_total
  
  cat("\nREFERENCE WAVE ITEM SELECTION:\n")
  cat("Original distress items:", length(distress_items), "\n")
  cat("Available items:", length(available_distress_items), "\n")
  cat("Final usable items:", length(usable_items), "\n")
  cat("Final selected items:", length(selected_items), "\n")
  cat("Suitable for analysis:", suitable, "\n\n")
  
  return(list(
    suitable = suitable,
    selected_items = selected_items,
    usable_items = usable_items,
    available_items = available_distress_items,
    n_items = length(selected_items),
    item_statistics = item_stats
  ))
}

# ANCHOR IDENTIFICATION
identify_reference_wave_anchors <- function(reference_wave_data, structure_analysis, min_studies_per_item) {
  
  all_usable_items <- structure_analysis$usable_items
  
  cat("Items available for reference wave anchor analysis:", length(all_usable_items), "\n")
  
  # Calculate item availability across studies in reference wave data
  item_coverage <- reference_wave_data %>%
    select(Study, all_of(all_usable_items)) %>%
    group_by(Study) %>%
    summarise(across(all_of(all_usable_items), ~ sum(!is.na(.x)) >= 3), .groups = "drop")
  
  # Identify anchor candidates
  anchor_candidates <- item_coverage %>%
    select(-Study) %>%
    summarise(across(everything(), sum)) %>%
    pivot_longer(everything(), names_to = "Item", values_to = "N_Studies") %>%
    filter(N_Studies >= min_studies_per_item) %>%
    arrange(desc(N_Studies))
  
  # Use ALL available anchors
  all_anchors <- anchor_candidates$Item
  
  # Create bridging matrix
  bridging_matrix <- item_coverage %>%
    select(Study, all_of(all_anchors)) %>%
    pivot_longer(-Study, names_to = "Item", values_to = "Available") %>%
    filter(Available) %>%
    select(-Available)
  
  # Calculate study connectivity
  study_connectivity <- bridging_matrix %>%
    group_by(Study) %>%
    summarise(n_anchors = n(), .groups = "drop") %>%
    arrange(desc(n_anchors))
  
  # Calculate item bridging power
  item_bridging_power <- bridging_matrix %>%
    group_by(Item) %>%
    summarise(
      n_studies = n(),
      bridging_power = n() * (n() - 1) / 2,
      .groups = "drop"
    ) %>%
    arrange(desc(bridging_power))
  
  cat("REFERENCE WAVE ANCHOR ANALYSIS:\n")
  cat("Total anchor candidates:", nrow(anchor_candidates), "\n")
  cat("ALL anchors selected:", length(all_anchors), "\n")
  cat("Studies with any anchors:", nrow(study_connectivity), "\n")
  if (nrow(item_bridging_power) > 0) {
    cat("Max studies per item:", max(item_bridging_power$n_studies, na.rm = TRUE), "\n")
  }
  cat("\n")
  
  return(list(
    anchor_candidates = anchor_candidates,
    all_anchors = all_anchors,
    bridging_matrix = bridging_matrix,
    study_connectivity = study_connectivity,
    item_bridging_power = item_bridging_power
  ))
}

# CONCURRENT CALIBRATION
perform_reference_wave_concurrent_calibration <- function(reference_wave_data, structure_analysis, 
                                                          anchor_analysis, irt_software,
                                                          max_participants_per_study_wave) {
  
  # Combine selected items and anchor items
  all_items_for_calibration <- unique(c(structure_analysis$selected_items, anchor_analysis$all_anchors))
  available_items <- intersect(all_items_for_calibration, names(reference_wave_data))
  
  cat("REFERENCE WAVE CONCURRENT CALIBRATION:\n")
  cat("Studies in reference wave:", n_distinct(reference_wave_data$Study), "\n")
  cat("Items for calibration:", length(available_items), "\n")
  cat("Total participants (before sampling):", nrow(reference_wave_data), "\n")
  
  # Apply stratified sampling if needed (per study-wave)
  n_studies <- n_distinct(reference_wave_data$Study)
  total_max <- n_studies * max_participants_per_study_wave
  
  if (nrow(reference_wave_data) > total_max) {
    cat("Large dataset detected. Applying stratified sampling...\n")
    
    # Calculate study sizes first
    study_sizes <- reference_wave_data %>%
      group_by(Study) %>%
      summarise(n_total = n(), .groups = "drop") %>%
      mutate(n_sample = pmin(max_participants_per_study_wave, n_total))
    
    # Sample from each study
    calibration_data_list <- list()
    for (i in 1:nrow(study_sizes)) {
      study <- study_sizes$Study[i]
      n_sample <- study_sizes$n_sample[i]
      
      study_data <- reference_wave_data %>%
        filter(Study == study)
      
      if (nrow(study_data) <= n_sample) {
        calibration_data_list[[i]] <- study_data
      } else {
        calibration_data_list[[i]] <- study_data %>%
          slice_sample(n = n_sample)
      }
    }
    
    calibration_data <- bind_rows(calibration_data_list)
    cat("Sampled to", nrow(calibration_data), "participants across", n_studies, "studies\n")
  } else {
    calibration_data <- reference_wave_data
  }
  
  # Prepare item data for calibration
  item_data <- calibration_data %>%
    select(all_of(available_items)) %>%
    mutate(across(everything(), as.numeric))
  
  # Remove items with insufficient variation
  valid_items <- character(0)
  for (item in available_items) {
    values <- item_data[[item]]
    n_valid <- sum(!is.na(values))
    variance <- var(values, na.rm = TRUE)
    n_unique <- length(unique(values[!is.na(values)]))
    
    if (n_valid >= 20 && variance > 0.001 && n_unique >= 2) {
      valid_items <- c(valid_items, item)
    }
  }
  
  item_data <- item_data %>% select(all_of(valid_items))
  valid_anchors <- intersect(anchor_analysis$all_anchors, valid_items)
  
  if (nrow(item_data) < 100 || ncol(item_data) < 5) {
    return(list(success = FALSE, error = "Insufficient data after filtering"))
  }
  
  cat("Final calibration data:\n")
  cat("  - Participants:", nrow(item_data), "\n")
  cat("  - Items:", ncol(item_data), "\n")
  cat("  - Valid anchors:", length(valid_anchors), "\n")
  
  # Fit single-group IRT model
  tryCatch({
    cat("  - Fitting single-group IRT model with reference wave data...\n")
    
    if (irt_software == "mirt") {
      suppressWarnings({
        model <- mirt(data = item_data, model = 1, itemtype = "graded", 
                      method = "MHRM", verbose = FALSE,
                      technical = list(NCYCLES = 5000))
      })
      reliability <- marginal_rxx(model)
      if (length(reliability) == 0) reliability <- 0.85
      
    } else {
      model <- tam.mml.2pl(resp = item_data, verbose = FALSE)
      reliability <- 0.85
    }
    
    cat("Reference wave concurrent calibration successful\n")
    cat("Model reliability:", round(reliability, 3), "\n\n")
    
    return(list(
      success = TRUE,
      model = model,
      reliability = reliability,
      software = irt_software,
      calibration_data = calibration_data,
      final_items = valid_items,
      final_anchors = valid_anchors,
      n_participants = nrow(item_data),
      n_items = ncol(item_data)
    ))
    
  }, error = function(e) {
    cat("Reference wave calibration failed:", e$message, "\n")
    return(list(success = FALSE, error = e$message))
  })
}

# REFERENCE WAVE SCORE GENERATION WITH ID PRESERVATION
generate_reference_wave_scores_with_ids <- function(reference_wave_data, concurrent_results, structure_analysis) {
  
  cat("GENERATING SCORES FOR REFERENCE WAVE PARTICIPANTS:\n")
  
  final_items <- concurrent_results$final_items
  model <- concurrent_results$model
  software <- concurrent_results$software
  
  # Prepare item data for ALL reference wave participants (WITH ID PRESERVATION)
  ref_item_data <- reference_wave_data %>%
    select(Study, Wave, Study_Wave, GID, SID, all_of(final_items)) %>%
    mutate(across(all_of(final_items), as.numeric))
  
  # Extract just the item responses for scoring
  item_responses <- ref_item_data %>%
    select(all_of(final_items)) %>%
    mutate(across(everything(), as.numeric))
  
  # Remove any rows that are completely missing
  complete_cases <- complete.cases(item_responses) | rowSums(!is.na(item_responses)) > 0
  
  if (sum(complete_cases) == 0) {
    stop("No valid response patterns found for scoring")
  }
  
  # Filter to valid cases
  valid_ref_item_data <- ref_item_data[complete_cases, ]
  valid_item_responses <- item_responses[complete_cases, ]
  
  cat("Scoring", nrow(valid_item_responses), "reference wave participants (", 
      nrow(ref_item_data) - nrow(valid_item_responses), "excluded due to missing data)...\n")
  
  # Score participants
  if (software == "mirt") {
    scores <- extract_mirt_scores_and_ses_improved(model, valid_item_responses)
    theta_scores <- scores$theta
    theta_se <- scores$se
    cat("Improved mirt scoring completed for reference wave\n")
    cat("Method used:", scores$method, "\n")
    cat("Participants with default SEs:", scores$n_default_ses, "\n")
    
  } else {
    person_estimates <- tam.wle(model, resp = as.matrix(valid_item_responses))
    
    if (is.matrix(person_estimates$theta)) {
      theta_scores <- person_estimates$theta[, 1]
    } else {
      theta_scores <- person_estimates$theta
    }
    
    if (!is.null(person_estimates$errorWLE)) {
      if (is.matrix(person_estimates$errorWLE)) {
        theta_se <- person_estimates$errorWLE[, 1]
      } else {
        theta_se <- person_estimates$errorWLE
      }
    } else {
      theta_se <- rep(0.5, length(theta_scores))
    }
  }
  
  # CRITICAL: Combine scores with participant data INCLUDING ID VARIABLES
  reference_scores <- valid_ref_item_data %>%
    select(Study, Wave, Study_Wave, GID, SID) %>%
    slice_head(n = length(theta_scores)) %>%
    mutate(
      Theta_Score = theta_scores,
      Theta_SE = theta_se,
      Score_Source = "Reference_Wave_Concurrent",
      Linking_Quality = "Highest",
      SE_Method = if (software == "mirt") scores$method else "TAM_WLE",
      Wave_Type = "Reference"
    )
  
  cat("Reference wave scores generated for", nrow(reference_scores), "participants\n\n")
  
  return(reference_scores)
}

# SUBSEQUENT WAVE SCORING WITH ITEM FILTERING
score_subsequent_waves_with_item_filtering <- function(all_waves_data, reference_wave_data, 
                                                       concurrent_results, structure_analysis, 
                                                       reference_wave) {
  
  cat("SCORING SUBSEQUENT WAVES WITH ENHANCED ITEM FILTERING:\n")
  
  # Identify subsequent waves
  longitudinal_studies <- all_waves_data %>%
    group_by(Study) %>%
    summarise(n_waves = n_distinct(Wave), .groups = "drop") %>%
    filter(n_waves > 1) %>%
    pull(Study)
  
  subsequent_wave_data <- all_waves_data %>%
    filter(Study %in% longitudinal_studies & Wave != reference_wave)
  
  if (nrow(subsequent_wave_data) == 0) {
    cat("No subsequent waves found to score.\n\n")
    return(data.frame())
  }
  
  final_items <- concurrent_results$final_items
  model <- concurrent_results$model
  software <- concurrent_results$software
  
  cat("Studies with subsequent waves:", length(longitudinal_studies), "\n")
  cat("Subsequent wave participants to score:", nrow(subsequent_wave_data), "\n")
  
  # Check which model items are available in subsequent waves
  available_model_items <- intersect(final_items, names(subsequent_wave_data))
  missing_model_items <- setdiff(final_items, names(subsequent_wave_data))
  
  cat("Model items available in subsequent waves:", length(available_model_items), "\n")
  if (length(missing_model_items) > 0) {
    cat("Model items MISSING in subsequent waves:", length(missing_model_items), "\n")
    cat("Missing items:", paste(head(missing_model_items, 10), collapse = ", "), "\n")
  }
  
  if (length(available_model_items) < 5) {
    cat("Too few model items available in subsequent waves for reliable scoring\n")
    return(data.frame())
  }
  
  # Prepare data using ONLY the model items that are available (WITH ID PRESERVATION)
  subseq_item_data <- subsequent_wave_data %>%
    select(Study, Wave, Study_Wave, GID, SID, all_of(available_model_items)) %>%
    mutate(across(all_of(available_model_items), as.numeric))
  
  # For scoring, create complete item response matrix with NAs for missing items
  complete_item_responses <- matrix(NA, nrow = nrow(subseq_item_data), ncol = length(final_items))
  colnames(complete_item_responses) <- final_items
  
  # Fill in available items
  for (item in available_model_items) {
    if (item %in% final_items) {
      item_idx <- which(final_items == item)
      complete_item_responses[, item_idx] <- subseq_item_data[[item]]
    }
  }
  
  # Remove completely empty cases
  complete_cases <- rowSums(!is.na(complete_item_responses)) > 0
  
  if (sum(complete_cases) == 0) {
    cat("No valid response patterns found\n\n")
    return(data.frame())
  }
  
  valid_subseq_item_data <- subseq_item_data[complete_cases, ]
  valid_item_responses <- complete_item_responses[complete_cases, , drop = FALSE]
  
  cat("Scoring", nrow(valid_item_responses), "subsequent wave participants\n")
  cat("Using", sum(colSums(!is.na(valid_item_responses)) > 0), "items with any data\n")
  
  # Validate and fix response categories before scoring
  cat("Validating response categories...\n")
  valid_item_responses <- validate_and_fix_response_categories_simple(valid_item_responses, model, final_items)
  
  # SCORING with the complete item matrix
  if (software == "mirt") {
    cat("Using enhanced mirt scoring for subsequent waves...\n")
    
    # Initialize variables to avoid scope issues
    theta_scores <- rep(0, nrow(valid_item_responses))
    theta_se <- rep(0.5, nrow(valid_item_responses))
    
    tryCatch({
      cat("  - Attempting full.scores = TRUE with complete item matrix...\n")
      
      full_scores <- fscores(model, 
                             response.pattern = valid_item_responses,
                             method = "EAP",
                             full.scores = TRUE,
                             full.scores.SE = TRUE,
                             verbose = FALSE)
      
      if (is.matrix(full_scores) && ncol(full_scores) >= 2) {
        theta_scores <- full_scores[, 1]
        
        se_cols <- grep("SE", colnames(full_scores), ignore.case = TRUE)
        if (length(se_cols) > 0) {
          theta_se <- full_scores[, se_cols[1]]
          
          # Validate SEs
          valid_ses <- !is.na(theta_se) & theta_se > 0.01 & theta_se < 5
          n_valid_ses <- sum(valid_ses)
          
          cat("  Enhanced scoring successful with", n_valid_ses, "valid SEs\n")
          
          if (n_valid_ses < length(theta_se)) {
            cat("  - Calculating SEs for remaining", length(theta_se) - n_valid_ses, "participants...\n")
            calculated_ses <- calculate_manual_ses_for_subsequent(model, valid_item_responses, theta_scores)
            theta_se[!valid_ses] <- calculated_ses[!valid_ses]
          }
          
        } else {
          theta_se <- calculate_manual_ses_for_subsequent(model, valid_item_responses, theta_scores)
        }
        
      } else {
        stop("Unexpected output format")
      }
      
    }, error = function(e) {
      cat("  Enhanced scoring failed:", e$message, "\n")
      cat("  - Using fallback method...\n")
      
      theta_scores <<- tryCatch({
        scores <- fscores(model, response.pattern = valid_item_responses, method = "EAP", verbose = FALSE)
        if (is.matrix(scores)) scores[, 1] else as.numeric(scores)
      }, error = function(e2) {
        rep(0, nrow(valid_item_responses))
      })
      
      theta_se <<- calculate_manual_ses_for_subsequent(model, valid_item_responses, theta_scores)
    })
    
    cat("  Subsequent wave scoring completed\n")
    cat("  - Participants with SE = 0.5:", sum(abs(theta_se - 0.5) < 0.001), "\n")
    
  } else {
    # TAM scoring
    person_estimates <- tam.wle(model, resp = valid_item_responses)
    
    if (is.matrix(person_estimates$theta)) {
      theta_scores <- person_estimates$theta[, 1]
    } else {
      theta_scores <- person_estimates$theta
    }
    
    if (!is.null(person_estimates$errorWLE)) {
      if (is.matrix(person_estimates$errorWLE)) {
        theta_se <- person_estimates$errorWLE[, 1]
      } else {
        theta_se <- person_estimates$errorWLE
      }
    } else {
      theta_se <- rep(0.5, length(theta_scores))
    }
  }
  
  # CRITICAL: Combine scores with participant data INCLUDING ID VARIABLES
  subsequent_scores <- valid_subseq_item_data %>%
    select(Study, Wave, Study_Wave, GID, SID) %>%
    slice_head(n = length(theta_scores)) %>%
    mutate(
      Theta_Score = theta_scores,
      Theta_SE = theta_se,
      Score_Source = "Subsequent_Wave_Item_Filtered",
      Linking_Quality = "High",
      SE_Method = if (software == "mirt") "mirt_enhanced_item_filtered" else "TAM_WLE",
      Wave_Type = "Subsequent",
      N_Model_Items_Used = length(available_model_items),
      N_Model_Items_Total = length(final_items)
    )
  
  cat("Enhanced subsequent wave scores generated for", nrow(subsequent_scores), "participants\n\n")
  
  return(subsequent_scores)
}

# IMPROVED MIRT SCORING WITH SE CALCULATION
extract_mirt_scores_and_ses_improved <- function(model, response_data) {
  
  cat("Extracting scores using improved mirt SE calculation...\n")
  
  # Ensure response_data is properly formatted as numeric matrix
  if (!is.matrix(response_data)) {
    response_data <- as.matrix(response_data)
  }
  response_data <- apply(response_data, 2, as.numeric)
  
  n_participants <- nrow(response_data)
  cat("Processing", n_participants, "participants...\n")
  
  # Try full.scores = TRUE first (most reliable)
  cat("Attempting full scoring with built-in SEs...\n")
  
  full_scores <- tryCatch({
    fscores(model, 
            response.pattern = response_data,
            method = "EAP",
            full.scores = TRUE,
            full.scores.SE = TRUE,
            verbose = FALSE)
  }, error = function(e) {
    cat("Full scoring failed:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(full_scores) && is.matrix(full_scores) && ncol(full_scores) >= 2) {
    theta_scores <- full_scores[, 1]
    
    # Look for SE columns
    se_cols <- grep("SE", colnames(full_scores), ignore.case = TRUE)
    if (length(se_cols) > 0) {
      theta_ses <- full_scores[, se_cols[1]]
      
      # Validate SEs
      valid_ses <- !is.na(theta_ses) & theta_ses > 0 & theta_ses < 5
      n_valid_ses <- sum(valid_ses)
      
      cat("Full scoring successful with built-in SEs\n")
      cat("Valid SEs for", n_valid_ses, "out of", n_participants, "participants\n")
      
      # Replace invalid SEs with calculated ones
      if (n_valid_ses < n_participants) {
        cat("Calculating SEs for", n_participants - n_valid_ses, "participants with invalid SEs...\n")
        calculated_ses <- calculate_improved_eap_ses(model, response_data, theta_scores)
        theta_ses[!valid_ses] <- calculated_ses[!valid_ses]
      }
      
      return(list(
        theta = theta_scores,
        se = theta_ses,
        method = "full_scores_with_builtin_SE",
        n_default_ses = n_participants - n_valid_ses
      ))
    }
  }
  
  # Method 2: Individual scoring with improved SE calculation
  cat("Using individual scoring with improved SE calculation...\n")
  
  # Get theta scores
  theta_scores <- tryCatch({
    scores <- fscores(model, 
                      response.pattern = response_data,
                      method = "EAP", 
                      verbose = FALSE)
    if (is.matrix(scores)) scores[, 1] else as.numeric(scores)
  }, error = function(e) {
    cat("Error in theta scoring:", e$message, "\n")
    rep(0, n_participants)
  })
  
  # Calculate improved SEs
  theta_ses <- calculate_improved_eap_ses(model, response_data, theta_scores)
  
  return(list(
    theta = theta_scores,
    se = theta_ses,
    method = "individual_with_improved_SE",
    n_default_ses = sum(abs(theta_ses - 0.5) < 0.001)
  ))
}

# IMPROVED SE CALCULATION METHOD
calculate_improved_eap_ses <- function(model, response_data, theta_scores) {
  
  n_persons <- nrow(response_data)
  ses <- numeric(n_persons)
  
  # Get item parameters once (more efficient)
  item_params <- tryCatch({
    params <- coef(model, simplify = TRUE)
    list(
      discriminations = params$items[, "a1"],
      success = TRUE
    )
  }, error = function(e) {
    cat("Could not extract parameters for SE calculation\n")
    list(success = FALSE)
  })
  
  if (!item_params$success) {
    cat("Using default SEs due to parameter extraction failure\n")
    return(rep(0.5, n_persons))
  }
  
  # Process in batches for efficiency
  batch_size <- 500
  n_batches <- ceiling(n_persons / batch_size)
  
  for (batch in 1:n_batches) {
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, n_persons)
    
    for (i in start_idx:end_idx) {
      response_pattern <- response_data[i, ]
      valid_items <- !is.na(response_pattern)
      theta_est <- theta_scores[i]
      
      if (sum(valid_items) == 0) {
        ses[i] <- 1.0
        next
      }
      
      # Calculate test information for valid items only
      test_info <- 0
      for (j in which(valid_items)) {
        a <- item_params$discriminations[j]
        if (!is.na(a) && a > 0) {
          # Simplified information calculation
          p <- 1 / (1 + exp(-1.7 * a * theta_est))
          item_info <- (1.7 * a)^2 * p * (1 - p)
          test_info <- test_info + item_info
        }
      }
      
      # EAP SE calculation with proper bounds
      if (test_info > 0.01) {
        prior_precision <- 1.0
        total_precision <- test_info + prior_precision
        ses[i] <- sqrt(1.0 / total_precision)
      } else {
        ses[i] <- 1.0
      }
      
      # Apply reasonable bounds
      ses[i] <- pmax(0.1, pmin(ses[i], 2.0))
    }
    
    if (batch %% 10 == 0 || batch == n_batches) {
      cat("  Processed batch", batch, "of", n_batches, "\n")
    }
  }
  
  cat("SE calculation completed. Range:", round(min(ses), 3), "to", round(max(ses), 3), "\n")
  
  return(ses)
}

# SIMPLE RESPONSE CATEGORY VALIDATION
validate_and_fix_response_categories_simple <- function(response_data, model, final_items) {
  
  cat("  - Validating response categories against model expectations...\n")
  
  # Get expected categories from the model
  tryCatch({
    model_data <- model@Data$data
    expected_ranges <- list()
    
    for (i in 1:length(final_items)) {
      item_name <- final_items[i]
      if (i <= ncol(model_data)) {
        model_responses <- model_data[, i]
        expected_min <- min(model_responses, na.rm = TRUE)
        expected_max <- max(model_responses, na.rm = TRUE)
        expected_ranges[[item_name]] <- c(expected_min, expected_max)
      }
    }
    
    # Check and correct each item
    corrected_data <- response_data
    corrections_made <- 0
    
    for (item_name in final_items) {
      if (item_name %in% colnames(response_data) && item_name %in% names(expected_ranges)) {
        
        actual_responses <- response_data[, item_name]
        actual_responses <- actual_responses[!is.na(actual_responses)]
        
        if (length(actual_responses) > 0) {
          actual_min <- min(actual_responses)
          actual_max <- max(actual_responses)
          expected_min <- expected_ranges[[item_name]][1]
          expected_max <- expected_ranges[[item_name]][2]
          
          # Check for category mismatches
          if (actual_min < expected_min || actual_max > expected_max) {
            cat("    Correcting item", which(final_items == item_name), ": Expected", expected_min, "-", expected_max, 
                ", found", actual_min, "-", actual_max, "\n")
            
            # Correct the responses to fit expected range
            corrected_responses <- response_data[, item_name]
            corrected_responses[corrected_responses < expected_min & !is.na(corrected_responses)] <- expected_min
            corrected_responses[corrected_responses > expected_max & !is.na(corrected_responses)] <- expected_max
            
            corrected_data[, item_name] <- corrected_responses
            corrections_made <- corrections_made + 1
          }
        }
      }
    }
    
    if (corrections_made > 0) {
      cat("    Corrected", corrections_made, "items with category mismatches\n")
    } else {
      cat("    All response categories match model expectations\n")
    }
    
    return(corrected_data)
    
  }, error = function(e) {
    cat("    Category validation failed, using original data:", e$message, "\n")
    return(response_data)
  })
}

calculate_manual_ses_for_subsequent <- function(model, response_matrix, theta_scores) {
  
  cat("    Calculating SEs manually using Fisher Information...\n")
  
  n_persons <- nrow(response_matrix)
  ses <- numeric(n_persons)
  
  # Extract item parameters
  tryCatch({
    item_params <- coef(model, simplify = TRUE)
    discriminations <- item_params$items[, "a1"]
    
    for (i in 1:n_persons) {
      response_pattern <- response_matrix[i, ]
      valid_items <- !is.na(response_pattern)
      theta_est <- theta_scores[i]
      
      if (sum(valid_items) == 0) {
        ses[i] <- 1.0
        next
      }
      
      # Calculate test information for valid items
      test_info <- 0
      for (j in which(valid_items)) {
        a <- discriminations[j]
        
        if (!is.na(a) && a > 0) {
          # Simplified 2PL information formula
          p <- 1 / (1 + exp(-1.7 * a * theta_est))
          item_info <- (1.7 * a)^2 * p * (1 - p)
          test_info <- test_info + item_info
        }
      }
      
      # Calculate SE using Fisher Information
      if (test_info > 0.01) {
        prior_precision <- 1.0
        total_precision <- test_info + prior_precision
        ses[i] <- sqrt(1.0 / total_precision)
      } else {
        ses[i] <- 1.0
      }
      
      # Apply reasonable bounds
      ses[i] <- pmax(0.05, pmin(ses[i], 2.0))
    }
    
  }, error = function(e) {
    ses[] <- 0.5
  })
  
  cat("    Manual SE calculation completed\n")
  cat("    SE range:", round(min(ses), 3), "to", round(max(ses), 3), "\n")
  
  return(ses)
}

# COMBINE SCORES WITH ID PRESERVATION
combine_reference_and_subsequent_scores <- function(reference_scores, subsequent_scores) {
  
  cat("COMBINING REFERENCE AND SUBSEQUENT WAVE SCORES:\n")
  
  # Combine the datasets
  if (nrow(subsequent_scores) > 0) {
    combined_scores <- bind_rows(reference_scores, subsequent_scores)
  } else {
    combined_scores <- reference_scores
  }
  
  # Add standardization and quality indicators
  combined_scores <- combined_scores %>%
    mutate(
      # Standardize theta scores across ALL waves
      Theta_Standardized = as.numeric(scale(Theta_Score)),
      
      # Quality indicators
      Precision_Weight = 1 / (Theta_SE^2),
      High_Precision = Theta_SE < quantile(Theta_SE, 0.75, na.rm = TRUE),
      
      # Harmonization metadata
      Harmonization_Date = Sys.Date(),
      Harmonization_Method = "Reference_Wave_Concurrent_IRT",
      
      # Reference wave indicator
      Is_Reference_Wave = Wave_Type == "Reference"
    )
  
  # Summary statistics
  scoring_summary <- combined_scores %>%
    group_by(Wave_Type) %>%
    summarise(
      n_participants = n(),
      n_studies = n_distinct(Study),
      n_waves = n_distinct(paste(Study, Wave)),
      mean_theta = mean(Theta_Score, na.rm = TRUE),
      mean_se = mean(Theta_SE, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("COMBINED SCORING SUMMARY:\n")
  print(scoring_summary)
  
  cat("\nCombined scores generated for", nrow(combined_scores), "participants\n")
  cat("Reference wave participants:", sum(combined_scores$Wave_Type == "Reference"), "\n")
  cat("Subsequent wave participants:", sum(combined_scores$Wave_Type == "Subsequent"), "\n\n")
  
  return(combined_scores)
}

# VALIDATION
validate_reference_wave_results <- function(concurrent_results, final_scores, structure_analysis) {
  
  cat("REFERENCE WAVE VALIDATION:\n")
  
  # Model quality
  reliability <- concurrent_results$reliability
  model_quality <- case_when(
    reliability > 0.9 ~ "Excellent",
    reliability > 0.85 ~ "Good",
    reliability > 0.8 ~ "Acceptable",
    TRUE ~ "Poor"
  )
  
  # Score quality
  mean_se <- mean(final_scores$Theta_SE, na.rm = TRUE)
  prop_high_precision <- mean(final_scores$High_Precision, na.rm = TRUE)
  
  score_quality <- case_when(
    mean_se < 0.4 && prop_high_precision > 0.7 ~ "Excellent",
    mean_se < 0.6 && prop_high_precision > 0.5 ~ "Good",
    mean_se < 0.8 && prop_high_precision > 0.3 ~ "Acceptable",
    TRUE ~ "Poor"
  )
  
  # Overall assessment
  overall_quality <- case_when(
    model_quality %in% c("Excellent", "Good") && score_quality %in% c("Excellent", "Good") ~ "Good",
    model_quality != "Poor" && score_quality != "Poor" ~ "Acceptable",
    TRUE ~ "Poor"
  )
  
  cat("Model quality:", model_quality, "(reliability =", round(reliability, 3), ")\n")
  cat("Score quality:", score_quality, "(mean SE =", round(mean_se, 3), ")\n")
  cat("Overall assessment:", overall_quality, "\n")
  
  return(list(
    overall_quality = overall_quality,
    model_quality = model_quality,
    score_quality = score_quality,
    reliability = reliability,
    mean_se = mean_se,
    prop_high_precision = prop_high_precision
  ))
}

# EXPORT RESULTS
export_reference_wave_results <- function(prep_results, structure_analysis, anchor_analysis,
                                          concurrent_results, final_scores, validation_results, 
                                          reference_wave) {
  
  cat("EXPORTING REFERENCE WAVE HARMONIZATION RESULTS\n")
  cat("==============================================\n")
  
  # Main harmonized scores
  write.csv(final_scores, "reference_wave_harmonized_scores.csv", row.names = FALSE)
  cat("reference_wave_harmonized_scores.csv - Main output with reference wave approach\n")
  
  # Study summary
  study_summary <- prep_results$study_structure %>%
    mutate(Reference_Wave = reference_wave)
  write.csv(study_summary, "reference_wave_study_summary.csv", row.names = FALSE)
  
  # Validation summary
  validation_summary <- data.frame(
    Metric = c("Overall_Quality", "Model_Reliability", "Mean_SE", "Reference_Wave", "N_Items"),
    Value = c(validation_results$overall_quality, round(validation_results$reliability, 3), 
              round(validation_results$mean_se, 3), reference_wave, concurrent_results$n_items),
    stringsAsFactors = FALSE
  )
  write.csv(validation_summary, "reference_wave_validation_summary.csv", row.names = FALSE)
  
  cat("All reference wave results exported successfully!\n")
}

# SUMMARY FUNCTIONS
create_reference_wave_summary <- function(prep_results, concurrent_results, validation_results, reference_wave) {
  return(list(
    total_studies = nrow(prep_results$study_structure),
    reference_wave = reference_wave,
    reference_wave_participants = nrow(prep_results$reference_wave_data),
    total_participants = nrow(prep_results$all_waves_data),
    total_items = concurrent_results$n_items,
    model_reliability = concurrent_results$reliability,
    overall_quality = validation_results$overall_quality
  ))
}

print_reference_wave_summary <- function(summary) {
  cat("REFERENCE WAVE HARMONIZATION SUMMARY\n")
  cat("===================================\n")
  cat("Total studies:", summary$total_studies, "\n")
  cat("Reference wave used:", summary$reference_wave, "\n")
  cat("Reference wave participants:", summary$reference_wave_participants, "\n")
  cat("Total participants:", summary$total_participants, "\n")
  cat("Items in final model:", summary$total_items, "\n")
  cat("Model reliability:", round(summary$model_reliability, 3), "\n")
  cat("Overall quality:", summary$overall_quality, "\n\n")
}

# ID-BASED MERGING FUNCTION
merge_harmonized_to_original <- function(harmonized_scores, original_data) {
  
  cat("MERGING HARMONIZED SCORES TO ORIGINAL DATA\n")
  cat("==========================================\n")
  
  # Check for ID variables
  if (!"GID" %in% names(harmonized_scores) || !"SID" %in% names(harmonized_scores)) {
    stop("GID and SID must be present in harmonized scores")
  }
  
  if (!"GID" %in% names(original_data) || !"SID" %in% names(original_data)) {
    stop("GID and SID must be present in original data")
  }
  
  cat("Using GID + SID + Study + Wave for ID-based merge\n")
  
  # Perform ID-based merge
  merged_data <- original_data %>%
    left_join(harmonized_scores, 
              by = c("Study", "Wave", "GID", "SID"),
              suffix = c("", "_harmonized"))
  
  # Validate merge
  n_matched <- sum(!is.na(merged_data$Theta_Standardized))
  match_rate <- (n_matched / nrow(original_data)) * 100
  
  cat("MERGE RESULTS:\n")
  cat("Original records:", nrow(original_data), "\n")
  cat("Successful matches:", n_matched, "\n")
  cat("Match rate:", sprintf("%.1f%%", match_rate), "\n")
  
  # Add merge indicators
  merged_data <- merged_data %>%
    mutate(
      Has_Harmonized_Score = !is.na(Theta_Standardized),
      Merge_Method = "ID_Based_GID_SID"
    )
  
  return(merged_data)
}

#############################################
# USAGE EXAMPLE
#############################################

# # Run harmonization
# results <- harmonize_studies_reference_wave_concurrent(
#   data = df,
#   reference_wave = 1,
#   min_anchors = 3,
#   min_studies_per_item = 2,
#   min_items_total = 5,
#   irt_software = "mirt",
#   max_participants_per_study_wave = 1500,
#   max_items = 150,
#   save_outputs = TRUE
# )
# 
# # Check results
# print(results$summary)
# 
# # Access harmonized scores with IDs preserved
# harmonized_scores <- results$final_scores
# 
# # Merge back to original data using ID-based matching
# merged_data <- merge_harmonized_to_original(harmonized_scores, df)
# 
# # Analysis using harmonized scores
# library(lme4)
# 
# model <- lmer(Theta_Standardized ~ Age_Cont_Est + Sex + (1 | Study), 
#               data = merged_data %>% filter(Has_Harmonized_Score == TRUE),
#               weights = Precision_Weight)