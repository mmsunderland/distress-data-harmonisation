# REFINED DIMENSIONALITY ANALYSIS FOR PARTIAL ANCHOR OVERLAP
# ============================================================
# Focus on:
# 1. Network connectivity analysis (anchor overlap structure)
# 2. Polychoric correlations between all distress items
# 3. Within-study unidimensional assessment
# 4. Within-study bifactor analysis and reliability

library(dplyr)
library(tidyr)
library(mirt)
library(psych)
library(ggplot2)
library(igraph)
library(polycor)
library(corrplot)
library(lavaan)  # Added for CFA with ordinal data

# ==============================================================================
# MAIN FUNCTION
# ==============================================================================

analyze_dimensionality_refined <- function(results,
                                           min_studies_per_anchor = 2,
                                           min_participants_for_study_analysis = 100,
                                           max_items_for_polychor = 50,
                                           save_outputs = TRUE) {
  
  cat("=== REFINED DIMENSIONALITY ANALYSIS ===\n")
  cat("=======================================\n\n")
  
  # Extract data
  calibration_data <- results$concurrent_calibration$calibration_data
  final_items <- results$concurrent_calibration$final_items
  anchor_items <- results$anchor_analysis$all_anchors
  
  cat("Dataset overview:\n")
  cat("  Total items:", length(final_items), "\n")
  cat("  Anchor items:", length(anchor_items), "\n")
  cat("  Study-specific items:", length(final_items) - length(anchor_items), "\n")
  cat("  Studies:", n_distinct(calibration_data$Study), "\n")
  cat("  Participants:", nrow(calibration_data), "\n\n")
  
  # ANALYSIS 1: Network Connectivity
  cat("ANALYSIS 1: NETWORK CONNECTIVITY\n")
  cat("================================\n")
  network_results <- analyze_anchor_network(
    calibration_data,
    anchor_items,
    min_studies_per_anchor
  )
  
  # ANALYSIS 2: Polychoric Correlations
  cat("\nANALYSIS 2: POLYCHORIC CORRELATIONS\n")
  cat("===================================\n")
  polychor_results <- analyze_polychoric_correlations(
    calibration_data,
    final_items,
    max_items_for_polychor
  )
  
  # ANALYSIS 3: Within-Study Unidimensional Assessment
  cat("\nANALYSIS 3: WITHIN-STUDY UNIDIMENSIONAL ASSESSMENT\n")
  cat("==================================================\n")
  within_study_unidim <- analyze_within_study_unidimensionality(
    calibration_data,
    final_items,
    min_participants_for_study_analysis
  )
  
  # ANALYSIS 4: Within-Study Bifactor Analysis
  cat("\nANALYSIS 4: WITHIN-STUDY BIFACTOR ANALYSIS\n")
  cat("==========================================\n")
  within_study_bifactor <- analyze_within_study_bifactor(
    calibration_data,
    final_items,
    min_participants_for_study_analysis
  )
  
  # Generate Summary
  cat("\nGENERATING SUMMARY\n")
  cat("==================\n")
  summary_results <- generate_refined_summary(
    network_results,
    polychor_results,
    within_study_unidim,
    within_study_bifactor
  )
  
  # Create Visualizations
  cat("\nGENERATING VISUALIZATIONS\n")
  cat("=========================\n")
  plots <- create_refined_plots(
    network_results,
    polychor_results,
    within_study_unidim,
    within_study_bifactor
  )
  
  # Export Results
  if (save_outputs) {
    export_refined_results(
      network_results,
      polychor_results,
      within_study_unidim,
      within_study_bifactor,
      summary_results,
      plots
    )
  }
  
  # Compile all results
  all_results <- list(
    network_analysis = network_results,
    polychoric_correlations = polychor_results,
    within_study_unidimensional = within_study_unidim,
    within_study_bifactor = within_study_bifactor,
    summary = summary_results,
    plots = plots
  )
  
  cat("\n✅ REFINED DIMENSIONALITY ANALYSIS COMPLETE!\n")
  cat("============================================\n")
  print_refined_summary(all_results)
  
  return(all_results)
}

# ==============================================================================
# ANALYSIS 1: NETWORK CONNECTIVITY
# ==============================================================================

analyze_anchor_network <- function(data, anchor_items, min_studies_per_anchor) {
  
  cat("Analyzing anchor item network structure...\n\n")
  
  # Create study-anchor coverage matrix
  study_anchor_matrix <- data %>%
    select(Study, all_of(anchor_items)) %>%
    group_by(Study) %>%
    summarise(
      across(all_of(anchor_items), 
             ~sum(!is.na(.x)) >= 3,
             .names = "{.col}"),
      n_participants = n(),
      .groups = "drop"
    )
  
  studies <- study_anchor_matrix$Study
  anchor_matrix <- as.matrix(study_anchor_matrix[, anchor_items])
  rownames(anchor_matrix) <- studies
  
  # Anchor statistics
  anchor_stats <- data.frame(
    Anchor = anchor_items,
    N_Studies = colSums(anchor_matrix),
    Prop_Studies = colMeans(anchor_matrix),
    Studies = sapply(anchor_items, function(a) {
      paste(studies[anchor_matrix[, a]], collapse = ", ")
    }),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(N_Studies))
  
  # Study statistics
  study_stats <- data.frame(
    Study = studies,
    N_Anchors = rowSums(anchor_matrix),
    N_Participants = study_anchor_matrix$n_participants,
    Anchors = sapply(studies, function(s) {
      paste(anchor_items[anchor_matrix[s, ]], collapse = ", ")
    }),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(N_Anchors))
  
  # Pairwise study overlap
  n_studies <- length(studies)
  pairwise_overlap <- matrix(0, n_studies, n_studies)
  rownames(pairwise_overlap) <- colnames(pairwise_overlap) <- studies
  
  study_pair_details <- list()
  
  for (i in 1:(n_studies-1)) {
    for (j in (i+1):n_studies) {
      shared_anchors <- anchor_items[anchor_matrix[i, ] & anchor_matrix[j, ]]
      n_shared <- length(shared_anchors)
      pairwise_overlap[i, j] <- n_shared
      pairwise_overlap[j, i] <- n_shared
      
      if (n_shared > 0) {
        pair_name <- paste(studies[i], studies[j], sep = "_")
        study_pair_details[[pair_name]] <- list(
          study1 = studies[i],
          study2 = studies[j],
          n_shared = n_shared,
          shared_anchors = shared_anchors
        )
      }
    }
  }
  
  # Create network graph
  edges <- data.frame()
  for (i in 1:(n_studies-1)) {
    for (j in (i+1):n_studies) {
      if (pairwise_overlap[i, j] > 0) {
        edges <- rbind(edges, data.frame(
          from = studies[i],
          to = studies[j],
          weight = pairwise_overlap[i, j]
        ))
      }
    }
  }
  
  if (nrow(edges) == 0) {
    cat("❌ ERROR: No study connections found!\n")
    cat("   Studies cannot be harmonized without shared anchor items.\n\n")
    
    return(list(
      connected = FALSE,
      n_components = n_studies,
      anchor_stats = anchor_stats,
      study_stats = study_stats,
      pairwise_overlap = pairwise_overlap
    ))
  }
  
  # Build igraph network
  network_graph <- graph_from_data_frame(edges, directed = FALSE, vertices = studies)
  
  # Network metrics
  is_connected <- is_connected(network_graph)
  n_components <- components(network_graph)$no
  graph_density <- edge_density(network_graph)
  
  # Additional metrics if connected
  diameter_val <- if (is_connected) diameter(network_graph) else NA
  mean_distance <- if (is_connected) mean_distance(network_graph) else NA
  
  # Study centrality
  study_degree <- degree(network_graph)
  study_betweenness <- betweenness(network_graph)
  
  study_centrality <- data.frame(
    Study = names(study_degree),
    Degree = study_degree,
    Betweenness = study_betweenness,
    N_Anchors = study_stats$N_Anchors[match(names(study_degree), study_stats$Study)],
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(Degree))
  
  # Component membership
  component_membership <- components(network_graph)$membership
  
  cat("NETWORK ANALYSIS RESULTS:\n")
  cat("  Studies:", n_studies, "\n")
  cat("  Anchor items:", length(anchor_items), "\n")
  cat("  Study connections:", nrow(edges), "\n")
  cat("  Fully connected:", is_connected)
  if (is_connected) {
    cat(" ✅\n")
  } else {
    cat(" 🚨\n")
  }
  cat("  Number of components:", n_components, "\n")
  cat("  Graph density:", round(graph_density, 3), "\n")
  if (!is.na(diameter_val)) {
    cat("  Network diameter:", diameter_val, "\n")
  }
  cat("\n")
  
  cat("HUB STUDIES (highest connectivity):\n")
  print(head(study_centrality, 5))
  cat("\n")
  
  if (!is_connected) {
    cat("⚠️  WARNING: Network is not fully connected!\n")
    cat("   Some studies cannot be compared directly.\n")
    cat("   Consider:\n")
    cat("   1. Adding more anchor items\n")
    cat("   2. Analyzing components separately\n")
    cat("   3. Using alternative linking methods\n\n")
  }
  
  return(list(
    connected = is_connected,
    n_components = n_components,
    graph_density = graph_density,
    diameter = diameter_val,
    mean_distance = mean_distance,
    network_graph = network_graph,
    anchor_stats = anchor_stats,
    study_stats = study_stats,
    study_centrality = study_centrality,
    pairwise_overlap = pairwise_overlap,
    study_pair_details = study_pair_details,
    component_membership = component_membership,
    n_studies = n_studies,
    n_anchors = length(anchor_items),
    anchor_matrix = anchor_matrix
  ))
}

# ==============================================================================
# ANALYSIS 2: POLYCHORIC CORRELATIONS
# ==============================================================================

analyze_polychoric_correlations <- function(data, all_items, max_items) {
  
  cat("Computing polychoric correlations between distress items...\n\n")
  
  # Select items for analysis
  if (length(all_items) > max_items) {
    cat("  Large number of items (", length(all_items), ").\n")
    cat("  Selecting top", max_items, "items by coverage for polychoric analysis...\n")
    
    item_coverage <- data %>%
      select(all_of(all_items)) %>%
      summarise(across(everything(), ~sum(!is.na(.x)))) %>%
      pivot_longer(everything(), names_to = "Item", values_to = "N_Valid") %>%
      arrange(desc(N_Valid))
    
    selected_items <- item_coverage$Item[1:max_items]
    
  } else {
    selected_items <- all_items
  }
  
  cat("  Items for polychoric analysis:", length(selected_items), "\n")
  
  # Prepare data
  item_data <- data %>%
    select(all_of(selected_items)) %>%
    mutate(across(everything(), as.numeric))
  
  # Initialize matrices
  n_items <- length(selected_items)
  polychor_matrix <- matrix(NA, n_items, n_items)
  rownames(polychor_matrix) <- colnames(polychor_matrix) <- selected_items
  diag(polychor_matrix) <- 1
  
  pairwise_n <- matrix(NA, n_items, n_items)
  rownames(pairwise_n) <- colnames(pairwise_n) <- selected_items
  
  # Compute polychoric correlations
  cat("  Computing pairwise polychoric correlations...\n")
  cat("  This may take several minutes for", n_items, "items\n")
  
  computed <- 0
  failed <- 0
  
  for (i in 1:(n_items-1)) {
    for (j in (i+1):n_items) {
      
      item1 <- selected_items[i]
      item2 <- selected_items[j]
      
      # Get valid cases
      valid_cases <- complete.cases(item_data[, c(item1, item2)])
      n_valid <- sum(valid_cases)
      pairwise_n[i, j] <- n_valid
      pairwise_n[j, i] <- n_valid
      
      # Compute polychoric if sufficient data
      if (n_valid >= 30) {
        
        tryCatch({
          
          x1 <- item_data[[item1]][valid_cases]
          x2 <- item_data[[item2]][valid_cases]
          
          # Check for sufficient variation
          if (length(unique(x1)) >= 2 && length(unique(x2)) >= 2) {
            
            polychor_result <- polychor(x1, x2, std.err = FALSE)
            
            if (!is.na(polychor_result) && abs(polychor_result) <= 1) {
              polychor_matrix[i, j] <- polychor_result
              polychor_matrix[j, i] <- polychor_result
              computed <- computed + 1
            } else {
              failed <- failed + 1
            }
          } else {
            failed <- failed + 1
          }
          
        }, error = function(e) {
          failed <<- failed + 1
        })
      } else {
        failed <- failed + 1
      }
    }
    
    if (i %% 10 == 0) {
      cat("    Progress:", round(i / length(selected_items) * 100), "%\n")
    }
  }
  
  cat("  ✅ Computation complete\n")
  cat("     Successful:", computed, "pairs\n")
  cat("     Failed/insufficient data:", failed, "pairs\n\n")
  
  # Summary statistics
  polychor_values <- polychor_matrix[upper.tri(polychor_matrix)]
  polychor_values <- polychor_values[!is.na(polychor_values)]
  
  cat("POLYCHORIC CORRELATION SUMMARY:\n")
  cat("  Valid correlations:", length(polychor_values), "\n")
  cat("  Mean:", round(mean(polychor_values, na.rm = TRUE), 3), "\n")
  cat("  Median:", round(median(polychor_values, na.rm = TRUE), 3), "\n")
  cat("  SD:", round(sd(polychor_values, na.rm = TRUE), 3), "\n")
  cat("  Range:", round(min(polychor_values, na.rm = TRUE), 3), "to",
      round(max(polychor_values, na.rm = TRUE), 3), "\n")
  
  # Percentiles
  quantiles <- quantile(polychor_values, probs = c(0.10, 0.25, 0.50, 0.75, 0.90), na.rm = TRUE)
  cat("  Percentiles:\n")
  cat("    10th:", round(quantiles[1], 3), "\n")
  cat("    25th:", round(quantiles[2], 3), "\n")
  cat("    50th:", round(quantiles[3], 3), "\n")
  cat("    75th:", round(quantiles[4], 3), "\n")
  cat("    90th:", round(quantiles[5], 3), "\n")
  
  # Check for problematic correlations
  n_negative <- sum(polychor_values < 0, na.rm = TRUE)
  n_weak <- sum(polychor_values < 0.10, na.rm = TRUE)
  n_strong <- sum(polychor_values > 0.70, na.rm = TRUE)
  
  cat("  Negative correlations:", n_negative, "\n")
  cat("  Weak correlations (<0.10):", n_weak, "\n")
  cat("  Strong correlations (>0.70):", n_strong, "\n\n")
  
  if (mean(polychor_values, na.rm = TRUE) > 0.30) {
    cat("  ✅ Average correlation suggests items measure related constructs\n\n")
  } else {
    cat("  ⚠️  Average correlation is low - items may be weakly related\n\n")
  }
  
  return(list(
    polychor_matrix = polychor_matrix,
    pairwise_n = pairwise_n,
    selected_items = selected_items,
    n_items = length(selected_items),
    summary_stats = list(
      mean = mean(polychor_values, na.rm = TRUE),
      median = median(polychor_values, na.rm = TRUE),
      sd = sd(polychor_values, na.rm = TRUE),
      min = min(polychor_values, na.rm = TRUE),
      max = max(polychor_values, na.rm = TRUE),
      quantiles = quantiles,
      n_negative = n_negative,
      n_weak = n_weak,
      n_strong = n_strong
    )
  ))
}

# ==============================================================================
# ANALYSIS 3: WITHIN-STUDY UNIDIMENSIONALITY (REVISED WITH LAVAAN CFA)
# ==============================================================================

analyze_within_study_unidimensionality <- function(data, all_items, min_participants) {
  
  cat("Assessing unidimensionality within each study using CFA for ordinal data...\n\n")
  
  studies <- unique(data$Study)
  within_study_results <- list()
  
  for (study in studies) {
    
    cat("Study:", study, "\n")
    
    # Get data for this study
    study_data <- data %>%
      filter(Study == study) %>%
      select(all_of(all_items)) %>%
      mutate(across(everything(), as.numeric))
    
    # Identify items with data in this study
    item_coverage <- colSums(!is.na(study_data))
    available_items <- names(item_coverage[item_coverage >= 10])
    
    cat("  Participants:", nrow(study_data), "\n")
    cat("  Items with data:", length(available_items), "\n")
    
    if (nrow(study_data) < min_participants) {
      cat("  ⚠️  Insufficient participants (need", min_participants, ")\n\n")
      within_study_results[[study]] <- list(
        study = study,
        sufficient_data = FALSE,
        n_participants = nrow(study_data),
        n_items = length(available_items)
      )
      next
    }
    
    if (length(available_items) < 5) {
      cat("  ⚠️  Insufficient items (need 5+)\n\n")
      within_study_results[[study]] <- list(
        study = study,
        sufficient_data = FALSE,
        n_participants = nrow(study_data),
        n_items = length(available_items)
      )
      next
    }
    
    # Get complete cases
    study_complete <- study_data %>%
      select(all_of(available_items)) %>%
      filter(complete.cases(.))
    
    cat("  Complete cases:", nrow(study_complete), "\n")
    
    if (nrow(study_complete) < 50) {
      cat("  ⚠️  Insufficient complete cases\n\n")
      within_study_results[[study]] <- list(
        study = study,
        sufficient_data = FALSE,
        n_participants = nrow(study_data),
        n_items = length(available_items),
        n_complete = nrow(study_complete)
      )
      next
    }
    
    # Unidimensional CFA assessment with ordinal indicators
    cat("  Running one-factor CFA with ordinal indicators (DWLS)...\n")
    
    # Build CFA model specification
    cfa_model <- paste0("distress =~ ", paste(available_items, collapse = " + "))
    
    # Fit CFA with ordinal indicators and DWLS estimator
    cfa_fit <- tryCatch({
      cfa(cfa_model, 
          data = study_complete, 
          ordered = available_items,  # Specify ordinal indicators
          estimator = "DWLS",         # Use DWLS for ordinal data
          std.lv = TRUE)
    }, error = function(e) NULL)
    
    if (!is.null(cfa_fit)) {
      # Extract fit indices
      fit_measures <- fitMeasures(cfa_fit, c("rmsea", "cfi", "tli", "chisq", "df", "pvalue"))
      
      cat("    RMSEA:", round(fit_measures["rmsea"], 3))
      if (fit_measures["rmsea"] < 0.08) cat(" ✅") else cat(" ⚠️")
      cat("\n")
      
      cat("    CFI:", round(fit_measures["cfi"], 3))
      if (fit_measures["cfi"] > 0.90) cat(" ✅") else cat(" ⚠️")
      cat("\n")
      
      cat("    TLI:", round(fit_measures["tli"], 3))
      if (fit_measures["tli"] > 0.90) cat(" ✅") else cat(" ⚠️")
      cat("\n")
      
      within_study_results[[study]] <- list(
        study = study,
        sufficient_data = TRUE,
        n_participants = nrow(study_data),
        n_items = length(available_items),
        n_complete = nrow(study_complete),
        available_items = available_items,
        unidim_rmsea = fit_measures["rmsea"],
        unidim_cfi = fit_measures["cfi"],
        unidim_tli = fit_measures["tli"],
        chisq = fit_measures["chisq"],
        df = fit_measures["df"],
        pvalue = fit_measures["pvalue"],
        cfa_fit = cfa_fit
      )
    } else {
      cat("    ❌ CFA failed to converge\n")
      within_study_results[[study]] <- list(
        study = study,
        sufficient_data = FALSE,
        n_participants = nrow(study_data),
        n_items = length(available_items),
        n_complete = nrow(study_complete),
        error = "CFA did not converge"
      )
    }
    
    cat("\n")
  }
  
  # Summary across studies
  successful_studies <- sum(sapply(within_study_results, function(x) 
    !is.null(x$sufficient_data) && x$sufficient_data))
  
  cat("WITHIN-STUDY UNIDIMENSIONALITY SUMMARY:\n")
  cat("  Studies analyzed:", length(studies), "\n")
  cat("  Successful analyses:", successful_studies, "\n")
  
  if (successful_studies > 0) {
    rmsea_values <- sapply(within_study_results, function(x) {
      if(!is.null(x$unidim_rmsea)) x$unidim_rmsea else NA
    })
    rmsea_values <- rmsea_values[!is.na(rmsea_values)]
    
    cfi_values <- sapply(within_study_results, function(x) {
      if(!is.null(x$unidim_cfi)) x$unidim_cfi else NA
    })
    cfi_values <- cfi_values[!is.na(cfi_values)]
    
    tli_values <- sapply(within_study_results, function(x) {
      if(!is.null(x$unidim_tli)) x$unidim_tli else NA
    })
    tli_values <- tli_values[!is.na(tli_values)]
    
    cat("  Mean RMSEA:", round(mean(rmsea_values, na.rm = TRUE), 3), "\n")
    cat("  Mean CFI:", round(mean(cfi_values, na.rm = TRUE), 3), "\n")
    cat("  Mean TLI:", round(mean(tli_values, na.rm = TRUE), 3), "\n")
    cat("  Studies with RMSEA < 0.08:", sum(rmsea_values < 0.08, na.rm = TRUE), 
        "out of", length(rmsea_values), "\n")
    cat("  Studies with CFI > 0.90:", sum(cfi_values > 0.90, na.rm = TRUE),
        "out of", length(cfi_values), "\n")
    
    if (mean(rmsea_values < 0.08, na.rm = TRUE) > 0.75) {
      cat("  ✅ Most studies show acceptable unidimensional fit\n")
    } else {
      cat("  ⚠️  Mixed unidimensionality across studies\n")
    }
  }
  
  cat("\n")
  
  return(list(
    study_results = within_study_results,
    n_studies = length(studies),
    n_successful = successful_studies
  ))
}

# ==============================================================================
# ANALYSIS 4: WITHIN-STUDY BIFACTOR ANALYSIS
# ==============================================================================
# Uses psych::omega() with poly=TRUE to compute polychoric correlations
# before bifactor rotation. This is appropriate for ordinal categorical data.
# The minres (minimum residual) factor method works well with polychoric correlations.

analyze_within_study_bifactor <- function(data, all_items, min_participants) {
  
  cat("Performing bifactor analysis within each study...\n\n")
  
  studies <- unique(data$Study)
  bifactor_results <- list()
  
  for (study in studies) {
    
    cat("Study:", study, "\n")
    
    # Get data for this study
    study_data <- data %>%
      filter(Study == study) %>%
      select(all_of(all_items)) %>%
      mutate(across(everything(), as.numeric))
    
    # Identify items with data
    item_coverage <- colSums(!is.na(study_data))
    available_items <- names(item_coverage[item_coverage >= 10])
    
    cat("  Items available:", length(available_items), "\n")
    
    if (nrow(study_data) < min_participants || length(available_items) < 6) {
      cat("  ⚠️  Insufficient data for bifactor analysis\n\n")
      bifactor_results[[study]] <- list(
        study = study,
        sufficient_data = FALSE,
        n_participants = nrow(study_data),
        n_items = length(available_items)
      )
      next
    }
    
    # Get complete cases
    study_complete <- study_data %>%
      select(all_of(available_items)) %>%
      filter(complete.cases(.))
    
    cat("  Complete cases:", nrow(study_complete), "\n")
    
    if (nrow(study_complete) < 100) {
      cat("  ⚠️  Insufficient complete cases for bifactor\n\n")
      bifactor_results[[study]] <- list(
        study = study,
        sufficient_data = FALSE,
        n_participants = nrow(study_data),
        n_items = length(available_items),
        n_complete = nrow(study_complete)
      )
      next
    }
    
    # Omega (bifactor) analysis with polychoric correlations for ordinal data
    cat("  Running omega (bifactor) analysis with polychoric correlations...\n")
    
    omega_result <- tryCatch({
      omega(study_complete, 
            nfactors = min(3, floor(ncol(study_complete) / 3)),
            fm = "minres",        # Minimum residual method works better with polychoric
            poly = TRUE,          # Use polychoric correlations for ordinal data
            rotate = "oblimin",
            plot = FALSE)
    }, error = function(e) NULL)
    
    if (!is.null(omega_result)) {
      omega_h <- omega_result$omega_h
      omega_total <- omega_result$omega.tot
      pev_general <- omega_result$omega.group[1, 1]
      
      cat("    Omega hierarchical (ωh):", round(omega_h, 3))
      if (omega_h > 0.70) {
        cat(" ✅\n")
      } else if (omega_h > 0.50) {
        cat(" ⚠️\n")
      } else {
        cat(" 🚨\n")
      }
      
      cat("    Omega total:", round(omega_total, 3), "\n")
      cat("    % variance by general:", round(pev_general * 100, 1), "%\n")
      
      bifactor_results[[study]] <- list(
        study = study,
        sufficient_data = TRUE,
        n_participants = nrow(study_data),
        n_items = length(available_items),
        n_complete = nrow(study_complete),
        available_items = available_items,
        omega_h = omega_h,
        omega_total = omega_total,
        pev_general = pev_general,
        omega_result = omega_result
      )
    } else {
      cat("    ❌ Bifactor analysis failed\n")
      bifactor_results[[study]] <- list(
        study = study,
        sufficient_data = FALSE,
        n_participants = nrow(study_data),
        n_items = length(available_items),
        n_complete = nrow(study_complete),
        error = "Bifactor failed"
      )
    }
    
    cat("\n")
  }
  
  # Summary across studies
  successful_studies <- sum(sapply(bifactor_results, function(x) 
    !is.null(x$sufficient_data) && x$sufficient_data))
  
  cat("WITHIN-STUDY BIFACTOR SUMMARY:\n")
  cat("  Studies analyzed:", length(studies), "\n")
  cat("  Successful analyses:", successful_studies, "\n")
  
  if (successful_studies > 0) {
    omega_h_values <- sapply(bifactor_results, function(x) {
      if(!is.null(x$omega_h)) x$omega_h else NA
    })
    omega_h_values <- omega_h_values[!is.na(omega_h_values)]
    
    cat("  Mean omega hierarchical:", round(mean(omega_h_values, na.rm = TRUE), 3), "\n")
    cat("  Studies with strong general (ωh > 0.70):", 
        sum(omega_h_values > 0.70, na.rm = TRUE), 
        "out of", length(omega_h_values), "\n")
    cat("  Studies with moderate+ general (ωh > 0.50):",
        sum(omega_h_values > 0.50, na.rm = TRUE),
        "out of", length(omega_h_values), "\n")
    
    if (mean(omega_h_values > 0.70, na.rm = TRUE) > 0.75) {
      cat("  ✅ Most studies have strong general factors\n")
    } else if (mean(omega_h_values > 0.50, na.rm = TRUE) > 0.75) {
      cat("  ⚠️  Most studies have moderate general factors\n")
    } else {
      cat("  🚨 Many studies have weak general factors\n")
    }
  }
  
  cat("\n")
  
  return(list(
    study_results = bifactor_results,
    n_studies = length(studies),
    n_successful = successful_studies
  ))
}

# ==============================================================================
# SUMMARY GENERATION
# ==============================================================================

generate_refined_summary <- function(network_results, polychor_results, 
                                     within_study_unidim, within_study_bifactor) {
  
  # Network summary
  network_summary <- list(
    fully_connected = network_results$connected,
    n_components = network_results$n_components,
    n_studies = network_results$n_studies,
    n_anchors = network_results$n_anchors,
    graph_density = network_results$graph_density,
    diameter = network_results$diameter
  )
  
  # Polychoric summary
  polychor_summary <- list(
    n_items = polychor_results$n_items,
    mean_correlation = polychor_results$summary_stats$mean,
    median_correlation = polychor_results$summary_stats$median,
    sd_correlation = polychor_results$summary_stats$sd
  )
  
  # Within-study unidimensional summary
  if (within_study_unidim$n_successful > 0) {
    rmsea_values <- sapply(within_study_unidim$study_results, function(x) {
      if(!is.null(x$unidim_rmsea)) x$unidim_rmsea else NA
    })
    rmsea_values <- rmsea_values[!is.na(rmsea_values)]
    
    cfi_values <- sapply(within_study_unidim$study_results, function(x) {
      if(!is.null(x$unidim_cfi)) x$unidim_cfi else NA
    })
    cfi_values <- cfi_values[!is.na(cfi_values)]
    
    tli_values <- sapply(within_study_unidim$study_results, function(x) {
      if(!is.null(x$unidim_tli)) x$unidim_tli else NA
    })
    tli_values <- tli_values[!is.na(tli_values)]
    
    unidim_summary <- list(
      n_studies_analyzed = within_study_unidim$n_successful,
      mean_rmsea = mean(rmsea_values, na.rm = TRUE),
      mean_cfi = mean(cfi_values, na.rm = TRUE),
      mean_tli = mean(tli_values, na.rm = TRUE),
      prop_good_fit = mean(rmsea_values < 0.08, na.rm = TRUE)
    )
  } else {
    unidim_summary <- list(
      n_studies_analyzed = 0,
      mean_rmsea = NA,
      mean_cfi = NA,
      mean_tli = NA,
      prop_good_fit = NA
    )
  }
  
  # Within-study bifactor summary
  if (within_study_bifactor$n_successful > 0) {
    omega_h_values <- sapply(within_study_bifactor$study_results, function(x) {
      if(!is.null(x$omega_h)) x$omega_h else NA
    })
    omega_h_values <- omega_h_values[!is.na(omega_h_values)]
    
    bifactor_summary <- list(
      n_studies_analyzed = within_study_bifactor$n_successful,
      mean_omega_h = mean(omega_h_values, na.rm = TRUE),
      prop_strong_general = mean(omega_h_values > 0.70, na.rm = TRUE),
      prop_moderate_general = mean(omega_h_values > 0.50, na.rm = TRUE)
    )
  } else {
    bifactor_summary <- list(
      n_studies_analyzed = 0,
      mean_omega_h = NA,
      prop_strong_general = NA,
      prop_moderate_general = NA
    )
  }
  
  return(list(
    network = network_summary,
    polychoric = polychor_summary,
    within_study_unidim = unidim_summary,
    within_study_bifactor = bifactor_summary
  ))
}

# ==============================================================================
# VISUALIZATION FUNCTIONS
# ==============================================================================

create_refined_plots <- function(network_results, polychor_results,
                                 within_study_unidim, within_study_bifactor) {
  
  plots <- list()
  
  # 1. Network diagram
  if (network_results$connected || network_results$n_components <= 3) {
    tryCatch({
      plots$network <- plot_network_diagram(network_results)
    }, error = function(e) {
      cat("  ⚠️  Network diagram failed:", e$message, "\n")
    })
  }
  
  # 2. Polychoric correlation heatmap
  tryCatch({
    plots$polychor_heatmap <- plot_polychoric_heatmap(polychor_results)
  }, error = function(e) {
    cat("  ⚠️  Polychoric heatmap failed:", e$message, "\n")
  })
  
  # 3. Polychoric correlation distribution
  tryCatch({
    plots$polychor_distribution <- plot_polychoric_distribution(polychor_results)
  }, error = function(e) {
    cat("  ⚠️  Polychoric distribution plot failed:", e$message, "\n")
  })
  
  # 4. Within-study unidimensional fit
  if (within_study_unidim$n_successful > 0) {
    tryCatch({
      plots$within_study_unidim <- plot_within_study_unidim(within_study_unidim)
    }, error = function(e) {
      cat("  ⚠️  Unidimensional fit plot failed:", e$message, "\n")
    })
  }
  
  # 5. Within-study omega hierarchical
  if (within_study_bifactor$n_successful > 0) {
    tryCatch({
      plots$within_study_omega <- plot_within_study_omega(within_study_bifactor)
    }, error = function(e) {
      cat("  ⚠️  Omega plot failed:", e$message, "\n")
    })
  }
  
  # 6. Anchor coverage by study
  tryCatch({
    plots$anchor_coverage <- plot_anchor_coverage(network_results)
  }, error = function(e) {
    cat("  ⚠️  Anchor coverage plot failed:", e$message, "\n")
  })
  
  return(plots)
}

plot_network_diagram <- function(network_results) {
  
  g <- network_results$network_graph
  
  # Layout
  layout_coords <- layout_with_fr(g)
  
  # Vertex data
  vertex_df <- data.frame(
    name = V(g)$name,
    x = layout_coords[, 1],
    y = layout_coords[, 2],
    degree = degree(g),
    betweenness = betweenness(g)
  )
  
  # Edge data
  edge_list <- as_data_frame(g, "edges")
  edge_df <- data.frame()
  
  for (i in 1:nrow(edge_list)) {
    from_vertex <- vertex_df[vertex_df$name == edge_list$from[i], ]
    to_vertex <- vertex_df[vertex_df$name == edge_list$to[i], ]
    
    edge_df <- rbind(edge_df, data.frame(
      x = from_vertex$x,
      y = from_vertex$y,
      xend = to_vertex$x,
      yend = to_vertex$y,
      weight = edge_list$weight[i]
    ))
  }
  
  p <- ggplot() +
    geom_segment(data = edge_df, 
                 aes(x = x, y = y, xend = xend, yend = yend, size = weight),
                 alpha = 0.4, color = "gray40") +
    geom_point(data = vertex_df, 
               aes(x = x, y = y, size = degree),
               color = "#1976D2", alpha = 0.7) +
    geom_text(data = vertex_df, 
              aes(x = x, y = y, label = name),
              vjust = -1.5, size = 3.5, fontface = "bold") +
    scale_size_continuous(range = c(2, 10), name = "Connections") +
    labs(title = "Study Network via Shared Anchor Items") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "right")
  return(p)
}

plot_polychoric_heatmap <- function(polychor_results) {
  
  cor_matrix <- polychor_results$polychor_matrix
  
  # Check for NA values and handle them
  na_counts <- rowSums(is.na(cor_matrix))
  max_na_allowed <- ncol(cor_matrix) * 0.5  # Allow up to 50% missing
  
  # Filter out items with too many missing correlations
  if (any(na_counts > max_na_allowed)) {
    cat("  Note: Removing", sum(na_counts > max_na_allowed), 
        "items with >50% missing correlations from heatmap\n")
    keep_items <- na_counts <= max_na_allowed
    cor_matrix <- cor_matrix[keep_items, keep_items]
  }
  
  # Check if we still have sufficient data for clustering
  remaining_na_prop <- sum(is.na(cor_matrix)) / (nrow(cor_matrix) * ncol(cor_matrix))
  
  # Decide on ordering method based on missing data
  if (remaining_na_prop > 0.1) {
    # Too many NAs for reliable clustering - use original order
    order_method <- "original"
    cat("  Note: Using original item order due to missing correlations (", 
        round(remaining_na_prop * 100, 1), "% missing)\n")
  } else if (remaining_na_prop > 0) {
    # Few NAs - use AOE (angular order of eigenvectors) which is more robust
    order_method <- "AOE"
    # Replace remaining NAs with 0 for visualization
    cor_matrix[is.na(cor_matrix)] <- 0
  } else {
    # No NAs - can use hierarchical clustering
    order_method <- "hclust"
  }
  
  # Create heatmap
  corrplot(cor_matrix,
           method = "color",
           type = "upper",
           order = order_method,
           tl.col = "black",
           tl.cex = 0.7,
           cl.cex = 0.8,
           title = "Polychoric Correlations Between Items",
           mar = c(0, 0, 2, 0),
           na.label = "square",  # Show NA values as squares
           na.label.col = "white")
  
  p <- recordPlot()
  return(p)
}

plot_polychoric_distribution <- function(polychor_results) {
  
  cor_values <- polychor_results$polychor_matrix[upper.tri(polychor_results$polychor_matrix)]
  cor_values <- cor_values[!is.na(cor_values)]
  
  plot_data <- data.frame(correlation = cor_values)
  
  p <- ggplot(plot_data, aes(x = correlation)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "black") +
    geom_vline(xintercept = mean(cor_values, na.rm = TRUE), 
               color = "red", linetype = "dashed", size = 1) +
    labs(title = "Distribution of Polychoric Correlations",
         x = "Correlation",
         y = "Frequency",
         subtitle = paste0("Mean = ", round(mean(cor_values, na.rm = TRUE), 3))) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12))
  
  return(p)
}

plot_within_study_unidim <- function(within_study_unidim) {
  
  # Extract fit statistics
  study_fits <- data.frame()
  for (study in names(within_study_unidim$study_results)) {
    result <- within_study_unidim$study_results[[study]]
    if (!is.null(result$sufficient_data) && result$sufficient_data) {
      study_fits <- rbind(study_fits, data.frame(
        Study = study,
        RMSEA = result$unidim_rmsea,
        CFI = result$unidim_cfi,
        TLI = result$unidim_tli,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  if (nrow(study_fits) == 0) {
    return(NULL)
  }
  
  # Reshape for plotting
  study_fits_long <- study_fits %>%
    pivot_longer(cols = c(RMSEA, CFI, TLI),
                 names_to = "Fit_Index",
                 values_to = "Value")
  
  p <- ggplot(study_fits_long, aes(x = Study, y = Value, fill = Fit_Index)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ Fit_Index, scales = "free_y") +
    geom_hline(data = data.frame(Fit_Index = "RMSEA", y = 0.08),
               aes(yintercept = y), linetype = "dashed", color = "red") +
    geom_hline(data = data.frame(Fit_Index = c("CFI", "TLI"), y = 0.90),
               aes(yintercept = y), linetype = "dashed", color = "red") +
    labs(title = "Within-Study Unidimensional Fit (CFA with Ordinal Indicators)",
         subtitle = "Dashed lines show acceptable fit thresholds",
         x = "Study",
         y = "Fit Index Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          legend.position = "none")
  
  return(p)
}

plot_within_study_omega <- function(within_study_bifactor) {
  
  # Extract omega values
  omega_data <- data.frame()
  for (study in names(within_study_bifactor$study_results)) {
    result <- within_study_bifactor$study_results[[study]]
    if (!is.null(result$sufficient_data) && result$sufficient_data) {
      omega_data <- rbind(omega_data, data.frame(
        Study = study,
        Omega_H = result$omega_h,
        Omega_Total = result$omega_total,
        PEV_General = result$pev_general,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  if (nrow(omega_data) == 0) {
    return(NULL)
  }
  
  omega_data_long <- omega_data %>%
    pivot_longer(cols = c(Omega_H, Omega_Total),
                 names_to = "Omega_Type",
                 values_to = "Value")
  
  p <- ggplot(omega_data_long, aes(x = reorder(Study, Value), y = Value, fill = Omega_Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = 0.70, linetype = "dashed", color = "yellow") +
    geom_hline(yintercept = 0.50, linetype = "dashed", color = "purple") +
    coord_flip() +
    labs(title = "Within-Study Bifactor Reliability",
         subtitle = "Yellow line: Strong general factor (ωh > 0.70) | Purple line: Moderate (ωh > 0.50)",
         x = "Study",
         y = "Omega Value",
         fill = "Reliability Type") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 9))
  
  return(p)
}

plot_anchor_coverage <- function(network_results) {
  
  anchor_mat <- network_results$anchor_matrix
  
  heatmap_long <- as.data.frame(anchor_mat) %>%
    mutate(Study = rownames(anchor_mat)) %>%
    pivot_longer(-Study, names_to = "Anchor", values_to = "Present")
  
  p <- ggplot(heatmap_long, 
                                  aes(x = Anchor, y = Study, fill = Present)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_manual(values = c("TRUE" = "#2E7D32", "FALSE" = "#E0E0E0"),
                      labels = c("TRUE" = "Present", "FALSE" = "Absent"),
                      name = NULL) +
    labs(title = "Anchor Item Coverage by Study",
         x = "Anchor Items",
         y = "Studies") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
          plot.title = element_text(face = "bold"))
  
  return(p)
}

# ==============================================================================
# EXPORT FUNCTIONS
# ==============================================================================

export_refined_results <- function(network_results, polychor_results,
                                   within_study_unidim, within_study_bifactor,
                                   summary_results, plots) {
  
  cat("Exporting refined dimensionality analysis results...\n")
  
  # 1. Network analysis results
  write.csv(network_results$anchor_stats, "refined_anchor_statistics.csv", row.names = FALSE)
  write.csv(network_results$study_stats, "refined_study_statistics.csv", row.names = FALSE)
  write.csv(network_results$study_centrality, "refined_study_centrality.csv", row.names = FALSE)
  
  # 2. Polychoric correlations
  write.csv(polychor_results$polychor_matrix, "refined_polychoric_matrix.csv", row.names = TRUE)
  
  # 3. Within-study unidimensional results (updated for CFA)
  unidim_df <- data.frame()
  for (study in names(within_study_unidim$study_results)) {
    result <- within_study_unidim$study_results[[study]]
    unidim_df <- rbind(unidim_df, data.frame(
      Study = study,
      Sufficient_Data = result$sufficient_data,
      N_Participants = result$n_participants,
      N_Items = result$n_items,
      N_Complete = if(!is.null(result$n_complete)) result$n_complete else NA,
      RMSEA = if(!is.null(result$unidim_rmsea)) result$unidim_rmsea else NA,
      CFI = if(!is.null(result$unidim_cfi)) result$unidim_cfi else NA,
      TLI = if(!is.null(result$unidim_tli)) result$unidim_tli else NA,
      stringsAsFactors = FALSE
    ))
  }
  write.csv(unidim_df, "refined_within_study_unidimensional.csv", row.names = FALSE)
  
  # 4. Within-study bifactor results
  bifactor_df <- data.frame()
  for (study in names(within_study_bifactor$study_results)) {
    result <- within_study_bifactor$study_results[[study]]
    bifactor_df <- rbind(bifactor_df, data.frame(
      Study = study,
      Sufficient_Data = result$sufficient_data,
      N_Participants = result$n_participants,
      N_Items = result$n_items,
      N_Complete = if(!is.null(result$n_complete)) result$n_complete else NA,
      Omega_H = if(!is.null(result$omega_h)) result$omega_h else NA,
      Omega_Total = if(!is.null(result$omega_total)) result$omega_total else NA,
      PEV_General = if(!is.null(result$pev_general)) result$pev_general else NA,
      stringsAsFactors = FALSE
    ))
  }
  write.csv(bifactor_df, "refined_within_study_bifactor.csv", row.names = FALSE)
  
  # 5. Overall summary
  summary_df <- data.frame(
    Category = c("Network", "Network", "Network", "Network",
                 "Polychoric", "Polychoric", "Polychoric",
                 "Within_Study_Unidim", "Within_Study_Unidim", "Within_Study_Unidim", "Within_Study_Unidim",
                 "Within_Study_Bifactor", "Within_Study_Bifactor", "Within_Study_Bifactor"),
    Metric = c("Fully_Connected", "N_Components", "N_Studies", "Graph_Density",
               "N_Items", "Mean_Correlation", "Median_Correlation",
               "N_Studies_Analyzed", "Mean_RMSEA", "Mean_CFI", "Prop_Good_Fit",
               "N_Studies_Analyzed", "Mean_Omega_H", "Prop_Strong_General"),
    Value = c(
      summary_results$network$fully_connected,
      summary_results$network$n_components,
      summary_results$network$n_studies,
      summary_results$network$graph_density,
      summary_results$polychoric$n_items,
      summary_results$polychoric$mean_correlation,
      summary_results$polychoric$median_correlation,
      summary_results$within_study_unidim$n_studies_analyzed,
      if(!is.null(summary_results$within_study_unidim$mean_rmsea)) 
        summary_results$within_study_unidim$mean_rmsea else NA,
      if(!is.null(summary_results$within_study_unidim$mean_cfi))
        summary_results$within_study_unidim$mean_cfi else NA,
      if(!is.null(summary_results$within_study_unidim$prop_good_fit))
        summary_results$within_study_unidim$prop_good_fit else NA,
      summary_results$within_study_bifactor$n_studies_analyzed,
      if(!is.null(summary_results$within_study_bifactor$mean_omega_h))
        summary_results$within_study_bifactor$mean_omega_h else NA,
      if(!is.null(summary_results$within_study_bifactor$prop_strong_general))
        summary_results$within_study_bifactor$prop_strong_general else NA
    ),
    stringsAsFactors = FALSE
  )
  write.csv(summary_df, "refined_overall_summary.csv", row.names = FALSE)
  
  # 6. Save plots with error handling
  if ("network" %in% names(plots) && !is.null(plots$network)) {
    tryCatch({
      ggsave("refined_network_diagram.png", plots$network, width = 10, height = 8)
    }, error = function(e) {
      cat("  ⚠️  Could not save network diagram:", e$message, "\n")
    })
  }
  if ("polychor_heatmap" %in% names(plots) && !is.null(plots$polychor_heatmap)) {
    tryCatch({
      png("refined_polychoric_heatmap.png", width = 1200, height = 1000)
      replayPlot(plots$polychor_heatmap)
      dev.off()
    }, error = function(e) {
      cat("  ⚠️  Could not save polychoric heatmap:", e$message, "\n")
    })
  }
  if ("polychor_distribution" %in% names(plots) && !is.null(plots$polychor_distribution)) {
    tryCatch({
      ggsave("refined_polychoric_distribution.png", plots$polychor_distribution, width = 8, height = 6)
    }, error = function(e) {
      cat("  ⚠️  Could not save polychoric distribution:", e$message, "\n")
    })
  }
  if ("within_study_unidim" %in% names(plots) && !is.null(plots$within_study_unidim)) {
    tryCatch({
      ggsave("refined_within_study_unidim.png", plots$within_study_unidim, width = 10, height = 8)
    }, error = function(e) {
      cat("  ⚠️  Could not save unidimensional plot:", e$message, "\n")
    })
  }
  if ("within_study_omega" %in% names(plots) && !is.null(plots$within_study_omega)) {
    tryCatch({
      ggsave("refined_within_study_omega.png", plots$within_study_omega, width = 10, height = 8)
    }, error = function(e) {
      cat("  ⚠️  Could not save omega plot:", e$message, "\n")
    })
  }
  if ("anchor_coverage" %in% names(plots) && !is.null(plots$anchor_coverage)) {
    tryCatch({
      ggsave("refined_anchor_coverage.png", plots$anchor_coverage, width = 12, height = 8)
    }, error = function(e) {
      cat("  ⚠️  Could not save anchor coverage plot:", e$message, "\n")
    })
  }
  
  cat("✅ Results exported (check console for any plot save warnings)!\n\n")
}

# ==============================================================================
# PRINT SUMMARY
# ==============================================================================

print_refined_summary <- function(results) {
  
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("REFINED DIMENSIONALITY ANALYSIS SUMMARY\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  summary <- results$summary
  
  # Network
  cat("1. NETWORK CONNECTIVITY\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  cat("   Fully connected:", summary$network$fully_connected)
  if (summary$network$fully_connected) {
    cat(" ✅\n")
  } else {
    cat(" 🚨\n")
  }
  cat("   Number of components:", summary$network$n_components, "\n")
  cat("   Studies:", summary$network$n_studies, "\n")
  cat("   Graph density:", round(summary$network$graph_density, 3), "\n")
  if (!is.na(summary$network$diameter)) {
    cat("   Network diameter:", summary$network$diameter, "\n")
  }
  cat("\n")
  
  # Polychoric
  cat("2. POLYCHORIC CORRELATIONS\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  cat("   Items analyzed:", summary$polychoric$n_items, "\n")
  cat("   Mean correlation:", round(summary$polychoric$mean_correlation, 3))
  if (summary$polychoric$mean_correlation > 0.30) {
    cat(" ✅\n")
  } else if (summary$polychoric$mean_correlation > 0.20) {
    cat(" ⚠️\n")
  } else {
    cat(" 🚨\n")
  }
  cat("   Median correlation:", round(summary$polychoric$median_correlation, 3), "\n")
  cat("   SD:", round(summary$polychoric$sd_correlation, 3), "\n")
  cat("\n")
  
  # Within-study unidimensional
  cat("3. WITHIN-STUDY UNIDIMENSIONALITY (CFA WITH ORDINAL INDICATORS)\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  cat("   Studies analyzed:", summary$within_study_unidim$n_studies_analyzed, "\n")
  if (summary$within_study_unidim$n_studies_analyzed > 0) {
    cat("   Mean RMSEA:", round(summary$within_study_unidim$mean_rmsea, 3), "\n")
    cat("   Mean CFI:", round(summary$within_study_unidim$mean_cfi, 3), "\n")
    cat("   Mean TLI:", round(summary$within_study_unidim$mean_tli, 3), "\n")
    cat("   Proportion with good fit (RMSEA < 0.08):", 
        round(summary$within_study_unidim$prop_good_fit, 2))
    if (summary$within_study_unidim$prop_good_fit > 0.75) {
      cat(" ✅\n")
    } else if (summary$within_study_unidim$prop_good_fit > 0.50) {
      cat(" ⚠️\n")
    } else {
      cat(" 🚨\n")
    }
  }
  cat("\n")
  
  # Within-study bifactor
  cat("4. WITHIN-STUDY BIFACTOR (GENERAL FACTOR STRENGTH)\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  cat("   Studies analyzed:", summary$within_study_bifactor$n_studies_analyzed, "\n")
  if (summary$within_study_bifactor$n_studies_analyzed > 0) {
    cat("   Mean omega hierarchical:", round(summary$within_study_bifactor$mean_omega_h, 3), "\n")
    cat("   Proportion with strong general (ωh > 0.70):", 
        round(summary$within_study_bifactor$prop_strong_general, 2))
    if (summary$within_study_bifactor$prop_strong_general > 0.75) {
      cat(" ✅\n")
    } else if (summary$within_study_bifactor$prop_moderate_general > 0.75) {
      cat(" ⚠️\n")
    } else {
      cat(" 🚨\n")
    }
    cat("   Proportion with moderate+ general (ωh > 0.50):",
        round(summary$within_study_bifactor$prop_moderate_general, 2), "\n")
  }
  cat("\n")
  
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("\n")
}

# ==============================================================================
# USAGE INSTRUCTIONS
# ==============================================================================

# After running harmonization:
# results <- harmonize_studies_reference_wave_concurrent(...)
#
# Run refined dimensionality analysis:
 dim_results <- analyze_dimensionality_refined(
   results = results,
   min_studies_per_anchor = 2,
   min_participants_for_study_analysis = 100,
   max_items_for_polychor = 70,
   save_outputs = TRUE
 )
#
# WHAT THIS PROVIDES:
# 
# 1. NETWORK CONNECTIVITY:
#    - Study network diagram showing how studies connect via anchors
#    - Checks if network is fully connected (CRITICAL)
#    - Identifies hub studies and network metrics
#
# 2. POLYCHORIC CORRELATIONS:
#    - All pairwise correlations between distress items
#    - Correlation matrix and heatmap
#    - Distribution of correlations
#
# 3. WITHIN-STUDY UNIDIMENSIONALITY:
#    - One-factor CFA with ordinal indicators (DWLS estimator)
#    - Fit indices (RMSEA, CFI, TLI) for each study
#    - Shows if items within each study are unidimensional
#
# 4. WITHIN-STUDY BIFACTOR:
#    - Omega hierarchical (general factor strength) for each study
#    - Uses polychoric correlations (appropriate for ordinal data)
#    - Omega total and % variance explained by general factor
#    - Identifies studies with weak general factors
#
# INTERPRETING RESULTS:
# - Network must be connected for harmonization to work
# - Mean polychoric > 0.30 suggests items measure related construct
# - Within-study ωh > 0.70 = strong general factor within that study
# - Consistency across studies supports harmonization validity
