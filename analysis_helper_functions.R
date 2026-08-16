# This R file contains analysis helper functions
library(nnet)


# Function to score the hill climb using glm
glm.bic <- function(node, parents, data, args) {

   if (length(parents) == 0){
     model = as.formula(paste(node, "~ 1"))
   }else{
     model = as.formula(paste(node, "~", paste(parents, collapse = "+")))
   }
  
   if (is.numeric(data[, node])){
     return(BIC(glm(model, data = data, family = "gaussian")) / 2)
    }else if (is.factor(data[, node])){
      if (length(unique(data[,node]))<=2){
        return(BIC(glm(model, data = data, family = "binomial")) / 2)
      }else{
        return(BIC(multinom(model, data = data))/2)
     }
   }
 }

# Function to get the tp, fp, and fn label per edge
apply_labels <- function(target_arcs, predicted_arcs, arcs = FALSE) {
  predicted_arcs <- predicted_arcs %>% mutate(Predicted = TRUE)
  target_arcs <- target_arcs %>% mutate(Label = TRUE)

  labeled_predictions <- predicted_arcs %>%
    full_join(target_arcs, by = c("Var1", "Var2")) %>%
    mutate(Label = replace_na(Label, FALSE), Predicted = replace_na(Predicted, FALSE))

  tp_edges <- labeled_predictions %>% filter(Predicted == TRUE & Label == TRUE)
  fp_edges <- labeled_predictions %>% filter(Predicted == TRUE & Label == FALSE)
  fn_edges <- labeled_predictions %>% filter(Predicted == FALSE & Label == TRUE)

  tp <- tp_edges %>% nrow()
  fp <- fp_edges %>% nrow()
  fn <- fn_edges %>% nrow()

  if (arcs == FALSE) {
    return(list(tp = tp, fp = fp, fn = fn))
  } else {
    # Return with data frames for fp and fn that have "from" and "to" columns
    tp_df <- tp_edges %>%
      select(Var1, Var2) %>%
      rename(from = Var1, to = Var2)
    fp_df <- fp_edges %>%
      select(Var1, Var2) %>%
      rename(from = Var1, to = Var2)
    fn_df <- fn_edges %>%
      select(Var1, Var2) %>%
      rename(from = Var1, to = Var2)

    return(list(tp = tp_df, fp = fp_df, fn = fn_df))
  }
}

# Scoring function to get performance metrics
score <- function(target_arcs, current_arcs, node_names, arcs = FALSE) {
  # target_arcs: data frame of valid edges from the edge_list
  # current_arcs: data frame of predicted edges
  # node_names: names of the nodes
  target_dag <- empty.graph(node_names)
  current_dag <- empty.graph(node_names)
  
  arcs(target_dag, check.cycles = FALSE) <- as.matrix(target_arcs[, c("Var1", "Var2")])
  arcs(current_dag, check.cycles = FALSE) <- as.matrix(current_arcs[, c("Var1", "Var2")])
  
  # r <- compare(target_dag, current_dag, arcs = FALSE)
  r <- apply_labels(target_arcs, current_arcs, arcs = FALSE)
  
  f1 <- (2 * r$tp) / (2 * r$tp + r$fp + r$fn)
  tpr <- r$tp / (r$tp + r$fn)
  fdr <- r$fp / (r$tp + r$fp)
  # hamming distance converts a DAG to a undirected graph first.
  hd <- hamming(current_dag, target_dag)
  shd <- shd(current_dag, target_dag, cpdag = FALSE)
  
  if (arcs == TRUE) {
    # r <- compare(target_dag, current_dag, arcs = TRUE)
    r <- apply_labels(target_arcs, current_arcs, arcs = TRUE)
  }
  
  return(c(r, f1=f1, tpr=tpr, fdr=fdr, hd=hd, shd=shd))
}

# visualization of causal model in relation to the expert dag
visualize <- function(target_arcs, current_arcs, node_names, main = "", sub = "", fontsize = 15.5, show_fn = TRUE) {
  # target_arcs: data frame of valid edges from the edge_list
  # current_arcs: data frame of predicted edges
  # node_names: names of the nodes
  node_mapping <- c(
    "Age" = "Age",
    "Alcohol" = "Alcohol",
    "Anxiety" = "Anxiety",
    "Catastrophizing" = "Cat",
    "Depression" = "Depr",
    "Education" = "Edu",
    "Fear_avoidance" = "FA",
    "Financial_level" = "Fin",
    "Obesity" = "Obesity",
    "Sex" = "Sex",
    "Sleep disturbance" = "Sleep",
    "Smoking" = "Smoking",
    "CCI" = "CCI",
    "PEG" = "PEG"
  )
  
  # Apply node mapping if provided
  if (!is.null(node_mapping)) {
    # Rename nodes in target_arcs
    target_arcs$Var1 <- ifelse(target_arcs$Var1 %in% names(node_mapping),
                                node_mapping[target_arcs$Var1],
                                target_arcs$Var1)
    target_arcs$Var2 <- ifelse(target_arcs$Var2 %in% names(node_mapping),
                                node_mapping[target_arcs$Var2],
                                target_arcs$Var2)

    # Rename nodes in current_arcs
    current_arcs$Var1 <- ifelse(current_arcs$Var1 %in% names(node_mapping),
                                 node_mapping[current_arcs$Var1],
                                 current_arcs$Var1)
    current_arcs$Var2 <- ifelse(current_arcs$Var2 %in% names(node_mapping),
                                 node_mapping[current_arcs$Var2],
                                 current_arcs$Var2)

    # Rename node_names
    node_names <- ifelse(node_names %in% names(node_mapping),
                         node_mapping[node_names],
                         node_names)
  }
  
  target_dag <- empty.graph(node_names)
  current_dag <- empty.graph(node_names)

  arcs(target_dag, check.cycles = FALSE) <- as.matrix(target_arcs[, c("Var1", "Var2")])
  arcs(current_dag, check.cycles = FALSE) <- as.matrix(current_arcs[, c("Var1", "Var2")])

  r <- compare(target_dag, current_dag, arcs = TRUE)

  # Extract FP and FN matrices
  fp_matrix <- as.matrix(r$fp[,c("from", "to")])
  fn_matrix <- as.matrix(r$fn[,c("from", "to")])

  if (show_fn == TRUE) {
    # Create combined DAG with current edges + false negative edges
    all_edges <- rbind(
      arcs(current_dag),                         # Current proto edges
      as.matrix(fn_matrix[,c("from","to")])      # Missing expert edges (FN)
    ) %>%
      as.data.frame()
  } else {
    # only the edges, not any fn
    all_edges <- arcs(current_dag) %>% 
      as.data.frame()
  }

  combined_dag <- empty.graph(node_names)
  arcs(combined_dag, check.cycles=FALSE) <- as.matrix(all_edges)

  # Create combined DAG plot
  combined_dag_plot <- graphviz.plot(combined_dag, layout="circo", render = FALSE, fontsize = fontsize, main = main, sub = sub)

  # Create edge names for both FP and FN matrices
  fp_edge_names <- paste(fp_matrix[,1], fp_matrix[,2], sep = "~")
  fn_edge_names <- paste(fn_matrix[,1], fn_matrix[,2], sep = "~")

  # Get existing edge names from combined DAG
  combined_edge_names <- edgeNames(combined_dag_plot)

  # Color logic: FP = red, FN = blue, TP = black
  edge_colors <- ifelse(combined_edge_names %in% fp_edge_names, "red",
                       ifelse(combined_edge_names %in% fn_edge_names, "blue", "black"))

  # Line type logic: FP and FN = dashed, TP = solid
  edge_lty <- ifelse(combined_edge_names %in% fp_edge_names, "dashed",
                    ifelse(combined_edge_names %in% fn_edge_names, "dashed", "solid"))

  # Set edge attributes using named vectors
  names(edge_colors) <- combined_edge_names
  names(edge_lty) <- combined_edge_names

  edgeRenderInfo(combined_dag_plot) <- list(
    col = edge_colors,
    lty = edge_lty
  )

  # Render the combined plot
  renderGraph(combined_dag_plot)

  # Add legend
#  legend("topright",
#         legend = c("True Positive", "False Positive", "False Negative"),
#         col = c("black", "red", "blue"),
#         lty = c("solid", "dashed", "dashed"),
#         lwd = 2,
#         bty = "n")

  # Return comparison results for further analysis
  # return(r)
}


# Proto-model directed edge set at a given strength threshold 
proto_at_threshold <- function(threshold) {
  f <- proto_df %>%
    filter(strength >= threshold, direction > 0) %>%
    filter(to != "Sex", to != "Age", from != "PEG")
  names(f)[1:2] <- c("Var1", "Var2")
  f %>%
    full_join(f, by = c("Var1" = "Var2", "Var2" = "Var1"), suffix = c(".fwd", ".rev")) %>%
    filter(replace_na(direction.fwd, 0) > replace_na(direction.rev, 0)) %>%
    transmute(Var1, Var2, source = "proto")
}

# Proto's continuous per-directed-edge weight = P(edge exists) * P(this direction | exists),
# used by the weighted-strength integration in instead of a pre-thresholded edge set
proto_weight_edges <- function() {
  proto_df %>%
    filter(to != "Sex", to != "Age", from != "PEG") %>%
    transmute(Var1 = from, Var2 = to, weight = strength * direction) %>%
    filter(weight > 0)
}

# Knowledge-system raw query data, loaded directly by backbone (model) and system
ks_cache <- list()
get_ks <- function(model, system) {
  key <- paste(model, system)
  if (is.null(ks_cache[[key]])) {
    ks_cache[[key]] <<- read_csv(file.path("data", model, paste0(system, ".csv")), show_col_types = FALSE) %>%
      select(Var1, Var2, Plausibility, Association, Temporality)
  }
  ks_cache[[key]]
}

# Directed edges a knowledge system "votes for" under a given query type
ks_vote_edges <- function(model, system, query_col, label = paste(system, query_col, sep = "_")) {
  df <- get_ks(model, system)
  df %>% filter(.data[[query_col]]) %>% transmute(Var1, Var2, source = label)
}

# A knowledge-system source as a weighted-edge table (weight_val applied to every TRUE row)
ks_weight_edges <- function(model, system, query_col, weight_val = 1) {
  df <- get_ks(model, system)
  df %>% filter(.data[[query_col]]) %>% transmute(Var1, Var2, weight = weight_val)
}

# Oracle weight: a source's own standalone precision (TP / (TP+FP)) against the expert
# graph, computed from its raw per-direction TRUE/FALSE calls and the per-row ground-truth
# Label already present in the raw query files. Upper-bound diagnostic only.
ks_oracle_precision <- function(model, system, query_col) {
  raw <- read_csv(file.path("data", model, paste0(system, ".csv")), show_col_types = FALSE)
  pos <- raw %>% filter(.data[[query_col]])
  tp <- sum(pos$Label, na.rm = TRUE)
  fp <- sum(!pos$Label, na.rm = TRUE)
  if (tp + fp == 0) return(0)
  tp / (tp + fp)
}

# Generic N-of-M vote fusion: retain an unordered pair only if >= min_sources distinct
# sources voted for it in some direction (1 = union, total #sources = intersection);
# resolve direction by vote count, breaking ties via the proto-model's direction oracle.
fuse_sources <- function(source_edge_dfs, min_sources = 1, tiebreak = c("proto", "drop")) {
  tiebreak <- match.arg(tiebreak)

  all_votes <- bind_rows(source_edge_dfs) %>% distinct(Var1, Var2, source)
  if (nrow(all_votes) == 0) return(tibble(Var1 = character(), Var2 = character()))

  all_votes <- all_votes %>%
    mutate(pair_id = map2_chr(Var1, Var2, ~ paste(sort(c(.x, .y)), collapse = "__")))

  pair_support <- all_votes %>% distinct(pair_id, source) %>% count(pair_id, name = "n_src")
  keep_pairs <- pair_support %>% filter(n_src >= min_sources) %>% pull(pair_id)

  dir_counts <- all_votes %>%
    filter(pair_id %in% keep_pairs) %>%
    count(pair_id, Var1, Var2, name = "votes")

  resolved <- dir_counts %>%
    group_by(pair_id) %>%
    mutate(top_votes = max(votes), n_top = sum(votes == top_votes)) %>%
    ungroup()

  clear <- resolved %>% filter(n_top == 1, votes == top_votes) %>% select(pair_id, Var1, Var2)
  tied <- resolved %>% filter(n_top > 1)

  if (nrow(tied) > 0 && tiebreak == "proto") {
    tb <- tied %>%
      left_join(proto_direction_pref, by = c("Var1", "Var2")) %>%
      mutate(direction = replace_na(direction, -1)) %>%
      group_by(pair_id) %>%
      filter(direction == max(direction)) %>%
      slice(1) %>%
      ungroup() %>%
      select(pair_id, Var1, Var2)
    clear <- bind_rows(clear, tb)
  }
  # tiebreak == "drop": tied pairs are simply excluded (already not in `clear`)

  clear %>% distinct(Var1, Var2)
}

# Weighted-strength fusion: sum weights per directed pair, keep pairs whose summed weight
# >= tau, resolve direction by comparing summed weight of the two directions (ties as above)
fuse_weighted <- function(weight_dfs, tau, tiebreak = c("proto", "drop")) {
  tiebreak <- match.arg(tiebreak)

  all_w <- bind_rows(weight_dfs)
  if (nrow(all_w) == 0) return(tibble(Var1 = character(), Var2 = character()))

  agg <- all_w %>%
    group_by(Var1, Var2) %>%
    summarise(total_weight = sum(weight), .groups = "drop") %>%
    filter(total_weight >= tau) %>%
    mutate(pair_id = map2_chr(Var1, Var2, ~ paste(sort(c(.x, .y)), collapse = "__")))

  if (nrow(agg) == 0) return(tibble(Var1 = character(), Var2 = character()))

  resolved <- agg %>%
    group_by(pair_id) %>%
    mutate(top_w = max(total_weight), n_top = sum(total_weight == top_w)) %>%
    ungroup()

  clear <- resolved %>% filter(n_top == 1, total_weight == top_w) %>% select(pair_id, Var1, Var2)
  tied <- resolved %>% filter(n_top > 1)

  if (nrow(tied) > 0 && tiebreak == "proto") {
    tb <- tied %>%
      left_join(proto_direction_pref, by = c("Var1", "Var2")) %>%
      mutate(direction = replace_na(direction, -1)) %>%
      group_by(pair_id) %>%
      filter(direction == max(direction)) %>%
      slice(1) %>%
      ungroup() %>%
      select(pair_id, Var1, Var2)
    clear <- bind_rows(clear, tb)
  }

  clear %>% distinct(Var1, Var2)
}

# Scoring wrapper around score()
score_edges <- function(edge_df) {
  r <- score(expert_edges, edge_df, node_names)
  tibble(n_edges = nrow(edge_df), tp = r$tp, fp = r$fp, fn = r$fn,
         f1 = r$f1, tpr = r$tpr, fdr = r$fdr, hd = r$hd, shd = r$shd)
}