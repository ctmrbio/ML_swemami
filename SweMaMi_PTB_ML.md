# Prepping the ML models here but I will run them directly on Gandalf


```R
library(tidyverse)
library(scutr)
library(MASS)
library(mixOmics) 
library(vctrs) 
library(caret)  
library(xgboost, Ckmeans.1d.dp) 
library(dplyr) 
library(ggplot2) 
library(stringi) 
library(arsenal) 
library(e1071) 
library(kernlab) 
library(nnet) 
library(randomForest) 
library(pROC) 
library(MLeval) 
library(skimr) 
library(RANN) 
library(cluster) 
library(mclust) 
library(tibble)
```

 

### prepping the files to use (will do this here)


```R
ancombc_input_files_csv <- list.files(path="/PATH", pattern=".csv", all.files=T, full.names=T)

for (files in ancombc_input_files_csv) {
    
    file_name <-  str_split(files,"/")
    file_name <- data.frame(file_name)
    file_name_2 <- str_split(file_name[10,],".csv")
    file_name_2 <- data.frame(file_name_2)
    name <- paste(file_name_2[1,],"",sep="")
    print(name)
    csv_file <- read.csv(files,sep=",")
    csv_file <- csv_file[,-c(1)]
    assign(name,csv_file)
}

metadata_Q1_for_ML <- read.csv("/PATH/metadata_Q1_for_ML.csv",sep=",")
metadata_Q2_for_ML <- read.csv("/PATH/metadata_Q2_for_ML.csv",sep=",")
```

    [1] "deltaCLR_F_species_metadata_ML"
    [1] "deltaCLR_F_V_species_metadata_ML"
    [1] "deltaCLR_V_species_metadata_ML"
    [1] "F1_case_control_clr"
    [1] "F1_cc_alpha_key"
    [1] "F1_clr_alpha_metadata"
    [1] "F2_case_control_clr"
    [1] "F2_cc_alpha_key"
    [1] "F2_clr_alpha_metadata"
    [1] "metadata_Q1_for_ML"
    [1] "metadata_Q2_for_ML"
    [1] "T1_all_species_metadata_var_imp_to_use"
    [1] "T1_metadata_alone"
    [1] "T2_all_species_metadata_var_imp_to_use"
    [1] "T2_metadata_alone"
    [1] "V1_case_control_clr"
    [1] "V1_cc_alpha_key"
    [1] "V1_clr_alpha_metadata"
    [1] "V1_F1_all_species_metadata"
    [1] "V1_F1_all_species"
    [1] "V1_F1_clr_alpha_metadata"
    [1] "V2_case_control_clr"
    [1] "V2_cc_alpha_key"
    [1] "V2_clr_alpha_metadata"
    [1] "V2_F2_all_species_metadata"
    [1] "V2_F2_all_species"
    [1] "V2_F2_clr_alpha_metadata"


### changing these to character, but might need to do it again, because it will go as integer again.


```R
metadata_Q2_for_ML[,c(1:38)] <- sapply(metadata_Q2_for_ML[,c(1:38)],function(x) as.character(x))

metadata_Q1_for_ML[,c(1:47)] <- sapply(metadata_Q1_for_ML[,c(1:47)],function(x) as.character(x))
```

### adding the _v and _f to the species so I know which they belong to once they're outputted by the algorithm.


```R

# Get the names of columns to be modified (excluding the first two columns)
columns_to_modify <- names(F2_case_control)[-c(1, 2)]

# Add "_f" to the selected column names
names(F2_case_control)[match(columns_to_modify, names(F2_case_control))] <- paste0(columns_to_modify, "_f")


# Get the names of columns to be modified (excluding the first two columns)
columns_to_modify <- names(F1_case_control)[-c(1, 2)]

# Add "_f" to the selected column names
names(F1_case_control)[match(columns_to_modify, names(F1_case_control))] <- paste0(columns_to_modify, "_f")

# Get the names of columns to be modified (excluding the first two columns)
columns_to_modify <- names(F1_F2_subtracted)[-c(1, 2)]

# Add "_f" to the selected column names
names(F1_F2_subtracted)[match(columns_to_modify, names(F1_F2_subtracted))] <- paste0(columns_to_modify, "_f")



# Get the names of columns to be modified (excluding the first two columns)
columns_to_modify <- names(V1_case_control)[-c(1, 2)]

# Add "_v" to the selected column names
names(V1_case_control)[match(columns_to_modify, names(V1_case_control))] <- paste0(columns_to_modify, "_v")

# Get the names of columns to be modified (excluding the first two columns)
columns_to_modify <- names(V2_case_control)[-c(1, 2)]

# Add "_v" to the selected column names
names(V2_case_control)[match(columns_to_modify, names(V2_case_control))] <- paste0(columns_to_modify, "_v")

# Get the names of columns to be modified (excluding the first two columns)
columns_to_modify <- names(V1_V2_subtracted)[-c(1, 2)]

# Add "_f" to the selected column names
names(V1_V2_subtracted)[match(columns_to_modify, names(V1_V2_subtracted))] <- paste0(columns_to_modify, "_f")


```

### combining the files. There are several combinations.

1) vaginal species + fecal species all T1


2) vaginal species + fecal species all + metadata Q1 T1

3) vaginal species + fecal species all T2


4) vaginal species + fecal species all + metadata Q2 T2



### Time point 1


```R
## vaginal species + fecal species all T1
V1_case_control_2 <- V1_case_control
F1_case_control_2 <- select(F1_case_control,-c("key"))
V1_F1_all_species <- merge(V1_case_control_2,F1_case_control_2,by="Studienummer") 
V1_F1_all_species <- na.omit(V1_F1_all_species)
V1_F1_all_species_2 <- select(V1_F1_all_species,-c("Studienummer","shannon_v","invsimp_v","richness_v","pielou_v","shannon_f","invsimp_f","richness_f","pielou_f"))
V1_case_control_2 <- select(V1_case_control,-c("Studienummer","shannon_v","invsimp_v","richness_v","pielou_v"))
V1_case_control_2 <- select(V1_case_control,-c("Studienummer","shannon_f","invsimp_f","richness_f","pielou_f"))

## vaginal species + fecal species all + metadata Q1 T1

V1_F1_all_species_metadata <- merge(V1_F1_all_species,metadata_Q1_for_ML,by="Studienummer") 
V1_F1_all_species_metadata <- na.omit(V1_F1_all_species_metadata)
V1_F1_all_species_metadata_2 <- select(V1_F1_all_species_metadata,-c("Studienummer","shannon_v","invsimp_v","pielou_v","richness_v","shannon_f","invsimp_f","richness_f","pielou_f"))

### metadata alone
T1_metadata_alone <- merge(V1_F1_all_species[,c(1:2)],metadata_Q1_for_ML,by="Studienummer")
T1_metadata_alone <- na.omit(T1_metadata_alone)
T1_metadata_alone_2 <- select(T1_metadata_alone,-c("Studienummer"))

```

### Time point 2


```R
## vaginal species + fecal species all T2
V2_case_control_2 <- V2_case_control
F2_case_control_2 <- select(F2_case_control,-c("key"))
V2_F2_all_species <- merge(V2_case_control_2,F2_case_control_2,by="Studienummer") 
V2_F2_all_species <- na.omit(V2_F2_all_species)
V2_F2_all_species_2 <- select(V2_F2_all_species,-c("Studienummer","shannon_v","invsimp_v","richness_v","pielou_v","shannon_f","invsimp_f","richness_f","pielou_f"))
V2_case_control_2 <- select(V2_case_control,-c("Studienummer","shannon_v","invsimp_v","richness_v","pielou_v"))
V2_case_control_2 <- select(V2_case_control,-c("Studienummer","shannon_f","invsimp_f","richness_f","pielou_f"))

## vaginal species + fecal species all + metadata Q2 T2

V2_F2_all_species_metadata <- merge(V2_F2_all_species,metadata_Q2_for_ML,by="Studienummer") 
V2_F2_all_species_metadata <- na.omit(V2_F2_all_species_metadata)
V2_F2_all_species_metadata_2 <- select(V2_F2_all_species_metadata,-c("Studienummer","shannon_v","invsimp_v","pielou_v","richness_v","shannon_f","invsimp_f","richness_f","pielou_f"))

### metadata alone
T2_metadata_alone <- merge(V2_F2_all_species[,c(1:2)],metadata_Q2_for_ML,by="Studienummer")
T2_metadata_alone <- na.omit(T2_metadata_alone)
T2_metadata_alone_2 <- Select(T2_metadata_alone,-c("Studienummer"))


### subtraction
T2_metadata_alone_2 <- select(T2_metadata_alone,-c("key"))
metadata_T1_T2 <- merge(T1_metadata_alone,T2_metadata_alone_2,by="Studienummer")
V1_V2_subtracted_2 <- select(V1_V2_subtracted,-c("key"))
F1_F2_subtracted_2 <- select(F1_F2_subtracted,-c("key"))
TP2_TP1_subtracted_species <- merge(V1_V2_subtracted_2,F1_F2_subtracted_2,by="Studienummer")
change_in_TP <- merge(metadata_T1_T2,TP2_TP1_subtracted_species,by="Studienummer")
change_in_TP_2 <- select(change_in_TP,-c("Studienummer"))

```

### Tuning


```R
### tuning methods
train.control_1 <- trainControl(method = "LOOCV", search ="random",summaryFunction=multiClassSummary,classProbs=T, savePredictions = T)
train.control_2 <- trainControl(method = "repeatedcv", repeats=5, search ="random",summaryFunction=multiClassSummary,classProbs=T, savePredictions = T)
train.control_3 <- trainControl(method = "boot", search ="random",summaryFunction=multiClassSummary,classProbs=T, savePredictions = T)

```

## function for ML

The analysis below was performed with all 3 training controls above. Alongside repeated with different methods (SVM, knn, nnet, rf) used within the prediction model. It was also done for two different pseudocounts to demonstrate that the addition of a pseudocount does not affect the results.

### Data extractions


```R

extract_rf_outputs <- function(
  model,
  prediction_probs,
  test_data,
  outcome_col = "key",
  positive_class = "Case"
) {
  
  # observed outcome
  obs <- factor(test_data[[outcome_col]])
  
  # predicted class from probabilities
  pred_class <- factor(
    ifelse(prediction_probs[[positive_class]] >= 0.5, positive_class, "Control"),
    levels = levels(obs)
  )
  
  # confusion matrix
  cm <- caret::confusionMatrix(
    data = pred_class,
    reference = obs,
    positive = positive_class
  )
  
  # ROC / AUC
  roc_obj <- pROC::roc(
    response = obs,
    predictor = prediction_probs[[positive_class]],
    levels = rev(levels(obs)),
    quiet = TRUE
  )
  
  auc_value <- as.numeric(pROC::auc(roc_obj))
  auc_ci <- as.numeric(pROC::ci.auc(roc_obj))
  
  # variable importance
  varimp <- caret::varImp(model)$importance %>%
    rownames_to_column("feature") %>%
    arrange(desc(Overall))
  
  # compact metrics
  metrics <- tibble(
    auc = auc_value,
    auc_ci_low = auc_ci[1],
    auc_ci_mid = auc_ci[2],
    auc_ci_high = auc_ci[3],
    sensitivity = cm$byClass["Sensitivity"],
    specificity = cm$byClass["Specificity"],
    ppv = cm$byClass["Pos Pred Value"],
    npv = cm$byClass["Neg Pred Value"],
    accuracy = cm$overall["Accuracy"],
    balanced_accuracy = cm$byClass["Balanced Accuracy"]
  )
  
  return(list(
    confusion_matrix = cm,
    roc = roc_obj,
    metrics = metrics,
    variable_importance = varimp
  ))
}
```

#### Ran twice with different pseudocount values


```R
run_rf_clr_iterations <- function(
  raw_df,
  outcome_col = "key",
  case_label = "Case",
  control_label = "Control",
  n_controls = 132,
  n_iter = 10,
  train_prop = 0.8,
  pseudocount = 1e-7,
  min_total_count = 100,
  output_dir,
  prefix = "T1_metadata_only",
  train_control,
  seed_range = 1000:5000,
  threshold = 0.5
) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  seed_list <- vector("list", n_iter)
  result_list <- vector("list", n_iter)
  metrics_list <- vector("list", n_iter)
  
  taxa_cols <- setdiff(colnames(raw_df), outcome_col)
  
  clr_transform <- function(count_df, pseudocount) {
    count_mat <- as.matrix(count_df)
    count_mat <- count_mat + pseudocount
    clr_mat <- compositions::clr(count_mat)
    clr_df <- as.data.frame(clr_mat)
    colnames(clr_df) <- colnames(count_df)
    clr_df
  }
  
  for (i in seq_len(n_iter)) {
    
    message("Running iteration: ", i)
    
    tryCatch({
      
      iter_seed <- sample(seed_range, 1)
      seed_list[[i]] <- iter_seed
      set.seed(iter_seed)
      
      control_subset <- raw_df %>%
        filter(.data[[outcome_col]] == control_label) %>%
        slice_sample(n = n_controls)
      
      case_set <- raw_df %>%
        filter(.data[[outcome_col]] == case_label)
      
      subsampled <- bind_rows(control_subset, case_set)
      
      indexes <- createDataPartition(
        subsampled[[outcome_col]],
        times = 1,
        p = train_prop,
        list = FALSE
      )
      
      train_raw <- subsampled[indexes, , drop = FALSE]
      test_raw  <- subsampled[-indexes, , drop = FALSE]
      
      train_counts <- train_raw[, taxa_cols, drop = FALSE]
      
      keep_features <- names(
        which(colSums(train_counts, na.rm = TRUE) > min_total_count)
      )
      
      train_counts_filtered <- train_raw[, keep_features, drop = FALSE]
      test_counts_filtered  <- test_raw[, keep_features, drop = FALSE]
      
      train_clr <- clr_transform(train_counts_filtered, pseudocount)
      test_clr  <- clr_transform(test_counts_filtered, pseudocount)
      
      train_dataset <- bind_cols(
        tibble(!!outcome_col := train_raw[[outcome_col]]),
        train_clr
      )
      
      test_dataset <- bind_cols(
        tibble(!!outcome_col := test_raw[[outcome_col]]),
        test_clr
      )
      
      train_dataset[[outcome_col]] <- factor(
        train_dataset[[outcome_col]],
        levels = c(control_label, case_label)
      )
      
      test_dataset[[outcome_col]] <- factor(
        test_dataset[[outcome_col]],
        levels = c(control_label, case_label)
      )
      
      rf_training <- train(
        as.formula(paste(outcome_col, "~ .")),
        data = train_dataset,
        method = "rf",
        tuneLength = 7,
        replace = FALSE,
        importance = TRUE,
        trControl = train_control
      )
      
      prediction_rfMod <- predict(
        rf_training,
        newdata = test_dataset,
        type = "prob"
      )
      
      extracted <- extract_rf_outputs(
        model = rf_training,
        prediction_probs = prediction_rfMod,
        test_data = test_dataset,
        outcome_col = outcome_col,
        positive_class = case_label,
        negative_class = control_label,
        threshold = threshold
      )
      
      metrics_list[[i]] <- extracted$metrics %>%
        mutate(
          iteration = i,
          seed = iter_seed,
          n_train = nrow(train_dataset),
          n_test = nrow(test_dataset),
          n_features_retained = length(keep_features),
          .before = 1
        )
      
      # File paths
      train_file <- file.path(output_dir, paste0(prefix, "_training_dataframe_", i, ".csv"))
      test_file <- file.path(output_dir, paste0(prefix, "_test_dataframe_", i, ".csv"))
      model_file <- file.path(output_dir, paste0(prefix, "_rf_model_", i, ".rds"))
      pred_file <- file.path(output_dir, paste0(prefix, "_predictions_", i, ".rds"))
      extracted_file <- file.path(output_dir, paste0(prefix, "_extracted_outputs_", i, ".rds"))
      varimp_file <- file.path(output_dir, paste0(prefix, "_variable_importance_", i, ".csv"))
      metrics_file <- file.path(output_dir, paste0(prefix, "_metrics_", i, ".csv"))
      features_file <- file.path(output_dir, paste0(prefix, "_features_retained_", i, ".rds"))
      
      write.csv(train_dataset, train_file, row.names = FALSE)
      write.csv(test_dataset, test_file, row.names = FALSE)
      saveRDS(rf_training, model_file)
      saveRDS(prediction_rfMod, pred_file)
      saveRDS(extracted, extracted_file)
      saveRDS(keep_features, features_file)
      
      write.csv(extracted$variable_importance, varimp_file, row.names = FALSE)
      write.csv(extracted$metrics, metrics_file, row.names = FALSE)
      
      result_list[[i]] <- list(
        iteration = i,
        seed = iter_seed,
        n_train = nrow(train_dataset),
        n_test = nrow(test_dataset),
        n_features_retained = length(keep_features),
        retained_features = keep_features,
        train_file = train_file,
        test_file = test_file,
        model_file = model_file,
        prediction_file = pred_file,
        extracted_file = extracted_file,
        varimp_file = varimp_file,
        metrics_file = metrics_file
      )
      
      message("Completed iteration: ", i)
      
    }, error = function(e) {
      message("Error in iteration ", i, ": ", conditionMessage(e))
      
      result_list[[i]] <- list(
        iteration = i,
        seed = seed_list[[i]],
        error = conditionMessage(e)
      )
    })
  }
  
  all_metrics <- bind_rows(metrics_list)
  
  seed_file <- file.path(output_dir, paste0(prefix, "_iteration_seeds.rds"))
  result_file <- file.path(output_dir, paste0(prefix, "_iteration_summary.rds"))
  all_metrics_file <- file.path(output_dir, paste0(prefix, "_all_metrics.csv"))
  
  saveRDS(seed_list, seed_file)
  saveRDS(result_list, result_file)
  write.csv(all_metrics, all_metrics_file, row.names = FALSE)
  
  list(
    seeds = seed_list,
    results = result_list,
    all_metrics = all_metrics,
    seed_file = seed_file,
    result_file = result_file,
    all_metrics_file = all_metrics_file
  )
}
```

### Metadata only


```R
rf_results_T1_metadata <- run_rf_clr_iterations(
  raw_df = T1_metadata_alone_2,
  outcome_col = "key",
  case_label = "Case",
  control_label = "Control",
  n_controls = 132,
  n_iter = 10,
  train_prop = 0.8,
  pseudocount = 1e-7,
  min_total_count = 100,
  output_dir = "/PATH/",
  prefix = "T1_metadata_alone_2",
  train_control = train.control_1
)

rf_results_T2_metadata <- run_rf_clr_iterations(
  raw_df = T2_metadata_alone_2,
  outcome_col = "key",
  case_label = "Case",
  control_label = "Control",
  n_controls = 83,
  n_iter = 10,
  train_prop = 0.8,
  pseudocount = 1e-7,
  min_total_count = 100,
  output_dir = "/PATH/",
  prefix = "T2_metadata_alone_2",
  train_control = train.control_1
)
```

#### Analyzing the output


```R
# Variable importance T1

varimp_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/T1_metadata_alone_2_extracted_outputs_", i, ".rds"))$variable_importance
})

varimp_all <- dplyr::bind_rows(varimp_list, .id = "iteration")

varimp_summary <- varimp_all %>%
  group_by(feature) %>%
  summarise(
    mean_importance = mean(relative_importance_score, na.rm = TRUE),
    sd_importance = sd(relative_importance_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_importance))


# CM, sensitivity, specificity



metrics_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/T1_metadata_alone_2_extracted_outputs_", i, ".rds"))$metrics
})

metrics_all <- dplyr::bind_rows(metrics_list, .id = "iteration")

metrics_summary_MD1 <- metrics_all %>%
  summarise(
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    sd_sensitivity = sd(sensitivity, na.rm = TRUE),
    
    mean_specificity = mean(specificity, na.rm = TRUE),
    sd_specificity = sd(specificity, na.rm = TRUE),
    
    mean_auc = mean(auc, na.rm = TRUE),
    sd_auc = sd(auc, na.rm = TRUE)
  )

cm_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/T1_metadata_alone_2_extracted_outputs_", i, ".rds"))$confusion_matrix
})

cm_tables <- lapply(cm_list, function(cm) {
  as.data.frame(cm$table)
})

cm_all <- dplyr::bind_rows(cm_tables, .id = "iteration")

cm_agg <- cm_all %>%
  group_by(Prediction, Reference) %>%
  summarise(n = sum(Freq), .groups = "drop")

tp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Case"]
tn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Control"]
fp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Control"]
fn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Case"]

sensitivity_MD1 <- tp / (tp + fn)
specificity_MD1 <- tn / (tn + fp)
```


```R
# Variable importance T2

varimp_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/T2_metadata_alone_2_extracted_outputs_", i, ".rds"))$variable_importance
})

varimp_all <- dplyr::bind_rows(varimp_list, .id = "iteration")

varimp_summary_MD2 <- varimp_all %>%
  group_by(feature) %>%
  summarise(
    mean_importance = mean(relative_importance_score, na.rm = TRUE),
    sd_importance = sd(relative_importance_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_importance))


# CM, sensitivity, specificity



metrics_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/T2_metadata_alone_2_extracted_outputs_", i, ".rds"))$metrics
})

metrics_all <- dplyr::bind_rows(metrics_list, .id = "iteration")

metrics_summary <- metrics_all %>%
  summarise(
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    sd_sensitivity = sd(sensitivity, na.rm = TRUE),
    
    mean_specificity = mean(specificity, na.rm = TRUE),
    sd_specificity = sd(specificity, na.rm = TRUE),
    
    mean_auc = mean(auc, na.rm = TRUE),
    sd_auc = sd(auc, na.rm = TRUE)
  )

cm_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/T2_metadata_alone_2_extracted_outputs_", i, ".rds"))$confusion_matrix
})

cm_tables <- lapply(cm_list, function(cm) {
  as.data.frame(cm$table)
})

cm_all <- dplyr::bind_rows(cm_tables, .id = "iteration")

cm_agg <- cm_all %>%
  group_by(Prediction, Reference) %>%
  summarise(n = sum(Freq), .groups = "drop")

tp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Case"]
tn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Control"]
fp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Control"]
fn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Case"]

sensitivity_MD2 <- tp / (tp + fn)
specificity_MD2 <- tn / (tn + fp)
```

### Vaginal T1 and T2


```R
#V1_case_control_2

rf_results_V1 <- run_rf_clr_iterations(
  raw_df = V1_case_control_2,
  outcome_col = "key",
  case_label = "Case",
  control_label = "Control",
  n_controls = 132,
  n_iter = 10,
  train_prop = 0.8,
  pseudocount = 1e-7,
  min_total_count = 100,
  output_dir = "/PATH/",
  prefix = "V1_case_control_2",
  train_control = train.control_1
)

rf_results_V2 <- run_rf_clr_iterations(
  raw_df = V2_case_control_2,
  outcome_col = "key",
  case_label = "Case",
  control_label = "Control",
  n_controls = 83,
  n_iter = 10,
  train_prop = 0.8,
  pseudocount = 1e-7,
  min_total_count = 100,
  output_dir = "/PATH/",
  prefix = "V2_case_control_2",
  train_control = train.control_1
)
```

#### analyzing the output


```R
# Variable importance V1

varimp_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/V1_case_control_2_extracted_outputs_", i, ".rds"))$variable_importance
})

varimp_all <- dplyr::bind_rows(varimp_list, .id = "iteration")

varimp_summary_V1 <- varimp_all %>%
  group_by(feature) %>%
  summarise(
    mean_importance = mean(relative_importance_score, na.rm = TRUE),
    sd_importance = sd(relative_importance_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_importance))


# CM, sensitivity, specificity



metrics_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/V1_case_control_2_extracted_outputs_", i, ".rds"))$metrics
})

metrics_all <- dplyr::bind_rows(metrics_list, .id = "iteration")

metrics_summary <- metrics_all %>%
  summarise(
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    sd_sensitivity = sd(sensitivity, na.rm = TRUE),
    
    mean_specificity = mean(specificity, na.rm = TRUE),
    sd_specificity = sd(specificity, na.rm = TRUE),
    
    mean_auc = mean(auc, na.rm = TRUE),
    sd_auc = sd(auc, na.rm = TRUE)
  )

cm_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/V1_case_control_2_extracted_outputs_", i, ".rds"))$confusion_matrix
})

cm_tables <- lapply(cm_list, function(cm) {
  as.data.frame(cm$table)
})

cm_all <- dplyr::bind_rows(cm_tables, .id = "iteration")

cm_agg <- cm_all %>%
  group_by(Prediction, Reference) %>%
  summarise(n = sum(Freq), .groups = "drop")

tp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Case"]
tn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Control"]
fp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Control"]
fn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Case"]

sensitivity_V1 <- tp / (tp + fn)
specificity_V1 <- tn / (tn + fp)
```


```R
# Variable importance V2

varimp_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/V2_case_control_2_extracted_outputs_", i, ".rds"))$variable_importance
})

varimp_all <- dplyr::bind_rows(varimp_list, .id = "iteration")

varimp_summary_V2 <- varimp_all %>%
  group_by(feature) %>%
  summarise(
    mean_importance = mean(relative_importance_score, na.rm = TRUE),
    sd_importance = sd(relative_importance_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_importance))


# CM, sensitivity, specificity



metrics_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/V2_case_control_2_extracted_outputs_", i, ".rds"))$metrics
})

metrics_all <- dplyr::bind_rows(metrics_list, .id = "iteration")

metrics_summary <- metrics_all %>%
  summarise(
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    sd_sensitivity = sd(sensitivity, na.rm = TRUE),
    
    mean_specificity = mean(specificity, na.rm = TRUE),
    sd_specificity = sd(specificity, na.rm = TRUE),
    
    mean_auc = mean(auc, na.rm = TRUE),
    sd_auc = sd(auc, na.rm = TRUE)
  )

cm_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/V2_case_control_2_extracted_outputs_", i, ".rds"))$confusion_matrix
})

cm_tables <- lapply(cm_list, function(cm) {
  as.data.frame(cm$table)
})

cm_all <- dplyr::bind_rows(cm_tables, .id = "iteration")

cm_agg <- cm_all %>%
  group_by(Prediction, Reference) %>%
  summarise(n = sum(Freq), .groups = "drop")

tp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Case"]
tn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Control"]
fp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Control"]
fn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Case"]

sensitivity_V2 <- tp / (tp + fn)
specificity_V2 <- tn / (tn + fp)
```

### Fecal T1 and T2


```R

rf_results_F1 <- run_rf_clr_iterations(
  raw_df = F1_case_control_2,
  outcome_col = "key",
  case_label = "Case",
  control_label = "Control",
  n_controls = 132,
  n_iter = 10,
  train_prop = 0.8,
  pseudocount = 1e-7,
  min_total_count = 100,
  output_dir = "/PATH/",
  prefix = "F1_case_control_2",
  train_control = train.control_1
)

rf_results_F2 <- run_rf_clr_iterations(
  raw_df = F2_case_control_2,
  outcome_col = "key",
  case_label = "Case",
  control_label = "Control",
  n_controls = 83,
  n_iter = 10,
  train_prop = 0.8,
  pseudocount = 1e-7,
  min_total_count = 100,
  output_dir = "/PATH/",
  prefix = "F2_case_control_2",
  train_control = train.control_1
)
```

#### Analyzing output


```R
# Variable importance F1

varimp_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/F1_case_control_2_extracted_outputs_", i, ".rds"))$variable_importance
})

varimp_all <- dplyr::bind_rows(varimp_list, .id = "iteration")

varimp_summary_F1 <- varimp_all %>%
  group_by(feature) %>%
  summarise(
    mean_importance = mean(relative_importance_score, na.rm = TRUE),
    sd_importance = sd(relative_importance_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_importance))


# CM, sensitivity, specificity



metrics_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/F1_case_control_2_extracted_outputs_", i, ".rds"))$metrics
})

metrics_all <- dplyr::bind_rows(metrics_list, .id = "iteration")

metrics_summary <- metrics_all %>%
  summarise(
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    sd_sensitivity = sd(sensitivity, na.rm = TRUE),
    
    mean_specificity = mean(specificity, na.rm = TRUE),
    sd_specificity = sd(specificity, na.rm = TRUE),
    
    mean_auc = mean(auc, na.rm = TRUE),
    sd_auc = sd(auc, na.rm = TRUE)
  )

cm_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/F1_case_control_2_extracted_outputs_", i, ".rds"))$confusion_matrix
})

cm_tables <- lapply(cm_list, function(cm) {
  as.data.frame(cm$table)
})

cm_all <- dplyr::bind_rows(cm_tables, .id = "iteration")

cm_agg <- cm_all %>%
  group_by(Prediction, Reference) %>%
  summarise(n = sum(Freq), .groups = "drop")

tp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Case"]
tn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Control"]
fp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Control"]
fn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Case"]

sensitivity_F1 <- tp / (tp + fn)
specificity_F1 <- tn / (tn + fp)
```


```R
# Variable importance F2

varimp_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/F2_case_control_2_extracted_outputs_", i, ".rds"))$variable_importance
})

varimp_all <- dplyr::bind_rows(varimp_list, .id = "iteration")

varimp_summary_F2 <- varimp_all %>%
  group_by(feature) %>%
  summarise(
    mean_importance = mean(relative_importance_score, na.rm = TRUE),
    sd_importance = sd(relative_importance_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_importance))


# CM, sensitivity, specificity



metrics_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/F2_case_control_2_extracted_outputs_", i, ".rds"))$metrics
})

metrics_all <- dplyr::bind_rows(metrics_list, .id = "iteration")

metrics_summary <- metrics_all %>%
  summarise(
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    sd_sensitivity = sd(sensitivity, na.rm = TRUE),
    
    mean_specificity = mean(specificity, na.rm = TRUE),
    sd_specificity = sd(specificity, na.rm = TRUE),
    
    mean_auc = mean(auc, na.rm = TRUE),
    sd_auc = sd(auc, na.rm = TRUE)
  )

cm_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/F2_case_control_2_extracted_outputs_", i, ".rds"))$confusion_matrix
})

cm_tables <- lapply(cm_list, function(cm) {
  as.data.frame(cm$table)
})

cm_all <- dplyr::bind_rows(cm_tables, .id = "iteration")

cm_agg <- cm_all %>%
  group_by(Prediction, Reference) %>%
  summarise(n = sum(Freq), .groups = "drop")

tp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Case"]
tn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Control"]
fp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Control"]
fn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Case"]

sensitivity_F2 <- tp / (tp + fn)
specificity_F2 <- tn / (tn + fp)
```

### Vagina and fecal T1 with and without metadata


```R

rf_results_V1_F1 <- run_rf_clr_iterations(
  raw_df = V1_F1_all_species_2,
  outcome_col = "key",
  case_label = "Case",
  control_label = "Control",
  n_controls = 132,
  n_iter = 10,
  train_prop = 0.8,
  pseudocount = 1e-7,
  min_total_count = 100,
  output_dir = "/PATH/",
  prefix = "V1_F1_all_species_2",
  train_control = train.control_1
)

rf_results_V1_F1_metadata <- run_rf_clr_iterations(
  raw_df = V1_F1_all_species_metadata_2,
  outcome_col = "key",
  case_label = "Case",
  control_label = "Control",
  n_controls = 83,
  n_iter = 10,
  train_prop = 0.8,
  pseudocount = 1e-7,
  min_total_count = 100,
  output_dir = "/PATH/",
  prefix = "V1_F1_all_species_metadata_2",
  train_control = train.control_1
)
```

#### Analyzing output


```R
# Variable importance V1 F1

varimp_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/V1_F1_all_species_2_extracted_outputs_", i, ".rds"))$variable_importance
})

varimp_all <- dplyr::bind_rows(varimp_list, .id = "iteration")

varimp_summary_V1_F1 <- varimp_all %>%
  group_by(feature) %>%
  summarise(
    mean_importance = mean(relative_importance_score, na.rm = TRUE),
    sd_importance = sd(relative_importance_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_importance))


# CM, sensitivity, specificity



metrics_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/V1_F1_all_species_2_extracted_outputs_", i, ".rds"))$metrics
})

metrics_all <- dplyr::bind_rows(metrics_list, .id = "iteration")

metrics_summary <- metrics_all %>%
  summarise(
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    sd_sensitivity = sd(sensitivity, na.rm = TRUE),
    
    mean_specificity = mean(specificity, na.rm = TRUE),
    sd_specificity = sd(specificity, na.rm = TRUE),
    
    mean_auc = mean(auc, na.rm = TRUE),
    sd_auc = sd(auc, na.rm = TRUE)
  )

cm_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/V1_F1_all_species_2_extracted_outputs_", i, ".rds"))$confusion_matrix
})

cm_tables <- lapply(cm_list, function(cm) {
  as.data.frame(cm$table)
})

cm_all <- dplyr::bind_rows(cm_tables, .id = "iteration")

cm_agg <- cm_all %>%
  group_by(Prediction, Reference) %>%
  summarise(n = sum(Freq), .groups = "drop")

tp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Case"]
tn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Control"]
fp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Control"]
fn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Case"]

sensitivity_V1_F1 <- tp / (tp + fn)
specificity_V1_F1 <- tn / (tn + fp)
```


```R
# Variable importance V1 F1 and metadata

varimp_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/V1_F1_all_species_metadata_2_extracted_outputs_", i, ".rds"))$variable_importance
})

varimp_all <- dplyr::bind_rows(varimp_list, .id = "iteration")

varimp_summary_V1_F1_metadata <- varimp_all %>%
  group_by(feature) %>%
  summarise(
    mean_importance = mean(relative_importance_score, na.rm = TRUE),
    sd_importance = sd(relative_importance_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_importance))


# CM, sensitivity, specificity



metrics_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/V1_F1_all_species_metadata_2_extracted_outputs_", i, ".rds"))$metrics
})

metrics_all <- dplyr::bind_rows(metrics_list, .id = "iteration")

metrics_summary <- metrics_all %>%
  summarise(
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    sd_sensitivity = sd(sensitivity, na.rm = TRUE),
    
    mean_specificity = mean(specificity, na.rm = TRUE),
    sd_specificity = sd(specificity, na.rm = TRUE),
    
    mean_auc = mean(auc, na.rm = TRUE),
    sd_auc = sd(auc, na.rm = TRUE)
  )

cm_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/V1_F1_all_species_metadata_2_extracted_outputs_", i, ".rds"))$confusion_matrix
})

cm_tables <- lapply(cm_list, function(cm) {
  as.data.frame(cm$table)
})

cm_all <- dplyr::bind_rows(cm_tables, .id = "iteration")

cm_agg <- cm_all %>%
  group_by(Prediction, Reference) %>%
  summarise(n = sum(Freq), .groups = "drop")

tp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Case"]
tn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Control"]
fp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Control"]
fn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Case"]

sensitivity_V1_F1_metadata <- tp / (tp + fn)
specificity_V1_F1_metadata <- tn / (tn + fp)
```

### Vaginal and fecal T2 with and without metadata


```R

rf_results_V2_F2 <- run_rf_clr_iterations(
  raw_df = V2_F2_all_species_2,
  outcome_col = "key",
  case_label = "Case",
  control_label = "Control",
  n_controls = 132,
  n_iter = 10,
  train_prop = 0.8,
  pseudocount = 1e-7,
  min_total_count = 100,
  output_dir = "/PATH/",
  prefix = "V2_F2_all_species_2",
  train_control = train.control_1
)

rf_results_V2_F2_metadata <- run_rf_clr_iterations(
  raw_df = V2_F2_all_species_metadata_2,
  outcome_col = "key",
  case_label = "Case",
  control_label = "Control",
  n_controls = 83,
  n_iter = 10,
  train_prop = 0.8,
  pseudocount = 1e-7,
  min_total_count = 100,
  output_dir = "/PATH/",
  prefix = "V2_F2_all_species_metadata_2",
  train_control = train.control_1
)

```

#### Analyzing output


```R
# Variable importance V2 F2

varimp_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/V2_F2_all_species_2_extracted_outputs_", i, ".rds"))$variable_importance
})

varimp_all <- dplyr::bind_rows(varimp_list, .id = "iteration")

varimp_summary_V2_F2 <- varimp_all %>%
  group_by(feature) %>%
  summarise(
    mean_importance = mean(relative_importance_score, na.rm = TRUE),
    sd_importance = sd(relative_importance_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_importance))


# CM, sensitivity, specificity



metrics_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/V2_F2_all_species_2_extracted_outputs_", i, ".rds"))$metrics
})

metrics_all <- dplyr::bind_rows(metrics_list, .id = "iteration")

metrics_summary <- metrics_all %>%
  summarise(
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    sd_sensitivity = sd(sensitivity, na.rm = TRUE),
    
    mean_specificity = mean(specificity, na.rm = TRUE),
    sd_specificity = sd(specificity, na.rm = TRUE),
    
    mean_auc = mean(auc, na.rm = TRUE),
    sd_auc = sd(auc, na.rm = TRUE)
  )

cm_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/V2_F2_all_species_2_extracted_outputs_", i, ".rds"))$confusion_matrix
})

cm_tables <- lapply(cm_list, function(cm) {
  as.data.frame(cm$table)
})

cm_all <- dplyr::bind_rows(cm_tables, .id = "iteration")

cm_agg <- cm_all %>%
  group_by(Prediction, Reference) %>%
  summarise(n = sum(Freq), .groups = "drop")

tp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Case"]
tn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Control"]
fp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Control"]
fn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Case"]

sensitivity_V2_F2 <- tp / (tp + fn)
specificity_V2_F2 <- tn / (tn + fp)
```


```R
# Variable importance V2 F2 and metadata

varimp_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/V2_F2_all_species_metadata_2_extracted_outputs_", i, ".rds"))$variable_importance
})

varimp_all <- dplyr::bind_rows(varimp_list, .id = "iteration")

varimp_summary_V2_F2_metadata <- varimp_all %>%
  group_by(feature) %>%
  summarise(
    mean_importance = mean(relative_importance_score, na.rm = TRUE),
    sd_importance = sd(relative_importance_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_importance))


# CM, sensitivity, specificity



metrics_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/V2_F2_all_species_metadata_2_extracted_outputs_", i, ".rds"))$metrics
})

metrics_all <- dplyr::bind_rows(metrics_list, .id = "iteration")

metrics_summary <- metrics_all %>%
  summarise(
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    sd_sensitivity = sd(sensitivity, na.rm = TRUE),
    
    mean_specificity = mean(specificity, na.rm = TRUE),
    sd_specificity = sd(specificity, na.rm = TRUE),
    
    mean_auc = mean(auc, na.rm = TRUE),
    sd_auc = sd(auc, na.rm = TRUE)
  )

cm_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/V2_F2_all_species_metadata_2_extracted_outputs_", i, ".rds"))$confusion_matrix
})

cm_tables <- lapply(cm_list, function(cm) {
  as.data.frame(cm$table)
})

cm_all <- dplyr::bind_rows(cm_tables, .id = "iteration")

cm_agg <- cm_all %>%
  group_by(Prediction, Reference) %>%
  summarise(n = sum(Freq), .groups = "drop")

tp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Case"]
tn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Control"]
fp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Control"]
fn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Case"]

sensitivity_V2_F2_metadata <- tp / (tp + fn)
specificity_V2_F2_metadata <- tn / (tn + fp)
```

### VarImp and alpha diversity
We found that the output from the final run -- V1/2 + F1/2 + metadata 1/2 was the best and therefore we continued the analysis using variable importance selection using the output from the above (V1/2 + F1/2 + metadata 1/2 for the remainder of the run.


```R
#varimp_summary_V1_F1_metadata
# V1_F1_all_species_metadata

# V1_F1_all_species_metadata (this df has the IDs so we can extract alpha diversity from it)

# timepoint 1
varimp_V1_F1_metadata <- V1_F1_all_species_metadata[,varimp_summary_V1_F1_metadata$feature]
varimp_V1_F1_metadata$key <- V1_F1_all_species_metadata$key
# with alpha
varimp_V1_F1_metadata_alpha <- V1_F1_all_species_metadata[,varimp_summary_V1_F1_metadata$feature]
varimp_V1_F1_metadata_alpha$key <- V1_F1_all_species_metadata$key
varimp_V1_F1_metadata_alpha$shannon_v <- V1_F1_all_species_metadata$shannon_v
varimp_V1_F1_metadata_alpha$invsimp_v <- V1_F1_all_species_metadata$invsimp_v
varimp_V1_F1_metadata_alpha$pielou_v <- V1_F1_all_species_metadata$pielou_v
varimp_V1_F1_metadata_alpha$richness_v <- V1_F1_all_species_metadata$richness_v
varimp_V1_F1_metadata_alpha$shannon_f <- V1_F1_all_species_metadata$shannon_f
varimp_V1_F1_metadata_alpha$invsimp_f <- V1_F1_all_species_metadata$invsimp_f
varimp_V1_F1_metadata_alpha$pielou_f <- V1_F1_all_species_metadata$pielou_f
varimp_V1_F1_metadata_alpha$richness_f <- V1_F1_all_species_metadata$richness_f

varimp_V1_F1_metadata_2 <- varimp_V1_F1_metadata[,-c("Studienummer")]
varimp_V1_F1_metadata_alpha_2 <- varimp_V1_F1_metadata_alpha[,-c("Studienummer")]

# variables
fecal_var_imp_T1 <- varimp_summary_V1_F1_metadata %>% filter(grepl("_f", feature)) #most important fecal species
vaginal_var_imp_T1 <- varimp_summary_V1_F1_metadata %>% filter(grepl("_v", feature)) #most important vaginal species
metadata_var_imp_T1 <- varimp_summary_V1_F1_metadata %>% filter(!grepl("k__", feature))#most important metadata

# Timepoint 2
varimp_V2_F2_metadata <- V2_F2_all_species_metadata_2[,varimp_summary_V2_F2_metadata$feature]
varimp_V2_F2_metadata$key <- V2_F2_all_species_metadata_2$key
# with alpha
varimp_V2_F2_metadata_alpha <- V2_F2_all_species_metadata_2[,varimp_summary_V2_F2_metadata$feature]
varimp_V2_F2_metadata_alpha$key <- V2_F2_all_species_metadata_2$key
varimp_V2_F2_metadata_alpha$shannon_v <- V2_F2_all_species_metadata$shannon_v
varimp_V2_F2_metadata_alpha$invsimp_v <- V2_F2_all_species_metadata$invsimp_v
varimp_V2_F2_metadata_alpha$pielou_v <- V2_F2_all_species_metadata$pielou_v
varimp_V2_F2_metadata_alpha$richness_v <- V2_F2_all_species_metadata$richness_v
varimp_V2_F2_metadata_alpha$shannon_f <- V2_F2_all_species_metadata$shannon_f
varimp_V2_F2_metadata_alpha$invsimp_f <- V2_F2_all_species_metadata$invsimp_f
varimp_V2_F2_metadata_alpha$pielou_f <- V2_F2_all_species_metadata$pielou_f
varimp_V2_F2_metadata_alpha$richness_f <- V2_F2_all_species_metadata$richness_f

varimp_V2_F2_metadata_2 <- varimp_V2_F2_metadata[,-c("Studienummer")]
varimp_V2_F2_metadata_alpha_2 <- varimp_V2_F2_metadata_alpha[,-c("Studienummer")]

#varibles
fecal_var_imp_T2 <- varimp_summary_V2_F2_metadata %>% filter(grepl("_f", feature)) #most important fecal species
vaginal_var_imp_T2 <- varimp_summary_V2_F2_metadata %>% filter(grepl("_v", feature)) #most important vaginal species
metadata_var_imp_T2 <- varimp_summary_V2_F2_metadata %>% filter(!grepl("k__", feature))#most important metadata

```

### Varimp without alpha TP1 and TP2


```R
rf_results_Varimp_TP1 <- run_rf_clr_iterations(
  raw_df = varimp_V1_F1_metadata_2,
  outcome_col = "key",
  case_label = "Case",
  control_label = "Control",
  n_controls = 132,
  n_iter = 10,
  train_prop = 0.8,
  pseudocount = 1e-7,
  min_total_count = 100,
  output_dir = "/PATH/",
  prefix = "results_Varimp_TP1",
  train_control = train.control_1
)

rf_results_Varimp_TP2 <- run_rf_clr_iterations(
  raw_df = varimp_V2_F2_metadata_2,
  outcome_col = "key",
  case_label = "Case",
  control_label = "Control",
  n_controls = 132,
  n_iter = 10,
  train_prop = 0.8,
  pseudocount = 1e-7,
  min_total_count = 100,
  output_dir = "/PATH/",
  prefix = "results_Varimp_TP2",
  train_control = train.control_1
)

```


```R
# CM, sensitivity, specificity



metrics_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/results_Varimp_TP1_extracted_outputs_", i, ".rds"))$metrics
})

metrics_all <- dplyr::bind_rows(metrics_list, .id = "iteration")

metrics_summary <- metrics_all %>%
  summarise(
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    sd_sensitivity = sd(sensitivity, na.rm = TRUE),
    
    mean_specificity = mean(specificity, na.rm = TRUE),
    sd_specificity = sd(specificity, na.rm = TRUE),
    
    mean_auc = mean(auc, na.rm = TRUE),
    sd_auc = sd(auc, na.rm = TRUE)
  )

cm_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/results_Varimp_TP1_extracted_outputs_", i, ".rds"))$confusion_matrix
})

cm_tables <- lapply(cm_list, function(cm) {
  as.data.frame(cm$table)
})

cm_all <- dplyr::bind_rows(cm_tables, .id = "iteration")

cm_agg <- cm_all %>%
  group_by(Prediction, Reference) %>%
  summarise(n = sum(Freq), .groups = "drop")

tp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Case"]
tn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Control"]
fp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Control"]
fn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Case"]

sensitivity_varimp_TP1 <- tp / (tp + fn)
specificity_varimp_TP1 <- tn / (tn + fp)
```


```R
# CM, sensitivity, specificity

metrics_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/results_Varimp_TP2_extracted_outputs_", i, ".rds"))$metrics
})

metrics_all <- dplyr::bind_rows(metrics_list, .id = "iteration")

metrics_summary <- metrics_all %>%
  summarise(
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    sd_sensitivity = sd(sensitivity, na.rm = TRUE),
    
    mean_specificity = mean(specificity, na.rm = TRUE),
    sd_specificity = sd(specificity, na.rm = TRUE),
    
    mean_auc = mean(auc, na.rm = TRUE),
    sd_auc = sd(auc, na.rm = TRUE)
  )

cm_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/results_Varimp_TP2_extracted_outputs_", i, ".rds"))$confusion_matrix
})

cm_tables <- lapply(cm_list, function(cm) {
  as.data.frame(cm$table)
})

cm_all <- dplyr::bind_rows(cm_tables, .id = "iteration")

cm_agg <- cm_all %>%
  group_by(Prediction, Reference) %>%
  summarise(n = sum(Freq), .groups = "drop")

tp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Case"]
tn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Control"]
fp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Control"]
fn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Case"]

sensitivity_varimp_TP2 <- tp / (tp + fn)
specificity_varimp_TP2 <- tn / (tn + fp)
```

### Varimp with alpha TP1 and TP2


```R
rf_results_Varimp_alpha <- run_rf_clr_iterations(
  raw_df = varimp_V1_F1_metadata_alpha_2,
  outcome_col = "key",
  case_label = "Case",
  control_label = "Control",
  n_controls = 132,
  n_iter = 10,
  train_prop = 0.8,
  pseudocount = 1e-7,
  min_total_count = 100,
  output_dir = "/PATH/",
  prefix = "results_Varimp_alpha_TP1",
  train_control = train.control_1
)

rf_results_Varimp_alpha <- run_rf_clr_iterations(
  raw_df = varimp_V2_F2_metadata_alpha_2,
  outcome_col = "key",
  case_label = "Case",
  control_label = "Control",
  n_controls = 132,
  n_iter = 10,
  train_prop = 0.8,
  pseudocount = 1e-7,
  min_total_count = 100,
  output_dir = "/PATH/",
  prefix = "results_Varimp_alpha_TP2",
  train_control = train.control_1
)
```

#### Analyzing output


```R
# CM, sensitivity, specificity varimp TP1 with alpha



metrics_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/results_Varimp_alpha_TP1_extracted_outputs_", i, ".rds"))$metrics
})

metrics_all <- dplyr::bind_rows(metrics_list, .id = "iteration")

metrics_summary <- metrics_all %>%
  summarise(
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    sd_sensitivity = sd(sensitivity, na.rm = TRUE),
    
    mean_specificity = mean(specificity, na.rm = TRUE),
    sd_specificity = sd(specificity, na.rm = TRUE),
    
    mean_auc = mean(auc, na.rm = TRUE),
    sd_auc = sd(auc, na.rm = TRUE)
  )

cm_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/results_Varimp_alpha_TP1_extracted_outputs_", i, ".rds"))$confusion_matrix
})

cm_tables <- lapply(cm_list, function(cm) {
  as.data.frame(cm$table)
})

cm_all <- dplyr::bind_rows(cm_tables, .id = "iteration")

cm_agg <- cm_all %>%
  group_by(Prediction, Reference) %>%
  summarise(n = sum(Freq), .groups = "drop")

tp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Case"]
tn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Control"]
fp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Control"]
fn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Case"]

sensitivity_varimp_alpha_TP1 <- tp / (tp + fn)
specificity_varimp_alpha_TP1 <- tn / (tn + fp)
```


```R
# CM, sensitivity, specificity varimp TP2 with alpha



metrics_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/results_Varimp_alpha_TP2_extracted_outputs_", i, ".rds"))$metrics
})

metrics_all <- dplyr::bind_rows(metrics_list, .id = "iteration")

metrics_summary <- metrics_all %>%
  summarise(
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    sd_sensitivity = sd(sensitivity, na.rm = TRUE),
    
    mean_specificity = mean(specificity, na.rm = TRUE),
    sd_specificity = sd(specificity, na.rm = TRUE),
    
    mean_auc = mean(auc, na.rm = TRUE),
    sd_auc = sd(auc, na.rm = TRUE)
  )

cm_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/results_Varimp_alpha_TP22_extracted_outputs_", i, ".rds"))$confusion_matrix
})

cm_tables <- lapply(cm_list, function(cm) {
  as.data.frame(cm$table)
})

cm_all <- dplyr::bind_rows(cm_tables, .id = "iteration")

cm_agg <- cm_all %>%
  group_by(Prediction, Reference) %>%
  summarise(n = sum(Freq), .groups = "drop")

tp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Case"]
tn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Control"]
fp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Control"]
fn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Case"]

sensitivity_varimp_alpha_TP2 <- tp / (tp + fn)
specificity_varimp_alpha_TP2 <- tn / (tn + fp)
```

### Change in timepoint species alongside metadata from both timepoints


```R
# change_in_TP_2
rf_results_Varimp_alpha <- run_rf_clr_iterations(
  raw_df = change_in_TP_2,
  outcome_col = "key",
  case_label = "Case",
  control_label = "Control",
  n_controls = 83,
  n_iter = 10,
  train_prop = 0.8,
  pseudocount = 1e-7,
  min_total_count = 100,
  output_dir = "/PATH/",
  prefix = "change_in_TP_2",
  train_control = train.control_1
)
```


```R
# CM, sensitivity, specificity change in TP and metadata



metrics_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/change_in_TP_2_extracted_outputs_", i, ".rds"))$metrics
})

metrics_all <- dplyr::bind_rows(metrics_list, .id = "iteration")

metrics_summary <- metrics_all %>%
  summarise(
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    sd_sensitivity = sd(sensitivity, na.rm = TRUE),
    
    mean_specificity = mean(specificity, na.rm = TRUE),
    sd_specificity = sd(specificity, na.rm = TRUE),
    
    mean_auc = mean(auc, na.rm = TRUE),
    sd_auc = sd(auc, na.rm = TRUE)
  )

cm_list <- lapply(1:10, function(i) {
  readRDS(paste0("/PATH/change_in_TP_2_extracted_outputs_", i, ".rds"))$confusion_matrix
})

cm_tables <- lapply(cm_list, function(cm) {
  as.data.frame(cm$table)
})

cm_all <- dplyr::bind_rows(cm_tables, .id = "iteration")

cm_agg <- cm_all %>%
  group_by(Prediction, Reference) %>%
  summarise(n = sum(Freq), .groups = "drop")

tp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Case"]
tn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Control"]
fp <- cm_agg$n[cm_agg$Prediction == "Case" & cm_agg$Reference == "Control"]
fn <- cm_agg$n[cm_agg$Prediction == "Control" & cm_agg$Reference == "Case"]

sensitivity_change_MD <- tp / (tp + fn)
specificity_change_MD <- tn / (tn + fp)
```

### plots -- varimp 


```R
#fecal_var_imp_T1

fecal_T1_varimp <- fecal_var_imp_T1[c(1:10),]
colnames(fecal_T1_varimp) <- c("vars","value")

F1_varimp_plot <- ggplot(fecal_T1_varimp, aes(x = value, y = reorder(vars, -value))) +
  geom_segment(aes(xend = 0, yend = reorder(vars, -value)), color = "grey") +  # Draw lines from x=0 to the value
  geom_point(color = "blue", size = 4) +  # Add the dots
  theme_minimal() +  # Clean theme
  labs(x = "VarImp scaled importance", y = "fecal species", title = "fecal species with top VIP variables TP1") +
  theme(panel.grid.major.y = element_blank(),  # Remove y gridlines
        panel.grid.minor = element_blank())    # Remove minor gridlines


#ggsave("/PATH/F1_varimp_plot.pdf", plot = F1_varimp_plot, width = 12, height = 12, dpi = 600)

#vaginal_var_imp_T1


vaginal_T1_varimp <- vaginal_var_imp_T1[c(1:10),]
colnames(vaginal_T1_varimp) <- c("vars","value")


V1_varimp_plot <- ggplot(vaginal_T1_varimp, aes(x = value, y = reorder(vars, -value))) +
  geom_segment(aes(xend = 0, yend = reorder(vars, -value)), color = "grey") +  # Draw lines from x=0 to the value
  geom_point(color = "blue", size = 4) +  # Add the dots
  theme_minimal() +  # Clean theme
  labs(x = "VarImp scaled importance", y = "vaginal species", title = "vaginal species with top VIP variables TP1") +
  theme(panel.grid.major.y = element_blank(),  # Remove y gridlines
        panel.grid.minor = element_blank())    # Remove minor gridlines

#ggsave("/PATH/V1_varimp_plot.pdf", plot = V1_varimp_plot, width = 12, height = 12, dpi = 600)

#metadata_T1_varimp
metadata_T1_varimp <- metadata_var_imp_T1[c(1:10),]
colnames(metadata_T1_varimp) <- c("vars","value")


metadata_T1_varimp_plot <- ggplot(metadata_T1_varimp, aes(x = value, y = reorder(vars, -value))) +
  geom_segment(aes(xend = 0, yend = reorder(vars, -value)), color = "grey") +  # Draw lines from x=0 to the value
  geom_point(color = "blue", size = 4) +  # Add the dots
  theme_minimal() +  # Clean theme
  labs(x = "VarImp scaled importance", y = "metadata", title = "metadata with top VIP variables TP1") +
  theme(panel.grid.major.y = element_blank(),  # Remove y gridlines
        panel.grid.minor = element_blank())    # Remove minor gridlines

#ggsave("/PATH/metadata_T1_varimp_plot.pdf", plot = metadata_T1_varimp_plot, width = 12, height = 12, dpi = 600)

```


```R
#fecal_var_imp_T2

fecal_T2_varimp <- fecal_var_imp_T2[c(1:10),]
colnames(fecal_T2_varimp) <- c("vars","value")

F2_varimp_plot <- ggplot(fecal_T2_varimp, aes(x = value, y = reorder(vars, -value))) +
  geom_segment(aes(xend = 0, yend = reorder(vars, -value)), color = "grey") +  # Draw lines from x=0 to the value
  geom_point(color = "blue", size = 4) +  # Add the dots
  theme_minimal() +  # Clean theme
  labs(x = "VarImp scaled importance", y = "fecal species", title = "fecal species with top VIP variables TP1") +
  theme(panel.grid.major.y = element_blank(),  # Remove y gridlines
        panel.grid.minor = element_blank())    # Remove minor gridlines


#ggsave("/PATH/F2_varimp_plot.pdf", plot = F2_varimp_plot, width = 12, height = 12, dpi = 600)

#vaginal_var_imp_T2


vaginal_T2_varimp <- vaginal_var_imp_T2[c(1:10),]
colnames(vaginal_T2_varimp) <- c("vars","value")


V2_varimp_plot <- ggplot(vaginal_T2_varimp, aes(x = value, y = reorder(vars, -value))) +
  geom_segment(aes(xend = 0, yend = reorder(vars, -value)), color = "grey") +  # Draw lines from x=0 to the value
  geom_point(color = "blue", size = 4) +  # Add the dots
  theme_minimal() +  # Clean theme
  labs(x = "VarImp scaled importance", y = "vaginal species", title = "vaginal species with top VIP variables TP1") +
  theme(panel.grid.major.y = element_blank(),  # Remove y gridlines
        panel.grid.minor = element_blank())    # Remove minor gridlines

#ggsave("/PATH/V2_varimp_plot.pdf", plot = V2_varimp_plot, width = 12, height = 12, dpi = 600)

#metadata_var_imp_T2
metadata_T2_varimp <- metadata_var_imp_T2[c(1:10),]
colnames(metadata_T2_varimp) <- c("vars","value")


metadata_T2_varimp_plot <- ggplot(metadata_T2_varimp, aes(x = value, y = reorder(vars, -value))) +
  geom_segment(aes(xend = 0, yend = reorder(vars, -value)), color = "grey") +  # Draw lines from x=0 to the value
  geom_point(color = "blue", size = 4) +  # Add the dots
  theme_minimal() +  # Clean theme
  labs(x = "VarImp scaled importance", y = "metadata", title = "metadata with top VIP variables TP1") +
  theme(panel.grid.major.y = element_blank(),  # Remove y gridlines
        panel.grid.minor = element_blank())    # Remove minor gridlines

#ggsave("/PATH/metadata_T2_varimp_plot.pdf", plot = metadata_T2_varimp_plot, width = 12, height = 12, dpi = 600)

```

### Plots sens/spec

##### Not denoised (original) TP1


```R
x_1 <- 1-specificity_V1_F1_metadata
y_1 <- sensitivity_V1_F1_metadata

x_2 <- 1-specificity_MD1
y_2 <- sensitivity_MD1


x_4 <- 1-specificity_V1
y_4 <- sensitivity_V1


x_3 <- 1-specificity_F1
y_3 <- sensitivity_F1

#### making figure

Figure <- ggplot() +
  geom_line(aes(x_1, y_1, color = "1. TP1 all taxa and metadata 0.77"), linewidth = 1) +
  geom_line(aes(x_2, y_2, color = "2. TP1 metadata only 0.74"), linewidth = 1) +
  geom_line(aes(x_4, y_4, color = "4. TP1 vaginal only 0.51"), linewidth = 1) +
  geom_line(aes(x_3, y_3, color = "3. TP1 fecal only 0.59"), linewidth = 1) +
  labs(title = "Initial AUC-ROC curves TP1", x = "False positive rate", y = "True positive rate") +
  scale_color_manual(values = c("1. TP1 all taxa and metadata 0.77" = "purple", "2. TP1 metadata only 0.74" = "red", "4. TP1 vaginal only 0.51" = "blue","3. TP1 fecal only 0.59" = "green")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  theme(axis.text.y = element_text(size = 18),
    axis.text.x = element_text(angle = 90, hjust = 1,size=18),# Adjust axis label size
    axis.title = element_text(size = 18),         # Adjust axis title size
    plot.title = element_text(size = 20, face = "bold"),  # Adjust plot title size
    panel.background = element_rect(fill = "white"),
       legend.text = element_text(size = 12))

#ggsave("/PATH/ML_TP1_1.pdf", plot = Figure, width = 12, height = 12, dpi = 600)

```

##### Denoised (varimp) TP1


```R
x_1 <- 1-specificity_varimp_alpha_TP1
y_1 <- sensitivity_varimp_alpha_TP1

x_2 <- 1-specificity_varimp_TP1
y_2 <- sensitivity_varimp_TP1

Figure <- ggplot() +
  geom_line(aes(x_1, y_1, color = "1. TP1 denoised with alpha 0.89"), linewidth = 1) +
  geom_line(aes(x_2, y_2, color = "2. TP1 denoised without alpha 0.86"), linewidth = 1) +

  labs(title = "Denoised AUC-ROC curves TP1", x = "False positive rate", y = "True positive rate") +
  scale_color_manual(values = c("1. TP1 denoised with alpha 0.89" = "purple", "2. TP1 denoised without alpha 0.86" = "red")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  theme(axis.text.y = element_text(size = 18),
    axis.text.x = element_text(angle = 90, hjust = 1,size=18),# Adjust axis label size
    axis.title = element_text(size = 18),         # Adjust axis title size
    plot.title = element_text(size = 20, face = "bold"),  # Adjust plot title size
    panel.background = element_rect(fill = "white"),
       legend.text = element_text(size = 12))

#ggsave("/PATH/denoised_ML_TP1_1.pdf", plot = Figure, width = 12, height = 12, dpi = 600)


```

##### Not denoised signal TP2


```R
x_1 <- 1-specificity_V2_F2_metadata
y_1 <- sensitivity_V2_F2_metadata

x_2 <- 1-specificity_MD2
y_2 <- sensitivity_MD2


x_4 <- 1-specificity_V2
y_4 <- sensitivity_V2


x_3 <- 1-specificity_F2
y_3 <- sensitivity_F2

#### making figure

Figure <- ggplot() +
  geom_line(aes(x_1, y_1, color = "1. TP2 all taxa and metadata 0.7"), linewidth = 1) +
  geom_line(aes(x_2, y_2, color = "2. TP2 metadata only 0.68"), linewidth = 1) +
  geom_line(aes(x_4, y_4, color = "4. TP2 vaginal only 0.54"), linewidth = 1) +
  geom_line(aes(x_3, y_3, color = "3. TP2 fecal only 0.59"), linewidth = 1) +
  labs(title = "Initial AUC-ROC curves TP1", x = "False positive rate", y = "True positive rate") +
  scale_color_manual(values = c("1. TP2 all taxa and metadata 0.7" = "red", "2. TP2 metadata only 0.68" = "orange", "4. TP2 vaginal only 0.54" = "blue","3. TP2 fecal only 0.59" = "green")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  theme(axis.text.y = element_text(size = 18),
    axis.text.x = element_text(angle = 90, hjust = 1,size=18),# Adjust axis label size
    axis.title = element_text(size = 18),         # Adjust axis title size
    plot.title = element_text(size = 20, face = "bold"),  # Adjust plot title size
    panel.background = element_rect(fill = "white"),
       legend.text = element_text(size = 12))

#ggsave("/PATH/ML_TP2_1.pdf", plot = Figure, width = 12, height = 12, dpi = 600)

```

##### Denoised TP2 and change between TPs


```R
x_1 <- 1-specificity_varimp_alpha_TP2
y_1 <- sensitivity_varimp_alpha_TP2

x_2 <- 1-specificity_varimp_TP2
y_2 <- sensitivity_varimp_TP2

x_3 <- 1-specificity_change_MD
y_3 <- sensitivity_change_MD

Figure <- ggplot() +
  geom_line(aes(x_1, y_1, color = "1. TP2 denoised with alpha 0.77"), linewidth = 1) +
  geom_line(aes(x_2, y_2, color = "2. TP2 denoised without alpha 0.75"), linewidth = 1) +
  geom_line(aes(x_3, y_3, color = "3. Change between TPs 0.82"), linewidth = 1) +
  labs(title = "Denoised AUC-ROC curves TP1", x = "False positive rate", y = "True positive rate") +
  scale_color_manual(values = c("1. TP2 denoised with alpha 0.77" = "red", "2. TP2 denoised without alpha 0.75" = "purple","3. Change between TPs 0.82" = "blue")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  theme(axis.text.y = element_text(size = 18),
    axis.text.x = element_text(angle = 90, hjust = 1,size=18),# Adjust axis label size
    axis.title = element_text(size = 18),         # Adjust axis title size
    plot.title = element_text(size = 20, face = "bold"),  # Adjust plot title size
    panel.background = element_rect(fill = "white"),
       legend.text = element_text(size = 12))

#ggsave("/PATH/denoised_ML_TP2_1.pdf", plot = Figure, width = 12, height = 12, dpi = 600)


```
