# ============================================================
#  2-layer MLP + Feature Attention (Keras)
#
# ============================================================

# -------------------- Libraries --------------------
library(plyr)
library(boot)
library(keras)
library(tfruns)
library(MLmetrics)
library(parallel)
library(snow)
library(reticulate)
library(tensorflow)
library(ggplot2)

set.seed(123)
tensorflow::tf$random$set_seed(123)

# -------------------- Inputs --------------------
file_path <- "/projects/apog/work/CNN_bin/miscellaneous/kept_genes.txt"

parallel_vars <- c("build_attention_mlp", "dropouts", "dropouts2", "seeds", 
                   "training", "training_target", "batchsizes", "Epochs", "j", "jj")

# Read the file, skipping the first row (header is irregular)
file_contents <- read.table(
  file_path,
  header = FALSE,
  skip = 1,
  col.names = c("Index", "Gene"),
  stringsAsFactors = FALSE
)$Gene

file_contents = file_contents[1:30]
cat("Number of genes:", length(file_contents), "\n")

# Reproducible seeds for tuning
set.seed(123)
seeds <- sample(1:10000, 30)
print(seeds)

df <- NULL
acc_nested_list <- list()
validation_split <- 0.3

# -------------------- Helper: rescale back to original --------------------
rescale <- function(x, min_val, max_val) {
  (x * (max_val - min_val)) + min_val
}

# -------------------- Helper: model builder with attention --------------------
# rate1 = dropout after Dense(100), rate2 = dropout after Dense(50)
build_attention_mlp <- function(input_dim, rate1, rate2) {
  input_layer <- layer_input(shape = c(input_dim), name = "input")
  
  # Feature-wise attention (per-sample): weights sum to 1 across features
  attention_weights <- input_layer %>%
    layer_dense(units = input_dim,
                activation = "softmax",
                name = "attention_weights")
  
  # Element-wise weighting of inputs by attention
  weighted_input <- layer_multiply(list(input_layer, attention_weights),
                                   name = "weighted_input")
  
  # 2-layer MLP
  dense1 <- weighted_input %>%
    layer_dense(units = 100, activation = "relu", name = "dense1") %>%
    layer_dropout(rate = rate1, name = "dropout1")
  
  dense2 <- dense1 %>%
    layer_dense(units = 50, activation = "relu", name = "dense2") %>%
    layer_dropout(rate = rate2, name = "dropout2")
  
  output_layer <- dense2 %>%
    layer_dense(units = 1, activation = "linear", name = "output")
  
  model <- keras_model(inputs = input_layer, outputs = output_layer)
  model %>% compile(
    loss = "mse",
    optimizer = optimizer_adam(),
    metrics = c("mse")
  )
  model
}

# -------------------- Main loop over genes --------------------
for (p in file_contents) {
  
  tryCatch({
    gene_name <- p
    start_time <- Sys.time()
    
    # ---------- Load per-gene matrix ----------
    gene_file <- paste0("/projects/apog/work/input/IHEC_Activity_1MB_hg38/", gene_name, ".txt.gz")
    if (file.exists(gene_file)) {
      cat("Loading:", gene_file, "\n")
      df <- read.table(gene_file, header = TRUE, sep = "\t")
    } else {
      cat("File not found for gene:", gene_name, "\n")
      next
    }
    
    cat("Processing gene:", gene_name, "\n")
    
    # ---------- Cluster for tuning ----------
    cl <- makeCluster(52)
    
    # ---------- Preprocess ----------
    df <- df[order(df$Sample), ]         # order by sample name
    rownames(df) <- NULL
    df <- df[, colSums(df != 0) > 0]     # drop all-zero columns
    df_samples <- df[, 1]
    df2 <- df[, 2:ncol(df)]
    df2 <- log2(df2 + 1)                 # log2 transform
    df3 <- cbind(df_samples, df2)
    
    # Partition lists
    test_samples <- "/projects/apog/work/CNN_bin/miscellaneous/partition0_test.csv"
    train_samples <- "/projects/apog/work/CNN_bin/miscellaneous/partition0_train.csv"
    
    test_sample2 <- read.csv(test_samples)
    test_sample3 <- test_sample2[, 1]
    test_data <- df3[df3$df_samples %in% test_sample3, ]
    
    train_samples2 <- read.csv(train_samples)
    train_samples3 <- train_samples2[, 1]
    train_data <- df3[df3$df_samples %in% train_samples3, ]
    
    # Keep only features + target (drop sample col)
    train_data <- train_data[, 2:ncol(train_data)]
    test_data  <- test_data[,  2:ncol(test_data)]
    train_traget_just_log_no_norm <- train_data[, ncol(train_data)]
    
    # Min-max normalization (fit on train, apply to test)
    train_min_values <- apply(train_data, 2, min)
    train_max_values <- apply(train_data, 2, max)
    
    normalize_min_max <- function(v) {
      if (max(v) == min(v)) return(rep(0, length(v)))
      (v - min(v)) / (max(v) - min(v))
    }
    
    train_data_normalized <- as.data.frame(apply(train_data, 2, normalize_min_max))
    
    min_max_normalize_test <- function(data, min_values, max_values) {
      normalized_data <- data.frame(
        lapply(names(data), function(col) {
          min_val <- min_values[col]
          max_val <- max_values[col]
          if (max_val == min_val) {
            rep(0, length(data[[col]]))
          } else {
            (data[[col]] - min_val) / (max_val - min_val)
          }
        })
      )
      colnames(normalized_data) <- colnames(data)
      normalized_data
    }
    
    test_data_normalized <- min_max_normalize_test(
      data = test_data, min_values = train_min_values, max_values = train_max_values
    )
    
    # Split X / y
    training        <- train_data_normalized[, 1:(ncol(train_data_normalized) - 1)]
    training_target <- train_data_normalized[, ncol(train_data_normalized)]
    test1           <- test_data_normalized[, 1:(ncol(test_data_normalized) - 1)]
    test1_target    <- test_data_normalized[, ncol(test_data_normalized)]
    
    cat("training class:", class(training), "target class:", class(training_target), "\n")
    
    training <- as.matrix(training)
    
    # ---------- Hyperparams ----------
    dropouts   <- c(0.2, 0.4)
    dropouts2  <- c(0.2, 0.4)
    batchsizes <- 32
    Epochs     <- 1200
    
    best_loss2 <- Inf
    b_hyper_dropouts <- NA
    b_hyper_dropouts2 <- NA
    best_seed <- NA
    
    # ---------- Tuning: 2 hidden layers with attention ----------
    cat("Two hidden layers mode with Attention\n")
    for (j in 1:length(dropouts)) {
      for (jj in 1:length(dropouts2)) {
        
        # Export data/vars to workers
        clusterExport(cl, parallel_vars, envir = environment())
        
        a2 <- parLapply(cl, seeds, function(seed) {
          require(keras)
          require(tensorflow)
          tensorflow::tf$random$set_seed(seed)
          
          input_dim <- ncol(training)
          model <- build_attention_mlp(
            input_dim = input_dim,
            rate1 = dropouts[j],
            rate2 = dropouts2[jj]
          )
          
          # Callbacks
          print_dot_callback <- callback_lambda(
            on_epoch_end = function(epoch, logs) {
              if (epoch %% 100 == 0) cat("\n")
              cat(".")
            }
          )
          early_stop <- callback_early_stopping(monitor = "val_loss", mode = "min", patience = 30)
          
          # Train
          history <- model %>% fit(
            x = training,
            y = training_target,
            epochs = Epochs,
            batch_size = batchsizes,
            shuffle = FALSE,            # NOTE: keras arg is 'shuffle'
            validation_split = 0.3,
            verbose = 0,
            callbacks = list(early_stop, print_dot_callback)
          )
          
          temp_loss <- mean(history$metrics$val_loss)
          
          list(
            loss = temp_loss,
            seed = seed
          )
        }) # parLapply
        
        # Track best config
        for (w in seq_along(a2)) {
          if (a2[[w]]$loss < best_loss2) {
            best_loss2        <- a2[[w]]$loss
            b_hyper_dropouts  <- dropouts[j]
            b_hyper_dropouts2 <- dropouts2[jj]
            best_seed         <- a2[[w]]$seed
          }
        }
        
      } # jj
    }   # j
    
    validation_err <- best_loss2
    cat("Best val_loss:", validation_err,
        "with dropout1:", b_hyper_dropouts, "dropout2:", b_hyper_dropouts2,
        "seed:", best_seed, "\n")
    
    # ---------- Re-train best model on full training data ----------
    tensorflow::tf$random$set_seed(best_seed)
    Best_model <- build_attention_mlp(
      input_dim = ncol(training),
      rate1 = b_hyper_dropouts,
      rate2 = b_hyper_dropouts2
    )
    
    early_stop <- callback_early_stopping(monitor = "val_loss", mode = "min", patience = 30)
    print_dot_callback <- callback_lambda(
      on_epoch_end = function(epoch, logs) {
        if (epoch %% 100 == 0) cat("\n")
        cat(".")
      }
    )
    
    Best_ModelFited <- Best_model %>% fit(
      x = training,
      y = training_target,
      epochs = Epochs,
      batch_size = batchsizes,
      shuffle = FALSE,
      validation_split = 0.3,
      verbose = 0,
      callbacks = list(early_stop, print_dot_callback)
    )
    
    # ---------- Save model ----------
    model_dir <- "/projects/apog/work/models/1MB/attention_mode/models/"
    dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
    name1 <- file.path(model_dir, paste0(gene_name, ".hdf5"))
    save_model_hdf5(Best_model, name1)
    cat("Best model is saved ->", name1, "\n")
    
    # ---------- Evaluate on test & train ----------
    min_origin_target <- min(train_traget_just_log_no_norm)
    max_origin_target <- max(train_traget_just_log_no_norm)
    
    test1 <- as.matrix(test1)
    predicted_test <- Best_model %>% predict(test1)
    predicted_test <- as.numeric(predicted_test)
    cat("test_statistics\n")
    MSE1 <- mean((test1_target - predicted_test)^2)
    cor1 <- suppressWarnings(cor(test1_target, predicted_test, method = "pearson"))
    scaled_predicted_target <- rescale(predicted_test, min_val = min_origin_target, max_val = max_origin_target)
    scaled_original_target  <- rescale(test1_target,  min_val = min_origin_target, max_val = max_origin_target)
    MSE_backscaled <- mean((scaled_original_target - scaled_predicted_target)^2)
    
    predicted_train <- Best_model %>% predict(training)
    predicted_train <- as.numeric(predicted_train)
    cat("train_statistics\n")
    MSE_train <- mean((training_target - predicted_train)^2)
    Cor_train <- suppressWarnings(cor(training_target, predicted_train, method = "pearson"))
    scaled_predicted_target_train <- rescale(predicted_train, min_val = min_origin_target, max_val = max_origin_target)
    scaled_original_target_train  <- rescale(training_target, min_val = min_origin_target, max_val = max_origin_target)
    MSE_backscaled_train <- mean((scaled_original_target_train - scaled_predicted_target_train)^2)
    
    end_time <- Sys.time()
    runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    accuracy1 <- data.frame(
      Name = gene_name,
      MSE_test_NN2 = MSE1,
      Cor_test_NN2 = cor1,
      MSE_S_test_NN2 = MSE_backscaled,
      MSE_train_NN2 = MSE_train,
      Cor_train_NN2 = Cor_train,
      MSE_S_train_NN2 = MSE_backscaled_train,
      val_error = validation_err,
      runtime_secs = runtime
    )
    
    acc_nested_list[[p]] <- accuracy1
    
    # ---------- Save aggregate stats ----------
    stat_dir <- "/projects/apog/work/models/1MB/attention_mode/stat"
    dir.create(stat_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(acc_nested_list, file.path(stat_dir, "NN2_attention_mode.RDS"))
    
    # ---------- Extract & save attention weights ----------
    # Build a small model to output the attention layer
    att_extractor <- keras_model(
      inputs = Best_model$input,
      outputs = get_layer(Best_model, "attention_weights")$output
    )
    att_matrix <- predict(att_extractor, training, batch_size = 128)  # n_samples x n_features
    
    # Global feature importance = mean attention over samples
    feature_importance <- colMeans(att_matrix)
    
    att_tbl <- data.frame(
      feature = colnames(training),
      mean_attention = feature_importance
    )
    
    att_dir <- "/projects/apog/work/models/1MB/attention_mode/weights"
    dir.create(att_dir, recursive = TRUE, showWarnings = FALSE)
    att_csv <- file.path(att_dir, paste0(gene_name, "_attention.csv"))
    write.csv(att_tbl, att_csv, row.names = FALSE)
    cat("Attention weights saved ->", att_csv, "\n")
    
    # ---------- Plot top-k features by attention ----------
    k <- min(30, length(feature_importance))
    top_idx <- order(feature_importance, decreasing = TRUE)[1:k]
    plot_df <- data.frame(
      feature = factor(colnames(training)[top_idx], levels = rev(colnames(training)[top_idx])),
      attention = feature_importance[top_idx]
    )
    
    p_plot <- ggplot(plot_df, aes(x = feature, y = attention)) +
      geom_col() +
      coord_flip() +
      labs(
        title = paste0(gene_name, ": Top ", k, " features by attention"),
        x = "Feature",
        y = "Mean attention weight"
      ) +
      theme_minimal()
    
    plot_dir <- "/projects/apog/work/models/1MB/attention_mode/plots"
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(filename = file.path(plot_dir, paste0(gene_name, "_attention_top", k, ".pdf")),
           plot = p_plot, width = 7, height = 6)
    
    cat("Saved attention plot\n")
    
    # ---------- Done with this gene ----------
    cat("Finished gene:", gene_name, "\n")
    
    # Stop cluster for this gene
    stopCluster(cl)
    
  }, error = function(e) {
    cat("ERROR in gene", p, ":", conditionMessage(e), "\n")
    # try to safely stop cluster if it exists
    try({ stopCluster(cl) }, silent = TRUE)
  })
  
} # end for file_contents

cat("All done.\n")
