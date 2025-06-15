# ===================================================================
# Plant Phenology Image Classification using Keras and TensorFlow in R
# ===================================================================
# This script trains a Convolutional Neural Network (CNN) to classify
# images of plants into four phenological states:
# 1. Flowering
# 2. Fruiting
# 3. Both flowering and fruiting
# 4. Neither
#
# The data pipeline is built using the 'tfdatasets' package for
# maximum performance and compatibility with TensorFlow.
# ===================================================================


# -------------------------------------------------------------------
# STEP 1: SETUP AND ENVIRONMENT PREPARATION
# -------------------------------------------------------------------

# --- Install necessary packages (if not already installed) ---
# install.packages("keras")
# install.packages("tensorflow")
# install.packages("tfdatasets")
# install.packages("magick") # For creating dummy data

# --- Load libraries ---
library(keras)
library(tensorflow)
library(tfdatasets)
library(magick)
Sys.setenv(CUDA_VISIBLE_DEVICES = "-1")

# --- (Optional) Force TensorFlow to use only the CPU ---
# Uncomment the line below if you encounter GPU/CUDA errors. This is the
# simplest way to bypass driver issues and get the code running.
# 

# --- (Optional but Recommended) Install Keras backend ---
# Run this command once if you haven't set up the Keras/TensorFlow backend.
# It creates an isolated Python environment for R to use.
# install_keras()


# -------------------------------------------------------------------
# STEP 2: LOAD DATA
# -------------------------------------------------------------------

setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")
image_dir<-"/media/huijieqiao/WD22T_11/invasive.plant.phenology/Data/Images"
image<-list.files(image_dir)
data<-readRDS("../Data/iNadata_ll.rda")

data$image.path<-sprintf("%s/%d.jpg", image_dir, data$id)
data$with.image<-file.exists(data$image.path)
data$image.name<-sprintf("%d.jpg", data$id)

data<-data[which(data$with.image==T),]



# -------------------------------------------------------------------
# STEP 3: DATA PREPARATION AND PREPROCESSING
# -------------------------------------------------------------------

# --- Convert text labels to numeric class IDs (0, 1, 2, 3) ---
phenology_classes<-unique(data$phenology)
phenology_data<-as.data.frame(data[, c("image.name", "phenology")])
phenology_data$phenology <- factor(phenology_data$phenology, levels = phenology_classes)
phenology_data$class_id <- as.numeric(phenology_data$phenology) - 1

# --- Split data into Training (70%), Validation (15%), and Test (15%) sets ---
set.seed(123) # for reproducibility
num_samples <- nrow(phenology_data)
all_indices <- 1:num_samples

train_indices <- sample(all_indices, size = round(0.8 * num_samples))
remaining_indices <- setdiff(all_indices, train_indices)
val_indices <- sample(remaining_indices, size = round(0.5 * length(remaining_indices))) # 50% of remaining
test_indices <- setdiff(remaining_indices, val_indices)

train_data <- phenology_data[train_indices, ]
val_data <- phenology_data[val_indices, ]
test_data <- phenology_data[test_indices, ]

cat(sprintf("Data split: %d training, %d validation, %d test images.\n",
            nrow(train_data), nrow(val_data), nrow(test_data)))


# -------------------------------------------------------------------
# STEP 4: CREATE TF.DATA PIPELINE
# -------------------------------------------------------------------
# This is the most robust way to feed data to a Keras model in R.

# --- Define constants and preprocessing functions ---
img_width <- 224
img_height <- 224
batch_size <- 32

data_augmentation_model <- keras_model_sequential(
  list(
    layer_random_flip(mode = "horizontal"),
    layer_random_rotation(factor = 0.2),
    layer_random_zoom(height_factor = 0.2, width_factor = 0.2),
    layer_random_contrast(factor = 0.2)
  )
)
# Function to decode and preprocess a single image
decode_and_resize <- function(path) {
  image_raw <- tf$io$read_file(path)
  # Use decode_jpeg or decode_png depending on your image format
  image <- tf$image$decode_jpeg(image_raw, channels = 3) 
  image <- tf$image$resize(image, size = shape(img_width, img_height))
  image <- image / 255.0 # Normalize pixel values to [0, 1]
  return(image)
}

# --- Convert R variables to TensorFlow constants for use in the graph ---
# This is crucial to avoid "not convertible to tensor" errors.
image_dir_tf <- tf$constant(image_dir)
filesep_tf <- tf$constant(.Platform$file.sep)
num_classes_tf <- tf$constant(length(phenology_classes), dtype = "int32")

# Main function to process a single element (image.name and class_id) of the dataset
# --- Create a function FACTORY that builds our mapping function ---
# This is a more advanced but safer pattern.
create_process_function <- function(image_filenames, image_class_ids) {
  
  # 1. Convert all necessary R data into TensorFlow constants INSIDE the factory
  image_dir_tf <- tf$constant(image_dir)
  filesep_tf <- tf$constant(.Platform$file.sep)
  num_classes_tf <- tf$constant(length(phenology_classes), dtype = "int32")
  
  # Convert the entire vectors of filenames and class_ids into TF constants
  filenames_tf <- tf$constant(image_filenames)
  class_ids_tf <- tf$constant(image_class_ids)
  
  # 2. The factory returns the ACTUAL function that `dataset_map` will use
  # This function takes only an INDEX as input.
  process_index <- function(index) {
    # Get the specific filename and class_id using the index
    filename <- tf$gather(filenames_tf, index)
    class_id <- tf$gather(class_ids_tf, index)
    
    # Now, build the path and process the image as before
    full_path <- tf$strings$join(list(image_dir_tf, filesep_tf, filename))
    image <- decode_and_resize(full_path)
    label <- tf$one_hot(tf$cast(class_id, dtype = "int32"), depth = num_classes_tf)
    
    return(list(image, label))
  }
  
  return(process_index)
}


# --- Create the training dataset pipeline ---
process_train_index <- create_process_function(train_data$image.name, train_data$class_id)

# b. Create a dataset of INDICES, from 0 to N-1
train_ds <- tf$data$Dataset$range(as.integer(nrow(train_data))) %>%
  dataset_map(process_train_index, num_parallel_calls = tf$data$AUTOTUNE) %>%
  dataset_shuffle(buffer_size = 1000) %>%
  dataset_batch(batch_size) %>%
  dataset_prefetch(buffer_size = tf$data$AUTOTUNE)

# 2. For the validation set
process_val_index <- create_process_function(val_data$image.name, val_data$class_id)
val_ds <- tf$data$Dataset$range(as.integer(nrow(val_data))) %>%
  dataset_map(process_val_index, num_parallel_calls = tf$data$AUTOTUNE) %>%
  dataset_batch(batch_size) %>%
  dataset_prefetch(buffer_size = tf$data$AUTOTUNE)

# 3. For the test set
process_test_index <- create_process_function(test_data$image.name, test_data$class_id)
test_ds <- tf$data$Dataset$range(as.integer(nrow(test_data))) %>%
  dataset_map(process_test_index, num_parallel_calls = tf$data$AUTOTUNE) %>%
  dataset_batch(batch_size) %>%
  dataset_prefetch(buffer_size = tf$data$AUTOTUNE)

cat("All tf.data pipelines created successfully using the robust factory pattern.\n")



# -------------------------------------------------------------------
# STEP 5: BUILD THE CNN MODEL
# -------------------------------------------------------------------
# 1. Define the input tensor. This is the entry point of your model.
# 1. Define the input tensor.
inputs <- layer_input(shape = c(img_width, img_height, 3))

# 2. Chain all layers starting from the input tensor.
#    This creates the output tensor in a single pipeline.
outputs <- inputs %>%
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = 'relu') %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  
  layer_conv_2d(filters = 64, kernel_size = c(3, 3), activation = 'relu') %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  
  layer_conv_2d(filters = 128, kernel_size = c(3, 3), activation = 'relu') %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  
  layer_flatten() %>%
  layer_dense(units = 512, activation = 'relu') %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = length(phenology_classes), activation = 'softmax')

# 3. Create the final model.
model <- keras_model(inputs = inputs, outputs = outputs)


# --- Print model summary (using the most robust method) ---
model$summary()

# -------------------------------------------------------------------
# STEP 6: COMPILE AND TRAIN THE MODEL
# -------------------------------------------------------------------


# --- Compile the model ---
model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_adam(),
  metrics = c('accuracy')
)
cat("Model compiled successfully using direct Python method.\n")

# --- Train the model ---
cat("Starting model training...\n")
history <- model %>% fit(
  train_ds,
  epochs = 50, # Adjust as needed
  validation_data = val_ds
)
cat("Model training complete.\n")

# --- Visualize training history ---
plot(history)


# -------------------------------------------------------------------
# STEP 7: EVALUATE THE MODEL ON THE TEST SET
# -------------------------------------------------------------------
# The test set was held out and not seen by the model during training.
# This provides an unbiased estimate of the model's performance.

cat("\nEvaluating model on the test set...\n")
evaluation_results <- model %>% evaluate(test_ds)

cat("\n--- Final Test Results ---\n")
cat(sprintf("Test Loss: %.4f\n", evaluation_results$loss))
cat(sprintf("Test Accuracy: %.2f%%\n", evaluation_results$accuracy * 100))


# -------------------------------------------------------------------
# STEP 8: DETAILED ANALYSIS AND PREDICTION
# -------------------------------------------------------------------

# --- Make predictions on the entire test set ---
cat("\nMaking predictions on the test set for detailed analysis...\n")
predictions_matrix <- model %>% predict(test_ds)

# Get the class index with the highest probability for each prediction
predicted_class_indices <- apply(predictions_matrix, 1, which.max)
predicted_labels <- phenology_classes[predicted_class_indices]

# Get the true labels from the test data frame
true_labels <- test_data$phenology

# --- Generate a confusion matrix ---
confusion_matrix <- table(Predicted = predicted_labels, True = true_labels)

cat("\n--- Confusion Matrix ---\n")
print(confusion_matrix)


# --- Function to predict a single, new image ---
predict_single_image <- function(model, image_path) {
  # 1. Load and preprocess the image using our TF function
  image_tensor <- decode_and_resize(tf$constant(image_path))
  
  # 2. Add a batch dimension (batch_size = 1)
  image_tensor <- tf$expand_dims(image_tensor, axis = 0)
  
  # 3. Predict
  pred_probs <- model %>% predict(image_tensor)
  
  # 4. Decode the prediction
  pred_index <- which.max(pred_probs[1,])
  pred_label <- phenology_classes[pred_index]
  confidence <- max(pred_probs[1,])
  
  return(list(label = pred_label, confidence = confidence, probabilities = pred_probs[1,]))
}

# --- Example of using the single prediction function ---
# (Replace with a real path to an image)
# test_image_path <- test_data$filename[1]
# full_test_image_path <- file.path(image_dir, test_image_path)
# 
# if (file.exists(full_test_image_path)) {
#   single_pred <- predict_single_image(model, full_test_image_path)
#   cat(sprintf("\nPrediction for '%s':\n  Label: %s\n  Confidence: %.2f%%\n",
#               basename(full_test_image_path),
#               single_pred$label,
#               single_pred$confidence * 100))
# } else {
#   cat("\nSingle test image not found. Skipping single prediction example.\n")
# }

# ======================= END OF SCRIPT =======================