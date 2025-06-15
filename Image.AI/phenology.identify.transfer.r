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


# --- (Optional) Force TensorFlow to use only the CPU ---
# Uncomment the line below if you encounter GPU/CUDA errors. This is the
# simplest way to bypass driver issues and get the code running.
# Sys.setenv(CUDA_VISIBLE_DEVICES = "-1")

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
    layer_random_flip("horizontal"),
    layer_random_rotation(0.2),
    layer_random_zoom(0.2),
    layer_random_contrast(0.2)
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
# 1. Load the pre-trained ResNet50V2 model without its top classification layer.
#    The weights are pre-trained on the ImageNet dataset.
base_model <- application_resnet50_v2(
  weights = "/home/huijieqiao/.keras/models/resnet50v2_weights_tf_dim_ordering_tf_kernels_notop.h5",
  include_top = FALSE, # Crucial: we provide our own classifier head.
  input_shape = c(img_width, img_height, 3)
)

# 2. Freeze the weights of the base model.
#    This prevents them from being updated during the initial training phase.
#    We are only training our new classifier head.
base_model$trainable <- FALSE

# 3. Create our new model on top of the base model using the Functional API.
#    Define the input matching the base model's input.
inputs <- layer_input(shape = c(img_width, img_height, 3))

# The input is first passed to the base_model.
# We set training=FALSE so that layers like BatchNormalization run in inference mode.
x <- base_model(inputs, training = FALSE)

# Add our custom classifier head on top of the base model's output.
x <- x %>%
  layer_global_average_pooling_2d() %>% # A good practice to reduce dimensions
  layer_dense(units = 256, activation = "relu") %>%
  layer_dropout(rate = 0.5) # Dropout is still very useful here

outputs <- x %>%
  layer_dense(units = length(phenology_classes), activation = 'softmax')

# 4. Create the final model object.
model <- keras_model(inputs = inputs, outputs = outputs)


# --- Print model summary ---
# You will see that most parameters are non-trainable, which is correct.
summary(model)


# ===================================================================
# STEP 6: COMPILE AND TRAIN THE MODEL (WITH EARLY STOPPING)
# ===================================================================

# --- Use the robust direct Python method call for compiling ---
model$compile(
  loss = 'categorical_crossentropy',
  # Use a slightly lower learning rate, which is common for transfer learning
  optimizer = optimizer_adam(learning_rate = 0.0005), 
  metrics = reticulate::r_to_py(c('accuracy'))
)

cat("Model compiled successfully using direct Python method.\n")


# --- Define an EarlyStopping callback ---
# This will monitor the validation loss and stop training if it
# stops improving, preventing overfitting.
early_stopping <- callback_early_stopping(
  monitor = "val_loss", # Monitor validation loss
  patience = 5,         # Stop if it doesn't improve for 5 consecutive epochs
  restore_best_weights = TRUE # IMPORTANT: Restores the model weights from the epoch with the best val_loss
)


# --- Train the model ---
# We can increase the number of epochs because EarlyStopping will handle it.
cat("Starting model training with transfer learning...\n")
history <- model$fit(
  x = train_ds,
  epochs = 50L, # Set a high number, early stopping will find the optimal point
  validation_data = val_ds,
  callbacks = list(early_stopping) # Add the callback here
)
cat("Model training complete.\n")

# --- Visualize training history ---
plot(history)
