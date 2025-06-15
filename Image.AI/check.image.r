library(magick)
library(data.table)
setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")
current_options <- magick::magick_options()
Sys.setenv(MAGICK_MEMORY_LIMIT = "1024GiB")

# Set max disk limit to 16 Gigabytes
Sys.setenv(MAGICK_DISK_LIMIT = "1024GiB")


image_dir<-"/media/huijieqiao/WD22T_11/invasive.plant.phenology/Data/Images"
image<-list.files(image_dir)
is_jpeg <- function(file_path) {
  info <- tryCatch({
    image_info(image_read(file_path))
  }, error = function(e) {
    cat(sprintf("Warning: Could not read file: %s\n", file_path))
    return(NULL)
  })
  
  if (is.null(info)) {
    return(-1)
  }
  return(ifelse(tolower(info$format) != "jpeg", 0, 1))
}
f<-image[1]
all<-list()
for (i in c(1:length(image))){
  f<-image[i]
  if (i %% 1e2==1){
    print(paste(i, length(image)))
  }
  is.jpg<-is_jpeg(sprintf("%s/%s", image_dir, f))
  item<-data.table(file=f, is.jpg=is.jpg)
  all[[length(all)+1]]<-item
}
all.df<-rbindlist(all)
no<-all.df[is.jpg==0]
for (i in c(1:nrow(no))){
  f<-no[i]$file
  f<-sprintf("%s/%s", image_dir, f)
  img <- image_read(f)
  
  info<-image_info(img)
  format<-tolower(info$format)
  new.f<-gsub("jpg", format, f)
  new.f<-gsub("/Images/", "/Images.Others/", new.f)
  file.rename(f, new.f)
  img <- image_flatten(img) 
  image_write(img, path = f, format = "jpg")
  
  
}

file_path<-"/media/huijieqiao/WD22T_11/invasive.plant.phenology/Data/Images/100176976.jpg"
