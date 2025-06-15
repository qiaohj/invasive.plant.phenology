
Sys.setenv(CUDA_VISIBLE_DEVICES = "0")
library(keras)
library(tensorflow)
gpus <- tf$config$list_physical_devices("GPU")
cat("Physical GPUs detected:\n")
print(gpus)

# 测试一下：简单运算看看在哪执行
tf$debugging$set_log_device_placement(TRUE)
a <- tf$constant(matrix(1, 2, 2))
b <- tf$constant(matrix(1, 2))
c <- tf$matmul(a, b)
print(c)

# 2. 创建一个简单的数据向量，就像您的 data$class_id 一样
# 这个向量包含0, 1, 2, 3，代表四个类别
test_ids <- c(0, 1, 2, 3, 0, 1, 2) 

# 3. 定义类别的总数
num_classes <- 4

# 4. 调用 to_categorical 函数 (使用我们之前确认的最稳健的写法)
# 我们用 tryCatch 包裹它，这样即使出错也不会崩溃，而是会打印信息
tryCatch({
  
  labels_test <- keras::to_categorical(test_ids, num_classes = num_classes)
  
  d# 如果成功，打印结果
  print("成功！to_categorical 函数运行正常。")
  print(labels_test)
  
}, error = function(e) {
  
  # 如果失败，打印错误信息
  print("失败了。问题仍然存在，即使在干净的环境中。")
  print(e)
  
})