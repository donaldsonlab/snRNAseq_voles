input_dir <- "/scratch/alpine/libr8020/shuffle_corrs/min1ct/"

my_files <- list.files(input_dir)
print(my_files)

all_corrs <- data.frame()
for (f in my_files) {
  df <- read.csv(paste0(input_dir, f))
  all_corrs <- rbind(all_corrs, df)
  
}

write.csv(all_corrs, "/scratch/alpine/libr8020/shuffle_corrs/combined_file/all_corrs_combined_min1ct.csv")
