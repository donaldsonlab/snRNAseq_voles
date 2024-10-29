## title: "Hotspot R analysis"
# Define a vector with all required packages
# (excluding Seurat and ComplexHeatmap for now)

required_packages <- c("dplyr",
                       "tidyr",
                       "ggplot2",
                       "corrr",
                       "stringr",
                       "Hmisc",
                       "ggpubr")


# Install any missing packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("corrr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("Hmisc"))
suppressPackageStartupMessages(library("ggpubr"))

#' This function reads a CSV file into a data frame.
#' @param file_name A string specifying the path to the CSV file to be read.

read_csv_file <- function(file_name) {
  tryCatch({
    data <- read.csv(file_name)
    return(data)
  },
  error = function(e) {
    if (grepl("cannot open the connection", e$message)) {
      warning(paste(" File not found:", file_name))
    } else {
      warning("An error occurred: ", e$message)
    }
    return(NA)
  }, finally = {
    message(" CSV read operation attempted.")
  }
  )
}


#' Process Metadata
#'
#' This function processes a metadata data frame by combining
#' 'sex' and 'SS_OS' columns into a new 'Group' column.
#'
#' @param metadata Data frame containing animal metadata to be processed.
#'
#' @return Returns the modified data frame with an additional 'Group' column,
#'         which is a concatenation of 'sex' and 'SS_OS'
#'         columns separated by an underscore.

process_metadata <- function(metadata) {
  tryCatch(
    {
      if (!is.data.frame(metadata)) {
        stop("metadata must be a data frame")
      }
      if (!("sex" %in% names(metadata)) || !("SS_OS" %in% names(metadata))) {
        stop("Both 'sex' and 'SS_OS' columns must exist")
      }
      metadata$Group <- paste(metadata$sex, metadata$SS_OS, sep = "_")
      return(metadata)
    },
    error = function(e) {
      warning("An error ocurred: ", e$message)
      return(NA)
    }
  )
}


#' Merge Two Data Frames
#'
#' This function merges two data frames based on a common column.
#'
#' @param df1 The first data frame to be merged.
#' @param df2 The second data frame to be merged.
#' @param common_column String specifying the column to merge on
#'
#' @return Returns a merged data frame based on the specified common column.

merge_dataframes <- function(df1, df2, common_column) {
  tryCatch(
    {
      if (!is.data.frame(df1)) {
        stop("Dataset 1 must be a dataframe")
      }
      if (!is.data.frame(df2)) {
        stop("Dataset 2 must be a dataframe")
      }
      if (!(common_column %in% names(df1)) ||
            !(common_column %in% names(df2))) {
        stop("Column to merge on must exist in both dataframes")
      }

      merged_data <- merge(df1, df2, by = common_column)
      return(merged_data)
    },
    error = function(e) {
      warning("An error ocurred: ", e$message)
    }
  )
}


#' Set Working Directory
#
#' @param directory_path A string specifying the directory path.
#'
#' @return There is no return value from this function.

set_working_directory <- function(directory_path) {
  tryCatch(
    {
      setwd(directory_path)
      cat("Working directory set to:", getwd(), "\n")
    },
    error = function(e) {
      if (grepl("cannot open the connection", e$message)) {
        warning(paste("Directory not found:", directory_path))
      } else {
        warning("An error occurred: ", e$message)
      }
      return(NA)
    }, finally = {
      message("Setting working directory attempted")
    }
  )
}


#' Set Column Comparisons
#'
#' This function generates a list of column comparisons based on provided
#' groups and sexes.
#' Each comparison is a pair of strings formed by concatenating elements from
#' 'sexes' and 'groups'.
#'
#' @param groups A vector of group names to be used in the comparisons.
#' @param sexes A vector of sexes to be used in the comparisons.
#'
#' @return Returns a list of string pairs. Each pair is created by
#' concatenating each sex with '_SS' and '_OS', and each group element.
#'
set_column_comparisons <- function(groups, sexes) {
  tryCatch({
    if (!is.vector(groups) || !is.vector(sexes)) {
      stop("Groups and sexes both need to be vectors")
    }
    comparisons <- list()
    for (sex in sexes) {
      for (group in groups) {
        comparisons <- c(comparisons, list(c(paste0(sex, "_SS"),
                                             paste0(sex, "_OS"))))
      }
    }
    return(comparisons)
  }, error = function(e) {
    warning("An error occurred: ", e$message)
    return(NULL)
  })
}


#' Set Group Levels in a Data Frame
#'
#' This function sets the 'Group' column to a specific order. .
#'
#' @param data .csv read in previous function that contains 'Group' data.
#'
#' @return Returns the data frame with the 'Group' column converted to a factor
#'         having levels "F_SS", "F_OS", "M_SS", "M_OS".
group_level <- function(data) {
  tryCatch(
           {
             if (!is.data.frame(data)) {
               stop("Data must be a dataframe")
             }
             if (!("Group" %in% names(data))) {
               stop("The 'Group' column does not exist in the data frame.")
             }
             data$Group <- factor(data$Group,
                                  levels = c("F_SS", "F_OS", "M_SS", "M_OS"))
             return(data)

           }, error = function(e) {
             warning("An error occurred: ", e$message)
             return(NA)
           })
}


#' Set Factor Levels for a Specific Column in a Dataset
#'
#' This function converts a specified column in a dataset into a factor with
#' defined levels, provided as a vector.
#'
#' @param dataset Data frame containing the column to be converted to a factor.
#' @param column_name Name of the column in the dataset as a string
#' @param levels_vector A vector specifying the levels for the factor.
#'
#' @return Returns the modified data frame with the specified column c
#' onverted to a factor with the provided levels.

set_factor_variables <- function(dataset, column_name, levels_vector) {
  tryCatch(
    {
      if (!(column_name %in% names(dataset))) {
        stop("Column does not exist in the dataframe")
      }
      if (!is.vector(levels_vector)) {
        stop("Levels must be a vector.")
      }
      dataset[[column_name]] <- factor(dataset[[column_name]],
                                       levels = levels_vector)
      return(data)

    }, error = function(e) {
      warning("An error occurred: ", e$message)
      return(NA)
    }
  )
}


#' Generate Column Names from a Specified Range in a Dataset
#'
#' This function extracts and returns the names of columns
#' in a dataset within a specified range.
#'
#' @param dataset Data frame from which column names are to be extracted.
#' @param first_column String representing the name of the starting column
#' @param last_column String representing the name of the ending column.
#'
#' @return Returns a vector of column names from the specified range.

generate_colnames <- function(dataset, first_column, last_column) {
  tryCatch({
    if (!(first_column %in% names(dataset)) || !(last_column %in%
                                                   names(dataset))) {
      stop("Both column names must exist in the specified dataset")
    }
    data_columns <- colnames(dataset)[first_column:last_column]
    return(data_columns)

    error <- function(e) {
      warning("An error occurred: ", e$message)
      return(NA)
    }
  })
}


#' Extract Module Column Names from a Dataset
#'
#' This function identifies and returns the names of all columns in the data
#' that contain the word 'Module'.
#'
#' @param dataset Data frame or tibble with desired columns.
#'
#' @return Returns a vector containing the names of all columns.
#'
get_modules <- function(dataset) {
  tryCatch(
           {
             if (!is.data.frame(dataset) && !("tbl_df" %in% class(dataset))) {
               stop("Input 'dataset' must be a data frame or tibble.")
             }
             mods <- select(dataset, contains("Module")) %>% colnames(.)
             return(mods)

             error <- function(e) {
               warning("An error occurred: ", e$message)
               return(NA)
             }
           })
}


#' Create and Optionally Save a Violin Plot
#'
#' This function creates a violin plot for a given dataset and module.
#' It compares groups based on predefined comparisons and
#' optionally saves the plot as a PNG file.
#'
#' @param data A data frame containing the data to be plotted.
#'The data frame must have a 'Group' column and a column for module scores.
#'
#' @param save A logical value indicating whether to save the plot.
#'             If TRUE, the plot is saved as a PNG file.
#'
#' @details The function performs the following operations:
#'          - Creates a violin plot that shows distribution of module scores
#'            across different groups.
#'          - Groups are "F_SS" vs "F_OS", "F_SS" vs "M_SS", "M_SS" vs "M_OS",
#'            and "F_OS" vs "M_OS". F = female, M = male, SS = same sex,
#'            OS = opposite sex.
#'          - The plot includes statistical comparison using the Wilcoxon test.
#'          - If 'save' is TRUE, the plot is saved as a PNG file.
#'          - Function assumes that global variable 'mod' holds name of module.
#'
violin <- function(data, save) {
  my_comparisons <- list(c("F_SS", "F_OS"),
                         c("F_SS", "M_SS"),
                         c("M_SS", "M_OS"),
                         c("F_OS", "M_OS"))
  plt <- ggplot(data, aes_string(x = "Group", y = mod, color = "Group",
                                 fill = "Group", alpha = 0.8)) +
    geom_violin(lwd = 0.5) +
    geom_point(position = position_dodge(width = 0.75),
               color = "slategrey", size = 2, alpha = 1) +
    # stat_compare_means(comparisons = my_comparisons,
                       # paired = FALSE, method = "wilcox.test") +

    scale_fill_manual(values = c("F_SS" = "mediumpurple",
                                 "F_OS" = "darkorchid4",
                                 "M_SS" = "lightseagreen",
                                 "M_OS" = "deepskyblue4")) +
    scale_color_manual(values = c("F_SS" = "mediumpurple",
                                  "F_OS" = "darkorchid4",
                                  "M_SS" = "lightseagreen",
                                  "M_OS" = "deepskyblue4")) +

    ylab("Module Score") +
    ggtitle(mod) +
    theme_classic() +
    theme(text = element_text(size = 40))
  print(plt)

  if (save == TRUE) {
    ggsave(paste(mod, "_violin.png"), plt,
           bg = "white",
           height = 8, width = 10,
           units = "in", device = "png")
  }
}


#' Merge Dataset and Metadata Based on Animal Column
#'
#' This function merges a primary dataset with a metadata dataset
#' based on a common 'animal' column.
#'
#' @param dataset Data frame representing the primary dataset,
#' containing module scores.
#' @param metadata Data frame representing the metadata.
#'
#' @return Returns a merged data frame.
#'
pairs <- function(dataset, metadata) {
  tryCatch({
    paired <- merge(dataset, metadata[c("animal", "pair")], on = "animal")
    return(paired)
  }, error = function(e) {
    cat("An error occurred: ", e$message, "\n")
    return(NULL)
  })
}


#' Create a Subset of a Dataset by Dropping Specified Columns
#'
#' This function takes a dataset and a vector of column names to drop,
#' then returns a subset of the dataset without those columns.
#'
#' @param dataset Data frame from which columns will be dropped.
#' @param drop_cols Character vector specifying names of the dropped columns.
#'
#' @return Returns subset of the input dataset with specified columns removed.
#'
cor_dataframe <- function(dataset, drop_cols) {
  if (!is.character(drop_cols)) {
    stop("'drop_cols' should be a character vector.")
  }
  if (ncol(dataset) == 0) {
    stop("'dataset' is empty.")
  }
  if (!all(drop_cols %in% names(dataset))) {
    stop("Some columns in 'drop_cols' do not exist in the dataset.")
  }
  tryCatch({
    core <- dataset[, !(names(dataset) %in% drop_cols)]
    return(core)
  }, error = function(e) {
    cat("An unexpected error occurred: ", e$message, "\n")
    return(NULL)
  })
}


#' Calculate Spearman Correlation Between Two Columns in a Dataset
#'
#' This function computes the Spearman correlation between two columns,
#' labeled 'O' and 'B', in a given dataset.
#'
#' @param data Data frame that must contain the columns 'O' and 'B'..
#'
#' @return Returns an object containing the Spearman correlation result
#' between the 'O' and 'B' columns of the dataset.

partner_correlation <- function(data) {
  tryCatch({
    cor <- rcorr(data$O, data$B, type = "spearman")
    return(cor)
  }, error = function(e) {
    cat("An unexpected error occurred: ", e$message, "\n")
    return(NULL)
  })
}


#' Extract the Second Element of the 'r' Component from a Correlation Object
#'
#' This function extracts the second element of the
#' 'r' (correlation coefficient) component from the correlation object
#' created in the previous function..
#'
#' @param cor Correlation object
#'
#' @return The second element of the 'r' component from the correlation object
#' as a numeric value.
#'
cor_r <- function(cor) {
  if (!is.list(cor) && !is.data.frame(cor)) {
    stop("Input must be a list or data frame.")
  }
  if (!"r" %in% names(cor)) {
    stop("Input does not have an 'r' component.")
  }
  if (length(cor$r) < 2) {
    stop("The 'r' component does not have at least two elements.")
  }
  if (is.null(cor$r[2]) || is.na(cor$r[2])) {
    stop("The second element of 'r' is NA or NULL.")
  }
  if (!is.numeric(cor$r)) {
    stop("The 'r' component is not numeric.")
  }
  return(cor$r[2])

  tryCatch({
    r <- cor$r[2]
    return(r)
  }, error = function(e) {
    cat("An error occurred: ", e$message, "\n")
  })
}

#' Extract the Second Element of the 'P' Component from a Correlation Object
#'
#' This function is designed to extract the second element of the 'P'
#' (p-value) component from correlation object created in a previous function.
#'
#' @param cor Correlation object
#'
#' @return Returns the second element of the 'P' component from the correlation
#' object as a numeric value.

cor_p <- function(cor) {
  if (!is.list(cor) && !is.data.frame(cor)) {
    stop("Input must be a list or data frame.")
  }
  if (!"P" %in% names(cor)) {
    stop("Input does not have a 'P' component.")
  }
  if (length(cor$P) < 2) {
    stop("The 'P' component does not have at least two elements.")
  }
  p <- cor$P[2]
  if (is.null(p) || is.na(p)) {
    stop("The second element of 'P' is NA or NULL.")
  }
  if (!is.numeric(p)) {
    if (is.character(p) || is.factor(p)) {
      p <- as.numeric(as.character(p))
    } else {
      stop("The 'P' component cannot be converted to numeric.")
    }
  }
  return(p)
}


#' Create and Save a Scatter Plot Showing Correlation Between Partners
#'
#' Generates a scatter plot that shows the correlation between partners.
#'
#' @param data Data frame containing the data to be plotted.
#' <ust contain columns named 'pair', 'color', 'O', 'B',
#' and variable representing 'mod'.
#'
#' @return Scatter plots as PNG files.
#' @details The function performs the following operations:
#'          - Transforms data into wide format using `pivot_wider` from 'tidyr'.
#'          - Calculates Spearman correlation coefficient and p-value between
#'            'O' and 'B' columns.
#'          - Creates scatter plot with points for each pair and linear model.
#'          - The plot title includes the correlation coefficient and p-value.
#'          - Plot is saved as PNG file with name based on the 'mod' variable.
corr_plot <- function(data, color = "black", font_size, save = FALSE) {
  data_wide <- tidyr::pivot_wider(data, id_cols = c("pair", "pair_type"),
                                  names_from = "color", values_from = mod)
  data <- data %>% filter(animal != "4918")
  min_val <- min(data[,mod])
  # print(min_val)
  max_val <- max(data[,mod])
  # print(max_val)
  
  cor <- partner_correlation(data_wide)
  
  r <- cor_r(cor)
  p <- cor_p(cor)
  
  lab <- paste("Rho = ", round(r, digits = 4),
               " p = ", round(p, digits = 4), sep = "")
  
  p <- ggplot(data_wide, aes(x = O, y = B, label = pair)) +
    geom_smooth(method = "lm", color = color, fill = color) +
    geom_point() +
    # scale_color_manual(values = c("FM" = "coral",
    #                               "FF"="slateblue",
    #                               "MM"="#1B9E77")) +
    ggtitle(paste0(mod, " ", lab)) +
    xlab("Partner 1 Avg. Gene Expr.") +
    ylab("Partner 2 Avg. Gene Expr.") +
    # coord_fixed(xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    theme(text = element_text(size = font_size))
  print(p)
  
  if (save == TRUE) {
    ggsave(paste(mod, "partner_corr_square.pdf"), p, bg = "white")
  }
  
  return(p)
}




corr_plot_bypair <- function(data, color = "black", font_size, save = FALSE) {
  data_wide <- tidyr::pivot_wider(data, id_cols = c("pair", "pair_type"),
                                  names_from = "color", values_from = mod)

  cor <- partner_correlation(data_wide)

  r <- cor_r(cor)
  p <- cor_p(cor)

  lab <- paste("R = ", round(r, digits = 4),
               " p = ", round(p, digits = 4), sep = "")

  p <- ggplot(data_wide, aes(x = O, y = B, label = pair)) +
    xlab("Partner 1") +
    ylab("Partner 2") +
    geom_smooth(method = "lm", color = color, fill = color) +
    geom_point(aes(color = pair_type)) +
    scale_color_manual(values = c("FM" = "coral",
                                  "FF"="slateblue",
                                  "MM"="#1B9E77")) +
    ggtitle(paste0(mod, " ", lab)) +
    theme(element_text(size = font_size)) +
    theme_classic()
  print(p)
  
  if (save == TRUE) {
    ggsave(paste(mod, "partner_corr_bypairtype.png"), p, bg = "white")
  }
  
  return(p)
}


#' Extract Unique Pair Values Excluding a Specified Value
#'
#' Extracts unique values from the 'pair' column in a dataset and
#' excludes a specified values.
#'
#' @param dataset Data frame that must contain a column named 'pair'.
#' @param exclusion Value to be excluded from the list of unique pairs.
#'
#' @return Returns a vector of unique values from the 'pair' column..
#'
use_pairs <- function(dataset, exclusion) {
  if (!"pair" %in% names(dataset)) {
    stop("Error: 'dataset' does not contain a 'pair' component.")
  }
  if (length(exclusion) != 1) {
    stop("Error: 'exclusion' should be a single value.")
  }
  pairs1 <- unique(dataset$pair)
  if (!exclusion %in% pairs1) {
    warning("Warning: 'exclusion' value not found in 'pairs1'.")
  }
  tryCatch({
    pairs2 <- pairs1[pairs1 != exclusion]
    return(pairs2)
  }, error = function(e) {
    cat("An error occurred: ", e$message, "\n")
  })
}


#' Create an Empty Data Frame for Storing Distance Data
#'
#' This function initializes an empty data frame for storing distance data
#' and associated attributes.
#'
#' @param pairdata A placeholder parameter for pair data
#' @param datatype A placeholder parameter for data type
#' @param realorfake A placeholder parameter for real or fake classification
#' @param distance A placeholder parameter for distance values
#' @param module A placeholder parameter for module information
#'
#' @return Returns an empty data frame with predefined columns
#'
make_dist_df <- function(pairdata, datatype, realorfake, distance, module) {
  tryCatch({
    all.dists <- data.frame(pairdata = character(), datatype = character(),
                            realorfake = character(), distance = numeric(),
                            module = character())
    return(all.dists)
  }, error = function(e) {
    cat("An error occurred: ", e$message, "\n")
  })
}


#' Generate Module Names in Order
#'
#' This function creates a sequence of module names,
#' numbered from 1 to a specified number.
#'
#' @param num_modules A numeric value specifying the number of modules to
#' generate names for, defined in a subsequent plotting function.
#'
#' @return Returns a vector of module names in the format 'Module.X',

generate_module_order <- function(num_modules) {
  if (!is.numeric(num_modules) || num_modules <= 0 ||
        num_modules != as.integer(num_modules)) {
    stop("Error: 'num_modules' must be a positive integer.")
  }
  tryCatch({
    module_names <- paste0("Module.", 1:num_modules)
    return(module_names)
  }, error = function(e) {
    cat("An error occurred: ", e$message, "\n")
  })
}


#' Calculate Rank Distances for Animal Pairs in a Dataset
#'
#' This function computes the rank distances between animal pairs in a given
#' dataset based on a specific module.
#' It assumes that 'mod' (a global variable) specifies the column
#' to be used for ranking.
#'
#' @param data A data frame containing the necessary columns for the
#' calculation, including 'animal', 'pair', and 'pair_type'.
#'
#' @return Returns a data frame containing
#' calculated rank distances for each animal pair.
#' @details Function performs the following operations:
#'          - Converts specified data to a matrix,
#'            computes ranks, and then calculates distances.
#'          - Sets row and column names of matrix based on animal names.
#'          - Checks for the presence of 'pair' and 'pair_type' columns.
#'          - Initializes a data frame to store distances.
#'          - Processes pairs and calculates distances.
get_rank_distances <- function(data) {
  if (!exists("mod")) {
    stop("'mod' is not defined.")
  }
  tryCatch({
    data_matrix <- as.matrix(data[, mod])
    rank_matrix <- rank(data_matrix)
    distance_matrix <- as.matrix(dist(rank_matrix))
  }, error = function(e) {
    stop("Error in matrix conversion or ranking: ", e$message)
  })
  tryCatch({
    rownames(distance_matrix) <- data$animal
    colnames(distance_matrix) <- rownames(distance_matrix)
  }, error = function(e) {
    stop("Error in setting row and column names: ", e$message)
  })
  if (!("pair" %in% names(data)) || !("pair_type" %in% names(data))) {
    stop("Both 'pair' and 'pair_type' columns must exist in 'data'.")
  }
  tryCatch({
    dist_df_rank <- make_dist_df(pair, type, fake_real, dist, Module)
  }, error = function(e) {
    stop("Error in initializing dist_df_rank: ", e$message)
  })
  tryCatch({
    pairs <- use_pairs(data, "4918x4967")
    for (pair in pairs) {
      ani1 <- str_split_i(pair, "x", 1)
      ani2 <- str_split_i(pair, "x", 2)

      if (!ani1 %in% rownames(distance_matrix) ||
            !ani2 %in% rownames(distance_matrix)) {
        next
      }
      distance <- distance_matrix[ani1, ani2]
      type <- data$pair_type[data$pair == pair][1]
      type <- paste(type, "real", sep = "_")
      mini_df <- data.frame("pair" = pair, "type" = type,
                            "fake_real" = "real",
                            "dist" = distance,
                            "Module" = mod)

      dist_df_rank <- rbind(dist_df_rank, mini_df)
    }
  }, error = function(e) {
    stop("Error in processing pairs: ", e$message)
  })

  return(dist_df_rank)
}


#' Generate and Save a Boxplot of Rank-Based Distances Across Modules
#'
#' This function creates a boxplot illustrating
#' rank-based distances for different modules in a dataset.
#' It highlights significant modules and draws a line
#' indicating the expected distance.
#'
#' @param data A data frame to be plotted
#'
#' @return Function generates and saves a boxplot as a PNG file.
#'
#' @details Function performs the following operations:
#'          - Sets order of modules using a generated sequence of module names.
#'          - Creates boxplot with points jittered for better visibility.
#'          - Highlights significant modules in a different color.
#'          - Adds horizontal line indicating the expected rank-based distance.
#'          - Saves plot as a PNG file named 'rank_based_partner_dists.png'.
plot_rank <- function(data, save = FALSE) {
  expected_dist <- mean(dist(1:38))
  num_modules <- 23
  mod_order <- generate_module_order(num_modules)
  data$Module <- factor(data$Module, levels = mod_order)
  plot <- ggplot(data, aes(x = Module, y = dist, fill = Module)) +
    geom_boxplot() + geom_point(position = position_jitter(w = 0.2, h = 0)) +
    geom_abline(slope = 0, intercept = expected_dist) + theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ylab("Rank-based distance") +
    scale_fill_manual(values = c("Module.1" = "gainsboro",
                                 "Module.2" = "gainsboro",
                                 "Module.3" = "gainsboro",
                                 "Module.4" = "gainsboro",
                                 "Module.5" = "gainsboro",
                                 "Module.6" = "deeppink4",
                                 #colored bc significant
                                 "Module.7" = "gainsboro",
                                 "Module.8" = "gainsboro",
                                 "Module.9" = "gainsboro",
                                 "Module.10" = "gainsboro",
                                 "Module.11" = "chartreuse4",
                                 #colored bc significant
                                 "Module.12" = "gainsboro",
                                 "Module.13" = "gainsboro",
                                 "Module.14" = "gainsboro",
                                 "Module.15" = "gainsboro",
                                 "Module.16" = "gainsboro",
                                 "Module.17" = "gainsboro",
                                 "Module.18" = "gainsboro",
                                 "Module.19" = "gainsboro",
                                 "Module.20" = "gainsboro",
                                 "Module.21" = "gainsboro",
                                 "Module.22" = "gainsboro",
                                 "Module.23" = "gainsboro"))
  print(plot)
  
  if (save == TRUE) {
    ggsave("rank_based_partner_dists.png", plot, width = 12, height = 8,
           units = "in", bg = "white", device = png)
  }

  return(plot)
}

plot_rank_group <- function(data, save = FALSE) {
  my_comparisons <- list(c("FM_real", "FF_real"), c("FM_real", "MM_real"), c("FF_real", "MM_real"))
  expected_dist <- mean(dist(1:38))
  plot <- ggplot(data, aes(x = type, y = dist, color = type, fill = type, alpha = 0.8)) +
    geom_violin(lwd = 0.5) +
    geom_point(color = "slategrey", alpha = 1) +
    geom_abline(slope = 0, intercept = expected_dist) +
    scale_color_manual(values = c("FM_real" = "coral",
                                  "FF_real"="slateblue",
                                  "MM_real"="#1B9E77")) +
    scale_fill_manual(values = c("FM_real" = "coral",
                                 "FF_real"="slateblue",
                                 "MM_real"="#1B9E77")) +
    ggtitle(mod) +
    ylab("Rank distance between partners") +
    # stat_compare_means(comparisons = my_comparisons, paired = FALSE, method = "wilcox.test") +
    theme_classic()
  
  print(plot)
  
  return(plot)
}

get_fake_pairs <- function(data) {
  #find all combos of not-true pairs
  all_animals <- data$animal
  all_animals <- all_animals[all_animals != c("4918")]
  pairs <- unique(data$pair)
  pairs <- pairs[pairs != "4918x4967"]
  l <- crossing(var1 = all_animals, var2= all_animals)
  l <- l %>% filter(!(var1==var2))
  
  #for not-true (fake) pairs - as a control!!!
  #remove true pairs from crossing matrix
  for (pair in pairs) {
    #find distance between pairs in matrix
    ani1 <- str_split_i(pair, "x", 1)
    ani2 <- str_split_i(pair, "x", 2)
    
    l <- l %>% filter(!(var1 == ani1 & var2 == ani2)) %>% filter(!(var1==ani2 & var2==ani1))
    
  }
  
  #remove reversed pairs
  fltr <- !duplicated(apply(l, 1, function(x) paste0(sort(x), collapse = "")))
  l <- l[fltr,]
  
  return(l)
}

get_euclidean_dists <- function(seurobj, mod_genes, data) {
  print(mod)
  mod_num <- strsplit(mod, "Module.")[[1]][2]
  mgs <- mod_genes$Gene[mod_genes$Module == mod_num]
  
  Idents(SCT_norm) <- "Ani"
  mod_expression <- AverageExpression(SCT_norm, features = mgs, assays = "SCT", slot = "counts")
  mod_expression <- mod_expression[1] %>% as.data.frame()
  mod_expression <- t(mod_expression)
  rownames(mod_expression) <- sub("SCT.", "", rownames(mod_expression))
  mod_expression <- mod_expression %>% as.data.frame()
  mod_expression <- mod_expression[!(row.names(mod_expression) %in% "4918"),]
  
  #get euclidean distances between all possible animal pairs
  mod_dist <- dist(mod_expression) %>% as.matrix()

  #distance between true pairs
  pairs <- unique(data$pair)
  pairs <- pairs[pairs!="4918x4967"]
  
  real_df <- data.frame(pair = character(), type = character(), fake_real = character(), dist = numeric())
  colnames(real_df) <- c("pair", "dist")
  for (pair in pairs) {
    ani1 <- str_split_i(pair, "x", 1)
    ani2 <- str_split_i(pair, "x", 2)
    
    distance <- mod_dist[ani1, ani2]
    type <- data$pair_type[data$pair==pair][1]
    type <- paste(type, "real", sep = "_")
    mini.df <- data.frame("pair" = pair, "type" = type, "fake_real" = "True Pair", "dist" = distance)
    
    real_df <- rbind(real_df, mini.df)
  }
  
  #for fake pairs
  l <- get_fake_pairs(merged_data)
  
  fake_list <- list()
  for (i in (1:nrow(l))) {
    row <- l[i,]
    ani1 <- row[1] %>% as.character()
    ani2 <- row[2] %>% as.character()
    pairname <- paste(ani1, ani2, sep = "x")
    distance <- mod_dist[ani1, ani2]
    revpair <- paste(ani2, ani1, sep = "x")
    if (revpair %in% names(fake_list)) {
    } else {
      fake_list[pairname] <- distance
    }
    
  }
  
  fake_list <- fake_list %>% unlist()
  
  fake_df <- data.frame(fake_list)
  fake_df$pair <- rownames(fake_df)
  fake_df$type <- "all_fake"
  fake_df$fake_real <- "Chance"
  colnames(fake_df)[1] ="dist"
  
  
  # modo_dists <- all.dists.df %>% filter(Module == "Module.6")
  wil <- wilcox.test(fake_df$dist, real_df$dist)
  print(wil)
  
  real_df$fr <- "True Pair"
  fake_df$fr <- "Chance"
  
  all_fr_dists <- rbind(real_df, fake_df)
  all_fr_dists$Module <- mod
  
  return(all_fr_dists)
}

plot_euclidean_dists <- function(euclidean_data, save = FALSE) {
  mod_df <- euclidean_data %>% filter(Module == mod)
  if (mod == "Module.6") {
    plt_color = "deeppink4"
  }
  else if (mod == "Module.11") {
    plt_color = "chartreuse4"
  }
  else {
    plt_color = "firebrick4"
  }
  p <- ggplot(mod_df, aes(x = fr, y = dist, color = fr, fill = fr, alpha = 0.8)) +
    geom_boxplot() +
    # geom_jitter(color = "slategrey", alpha = 0.8, width = 0.2, height = 0) +
    stat_compare_means() +
    scale_fill_manual(values = c("Chance" = "grey69",
                                 "True Pair" = plt_color)) +
    scale_color_manual(values = c("Chance" = "grey69",
                                  "True Pair" = plt_color)) +
    ggtitle(paste(mod, "Euclidean Distance")) +
    xlab("") +
    ylab("Euclidean Distance") +
    theme_classic()
  
  print(p)
  
  if (save == TRUE) {
    ggsave(paste(mod, "euclidean dist.pdf"), bg = "white", width = 8, height = 8, units = "in", device = pdf)
  }
  
  return(p)
  

}
  
plot_euclidean_dists_pairtype <- function(euclidean_data, save = FALSE) {
  mod_df <- euclidean_data %>% filter(Module == mod & fake_real == "True Pair")
  my_comparisons <- list(c("FM_real", "FF_real"), c("FM_real", "MM_real"), c("FF_real", "MM_real"))
  p <- ggplot(mod_df, aes(x = type, y = dist, color = type, fill = type, alpha = 0.8)) +
    geom_violin() +
    geom_jitter(color = "slategrey", alpha = 1, width = 0.2, height = 0) +
    # stat_compare_means(comparisons = my_comparisons) +
    scale_color_manual(values = c("FM_real" = "coral",
                                  "FF_real"="slateblue",
                                  "MM_real"="#1B9E77")) +
    scale_fill_manual(values = c("FM_real" = "coral",
                                 "FF_real"="slateblue",
                                 "MM_real"="#1B9E77")) +
    ggtitle(paste(mod, "Euclidean Distance Between Partners")) +
    xlab("Pair Type") +
    ylab("Euclidean Distance") +
    theme_classic()
  
  print(p)
  return(p)
}

