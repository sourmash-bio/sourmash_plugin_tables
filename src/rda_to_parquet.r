#! /usr/bin/env Rscript

# wget https://github.com/waldronlab/curatedMetagenomicData/raw/refs/heads/devel/data/sampleMetadata.rda

library(argparse)

if (!requireNamespace("polars", quietly = TRUE)) {
    cat("Polars package not found. Installing...\n")
    # Set CRAN mirror (if not already set)
    options(repos = c(CRAN = "https://cran.rstudio.com"))
    install.packages("codetools")
    Sys.setenv(NOT_CRAN = "true")
    install.packages("polars", repos = "https://community.r-multiverse.org")
}

library(polars)

p <- ArgumentParser(description = "Load an .rda file and save its contents as CSV or Parquet.")

p$add_argument("input_file", help = "Path to the .rda file to load")
p$add_argument("output_file", help = "Path to save the output file. Must be .csv or .parquet")
p$add_argument("-v", "--verbose", action="store_true", help = "Display full Dataframe")

args <- p$parse_args()

# Validate all file extension for 
input_ext <- tools::file_ext(args$input_file)
if (!input_ext == "rda") {
    stop("Input file must have .rda extension.")
} else { cat(sprintf("Processing %s input file.\n", args$input_file)) }

data_env <- new.env()  # Create a new environment to store the loaded objects
load(args$input_file, envir = data_env)

# Check for exactly one object in the .rda file
objects <- ls(data_env)
if (length(objects) != 1) {
    stop(sprintf("Expected exactly one object in the .rda file, but found %d: %s", length(objects), paste(objects, collapse = ", ")))
} else { cat(sprintf("Successfully loaded %s input file. \n", paste(objects))) }

data <- data_env[[objects[1]]]

# Ensure data is a data frame, and convert if necessary
if (!inherits(data, "data.frame")) {
    if (inherits(data, "list")) {
        data <- as.data.frame(data)
    } else {
        stop("The loaded object is neither a data frame nor a convertible list.")
    }
}

# Check the structure of the data for any unsupported types
#str(data)

output_ext <- tools::file_ext(args$output_file)
if (!output_ext %in% c("csv", "parquet")) {
    stop("Output file must have either .csv or .parquet extension.")
} else { cat(sprintf("Saving to %s output file. \n", args$output_file)) }

pl_df <- polars::as_polars_df(data)

if (args$verbose) {
    Sys.setenv(POLARS_FMT_MAX_COLS = -1)
    cat("Verbose output:\n")
    cat("DataFrame:\n")
    print(pl_df)
    cat("\nDataFrame Description:\n")
    print(pl_df$describe())
    cat("\nDataFrame Dtypes:\n")
    print(pl_df$dtypes)
    cat("\nDataFrame Column Names:\n")
    print(pl_df$columns)
} else {
    cat("DataFrame:\n")
    print(pl_df)
    cat("\nDataFrame Description:\n")
    print(pl_df$describe())
    cat("\nDataFrame Column Names:\n")
    print(pl_df$columns)
}

# Get unique values from the 'disease' column
unique_diseases <- pl_df$disease$unique()

# Print the unique values
print(unique_diseases)

# Filter the DataFrame by strings in the 'disease' column
#filtered_df <- pl_df$filter(pl_df$disease$is_in(c("CRC", "adenoma")))


if (output_ext == "csv") {
    pl_df$write_csv(args$output_file)
} else if (output_ext == "parquet") {
    pl_df$write_parquet(args$output_file)
}

cat(sprintf("File saved as %s.\n", args$output_file))
