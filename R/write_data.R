#' Write TiNDA data
#'
#' Function to write TiNDA data to a tsv file with TiN_Class information
#'
#' @param tinda_object TiNDA object produced by the main TiNDA function
#' @param file_name Name of the file
#'
#' @importFrom readr write_tsv
#' @export
write_data <- function(tinda_object, file_name){
  write_tsv(tinda_object$data, file_name)
}

#' Run TiNDA pipeline from input file to output file
#'
#' Convenience function that reads input data, runs TiNDA analysis,
#' and writes results to output file in a single call.
#'
#' @param input_file Path to input TSV file with variant data
#' @param output_file Path to output TSV file for results
#' @param sample_name Sample name for the analysis title. If NULL, uses input file basename.
#' @param data_source Data source type ('WGS' or 'WES'). Default: 'WGS'
#' @param ... Additional arguments passed to TiNDA()
#'
#' @return TiNDA object (invisibly) and writes results to output_file
#'
#' @examples
#' \dontrun{
#' run_pipeline("variants.tsv", "results.tsv", sample_name = "sample_1")
#' }
#'
#' @importFrom readr read_tsv
#' @export
run_pipeline <- function(input_file, output_file, sample_name = NULL,
                         data_source = "WGS", ...) {
  # Read input data
  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file)
  }

  tbl <- read_tsv(input_file, show_col_types = FALSE)

  # Set sample name from file if not provided
  if (is.null(sample_name)) {
    sample_name <- tools::file_path_sans_ext(basename(input_file))
  }

  # Run TiNDA analysis
  result <- TiNDA(tbl, sample_name = sample_name, data_source = data_source, ...)

  # Write results
  write_data(result, output_file)

  cat("Results written to:", output_file, "\n")

  invisible(result)
}
