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
