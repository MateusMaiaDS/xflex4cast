#' @importFrom data.table data.table :=
#' @importFrom stats predict quantile sigma
#' @importFrom utils data flush.console
NULL

#' A simple progress bar for package
#'
#' @param current numeric for the current iteration.
#' @param total total number of iterations.
#' @param bar_length length of the bart displayed in the console.
#'
#' @details Internal use.
#'
sr_progress_bar <- function(current, total, bar_length = 50) {
  # Calculate the percentage of completion
  progress <- current / total
  percent <- round(progress * 100)

  # Calculate the number of characters for the progress bar
  num_hashes <- round(progress * bar_length)
  num_spaces <- bar_length - num_hashes

  # Create the progress bar string
  progress_bar <- paste0("[", paste(rep("#", num_hashes), collapse = ""), paste(rep(" ", num_spaces), collapse = ""), "]")

  # Display the progress bar with the percentage and message
  cat(sprintf("\r%s %d%% - Running the %d-th out %d observation", progress_bar, percent, current, total))

  # Flush the console output
  flush.console()
}

#' Convert prediction table into day-by-quantile matrix
#'
#' @param dt A data.table with columns: quantile, valid_time, value
#' @return A numeric matrix with rows = valid_time and columns = quantile
quantile_matrix <- function(dt) {

  x <- copy(dt)

  # Cast with rows = valid_time and columns = quantile
  wide <- dcast(
    x,
    valid_time ~ quantile,
    value.var = "value"
  )

  # Extract matrix
  vt <- wide$valid_time
  M <- as.matrix(wide[, -1])
  rownames(M) <- as.character(vt)

  return(M)
}
