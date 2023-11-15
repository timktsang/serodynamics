#' Example input data
#'
#' This is an example of the input data used in the \code{serodynamics} function. This data frame illustrates the format of the input data.
#' @docType data
#' @usage data(inputdata)
#' @format A data frame with 9 variables, where each row represents an individual:
#' \describe{
#'   \item{age_group}{0: children, 1: adults, 2: older adults}
#'   \item{start_time}{start of follow-up}
#'   \item{end_time}{end of follow-up}
#'   \item{time1}{date of first serum collection}
#'   \item{time2}{date of second serum collection}
#'   \item{time3}{date of third serum collection}
#'   \item{HAI_titer_1}{HAI titer for first serum collection}
#'   \item{HAI_titer_2}{HAI titer for second serum collection}
#'   \item{HAI_titer_3}{HAI titer for third serum collection}
#' }
#' @family inputdata
"inputdata"

#' Example flu activity data
#'
#' This is an example of the flu activity data used in the \code{serodynamics} function. This data frame specifies the format of the flu activity data.
#' @docType data
#' @usage data(flu_activity)
#' @format A data frame with 1 variable, where each row represents a date, and it should match the date in the input data:
#' \describe{
#'   \item{h1.activity}{This is the influenza activity from surveillance data. It can be on a relative scale, as the model includes a scale parameter to estimate infection probability.}
#'}
#' @family example_data
"flu_activity"
