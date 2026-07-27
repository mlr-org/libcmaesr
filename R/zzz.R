#' @import checkmate
#' @import mlr3misc
#' @importFrom stats runif
"_PACKAGE"


#' @title CMAES Algorithm Names
#' @description A vector of strings containing the names of the CMAES variants.
#' See \url{https://cma-es.github.io/libcmaes/doc/html/classlibcmaes_1_1CMAParameters.html} for details.
#' @return A `character` vector of algorithm names.
#' @examples
#' cmaes_algos
#' @name cmaes_algos
#' @export
cmaes_algos = c(
  "cmaes",
  "ipop",
  "bipop",
  "acmaes",
  "aipop",
  "abipop",
  "sepcmaes",
  "sepipop",
  "sepbipop",
  "sepacmaes",
  "sepaipop",
  "sepabipop",
  "vdcma",
  "vdipopcma",
  "vdbipopcma"
)
