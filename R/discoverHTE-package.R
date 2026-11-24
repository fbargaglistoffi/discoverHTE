#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @rawNamespace import(stats, except = D)
#' @import collapse
#' @import ggplot2
#' @importFrom splines ns
## usethis namespace: end
NULL

# `.c(Y, A, X) %=% process_formula(...)` is collapse's multiple assignment. It creates
# these bindings at run time, which codetools cannot see, so `R CMD check` reports them
# as undefined globals in every method_*() closure.
utils::globalVariables(c("Y", "A", "X"))
