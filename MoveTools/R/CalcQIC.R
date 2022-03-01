#' Calculate Quasi-likelihood under the independence model criterion (QIC)
#'
#' Takes a fit clogit() model with a cluster term and provides the QIC based on methods of Pan et al. 2001. Written by Mathieu Basille and Jerod Merkle. Last updated January 2022.
#'
#' @param mod A fitted clogit() model including a cluster term from the survival package.
#' @param details If TRUE, provides all the details of QIC calculation
#'
#' @return   Returns a dataframe providing the QIC results.
#'
#' @examples
#' # none

#' @export


CalcQIC <- function(mod,
                    details = FALSE
                    ){
  if (!exists("naive.var", mod))
    stop("QIC can be computed only if robust variances are estimated.")
  trace <- sum(diag(solve(mod$naive.var) %*% mod$var))
  quasi <- mod$loglik[2]
  if (details)
    return(data.frame(QICR = -2 * quasi + 2 * trace, QICI = -2 *
                        quasi + 2 * length(mod$coefficients), QuasiLL = quasi,
                      n = mod$n, nevent = mod$nevent, K = length(mod$coefficients),
                      Trace = trace))
  else return(-2 * quasi + 2 * trace)
}
