#' Estimate an IV-optimal individualized treatment rule
#'
#' \code{IV_PILE} estimates an IV-optimal individualized treatment
#' rule given a dataset with estimated partial identification intervals
#' for each instance.
#'
#'
#' @param dt A dataframe whose first column is a binary IV 'Z', followed
#'           by q columns of observed covariates, a binary
#'           treatment indicator 'A', a binary outcome 'Y',
#'           lower endpoint of the partial identification interval 'L',
#'           and upper endpoint of the partial identification interval 'U'.
#'           The dataset has q+5 columns in total.
#'
#' @param kernel The kernel used in the weighted SVM algorithm. The user
#'               may choose between 'linear' (linear kernel) and
#'               'radial' (Gaussian RBF kernel).
#'
#' @param C Cost of violating the constraint. This is the parameter C in
#'          the Lagrange formulation.
#'
#' @param sig Sigma in the Gaussian RBF kernel. Default is set to
#'            1/dimension of covariates, i.e., 1/q. This parameter
#'            is not relevant for linear kernel.
#'
#' @return An object of the type \code{wsvm}, inheriting from \code{svm}.
#'
#' @examples
#' \dontrun{
#' # It is necessary to install the package locClass in order
#' # to run the following code.
#'
#' attach(dt_Rouse)
#' # Construct an IV out of differential distance to two-year versus
#' # four-year college. Z = 1 if the subject lives not farther from
#' # a 4-year college compared to a 2-year college.
#' Z = (dist4yr <= dist2yr) + 0
#'
#' # Treatment A = 1 if the subject attends a 4-year college and 0
#' # otherwise.
#' A = 1 - twoyr
#'
#' # Outcome Y = 1 if the subject obtained a bachelor's degree
#' Y = (educ86 >= 16) + 0
#'
#' # Prepare the dataset
#' dt = data.frame(Z, female, black, hispanic, bytest, dadsome,
#'      dadcoll, momsome, momcoll, fincome, fincmiss, A, Y)
#'
#' # Estimate the Balke-Pearl bound by estimating each constituent
#' # conditional probability p(Y = y, A = a | Z, X) with a multinomial
#' # regression.
#' dt_with_BP_bound_multinom = estimate_BP_bound(dt, method = 'multinom')
#'
#' # Estimate the IV-optimal individualized treatment rule using a
#' # linear kernel, under the putative IV and the Balke-Pearl bound.
#'
#'
#' iv_itr_BP_linear = IV_PILE(dt_with_BP_bound_multinom, kernel = 'linear')
#'}
#'
#' @importFrom rlang .data
#' @export
#'
#'
#'
IV_PILE <- function(dt, kernel = 'linear',
                    C = 1, sig = 1/(ncol(dt) - 5)){

  if (!requireNamespace("locClass", quietly = TRUE)) {
    stop("Package \"locClass\" needed for this function to work.
         Please install it.",
         call. = FALSE)
  }

  n = dim(dt)[1]

  # Create labels
  ind_p = which(dt$L > 0)
  ind_n = which(dt$U < 0)
  ind_sp = which((dt$L <= 0) & (dt$U >= 0) & (abs(dt$U) >= abs(dt$L)))
  ind_sn = which((dt$L <= 0) & (dt$U >= 0) & (abs(dt$U) < abs(dt$L)))
  ind_semi = c(ind_sp, ind_sn)

  labels = numeric(n)
  labels[ind_p] = 1
  labels[ind_sp] = 1
  labels = as.factor(labels > 0)

  # Create weights
  weights = numeric(n)
  weights[ind_p] = abs(dt$U[ind_p])
  weights[ind_n] = abs(dt$U[ind_n])
  weights[ind_semi] = abs(abs(dt$L[ind_semi]) - abs(dt$U[ind_semi]))

  # Run weighted SVM
  dt_md = dplyr::select(dt, -.data$Z, -.data$A, -.data$Y, -.data$L, -.data$U)
  iv_pile = locClass::wsvm(labels ~ ., data = dt_md,
                           type = 'C-classification',
                           kernel = kernel,
                           case.weights = weights,
                           cost = C,
                           gamma = sig,
                           fitted = TRUE)
  return(iv_pile)
}
