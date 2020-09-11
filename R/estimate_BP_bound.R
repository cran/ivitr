#' Estimate the Balke-Pearl bound for each instance in a dataset
#'
#' \code{estimate_BP_bound} estimates the Balke-Pearl bound for
#' each instance in the input dataset with a binary IV, observed
#' covariates, a binary treatment indicator, and a binary outcome.
#'
#'
#' @param dt A dataframe whose first column is a binary IV 'Z', followed
#'           by q columns of observed covariates, followed by a binary
#'           treatment indicator 'A', and finally followed by a binary
#'           outcome 'Y'. The dataset has q+3 columns in total.
#'
#' @param method A character string indicator the method used to estimate
#'               each constituent conditional probability of the
#'               Balke-Pearl bound. Users can choose to fit multinomial
#'               regression by setting method = 'multinom', and random
#'               forest by setting method = 'rf'.
#'
#' @param nodesize Node size to be used in a random forest algorithm if
#'                 method is set to 'rf'. The default value is set to 5.
#'
#' @return The original dataframe with two additional columns: L and U.
#'         L indicates the Balke-Pearl lower bound and U is the Balke-Pearl
#'         upper bound.
#'
#' @examples
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
#' # Calculate the Balke-Pearl bound by estimating each constituent
#' # conditional probability p(Y = y, A = a | Z, X) with a random
#' # forest.
#' dt_with_BP_bound_rf = estimate_BP_bound(dt, method = 'rf', nodesize = 5)
#'
#' # Calculate the Balke-Pearl bound by estimating each constituent
#' # conditional probability p(Y = y, A = a | Z, X) with a multinomial
#' # regression.
#' dt_with_BP_bound_multinom = estimate_BP_bound(dt, method = 'multinom')
#'
#' @importFrom stats predict
#' @importFrom rlang .data
#' @export
#'

estimate_BP_bound <- function(dt, method = 'rf', nodesize = 5){

  # Re-code 2*2 = 4 different (A, Y) combinations as 4 categorical
  # outcomes.
  dt$Outcome = 1*(dt$Y == 0 & dt$A == 0) + 2*(dt$Y == 0 & dt$A == 1) +
    3*(dt$Y == 1 & dt$A == 0) + 4*(dt$Y == 1 & dt$A == 1)

  # Prepare a dataset for estimating conditional probabilities
  dt_md = dplyr::select(dt, -.data$A, -.data$Y)

  # Prepare two datasets for predicting P(Y = y, A = a | Z = 0, X)
  # and P(Y = y, A = a | Z = 1, X)
  dt_predict_z_0 = dt_md
  dt_predict_z_0$Z = 0
  dt_predict_z_1 = dt_md
  dt_predict_z_1$Z = 1

  # Use the random forest method
  if (method == 'rf'){
    md = randomForest::randomForest(as.factor(Outcome) ~.,
                                    data = dt_md, nodesize = nodesize)

    # Predict P(Y = y, A = a | Z = 0, X)
    prob_res_z_0 = predict(md, newdata = dt_predict_z_0, type = 'prob')

    # Predict P(Y = y, A = a | Z = 1, X)
    prob_res_z_1 = predict(md, newdata = dt_predict_z_1, type = 'prob')
  }
  else if (method == 'multinom'){
    # Fit a multinomial logit regression
    md = nnet::multinom(Outcome ~., data = dt_md, trace = FALSE)
    prob_res_z_0 = predict(md, newdata = dt_predict_z_0,  'probs')
    prob_res_z_1 = predict(md, newdata = dt_predict_z_1,  'probs')
  } else
  return('Need to specify one of the following two methods:
                rf or multinom')

  # Balke-Pearl lower bound
  l_1 = prob_res_z_0[,1] + prob_res_z_1[,4] - 1
  l_2 = prob_res_z_1[,1] + prob_res_z_1[,4] - 1
  l_3 = prob_res_z_0[,4] + prob_res_z_1[,1] - 1
  l_4 = prob_res_z_0[,1] + prob_res_z_0[,4] - 1
  l_5 = 2*prob_res_z_0[,1] + prob_res_z_0[,4] + prob_res_z_1[,3] + prob_res_z_1[,4] - 2
  l_6 = prob_res_z_0[,1] + 2*prob_res_z_0[,4] + prob_res_z_1[,1] + prob_res_z_1[,2] - 2
  l_7 = prob_res_z_0[,3] + prob_res_z_0[,4] + 2*prob_res_z_1[,1] + prob_res_z_1[,4] - 2
  l_8 = prob_res_z_0[,1] + prob_res_z_0[,2] + prob_res_z_1[,1] + 2*prob_res_z_1[,4] - 2

  bk_lower_bound = pmax(l_1, l_2, l_3, l_4, l_5, l_6, l_7, l_8)

  # Balke-Pearl upper bound
  u_1 = 1 - prob_res_z_0[,3] - prob_res_z_1[,2]
  u_2 = 1 - prob_res_z_0[,2] - prob_res_z_1[,3]
  u_3 = 1 - prob_res_z_0[,2] - prob_res_z_0[,3]
  u_4 = 1 - prob_res_z_1[,2] - prob_res_z_1[,3]
  u_5 = 2 - 2*prob_res_z_0[,2] - prob_res_z_0[,3] - prob_res_z_1[,3] - prob_res_z_1[,4]
  u_6 = 2 - prob_res_z_0[,2] - 2*prob_res_z_0[,3] - prob_res_z_1[,1] - prob_res_z_1[,2]
  u_7 = 2 - prob_res_z_0[,3] - prob_res_z_0[,4] - 2*prob_res_z_1[,2] - prob_res_z_1[,3]
  u_8 = 2 - prob_res_z_0[,1] - prob_res_z_0[,2] - prob_res_z_1[,2] - 2*prob_res_z_1[,3]

  bk_upper_bound = pmin(u_1, u_2, u_3, u_4, u_5, u_6, u_7, u_8)

  dt$L = bk_lower_bound
  dt$U = bk_upper_bound
  dt = dplyr::select(dt, -.data$Outcome)

  return(dt)
}
