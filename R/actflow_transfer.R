transfer_function <- function(activity, transfer = c("linear", "relu", "sigmoid", "logit"), threshold = 0, a = 1) {
  transfer <- match.arg(transfer)
  if (transfer == "linear") return(activity)
  if (transfer == "relu") return(activity * (activity > threshold))
  if (transfer == "sigmoid") return(1 / (1 + exp(-activity)))
  if (transfer == "logit") return((1 / a) * log(activity / (1 - activity)))
  activity
}
