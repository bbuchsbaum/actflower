pearson_vec <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- stats::complete.cases(x, y)
  if (sum(ok) < 2) {
    return(NA_real_)
  }
  xx <- x[ok]
  yy <- y[ok]
  if (stats::sd(xx) == 0 || stats::sd(yy) == 0) {
    return(NA_real_)
  }
  stats::cor(xx, yy)
}

r2_vec <- function(y_true, y_pred) {
  y_true <- as.numeric(y_true)
  y_pred <- as.numeric(y_pred)
  ok <- stats::complete.cases(y_true, y_pred)
  if (sum(ok) < 2) {
    return(NA_real_)
  }
  yt <- y_true[ok]
  yp <- y_pred[ok]
  ss_res <- sum((yt - yp)^2)
  ss_tot <- sum((yt - mean(yt))^2)
  if (ss_tot == 0) {
    return(NA_real_)
  }
  1 - ss_res / ss_tot
}

mae_vec <- function(y_true, y_pred) {
  y_true <- as.numeric(y_true)
  y_pred <- as.numeric(y_pred)
  ok <- stats::complete.cases(y_true, y_pred)
  if (sum(ok) < 1) {
    return(NA_real_)
  }
  mean(abs(y_true[ok] - y_pred[ok]))
}
