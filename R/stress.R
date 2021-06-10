#' Fits a linear function to describe a stress response
stress_linear <- function() {
  function(srel) {
    1 - srel
  }
}

#' Fits a convex function to describe a stress response
stress_convex <- function(fshape) {
  function(srel) {
    1 - (exp(fshape * srel) - 1) / (exp(fshape) - 1)
  }
}

#' Fits a logistic function to describe a stress response
stress_logistic <- function(sn, sx) {
  r <- -log(2 * sn * (sx - 0.5) / (sx - sn)) * 2
  kmax <- (sn * sx) / (sn + (sx - sn) * exp(-r))
  kmin <- sn
  function(srel) {
    ks <- (sn * sx) / (sn + (sx - sn) * exp(-r * (1 - srel)))
    (ks - kmin) / (kmax - kmin)
  }
}

#' Stress response based on a spline function
stress_spline <- function(known_srel, known_ks, method = "hyman", ties = mean) {
  fun <- splinefun(x = known_srel, y = known_ks, method = method, ties = ties)
  function(srel) {
    fun(x = srel)
  }
}
