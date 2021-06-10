day_mth <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
mid_day <- c(15, 43, 74, 104, 135, 165, 196, 227, 257, 288, 318, 349)

#' Calcula a evapotranspiração potencial (ETP) por meio do método de Thornthwaite e Mather
#' @param  temp a vector containing 12 values for air temperature or a raster::RasterStack
#' with 12 layers of air temperature.
#' @param ... other arguments
#' @param lat latitude
#' @export
thornthwaite_pet <- function(temp, ...) {
  UseMethod("thornthwaite_pet", temp)
}

#' @export
thornthwaite_pet.default <- function(temp, lat) {
  stopifnot(length(temp) == 12)
  heat_index = sum((temp / 5) ** 1.514)
  alpha <- 675e-9 * heat_index ** 3 - 771e-7 * heat_index ** 2 + 1792e-5 * heat_index + 0.49239
  pet <- 16 * ((10 * temp / heat_index) ** alpha)
  dl <- day_length(mid_day, lat)
  pet * day_mth * dl / 360
}

#' @export
thornthwaite_pet.RasterStack <- function(temp) {
  raster::calc(raster::stack(temp, raster::init(temp, "y")), fun = function(x){
    temp <- x[1:12]
    lat <- x[13]
    thornthwaite_pet.default(temp, lat)
  })
}

#' Calculates the water balance using the method of Thornthwaite & Mather
#' @description Calculates the water balance using the method of Thornthwaite & Mather
#' @param p precipitation
#' @export
water_balance <- function(p, pet, taw, ..., aw_initial = taw, max_iter = 256, max_error = 0, sequential = FALSE) {
  UseMethod("water_balance", p)
}

#' Calculates the water balance using the method proposed by Thornthwaite & Mather.
#' @param p precipitation
#' @param pet potential (maximum) evapotranspiration
#' @param taw total available soil water
#' @param aw_initial soil water content at the beginning of calculations.
#' @param max_iter the maximum number of iterations to perform until it stops (see details). Only valid if sequential is set to FALSE.
#' @param max_error the maximum allowed difference between the previous and the current calculation (see details). Only valid if sequential is set to FALSE.
#' @param sequential whether it will calculate a sequential or a cyclical water balance
#' @details The function iterates over each value of p-pet calculating the available
#' soil water (aw) at each step. If sequential is set FALSE and the last value of p-pet is reached,
#' the function returns to
#' the first value of p-pet and calculates a new value for aw. The new aw is compared to the aw calculated previously.
#' It will stop only when the new aw is equal to the old aw, or when the number of iterations exceeds max_iter
#' or the difference between the new aw and the old aw is lesser than max_error. If sequential is set TRUE,
#' then max_error and max_iter values are irrelevant.
#' @export
water_balance.default <- function(p, pet, taw, aw_initial = taw, max_iter = 256, max_error = 0, sequential = FALSE) {
  stopifnot(length(p) == length(pet))

  len <- length(p)

  dif <- p - pet
  aw = numeric(0)
  aw_var = numeric(0)

  i <- 0
  error <- Inf

  previous_aw <- aw_initial

  finished <- len == 0

  while (!finished) {
    i <- i + 1
    n <- if (sequential) i else 1 + ((i - 1) %% len)
    y <- ifelse(dif[[n]] >= 0, min(previous_aw + dif[[n]], taw), previous_aw * exp(dif[[n]] / taw))
    if (!sequential & i > length(p)) {
      error <- abs(y - aw[[n]])
    }
    if ((sequential & (i == len)) | (!sequential & (error <= max_error | i == max_iter))) {
      finished <- TRUE
    }
    aw[[n]] <- y
    aw_var[[n]] <- aw[[n]] - previous_aw
    previous_aw <- aw[[n]]
  }

  ret <- list(
    taw = taw,
    seq = sequential,
    wb = data.frame(
      step = 1:len,
      p = p,
      pet = pet,
      dif = dif,
      aw = aw,
      aw_var = aw_var,
      et = pmin(pet, p + abs(aw_var)),
      surplus = dif - aw_var
    )
  )
  class(ret) <- "WaterBalance"
  ret
}

#' @export
water_balance.RasterStack <- function(p, pet, taw, aw_initial = taw, max_iter = 256, max_error = 0, sequential = FALSE, verbose = 1) {
  len <- raster::nlayers(p)
  max_taw <- raster::maxValue(taw)

  dif <- p - pet
  aw = list()
  aw_var = list()

  i <- 0
  error <- Inf

  previous_aw <- aw_initial

  finished <- len == 0

  while (!finished) {
    i <- i + 1
    n <- if (sequential) i else 1 + ((i - 1) %% len)
    y <- (dif[[n]] >= 0) * min(previous_aw + dif[[n]], taw) + (dif[[n]] < 0) * previous_aw * exp(dif[[n]] / taw)
    if (!sequential & i > len) {
      error <- abs(raster::minValue(y - aw[[n]]))
    }
    if ((sequential & (i == len)) | (!sequential & (error <= max_error | i == max_iter))) {
      finished <- TRUE
    }

    if (verbose > 0) {
      if (sequential) {
        cat(paste0("mode: sequential; iteration: ", i, "\n"))
      } else {
        error_label <- if (error == Inf) "NA" else error
        cat(paste0("mode: cyclical; iteration: ", i, "; current error: ", error_label, "\n"))
      }
      if (verbose > 1){
        raster::plot(y, zlim = c(0, max_taw))
      }
    }

    aw[[n]] <- y
    aw_var[[n]] <- aw[[n]] - previous_aw
    previous_aw <- aw[[n]]
  }

  aw <- raster::stack(aw)
  aw_var <- raster::stack(aw_var)

  ret <- list(
    taw = taw,
    seq = sequential,
    wb = list(
      step = 1:len,
      p = p,
      pet = pet,
      dif = dif,
      aw = aw,
      aw_var = aw_var,
      et = raster::stack(lapply(1:len, function(n) min(pet[[n]], p[[n]] + abs(aw_var[[n]])))),
      surplus = dif - aw_var
    )
  )
  class(ret) <- c("WaterBalanceRaster", "WaterBalance")
  ret
}
