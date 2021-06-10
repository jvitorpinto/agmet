ext_rad_raster <- function(t0, t1, j, rst) {
  sc <- equation_of_time(j)
  tm <- (ifelse(t0 > t1, (24 + t1 + t0), (t0 + t1)) / 2) %% 24
  tdif <- ifelse(t0 > t1, 24 + t1 - t0, t1 - t0)
}
