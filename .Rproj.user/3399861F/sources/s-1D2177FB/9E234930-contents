limit_range <- function(x, x_min, x_max) {
  UseMethod("limit_range", x)
}

limit_range.default <- function(x, x_min, x_max) {
  pmax(pmin(x, x_max), x_min)
}

limit_range.RasterLayer <- function(x, x_min, x_max) {
  max(min(x, x_max), x_min)
}

limit_angle <- function(theta) {
  ((theta + pi) %% (2 * pi)) - pi
}

kelvin_to_celsius <- function(temp_kelvin) {
  temp_kelvin - 273.15
}

celsius_to_kelvin <- function(temp_celsius) {
  temp_celsius + 273.15
}

deg_to_rad <- function(deg) {
  deg * pi / 180
}

rad_to_deg <- function(rad) {
  rad * 180 / pi
}

to_deg_west <- function(lon) {
  ifelse(lon < 0, abs(lon), 360 - lon)
}

central_tz_lon <- function(lon) {
  round(ifelse(lon < 0, abs(lon), 360 - lon) / 15) * 15
}

lon_to_tz <- function(lon) {
  paste0("Etc/GMT", sprintf("%+i", -round(lon/15)))
}
