#' Calculates the equation of time using
#' @description
#' Calculates the seasonal correction of solar time.
#' @usage metric_equation_of_time(dt)
#' @param dt date and time (timezone aware).
#' @export
metric_equation_of_time <- function(dt) {
  dt <- lubridate::with_tz(dt, 'UTC')
  b <- 2 * pi * (lubridate::yday(dt)-81) / 364
  0.1645*sin(2*b) - 0.1255*cos(b) - 0.025*sin(b)
}

metric_solar_declination <- function(dt) {
  dt <- lubridate::with_tz(dt, 'UTC')
  0.409 * sin(2 * pi * lubridate::yday(dt) / 365 - 1.39)
}

metric_relative_distance_earth_sun <- function(dt) {
  1 + 0.033 * cos(2 * pi * lubridate::yday(lubridate::with_tz(dt, 'UTC')) / 365)
}

metric_hour_angle <- function(lon, dt) {
  dt <- with_tz(dt, 'UTC')
  dt <- as.numeric(difftime(dt, floor_date(dt, 'days'), units = 'hours'))
  sz <- metric_equation_of_time(dt)
  res <- (((dt + lon / 15 + sz) - 12) * pi / 12)
  ((res+pi) %% (2*pi)) - pi
}

metric_solar_incidence_angle <- function(terrain, datetime) {
  phi = terra::init(terrain, 'y') * pi / 180
  lon = terra::init(terrain, 'x')
  delta = metric_solar_declination(datetime)
  omega = metric_hour_angle(lon, datetime)
  slope = terra::terrain(terrain, v = 'slope', neighbors = 8, unit = 'radians')
  aspect = terra::terrain(terrain, v = 'aspect', neighbors = 8, unit = 'radians')

  sd = sin(delta);  cd = cos(delta)
  sp = sin(phi);    cp = cos(phi)
  so = sin(omega);  co = cos(omega)
  ss = sin(slope);  cs = cos(slope)
  sa = sin(aspect); ca = cos(aspect)

  res = sd*sp*cs + sd*cp*ss*ca + cd*cp*cs*co - cd*sp*ss*ca*co - cd*sa*ss*so
  acos(res)
}

# Calculates the precipitable water on the atmosphere.
# patm:
#  near-surface atmospheric pressure [kPa]
# ea:
#  near surface vapor pressure [kPa]
metric_precipitable_water <- function(patm, ea) {
  0.14 * ea * patm + 2.1
}

# Calculates the clearness index for direct beam radiation.
# patm:
#  atmospheric pressure [kPa]. May be a raster or a number.
# theta:
#  solar incidence angle over a horizontal surface.
# kt:
#  unitless turbidity coefficient. kt = 1 for clean air and kt = 0.5 for extremely polluted, turbid or dusty air. May be a raster or a number.
# pw:
#  precipitable water in the atmosphere [mm]. May be calculated using metric_precipitable_water. May be a raster or a number.
# --------------------
# Reference:
# Allen; Walter; Elliott; Howell; Itenfisu; Jensen; Snyder. The ASCE Standardized Reference Evapotranspiration Equation.
# Available on: https://doi.org/10.1061/9780784408056.
metric_kb <- function(patm, theta, pw, kt = 1){
  cos_theta = cos(theta)
  p1 = -0.00146 * patm / (kt * cos_theta)
  p2 = -0.075 * (pw / cos_theta) ^ 0.4
  0.98 * exp(p1 + p2)
}

# Calculates the clearness index for diffuse beam radiation.
# kb - number or raster
#
# --------------------
# Reference:
# Allen; Walter; Elliott; Howell; Itenfisu; Jensen; Snyder. The ASCE Standardized Reference Evapotranspiration Equation.
# Available on: https://doi.org/10.1061/9780784408056.
metric_kd <- function(kb){
  terra::ifel(kb >= 0.15, 0.35 - 0.36*kb, 0.18 + 0.82*kb)
}

metric_broadband_atmospheric_transmissivity <- function(patm, theta, pw, kt = 1) {
  kb = metric_kb(patm, theta, pw, kt)
  kd = metric_kd(kb)
  (kb + kd)
}
