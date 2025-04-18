
#' Calculates the radiant emittance of a body as a function of its temperature
#' @description Calculates the radiant emittance of a body as a function of its
#'  absolute temperature.
#' @param temp absolute temperature (K)
#' @param epsilon body's emissivity (adimensional between 0 and 1). A value of 1
#'  calculates black body radiant emittance.
#' @export
stefan_boltzmann_law <- function(temp, epsilon = 1.) {
  epsilon * 5.670374419184429453e-8 * temp ** 4
}

#' Calcula a declinação solar
#' @description Calcula a declinação solar. Declinação solar (δ) é o ângulo
#' formado entre o Sol e o plano do equador terrestre.
#' @usage solar_declination(j)
#' @param j dia do ano (1 a 365 ou 366), objeto Date ou POSIXt
#' @details Ao longo do ano, a
#' declinação solar varia de -23.44° a 23.44° (em radianos: -0.4091052 rad a
#' +0.4091052 rad), com valores negativos ocorrendo no verão do hemisfério sul e
#' valores positivos no verão do hemisfério norte. A declinação solar é igual a
#' 0 durante os equinócios.
#' @export
solar_declination <- function(j) {
  UseMethod("solar_declination", j)
}

#' @export
solar_declination.default <- function(j) {
  0.409105176667471 * sin((2 * pi * j / 365) - 1.39)
}

#' @export
solar_declination.Date <- function(j) {
  solar_declination.default(lubridate::yday(j))
}

#' @export
solar_declination.POSIXt <- function(j) {
  gmt_date_time <- lubridate::with_tz(j, tzone = "GMT")
  solar_declination.default(lubridate::yday(gmt_date_time))
}

#' Calcula a distância inversa relativa entre a Terra e Sol
#' @description A distância relativa inversa entre a Terra e o Sol é um número
#' variando de 0,967 UA no afélio a 1,033 UA no periélio.
#' @details A Terra orbita o Sol em uma trajetória aproximadamente elíptica. O momento em
#' que a Terra está mais próxima do Sol é chamado periélio, e o momento onde a distância
#' é máxima é chamado de afélio. Quanto mais próximo do limite superior (1,033),
#' maior a proximidade entre a Terra e o Sol.
#' @param j numeric representando dia do ano, entre 1 e 365, Date ou POSIXt.
#' @export
inv_rel_dist_sun <- function(j) {
  UseMethod("inv_rel_dist_sun", j)
}

#' @export
inv_rel_dist_sun.default <- function(j) {
  1 + 0.033 * cos(2 * pi * j / 365)
}

#' @export
inv_rel_dist_sun.Date <- function(date) {
  inv_rel_dist_sun.default(lubridate::yday(date))
}

#' @export
inv_rel_dist_sun.POSIXt <- function(date_time) {
  gmt_date <- lubridate::with_tz(date_time, "GMT")
  inv_rel_dist_sun.default(lubridate::yday(gmt_date))
}

#' Calculates the water vapor saturation pressure as a function of temperature.
#' @param temp temperature in Celsius (ºC). Rasters are also supported.
#' @return a vector of the same length as temp with the calculated
#' water vapor saturation pressure in kPa for each temperature.
#' @examples
#' > tetens_equation(seq(5, 30, by = 5))
#' [1] 0.872311 1.227963 1.705346 2.338281 3.167778 4.243065
#' @export
tetens_equation <- function(temp) {
  0.6108 * exp(17.27 * temp / (237.3 + temp))
}

#' @export
latent_heat <- function(temp) {
  l0 <- 2.501e6
  cpcv <- 2.180e3
  t0 <- 273.16
  l0 - cpcv * (temp - t0)
}

#' Calculates the water vapor saturation pressure as a function of temperature.
#' @param temp temperature in Kelvin (K).
#' @return a vector of the same length as temp with the calculated
#' water vapor saturation pressure in Pa for each temperature.
#' @export
ambaum_equation <- function(temp) {
  l0 <- 2.501e6
  cpcv <- 2.180e3
  t0 <- 273.16
  rv <- 461.52
  lambda <- latent_heat(temp)
  611.655 * exp(l0 / (rv * t0) - lambda / (rv * temp)) * (t0 / temp) ^ (cpcv / rv)
}

#' Calculates the slope of the equation of Tetens
#' @param  temp temperature in Kelvin (K).
#' @return a vector with the same length as temp with the calculated
#' slope in Pa/K.
#' @export
ambaum_slope <- function(temp) {
  lambda <- latent_heat(temp)
  rv <- 461.52
  es <- ambaum_equation(temp)
  lambda * es / (rv * temp * temp)
}

#' Calculates the psychrometric constant
#' @param patm atmospheric pressure (Pa).
#' @return Psychrometric constant (Pa / K)
#' @export
psy_const <- function(patm) {
  6.65e-4 * patm
}

#' Estimates wind speed at height z2 based on wind speed measued at
#' height z1, given the height of the zero plane displacement and the roughness
#' length of the surface.
#'
#' For a standardized FAO Peman Monteith ET0 surface
#' you can use the default values for d and z0. If the wind speed is measured
#' at a standard weather station, which measures wind speed
#' at a 10m height, you can use the default value for z1.
#' @param u1 Wind speed [m/s] measured at height z1.
#' @param z2 Height z2 [m].
#' @param z1 Height z1 [m]. The default is 10.
#' @param d Zero plane displacement height. If not set, a default value = 0.08 will be set, which
#' is the zero plane displacement height estimated for a 0.12m height uniform crop.
#' @param z0 Roughness length. If not set, a default value of 0.01476 will be set, which
#' corresponds to the roughness length of the standardized FAO ET0 Penman-Monteith
#' equation.
#' @return Wind speed at height z2.
#' @export
log_wind_profile <- function(u1, z2 = 2, z1 = 10, d = 0.084, z0 = 0.01476) {
  u1 * log((z2 - d) / z0) / log((z1 - d) / z0)
}

#' Calculates the slope of the equation of Tetens
#' @param  temp temperature in Celsius (ºC). Rasters are also supported.
#' @return a vector of the same length as temp with the calculated
#' slope in kPa/ºC.
#' @export
tetens_slope <- function(temp) {
  4098 * tetens_equation(temp) / ((237.3 + temp) ** 2)
}

#' @export
sunset_hour_angle <- function(j, lat) {
  delta <- solar_declination(j)
  phi <- deg_to_rad(lat)
  x <- limit_range(-tan(phi) * tan(delta), -1, 1)
  acos(x)
}

#' @export
equation_of_time <- function(j) {
  b <- 2 * pi * (j - 81) / 364
  0.1645 * sin(2 * b) - 0.1255 * cos(b) - 0.025 * sin(b)
}

# Calculates the hour angle
#' @export
hour_angle <- function(dt, lon) {
  UseMethod("hour_angle", dt)
}

#' @export
hour_angle.default <- function(dt, lon, j) {
  sc <- equation_of_time(j)
  ret <- (dt + lon / 15 + sc - 12) * pi / 12
  return(ret)
}

#' @export
hour_angle.POSIXt <- function(dt, lon) {
  dt <- lubridate::with_tz(dt, tzone = "GMT")
  t <- as.numeric(lubridate::as.duration(dt - lubridate::floor_date(dt, unit = "days"))) / 3600
  ret <- hour_angle.default(t, lon, lubridate::yday(dt))
  return(ret)
}

#' Calculates atmospheric pressure at a given altitude above sea level.
#' @param z Altitude above sea level [m].
#' @param temp Meam atmospheric temperature [K]. The default is 293.15.
#' @param lb Temperature lapse rate [K / m] (i.e. how many Kelvin the
#'  temperature of air decreases with a 1 m increase in altitude). The default
#'  is 5e-3 K / m.
#' @return Atmospheric pressure at altitude z.
#' @export
atmospheric_pressure <- function(z, temp = 293.15, lb = 6.5e-3) {
  p0 <- 101325 # standard atmospheric pressure at sea level.
  g <- 9.80665 # standard gravity on Earth's surface.
  rd <- 287.058 # specific gas constant for dry air.
  power <- -g / (rd * lb)
  return(p0 * ((temp + lb * z) / temp) ** power)
}

#' ext_rad
#' @description calculates the exo-atmospheric radiation, in MJ/m^2, for time intervals
#'  less than a whole day. It can be used to calculate extraterrestrial radiations
#'  for a whole day, but computations are a lot more expensive than using ext_rad_day.
#' @aliases
#'  ext_rad.default(t0, t1, lat, lon, j)
#'  ext_rad.POSIXct(t0, t1, lat, lon)
#' @param t0,t1 time at the beginning and end of the period.
#' @param lat latitude in degres from -180 to 180 (negative values are west of Greenwich)
#' @param lon longitude in degres from -90 to 90 (negative values are south of Equator)
#' @param ... other variables
#'  @param j day of the year (1 to 366) (only valid if t0 and t1 are not Date or POSIXct objects)
#' @examples
#' t0 <- 0:23
#' t1 <- (t0 + 1) %% 24
#' ext_rad(t0, t1, -30, -48, j = 90)
#' @references ALLEN, R. G. et al. Crop evapotranspiration - Guidelines for
#'  computing crop water requirements - FAO Irrigation and drainage paper 56. FAO: Rome, 1998.
#' @export
ext_rad <- function(t0, t1, lat, lon, ...) {
  UseMethod("ext_rad", t0)
}

#' @export
ext_rad.default <- function(t0, t1, lat, lon, j) {
  # it can be better!
  # calculates hour angle range (omega). Omega usually varies from
  # omega - omega_dif / 2 to omega + omega_dif / 2, except during sunrise and
  # sunset.
  omega_dif <- ifelse(t0 > t1, 24 + t1 - t0, t1 - t0) * pi / 12
  # calculates the the time between t0 and t1
  tm <- (ifelse(t0 > t1, (24 + t1 + t0), (t0 + t1)) / 2) %% 24

  omega <- hour_angle.default(tm, lon, j = j)
  delta <- solar_declination(j)
  phi <- deg_to_rad(lat)
  dr <- inv_rel_dist_sun(j)
  omega_s <- sunset_hour_angle(j, lat)

  # fd is TRUE if
  fd <- -tan(phi) * tan(delta) < -1

  omega_1 <- omega - omega_dif / 2
  omega_1 <- ifelse(fd, omega_1, limit_range(omega_1, x_min = -omega_s, x_max = omega_s))
  omega_2 <- omega + omega_dif / 2
  omega_2 <- ifelse(fd, omega_2, limit_range(omega_2, x_min = -omega_s, x_max = omega_s))

  omega_dif <- ifelse(omega_2 < omega_1, 2 * pi + omega_2 - omega_1, omega_2 - omega_1)

  (omega_dif * sin(phi) * sin(delta) + cos(phi) * cos(delta) * (sin(omega_2) - sin(omega_1))) * dr * 59.04/pi
}

#' @export
ext_rad.POSIXt <- function(t0, t1, lat, lon, j = NULL) {
  t0 <- do.call(c, Map(function(dt) lubridate::with_tz(dt, 'GMT'), t0))
  t1 <- do.call(c, Map(function(dt) lubridate::with_tz(dt, 'GMT'), t1))
  tm <- do.call(c, Map(function(x, y) mean(c(x, y)), t0, t1))

  ext_rad.default(
    (as.numeric(t0) - as.numeric(lubridate::floor_date(t0, unit = "days"))) / 3600,
    (as.numeric(t1)  - as.numeric(lubridate::floor_date(t0, unit = "days"))) / 3600,
    lat, lon,
    lubridate::yday(tm)
  )
}

#' Calcula a radiação solar exoatmosférica
#' @description Calcula a radiação solar extraterrestre para o dia e a latitude informadas
#' @param j dia do ano (de 1 a 366) ou Date ou POSIXct.
#' @param lat latitude, em graus decimais (valores negativos para o hemisfério sul)
#' @export
ext_rad_day <- function(j, lat) {
  dr <- inv_rel_dist_sun(j)
  delta <- solar_declination(j)
  phi <- deg_to_rad(lat)
  omega_s <- sunset_hour_angle(j, lat)
  gsc <- 1.360e3
  secs <- 8.64e4
  secs * gsc * dr * (omega_s * sin(phi) * sin(delta) + cos(phi) * cos(delta) * sin(omega_s)) / pi
}

#' Calcula a duração do dia
#' @description Calcula a duração do dia para a data informada e a latitude
#' @param j dia do ano (inteiro de 1 a 366), Date ou POSIXct
#' @param lat latitude, em graus decimais (valores negativos para hemisfério sul)
#' @return duração do dia, em horas decimais.
#' @export
day_length <- function(j, lat) {
  sunset_hour_angle(j, lat) * 24 / pi
}

#' @export
solar_elevation_angle <- function(dt, lat, lon){
  ha <- hour_angle(dt, lon)
  phi <- lat * pi / 180
  declination <- solar_declination(dt)
  return(asin((cos(phi) * cos(declination) * cos(ha)) + (sin(phi) * sin(declination))))
}


#' @export
solar_zenith_angle <- function(dt, lat, lon){
  ha <- hour_angle(dt, lon)
  phi <- lat * pi / 180
  declination <- solar_declination(dt)
  return(acos((cos(phi) * cos(declination) * cos(ha)) + (sin(phi) * sin(declination))))
}

#' Calcula o saldo de radiação de ondas longas
#' @description Calcula a radiação de ondas longas conforme metodologia
#' proposta por Allen et al. (1998)
#' @param tmax Temperatura máxima do ar medida dentro do intervalo de tempo do cálculo [ºC]
#' @param tmin Temperatura mínima do ar medida dentro do intervalo de tempo do cálculo [ºC]
#' @param e Pressão parcial do vapor [kPa] média dentro do intervalo de tempo do cálculo.
#' @param rs Radiação solar de onda curta incidente sobre uma superfície horizontal, durante o intervalo de tempo do cálculo [MJ/m²]
#' @param rs0 Radiação solar de onda curta em condições de céu claro [MJ/m²]
#' @param hourly TRUE para intervalo de tempo de 1 hora, FALSE caso contrário.
#' @references
#' ALLEN, R. G.; PEREIRA, L. S.; RAES, D.; SMITH, M. FAO Irrigation and Drainage Paper no 56. FAO: Rome, 1998.
#' @export
net_longwave <- function(tmax, tmin, e, rs, rs0) {
  tmaxk <- celsius_to_kelvin(tmax)
  tmink <- celsius_to_kelvin(tmin)
  sigma <- 4.903e-9
  radRatio <- ifelse(rs0 == 0, 0, pmin(rs / rs0, 1))
  ret <- -sigma * (tmaxk ^ 4 + tmink ^ 4) / 2
  ret <- ret * (0.34 - 0.14 * sqrt(e))
  ret * ((1.35 * radRatio) - 0.35)
}

#' Calculates net shortwave radiation
#' @description Calculates net shortwave radiation, which is the difference between
#' incoming shortwave radiation and reflected shortwave radiation over a surface.
#' @param rs global radiation (W/m², MJ/m²/d, kJ/m², etc...)
#' @param albedo surface albedo (adimensional) between 0 and 1
#' @details O parâmetro albedo possui valor padrão de 0.23. Este é o albedo para a cultura
#'  de referência, utilizado no cálculo da evapotranspiração de referência. Para maiores detalhes
#'  veja a seção 'details' da função \code{\link{fao_et0}}.
#' @export
net_shortwave <- function(rs, albedo = 0.23) {
  (1 - albedo) * rs
}

#' Calcula a evapotranspiração de referência
#' @description Calcula a evapotranspiração de referência (ET0, mm) pelo método de Penman-Monteith
#'  padronizado pela FAO (Allen et al., 1998) com ou sem as alterações propostas pela American Society of Civil Engineers (ASCE, 2005)
#' @param temp temperatura média do ar (ºC) (ver detalhes abaixo)
#' @param rn saldo de radiação [MJ/m²]
#' @param g fluxo de calor no solo [MJ/m²] (ver detalhes)
#' @param gamma coeficiente psicrométrico (ver detalhes)
#' @param es pressão de saturação do vapor de água (kPa) (ver detalhes)
#' @param e pressão parcial do vapor de água (kPa) (ver detalhes)
#' @param hourly FALSE para calcular a evapotranspiração diária (mm/dia) (valor padrão), ou TRUE para
#'  calcular a evapotranspiração horária (mm/h)
#' @param mode um dos seguintes valores possíveis 'FAO', 'ASCE-ET0', 'ASCE-ETR' (ver detalhes)
#' @param day válido somente no cálculo da evapotranspiração horária. Um vetor booleano informando para
#'  cada registro se é dia ou noite (TRUE para dia e FALSE para noite).
#' @details Para fins de padronização, no cálculo para intervalos de tempo de 1 dia ou maiores, deve-se definir \code{temp}
#'  como a média entre a temperatura máxima
#'  e a temperatura mínima do ar medidas no período de tempo em que se deseja calcular a evapotranspiração,
#'  e não como a média de todas as medições de temperatura no período (Allen et al., 1998).\cr
#'  Caso não se disponha de medições diretas do parâmetro \code{rn}, pode-se estimar o seu
#'  valor por meio das funções \code{\link{net_shortwave}} e \code{\link{net_longwave}}.\cr
#'  Para o cálculo da evapotranspiração
#'  de referência, o albedo da cultura de referência deve ser fixado em 0.23 tanto para a grama quanto para a alfafa (Allen et al., 2005).
#'  Caso não se disponha de medições diretas do parâmetro \code{g}, Allen et al. (1998) recomendam
#'  atribuir o valor de 0 para cálculos de evapotranspiração diária ou \code{0.1 * rn} durante o dia
#'  e \code{0.5 * rn} durante a noite, para cálculos de evapotranspiração horária tendo como referência a grama.
#'  Caso a referência adotada seja a alfafa e não haja medições diretas do fluxo de calor no solo, deve-se adotar
#'  \code{0.04 * rn} durante o dia e \code{0.2 * rn} durante a noite (Allen et al., 2005)\cr
#'  O valor de \code{es} pode ser calculado por meio da função \code{\link{tetens_equation}}. Recomenda-se
#'  definir es como a média entre os valores retornados por \code{\link{tetens_equation}} para
#'  a temperatura mínima e para a temperatura máxima no período, para cálculos realizados
#'  no intervalo de tempo de 1 dia ou maior.\cr
#'  O parâmetro \code{mode} pode ser definido como 'FAO', 'ASCE-ET0' ou 'ASCE-ETR'. O valor padrão é 'ASCE-ET0'. Se definido
#'  como 'FAO', o cálculo é realizado com os parâmetros definidos em Allen et al. (1998). Se definido como 'ASCE-ET0',
#'  o cálculo considera os parâmetros definidos para a grama de referência conforme Allen et al. (2006). Se definido
#'  como 'ASCE-ETR', os parâmetros utilizados correspondem ao da alfafa como cultura de referência.
#' @references
#' ALLEN, R. G.; PEREIRA, L. S.; RAES, D.; SMITH, M. FAO Irrigation and Drainage Paper no 56. FAO: Rome, 1998.\cr
#' ALLEN, et al.. A recommendation on standardized surface resistance for hourly calculation of reference ETo
#'  by the FAO56 Penman-Monteith method. Agricultural Water Management, 81, p. 1-22, 2006.
#'  Available from: <https://www.sciencedirect.com/science/article/abs/pii/S037837740500154X>
#' @export
fao_et0 <- function(temp, rn, g, gamma, u2, es, e, hourly = FALSE, mode = 'ASCE-ET0', day = NA) {

  delta <- tetens_slope(temp)

  if (hourly) {
    # for hourly calculations
    cn <- switch(mode,
                 'FAO' = 37,
                 'ASCE-ET0' = 37,
                 'ASCE-ETR' = 66)
    cd <- switch(mode,
                 'FAO' = 0.34,
                 'ASCE-ET0' = ifelse(day, 0.24, 0.96),
                 'ASCE-ETR' = ifelse(day, 0.25, 1.7))
  } else {
    # for daily calculations
    cn <- switch(mode,
                 'FAO' = 900,
                 'ASCE-ET0' = 900,
                 'ASCE-ETR' = 1600)
    cd <- switch(mode,
                 'FAO' = 0.34,
                 'ASCE-ET0' = 0.34,
                 'ASCE-ETR' = 0.38)
  }

  ((0.408*delta*(rn-g)) + (gamma*cn*u2*(es-e) / (temp+273))) / (delta+gamma * (1 + cd*u2))
}
