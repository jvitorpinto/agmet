# default name for variables
# - dt: datetime
# - temp: temperature
# - red: red band reflectance
# - nir: near-infrared band reflectance
# - epsilon: emissivity
# - theta: solar incidence angle
# - lat: latitude in degrees
# - lon: longitude in degrees

#' @description
#' Returns a pan-sharpened RGB image using red, green, blue and panchromatic reflectances
#' @param r red band reflectance
#' @param g green band reflectance
#' @param b blue band reflectance
#' @param p panchromatic band reflectance
#' @param min minimum reflectance value (i.e. it will correspond to R, G, or B 0x00)
#' @param max maximum reflectance value (i.e. it wiil correspond to R, G, or B 0xFF)
#' @param method resampling method passed to the `function terra::resample`
metric_pan_sharpening <- function(r, g, b, p, min = 0, max = 1, method = 'near') {
  rgb <- c(r, g, b)
  rgb <- 255 * terra::clamp((rgb - min) / (max - min), 0, 1)
  terra::RGB(rgb) <- 1:3
  val <- terra::clamp((p - min) / (max - min), 0, 1)
  hs <- terra::colorize(rgb, 'hsv')[[c('hue', 'saturation')]]
  hs <- terra::resample(hs, p, method = 'near')
  hsv <- c(hs, val)
  hsv
  terra::RGB(hsv, type = 'hsv') <- 1:3
  terra::colorize(hsv, 'rgb')
}

earth_sun_distance <- function(dt) {
  eccentricity <- 0.0167
  orbital_period <- 365.25
}

#' Calculates the square of the relative distance Earth-Sun.
metric_relative_distance_earth_sun <- function(dt) {
  j <- lubridate::yday(lubridate::with_tz(dt, 'UTC'))
  1 / (1 + 0.033*cos(2*pi*j/365))
}

#' Calculates the hour angle at a given time based on the longitude.
#' @param lon SpatRaster containing the longitude of each raster cell.
#' @export
metric_hour_angle <- function(lon, dt) {
  dt <- lubridate::with_tz(dt, 'UTC')
  dt <- as.numeric(difftime(dt, lubridate::floor_date(dt, 'days'), units = 'hours'))
  sz <- equation_of_time(dt)
  res <- (((dt + lon / 15 + sz) - 12) * pi / 12)
  ((res+pi) %% (2*pi)) - pi
}

#' Calculates the solar incidence angle over a horizontal surface.
#' @export
metric_horizontal_incidence_angle <- function(lat, lon, datetime) {
  phi = lat * pi / 180
  delta = solar_declination(datetime)
  omega = metric_hour_angle(lon, datetime)

  sd = sin(delta);  cd = cos(delta)
  sp = sin(phi);    cp = cos(phi)
  so = sin(omega);  co = cos(omega)

  res = sd*sp + cd*cp*co
  acos(res)
}

#' @description
#' Calculates the solar incidence angle over a sloped surface.
#'
#' @param phi latitude in radians as SpatRaster.
#' @param slope terrain slope in radians as SpatRaster.
#' @param aspect terrain aspect in radians as SpatRaster.
#' @param delta solar declination as numeric.
#' @param omega hour angle as numeric.
#'
#' @details
#' The parameter `phi` can be obtained by passing a raster with the CRS correctly
#' set to the function `metric_get_coordinates`, and then multiplying the result
#' by 180/pi.
#'
#' The parameters `aspect` and `slope` can be obtained with the function terra::terrain.
#'
#' The parameter `delta` can be obtained with the function `solar_declination`.
#'
#' The parameter `omega` can be obtained with the function `metric_hour_angle`,
#'
metric_solar_incidence_angle <- function(phi, slope, aspect, delta, omega) {
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

#' @description
#' Calculates the broadband atmospheric transmissivity.
#' @param theta horizontal incidence angle calculated by `metric_horizontal_incidence_angle`
#' @param ea vapor pressure [kPa]
#' @param patm atmospheric pressure [kPa], can be calculated from the altitude using `atmospheric_pressure`
#' @param kt atmospheric turbidity.
#' @export
metric_broadband_atmospheric_transmissivity <- function(patm, theta, ea, kt = 1) {
  kb = metric_kb(patm, theta, ea, kt)
  kd = metric_kd(kb)
  (kb + kd)
}


metric_exoatmospheric_irradiance <- function(terrain, pw, kt, dt) {
  squared_earth_sun_distance <- metric_relative_distance_earth_sun(dt)
  cos_theta <- cos(metric_solar_incidence_angle(terrain, dt))
  theta_hor <- metric_horizontal_incidence_angle(init(terrain, 'y'), init(terrain, 'x'), dt)
  metric_broadband_atmospheric_transmissivity(patm, theta_hor, pw, kt)
}

#' Metric NDVI
#' @description
#' Calculates the Normalized Difference Vegetation Index.
#' @param red reflectance on the red band.
#' @param nir reflectance on the near-infrared band.
#' @details
#' The normalizeed difference vegetation Index is calculated as
#' NDVI = (nir - red) / (nir + red).
#' Which results in a number between -1 and 1.
metric_ndvi <- function(red, nir) {
  (nir - red) / (nir + red)
}

#' Calculates the albedo from a landsat 8 image using the procedure described by
#' Silva et al. (2016).
#' @description
#' Calculates the broadband albedo through the integration of narrow band reflectance
#' for bands 2, 3, 4, 5, 6, and 7 of Landsat 8.
#'
#' @param b2 Landsat 8 top-of-atmosphere reflectance at band 2
#' @param b3 Landsat 8 top-of-atmosphere reflectance at band 3
#' @param b4 Landsat 8 top-of-atmosphere reflectance at band 4
#' @param b5 Landsat 8 top-of-atmosphere reflectance at band 5
#' @param b6 Landsat 8 top-of-atmosphere reflectance at band 6
#' @param b7 Landsat 8 top-of-atmosphere reflectance at band 7
#' @param alpha_atm atmospheric albedo (see details)
#' @param tau atmospheric transmissivity
#'
#' @details
#' The atmospheric albedo can be obtained from a radiative transfer model, such as
#' MODTRAN. It is usually a value between 0.025 and 0.040. Here we use a 'default'
#' value of 0.03, but we recommend you to adapt the value according to your needs.
#'
metric_albedo_l8 <- function(b2, b3, b4, b5, b6, b7, tau, alpha_atm = 0.03) {
  w2 <- 0.300
  w3 <- 0.277
  w4 <- 0.233
  w5 <- 0.143
  w6 <- 0.036
  w7 <- 0.012

  alpha_toa <- w2*b2 + w3*b3 + w4*b4 + w5*b5 + w6*b6 + w7*b7 / (w2 + w3 + w4 + w5 + w6 + w7)

  terra::clamp((alpha_toa - alpha_atm) / (tau * tau), 0, 1)
}

metric_albedo_l9 <- metric_albedo_l8

#' Retrieves the longitude and latitude of each raster cell in a new raster, regardless
#' of the coordinate reference system of the source image.
#' @param src source image
metric_get_coordinates <- function(src) {
  coords <- c(terra::init(src, 'x'), terra::init(src, 'y'))
  names(coords) <- c('x', 'y')
  terra::values(coords) <- sf::sf_project(from = sf::st_crs(src), to = sf::st_crs(4326), coords)
  coords
}

#' @export
metric_terrain <- function(terrain) {
  slope <- terra::terrain(terrain, v = 'slope', neighbors = 8, unit = 'radians')
  aspect <- terra::terrain(terrain, v = 'aspect', neighbors = 8, unit = 'radians')
  coords <- metric_get_coordinates(terrain)
  coords_rad <- coords * pi/180
  res <- c(terrain, slope, aspect, coords, coords_rad)
  names(res) <- c('ALTITUDE', 'SLOPE', 'ASPECT', 'LON_DEG', 'LAT_DEG', 'LON_RAD', 'LAT_RAD')
  res
}

#' Calculates the soil adjusted vegetation index (SAVI).
#' @param red top-of-atmosphere reflectance on the red band.
#' @param nir top-of-atmosphere reflectance on the near-infrared band.
#' @param const numeric constant usually set to 0.1 for METRIC applications.
metric_savi <- function(red, nir, const = 0.1) {
  (1 + const) * (nir - red) / (const + nir + red)
}

#' Calculates the leaf area index based on SAVI.
#' @param savi soil adjusted vegetation index.
#' @export
metric_leaf_area_index <- function(x) {
  UseMethod('metric_leaf_area_index', x)
}

#' @export
metric_leaf_area_index.default <- function(savi) {
  savi <- limit_range(savi, 0.1, 0.69)
  limit_range(-log((0.69 - savi) / 0.59) / 0.91, 0, 6)
}

#' @export
metric_leaf_area_index.SpatRaster <- function(savi) {
  savi <- terra::clamp(savi, 0.1, 0.69)
  terra::clamp(-log((0.69 - savi) / 0.59) / 0.91, 0, 6)
}

#' @export
metric_leaf_area_index.Metric <- function(x) {
  metric_leaf_area_index.SpatRaster(x$rast$SAVI)
}

metric_broadband_surface_emissivity <- function(lai) {
  min(0.95 + 0.01 * lai, 0.98)
}

metric_effective_atmospheric_emissivity <- function(tau) {
  0.85 * (-log(tau)) ^ 0.09
}

#===============================================================================
# Landsat 8 specific functions
#===============================================================================

extract_param <- function(data, param, type = as.character) {
  filtered_data <- data[stringr::str_detect(names(data), paste0('^', param, '_'))]
  names <- stringr::str_remove(names(filtered_data), paste0('^', param, '_'))
  res <- data.frame(
    name = names,
    value = type(filtered_data)
  )
  names(res)[2] <- tolower(param)
  res
}

landsat8_level1_radiometric_rescaling <- function(x) {
  res <- list(
    extract_param(x, 'RADIANCE_MULT', as.numeric),
    extract_param(x, 'RADIANCE_ADD', as.numeric),
    extract_param(x, 'REFLECTANCE_MULT', as.numeric),
    extract_param(x, 'REFLECTANCE_ADD', as.numeric)
  )
  purrr::reduce(res, merge, all = T)
}

landsat8_level1_thermal_constants <- function(x) {
  res <- list(
    extract_param(x, 'K1_CONSTANT', as.numeric),
    extract_param(x, 'K2_CONSTANT', as.numeric)
  )
  purrr::reduce(res, merge, all = T)
}

landsat8_load_projection_attributes <- function(x) {
  proj_crs <- sf::st_crs(32600 + as.integer(x$UTM_ZONE))
  proj_pts <- c(
    x$CORNER_LL_PROJECTION_X_PRODUCT, x$CORNER_LL_PROJECTION_Y_PRODUCT,
    x$CORNER_UL_PROJECTION_X_PRODUCT, x$CORNER_UL_PROJECTION_Y_PRODUCT,
    x$CORNER_UR_PROJECTION_X_PRODUCT, x$CORNER_UR_PROJECTION_Y_PRODUCT,
    x$CORNER_LR_PROJECTION_X_PRODUCT, x$CORNER_LR_PROJECTION_Y_PRODUCT,
    x$CORNER_LL_PROJECTION_X_PRODUCT, x$CORNER_LL_PROJECTION_Y_PRODUCT
  )
  proj_pts <- matrix(as.numeric(proj_pts), ncol = 2, byrow = T)
  proj_pts <- sf::st_sfc(sf::st_polygon(list(proj_pts)), crs = proj_crs)

  # do the same for latitude/longitude coordinates.
  coor_crs <- sf::st_crs(4326)
  coor_pts <- c(
    x$CORNER_LL_LON_PRODUCT, x$CORNER_LL_LAT_PRODUCT,
    x$CORNER_UL_LON_PRODUCT, x$CORNER_UL_LAT_PRODUCT,
    x$CORNER_UR_LON_PRODUCT, x$CORNER_UR_LAT_PRODUCT,
    x$CORNER_LR_LON_PRODUCT, x$CORNER_LR_LAT_PRODUCT,
    x$CORNER_LL_LON_PRODUCT, x$CORNER_LL_LAT_PRODUCT
  )
  coor_pts <- matrix(as.numeric(coor_pts), ncol = 2, byrow = T)
  coor_pts <- sf::st_sfc(sf::st_polygon(list(coor_pts)), crs = coor_crs)

  list(
    projection = proj_pts,
    lonlat = coor_pts
  )
}

#' Reads landsat metadata from a json file and converts it to a more accessible
#' format that R can understand.
#' @description
#' A short description...
#' @export
landsat8_read_meta <- function(filename) {
  meta <- jsonlite::read_json(filename)

  res <- list()
  res$filename <- filename
  res$processing_level <- meta$LANDSAT_METADATA_FILE$PRODUCT_CONTENTS$PROCESSING_LEVEL

  filenames = extract_param(meta$LANDSAT_METADATA_FILE$PRODUCT_CONTENTS, 'FILE_NAME')
  datatypes =  extract_param(meta$LANDSAT_METADATA_FILE$PRODUCT_CONTENTS, 'DATA_TYPE')
  res$files <- merge(filenames, datatypes, all = T)

  res$datetime = lubridate::as_datetime(paste(meta$LANDSAT_METADATA_FILE$IMAGE_ATTRIBUTES$DATE_ACQUIRED, meta$LANDSAT_METADATA_FILE$IMAGE_ATTRIBUTES$SCENE_CENTER_TIME, sep = 'T'))

  if (res$processing_level == 'L2SP') {
    level2_surface_reflectance_parameters <- meta$LANDSAT_METADATA_FILE$LEVEL2_SURFACE_REFLECTANCE_PARAMETERS
    level2_surface_reflectance_parameters <- list(
      extract_param(level2_surface_reflectance_parameters, 'REFLECTANCE_MAXIMUM', as.numeric),
      extract_param(level2_surface_reflectance_parameters, 'REFLECTANCE_MINIMUM', as.numeric),
      extract_param(level2_surface_reflectance_parameters, 'QUANTIZE_CAL_MAX', as.numeric),
      extract_param(level2_surface_reflectance_parameters, 'QUANTIZE_CAL_MIN', as.numeric),
      extract_param(level2_surface_reflectance_parameters, 'REFLECTANCE_MULT', as.numeric),
      extract_param(level2_surface_reflectance_parameters, 'REFLECTANCE_ADD', as.numeric)
    )
    res$surface_reflectance_params <- purrr::reduce(level2_surface_reflectance_parameters, merge, all = T)

    surface_temperature_parameters <- meta$LANDSAT_METADATA_FILE$LEVEL2_SURFACE_TEMPERATURE_PARAMETERS
    surface_temperature_parameters <- list(
      extract_param(surface_temperature_parameters, 'TEMPERATURE_MAXIMUM', as.numeric),
      extract_param(surface_temperature_parameters, 'TEMPERATURE_MINIMUM', as.numeric),
      extract_param(surface_temperature_parameters, 'QUANTIZE_CAL_MAXIMUM', as.numeric),
      extract_param(surface_temperature_parameters, 'QUANTIZE_CAL_MINIMUM', as.numeric),
      extract_param(surface_temperature_parameters, 'TEMPERATURE_MULT', as.numeric),
      extract_param(surface_temperature_parameters, 'TEMPERATURE_ADD', as.numeric)
    )
    res$surface_temperature_params <- purrr::reduce(surface_temperature_parameters, merge, all = T)
  }

  if (res$processing_level == 'L1TP') {
    res$level1_radiometric_rescaling = landsat8_level1_radiometric_rescaling(meta$LANDSAT_METADATA_FILE$LEVEL1_RADIOMETRIC_RESCALING)
    res$level1_thermal_constants = landsat8_level1_thermal_constants(meta$LANDSAT_METADATA_FILE$LEVEL1_THERMAL_CONSTANTS)
  }

  res$coords <- landsat8_load_projection_attributes(meta$LANDSAT_METADATA_FILE$PROJECTION_ATTRIBUTES)

  res
}

landsat8_read_reflectance_file <- function(dir, filename, add, mult, aoi) {
  r <- terra::rast(paste(dir, filename, sep = '/'))
  r <- terra::crop(r, sf::st_transform(aoi, sf::st_crs(r))) * mult + add
  r
}

#' @export
landsat8_level2 <- function(meta, aoi) {
  dir <- dirname(meta$filename)

  ref_files <- merge(meta$files, meta$surface_reflectance_params)
  ref <- sapply(split(ref_files, 1:nrow(ref_files)), function(x){
    landsat8_read_reflectance_file(dir, x$file_name, x$reflectance_add, x$reflectance_mult, aoi)
  })
  ref <- terra::rast(ref)
  names(ref) <- ref_files$name

  temp_files <- merge(meta$files, meta$surface_temperature_params)
  temp <- sapply(split(temp_files, 1:nrow(temp_files)), function(x){
    r <- terra::rast(paste(dir, x$file_name, sep = '/'))
    r <- terra::crop(r, sf::st_transform(aoi, sf::st_crs(r))) * x$temperature_mult + x$temperature_add
    r
  })
  temp <- terra::rast(temp)
  names(temp) <- temp_files$name

  list(
    meta = meta,
    datetime = meta$datetime,
    rast = c(ref, temp)
  )
}

#' @export
landsat8_level1 <- function(meta, aoi) {
  dir <- dirname(meta$filename)
  ref_files <- merge(meta$files, meta$level1_radiometric_rescaling)[,c('name', 'file_name', 'reflectance_add', 'reflectance_mult')]
  ref_files <- na.omit(ref_files)
  ref <- sapply(split(ref_files, 1:nrow(ref_files)), function(x){
    r <- terra::rast(paste(dir, x$file_name, sep = '/'))
    r <- terra::crop(r, sf::st_transform(aoi, sf::st_crs(r))) * x$reflectance_mult + x$reflectance_add
    r
  })
  names(ref) <- ref_files$name
  pan = terra::rast(ref[ref_files$name == 'BAND_8'])
  names(pan) <- 'BAND_8'
  list(
    meta = meta,
    ref = terra::rast(ref[ref_files$name != 'BAND_8']),
    pan = pan
  )
}

#' @description
#' Loads Landsat 8 raster data and delivers it with the input format
#' of METRIC.
#' @param x metadata filename.
#' @export
metric_load_landsat_8 <- function(toa, sr, tau) {
  rgb <- terra::as.int(terra::round(c(
    255 * terra::clamp(sr$rast$BAND_4 / 0.15, 0, 1),
    255 * terra::clamp(sr$rast$BAND_3 / 0.15, 0, 1),
    255 * terra::clamp(sr$rast$BAND_2 / 0.15, 0, 1)
  )))

  ndvi <- metric_ndvi(sr$rast$BAND_4, sr$rast$BAND_5)
  savi <- metric_savi(toa$ref$BAND_4, toa$ref$BAND_5, 0.1)

  # weather data is needed to calculate the broadband albedo.
  # alpha <- metric_albedo_l8(
  #   toa$ref$BAND_2,
  #   toa$ref$BAND_3,
  #   toa$ref$BAND_4,
  #   toa$REF$BAND_5,
  #   toa$REF$BAND_6,
  #   toa$REF$BAND_7,
  #   tau
  # )
  alpha <- metric_albedo_l8(
    toa$ref$BAND_2,
    toa$ref$BAND_3,
    toa$ref$BAND_4,
    toa$ref$BAND_5,
    toa$ref$BAND_6,
    toa$ref$BAND_7,
    tau
  )

  metric_input_raster(rgb, ndvi, savi, alpha, sr$rast$BAND_ST_B10)
}

#' @export
metric_input_raster <- function(rgb, ndvi, savi, alpha, temp) {
  rast <- c(rgb, ndvi, savi, alpha, temp)
  names(rast) <- c('R', 'G', 'B', 'NDVI', 'SAVI', 'ALBEDO', 'TEMP')
  rast
}

#' Creates a object that represents a model with its data.
#' @param rast raster data. Must contain, at least, NDVI, SAVI, albedo and surface temperature.
#' @export
metric <- function(rast, terrain, dt, weather_station = NULL) {
  model <- list(
    rast = c(terrain, rast),
    datetime = dt
  )
  class(model) <- c('Metric')
  model
}


# testar se função funciona com raster
fun_psy200 <- function(l) {
  x <- ((1 - 16*200/l)^0.25)
  ifelse(l==0,0,ifelse(l>0,-5*2/l,2*log((1+x)/2)+log((1+x*x)/2)-2*atan(x)+0.5*pi))
}

# mesma coisa
fun_psy <- function(z, l) {
  x<-((1-16*z/l)^0.25)
  ifelse(l==0,0,ifelse(l>0,-5*z/l,2*log((1+x*x)/2)))
}

# mesma coisa
fun_psy200_r <- function(l) {
  x <- ((1 - 16*200/l)^0.25)
  terra::ifel(l==0,0,terra::ifel(l>0,-5*2/l,2*log((1+x)/2)+log((1+x*x)/2)-2*atan(x)+0.5*pi))
}

# mesma coisa
fun_psy_r <- function(z, l) {
  x<-((1-16*z/l)^0.25)
  terra::ifel(l==0,0,terra::ifel(l>0,-5*z/l,2*log((1+x*x)/2)))
}


#' @export
metric_run <- function(model, anchor_points, weather_data) {
  # Calculates leaf area index
  lai <- metric_leaf_area_index(model$rast$SAVI)
  names(lai) <- 'LEAF_AREA_INDEX'

  epsilon_0 <- metric_broadband_surface_emissivity(lai)

  lw_out <- stefan_boltzmann_law(model$rast$TEMP, epsilon_0)

  epsilon_a <- metric_effective_atmospheric_emissivity(tau)

  lw_in <- stefan_boltzmann_law(model$rast$TEMP, epsilon_a)

  theta_rel <- metric_solar_incidence_angle(
    model$rast$LAT_RAD,
    model$rast$SLOPE,
    model$rast$ASPECT,
    solar_declination(model$datetime),
    metric_hour_angle(model$rast$LON_DEG, model$datetime)
  )

  theta_hor <- metric_horizontal_incidence_angle(
    model$rast$LAT_RAD,
    model$rast$LON_DEG,
    model$datetime
  )

  dr <- metric_relative_distance_earth_sun(model$datetime)

  sw_in <- 1367 * cos(theta_rel) * tau / dr

  sw_out <- model$rast$ALBEDO * sw_in

  rn <- (sw_in - sw_out) + (lw_in - lw_out)
  names(rn) <- c('NET_RADIATION')

  g <- rn * (model$rast$TEMP - 273.15) * (0.0038 + 0.0074 * model$rast$ALBEDO) * (1 - 0.98 * model$rast$NDVI ^ 4)
  names(g) <- c('SOIL_HEAT_FLUX')

  # ap (anchor_points) é um objeto SpatVector contendo todos os dados do raster
  # do modelo para o local dos pontos de ancoragem
  ap <- terra::extract(c(model$rast, lai, rn, g), anchor_points, bind = T)
  ap$MOMENTUM_ROUGHNESS_LENGTH <- 0.018 * ap$LEAF_AREA_INDEX # mudar aqui para calcular um raster

  # momentum roughness length (eq 33, 34a or 34b)
  w <- weather_data
  w <- sf::st_as_sf(terra::extract(lai, w, bind=T)) # mudar aqui para extrair do raster acima
  w$MOMENTUM_ROUGHNESS_LENGTH <- 0.018 * w$LEAF_AREA_INDEX # retirar esta linha

  # heights
  z1 = 0.1
  z2 = 2.0

  # wind speed at 200 m (eq 32)
  u200 <- w$u * log(200 / w$MOMENTUM_ROUGHNESS_LENGTH) / log(w$z / w$MOMENTUM_ROUGHNESS_LENGTH) # permanence inalterado

  # friction velocity u* (eq 31)
  ap$FRICTION_VELOCITY <- 0.41 * u200 / log(200 / ap$MOMENTUM_ROUGHNESS_LENGTH) # mudar para calcular um raster

  # aerodynamic resistance (eq 28)
  ap$AERODYNAMIC_RESISTANCE <- log(z2/z1) / (0.41*ap$FRICTION_VELOCITY) # mudar para calcular um raster

  # DT é o gradiente de temperatura
  ap$DT = 0
  old_dt <- ap$DT + 1 # aqui se define um valor para old_dt diferente de dt para não comprometer o loop a seguir
  ap$LE <- ap$value

  n <- 0
  while (max(abs(ap$DT - old_dt)) > 1e-6 & n < 100) {
    n <- n + 1
    # air density
    ap$AIR_DENSITY = atmospheric_pressure(ap$ALTITUDE) / (1.01 * (ap$TEMP - ap$DT) * 287)

    cp <- 1307
    ap$H = ap$NET_RADIATION - ap$SOIL_HEAT_FLUX - ap$value
    old_dt <- ap$DT
    ap$DT <- ap$H * ap$AERODYNAMIC_RESISTANCE / (ap$AIR_DENSITY * cp)

    #--------------------
    # monin obukov
    ap$L <- -ap$AIR_DENSITY * cp * (ap$FRICTION_VELOCITY ^ 3) * ap$TEMP / (0.41 * 9.80665 * ap$H)
    psy_200 <- fun_psy200(ap$L)
    psy_z2 <- fun_psy(z2, ap$L)
    psy_z1 <- fun_psy(z1, ap$L)

    #ap_rec[[length(ap_rec) + 1]] <- ap
    ap$FRICTION_VELOCITY <- u200 * 0.41 / (log(200 / ap$MOMENTUM_ROUGHNESS_LENGTH) - psy_200)
    ap$AERODYNAMIC_RESISTANCE <- (log(z2 / z1) - psy_z2 + psy_z1) / (0.41 * ap$FRICTION_VELOCITY)
  }

  ap$TEMP_DATUM <- ap$TEMP - ap$ALTITUDE * 6.5/1000 # lapse-rate de 6.5°C/1000 m
  m <- lm(DT ~ TEMP_DATUM, data = ap)

  #===============================================================================
  # Após o loop
  #===============================================================================

  temp_datum <- model$rast$TEMP - model$rast$ALTITUDE * 6.5 / 1000

  a <- coef(m)[1]
  b <- coef(m)[2]

  # dt for the whole image
  dt <- a + b * temp_datum

  # air density
  rho <- atmospheric_pressure(model$rast$ALTITUDE) / (1.01 * (model$rast$TEMP - dt) * 287)

  roughness_length <- 0.018 * lai

  #===============================
  # new loop
  #===============================
  uf <- 0.41 * u200 / log(200 / roughness_length)
  r_ah <- log(z2 / z1) / (uf * 0.41)

  old_h <- 0
  #===================================
  # starts here
  #==================================
  for (i in 1:10) {
    h <- rho * cp * dt / r_ah
    l <- -rho * cp * (uf ^ 3) * model$rast$TEMP / (0.41 * 9.80665 * h)
    psy_200 <- fun_psy200_r(l)
    psy_z2 <- fun_psy_r(z2, l)
    psy_z1 <- fun_psy_r(z1, l)

    uf <- u200 * 0.41 / (log(200 / roughness_length) - psy_200)
    r_ah <- (log(z2 / z1) - psy_z2 + psy_z1) / (0.41 * uf)
    old_h <- h
  }

  #------------------------------------

  le <- rn - g - h

  #-------------------------------------------------------------------------------

  model$rast$NET_RADIATION <- rn
  model$rast$SOIL_HEAT_FLUX <- g
  model$rast$SENSIBLE_HEAT_FLUX <- h
  model$rast$LATENT_HEAT_FLUX <- le

  model
}
