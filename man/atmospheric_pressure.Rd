% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formulae.R
\name{atmospheric_pressure}
\alias{atmospheric_pressure}
\title{Calculates atmospheric pressure at a given altitude above sea level.}
\usage{
atmospheric_pressure(z, temp = 293.15, lb = 0.0065)
}
\arguments{
\item{z}{Altitude above sea level [m].}

\item{temp}{Meam atmospheric temperature [K]. The default is 293.15.}

\item{lb}{Temperature lapse rate [K / m] (i.e. how many Kelvin the
temperature of air decreases with a 1 m increase in altitude). The default
is 5e-3 K / m.}
}
\value{
Atmospheric pressure at altitude z.
}
\description{
Calculates atmospheric pressure at a given altitude above sea level.
}
