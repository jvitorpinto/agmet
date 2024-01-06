
Agmet is a collection of functions to calculate solar position,
reference evapotranspiration, and soil water balance.

The following document describes the procedures to install Agmet
directly from GitHub as well as some of the functionalities embedded in
this package.

## Package’s installation

To install Agmet in your computer with R already installed do the
following

``` r
install.packages('devtools')
```

    ## Installing package into '/home/jvitorpinto/R/x86_64-pc-linux-gnu-library/4.3'
    ## (as 'lib' is unspecified)

``` r
devtools::install_github('https://github.com/jvitorpinto/agmet.git')
```

    ## Skipping install of 'agmet' from a github remote, the SHA1 (48fb206e) has not changed since last install.
    ##   Use `force = TRUE` to force installation

## Using Agmet

The function `stefan_boltzmann_law` calculates the amount of radiation
emitted by a black body as a function of its absolute temperature.

``` r
temp <- seq(273.15, 473.15, by = 0.1)
y <- stefan_boltzmann_law(temp)

plot(temp, y, type = 'l', xlab = expression(italic(T) ~ ('K')), ylab = expression('Black body emittance' ~ (W ~ m^-2)))
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->
