\name{tetens}
\alias{tetens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates water vapour saturation pressure as a function of temperature.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{tetens(temp)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{temp}{Temperature in celsius (ºC)}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the water vapour saturation pressure in kPa.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{João Vitor de Nóvoa Pinto}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (temp)
{
  0.6108 * exp(17.27 * temp / (237.3 + temp))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
