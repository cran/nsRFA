% 2005-09-23, Alberto Viglione
%
\name{FEH1000}
\docType{data}
\alias{FEH1000}
\alias{cd}
\alias{am}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Data-sample}
\description{
Flood Estimation Handbook flood peak data CD-ROM.
}
\usage{
 data(FEH1000)
}
\format{Data.frames:

\code{am} is the data.frame of the annual maximum flows with 3 columns: 
 \code{number}, the code of the station; 
 \code{date}, date of the annual maximum; 
 \code{year}, year of the annual maximum (we consider hydrologic year: 1 october - 30 september);
 \code{am}, the value of the annual maximum flow [m3/s].

\code{cd} is the data.frame of parameters of 1000 catchements with 24 columns: 
 \code{number}, the code of the station;
 \code{nominal_area}, catchment drainage area [km2];
 \code{nominal_ngr_x}, basin outflow coordinates [m];
 \code{nominal_ngr_y}, basin outflow coordinates [m];
 \code{ihdtm_ngr_x}, basin outflow coordinates by Institute of Hydrology digital terrain model [m]; 
 \code{ihdtm_ngr_y}, basin outflow coordinates by Institute of Hydrology digital terrain model [m];
 \code{dtm_area}, catchment drainage area [km2] derived by CEH using their DTM (IHDTM);
 \code{saar4170}, standard average annual rainfall 1941-1970 [mm]; 
 \code{bfihost}, baseflow index derived from HOST soils data;
 \code{sprhost}, standard percentage runoff derived from HOST soils data;
 \code{farl}, index of flood attenuation due to reservoirs and lakes; 
 \code{saar}, standard average annual rainfall 1961-1990 [mm];
 \code{rmed_1d}, median annual maximum 1-day rainfall [mm];
 \code{rmed_2d}, median annual maximum 2-days rainfall [mm];
 \code{rmed_1h}, median annual maximum 1-hour rainfall [mm];
 \code{smdbar}, mean SMD (soil moisture deficit) for the period 1961-1990 calculated from MORECS month-end values [mm];
 \code{propwet}, proportion of time when soil moisture deficit <=6 mm during 1961-90, defined using MORECS;
 \code{ldp}, longest drainage path [km], defined by recording the greatest distance from a catchment node to the defined outlet;
 \code{dplbar}, mean drainage path length [km];
 \code{altbar}, mean catchment altitude [m];
 \code{dpsbar}, mean catchement slope [m/km];
 \code{aspbar}, index representing the dominant aspect of catchment slopes (its values increase clockwise from zero to 360, starting from the north). Mean direction of all inter-nodal slopes with north being zero;
 \code{aspvar}, index describing the invariability in aspect of catchment slopes. Values close to one when all slopes face a similar direction;
 \code{urbext1990}, extent of urban and suburban land cover in 1990 [fraction].

}
\source{http://www.environment-agency.gov.uk/hiflowsuk/}
\references{
Institute of Hydrology (1999) Flood Estimation Handbook, Institute of Hydrology, Oxford.

http://www.environment-agency.gov.uk/hiflowsuk/

Hosking, J.R.M. and Wallis, J.R. (1997) Regional Frequency Analysis: an approach based on L-moments, Cambridge Un
iversity Press, Cambridge, UK.
}

\examples{
data(FEH1000)
names(cd)
am[1:20,]
}
\keyword{datasets}
