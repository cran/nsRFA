% 2005-09-23, Alberto Viglione
%
\name{hydroSIMN}
\docType{data}
\alias{hydroSIMN}
\alias{annualflows}
\alias{parameters}
\alias{meanmonthlyflows}
\alias{monthlyflows}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Data-sample}
\description{
SIMN (Servizio Idrografico e Mareografico Nazionale) flow data samples and catchment parameters.
}
\usage{
 data(hydroSIMN)
}
\format{Data.frames:

\code{annualflows} is the data.frame of the annual flows with 3 columns: 
\code{cod}, the code of the station; 
\code{anno}, the year; 
\code{dato}, the value of the annual flow [mm].

\code{parameters} is the data.frame of parameters of 47 catchements with 16 columns: 
\code{cod}, the code of the station;
\code{Dm}, the mean annual streamflow [mm] as reported in the \sQuote{Pubblicazione n. 17};
\code{Am}, the mean annual rainfall [mm] as reported in the \sQuote{Pubblicazione n. 17};
\code{S}, area of the plane projection of the drainage basin [km2];
\code{Hm}, mean elevation of the drainage basin [m a.s.l.]; 
\code{Pm}, mean slope of the basin [\%]:
 \deqn{P_m=arctg(2(H_{med} - H_{min})/\sqrt{S})} 
 where \eqn{S} is the basin area, \eqn{H_{med}} the median elevation and \eqn{H_{min}} the elevation of the closing section. 
 \code{Pm} is a slope measure of a square equivalent basin, and does not account for basin shape;
\code{LLDP},  length of the longest drainage path [km]. 
 The longest drainage path is the longest path between the basin outlet and the most distant point on the basin border, following drainage directions. 
 Actually the longest drainage path corresponds to the main stream plus the path on the hillslope that connects the stream source to the basin border;
\code{PLDP}, slope of the longest drainage path [\%]. Average of the slope values associated to each pixel in the longest drainage path;
\code{S2000}, area above 2000 m a.s.l. [\%];
\code{EST}, \sQuote{easting}, sine of the angle between the segment connecting the center of mass and the outlet of the basin and the north. 
 \code{EST} is 1 if the basin is oriented eastward, -1 if it is oriented westward;
\code{NORD}, \sQuote{northing}, cosine of the angle between the segment connecting the center of mass and the outlet of the basin and the north.
 \code{NORD} is 1 if the basin is oriented northward, -1 if it is oriented southward;
\code{Rc}, circularity ratio Rc.
 Ratio between the basin area and the area of a circle having the same perimeter:
 \deqn{R_c = \frac{4 \pi S}{P^2}}
 where \eqn{P} is the watershed perimeter;
\code{Xbar}, longitude [deg] of the centroid of the plane projection of the drainage basin;
\code{Ybar}, latitude [deg] of the centroid of the plane projection of the drainage basin;
\code{IT}, Thornthwaite index:  a global moisture index that can be estimated, in its simplest form, as the ratio
 \deqn{I_T=\frac{A_m - ET_0}{ET_0}}
 where \eqn{ET_0} is the mean annual potential evapotranspiration on the basin;
\code{IB}, Budyko index: a radiational aridity index expressed as
 \deqn{I_B=\frac{R_n}{\lambda A_m}}
 where \eqn{R_n} is the mean annual net radiation and \eqn{\lambda} is the latent vaporization heat. 
 Values assumed by \eqn{I_B} are lower than 1 for humid regions and greater than 1 in arid regions.

\code{meanmonthlyflows} is the data.frame of the mean monthly streamflows [mm] as reported in the \sQuote{Pubblicazione n. 17}.
It has 13 columns because the first one, \code{cod}, is the code of the station.

\code{monthlyflows} is the data.frame of the monthly streamflows [mm] with 4 columns:
\code{cod}, the code of the station;
\code{anno}, the year;
\code{mese}, the month;
\code{dato}, the value of the annual flow [mm].

}
%\source{}

%\references{}
%\author{Alberto Viglione}
\note{For information on the package and the Author, and for all the references, see \code{\link{nsRFA}}.}

\examples{
data(hydroSIMN)
annualflows
summary(annualflows)
x <- annualflows["dato"][,]
cod <- annualflows["cod"][,] 
split(x,cod)
sapply(split(x,cod),mean)
sapply(split(x,cod),median)
sapply(split(x,cod),quantile)
sapply(split(x,cod),Lmoments)

parameters
}
\keyword{datasets}
