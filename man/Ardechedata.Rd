% 2014-08-30, Alberto Viglione
%
\name{Ardechedata}
\docType{data}
\alias{Ardechedata}
\alias{Ardeche_areas}
\alias{Ardeche_ungauged_extremes}
\alias{Beauvene_cont}
\alias{Chambonas_cont}
\alias{SaintLaurent_cont}
\alias{SaintMartin_cont}
\alias{SaintMartin_hist}
\alias{Vogue_cont}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Data-sample}
\description{
Systematic flood data and historical flood data for the Ard\`eche region (France) as described in: \\
Naulet, R. (2002). Utilisation de l'information des crues historiques pour une meilleure pr\'ed\'etermination du risque d'inondation. Application au bassin de l'Ard\`eche \`a Vallon Pont-d'Arc et St-Martin d'Ard\`eche. PhD thesis at CEMAGREF, Lyon, and at the Universit\'e du Qu\'ebec, pp. 322. (http://www.lthe.fr/OHM-CV/Documents/theses/these\_naulet.pdf) \\
and \\
Nguyen, C.C. (2012). Am\'elioration des approches Bay\'esiennes MCMC pour l'analyse r\'egionale des crues (Improvement of BayesianMCMC approaches for regional flood frequency analyses). PhD thesis at the University of Nantes, pp. 192. (http://archive.bu.univ-nantes.fr/pollux/show.action?id=178da900-8dd8-491b-8720-011608382b98)
}
\usage{
 data(Ardechedata)
}
\format{
\code{Ardeche_areas}, areas (km2) of the gauged and ungauged catchments in the Ard\`eche region (France); 
\code{Ardeche_ungauged_extremes}, flood peaks (m3/s) reconstructed in ungauged catchments and number of years for which the peak was not exceeded;
\code{Beauvene_cont}, sistematic flood peaks (m3/s) recorded at one station;
\code{Chambonas_cont}, sistematic flood peaks (m3/s) recorded at one station;
\code{SaintLaurent_cont}, sistematic flood peaks (m3/s) recorded at one station;
\code{SaintMartin_cont}, sistematic flood peaks (m3/s) recorded at one station;
\code{SaintMartin_hist}, values for historical peaks (m3/s) for one station and for flood perception thresholds (m3/s) non exceeded in the periods indicated;
\code{Vogue_cont}, sistematic flood peaks (m3/s) recorded at one station.
}
%\source{}
%\references{}
%\author{Olivier Paiastre}
%\note{}

\examples{
data(Ardechedata)
SaintMartin_cont
SaintMartin_hist
}
\keyword{datasets}
