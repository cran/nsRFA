% 2005-09-23, Alberto Viglione
%
\name{exampleRFA}
\alias{exampleRFA}
\alias{exampleRFA01}
\alias{exampleRFA02}
\alias{exampleRFA03}
\alias{exampleRFA04}
\alias{exampleRFA05}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Example of package functionality}
\description{
 Examples of Regional Frequency Analysis.
}
\usage{
 exampleRFA01 ()
 exampleRFA02 ()
 exampleRFA03 ()
 exampleRFA04 ()
 exampleRFA05 ()
}
%- maybe also `usage' for other objects documented here.
%\arguments{\item{}{}}
}
%\details{}
\value{\code{exampleRFA01} and \code{exampleRFA02} are about the Regional Frequency Analysis of annual flows in Piemonte and Valle d'Aosta.
\code{exampleRFA01} treats the regionalization of the index-flow.
\code{exampleRFA02} treats the regionalization of the growth-curve.

\code{exampleRFA03} deals with the formation of pooling groups using the Region of Influence (ROI) approach with the FEH data.

\code{exampleRFA04} deals with the estimation of the index-flood with the FEH data.

\code{exampleRFA05} deals with the determination of the classification variables with the FEH data.
}
\references{
Viglione A., Claps P., Laio F. (2006) Utilizzo di criteri di prossimit\`a nell'analisi regionale del deflusso annuo, XXX Convegno di Idraulica e Costruzioni Idrauliche - IDRA 2006, Roma, 10-15 Settembre 2006.

Viglione A., Claps P., Laio F. (2006) Water resources assessment and management under water scarcity scenarios, chapter Mean annual runoff estimation in North-Western Italy. CDSU, Milan.

Viglione A. (2007) Metodi statistici non-supervised per la stima di grandezze idrologiche in siti non strumentati, PhD thesis
, In press.

Institute of Hydrology (1999) Flood Estimation Handbook, Institute of Hydrology, Oxford.

http://www.environment-agency.gov.uk/hiflowsuk/

Hosking, J.R.M. and Wallis, J.R. (1997) Regional Frequency Analysis: an approach based on L-moments, Cambridge University Press, Cambridge, UK.

}

\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

%\note{.}

%\seealso{}

\examples{
exampleRFA01()
#exampleRFA02()          # it takes some time
#exampleRFA03()          # it takes some time
exampleRFA04()
exampleRFA05()
}
\keyword{manip}
