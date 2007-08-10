% 2005-09-23, Alberto Viglione
%
\name{REGRDIAGNOSTICS}
\alias{REGRDIAGNOSTICS}
\alias{R2.lm}
\alias{prt.lm}
\alias{mantel.lm}
\alias{vif.lm}
\alias{RMSE.lm}
\alias{MAE.lm}
\alias{predinterval.lm}
\alias{jackknife1.lm}
\alias{RMSEjk.lm}
\alias{MAEjk.lm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Diagnostics of regressions}
\description{
 Diagnostics of the output of \code{lm}, that is used to fit linear models. 
}
\usage{
 R2.lm (x)
 prt.lm (x)
 mantel.lm (x, Nperm = 1000)
 vif.lm (x)
 RMSE.lm (x) 
 MAE.lm (x)
 predinterval.lm (x, level = 0.95)
 jackknife1.lm (x)
 RMSEjk.lm (x)
 MAEjk.lm (x)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{object of class ``lm'' (output of `lm')}
 \item{Nperm}{number of permutations}
 \item{level}{significance level}
}
%\details{}
\value{
 \code{R2.lm} returns the coefficient of determination \eqn{R^2} and the adjusted coefficient of determination \eqn{R^2_{adj}} of the regression.

 \code{prt.lm} returns the probability \eqn{Pr(>|t|)} of the significance test (Student t) of the independent variables.
 If the value is 0.06 for a regressor, its coefficient is not significantly different from 0 for a test with significance level of 5\%. 

 \code{mantel.lm} returns the probability \eqn{P} of the Mantel test on every variable conditionated to the others.
 It substitutes \code{prt.lm} when dealing with distance matrices.
 If \eqn{P} is, for example, 0.92, the variable should be considered significant with significance level greater of 8\%.

 \code{vif.lm} returns the variance inflation factors (VIF) of the independent values of the regression. If \eqn{VIF > 5} (or 10) there is a problem of multicollinearity.

 \code{RMSE.lm} returns the root mean squared error of the regression.

 \code{MAE.lm} returns the mean absolute error of the regression.

 \code{predinterval.lm} returns the prediction intervals at a specified \code{level} in correspondence to the fitted data.

 \code{jackknife1.lm} returns predicted values by a jackknife (cross-validation) procedure. 
 The procedure (remove 1 observation, fit the model, estimate in the removed point) is repeated for all the points.

 \code{RMSEjk.lm} returns the root mean squared error of the cross-validation (performed with \code{jackknife1.lm}).

 \code{MAEjk.lm} returns the mean absolute error of the cross-validation (performed with \code{jackknife1.lm}).
}
\references{
 Viglione A., Claps P., Laio F. (2006) Utilizzo di criteri di prossimit\`a nell'analisi regionale del deflusso annuo, XXX Convegno di Idraulica e Costruzioni Idrauliche - IDRA 2006, Roma, 10-15 Settembre 2006.

 Viglione A., Claps P., Laio F. (2006) Water resources assessment and management under water scarcity scenarios, chapter Meanannual runoff estimation in North-Western Italy. CDSU, Milan.

}
\author{Alberto Viglione, e-mail: \email{alviglio@tiscali.it}.}

%\note{.}

\seealso{\code{\link{lm}}, \code{\link{summary.lm}}, \code{\link{predict.lm}}}

\examples{
data(hydroSIMN)

D <- annualflows["dato"][,]
cod <- annualflows["cod"][,]

#Dm <- tapply(D,cod,mean)
#datregr <- cbind(Dm,parameters)
datregr <- parameters
regr0 <- lm(Dm ~ .,datregr); summary(regr0)
regr1 <- lm(Dm ~ Am + Hm + Ybar,datregr); summary(regr1)

R2.lm(regr0)
R2.lm(regr1)

prt.lm(regr0)
prt.lm(regr1)

vif.lm(regr0)
vif.lm(regr1)

RMSE.lm(regr0)
RMSE.lm(regr1)

MAE.lm(regr0)
MAE.lm(regr1)

predinterval.lm(regr0)

jackknife1.lm(regr0)
jackknife1.lm(regr1)

RMSEjk.lm(regr0)
RMSEjk.lm(regr1)

MAEjk.lm(regr0)
MAEjk.lm(regr1)

# mantel test on distance matrices
#Y <- AD.dist(D,cod)             # it takes some time
#X <- data.frame(apply(datregr[,c("Hm","Ybar")],2,dist))
#dati <- cbind(as.numeric(Y),X)
#modello <- lm(Y ~ Hm + Ybar, dati)
#mantel.lm(modello, Nperm=100)

}
\keyword{regression}
