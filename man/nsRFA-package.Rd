\name{nsRFA-package}
\alias{nsRFA-package}
\alias{nsRFA}
\docType{package}
\title{
Non-supervised Regional Frequency Analysis
}
\description{
The estimation of hydrological variables in ungauged basins is a very important topic for many purposes, from research to engineering applications and water management (see the PUB project, Sivapalan et al., 2003). 
Regardless of the method used to perform such estimation, the underlying idea is to transfer the hydrological information from gauged to ungauged sites. 
When observations of the same variable at different measuring sites are available and many data samples are used together as source of information, the methods are called regional methods. 
The well known regional frequency analysis (e.g. Hosking and Wallis, 1997), where the interest is in the assessment of the frequency of hydrological events, belong to this class of methods. 
In literature, the main studied variable is the flood peak and the most used regional approach is the index-flood method of Dalrymple (1960), in which it is implicitly assumed that the flood frequency distribution for different sites belonging to a
homogeneous region is the same except for a site-specific scale factor, the index-flood (see Hosking and Wallis, 1997, for details). 
Hence, the estimation of the flood frequency distribution for an ungauged site can be divided into two parts:
estimation of the index-flood (more in general, the index-value) through linear/non-linear relations with climatic and basin descriptors;
estimation of the adimentional flood frequency distribution, the growth curve, assigning the ungauged basin to one homogeneous region.

\code{nsRFA} is a collection of statistical tools for objective (non-supervised) applications of the Regional Frequency Analysis methods in hydrology.
This does not mean that Regional Frequency Analysis should be non-supervised.
These tools are addressed to experts, to help their expert judgement.
The functions in \code{nsRFA} allow the application of the index-flood method in the following points:
\tabular{ll}{
- \tab regionalization of the index-value;\cr
- \tab formation of homogeneous regions for the growth curves;\cr
- \tab fit of a distribution function to the empirical growth curve of each region;\cr
}



\bold{Regionalization of the index-value}

The index-value can be either the sample mean (e.g. Hosking and Wallis, 1997) or the sample median (e.g. Robson and Reed, 1999) or another scale parameter.
Many methodological approaches are available for the index-value estimation, and their differences can be related to the amount of information available.
Excluding direct methods, that use information provided by flow data available at the station of interest, regional estimation methods require ancillary hydrological and physical information.
Those methods can be divided in two classes: the multiregressive approach and the hydrological simulation approach.
For both of them, the best estimator is the one that optimizes some criterion, such as the minimum error, the minimum variance or the maximum efficiency.
Due to its simplicity, the most frequently used method is the multiregressive approach (see e.g. Kottegoda & Rosso, 1998; Viglione et al., 2007a), that relates the index-flow to catchment characteristics, such as climatic indices, geologic and morphologic parameters, land cover type, etc., through linear or non-linear equations.

R provides the functions \code{\link{lm}} and \code{\link{nls}} for linear and non-linear regressions (package \code{\link{stats}}).
With the package \code{nsRFA}, a tool to select the best linear regressions given a set of candidate descriptors, \code{\link{bestlm}}, is provided.
In \code{\link{REGRDIAGNOSTICS}} several functions are provided to analyze the output of \code{\link{lm}}, such as: 
the coefficient of determination (classic and adjusted);
the Student t significance test;
the variance inflation factor (VIF);
the root mean squared error (RMSE);
the mean absolute error (MAE);
the prediction intervals;
predicted values by a jackknife (cross-validation) procedure.
The functions in \code{\link{DIAGNOSTICS}} provide more general diagnostics of model results (that can be also non-linear), comparing estimated values with observed values.

More details are provided in vignettes:
\tabular{ll}{
\code{nsRFA_ex01} \tab How to use the package nsRFA (example 1):\cr
                  \tab Regional frequency analysis of the annual flows in Piemonte\cr
                  \tab and Valle d'Aosta\cr
}
that can be accessed via \code{vignette("nsRFA_ex01", package="nsRFA")}.



\bold{Formation of homogeneous regions for the growth curves}

Different techniques exist, for example those that lead to the formation of fixed regions through cluster analysis (Hosking and Wallis, 1997, Viglione, 2007), or those based on the method of the region of influence (ROI, Burn, 1990).
The regional procedure can be divided into two parts: the formation of regions and the assignment of an ungauged site to one
of them.
Concerning the first point, the sites are grouped according to their similarity in terms of those basin descriptors that are assumed to explain the shape of the growth curve.
This shape is usually quantified in a parametric way.
For instance, the coefficient of variation (CV) or the L-CV of the curve can be used for this purpose.
The package \code{nsRFA} provide the functions in \code{\link{moments}} and \code{\link{Lmoments}} to calculate sample moments and L-moments.
Subsequently, the selected parameter is related with basin descriptors through a linear or a more complex model.
A regression analysis is performed with different combination of descriptors, and descriptors that are strongly related with the parameter are used to group sites in regions.
The same tools used for the regionalization of the index value, i.e. \code{\link{bestlm}}, \code{\link{REGRDIAGNOSTICS}} and \code{\link{DIAGNOSTICS}}, can be used if the parametric method is chosen.

\code{nsRFA} also provide a non-parametric approach that considers the dimensionless growth curve as a whole (see, Viglione et al., 2006; Viglione, 2007). 
The multiregressive approach can still be used if we reason in terms of (dis)similarity between pairs of basins in the following way:
(1) for each couple of stations, a dissimilarity index between non-dimensional curves is calculated using a quantitatively predefined metric, for example using the Anderson-Darling statistic (\code{\link{A2}}), and organising the distances in a matrix with \code{\link{AD.dist}};
(2) for each basin descriptor, the absolute value (or another metric) of the difference between its measure in two basins is used as distance between them, using \code{\link{dist}} of the package \code{\link{stats}} to obtain distance matrices;
(4) a multiregressive approach (\code{\link{bestlm}}, \code{\link{lm}}) is applied considering the matrices as variables and the basin descriptors associated to the best regression are chosen;
(5) the significance of the linear relationship between distance matrices is assessed through the Mantel test with \code{\link{mantel.lm}}.

In the suitable-descriptor's space, stations with similar descriptor values can be grouped into disjoint regions through a cluster analysis (using functions in \code{\link{traceWminim}}) or the ROI method can be used adapting a region to the ungauged basin (\code{\link{roi}}).
In both cases, the homogeneity of the regions can be assessed with the functions in \code{\link{HOMTESTS}}, where the Hosking and Wallis heterogeneity measures (\code{\link{HW.tests}}, see Hosking and Wallis, 1997) and the Anderson-Darling homogeneity test (\code{\link{ADbootstrap.test}}, see Viglione et al., 2007b) are provided.

More details are provided in vignettes:
\tabular{ll}{
\code{nsRFA_ex01} \tab How to use the package nsRFA (example 1):\cr
                  \tab Regional frequency analysis of the annual flows in Piemonte\cr
                  \tab and Valle d'Aosta\cr
\code{nsRFA_ex02} \tab How to use the package nsRFA (example 2):\cr
                  \tab Region-Of-Influence approach, some FEH examples\cr
}
that can be accessed via \code{vignette("nsRFA_ex01", package="nsRFA")}.





\bold{Fit of a distribution function to the empirical growth curve of each region}

Once an homogeneous region is defined, the empirical growth curves can be pooled together and a probability distribution can be fitted to the pooled sample.
The choice of the best distribution can be assisted by a Model Selection Criteria with \code{\link{MSClaio2008}} (see, Laio et al., 2008).
The parameters of the selected distribution can be estimated using the method of moments (\code{\link{moment_estimation}}), L-moments (\code{\link{par.GEV}}, \code{\link{par.genpar}}, \code{\link{par.gamma}}, ...) or maximum-likelihood (\code{\link{MLlaio2004}}).
Goodness-of-fit tests are also available: the Anderson-Darling goodness of fit test with \code{\link{GOFlaio2004}} (Laio. 2004), and Monte-Carlo based tests with \code{\link{GOFmontecarlo}}.
Confidence intervals for the fitted distribution can be calculated with a Markov Chain Monte Carlo algorithm, using \code{\link{BayesianMCMC}}.

More details are provided in vignettes:
\tabular{ll}{
\code{nsRFA_ex01} \tab How to use the package nsRFA (example 1):\cr
                  \tab Regional frequency analysis of the annual flows in Piemonte\cr
                  \tab and Valle d'Aosta\cr
\code{MSClaio2008} \tab Model selection techniques for the frequency analysis\cr 
                   \tab of hydrological extremes: the MSClaio2008 R function\cr
}
that can be accessed via \code{vignette("nsRFA_ex01", package="nsRFA")}.





% ------------------------------------------------------------ %

\bold{Other functions}

\code{\link{varLmoments}} provides distribution-free unbiased estimators of the variances and covariances of sample L-moments, as described in Elamir and Seheult (2004).

More details are provided in vignettes:
\tabular{ll}{
\code{Fig1ElamirSeheult} \tab Figure 1 in Elamir and Seheult (2004)\cr
}


}	% end of description
\details{
\tabular{ll}{
Package: \tab nsRFA\cr
Version: \tab 0.7\cr
}
The package provides several tools for Regional Frequency Analysis of hydrological variables.
The first version dates to 2006 and was developed in Turin at the Politecnico by Alberto Viglione.

For a complete list of the functions, use \code{library(help="nsRFA").}

\bold{Main changes in version 0.7}
\tabular{ll}{
0.7-2: \tab ???;\cr
0.7-1: \tab refinement of function \code{\link{BayesianMCMC}};\cr
0.7-0: \tab plotting position for historical information in \code{\link{DISTPLOTS}};\cr
}
\bold{Main changes in version 0.6}
\tabular{ll}{
0.6-9: \tab new vignette \code{Fig11GriffisStedinger};\cr
0.6-8: \tab exponential and Gumbel distributions added in \code{\link{GOFmontecarlo}};\cr
0.6-6: \tab some plotting position/probability plots have been added in \code{\link{DISTPLOTS}};\cr
0.6-4: \tab refinement of function \code{\link{BayesianMCMC}};\cr
0.6-2: \tab new vignette \code{nsRFA_ex02};\cr
0.6-2: \tab refinement of function \code{\link{BayesianMCMC}};\cr
0.6-0: \tab new vignette \code{nsRFA_ex01};\cr
0.6-0: \tab new function \code{\link{bestlm}};\cr
0.6-0: \tab the plotting position/probability plots in \code{\link{DISTPLOTS}} have been reshaped;\cr
0.6-0: \tab this list of changes has been added;\cr
}
}
\author{
Alberto Viglione

Maintainer: Alberto Viglione <viglione@hydro.tuwien.ac.at>
}
\references{
All the manual references are listed here:

Beirlant, J., Goegebeur, Y., Teugels, J., Segers, J., 2004. Statistics of 
  Extremes: Theory and Applications. John Wiley and Sons Inc., 490 pp.

Burn, D.H., 1990. Evaluation of regional flood frequency analysis with a
  region of influence approach. Water Resources Research 26(10), 2257-2265.

Castellarin, A., Burn, D.H., Brath, A., 2001. Assessing the effectiveness of
  hydrological similarity measures for flood frequency analysis. Journal of
  Hydrology 241, 270-285.

Chowdhury, J.U., Stedinger, J.R., Lu, L.H., Jul. 1991. Goodness-of-fit tests
  for regional generalized extreme value flood distributions. Water Resources
  Research 27(7), 1765-1776.

D'Agostino, R., Stephens, M., 1986. Goodness-of-fit techniques. Vol.68 of
  Statistics: textBooks and monographs. Department of Statistics, Southern
  Methodist University, Dallas, Texas.

Dalrymple, T., 1960. Flood frequency analyses. Vol. 1543-A of Water Supply
  Paper. U.S. Geological Survey, Reston, Va.

Durbin, J., Knott, M., 1971. Components of Cramer-von Mises statistics. 
  London School of Economics and Political Science, pp. 290-307.

El Adlouni, S., Bob\'ee, B., Ouarda, T.B.M.J., 2008. On the tails of extreme 
  event distributions in hydrology. Journal of Hydrology 355(1-4), 16-33.

Elamir, E. A.H., Seheult, A.H., 2004. Exact variance structure of sample
  l-moments. Journal of Statistical Planning and Inference 124, 337-359.

Everitt, B., 1974. Cluster Analysis. Social Science Research Council. 
  Halsted Press, New York.

Fill, H., Stedinger, J., 1998. Using regional regression within index flood
  procedures and an empirical bayesian estimator. Journal of Hydrology
  210(1-4), 128-145.

Greenwood, J., Landwehr, J., Matalas, N., Wallis, J., 1979. Probability
  weighted moments: Definition and relation to parameters of several
  distributions expressible in inverse form. Water Resources Research 15,
  1049-1054.

Hosking, J., 1986. The theory of probability weighted moments. Tech. Rep.
  RC12210, IBM Research, Yorktown Heights, NY.

Hosking, J., 1990. L-moments: analysis and estimation of distributions using
  linear combinations of order statistics. J. Royal Statistical Soc. 52,
  105-124.

Hosking, J., Wallis, J., 1993. Some statistics useful in regional frequency
  analysis. Water Resources Research 29(2), 271-281.

Hosking, J., Wallis, J., 1997. Regional Frequency Analysis: An Approach Based
  on L-Moments. Cambridge University Press.

Institute of Hydrology, 1999. Flood Estimation Handbook, Institute of 
  Hydrology, Oxford.

Kendall, M., Stuart, A., 1961-1979. The Advanced Theory of Statistics.
  Charles Griffin \& Company Limited, London.

Kottegoda, N.T., Rosso, R., 1997. Statistics, Probability, and Reliability for
  Civil and Environmental Engineers, international Edition. McGraw-Hill
  Companies.

Laio, F., 2004. Cramer-von Mises and Anderson-Darling goodness of fit tests for 
  extreme value distributions with unknown parameters, Water Resour. Res., 
  40, W09308, doi:10.1029/2004WR003204.

Laio, F., Tamea, S., 2007. Verification tools for probabilistic forecast of
  continuous hydrological variables. Hydrology and Earth System Sciences 11,
  1267-1277.

Laio, F., Di Baldassarre G., Montanari A., 2008. Model selection techniques 
  for the frequency analysis of hydrological extremes, Water Resour. Res., 
  Under Revision.

Regione Piemonte, 2004. Piano di tutela delle acque. Direzione Pianificazione 
  Risorse Idriche.

Robson, A., Reed, D., 1999. Statistical procedures for flood frequency
  estimation. In: Flood Estimation HandBook. Vol.~3. Institute of Hydrology
  Crowmarsh Gifford, Wallingford, Oxfordshire.

Sankarasubramanian, A., Srinivasan, K., 1999. Investigation and comparison of
  sampling properties of l-moments and conventional moments. Journal of
  Hydrology 218, 13-34.

Scholz, F., Stephens, M., 1987. K-sample Anderson-Darling tests. Journal of 
  American Statistical Association, 82(399), pp. 918-924.

Sivapalan, M., Takeuchi, K., Franks, S.W., Gupta, V.K., Karambiri, H., Lakshmi, V.,
  Liang, X., McDonnell, J.J., Mendiondo, E.M., O'Connell, P.E., Oki, T., Pomeroy, J.W,
  Schertzer, D., Uhlenbrook, S., Zehe, E., 2003. IAHS decade on Predictions in 
  Ungauged Basins (PUB), 2003-2012: Shaping an exciting future for the hydrological 
  sciences, Hydrological Sciences - Journal - des Sciences Hydrologiques, 48(6), 
  857-880.

Smouse, P.E., Long, J.C., Sokal, R.R., 1986. Multiple regression and correlation 
  extensions of the Mantel test of matrix correspondence. Systematic Zoology, 
  35(4), 627-632.

Stedinger, J., Lu, L., 1995. Appraisal of regional and index flood quantile
  estimators. Stochastic Hydrology and Hydraulics 9(1), 49-75.

Stedinger, J.R., Vogel, R.M. and Foufoula-Georgiou, E., 1993. Frequency analysis of extreme events. 
  In David R. Maidment, editor, Hand-Book of Hydrology, chapter 18. McGraw-Hill Companies, 
  international edition.

Viglione, A., Claps, P., Laio, F., 2006. Utilizzo di criteri di prossimit\`a 
  nell'analisi regionale del deflusso annuo, XXX Convegno di Idraulica e 
  Costruzioni Idrauliche - IDRA 2006, Roma, 10-15 Settembre 2006.

Viglione, A., 2007. Metodi statistici non-supervised per la stima di grandezze 
  idrologiche in siti non strumentati, PhD thesis at the Politecnico of Turin,
  February 2007.

Viglione, A., Claps, P., Laio, F., 2007a. Mean annual runoff estimation in 
  North-Western Italy, In: G. La Loggia (Ed.) Water resources assessment 
  and management under water scarcity scenarios, CDSU Publ. Milano. 

Viglione, A., Laio, F., Claps, P., 2007b. A comparison of homogeneity tests
  for regional frequency analysis. Water Resources Research 43~(3).

Vogel, R., Fennessey, N., 1993. L moment diagrams should replace product moment
  diagrams. Water Resources Research 29(6), 1745-1752.

Vogel, R., Wilson, I., 1996. Probability distribution of annual maximum, mean,
  and minimum streamflows in the united states. Journal of Hydrologic
  Engineering 1(2), 69-76.

Ward, J., 1963. Hierarchical grouping to optimize an objective function, Journal 
  of the American Statistical Association, 58, pp. 236-244.
}
\keyword{package}
\keyword{htest}
\keyword{distribution}
\keyword{univar}
\keyword{hplot}
%\seealso{}
%\examples{}
