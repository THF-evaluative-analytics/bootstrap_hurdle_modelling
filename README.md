# Bootstrapping effect estimates from Hurdle models

#### Project Status: [In progess]


## Project Description

- Much of the IAU's analytical work focusses on quantifying the impact of selected health-care improvement programmes on key outcome measures evidenced by observational studies.  Conventional statistical outcome models may fall short of adequately capturing relevant features of the distribution across comparator groups of selected outcomes (notably length of hospital stay), thereby casting doubts on the reliability of ensuing effect estimates.  Two-part models (e.g. of the Hurdle variety as illustrated in Cameron and Trivedi, 2013, Sec. 4.5) are especially suited to describe discrete distributions exhibiting a mode at 0 which is too marked to be satisfactorily accommodated by standard statistical models.  On the other hand, when a Hurdle model is utilised to fit count outcomes in the context of a programme evaluation, separate sets of inferences around impact are obtained: one for strictly positive outcome counts, and another for the ocurrence of null instead of positive counts.  A complication though arises in that a single measure of intervention impact, which isn't available in standard software implementations, is normally desired to usefully inform policy-making.

- This project aims at providing stand-alone R code implementing a bootstrap approach (Efron and Tibshirani, 1993) for pooling effect estimates from the positive and zero parts of a Hurdle model fitted to a rectangular data-set including a binary intervention indicator.

- Provided code will fit a user-defined Hurdle model to a pre-specified count outcome, a set of covariates and a binary intervention assignment indicator, bootstrap resulting effect estimates to obtain corresponding confidence intervals and a p-value for the null hypothesis of no intervention impact.  All effect estimates are averaged across the distribution of all covariates to remove their dependency on any covariate value, thus ensuring their meaningful interpretation.

- GLMs for the zero (modelled via Logistic regression) and positive (represented by either a Poisson, Negative Binomial or Geometric regression) parts making up the Hurdle model are separately fitted via maximum likelihood; 95% confidence intervals around effect estimates are subsequently built via either accelerated bias-corrected or percentile-based boostrap.


## Data source

Data are extracted, collated and processed from a variety of sources (chiefly but not limited to SUS) by data managers and analysts in the IAU.  Due to data-sharing agreements stipulated by the IAU, none of the final data-sets the code is intended for use on are publicly available.


## How does it work?

(i)  A rectangular data-set (one row per case, one column per variable) inclusive of an individual-level intervention assignment binary indicator.

What you need to do to reproduce the analysis or re-use the code on your local machine.  

### Requirements

This script was written in R version 3.6.2 and Emacs version 26.3. 

The following R packages (available on CRAN) are needed: 

* [boot] (https://cran.r-project.org/package=boot)
* [pscl] (https://cran.r-project.org/package=pscl)


### Getting started

Describe the way in which the code can be used. 


## Useful references

--  A. C. Cameron and P. K. Trivedi.  Regression Analysis of Count Data.  Econometric Society Monographs.  Cambridge University Press, NY, 2nd edition, 2013.  doi:  10.1007/978-1-4757-3692-2

--  B. Efron and R. J. Tibshirani.  An Introduction to the Bootstrap.  Monographs on Statistics and Applied Probability.  Springer Science+Business Media, New York, NY, 1993.  doi:  10.1007/978-1-4899-4541-9

--  H. Herwartz, N. Klein, and C. Strumann.  Modelling Hospital Admission  and  Length  of  Stay  by  Means  of  Generalised  Count  Data  Models.  Journal of Applied Econometrics, 31(6):1159–1182, 2016.  doi:10.1002/jae.2454


## Authors

* Stefano Conti - [e. stefano.conti@nhs.net] - [https://github.com/sconti555]


## License

This project is licensed under the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl-3.0.html).


## Acknowledgments

* Current and past colleagues in the IAU for their existing coding work, encouragement, constructive views and critical comments.
* CRAN and GNU for allowing work on this project to be carried out exclusively from freely available, open-source software.
