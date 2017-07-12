# phenor

The phenor R package is part of a phenology modelling framework in R leveraging measurements of vegetation phenology from the PhenoCam network and existing retrospective and projected climate data (e.g. Daymet, CMIP5, Berkeley Earth and others) in support of model comparisons and model development.

The package is curently focussed on contiguous US and relies heavily on Daymet climate data queried using [daymetr](https://github.com/khufkens/daymetr) for model optimization. Phenological model validation data is primarily derived from PhenoCam time series (i.e. through [phenocamr](https://github.com/khufkens/phenocamr)). But functions are included to ingest different phenology data sources such as the MODIS MCD12Q2 phenology product using the MODISTools R package. We refer to Hufkens et al. (2017, see below) for an in depth description and a worked example of the phenor R package.

For full transparency the code used to generate the referenced publication is provided in the analysis folder of the package. This folder is not compiled into the R package using the installation instructions below. To run these scripts please clone the github repository. Keep in mind that some of the scripts will take a significant amount of time to finish. As such, some data generated for the manuscript is included in the R package. Some scripts to generate model comparison figures and summary statistics therefore rely on these precompiled datasets rather than clean runs.

## Installation

To install the toolbox in R run the following commands in a R terminal

```R
if(!require(devtools)){install.package(devtools)}
devtools::install_github("khufkens/phenor")
library(phenor)
```
Additional PhenoCam data (Richarson et al. 2017) can be downloaded as:

```
git clone https://github.com/khufkens/phenocam_dataset.git
```

## Use

Example code below shows that in few lines a modelling exercise can be set up. You can either download your own data using the phenocamr pacakge and format them correctly using **format_phenocam()**

```R
# The command below downloads all time series for deciduous broadleaf
# data at the bartlett PhenoCam site and estimates the
# phenophases.
download_phenocam(vegetation = "DB",
                  site = "bartlett",
                  phenophase = TRUE)

# process phenocam transition files into a consistent format
phenocam_data = format_phenocam("/foo/bar/transition_dates/")
```

Alternatively you can use the included data which consists of 370 site years of deciduous broadleaf forest sites as included in PhenoCam 1.0 dataset (Richardson et al. 2017). The file is loaded using a standard **data()** call or reference them directly (e.g. *phenocam_DB*).

```R
# load the included data using
data("phenocam_DB")
```

### Model development, and parameter estimation

The gathered data can now be used in model validation. Currently 17 models as described by Basler (2016) are provided in the pacakge. These models include: null, LIN, TT, TTs, PTT, PTTs, M1, M1s, PA, Pab, SM1, SM1b, SQ, SQb, UN, UM1, PM1 and PM1b models as described in Basler (2016). Parameter values associated with the models are provided in a separate file included in the package.

```R
# comma separated parameter file inlcuded in the package
# for all the included models, this file is used by default
# in all optimization routines
path = sprintf("%s/extdata/parameter_ranges.csv",path.package("phenor"))
par_ranges = read.table(path,
                        header = TRUE,
                        sep = ",")
```

Your own model development can be done by creating similar functions which take in the data as described by the functions above. Below is the function to optimize the model parameters for the *TT* (thermal time) model using Generalized Simulated Annealing (GenSA). Uppper and lower constraints to the parameter space have to be provided, and in case of GenSA initial parameters are estimated when par = NULL.

```R
# optimize model parameters
set.seed(1234)
optim.par = optimize_parameters(par = NULL,
                          data = phenocam_data,
                          cost = rmse,
                          model = "TT",
                          method = "GenSA",
                          lower = c(1,-5,0),
                          upper = c(365,10,2000))
```

After a few minutes optimal parameters are provided. We can use these now to calculate our model values and some accuracy metrics.

```R
# now run the model for all data in the nested list using the estimated parameters
modelled = estimate_phenology(data = data,
                              par = optim.par$par)
```

Easy model validation can be achieved using the **model_validation()** and **model_comparison()** functions. Which either allow for quick screening for model development, or the comparison of a suite of models (using different starting parameters).

## References

Hufkens K., Basler J. D., Milliman T. Melaas E., Richardson A.D. 2017 An integrated phenology modelling framework in R: Phenology modelling with phenor. In Review

Richardson, A.D., Hufkens, K., Milliman, T., Aubrecht, D.M., Chen, M., Gray, J.M., Johnston, M.R., Keenan, T.F., Klosterman, S.T., Kosmala, M., Melaas, E.K., Friedl, M.A., Frolking, S. 2017. Tracking vegetation phenology across diverse North American biomes using PhenoCam imagery. In Review.

## Acknowledgements

This project was is supported by the National Science Foundation’s Macro-system Biology Program (award EF-1065029).
