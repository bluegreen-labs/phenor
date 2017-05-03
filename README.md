# phenor

Here I combine three of my R packages to facilitate fast and reproducible phenology model development using both PhenoCam and MODIS phenology (MCD12Q2) data.

I combine the query interface in [daymetr](https://github.com/khufkens/daymetr), to access DAYMET climatological data. The [phenocamr](https://github.com/khufkens/phenocamr) package provides easy access to raw PhenoCam Gcc time series and estimated phenological phases. MODIS data is provided through single pixel querying of the ORNL MODIS DAAC using the MODISTools R package. Gathered or processed data is subsequently formated and fit to various phenology models using the phenor package.

## Installation

To install the toolbox in R run the following commands in a R terminal

```R
require(devtools) # load the devtools library
install_bitbucket("khufkens/phenor") # install the package
```

## Use

Example code below shows that in few lines a modelling exercise can be set up. You can either download your own data using the phenocamr pacakge and process / format them correctly using format_phenocam()

```R
# # downloads of PhenoCam data are facilitated using the functions
# # using the phenophase = TRUE flags, transition dates are estimated
# # on the go. The command below downloads all deciduous broadleaf sites
# # and estimates their respective phenophases.
# download.phenocam(vegetation = "DB", phenophase = TRUE)

# process phenocam transition files into a consistent format
phenocam_data = format_phenocam("/foo/bar/transition_dates/")

```

Alternatively you can use the included data which consists of 370 site years of deciduous broadleaf forest sites as included in PhenoCam 1.0 dataset (Richardson et al. 2017). The file is loaded using a standard data() call. The dataset is named *phenocam_DB*.

```R
# load the included data using
data("phenocam_DB")
```

# Model development, and parameter estimation

The gathered data can now be used in model validation. Currently 17 models as described by Basler (2016) are provided in the pacakge. These models include: null, LIN, TT, TTs, PTT, PTTs, M1, M1s, PA, Pab, SM1, SM1b, SQ, SQb, UN, UM1, PM1 and PM1b models as described in Basler (2016). Parameter values associated with the models are provided in a separate file included in the package.

```R
# comma separated parameter file inlcuded in the package
path = sprintf("%s/extdata/parameter_ranges.csv",path.package("phenor"))

par_ranges = read.table(path,
                        header = TRUE,
                        sep = ",")
```

Your own model development can be done by creating similar functions which take in the data as described by the functions above. Below is the function to optimize the model parameters for the *TT* (thermal time) model using Generalized Simulated Annealing (GenSA). Uppper and lower constraints to the parameter space have to be provided, and in case of GenSA initial parameters are estimated when par = NULL.

```R
# optimize model parameters
set.seed(1234)
optim.par = optimize.parameters(par = NULL,
                          data = phenocam_data,
                          cost = cost.function,
                          model = "TT",
                          method = "GenSA",
                          lower = c(0,-20,0,0,0,0,1),
                          upper = c(10,10,100,10,10,1000,365))
```

After a few minutes optimal parameters are provided. We can use these now to calculate our model values and some accuracy metrics.

```R
# now run the model for all data in the nested list using the estimated parameters
modelled = estimate_phenology(data = data, par = optim.par$par)
```

Easy model validation can be achieved using the model_validation() and model_comparison() functions. Which either allow for quick screening for model development, or the comparison of a suite of models (using different starting parameters).

### Dependencies

The code depends on the following R packages: devtools, GenSA, rgenoud, zoo, daymetr and phenocamr
