# PhenoR

A toolbox to model phenological responses with the standard spring warming and sequential models as described in Melaas et al. 2013 (AFM).

## Installation

To install the toolbox in R run the following commands in a R terminal

	require(devtools) # load the devtools library
	install_bitbucket("khufkens/PhenoR") # install the package
	
## Use

# Fast PhenoCam & MODIS model integration.

Here I combine three of my R packages to facilitate fast and reproducible phenology model development using both PhenoCam and MODIS phenology (MCD12Q2) data.

I combine the query interface in daymetr, to access DAYMET climatological data. The phenocamr package provides easy access to raw PhenoCam Gcc time series and estimated phenological phases. MODIS data is provided through single pixel querying of the ORNL MODIS DAAC. Gathered and / or processed data is subsequently formated and fit to various models using the phenor package.

Example code below shows that in few lines a modelling exercise can be set up.

```R
# process phenocam transition files into a consistent format
phenocam_data = process.phenocam("~/Dropbox/Data/tmp/transition_dates/")

# # downloads of PhenoCam data are facilitated using the functions
# # using the phenophase = TRUE flags, transition dates are estimated
# # on the go. The command below downloads all deciduous broadleaf sites
# # and estimates their respective phenophases.
# download.phenocam(vegetation = "DB", phenophase = TRUE)
```

# Model development, and parameter estimation

The gathered data can now be used in model validation. Out of the box 8 models are provided as described in the Eli's paper (Melaas et al. 2016, Global Change Biology paper). All models are listed in the phenology.models.r file and provided by loading the package. Your own model development can be done by creating similar functions which take in the data as described by the functions above.

Also consider that formats are consistent between both process.*() models and one can mix data streams, or append to them as one sees fit.

Below is the function to optimize the model parameters for the *SEQ1.3* model using Generalized Simulated Annealing (GenSA). Uppper and lower constraints to the parameter space have to be provided, and in case of GenSA initial parameters are estimated when par = NULL.

The sites used are those used in Eli's paper, there is a slightly higher number of site years due to longer time series used.


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

## Notes

All models included in the Melaas et al. 2013 study are included. As of now, these are untested. No guarantees are provided as to the quality of the functions.

### Dependencies

The code depends on the following R packages: devtools, GenSA
