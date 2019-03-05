#' Model cross valiadation
#'
#' Faciliates model cross-validation using k-fold or leave-one-out cross
#' validation.
#'
#' @param k perform k-fold cross validation
#' (set k=0 to perform a leave-one-out cross-validation)
#' @param cv_seed a vector with random seeds for cross validation
#' (to shuffle the input data)
#' @param cv_set  (optional) use custom sets cross validation
#' (requires a list of vectors the indices of the data used for validation
#' (hold-out), the remaining data will be used for calibration)
#' @param ncores  number of cores to use for crossvalidation
#' (requires 'doParallel' and 'foreach'  library)
#' @param random_seeds a vector with random seeds for optimization replicates
#' @param models list of models to compare
#' @param data which standard or custom dataset to use
#' @param method which optimization method to use, GenSA or rgenoud
#' (default = GenSA)
#' @param control additional optimization control parameters
#' (default = list(max.call = 5000, temperature = 10000))
#' @param par_ranges location of the parameter ranges of the models
#' @keywords phenology, model, validation, comparison
#' @importFrom foreach %dopar%
#' @export
#' @examples
#'
#' \dontrun{
#'  pr_cross_validate(k=10,CV_seed = 123, models =  c("LIN","TT","PTT") ,
#'  data = phenocam_DB[c(1:10)])
#'
#' # will return the data of each validation and some summarizing statistics
#' }

pr_cross_validate <- function(
  k = 10,
  cv_seed = 1234,
  cv_set = NULL,
  ncores = 1,
  random_seeds = c(1, 12, 40),
  models = c(
    "LIN",
    "TT",
    "TTs",
    "PTT",
    "PTTs",
    "M1",
    "M1s",
    "AT",
    "SQ",
    "SQb",
    "SM1",
    "SM1b",
    "PA",
    "PAb",
    "PM1",
    "PM1b",
    "UN",
    "UM1",
    "SGSI",
    "AGSI"
  ),
  data = phenor::phenocam_DB,
  method = "GenSA",
  control = list(max.call = 5000, temperature = 10000),
  par_ranges = sprintf("%s/extdata/parameter_ranges.csv",
                       path.package("phenor"))
) {

  # get a specified subset out of flat data set
  get_cv_subset <- function(data,ids){
    CV_subset <- list(
      site = data$site[ids, drop = FALSE],
      location = data$location[, ids, drop = FALSE],
      doy = data$doy,
      transition_dates = data$transition_dates[ids, drop = FALSE],
      ltm = data$ltm[, ids, drop = FALSE],
      Ti = data$Ti[, ids, drop = FALSE],
      Tmini = data$Tmini[, ids, drop = FALSE],
      Tmaxi = data$Tmaxi[, ids, drop = FALSE],
      Li = data$Li[, ids, drop = FALSE],
      Pi = data$Pi[, ids, drop = FALSE],
      VPDi = data$VPDi[, ids, drop = FALSE],
      georeferencing = data$georeferencing[ids, drop = FALSE]
    )
    return (CV_subset)
  }

  # Train and validate a fold in k-fold (or LOO) CV
  fold <- function(i){
    print(sprintf(" %s/%s",i,k))
    val_data<-get_cv_subset(data,cv_set[[i]])
    train_data<-get_cv_subset(data,c(1:n) [-cv_set[[i]]])

    # Estimate Parameters on training dataset
    CVfold<-model_comparison(
      random_seeds,
      models,
      train_data,
      method,
      control,
      par_ranges)
    CVfold$validation_measured<-val_data$transition_dates

    # estimate model output for validation data using the optimized parameters
    for (j in 1:length(models)){
      vp<-matrix(NA,ncol = length(val_data$transition_dates),
                 nrow=length(random_seeds))
      for (l in 1:length(random_seeds)){
        out<-estimate_phenology(data = val_data,model = models[j],
                                par = CVfold$modelled[[j]]$parameters[l,])
        vp[l,] <- out
      }
      cv_fold$modelled[[j]]$validation_predicted_values<-vp
    }
    return(cv_fold)
  }

  # convert to a flat format for speed
  # The call to model_copmarison call further down will call this again;
  # but already flat datasets are ignored by the function
  data <- flat_format(data)
  cvindex <- c(1:length(data$transition_date))

  # Prepare training and validation dataset indices ...
  n <- length(cvindex)

  if (is.null(cv_set)) {

    # Create the CV sets here
    if (k == 0) {
      k <- n
    }

    vc <- (n - n %% k) / k

    if (n %% k > 0) {
      print (
        sprintf(
          "%s-fold CV: Warning: %s random value(s) dropped while creating %s sets of equal size (%s values each)",
          k,
          n %% k,
          k,
          vc
        )
      )
    }
    if (k == n) {
      # LEAVE-ONE-OUT
      rnd_set <- c(1:n)
      print ("LEAVE-ONE_OUT CV")
    } else {
      # K-FOLD
      set.seed(cv_seed)
      rnd_set <- sample(cvindex)
      cat(sprintf("%s-FOLD CV \n", k))
    }
    cv <- rnd_set[1:(n - n %% k)]
    cv_set = list()
    for (i in 1:k) {
      vci <- c(1:vc) + vc * (i - 1)
      cv_set[[i]] <- cv[vci]
    }
  } else {
    #..or use the ones provided by the user

    if (length(cv_set) != k) {
      stop("Invalid parameter cv_set provided")
    }
  }

  # Perform CV
  if (ncores == 1){ # Synchronous mode

    results <- list()
    for (i in 1:k){
      results[[i]] = fold(i)
      }

  } else { # Parallel mode

    cat(sprintf("Using %s cores to calculate the cross validation (%s sets).\n
                  This may still take a while...\n",ncores,k))

    # Use n cores for parallel CV
    cl<- parallel::makeCluster(ncores)

    doParallel::registerDoParallel()
    on.exit(parallel::stopCluster(cl))
    results <- foreach::foreach(
      i=1:k,
      .export=c("model_comparison",
                "estimate_phenology")) %dopar% {fold(i)}
  }
  names(results)<-c(1:k)


  #Caluclate some basic stats
  cv_stat=list()
  for (m in 1:length(models)) {

    tmp_trn <- matrix(NA, ncol = length(random_seeds), nrow = k)
    tmp_val <- matrix(NA, ncol = length(random_seeds), nrow = k)

    for (i in 1:k) {
      for (j in 1:length(random_seeds)) {
        out <- results[[i]]$modelled[[m]]$predicted_values[j, ]
        val <- results[[i]]$measured
        tmp_trn[i, j] = sqrt(mean((val - out) ^ 2, na.rm = T))
        out <- results[[i]]$modelled[[m]]$validation_predicted_values[j, ]
        val <- results[[i]]$validation_measured
        tmp_val[i, j] = sqrt(mean((val - out) ^ 2, na.rm = T))
      }
    }
    mrmse_trn <- mean(tmp_trn)
    mrmse_val <- mean(tmp_val)
    cv_stat[[m]] <- list(training_errors = tmp_trn,
                         validation_errors = tmp_val,
                         mean_training_error = mrmse_trn,
                         mean_validation_error = mrmse_val)
  }

  names(cv_stat) <- models

  # combine and return
  results <- list(results, cv_stat)
  names(results) <- c("cv_predictions", "cv_statistics")
  return (results)
}

