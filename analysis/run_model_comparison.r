# Run model validation comparison for the different datasets
# Save the resulting data as a rda file!
library(phenor)

# 12 random seeds (hard coded, but can be generated using code below)
#seeds = sample(1:1000, 12)
seeds = c(204,680,364,350,
          62,481,397,17,
          125,395,500,325)

# GenSA control parameters
con = list(max.call = 40000,
           temperature = 20000)

# phenor datasets to evaluate
datasets = c("phenocam_DB",
             "phenocam_EN",
             "phenocam_GR")

# create output matrix
comparison = list()

# run comparison for the standard data sets
for (i in 1:length(datasets)){

  # provide feedback
  cat(sprintf("processing: %s\n",datasets[i]))

  # subset data and optimize
  data(list = datasets[i])
  comparison[i] = list(model_comparison(dataset = datasets[i],
                                        random_seeds = seeds,
                                        control = con))
}

# give the list items proper names
names(comparison) = datasets

# save the data
saveRDS(comparison,"~/comparison.rds")

# run a similar analysis for the sites
# as mentioned in Melaas et al. 2016
# site are listed below
sites_melaas = c(
  "harvard",
  "bartlettir",
  "acadia",
  "mammothcave",
  "nationalcapital",
  "dollysods",
  "smokylook",
  "upperbuffalo",
  "boundarywaters",
  "groundhog",
  "umichbiological2",
  "queens"
)

# create subset based upon the sites listed above
phenocam_DB_subset = phenocam_DB[sites_melaas]

# run the comparison
melaas = model_comparison(dataset = phenocam_DB_subset,
                          random_seeds = seeds,
                          control = con)

# save the data
saveRDS(melaas,"~/melaas.rds")

