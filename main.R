## ---------------------------
##
## Script name: main.R
##
## Purpose of script: Running infection forecasting model
##
## Author: Mr. Ashleigh C. Myall
##
## Date Created: 04-09-2021
##
## Copyright (c) Ashleigh C. Myall 2020
## Email: a.myall19@imperial.ac.uk
##
## ---------------------------
##
## Notes: This repo provides an implementable example of the model proposed in 
##        Myall et al. 2021 [link]. R/ is the folder for scripts that contain R 
##        functions. All functions are documented with roxygen2 syntax.
##
## ---------------------------
##
## Todo: -add paper link
##
## ---------------------------

################################################################################

### Load Libs

################################################################################

library(dplyr)
library(tidyr)
library(igraph)
library(tidyverse)
library(readr)
library(caret)

################################################################################

###
###                               Functions
###

################################################################################

source("./R/functions.R")


################################################################################

##########
##########                       Main Work Flow 
##########

################################################################################



# ------------------------------------------------------------------------------

##
##  0. Example data
##

# The framework is flexible and can integrate a wide range of variables. However, 
# the main functionality relies on patient pathways, with test results (date and 
# result) that detail when and where someone was identified as infected. We provide 
# examples of both data which can be called using getExamplePathways() and 
# getExampleTesting().

#  Read in example patient pathways
pathways = read.csv("data/examplePathways.csv")

# Read in example patient testing results
tests = read.csv("data/exampleTests.csv")

# Join pathways with test results by patient identifier
pathwaysWithTests = left_join(pathways,tests, by = c("Ptnumber"))

# -------

# Pathways are stored in a long format, specifying a patient identifier in pathways$Ptnumber, 
# their location and time at that location in pathways$location, and pathways$t 
# respectively, and then the date of their first positive test result in pathways$posTestResDt.

# In addition, data in data/staticVars.csv can contain fixed attributes of a patient 
# used for prediction, and data in data/contextualVars.csv includes background 
# statistics that change over time, which can also be integrated into the prediction 
# framework.

# Read in example patient static variables
staticVars = read.csv("data/exampleStaticVars.csv")

# Read in example contextual variables
contextualVars = read.csv("data/exampleContextualVars.csv")




# ------------------------------------------------------------------------------

##
##  1. Pre-process
##

# The function preProRollingWind() is a wrapper function that is applied over a 
# sliding time window of length feature_n. In summary, the function: (i) splits 
# the data into windows, (ii) constructs a contact network, and then centrality 
# of each patient across it, (iii) derives the background contextual variable for 
# the window, and (iii) joins the patient statistics within the window with static 
# variables.



# Run pre-proccesing over a sliding time window
preProRollingWind(pathwaysWithTests,  # Patient pathways with tests (when a patient become positive)
                  staticVars,         # Static variables (i.e. age, gender, ...)
                  contextualVars,     # Background contextual variables (hospital infection numbers)
                  feature_n = 14,     # Time window size to compute variables over
                  prediction_n = 7)   # Forecasting horizon

# -------

# Pre-processing data over extended periods or with many individuals can increase 
# computational expenditure. We advise the user to experiment with the parallel 
# implementation of the pre-processing function, preProRollingWindPar().

# Run pre-proccesing in parallel
#preProRollingWindPar(..., cores = 10)



# ------------------------------------------------------------------------------

##
##  2. Clean data
##

# Data pre-processing files are read in, aggregated, and cleaned, to produce final 
# data suitable for statistical analysis.

# Loads and aggregates saved files (fixed and static variables are read in and joined by patient-timecode IDs)
stat.df = loadPreData()

# Clean (add infection labels)
stat.df.clean = cleanStatData(stat.df)

# Prepare modelling dataset (Scale data, remove redundency, under-sample, split into train/test)
trainTestData.l = prepModData(stat.df.clean)


# ------------------------------------------------------------------------------

##
##  3. Statistical analysis
##

# The cleaned datasets are finally analysed following the paper by (i) performing 
# a univariate analysis over variables grouped and averaged across patients and 
# (ii) fitting a model to predict disease. The prediction model is implemented 
# using caret and thus can be flexibly changed and compared.

# Univariate analysis
uniVarAnalysis(trainTestData.l)

# Run prediction model
caret::train(Infection ~., data = trainTestData.l$train,method='rf')




