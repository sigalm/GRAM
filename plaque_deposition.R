
# Set directory
setwd("~/Documents/GitHub/GRAM")

# Load libraries
library(dplyr)
library(reshape2)
library(ggplot2)
source("amyloid_model.R")


model1 <- my_model(n_i=10,n_t=20)

model1$model_plot

model1$model_output


# increase n_i.
# increase sd - experiment with higher variability across time. 
# what's the test performance of pet scans compared to autopsy-determined pathology?

# Is it true that amyloid accumulates in this steady rate?
# Relative importance of amyloid (outside neurons) vs tau (inside neurons)
# Pivot toward protraying the tau process instead of amyloid.
# best way is to look at brain after death and to Braak staging.
# can do tau PET while alive, and compare to autopsy result Braak
# Tau PET is a good but imperfect window into the autopsy-measured Braak
# How good are CSF and plasma biomarkers against tau PET?

# Portray Braak histopathological development of disease
# Layer on tau PET to see that

# Prevalence / absence of pathology
# If present, how does it evolve?
# How do different tests perform at different stages of its evolution?
# Tau PET can miss early stages of tau deposition -- 
## plasma tau might do better for early detection
# treatments don't clear tau, it's not a measure that corresponds to treatment

# How does Braak translate to fundamental equations?
# a test like tau pet is not dichotomous -- test performance to allow for multiple stages of results
# How do you characterize test performance
# Interpretation of tau PET



# 1)	Divide individuals into have or don't have AD pathology
# 
# 2)	For individuals that do have AD path, portray evolution over time ... onset, some quant measure of (growing) severity ... check path literature eg concentration of amyloid plaques.
# 
# a.	PET tau is better than plaques (can use plasma). Tracks better with trajectory of disease and clinical severity. What you’re really trying to understand is the density of the plaques. They’ll stage patients after death. 
# 
# 3)	Portray test perf (sens / spec) of specific tests for specified reference, e.g., tau PET vs actual pathology, or plasma biomarker (specific one) vs. tau PET. Thus, as underlying path gets more severe, the tests are more likely to detect.
# 
# 4)	For specified prevalence of AD path, translate sens & spec to PPV and NPV.
# 
# 5)	MCI or dementia determination by which cognitive function assessment/s? See issue 6a below.
