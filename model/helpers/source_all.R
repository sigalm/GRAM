# Parameter blocks that are deliberately not in model/setup.R: setup.R describes
# natural history, these describe what is done to it.
source("model/test_properties.R")
source("model/intervention_properties.R")

source("model/helpers/run_wrappers.R")
source("model/helpers/output_formatters.R")
lapply(list.files("model/helpers", pattern = "_helpers\\.R$", full.names = TRUE), source) 
lapply(list.files("model/modules", pattern = "^module_.*\\.R$", full.names = TRUE), source)
lapply(list.files("model/config", pattern = "_config\\.R$", full.names = TRUE), source)
