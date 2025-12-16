#========== Shell for Disease Progression Model ================
#   Goal: Make an initial structure for disease progression model
# 	    (Total of 6 health states:  no-MCI (cognitively intact)
#				         MCI, 
#                                                      mild dementia, 
#                                                      moderate dementia, 
#                                                      severe dementia, and
#                                                      death.
#   	   Numbers here are arbitrary. 
#	   The purpose of this exercise is to build a structure, not to make things perfect.
#----------------------------------------------------------------------------------------------------------------------
#   Assumptions: (1) For now, the cycle is 1 year. We have 50 cycles, starting from age 50.
#   		          (Note that each patient can have more than 1 MMSE per year. 
#   		           For now, assume there is only 1 MMSE observation per year.)
#                         
#   		   (2) The transition probability from the cognitively intact status to a cognitively impaired status
#			is a logistic function of “age, education attainment, race/ethnicity, obesity”, 
#			which estimate log OR
#
# 		   (3) If one develops a cognitive impairment, then update (or, lower) MMSE using the equation 
#		         of “rate of MMSE change” from (Getsios, 2010): 
# 	                    “Rate of MMSE change  
#			             = f(previous MMSE, time since MCI onset, age at MCI onset)” 
#----------------------------------------------------------------------------------------------------------------------
#   Steps : (1) Import data 
#               (2) Generate new variables 
#   	     (3) Set the coefficient values and initialize parameters.
# 	     (4) Make matrices and vectors, and initialize their first elements
#	     (5) Run simulation loops
#               (6) Make a stacked area chart
#	     (7) Make graphs to analyze disease progression over time
#----------------------------------------------------------------------------------------------------------------------
#   Revisions : 2.Feb.2024: 	(1) age variable fix in MMSE equation
#				(2) sample size increase: 30 => 300
#				(3) random seed, instead of one that depends on i and t
#		8.Feb.2024:   Mortality change
#		12.Feb.2024: Added “age interactions” in the MCI transition function
#		21.Feb.2024: Relative Risk of Severe Dementia made Age-Specific
#===================================================================
#---------------------------------------------  Directory Prep ---------------------------------------------------
#-- Check the current directory 
getwd()
#-- Set the current directory to the location where my dataset is located
setwd("/Users/geenakim/Desktop/UCSF - Dementia/Men")
#----------------------------------------------------------------------------------------------------------------------
#------- Import hypothetical data to describe population
dementia_data <- read.csv("dementia_hypothetical_data - Data_initially_all_starts_at_50.csv")
#------- Import life table data to get the probability of dying next year by age
life_table_data <- read.csv("table_prob_die_next_year.csv", header = TRUE)
#------- Generate new variables 
#-- Initialize the variable to NA. To be replaced with actual value when MCI occurs
dementia_data$age_at_mci_onset <- NA				
#-- Ordinally encode the "edu" variable, i.e. generate a new variable with numbers instead of words. Nested “ifelse” approach. 
dementia_data$edu_num <- ifelse(dementia_data$edu == "college_more", 1, ifelse(dementia_data$edu == "high_ged", 2, 3))
#-- Generate a new variable "hisp_black" using race_eth, i.e. generate a new variable assigning 1 to hispanic or black and 0 otherwise.
dementia_data$hisp_black <- ifelse(dementia_data$race_eth %in% c("black", "hisp"), 1, 0)
#-- Generate a new variable "obesity_num", converting y/n to 1/0
dementia_data$obesity_num <- ifelse(dementia_data$obesity == "n", 0, 1)

#-------- Call the “plot_dementia_progress” function 
source("/Users/geenakim/Desktop/UCSF - Dementia/Men/plot_dementia_progress_24Mar2024.R")
plot_dementia_progress(-5.0, 0.001, -0.0001, 0.00005 , 0.001, 0.001, 0.001, -0.002, -0.002, -0.002)
