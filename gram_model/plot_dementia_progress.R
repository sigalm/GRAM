#=================== Modeling Disease Progression: dementia_disease_progression_plot.R  ======================
#   Function of plotting dementia_disease_progression as a function of coefficients of MCI transition function
#=============================================================================================================

##-------------  Define a function to plot dementia progression —------------##
plot_dementia_progress  <- function(b0, b1, b11, b12, b2, b3, b4, b1edu, b1hisp, b1ob) {
  #- Ultimately, we will calibrate these coefficients.
  # b0: 		intercept 
  # b1: 		coefficient for “age-50” 
  # b11: 		coefficient for (age -50)^2
  # b12: 		coefficient for (age -50)^3
  # b2: 		coefficient for low education
  # b3: 		coefficient for hispanic or black
  # b4: 		coefficient for obesity status
  # b1edu: 	coefficient for education * (age -50)
  # b1hisp: coefficient for Hispanic/Black * (age -50)
  # b1ob: 	coefficient for obesity *(age -50)
  
  # ===============SET UP LOGS/RECORDS ==================
  start_time <- Sys.time()
  
  timestamp <- format(start_time, "%Y%m%d_%H%M")
  log_file_name <- paste0("log_", timestamp, ".txt")
  log_file <- file(log_file_name, open = "a")
  
  log_output <- function(level, message, file) {
    cat(paste(Sys.time(), paste0("[", level, "]"), message, "\n"), file = file, append = TRUE)
  }
  
  #------- Initialize model parameters
  n_i <- 3000   			 # number of individuals. 
  n_t <- 50			 # number of maximum cycles per individual
  
  log_output(100, paste0("Starting simulation for n_i = ", n_i, 
                         " and n_t = ", n_t, "."), log_file)
  
  #------- Make matrices of “n_i by n_t + 1” to keep track of variable values by cycle for each individual
  #-- The rows are labeled as “indiv 1”, “indiv 2”, …, “indiv n_i”
  #-- The columns are labeled as “cycle 0”, “cycle 1”, …, “cycle n_t”
  m_age <-   
    m_alive <- 						
    m_mci_flag <- 
    m_mmse <- 
    m_delta_mmse <- 
    m_time_since_mci_onset <- matrix(nrow = n_i, ncol = n_t + 1, dimnames = list(paste("indiv", 1:n_i), paste("cycle", 0:n_t)))
  m_seed 		      <- matrix(nrow = n_i, ncol = n_t, dimnames = list(paste("indiv", 1:n_i), paste("cycle", 1:n_t)))
  m_p_dead_next_yr_all_cause <- as.matrix(life_table_data) # this matrix(age, qx) contains the probability of dying next year for each age
  
  #------- Initialize the matrix elements for cycle 0 of each individual i.
  m_age[ ,1] <- 50     				  # in this initial setup, all enter simulation at age 50. 
  m_alive[ ,1] <- 1				      # all enter simulation alive. 
  m_mci_flag[ ,1] <- 0				  # in this initial setup, none has any cognitive impairment when entering simulation.
  m_mmse[ ,1] <- 30				      # in this initial setup, all has MMSE = 30 when entering simulation.
  m_time_since_mci_onset[ ,1]  <- 0                      # in this initial setup, time since MCI onset is 0  when entering simulation.
  
  #------- prep for setting random seed number
  set.seed(20240202)					  # set seed using a fixed seed number to generate random number from a uniform distribution
  num_of_seed_val  <- n_i * n_t			  # number of seed values
  v_random_num <- runif(num_of_seed_val)*10^7	  # make a vector populated with random numbers
  #-- store the random seeds in a matrix 
  m_rand_unif <- matrix(v_random_num, nrow = n_i, ncol = n_t, dimnames = list(paste("indiv", 1:n_i), paste("cycle", 1:n_t)))
  
  # ================= SIMULATION LOOP: START ===================#
  
  #------- Initialization of vectors: Make vectors with “n_t+1” elements to trace the number of people by health state category, then zero out the values.
  v_total_pop_alive  <-
    v_total_pop_dead  <-
    v_pop_cog_imp_none <-
    v_pop_cog_imp_mci <-  
    v_pop_cog_imp_mild_dementia <-  
    v_pop_cog_imp_mod_dementia <-  
    v_pop_cog_imp_severe_dementia  <- rep(0, n_t+1)    
  
  #----- 1st LOOP —-# Calculate the probability of dying next year for those cognitively intact at cycle t=1.
  t <- 1
    for (i in 1:n_i) {	
      #-- Calculate the probability of dying next year for those cognitively intact at cycle t.
      #- Note that the CDC data for the probability of dying next year does not account for the severity of dementia. 
      #- Thus, we adjust the probability of dying for the severity of dementia based on the Relative Risks (RRs) with the Cognitively Intact being its reference group. The RRs are from Anderson 2010 https://pubmed.ncbi.nlm.nih.gov/20110702/ 
      #- We take the following steps (1), (2), (3), and (4) to calculate the intended probability. 
      
      if (m_alive[i,t] == 1) { 
        #-- (1) Count the number of people in each cognitive impairment state (none, mci, mild dementia, moderate dementia, severe dementia) in cycle t. (The current criterion for categorization may change.). We need this count to calculate the probability of dying next year by disease state. 
          #- mmse = 30: 	No cognitive impairment
          #- mmse = 26~29: 	MCI 
          #- mmse = 21~25: 	Mild dementia
          #- mmse = 11-20: 	Moderate dementia
          #- mmse = 0-10: 	Severe dementia
        v_pop_cog_imp_none[t]  <-  v_pop_cog_imp_none[t] + (m_mmse[i, t] == 30)
        v_pop_cog_imp_mci[t]   <-  v_pop_cog_imp_mci[t] + (m_mmse[i, t] >= 26 & m_mmse[i, t] < 30)
        v_pop_cog_imp_mild_dementia[t] <- v_pop_cog_imp_mild_dementia[t] +  (m_mmse[i, t] >= 21 & m_mmse[i, t] < 26)
        v_pop_cog_imp_mod_dementia[t]  <- v_pop_cog_imp_mod_dementia[t] +  (m_mmse[i, t] >= 11 & m_mmse[i, t] < 21)
        v_pop_cog_imp_severe_dementia[t] <- v_pop_cog_imp_severe_dementia[t] + (m_mmse[i, t] >= 0 & m_mmse[i, t] < 11)
      } else {m_alive[i,t+1] <- 0                                                          # if dead, skip the “i” loop after defining mortality status at t+1
        }	
    }
    #-- (2) Count the number of people alive/dead in cycle t
    v_total_pop_alive[t] <- v_pop_cog_imp_none[t] + v_pop_cog_imp_mci[t] + v_pop_cog_imp_mild_dementia[t] + v_pop_cog_imp_mod_dementia[t]+ v_pop_cog_imp_severe_dementia[t] 	
    v_total_pop_dead[t]  <- n_i - sum(m_alive[, t]) 
    
  #------- 2nd LOOP —------#
for (t in 1:(n_t-1)) {                    # start loop for time cycles: “t in 1:n_t” means “from t=1 to t=n_t”

  #-- Record t loop
  log_output(3, sprintf("Simulating cycle t=%s/%s", t, n_t), log_file) 
  
      for (i in 1:n_i) {						      # start loop for each individual
      
     # set.seed(as.integer(round(m_rand_unif[i, t])))
      set.seed(round(m_rand_unif[i,t]))	      # set random seed 
      m_age[i,t+1] <- m_age[i,t]+1			      # update age 
      
      # m_alive[i, t + 1] <- 1				        # just to test what happens when all live through 100
      
      if (m_alive[i,t] == 1) { 
        
        #-- (3) Then, we calculate the probability of the reference group (with the cognitively intact people) dying next year. 
        
        #-- Numerator: v_total_pop_alive[t] * all cause mortality m_p_dead_next_yr_all_cause[age_index, 2]
        age_index <- m_age[i,t+1] - 50 
        numerator <- v_total_pop_alive[t] * m_p_dead_next_yr_all_cause[age_index, 2]   # same for all individual i
        
        #-- Denominator: depends on the Relative Risk (RR). 
        ## RR for MCI: 			          1.82,
        ##	      mild dementia:		  2.92
        ##	      moderate dementia:	3.85
        ##	      severe dementia:		9.52	 
        # denom <-  	v_pop_cog_imp_none[t] + v_pop_cog_imp_mci[t]*1.82 + v_pop_cog_imp_mild_dementia[t]*2.92 + v_pop_cog_imp_mod_dementia[t]*3.85 + v_pop_cog_imp_severe_dementia[t]*9.52    		# this is without age-adjustment => commented out
        #-- (3-1). Here, we make the RR for severe dementia “age-specific.”
        # The Relative Risks (RRs) from the Anderson paper are not age-specific. Considering that mortality is generally very low at an early age, it is likely that individuals with severe dementia will have a higher mortality than the average for that age group. However, at a later age, when mortality is generally higher, the mortality of those with severe dementia may not be as high as at an earlier age. 
        ## First, try RR of 	15 for those aged 50-64
        ## 	        		      9  for those aged 65-79
        ##			              2  for ages 80+
        
        rr_sev_dement <- 	15*(m_age[i,t]>=50 & m_age[i,t]<65) + 9*(m_age[i,t]>=65 & m_age[i,t]<80) + 2*(m_age[i,t]>=80)
        
        #- denom differ by individual i
        denom <-  		v_pop_cog_imp_none[t] + v_pop_cog_imp_mci[t]*1.82 + v_pop_cog_imp_mild_dementia[t]*2.92 + v_pop_cog_imp_mod_dementia[t]*3.85 + v_pop_cog_imp_severe_dementia[t] *rr_sev_dement  
  
        if (any(is.na(denom))) {
          print("denom contains NA")
          print(paste("i:", i, "t:", t))
        } 
              
        #-- (4) Calculate the mortality for the cognitively intact then the mortality by severity of cognitive impairment using inputs calculated above. 
        #-  differ by individual i
        m_p_dead_next_yr_cog_intact  		  <- numerator / denom   	
        m_p_dead_next_yr_mci  		        <- min(1,m_p_dead_next_yr_cog_intact*1.82) 
        m_p_dead_next_yr_mild_dementia    <- min(1,m_p_dead_next_yr_cog_intact*2.92) 
        m_p_dead_next_yr_mod_dementia     <- min(1,m_p_dead_next_yr_cog_intact*3.85) 
        #m_p_dead_next_yr_severe_dementia <- min(1,m_p_dead_next_yr_cog_intact*9.52) 
        # this is commented out due to age-adjustment for the mortality of the severely demented
        m_p_dead_next_yr_severe_dementia  <- min(1,m_p_dead_next_yr_cog_intact*rr_sev_dement) 
        
        #- Determine whether this person is alive or dead at t+1 based on the Bernoulli distribution - a special case of binomial - using the probability of dying next year calculated above.  
        
        if (m_mmse[i, t] == 30) {
          m_alive[i, t+1] <- rbinom(n = 1, size = 1, prob = (1 - m_p_dead_next_yr_cog_intact))
        } else if (m_mmse[i, t] >= 26 & m_mmse[i, t] < 30) {
          m_alive[i, t+1] <- rbinom(n = 1, size = 1, prob = (1 - m_p_dead_next_yr_mci))
        } else if (m_mmse[i, t] >= 21 & m_mmse[i, t] < 26) {
          m_alive[i, t+1] <- rbinom(n = 1, size = 1, prob = (1 - m_p_dead_next_yr_mild_dementia))
        } else if (m_mmse[i, t] >= 11 & m_mmse[i, t] < 21) {
          m_alive[i, t+1] <- rbinom(n = 1, size = 1, prob = (1 - m_p_dead_next_yr_mod_dementia))
        } else if (m_mmse[i, t] >= 0 & m_mmse[i, t] < 11) {
          m_alive[i, t+1] <- rbinom(n = 1, size = 1, prob = (1 - m_p_dead_next_yr_severe_dementia))
        }
      } else {m_alive[i,t+1] <- 0                                                          # if dead, skip the “i” loop after defining mortality status at t+1
              next}	
  
      #-- Calculate the MCI transition probability for individual i in cycle t
      x <- b0 + b1*(m_age[i,t]-50) + b2*dementia_data$edu_num[i] + b3*dementia_data$hisp_black[i] + b4*dementia_data$obesity_num[i] + b11*(m_age[i,t]-50)^2 + b12*(m_age[i,t] -50)^3 + b1edu * (m_age[i,t] -50) * dementia_data$edu_num[i]  + b1hisp* (m_age[i,t] -50) * dementia_data$hisp_black[i] + b1ob * (m_age[i,t] -50) * dementia_data$obesity_num[i] 
      m_p_none_to_mci <- exp(x) / (1 + exp(x)) 		 	# logistic function
      
      #-- Assign next cycle’s flag for being cognitively impaired (MCI or worse) 
      #- if this MCI flag at t is not 1, then assign MCI flag at t+1 based on the Bernoulli distribution - a special case of binomial - using the calculated transition probability
      #- If the MCI flag at t is 1, then the MCI flag at t+1 is always 1. (MCI does not get better.) 
      
      # if (m_mci_flag[i, t] != 1) {
      #   m_mci_flag[i,t+1] <- rbinom(n = 1, size = 1, prob = m_p_none_to_mci)
      # } else {m_mci_flag[i,t+1] <- 1}	
      
      # Set the MCI_flag value at t+1, depending on the MCI_flag value at t 
      if (m_mci_flag[i, t] == 1) {
        m_mci_flag[i, (t+1):ncol(m_mci_flag)] <- 1  # Set t+1, t+2, etc. to 1 once m_mci_flag[i, t] is 1
      } else {       # if (m_mci_flag[i, t] != 1)
        m_mci_flag[i,t+1] <- rbinom(n = 1, size = 1, prob = m_p_none_to_mci)
      }
    

#-- Define “age_at_mci_onset[i]” when the MCI flag turns on at t+1
if (m_mci_flag[i, t] != 1 & m_mci_flag[i, t + 1] == 1) {
 	dementia_data$age_at_mci_onset[i] <-  m_age[i, t + 1]
}
#-- If the MCI flag at t+1 = 1 then update (lower) MMSE using the equation for rate of MMSE change (Getsios, 2010)
#- applying the rate of change equation can yield negative MMSE value, in which case, MMSE is set to 0.
#- if the MCI flag at t+1 is not 1, then MMSE is set to 30.

if (m_mci_flag[i, t + 1] == 1) {
m_time_since_mci_onset[i, t + 1] <-  m_age[i, t + 1] - dementia_data$age_at_mci_onset[i] + 1
m_delta_mmse[i, t + 1] <- -5.4663 - 0.4299 * min(m_mmse[i, t], 9) - 0.0042 * max(0, min(m_mmse[i, t] - 9, 9)) + 0.1415 * max(0, min(m_mmse[i, t] - 18, 12)) - 0.0791 * (m_mmse[i, 1] - 30) / 1 + 0.0747 * 80		
# no one entered the data with MCI, so the denominator of the term “-0.0791 * (m_mmse[i, 1] - 30) / 1” is 1
# m_delta_mmse[i, t + 1] <- -5.4663 - 0.4299 * min(m_mmse[i, t], 9) - 0.0042 * max(0, min(m_mmse[i, t] - 9, 9)) + 0.1415 * max(0, min(m_mmse[i, t] - 18, 12)) - 0.0791 * (m_mmse[i, 1] - 30) / 1 + 0.0747 * dementia_data$age_at_mci_onset[i]    
m_mmse[i, t + 1] <- round(min(30, max(0, m_delta_mmse[i, t + 1] + m_mmse[i, t])))
} else {	m_mmse[i, t + 1] <- 30
m_time_since_mci_onset[i,t + 1]  <- 0  
}

      
      #-- Update vectors of population by severity of cognitive impairment for cycle t+1
        v_pop_cog_imp_none[t+1]  <-  v_pop_cog_imp_none[t+1] + (m_mmse[i,t+1] == 30)
        v_pop_cog_imp_mci[t+1]   <-  v_pop_cog_imp_mci[t+1] + (m_mmse[i,t+1] >= 26 & m_mmse[i,t+1] < 30)
        v_pop_cog_imp_mild_dementia[t+1] <- v_pop_cog_imp_mild_dementia[t+1] + (m_mmse[i,t+1] >= 21 & m_mmse[i,t+1] < 26)
        v_pop_cog_imp_mod_dementia[t+1]  <- v_pop_cog_imp_mod_dementia[t+1] + (m_mmse[i,t+1] >= 11 & m_mmse[i,t+1] < 21)
        v_pop_cog_imp_severe_dementia[t+1] <- v_pop_cog_imp_severe_dementia[t+1] + (m_mmse[i,t+1] >= 0 & m_mmse[i,t+1] < 11)
 
      #-- Count the number of people alive/dead in cycle t+1
        v_total_pop_alive[t+1] <- v_pop_cog_imp_none[t+1] + v_pop_cog_imp_mci[t+1] + v_pop_cog_imp_mild_dementia[t+1] + v_pop_cog_imp_mod_dementia[t+1] + v_pop_cog_imp_severe_dementia[t+1] 	
        v_total_pop_dead[t+1]  <- n_i - v_total_pop_alive[t+1]
        
      # -- Return the modified dementia_data invisibly. Otherwise, no update is made in the dementia_data.
      invisible({
        assign("dementia_data", dementia_data, envir = .GlobalEnv)
        assign("m_alive", m_alive, envir = .GlobalEnv)
        assign("m_mci_flag", m_mci_flag, envir = .GlobalEnv)
        assign("m_mmse", m_mmse, envir = .GlobalEnv)
        assign("m_delta_mmse", m_delta_mmse, envir = .GlobalEnv)
        assign("m_time_since_mci_onset", m_time_since_mci_onset, envir = .GlobalEnv)
      })
   
    }  # close loop for individuals (i)		
  }  # close loop for cycles (t)
  
 # v_total_pop_dead[n_t+1] <- n_i           # All die by age 100   
  #v_total_pop_dead[n_t+1] <- 0						  # Just to test how people would behave when all live by age 100   
  
  # ================= SIMULATION LOOP: END ===================#
  log_output(100, "Simulation complete. Preparing results.", log_file)  
    
  #####-----Plot a stacked area chart of cognitive impairment progression ----#####
  
  #-- make a trace data for cognitive impairment states and the death state 
  library(ggplot2)
  library(tidyr)
  trace <- data.frame(   age = 50:100,
                         cog_intact 		= v_pop_cog_imp_none,
                         MCI 			      = v_pop_cog_imp_mci,
                         mild_dementia 	= v_pop_cog_imp_mild_dementia,
                         mod_dementia 	= v_pop_cog_imp_mod_dementia,
                         severe_dementia 	= v_pop_cog_imp_severe_dementia,
                         death 			    = v_total_pop_dead
  )
  
  #-- reshape the trace data from wide to long 
  trace_long 		    <- pivot_longer(trace, -age, names_to = "state", values_to = "value")
  desired_order 		<- c("death", "severe_dementia", "mod_dementia", "mild_dementia", "MCI", "cog_intact")
  trace_long$state 	<- factor(trace_long$state, levels = desired_order)
  
  #-- plot a stacked area chart
  # Define color palette
  my_colors <- c("cog_intact" 			= "pink",
                 "MCI" 			        = "yellow",
                 "mild_dementia" 		= "green",
                 "mod_dementia" 		= "purple",
                 "severe_dementia" 	= "orange",
                 "death" 			      = "white"
  )
  
  plot_dementia_progression <- ggplot(trace_long, aes(x = age, y = value, group=state, fill = state, order = dplyr::desc(state))) +
    geom_area(alpha = .7, color = "black", linewidth = 0.3) +
    geom_line(position = "stack", linewidth = .4) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +  # Add horizontal line at y = 0
    geom_vline(xintercept = c(50, 100), color = "black", linewidth = 0.3) +  
      scale_fill_manual(values = my_colors) +
      labs(title = "Cognitive Impairment Progression",
         x = "Age",
         y = "Number of People") +
    theme(	plot.title = element_text(size = 10, face = "bold"),
           legend.position = "bottom", 
           legend.key.size = unit(0.4, "cm"),
           axis.text = element_text(size = 8),    # Adjust x- and y-axis text size
           axis.title = element_text(size = 10),  # Adjust x- and y-axis title size
    ) +
    guides(fill = guide_legend(reverse = TRUE))
 
# =============== Record Total run time & Close the log file ==================
end_time <- Sys.time()
log_output(100, paste0("Results saved. Total run time: ", 
                    difftime(end_time, start_time)), log_file)

  if (!is.null(log_file)) {
      close(log_file)
      }
  
return(plot_dementia_progression)
}
#############################################################################

