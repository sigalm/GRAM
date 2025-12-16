######################################## GRAM HELPER FUNCTIONS: PLOTTERS ########################################

# generate figures
f.make_figures <- function(l.out, l.inputs) {
  
  # Reshape syndrome data to long format for use with ggplot.
  
  #### True (unobserved) status
  t.trace_syndrome_true <- as.data.frame(l.out$state_trace[ ,c("healthy", "mci", "mil", "mod", "sev", "dth")],) %>%
    mutate(year = (1:l.inputs[["n.cycle"]]) + l.inputs[["AGE_start_mean"]]-1,
           mod_sev = mod + sev) %>%
    select(-mod, -sev) 
  t.trace_syndrome_true_long <- t.trace_syndrome_true %>%
    pivot_longer(cols = -year, names_to = "syndrome", values_to = "proportion") %>%
    mutate(syndrome = fct_rev(factor(syndrome, levels = c("healthy", "mci", "mil", "mod_sev", "dth"))))
  
  fig.progression_true <- ggplot(t.trace_syndrome_true_long, aes(x = year, y = proportion, fill = syndrome)) +
    geom_area(alpha = 0.8, position = "stack") + 
    geom_path(aes(group = syndrome), position = "stack", color = "black", linewidth = 0.5) + 
    scale_x_continuous(breaks =  seq(min(t.trace_syndrome_true_long$year), max(t.trace_syndrome_true_long$year), by = 5), minor_breaks = NULL) +
    scale_fill_manual(values = c("white","purple","green","yellow","pink"),
                      labels = c("Death","Moderate to severe dementia","Mild dementia","MCI","Cognitively intact"),
                      guide = guide_legend(override.aes = list(colour = "black", size = 0.5))) +
    labs(title = "GRAM: Progression of Cognitive Impairment",
         subtitle = paste0(l.inputs[["scenario"]], "\n(N = ", l.inputs[["n.ind"]], " individuals)"),
         x = "Age",
         y = "Proportion of Population",
         fill = "True Cognitive Status") +
    theme(axis.text = element_text(size = 14),  # Increase axis tick font size
          axis.title = element_text(size = 16),
          panel.grid.major = element_line(color = "gray50", linewidth = 0.8),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          title = element_text(size = 16))
  
  return(fig.progression_true)
}

f.make_histogram <- function(values, lab) {
  values <- data.frame(values = values)
  h <- ggplot(data = values, aes(x = values)) +
    geom_histogram(aes(y = after_stat(count) / sum(after_stat(count)) * 100),
                   binwidth = 5) +
    labs(title = lab,
         x = lab,
         y = "Percentage")
  
  return(h)
}
