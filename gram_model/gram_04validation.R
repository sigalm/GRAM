# gram_04validation.R
# Quick validation checks for GRAM simulation results
# Assumes a.out is available in the environment (output from f.run or similar)

# 1. Inspect people who clear memloss and have CDR progression
validate_memloss_clearance <- function(a.out) {
  n_cycles <- dim(a.out)[1]
  n_ind <- dim(a.out)[3]
  memloss <- a.out[, "MEMLOSS", ]
  cdr <- a.out[, "CDR", ]
  res <- data.frame(id = integer(), cycles_in_memloss = integer(), cdr_before = numeric(), cdr_after = numeric())

  for (i in 1:n_ind) {
    memloss_i <- memloss[, i]
    cdr_i <- cdr[, i]
    clear_cycle <- which(diff(memloss_i) == -1) + 1  # cycle where memloss goes from 1 to 0
    if (length(clear_cycle) > 0) {
      for (cycle in clear_cycle) {
        # Only consider if CDR increases after memloss clears
        cdr_before <- cdr_i[cycle - 1]
        cdr_after <- cdr_i[cycle:min(cycle + 2, n_cycles)]
        if (any(cdr_after > cdr_before, na.rm = TRUE)) {
          # Calculate how many cycles memloss was 1 before clearing
          start_memloss <- max(which(memloss_i[1:(cycle - 1)] == 0)) + 1
          if (is.finite(start_memloss)) {
            cycles_in_memloss <- (cycle - 1) - start_memloss + 1
            res <- rbind(res, data.frame(id = i, cycles_in_memloss = cycles_in_memloss, cdr_before = cdr_before, cdr_after = min(cdr_after[cdr_after > cdr_before], na.rm = TRUE)))
          }
        }
      }
    }
  }
  print("\nValidation 1: MEMLOSS duration before clearance and CDR progression:")
  if (nrow(res) > 0) {
    print(summary(res$cycles_in_memloss))
    print("If only 1, that's a problem. Should see >1 too.")
  } else {
    print("No subjects found who clear MEMLOSS and have CDR progression.")
  }
}


# 3. Anyone who becomes syn 0.5 stays for exactly two cycles, then becomes 1
validate_syn_half_duration <- function(a.out) {
  n_cycles <- dim(a.out)[1]
  n_ind <- dim(a.out)[3]
  syn <- a.out[, "SYN", ]
  bad_cases <- data.frame(id = integer(), cycle = integer(), duration = integer())
  for (i in 1:n_ind) {
    syn_i <- syn[, i]
    half_cycles <- which(syn_i == 0.5)
    if (length(half_cycles) > 0) {
      # Find runs of consecutive 0.5
      rle_half <- rle(syn_i == 0.5)
      idx <- which(rle_half$values)
      pos <- cumsum(rle_half$lengths)
      for (j in idx) {
        end_cycle <- pos[j]
        start_cycle <- end_cycle - rle_half$lengths[j] + 1
        duration <- end_cycle - start_cycle + 1
        # Check that after two cycles, syn becomes 1, or allow shorter duration if next value is NA (death)
        next_syn <- if (end_cycle < n_cycles) syn_i[end_cycle + 1] else NA
        # Only flag if duration != 2 AND the next value is not NA (death) or 1
        if ((duration != 2 && (is.na(next_syn) || next_syn != 1)) || (duration == 2 && !(is.na(next_syn) || next_syn == 1))) {
          bad_cases <- rbind(bad_cases, data.frame(id = i, cycle = start_cycle, duration = duration))
        }
      }
    }
  }
  print("\nValidation 3: SYN==0.5 duration check:")
  if (nrow(bad_cases) > 0) {
    print("Subjects with incorrect SYN==0.5 duration or transition:")
    print(bad_cases)
    # Print SYN trajectories for first 5 flagged cases
    n_debug <- min(5, nrow(bad_cases))
    if (n_debug > 0) {
      print("Debug: SYN trajectories for first flagged cases:")
      for (k in 1:n_debug) {
        id <- bad_cases$id[k]
        cycle <- bad_cases$cycle[k]
        syn_i <- a.out[, "SYN", id]
        print(paste0("ID ", id, " (cycle ", cycle, "): ", paste(round(syn_i, 2), collapse=",")))
      }
    }
  } else {
    print("PASS: All subjects with SYN==0.5 stay for 2 cycles, then become 1.")
  }
}

# Example usage (assuming a.out is loaded):
validate_memloss_clearance(sim_calib$sim$output)
validate_syn_half_duration(sim_calib$sim$output)
