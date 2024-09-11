# Load required libraries
library(dplyr)
library(tidyverse)

simulate_population <- function(n) {
  data.frame(
    id = 1:n,
    age = sample(18:90, n, replace = TRUE),
    chronic_disease = sample(c(0, 1), n, replace = TRUE),
    race_ethnicity = sample(c("White", "Black", "Hispanic", "Other"), n, replace = TRUE),
    area_residence = sample(c("Urban", "Rural"), n, replace = TRUE)
  )
}

simulate_engagement <- function(pop, intervention, main_effect_size, hte_effect_size, noise_scale = 0.01) {
  base_prob <- 0.05
  
  intervention_effect <- numeric(nrow(pop))
  hte <- numeric(nrow(pop))
  
  empathetic_mask <- intervention == "empathetic"
  factual_mask <- intervention == "factual"
  weekday_mask <- intervention == "weekday"
  weekend_mask <- intervention == "weekend"
  human_mask <- intervention == "human"
  llm_mask <- intervention == "LLM"
  
  intervention_effect[empathetic_mask | factual_mask] <- main_effect_size
  hte[empathetic_mask] <- ifelse(pop$age[empathetic_mask] > 50, hte_effect_size, 0)
  hte[factual_mask] <- ifelse(pop$age[factual_mask] <= 50, hte_effect_size, 0)
  
  intervention_effect[weekday_mask | weekend_mask] <- main_effect_size
  hte[weekday_mask] <- ifelse(pop$race_ethnicity[weekday_mask] %in% c("Black", "Hispanic"), hte_effect_size, 0)
  hte[weekend_mask] <- ifelse(!(pop$race_ethnicity[weekend_mask] %in% c("Black", "Hispanic")), hte_effect_size, 0)
  
  intervention_effect[human_mask | llm_mask] <- main_effect_size
  hte[human_mask] <- ifelse(pop$chronic_disease[human_mask] == 1, hte_effect_size, 0)
  hte[llm_mask] <- ifelse(pop$chronic_disease[llm_mask] == 0, hte_effect_size, 0)
  
  noise <- rnorm(nrow(pop), mean = 0, sd = noise_scale)
  
  prob <- pmin(pmax(base_prob + intervention_effect + hte + noise, 0), 1)
  engagement <- rbinom(nrow(pop), 1, prob)
  
  pop$engagement <- engagement
  return(pop)
}

simulate_ab_tests <- function(n, main_effect_size, hte_effect_size) {
  pop <- simulate_population(round(n/3))
  
  pop_test1 <- pop
  pop_test1$intervention <- sample(c("empathetic", "factual"), round(n/3), replace = TRUE)
  pop_test1 <- simulate_engagement(pop_test1, pop_test1$intervention, main_effect_size, hte_effect_size)
  
  pop_test2 <- pop
  pop_test2$intervention <- sample(c("weekday", "weekend"), round(n/3), replace = TRUE)
  pop_test2 <- simulate_engagement(pop_test2, pop_test2$intervention, main_effect_size, hte_effect_size)
  
  pop_test3 <- pop
  pop_test3$intervention <- sample(c("LLM", "human"), round(n/3), replace = TRUE)
  pop_test3 <- simulate_engagement(pop_test3, pop_test3$intervention, main_effect_size, hte_effect_size)
  
  return(list(pop_test1, pop_test2, pop_test3))
}

simulate_smart <- function(n, main_effect_size, hte_effect_size) {
  pop <- simulate_population(n)
  
  pop$intervention1 <- factor(sample(c("empathetic", "factual"), n, replace = TRUE))
  pop <- simulate_engagement(pop, as.character(pop$intervention1), main_effect_size, hte_effect_size)
  
  non_responders <- pop[pop$engagement == 0, ]
  if (nrow(non_responders) > 0) {
    non_responders$intervention2 <- factor(sample(c("weekday", "weekend"), nrow(non_responders), replace = TRUE))
    non_responders <- simulate_engagement(non_responders, as.character(non_responders$intervention2), main_effect_size, hte_effect_size)
    pop[pop$engagement == 0, names(non_responders)] <- non_responders
  } else {
    pop$intervention2 <- factor(c("weekday", "weekend")[1])
  }
  
  non_responders <- pop[pop$engagement == 0, ]
  if (nrow(non_responders) > 0) {
    non_responders$intervention3 <- factor(sample(c("LLM", "human"), nrow(non_responders), replace = TRUE))
    non_responders <- simulate_engagement(non_responders, as.character(non_responders$intervention3), main_effect_size, hte_effect_size)
    pop[pop$engagement == 0, names(non_responders)] <- non_responders
  } else {
    pop$intervention3 <- factor(c("LLM", "human")[1])
  }
  
  return(pop)
}

calculate_cost_effectiveness <- function(n, method, main_effect_size, hte_effect_size, annual_participants = 10000) {
  fixed_costs <- 89750 # Annual fixed costs
  fixed_cost_per_participant <- fixed_costs / annual_participants
  
  if (method == "ab") {
    pop_list <- simulate_ab_tests(n, main_effect_size, hte_effect_size)
    
    cost <- sum(sapply(pop_list, function(pop) {
      sum(pop$engagement * ifelse(pop$intervention == "empathetic", 1.23,
                                  ifelse(pop$intervention == "factual", 0.63,
                                         ifelse(pop$intervention == "weekday" | pop$intervention == "weekend", 0.03,
                                                ifelse(pop$intervention == "LLM", 1.23, 5.63)))))
    }))
    
    effectiveness <- sum(sapply(pop_list, function(pop) sum(pop$engagement))) / (3 * n)
    
  } else { # SMART design
    pop <- simulate_smart(n, main_effect_size, hte_effect_size)
    
    cost <- sum(pop$engagement * 1.23) # Stage 1
    non_responders <- pop[pop$engagement == 0, ]
    if (nrow(non_responders) > 0) {
      cost <- cost + nrow(non_responders) * 0.03 # Stage 2
      non_responders <- non_responders[non_responders$engagement == 0, ]
      if (nrow(non_responders) > 0) {
        cost <- cost + sum(non_responders$engagement * ifelse(non_responders$intervention3 == "LLM", 1.23, 5.63)) # Stage 3
      }
    }
    
    effectiveness <- sum(pop$engagement) / n
  }
  
  total_cost <- (fixed_cost_per_participant * n) + cost
  
  if (effectiveness == 0) {
    return(list(cost = total_cost, effectiveness = effectiveness, cost_effectiveness = Inf))
  } else {
    return(list(cost = total_cost, effectiveness = effectiveness, cost_effectiveness = total_cost / (effectiveness * n)))
  }
}

calculate_net_benefit <- function(pop, threshold) {
  TP <- sum(pop$engagement == 1)
  FP <- sum(pop$engagement == 0)
  N <- nrow(pop)
  
  NB <- (TP/N) - (FP/N) * (threshold / (1 - threshold))
  return(NB)
}

calculate_power_and_fpr <- function(n, method, main_effect_size, hte_effect_size, total_tests = 100) {
  true_positive_age <- 0
  true_positive_race <- 0
  true_positive_disease <- 0
  false_positive_age <- 0
  false_positive_race <- 0
  false_positive_disease <- 0
  
  for (i in 1:total_tests) {
    if (method == "ab") {
      pop_list <- simulate_ab_tests(n, main_effect_size, hte_effect_size)
      
      for (j in 1:length(pop_list)) {
        pop <- pop_list[[j]]
        model <- glm(engagement ~ intervention + age + chronic_disease + race_ethnicity + area_residence +
                       intervention:age + intervention:chronic_disease + intervention:race_ethnicity,
                     data = pop, family = binomial, method = "glm.fit", epsilon = 1e-8)
        
        p_values <- coef(summary(model))[, "Pr(>|z|)"]
        true_positive_age <- true_positive_age + any(p_values[grepl("intervention(empathetic|factual):age", names(p_values))] < 0.05, na.rm = TRUE)
        true_positive_race <- true_positive_race + any(p_values[grepl("intervention(weekday|weekend):race_ethnicity", names(p_values))] < 0.05, na.rm = TRUE)
        true_positive_disease <- true_positive_disease + any(p_values[grepl("intervention(LLM|human):chronic_disease", names(p_values))] < 0.05, na.rm = TRUE)
        false_positive_age <- false_positive_age + any(p_values[grepl("intervention(weekday|weekend|LLM|human):age", names(p_values))] < 0.05, na.rm = TRUE)
        false_positive_race <- false_positive_race + any(p_values[grepl("intervention(empathetic|factual|LLM|human):race_ethnicity", names(p_values))] < 0.05, na.rm = TRUE)
        false_positive_disease <- false_positive_disease + any(p_values[grepl("intervention(empathetic|factual|weekday|weekend):chronic_disease", names(p_values))] < 0.05, na.rm = TRUE)
      }
      
    } else {
      pop <- simulate_smart(n, main_effect_size, hte_effect_size)
      
      model <- glm(engagement ~ intervention1 + intervention2 + intervention3 +
                     age + chronic_disease + race_ethnicity + area_residence +
                     intervention1:age + intervention2:age + intervention3:age +
                     intervention1:race_ethnicity + intervention2:race_ethnicity + intervention3:race_ethnicity +
                     intervention1:chronic_disease + intervention2:chronic_disease + intervention3:chronic_disease,
                   data = pop, family = binomial, method = "glm.fit", epsilon = 1e-8)
      
      p_values <- coef(summary(model))[, "Pr(>|z|)"]
      true_positive_age <- true_positive_age + any(p_values[grepl("intervention1(empathetic|factual):age", names(p_values))] < 0.05, na.rm = TRUE)
      true_positive_race <- true_positive_race + any(p_values[grepl("intervention2(weekday|weekend):race_ethnicity", names(p_values))] < 0.05, na.rm = TRUE)
      true_positive_disease <- true_positive_disease + any(p_values[grepl("intervention3(LLM|human):chronic_disease", names(p_values))] < 0.05, na.rm = TRUE)
      false_positive_age <- false_positive_age + any(p_values[grepl("intervention2(weekday|weekend):age|intervention3(LLM|human):age", names(p_values))] < 0.05, na.rm = TRUE)
      false_positive_race <- false_positive_race + any(p_values[grepl("intervention1(empathetic|factual):race_ethnicity|intervention3(LLM|human):race_ethnicity", names(p_values))] < 0.05, na.rm = TRUE)
      false_positive_disease <- false_positive_disease + any(p_values[grepl("intervention1(empathetic|factual):chronic_disease|intervention2(weekday|weekend):chronic_disease", names(p_values))] < 0.05, na.rm = TRUE)
    }
  }
  
  if (method == "ab") {
    power_age <- true_positive_age / (total_tests * length(pop_list))
    power_race <- true_positive_race / (total_tests * length(pop_list))
    power_disease <- true_positive_disease / (total_tests * length(pop_list))
    fpr_age <- false_positive_age / (total_tests * length(pop_list)/3)
    fpr_race <- false_positive_race / (total_tests * length(pop_list)/3)
    fpr_disease <- false_positive_disease / (total_tests * length(pop_list)/3)
  } else {
    power_age <- true_positive_age / total_tests
    power_race <- true_positive_race / total_tests
    power_disease <- true_positive_disease / total_tests
    fpr_age <- false_positive_age / total_tests
    fpr_race <- false_positive_race / total_tests
    fpr_disease <- false_positive_disease / total_tests
  }
  
  return(list(power_age = power_age, power_race = power_race, power_disease = power_disease,
              fpr_age = fpr_age, fpr_race = fpr_race, fpr_disease = fpr_disease))
}

run_simulations <- function(n, method, main_effect_size, hte_effect_size, num_bootstraps = 10, total_tests = 100) {
  results <- replicate(num_bootstraps, {
    power_and_fpr <- calculate_power_and_fpr(n, method, main_effect_size, hte_effect_size, total_tests)
    cost_effectiveness <- calculate_cost_effectiveness(n, method, main_effect_size, hte_effect_size)
    
    # Calculate net benefit for a range of thresholds
    if (method == "ab") {
      pop_list <- simulate_ab_tests(n, main_effect_size, hte_effect_size)
      pop <- do.call(rbind, pop_list)
    } else {
      pop <- simulate_smart(n, main_effect_size, hte_effect_size)
    }
    net_benefits <- sapply(seq(0.05, 0.30, by = 0.01), function(threshold) calculate_net_benefit(pop, threshold))
    
    c(power_and_fpr$power_age, power_and_fpr$power_race, power_and_fpr$power_disease,
      power_and_fpr$fpr_age, power_and_fpr$fpr_race, power_and_fpr$fpr_disease,
      cost_effectiveness$cost_effectiveness, net_benefits)
  })
  
  power_age_ci <- quantile(results[1, ], c(0.025, 0.975), na.rm = TRUE)
  power_race_ci <- quantile(results[2, ], c(0.025, 0.975), na.rm = TRUE)
  power_disease_ci <- quantile(results[3, ], c(0.025, 0.975), na.rm = TRUE)
  fpr_age_ci <- quantile(results[4, ], c(0.025, 0.975), na.rm = TRUE)
  fpr_race_ci <- quantile(results[5, ], c(0.025, 0.975), na.rm = TRUE)
  fpr_disease_ci <- quantile(results[6, ], c(0.025, 0.975), na.rm = TRUE)
  ce_ci <- quantile(results[7, ], c(0.025, 0.975), na.rm = TRUE)
  
  data.frame(
    method = method,
    n = n,
    main_effect_size = main_effect_size,
    hte_effect_size = hte_effect_size,
    power_age = mean(results[1, ], na.rm = TRUE),
    power_age_lower = power_age_ci[1],
    power_age_upper = power_age_ci[2],
    power_race = mean(results[2, ], na.rm = TRUE),
    power_race_lower = power_race_ci[1],
    power_race_upper = power_race_ci[2],
    power_disease = mean(results[3, ], na.rm = TRUE),
    power_disease_lower = power_disease_ci[1],
    power_disease_upper = power_disease_ci[2],
    fpr_age = mean(results[4, ], na.rm = TRUE),
    fpr_age_lower = fpr_age_ci[1],
    fpr_age_upper = fpr_age_ci[2],
    fpr_race = mean(results[5, ], na.rm = TRUE),
    fpr_race_lower = fpr_race_ci[1],
    fpr_race_upper = fpr_race_ci[2],
    fpr_disease = mean(results[6, ], na.rm = TRUE),
    fpr_disease_lower = fpr_disease_ci[1],
    fpr_disease_upper = fpr_disease_ci[2],
    cost_effectiveness = mean(results[7, ], na.rm = TRUE),
    ce_lower = ce_ci[1],
    ce_upper = ce_ci[2],
    net_benefit = colMeans(results[8:33, ], na.rm = TRUE)
  )
}

plot_results <- function(min_n = 100, max_n = 1000, step_size = 100, main_effect_sizes, hte_effect_sizes) {
  sample_sizes <- seq(min_n, max_n, by = step_size)
  results <- bind_rows(
    lapply(sample_sizes, function(n) {
      bind_rows(
        lapply(main_effect_sizes, function(main_effect_size) {
          lapply(hte_effect_sizes, function(hte_effect_size) {
            bind_rows(
              run_simulations(n, "ab", main_effect_size, hte_effect_size),
              run_simulations(n, "smart", main_effect_size, hte_effect_size)
            )
          }) %>% bind_rows()
        }) %>% bind_rows()
      )
    }) %>% bind_rows()
  )
  
  power_age_plot <- ggplot(results, aes(x = n, y = power_age, color = method)) +
    geom_line() +
    geom_ribbon(aes(ymin = power_age_lower, ymax = power_age_upper, fill = method), alpha = 0.2) + ylim(0, 1) +
    labs(x = "Sample Size", y = "Power (Age HTE)", title = "Power vs Sample Size (Age HTE)") +
    theme_minimal() +
    facet_grid(main_effect_size ~ hte_effect_size, labeller = label_both)
  
  power_race_plot <- ggplot(results, aes(x = n, y = power_race, color = method)) +
    geom_line() +
    geom_ribbon(aes(ymin = power_race_lower, ymax = power_race_upper, fill = method), alpha = 0.2) + ylim(0, 1) +
    labs(x = "Sample Size", y = "Power (Race/Ethnicity HTE)", title = "Power vs Sample Size (Race/Ethnicity HTE)") +
    theme_minimal() +
    facet_grid(main_effect_size ~ hte_effect_size, labeller = label_both)
  
  power_disease_plot <- ggplot(results, aes(x = n, y = power_disease, color = method)) +
    geom_line() +
    geom_ribbon(aes(ymin = power_disease_lower, ymax = power_disease_upper, fill = method), alpha = 0.2) + ylim(0, 1) +
    labs(x = "Sample Size", y = "Power (Chronic Disease HTE)", title = "Power vs Sample Size (Chronic Disease HTE)") +
    theme_minimal() +
    facet_grid(main_effect_size ~ hte_effect_size, labeller = label_both)
  
  fpr_age_plot <- ggplot(results, aes(x = n, y = fpr_age, color = method)) +
    geom_line() +
    geom_ribbon(aes(ymin = fpr_age_lower, ymax = fpr_age_upper, fill = method), alpha = 0.2) + ylim(0, 1) +
    labs(x = "Sample Size", y = "False Positive Rate (Age HTE)", title = "False Positive Rate vs Sample Size (Age HTE)") +
    theme_minimal() +
    facet_grid(main_effect_size ~ hte_effect_size, labeller = label_both)
  
  fpr_race_plot <- ggplot(results, aes(x = n, y = fpr_race, color = method)) +
    geom_line() +
    geom_ribbon(aes(ymin = fpr_race_lower, ymax = fpr_race_upper, fill = method), alpha = 0.2) + ylim(0, 1) +
    labs(x = "Sample Size", y = "False Positive Rate (Race/Ethnicity HTE)", title = "False Positive Rate vs Sample Size (Race/Ethnicity HTE)") +
    theme_minimal() +
    facet_grid(main_effect_size ~ hte_effect_size, labeller = label_both)
  
  fpr_disease_plot <- ggplot(results, aes(x = n, y = fpr_disease, color = method)) +
    geom_line() +
    geom_ribbon(aes(ymin = fpr_disease_lower, ymax = fpr_disease_upper, fill = method), alpha = 0.2) + ylim(0, 1) +
    labs(x = "Sample Size", y = "False Positive Rate (Chronic Disease HTE)", title = "False Positive Rate vs Sample Size (Chronic Disease HTE)") +
    theme_minimal() +
    facet_grid(main_effect_size ~ hte_effect_size, labeller = label_both)
  
  ce_plot <- ggplot(results, aes(x = n, y = cost_effectiveness, color = method)) +
    geom_line() +
    geom_ribbon(aes(ymin = ce_lower, ymax = ce_upper, fill = method), alpha = 0.2) +
    labs(x = "Sample Size", y = "Cost-Effectiveness ($/patient engaged)", 
         title = "Cost-Effectiveness vs Sample Size") +
    theme_minimal() +
    facet_grid(main_effect_size ~ hte_effect_size, labeller = label_both)
  
  nb_data <- results %>%
    pivot_longer(cols = starts_with("net_benefit"), 
                 names_to = "threshold", 
                 values_to = "net_benefit") %>%
    mutate(threshold = as.numeric(sub("net_benefit", "", threshold)) / 100)
  
  nb_plot <- ggplot(nb_data, aes(x = threshold, y = net_benefit, color = method, group = interaction(method, n, main_effect_size, hte_effect_size))) +
    geom_line() +
    labs(x = "Threshold Probability", y = "Net Benefit", 
         title = "Decision Curve Analysis") +
    theme_minimal() +
    facet_grid(main_effect_size ~ hte_effect_size, labeller = label_both) +
    scale_x_continuous(limits = c(0.05, 0.30), breaks = seq(0.05, 0.30, by = 0.05))
  
  list(power_age_plot = power_age_plot, 
       power_race_plot = power_race_plot,
       power_disease_plot = power_disease_plot, 
       fpr_age_plot = fpr_age_plot,
       fpr_race_plot = fpr_race_plot, 
       fpr_disease_plot = fpr_disease_plot,
       ce_plot = ce_plot,
       nb_plot = nb_plot)
}

# Run simulations and generate plots
main_effect_sizes <- c(0.1, 0.2, 0.3)
hte_effect_sizes <- c(0.05, 0.1, 0.15)
plots <- plot_results(min_n = 100, max_n = 1000, step_size = 100, main_effect_sizes, hte_effect_sizes)

# Save plots
ggsave("power_age_plot.png", plot = plots$power_age_plot, width = 8, height = 6, dpi = 300)
ggsave("power_race_plot.png", plot = plots$power_race_plot, width = 8, height = 6, dpi = 300)
ggsave("power_disease_plot.png", plot = plots$power_disease_plot, width = 8, height = 6, dpi = 300)
ggsave("fpr_age_plot.png", plot = plots$fpr_age_plot, width = 8, height = 6, dpi = 300)
ggsave("fpr_race_plot.png", plot = plots$fpr_race_plot, width = 8, height = 6, dpi = 300)
ggsave("fpr_disease_plot.png", plot = plots$fpr_disease_plot, width = 8, height = 6, dpi = 300)
ggsave("ce_plot.png", plot = plots$ce_plot, width = 8, height = 6, dpi = 300)
ggsave("nb_plot.png", plot = plots$nb_plot, width = 8, height = 6, dpi = 300)

# Generate results table
results_table <- bind_rows(
  lapply(main_effect_sizes, function(main_effect_size) {
    lapply(hte_effect_sizes, function(hte_effect_size) {
      bind_rows(
        run_simulations(1000, "ab", main_effect_size, hte_effect_size),
        run_simulations(1000, "smart", main_effect_size, hte_effect_size)
      )
    }) %>% bind_rows()
  }) %>% bind_rows()
)
print(results_table)
write.csv(results_table, file = "results_updated.csv", row.names = FALSE)
