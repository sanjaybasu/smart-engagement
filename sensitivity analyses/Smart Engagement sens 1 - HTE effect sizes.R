# By running this updated code, you can compare the performance of the A/B test and SMART trial under different disease HTE effect sizes while keeping the age and race HTE effect sizes constant. This will help you assess whether the SMART trial's advantage in detecting the disease HTE becomes more pronounced as the effect size increases, providing insight into the targeted population factor (answer 2).

library(dplyr)
library(tidyverse)
library(parallel)

simulate_population <- function(n) {
  data.frame(
    id = 1:n,
    age = sample(18:90, n, replace = TRUE),
    chronic_disease = sample(c(0, 1), n, replace = TRUE),
    race_ethnicity = sample(c("White", "Black", "Hispanic", "Other"), n, replace = TRUE),
    area_residence = sample(c("Urban", "Rural"), n, replace = TRUE)
  )
}

simulate_engagement <- function(pop, intervention, main_effect_size, hte_effect_sizes, noise_scale = 0.01) {
  base_prob <- 0.1
  
  intervention_effect <- numeric(nrow(pop))
  hte <- numeric(nrow(pop))
  
  empathetic_mask <- intervention == "empathetic"
  factual_mask <- intervention == "factual"
  weekday_mask <- intervention == "weekday"
  weekend_mask <- intervention == "weekend"
  human_mask <- intervention == "human"
  llm_mask <- intervention == "LLM"
  
  intervention_effect[empathetic_mask | factual_mask] <- main_effect_size
  hte[empathetic_mask] <- ifelse(pop$age[empathetic_mask] > 50, hte_effect_sizes[1], 0)
  hte[factual_mask] <- ifelse(pop$age[factual_mask] <= 50, hte_effect_sizes[1], 0)
  
  intervention_effect[weekday_mask | weekend_mask] <- main_effect_size
  hte[weekday_mask] <- ifelse(pop$race_ethnicity[weekday_mask] %in% c("Black", "Hispanic"), hte_effect_sizes[2], 0)
  hte[weekend_mask] <- ifelse(!(pop$race_ethnicity[weekend_mask] %in% c("Black", "Hispanic")), hte_effect_sizes[2], 0)
  
  intervention_effect[human_mask | llm_mask] <- main_effect_size
  hte[human_mask] <- ifelse(pop$chronic_disease[human_mask] == 1, hte_effect_sizes[3], 0)
  hte[llm_mask] <- ifelse(pop$chronic_disease[llm_mask] == 0, hte_effect_sizes[3], 0)
  
  noise <- rnorm(nrow(pop), mean = 0, sd = noise_scale)
  
  prob <- pmin(pmax(base_prob + intervention_effect + hte + noise, 0), 1)
  engagement <- rbinom(nrow(pop), 1, prob)
  
  pop$engagement <- engagement
  return(pop)
}

simulate_ab_tests <- function(n, main_effect_size, hte_effect_sizes) {
  pop <- simulate_population(round(n/3))
  
  pop_test1 <- pop
  pop_test1$intervention <- sample(c("empathetic", "factual"), round(n/3), replace = TRUE)
  pop_test1 <- simulate_engagement(pop_test1, pop_test1$intervention, main_effect_size, hte_effect_sizes)
  
  pop_test2 <- pop
  pop_test2$intervention <- sample(c("weekday", "weekend"), round(n/3), replace = TRUE)
  pop_test2 <- simulate_engagement(pop_test2, pop_test2$intervention, main_effect_size, hte_effect_sizes)
  
  pop_test3 <- pop
  pop_test3$intervention <- sample(c("LLM", "human"), round(n/3), replace = TRUE)
  pop_test3 <- simulate_engagement(pop_test3, pop_test3$intervention, main_effect_size, hte_effect_sizes)
  
  return(list(pop_test1, pop_test2, pop_test3))
}

simulate_smart <- function(n, main_effect_size, hte_effect_sizes) {
  pop <- simulate_population(n)
  
  pop$intervention1 <- factor(sample(c("empathetic", "factual"), n, replace = TRUE))
  pop <- simulate_engagement(pop, as.character(pop$intervention1), main_effect_size, hte_effect_sizes)
  
  non_responders <- pop[pop$engagement == 0, ]
  if (nrow(non_responders) > 0) {
    non_responders$intervention2 <- factor(sample(c("weekday", "weekend"), nrow(non_responders), replace = TRUE))
    non_responders <- simulate_engagement(non_responders, as.character(non_responders$intervention2), main_effect_size, hte_effect_sizes)
    pop[pop$engagement == 0, names(non_responders)] <- non_responders
  } else {
    pop$intervention2 <- factor(c("weekday", "weekend")[1])
  }
  
  non_responders <- pop[pop$engagement == 0, ]
  if (nrow(non_responders) > 0) {
    non_responders$intervention3 <- factor(sample(c("LLM", "human"), nrow(non_responders), replace = TRUE))
    non_responders <- simulate_engagement(non_responders, as.character(non_responders$intervention3), main_effect_size, hte_effect_sizes)
    pop[pop$engagement == 0, names(non_responders)] <- non_responders
  } else {
    pop$intervention3 <- factor(c("LLM", "human")[1])
  }
  
  return(pop)
}

calculate_cost_effectiveness <- function(n, method, main_effect_size, hte_effect_sizes) {
  pop <- simulate_population(n)
  
  if (method == "ab") {
    pop_empathetic <- simulate_engagement(pop, "empathetic", main_effect_size, hte_effect_sizes)
    pop_factual <- simulate_engagement(pop, "factual", main_effect_size, hte_effect_sizes)
    pop_weekday <- simulate_engagement(pop, "weekday", main_effect_size, hte_effect_sizes)
    pop_weekend <- simulate_engagement(pop, "weekend", main_effect_size, hte_effect_sizes)
    pop_llm <- simulate_engagement(pop, "LLM", main_effect_size, hte_effect_sizes)
    pop_human <- simulate_engagement(pop, "human", main_effect_size, hte_effect_sizes)
    
    cost <- sum(pop_empathetic$engagement) * 1.13 + 
      sum(pop_factual$engagement) * 0.60 +
      sum(pop_weekday$engagement) * 0.03 +
      sum(pop_weekend$engagement) * 0.10 +
      sum(pop_llm$engagement) * 0.60 +
      sum(pop_human$engagement) * 5.60
    
    effectiveness <- length(unique(c(
      pop_empathetic$id[pop_empathetic$engagement == 1],
      pop_factual$id[pop_factual$engagement == 1],
      pop_weekday$id[pop_weekday$engagement == 1],
      pop_weekend$id[pop_weekend$engagement == 1],
      pop_llm$id[pop_llm$engagement == 1],
      pop_human$id[pop_human$engagement == 1]
    ))) / n
    
  } else {
    pop_empathetic <- simulate_engagement(pop, "empathetic", main_effect_size, hte_effect_sizes)
    non_responders_weekday <- pop[pop_empathetic$engagement == 0, ]
    non_responders_weekend <- pop[pop_empathetic$engagement == 0, ]
    non_responders_human <- pop[pop_empathetic$engagement == 0, ]
    
    if (nrow(non_responders_weekday) > 0) {
      non_responders_weekday <- simulate_engagement(non_responders_weekday, "weekday", main_effect_size, hte_effect_sizes)
    }
    
    if (nrow(non_responders_weekend) > 0) {
      non_responders_weekend <- simulate_engagement(non_responders_weekend, "weekend", main_effect_size, hte_effect_sizes)
    }
    
    if (nrow(non_responders_human) > 0) {
      non_responders_human <- simulate_engagement(non_responders_human, "human", main_effect_size, hte_effect_sizes)
    }
    
    cost <- sum(pop_empathetic$engagement) * 1.13 +
      sum(non_responders_weekday$engagement) * 0.03 +
      sum(non_responders_weekend$engagement) * 0.10 +
      sum(non_responders_human$engagement) * 5.60
    
    effectiveness <- length(unique(c(
      pop_empathetic$id[pop_empathetic$engagement == 1],
      non_responders_weekday$id[non_responders_weekday$engagement == 1],
      non_responders_weekend$id[non_responders_weekend$engagement == 1],
      non_responders_human$id[non_responders_human$engagement == 1]
    ))) / n
  }
  
  if (effectiveness == 0) {
    return(list(cost = cost, effectiveness = effectiveness, cost_effectiveness = Inf))
  } else {
    return(list(cost = cost, effectiveness = effectiveness, cost_effectiveness = cost / effectiveness))
  }
}

calculate_power_and_fpr <- function(n, method, main_effect_size, hte_effect_sizes, total_tests = 100) {
  true_positive_age <- 0
  true_positive_race <- 0
  true_positive_disease <- 0
  false_positive_age <- 0
  false_positive_race <- 0
  false_positive_disease <- 0
  
  for (i in 1:total_tests) {
    if (method == "ab") {
      pop_list <- simulate_ab_tests(n, main_effect_size, hte_effect_sizes)
      
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
      pop <- simulate_smart(n, main_effect_size, hte_effect_sizes)
      
      model <- glm(engagement ~ intervention1 + intervention2 + intervention3 +
                     age + chronic_disease + race_ethnicity + area_residence +
                     intervention1:age + intervention2:race_ethnicity + intervention3:chronic_disease,
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

run_simulations <- function(n, method, main_effect_size, hte_effect_sizes, num_bootstraps = 10, total_tests = 100) {
  results <- replicate(num_bootstraps, {
    power_and_fpr <- calculate_power_and_fpr(n, method, main_effect_size, hte_effect_sizes, total_tests)
    cost_effectiveness <- calculate_cost_effectiveness(n, method, main_effect_size, hte_effect_sizes)
    c(power_and_fpr$power_age, power_and_fpr$power_race, power_and_fpr$power_disease,
      power_and_fpr$fpr_age, power_and_fpr$fpr_race, power_and_fpr$fpr_disease,
      cost_effectiveness$cost_effectiveness)
  })
  
  power_age_ci <- quantile(results[1, ], c(0.025, 0.975), na.rm = TRUE)
  power_race_ci <- quantile(results[2, ], c(0.025, 0.975), na.rm = TRUE)
  power_disease_ci <- quantile(results[3, ], c(0.025, 0.975), na.rm = TRUE)
  fpr_age_ci <- quantile(results[4, ], c(0.025, 0.975), na.rm = TRUE)
  fpr_race_ci <- quantile(results[5, ], c(0.025, 0.975), na.rm = TRUE)
  fpr_disease_ci <- quantile(results[6, ], c(0.025, 0.975), na.rm = TRUE)
  ce_ci <- quantile(results[7, ], c(0.025, 0.975), na.rm = TRUE)
  
  # Separate HTE effect sizes in the results table
  data.frame(
    method = method,
    n = n,
    main_effect_size = main_effect_size,
    hte_effect_size_age = hte_effect_sizes[1],  # Separate HTE effect sizes
    hte_effect_size_race = hte_effect_sizes[2],
    hte_effect_size_disease = hte_effect_sizes[3],
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
    ce_upper = ce_ci[2]
  )
}

plot_results <- function(results_df) {
  # Correctly display and save plots
  
  # Power Plots
  power_age_plot <- ggplot(results_df, aes(x = as.factor(hte_effect_size_disease), y = power_age, color = method)) +
    geom_line(aes(group = method)) + 
    geom_ribbon(aes(ymin = power_age_lower, ymax = power_age_upper, fill = method), alpha = 0.2) +
    labs(x = "Disease HTE Effect Size", y = "Power (Age HTE)", title = "Power vs Disease HTE Effect Size (Age HTE)") +
    theme_minimal()
  
  power_race_plot <- ggplot(results_df, aes(x = as.factor(hte_effect_size_disease), y = power_race, color = method)) +
    geom_line(aes(group = method)) + 
    geom_ribbon(aes(ymin = power_race_lower, ymax = power_race_upper, fill = method), alpha = 0.2) +
    labs(x = "Disease HTE Effect Size", y = "Power (Race/Ethnicity HTE)", title = "Power vs Disease HTE Effect Size (Race/Ethnicity HTE)") +
    theme_minimal()
  
  power_disease_plot <- ggplot(results_df, aes(x = as.factor(hte_effect_size_disease), y = power_disease, color = method)) +
    geom_line(aes(group = method)) + 
    geom_ribbon(aes(ymin = power_disease_lower, ymax = power_disease_upper, fill = method), alpha = 0.2) +
    labs(x = "Disease HTE Effect Size", y = "Power (Chronic Disease HTE)", title = "Power vs Disease HTE Effect Size (Chronic Disease HTE)") +
    theme_minimal()
  
  # False Positive Rate Plots
  fpr_age_plot <- ggplot(results_df, aes(x = as.factor(hte_effect_size_disease), y = fpr_age, color = method)) +
    geom_line(aes(group = method)) +
    geom_ribbon(aes(ymin = fpr_age_lower, ymax = fpr_age_upper, fill = method), alpha = 0.2) +
    labs(x = "Disease HTE Effect Size", y = "False Positive Rate (Age HTE)", title = "False Positive Rate vs Disease HTE Effect Size (Age HTE)") +
    theme_minimal()
  
  fpr_race_plot <- ggplot(results_df, aes(x = as.factor(hte_effect_size_disease), y = fpr_race, color = method)) +
    geom_line(aes(group = method)) +
    geom_ribbon(aes(ymin = fpr_race_lower, ymax = fpr_race_upper, fill = method), alpha = 0.2) +
    labs(x = "Disease HTE Effect Size", y = "False Positive Rate (Race/Ethnicity HTE)", title = "False Positive Rate vs Disease HTE Effect Size (Race/Ethnicity HTE)") +
    theme_minimal()
  
  fpr_disease_plot <- ggplot(results_df, aes(x = as.factor(hte_effect_size_disease), y = fpr_disease, color = method)) +
    geom_line(aes(group = method)) +
    geom_ribbon(aes(ymin = fpr_disease_lower, ymax = fpr_disease_upper, fill = method), alpha = 0.2) +
    labs(x = "Disease HTE Effect Size", y = "False Positive Rate (Chronic Disease HTE)", title = "False Positive Rate vs Disease HTE Effect Size (Chronic Disease HTE)") +
    theme_minimal()
  
  # Cost-Effectiveness Plot
  ce_plot <- ggplot(results_df, aes(x = as.factor(hte_effect_size_disease), y = cost_effectiveness, color = method)) +
    geom_line(aes(group = method)) +
    geom_ribbon(aes(ymin = ce_lower, ymax = ce_upper, fill = method), alpha = 0.2) +
    labs(x = "Disease HTE Effect Size", y = "Cost-Effectiveness", title = "Cost-Effectiveness vs Disease HTE Effect Size") +
    theme_minimal()
  
  # Print plots
  print(power_age_plot)
  print(power_race_plot)
  print(power_disease_plot)
  print(fpr_age_plot)
  print(fpr_race_plot)
  print(fpr_disease_plot)
  print(ce_plot)
  
  # Save plots
  ggsave("power_age_plot1.png", plot = power_age_plot, width = 8, height = 6, dpi = 300)
  ggsave("power_race_plot1.png", plot = power_race_plot, width = 8, height = 6, dpi = 300)
  ggsave("power_disease_plot1.png", plot = power_disease_plot, width = 8, height = 6, dpi = 300)
  ggsave("fpr_age_plot1.png", plot = fpr_age_plot, width = 8, height = 6, dpi = 300)
  ggsave("fpr_race_plot1.png", plot = fpr_race_plot, width = 8, height = 6, dpi = 300)
  ggsave("fpr_disease_plot1.png", plot = fpr_disease_plot, width = 8, height = 6, dpi = 300)
  ggsave("ce_plot1.png", plot = ce_plot, width = 8, height = 6, dpi = 300)
}

# --- Main Analysis ---

main_effect_size <- 0.2
hte_effect_sizes_list <- list(
  c(0.05, 0.05, 0.05),
  c(0.05, 0.05, 0.10),
  c(0.05, 0.05, 0.15)
)

n <- 5000

# Generate results table
results_table <- bind_rows(
  lapply(hte_effect_sizes_list, function(hte_effect_sizes) {
    bind_rows(
      run_simulations(n, "ab", main_effect_size, hte_effect_sizes),
      run_simulations(n, "smart", main_effect_size, hte_effect_sizes)
    )
  })
)

# Generate and display plots
plot_results(results_table)

print(results_table)
write.csv(results_table, file = "results_sens1.csv", row.names = FALSE)
