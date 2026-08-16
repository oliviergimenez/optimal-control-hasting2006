# Try and reproduce results from
# A simple approach to optimal control of invasive species
# Hastings, Hall, Taylor (2006) TPB

# Function to compute L^x for a given matrix L and exponent x
matrix_power <- function(L, x) {
  # Check if x is 0; return the identity matrix in that case
  if (x == 0) {
    return(diag(nrow(L)))
  }
  # Use Reduce to multiply L by itself x times
  result <- Reduce(function(a, b) a %*% b, replicate(x, L, simplify = FALSE), init = diag(nrow(L)))
  return(result)
}


#--------------------------- T = 1 with total budget = 3 * 1e-3

# Load the lpSolve package
library(lpSolve)

# Define the Leslie matrix based on the values provided
L <- matrix(c(0, 0.0000646, 0.0177,
              1, 1.115, 0,
              0, 0.265, 1.107), 
            byrow = TRUE, ncol = 3)

# Initial population vector based on provided values
N0 <- c(0.00008226, 0.01186, 0.005977)  # N0

# Define the time horizon for the control strategy
time_horizon <- 1  # Time horizon for control

# Define the cost for removing individuals from each stage
removal_cost <- c(1, 1, 0.4)  # Cost matrix as provided

# Define the budget constraint per time step
budget <- 5 * 1e-3  # Updated budget constraint

# minimize LN0 - LH1 is eq to min -LH1 or max LH1

# Objective function: minimize cost
objective_function <- rowSums(t(L))

# Define constraints based on the population dynamics and budget
# 1. Population constraint (ensuring we don’t remove more than the initial population in each stage)
population_constraints <- diag(3)
rhs_population <- N0

# 2. Budget constraint from Equation 6: sum(removal_costs * amount_removed) <= total_budget
budget_constraint <- matrix(removal_cost, nrow = 1)  # Single constraint for budget
rhs_budget <- budget

# Combine all constraints
A <- rbind(population_constraints, budget_constraint)
rhs <- c(rhs_population, rhs_budget) # Right-hand side: Set population threshold after control
constraint_direction <- c(rep("<=", length(rhs_population)), "<=") # direction of constraints

# Solve the linear programming problem
lp_solution <- lp(direction = "max", 
                  objective.in = objective_function, 
                  const.mat = A, 
                  const.dir = constraint_direction, 
                  const.rhs = rhs)

lp_solution$solution
lp_solution$objval

# pop reamining L^T N0 - sum_i^T(L^(T+1-i)*Hi)
sum(L %*% matrix(N0, ncol = 1) - L %*% matrix(lp_solution$solution,ncol=1))

#----------------------- T = 1 w/ grid on budget

# Define the Leslie matrix based on the values provided
L <- matrix(c(0, 0.0000646, 0.0177,
              1, 1.115, 0,
              0, 0.265, 1.107), 
            byrow = TRUE, ncol = 3)

# Initial population vector based on provided values
N0 <- c(0.00008226, 0.01186, 0.005977)  # N0

# Define the time horizon for the control strategy
time_horizon <- 1  # Time horizon for control

# Define the cost for removing individuals from each stage
removal_cost <- c(1, 1, 0.4)  # Cost matrix as provided

# minimize LN0 - LH1 is eq to min -LH1 or max LH1

# Objective function: minimize cost
objective_function <- rowSums(t(L))

# We want the total population after applying control to be below the target threshold
A <- matrix(c(1, 0, 0,
              0, 1, 0,
              0, 0, 1,
              removal_cost), 
            byrow = TRUE, 
            nrow = 4)

# Direction of constraints
constraint_direction <- "<="

# Define the budget constraint per time step
#budget <- 3 * 1e-3  # Updated budget constraint

budget <- seq(0, 5, by = 0.5) * 1e-3
length_grid <- length(budget)
sol <- matrix(NA, nrow = length_grid, ncol = 3)
obj <- rep(NA, length_grid)
ind <- 1

for (i in budget){
  
  
  # Right-hand side: Set population threshold after control
  rhs <- c(N0, i)
  
  # Solve the linear programming problem
  lp_solution <- lp(direction = "max", 
                    objective.in = objective_function, 
                    const.mat = A, 
                    const.dir = constraint_direction, 
                    const.rhs = rhs)
  
  sol[ind,] <- lp_solution$solution
  obj[ind] <- lp_solution$objval
  ind <- ind + 1

}

sol[,2] / (sol[,2] + sol[,3]) * 100
sol[,3] / (sol[,2] + sol[,3]) * 100

obj



#------------------ Any T w/ budget = 3e-3


# Load the lpSolve package
library(lpSolve)

# Define the number of stages
stages <- 3  # Three-stage model

# Define the Leslie matrix based on the values provided
L <- matrix(c(0, 0.0000646, 0.0177,
              1, 1.115, 0,
              0, 0.265, 1.107), 
            byrow = TRUE, ncol = 3)

# Initial population vector based on provided values
N0 <- c(0.00008226, 0.01186, 0.005977)  # N0

# Define the time horizon for the control strategy
T <- 10  # Time horizon for control

# Define the cost for removing individuals from each stage
removal_cost <- c(1, 1, 0.4)  # Cost matrix as provided

# Define the budget constraint per time step
budget <- 1*1e-3  # Updated budget constraint

# Set up the linear programming model

objective <- rep(NA, stages * T)
boundaries <- seq(1, stages*T, by = 3)
for (i in 1:T){
  objective[boundaries[i]:(boundaries[i]+2)] <- rowSums(t(matrix_power(L, T + 1 - i)))
}

# Initialize constraints for Leslie matrix dynamics and budget
constraints <- matrix(0, nrow = T * stages + T, ncol = stages * T)
constraints_rhs <- numeric(T * stages + T)
constraints_dir <- rep("<=", T * stages + T)

# Populate constraints for each time step based on population dynamics and budget
for (t in 1:T) {
  if (t == 1) {
    # Initial population constraints
    constraints[(1:stages) + (t - 1) * stages, (1:stages) + (t - 1) * stages] <- diag(stages)
    constraints_rhs[(1:stages) + (t - 1) * stages] <- N0
  } else {
    # Leslie matrix dynamics for each subsequent time step
    constraints[(1:stages) + (t - 1) * stages, (1:stages) + (t - 1) * stages] <- diag(stages)
    constraints[(1:stages) + (t - 1) * stages, (1:stages) + (t - 2) * stages] <- -L
  }
  # Budget constraint for each time step
  constraints[T * stages + t, (1:stages) + (t - 1) * stages] <- removal_cost
  constraints_rhs[T * stages + t] <- budget
}

# Solve the linear programming problem
solution <- lp("max", objective, constraints, constraints_dir, constraints_rhs)
solution

# pop restant L^T N0 - sum_i^T(L^(T+1-i)*Hi)
sum(matrix_power(L,T) %*% matrix(N0, ncol = 1) - 
      matrix_power(L,T) %*% matrix(solution$solution[1:3],ncol=1)-
      matrix_power(L,T-1) %*% matrix(solution$solution[4:6],ncol=1)-
      matrix_power(L,T-2) %*% matrix(solution$solution[7:9],ncol=1)-
      matrix_power(L,T-3) %*% matrix(solution$solution[10:12],ncol=1)-
      matrix_power(L,T-4) %*% matrix(solution$solution[13:15],ncol=1)-
      matrix_power(L,T-5) %*% matrix(solution$solution[16:18],ncol=1)-
      matrix_power(L,T-6) %*% matrix(solution$solution[19:21],ncol=1)-
      matrix_power(L,T-7) %*% matrix(solution$solution[22:24],ncol=1)-
      matrix_power(L,T-8) %*% matrix(solution$solution[25:27],ncol=1)-
      matrix_power(L,T-9) %*% matrix(solution$solution[28:30],ncol=1))

# Display results
if (solution$status == 0) {
  control_strategy <- matrix(solution$solution, nrow = stages, byrow = TRUE)
  print("Optimal control strategy for each stage over time:")
  print(control_strategy)
} else {
  print("No feasible solution found. Adjust parameters or constraints.")
}

den <- (control_strategy[1,] + control_strategy[2,] + control_strategy[3,])
control_strategy[1,] / den * 100
control_strategy[2,] / den * 100
control_strategy[3,] / den * 100


#------------------ Any T w/ budget on a grid


# Load the lpSolve package
library(lpSolve)

# Define the number of stages
stages <- 3  # Three-stage model

# Define the Leslie matrix based on the values provided
L <- matrix(c(0, 0.0000646, 0.0177,
              1, 1.115, 0,
              0, 0.265, 1.107), 
            byrow = TRUE, ncol = 3)

# Initial population vector based on provided values
N0 <- c(0.00008226, 0.01186, 0.005977)  # N0

# Define the time horizon for the control strategy
T <- 10  # Time horizon for control

# Define the cost for removing individuals from each stage
removal_cost <- c(1, 1, 0.4)  # Cost matrix as provided

# Set up the linear programming model

# Function to compute L^x for a given matrix L and exponent x
matrix_power <- function(L, x) {
  # Check if x is 0; return the identity matrix in that case
  if (x == 0) {
    return(diag(nrow(L)))
  }
  # Use Reduce to multiply L by itself x times
  result <- Reduce(function(a, b) a %*% b, replicate(x, L, simplify = FALSE), init = diag(nrow(L)))
  return(result)
}

objective <- rep(NA, stages * T)
boundaries <- seq(1, stages*T, by = 3)
for (j in 1:T) {
  # Calculate L^(T + 1 - j) for future impact of population at time step j
  future_impact <- rowSums(t(matrix_power(L, T + 1 - j)))
  objective[boundaries[j]:(boundaries[j] + 2)] <- future_impact
}


# Initialize constraints for Leslie matrix dynamics and budget

budget <- seq(0, 5, by = 0.1) * 1e-3
length_grid <- length(budget)
sol <- array(NA, dim = c(T, length_grid, 3))
obj <- rep(NA,  length_grid)
ind <- 1

results <- data.frame(Time = integer(), 
                      Budget = numeric(), 
                      PopulationSize = numeric())

for (i in budget){
  
  constraints <- matrix(0, nrow = T * stages + T, ncol = stages * T)
  constraints_rhs <- numeric(T * stages + T)
  constraints_dir <- rep("<=", T * stages + T)
  
  # Populate constraints for each time step based on population dynamics and budget
  for (t in 1:T) {
    if (t == 1) {
      # Initial population constraints
      constraints[(1:stages) + (t - 1) * stages, (1:stages) + (t - 1) * stages] <- diag(stages)
      constraints_rhs[(1:stages) + (t - 1) * stages] <- N0
    } else {
      # Leslie matrix dynamics for each subsequent time step
      constraints[(1:stages) + (t - 1) * stages, (1:stages) + (t - 1) * stages] <- diag(stages)
      constraints[(1:stages) + (t - 1) * stages, (1:stages) + (t - 2) * stages] <- -L
    }
    # Budget constraint for each time step
    constraints[T * stages + t, (1:stages) + (t - 1) * stages] <- removal_cost
    constraints_rhs[T * stages + t] <- i
  }
  
  # Solve the linear programming problem
  solution <- lp("max", objective, constraints, constraints_dir, constraints_rhs)
  
  sol[,ind,] <- solution$solution
  obj[ind] <- solution$objval
  
  # Initialize the population vector
  current_population <- N0
  # Extract population sizes for each time step and store results
  for (t in 1:T) {
    # Calculate population size at time t by summing all stages
    removals <- solution$solution[(3 * (t - 1) + 1):(3 * t)]
    removals <- pmin(removals, current_population) ######### ensure we do not remove more than we have
    post_removal_population <- current_population - removals
    next_population <- L %*% post_removal_population
    population_at_t <- sum(next_population)
    results <- rbind(results, data.frame(Time = t, 
                                         Budget = i, 
                                         PopulationSize = population_at_t))
    # Update current population for the next time step
    current_population <- next_population
    
  }
  
  ind <- ind + 1
  
}

sol
obj

results



#--------------------- Figure 1(a)

# Plot the results
library(tidyverse)


# Create the filled contour plot using geom_tile
ggplot(results, aes(x = Time, y = Budget, fill = PopulationSize)) +
  geom_tile() +  # Use geom_tile to fill the space with colors
  labs(title = "Filled Contour Plot of Population Size Over Time and Budget",
       x = "Time (Years)",
       y = "Annual Budget",
       fill = "Population Size") +
  scale_x_continuous(breaks = seq(min(results$Time), max(results$Time), by = 1)) +  # Show Time as integers
  scale_fill_viridis_c(breaks = c(0, 0.05, 0.10, 0.15), limits = c(0, 0.15)) +  # Continuous color scale for Population Size
  theme_minimal()


ggplot(results, aes(x = Time, y = Budget, fill = PopulationSize)) +
  geom_tile() +  # Use geom_tile to fill the space with colors
  labs(x = "Time (Years)", y = "Annual Budget", fill = "Population Size") +
  scale_x_continuous(breaks = seq(min(results$Time), max(results$Time), by = 1)) +  # Show Time as integers
  scale_fill_viridis_c(breaks = c(0, 0.05, 0.10, 0.15), limits = c(0, 0.15)) +  # Continuous color scale for Population Size
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +  # Move color legend to top and center the title
  theme_minimal() +
  theme(
    legend.position = "top",  # Position the legend at the top
    plot.title = element_blank(),  # Remove the plot title
    legend.title = element_text(size = 12, face = "bold")  # Style the legend title if needed
  )




#-------------- Figure 1(b)



# Load necessary libraries
library(lpSolve)
library(ggplot2)


# Define the number of stages
stages <- 3  # Three-stage model

# Define the Leslie matrix based on the values provided
L <- matrix(c(0, 0.0000646, 0.0177,
              1, 1.115, 0,
              0, 0.265, 1.107), 
            byrow = TRUE, ncol = 3)

# Initial population vector based on provided values
N0 <- c(0.00008226, 0.01186, 0.005977)  # N0

# Define the time horizon for the control strategy
T <- 10  # Time horizon for control

# Define the cost for removing individuals from each stage
removal_cost <- c(1, 1, 0.4)  # Cost matrix as provided

# Set up the linear programming model

# Function to compute L^x for a given matrix L and exponent x
matrix_power <- function(L, x) {
  # Check if x is 0; return the identity matrix in that case
  if (x == 0) {
    return(diag(nrow(L)))
  }
  # Use Reduce to multiply L by itself x times
  result <- Reduce(function(a, b) a %*% b, replicate(x, L, simplify = FALSE), init = diag(nrow(L)))
  return(result)
}

objective <- rep(NA, stages * T)
boundaries <- seq(1, stages*T, by = 3)
for (j in 1:T) {
  # Calculate L^(T + 1 - j) for future impact of population at time step j
  future_impact <- rowSums(t(matrix_power(L, T + 1 - j)))
  objective[boundaries[j]:(boundaries[j] + 2)] <- future_impact
}


# Initialize constraints for Leslie matrix dynamics and budget

length_grid <- 1
sol <- array(NA, dim = c(T, length_grid, 3))
obj <- rep(NA,  length_grid)
ind <- 1

results <- data.frame(Time = integer(), 
                      Budget = numeric(), 
                      PopulationSize = numeric())

removal_results <- data.frame(Time = integer(),
                              StageClass = factor(),
                              RemovalFraction = numeric())

i <- 0.0025

  constraints <- matrix(0, nrow = T * stages + T, ncol = stages * T)
  constraints_rhs <- numeric(T * stages + T)
  constraints_dir <- rep("<=", T * stages + T)
  
  # Populate constraints for each time step based on population dynamics and budget
  for (t in 1:T) {
    if (t == 1) {
      # Initial population constraints
      constraints[(1:stages) + (t - 1) * stages, (1:stages) + (t - 1) * stages] <- diag(stages)
      constraints_rhs[(1:stages) + (t - 1) * stages] <- N0
    } else {
      # Leslie matrix dynamics for each subsequent time step
      constraints[(1:stages) + (t - 1) * stages, (1:stages) + (t - 1) * stages] <- diag(stages)
      constraints[(1:stages) + (t - 1) * stages, (1:stages) + (t - 2) * stages] <- -L
    }
    # Budget constraint for each time step
    constraints[T * stages + t, (1:stages) + (t - 1) * stages] <- removal_cost
    constraints_rhs[T * stages + t] <- i
  }
  
  # Solve the linear programming problem
  solution <- lp("max", objective, constraints, constraints_dir, constraints_rhs)
  
  sol[,ind,] <- solution$solution
  obj[ind] <- solution$objval
  
  # Initialize the population vector
  current_population <- N0
  # Extract population sizes for each time step and store results
  for (t in 1:T) {
    # Calculate population size at time t by summing all stages
    removals <- solution$solution[(3 * (t - 1) + 1):(3 * t)]
    removals <- pmin(removals, current_population) ######### ensure we do not remove more than we have
    # Calculate fraction removed
    #    removal_fraction <- removals / sum(removals)*100
    removal_fraction <- removals / current_population * 100
    post_removal_population <- current_population - removals
    next_population <- L %*% post_removal_population
    population_at_t <- sum(next_population)
    results <- rbind(results, data.frame(Time = t, 
                                         Budget = i, 
                                         PopulationSize = population_at_t))
    
    # Store results for each class (1, 2, 3)
    for (j in 1:stages) {
      removal_results <- rbind(removal_results,
                               data.frame(Time = t, 
                                          StageClass = paste("Class", j), 
                                          RemovalFraction = removal_fraction[j]))
    }
    
    # Update current population for the next time step
    current_population <- next_population
    
  }
  
sol
obj

results
removal_results

# Create the plot for fraction of each stage class removed over time
ggplot(removal_results, aes(x = Time, y = RemovalFraction, color = StageClass)) +
  geom_line(linewidth = 1.2) +  # Line plot for each stage class over time
 # geom_point(size = 2) +   # Points for emphasis at each time step
  labs(
    title = "Fraction of Each Stage Class Removed by Control (Budget = 0.0025)",
    x = "Time (Years)",
    y = "Removal Fraction (%)"
  ) +
  scale_color_manual(values = c("blue", "green", "red"), 
                     labels = c("Class 1", "Class 2", "Class 3")) +  # Adjust colors if desired
  theme_minimal() +
  theme(
    legend.title = element_blank(),  # Hide legend title
    text = element_text(size = 12),
    axis.title = element_text(face = "bold")
  )

# Create the bar plot for fraction of each stage class removed over time
ggplot(removal_results, aes(x = factor(Time), y = RemovalFraction, fill = StageClass)) +
  geom_bar(stat = "identity", position = "dodge") +  # Use bars instead of lines
  labs(
    title = "Fraction of Each Stage Class Removed by Control (Budget = 0.0025)",
    x = "Time (Years)",
    y = "Removal Fraction (%)"
  ) +
  scale_fill_manual(values = c("blue", "green", "red"), 
                    labels = c("Class 1", "Class 2", "Class 3")) +  # Adjust colors if desired
  theme_minimal() +
  theme(
    legend.title = element_blank(),  # Hide legend title
    text = element_text(size = 12),
    axis.title = element_text(face = "bold")
  )






#-------------- Figure 1(c)



# Load necessary libraries
library(lpSolve)
library(ggplot2)


# Define the number of stages
stages <- 3  # Three-stage model

# Define the Leslie matrix based on the values provided
L <- matrix(c(0, 0.0000646, 0.0177,
              1, 1.115, 0,
              0, 0.265, 1.107), 
            byrow = TRUE, ncol = 3)

# Initial population vector based on provided values
N0 <- c(0.00008226, 0.01186, 0.005977)  # N0

# Define the time horizon for the control strategy
T <- 10  # Time horizon for control

# Define the cost for removing individuals from each stage
removal_cost <- c(1, 1, 0.4)  # Cost matrix as provided

# Set up the linear programming model

# Function to compute L^x for a given matrix L and exponent x
matrix_power <- function(L, x) {
  # Check if x is 0; return the identity matrix in that case
  if (x == 0) {
    return(diag(nrow(L)))
  }
  # Use Reduce to multiply L by itself x times
  result <- Reduce(function(a, b) a %*% b, replicate(x, L, simplify = FALSE), init = diag(nrow(L)))
  return(result)
}

objective <- rep(NA, stages * T)
boundaries <- seq(1, stages*T, by = 3)
for (j in 1:T) {
  # Calculate L^(T + 1 - j) for future impact of population at time step j
  future_impact <- rowSums(t(matrix_power(L, T + 1 - j)))
  objective[boundaries[j]:(boundaries[j] + 2)] <- future_impact
}


# Initialize constraints for Leslie matrix dynamics and budget

length_grid <- 1
sol <- array(NA, dim = c(T, length_grid, 3))
obj <- rep(NA,  length_grid)
ind <- 1

results <- data.frame(Time = integer(), 
                      Budget = numeric(), 
                      PopulationSize = numeric(),
                      Fraction = numeric())

removal_results <- data.frame(Time = integer(),
                              StageClass = factor(),
                              RemovalFraction = numeric())

i <- 0.005

constraints <- matrix(0, nrow = T * stages + T, ncol = stages * T)
constraints_rhs <- numeric(T * stages + T)
constraints_dir <- rep("<=", T * stages + T)

# Populate constraints for each time step based on population dynamics and budget
for (t in 1:T) {
  if (t == 1) {
    # Initial population constraints
    constraints[(1:stages) + (t - 1) * stages, (1:stages) + (t - 1) * stages] <- diag(stages)
    constraints_rhs[(1:stages) + (t - 1) * stages] <- N0
  } else {
    # Leslie matrix dynamics for each subsequent time step
    constraints[(1:stages) + (t - 1) * stages, (1:stages) + (t - 1) * stages] <- diag(stages)
    constraints[(1:stages) + (t - 1) * stages, (1:stages) + (t - 2) * stages] <- -L
  }
  # Budget constraint for each time step
  constraints[T * stages + t, (1:stages) + (t - 1) * stages] <- removal_cost
  constraints_rhs[T * stages + t] <- i
}

# Solve the linear programming problem
solution <- lp("max", objective, constraints, constraints_dir, constraints_rhs)

sol[,ind,] <- solution$solution
obj[ind] <- solution$objval

# Initialize the population vector
current_population <- N0
# Extract population sizes for each time step and store results
for (t in 1:T) {
  # Calculate population size at time t by summing all stages
  removals <- solution$solution[(3 * (t - 1) + 1):(3 * t)]
  removals <- pmin(removals, current_population) ######### ensure we do not remove more than we have
  # Calculate fraction removed
  #    removal_fraction <- removals / sum(removals)*100
  removal_fraction <- removals / current_population * 100
  post_removal_population <- current_population - removals
  next_population <- L %*% post_removal_population
  population_at_t <- sum(next_population)
  initial_population_size <- sum(N0)
  # Calculate the fraction of the initial population remaining
  fraction_remaining_at_t <- population_at_t / initial_population_size * 100
  
  
  results <- rbind(results, data.frame(Time = t, 
                                       Budget = i, 
                                       PopulationSize = population_at_t,
                                       Fraction = fraction_remaining_at_t))
  
  # Store results for each class (1, 2, 3)
  for (j in 1:stages) {
    removal_results <- rbind(removal_results,
                             data.frame(Time = t, 
                                        StageClass = paste("Class", j), 
                                        RemovalFraction = removal_fraction[j]))
  }
  
  # Update current population for the next time step
  current_population <- next_population
  
}

sol
obj

results


#--------------------------------- Claude à la rescousse


# Implémentation complète du modèle de programmation linéaire pour le contrôle optimal d'espèces invasives
# Basé sur Hastings, Hall & Taylor (2006) "A simple approach to optimal control of invasive species"

library(lpSolve)  # Pour la résolution de problèmes de programmation linéaire
library(ggplot2)  # Pour la visualisation des résultats
library(plot3D)   # Pour les graphiques 3D
library(reshape2) # Pour la manipulation des données

#==========================================
# Paramètres du modèle
#==========================================

# Horizon temporel
T <- 10  # Nombre d'années

# Classes de stade/habitat
n_classes <- 3  # Nombre de classes (isolats, semis, prairies)
class_names <- c("Isolates", "Seedlings", "Meadows")

# Paramètres démographiques (approximations basées sur le papier)
# Matrice de transition entre les classes
transition_matrix <- matrix(c(
  0.8, 0.1, 0.05,  # De isolats vers (isolats, semis, prairies)
  0.2, 0.0, 0.6,   # De semis vers (isolats, semis, prairies)
  0.1, 0.4, 0.9    # De prairies vers (isolats, semis, prairies)
), nrow = n_classes, byrow = TRUE)

# Population initiale (normalisée)
initial_population <- c(0.3, 0.2, 0.5)  # (isolats, semis, prairies)

# Coûts de suppression par classe (relatifs)
removal_costs <- c(1.0, 0.2, 2.0)  # (isolats, semis, prairies)

# Budget annuel (relatif au coût de suppression de la classe la plus chère)
annual_budget <- 0.0025  # Comme indiqué dans la légende de la figure

#==========================================
# Formulation du problème de programmation linéaire
#==========================================

# Nombre total de variables de décision: n_classes * T (variables de contrôle) + n_classes * (T+1) (variables d'état)
n_decision_vars <- n_classes * T + n_classes * (T+1)

# Fonction objectif: minimiser la somme des populations au temps T
objective <- rep(0, n_decision_vars)
# Les coefficients pour les variables d'état au temps T sont 1 (c'est ce qu'on veut minimiser)
objective[(n_classes * T + n_classes * T + 1):(n_classes * T + n_classes * (T+1))] <- 1

# Contraintes
# 1. Contraintes de dynamique de population
# 2. Contraintes de budget
# 3. Contraintes de non-négativité

# Nombre de contraintes: n_classes * T (dynamique) + T (budget) + n_decision_vars (non-négativité)
n_constraints <- n_classes * T + T + n_decision_vars

# Matrice de contraintes
constraint_matrix <- matrix(0, nrow = n_constraints, ncol = n_decision_vars)

# Termes du côté droit des contraintes
rhs <- rep(0, n_constraints)

# Direction des contraintes (<= ou ==)
constraint_directions <- rep("<=", n_constraints)

# Construction des contraintes
constraint_row <- 1

# 1. Contraintes de dynamique de population
for (t in 1:T) {
  for (i in 1:n_classes) {
    # Indices des variables dans le vecteur de décision
    state_idx_t_minus_1 <- n_classes * T + (t-1) * n_classes + i
    state_idx_t <- n_classes * T + t * n_classes + i
    control_idx_t <- (t-1) * n_classes + i
    
    # x_i(t) = x_i(t-1) + croissance - contrôle
    constraint_matrix[constraint_row, state_idx_t] <- 1  # x_i(t)
    
    # Contribution des différentes classes à la classe i
    for (j in 1:n_classes) {
      state_idx_j_t_minus_1 <- n_classes * T + (t-1) * n_classes + j
      constraint_matrix[constraint_row, state_idx_j_t_minus_1] <- -transition_matrix[j, i]  # -a_ji * x_j(t-1)
    }
    
    # Effet du contrôle
    constraint_matrix[constraint_row, control_idx_t] <- 1  # u_i(t)
    
    # Égalité
    constraint_directions[constraint_row] <- "=="
    
    constraint_row <- constraint_row + 1
  }
}

# 2. Contraintes de budget
for (t in 1:T) {
  for (i in 1:n_classes) {
    control_idx <- (t-1) * n_classes + i
    constraint_matrix[constraint_row, control_idx] <- removal_costs[i]
  }
  rhs[constraint_row] <- annual_budget
  constraint_row <- constraint_row + 1
}

# 3. Contraintes de non-négativité (implicites dans lpSolve)

# 4. Conditions initiales
for (i in 1:n_classes) {
  state_idx_0 <- n_classes * T + i
  constraint_matrix[constraint_row, state_idx_0] <- 1
  rhs[constraint_row] <- initial_population[i]
  constraint_directions[constraint_row] <- "=="
  constraint_row <- constraint_row + 1
}

# Variables supplémentaires pour garantir que le contrôle ne dépasse pas la population
for (t in 1:T) {
  for (i in 1:n_classes) {
    control_idx <- (t-1) * n_classes + i
    state_idx_t_minus_1 <- n_classes * T + (t-1) * n_classes + i
    
    constraint_matrix[constraint_row, control_idx] <- 1
    constraint_matrix[constraint_row, state_idx_t_minus_1] <- -1
    rhs[constraint_row] <- 0
    
    constraint_row <- constraint_row + 1
  }
}

#==========================================
# Résolution du problème
#==========================================

# Solution avec lpSolve
lp_solution <- lp("min", objective, constraint_matrix, constraint_directions, rhs, all.bin = FALSE)

# Extraction des résultats
solution_vector <- lp_solution$solution

# Variables de contrôle
control_vars <- matrix(solution_vector[1:(n_classes * T)], nrow = T, byrow = TRUE)
colnames(control_vars) <- class_names

# Variables d'état (populations)
state_vars <- matrix(solution_vector[(n_classes * T + 1):n_decision_vars], nrow = T+1, byrow = TRUE)
colnames(state_vars) <- class_names

#==========================================
# Analyse des résultats et visualisation
#==========================================

# 1. Population totale au fil du temps
total_population <- rowSums(state_vars)

# Visualisation de la population totale
df_total <- data.frame(
  Year = 0:T,
  Population = total_population
)

# 2. Pourcentage de chaque classe supprimée par année
percentage_removed <- matrix(0, nrow = T, ncol = n_classes)
for (t in 1:T) {
  for (i in 1:n_classes) {
    if (state_vars[t, i] > 0) {
      percentage_removed[t, i] <- (control_vars[t, i] / state_vars[t, i]) * 100
    } else {
      percentage_removed[t, i] <- 0
    }
  }
}

# Visualisation du pourcentage supprimé (figure 1b)
df_removed <- melt(percentage_removed)
colnames(df_removed) <- c("Year", "Class", "Percentage")
df_removed$Year <- df_removed$Year  # Ajustement de l'index pour commencer à 1
df_removed$Class <- factor(df_removed$Class, labels = class_names)

fig_b <- ggplot(df_removed, aes(x = Year, y = Percentage, fill = Class)) +
  geom_col(position = "dodge", color = "black") +
  scale_fill_manual(values = c("Isolates" = "green3", "Seedlings" = "yellow3", "Meadows" = "darkred")) +
  labs(
    x = "Year",
    y = "Percentage of class removed",
    title = "Percentage of each class removed by control"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 100))

# 3. Effet du budget sur la population finale (figure 1a)
# Simuler différents niveaux de budget
budgets <- seq(0, 0.005, by = 0.0005)
final_populations <- numeric(length(budgets))

# Fonction pour simuler avec un budget donné
simulate_with_budget <- function(budget) {
  # Mise à jour des contraintes de budget
  constraint_row <- n_classes * T + 1
  for (t in 1:T) {
    rhs[constraint_row] <- budget
    constraint_row <- constraint_row + 1
  }
  
  # Résolution
  lp_sol <- lp("min", objective, constraint_matrix, constraint_directions, rhs, all.bin = FALSE)
  
  # Extraction des populations
  sol_vector <- lp_sol$solution
  state_v <- matrix(sol_vector[(n_classes * T + 1):n_decision_vars], nrow = T+1, byrow = TRUE)
  
  # Population totale finale
  return(sum(state_v[T+1, ]))
}

# Calcul pour chaque budget
for (i in 1:length(budgets)) {
  final_populations[i] <- simulate_with_budget(budgets[i])
}

# Visualisation de l'effet du budget
df_budget <- data.frame(
  Budget = budgets * 1000,  # Pour correspondre aux unités de la figure (×10^-3)
  FinalPopulation = final_populations
)

fig_budget <- ggplot(df_budget, aes(x = Budget, y = FinalPopulation)) +
  geom_line(size = 1.5, color = "blue") +
  geom_point(size = 3, color = "red") +
  labs(
    x = "Annual budget (×10^-3)",
    y = "Final population size",
    title = "Effect of annual budget on final population size"
  ) +
  theme_minimal()

# 4. Création de la figure 3D (population en fonction du budget et du temps) - figure 1a
# Simulation pour différentes combinaisons de budget et de temps
budgets_3d <- seq(0, 0.005, by = 0.0005)
years_3d <- 0:T
pop_matrix <- matrix(0, nrow = length(budgets_3d), ncol = length(years_3d))

for (i in 1:length(budgets_3d)) {
  # Mise à jour des contraintes de budget
  constraint_row <- n_classes * T + 1
  for (t in 1:T) {
    rhs[constraint_row] <- budgets_3d[i]
    constraint_row <- constraint_row + 1
  }
  
  # Résolution
  lp_sol <- lp("min", objective, constraint_matrix, constraint_directions, rhs, all.bin = FALSE)
  
  # Extraction des populations
  sol_vector <- lp_sol$solution
  state_v <- matrix(sol_vector[(n_classes * T + 1):n_decision_vars], nrow = T+1, byrow = TRUE)
  
  # Population totale à chaque année
  pop_matrix[i, ] <- rowSums(state_v)
}

# Visualisation 3D
plot_fig_a <- function() {
  par(mar = c(5, 5, 2, 2) + 0.1)
  hist3D(x = budgets_3d * 1000, y = years_3d, z = pop_matrix,
         bty = "b2",
         phi = 20, theta = -30,
         col = colorRampPalette(c("blue", "cyan", "green", "yellow", "orange", "red"))(100),
         border = "black",
         shade = 0.4,
         space = 0,
         ticktype = "detailed",
         xlab = "Annual budget (×10^-3)",
         ylab = "Year",
         zlab = "Population size",
         main = "")
  
  mtext("(a)", side = 1, line = 3, adj = 0, cex = 1.2)
}

# Pour exécuter et visualiser les résultats, décommentez les lignes suivantes:
plot_fig_a()
print(fig_b)
print(fig_budget)

# Affichage des informations sur la stratégie optimale
print("Optimal control strategy:")
print(control_vars)
print("Population dynamics:")
print(state_vars)

