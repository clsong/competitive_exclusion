library(here)
library(tidyverse)
library(mvtnorm)
library(DiagrammeR)
library(RColorBrewer)
library(latex2exp)
library(viridis)
library(tictoc)
library(furrr)
library(purrrgress)
library(progressr)
library(ggridges)

get_alpha <- function(density_form) {
  # a; species interactions ------------------------------------------------------------
  get_competition_coefficient <- function(species_name) {
    here("data", paste0(species_name, ".txt")) %>%
      read_table2(col_names = FALSE) %>%
      mutate(X1 = X1 * (1 - X2) / X2) %>% # mean of negative binomial distribution
      select(-X2) %>%
      magrittr::set_colnames(c("lambda", "SP", "EG", "BH", "BD", "AB", "leaf")) %>%
      mutate(sample = row_number()) %>%
      select(-leaf, -lambda) %>%
      gather(species_j, interaction, -sample) %>%
      mutate(species_i = species_name) %>%
      select(sample, species_i, species_j, interaction)
  }

  alpha <- c("AB", "BD", "BH", "EG", "SP") %>%
    map_dfr(~ get_competition_coefficient(.x))

  # lambda; intrinsic growth rates ------------------------------------------------------------
  get_lambda <- function(species_name) {
    here("data", paste0(species_name, ".txt")) %>%
      read_table2(col_names = FALSE) %>%
      mutate(X1 = X1 * (1 - X2) / X2) %>% # mean of negative binomial distribution
      select(X1) %>%
      rename(lambda = X1) %>%
      mutate(
        sample = row_number(),
        species_j = species_name
      )
  }

  lambda <- c("AB", "BD", "BH", "EG", "SP") %>%
    map_dfr(~ get_lambda(.x))

  # omega; adult survival ----------------------------------------------------------
  EGAdultSurv <- read_csv("data/EGAdultSurv.txt", col_names = FALSE) %>%
    rename(omega = X1) %>%
    mutate(
      species_j = "EG",
      sample = row_number()
    )

  SPAdultSurv <- read_csv("data/SPAdultSurv.txt", col_names = FALSE) %>%
    rename(omega = X1) %>%
    mutate(
      species_j = "SP",
      sample = row_number()
    )

  AdultSurv <- bind_rows(EGAdultSurv, SPAdultSurv)

  # g; germination -------------------------------------------------------------
  germMat <- read_table2("data/germMat.txt", col_types = cols(X6 = col_skip())) %>%
    gather(species_j, g)

  # v; transition -------------------------------------------------------------------
  perTrans <- read_table2("data/perTrans.txt", col_names = FALSE) %>%
    select(v = X4) %>%
    mutate(
      sample = row_number()
    )

  # all ---------------------------------------------------------------------
  alpha %>%
    left_join(germMat, by = "species_j") %>%
    left_join(lambda, by = c("sample", "species_j")) %>%
    left_join(AdultSurv, by = c("sample", "species_j")) %>%
    left_join(perTrans, by = c("sample")) %>%
    mutate(type = if_else(species_j %in% c("AB", "BD", "BH"), "Annual", "Perennial")) %>%
    mutate(v = ifelse(type == "Perennial", v, NA)) %>%
    mutate(alpha_effective = case_when(
      type == "Annual" ~ interaction * g,
      type == "Perennial" & density_form == 0 ~ interaction * g,
      type == "Perennial" & density_form == 1 ~ interaction * g * (1 + v / (1 - omega)),
      type == "Perennial" & density_form == 2 ~ interaction * g * (1 + sqrt(v / (lambda * (1 - omega))))
      # stage * value * (g + g * v / (1 - omega))) + (1-stage) * value * g
    )) %>%
    mutate(density_form = density_form) %>%
    select(sample, species_i, species_j, alpha_effective, interaction, density_form)
}

get_r <- function(density_form) {
  # lambda; species interactions ------------------------------------------------------------
  get_lamabda <- function(species_name) {
    here("data", paste0(species_name, ".txt")) %>%
      read_table2(col_names = FALSE) %>%
      mutate(X1 = X1 * (1 - X2) / X2) %>%
      select(X1) %>%
      rename(lambda = X1) %>%
      mutate(
        sample = row_number(),
        species_i = species_name
      )
  }

  r <- c("AB", "BD", "BH", "EG", "SP") %>%
    map_dfr(~ get_lamabda(.x))

  # omega; adult survival ----------------------------------------------------------
  EGAdultSurv <- read_csv("data/EGAdultSurv.txt", col_names = FALSE) %>%
    rename(omega = X1) %>%
    mutate(
      species_i = "EG",
      sample = row_number()
    )

  SPAdultSurv <- read_csv("data/SPAdultSurv.txt", col_names = FALSE) %>%
    rename(omega = X1) %>%
    mutate(
      species_i = "SP",
      sample = row_number()
    )

  AdultSurv <- bind_rows(EGAdultSurv, SPAdultSurv)

  # g; germination -------------------------------------------------------------
  germMat <- read_table2("data/germMat.txt", col_types = cols(X6 = col_skip())) %>%
    gather(species_i, g)

  # v; transition -------------------------------------------------------------------
  perTrans <-
    read_table2("data/perTrans.txt", col_names = FALSE) %>%
    select(v = X4) %>%
    mutate(
      sample = row_number()
    )


  # combine data ------------------------------------------------------------
  r %>%
    left_join(germMat, by = "species_i") %>%
    left_join(AdultSurv, by = c("sample", "species_i")) %>%
    left_join(perTrans, by = c("sample")) %>%
    mutate(type = if_else(species_i %in% c("AB", "BD", "BH"), "Annual", "Perennial")) %>%
    distinct() %>%
    mutate(r_effective = case_when(
      type == "Annual" ~ lambda - 1,
      type == "Perennial" & density_form == 0 ~ lambda - 1,
      type == "Perennial" & density_form == 1 ~ lambda * v / (1 - omega) - 1,
      type == "Perennial" & density_form == 2 ~ sqrt(lambda * v / (1 - omega)) - 1
    )) %>%
    mutate(v = ifelse(type == "Perennial", v, NA)) %>%
    mutate(density_form = density_form)
}

calculate_Omega <- function(alpha) {
  S <- nrow(alpha)
  omega <- function(S, Sigma) {
    m <- matrix(0, S, 1)
    a <- matrix(0, S, 1)
    b <- matrix(Inf, S, 1)
    d <- pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
    out <- d[1]^(1 / S)
    return(out)
  }
  f <- function(m) class(try(solve(t(m) %*% m), silent = T)) == "matrix"
  if (f(alpha) == FALSE) {
    return(0)
  }
  else {
    Sigma <- solve(t(alpha) %*% alpha)
    return(omega(S, Sigma))
  }
}

determine_omega_feasibility <- function(data, data_type = "empirical") {
  sample_points <- function(S, steps = 50) {
    fixed <- 1

    points <- matrix(0, nrow = 1, ncol = S)
    colnames(points) <- paste("a_", 1:S)

    for (i in c(1:((steps - 2) / (S - 2)))) {
      fixed <- i
      up <- c(1:(steps - fixed * (S - 2) - 1))
      down <- c(steps - fixed * (S - 2) - up)
      len_row <- (length(up))

      make_points <- cbind(matrix(up, ncol = 1, nrow = len_row), matrix(down, ncol = 1, nrow = len_row), matrix(fixed, ncol = (S - 2), nrow = len_row))
      points <- rbind(points, make_points)
    }

    points <- points / steps

    initial_points <- nrow(points)

    orders <- combinat::permn(1:S)
    orders <- orders[-1]

    for (val in orders) {
      new_points <- points[1:initial_points, val]
      points <- rbind(points, new_points)
    }

    points <- unique(points,
      incomparables = FALSE, MARGIN = 1,
      fromLast = FALSE
    )
    points <- points[-1, ]
    return(points)
  }

  check_feasibility <- function(alpha, r) {
    abundance <- solve(alpha, r)
    if (sum(abundance < 0) == 0) {
      return(1)
    } else {
      return(0)
    }
  }

  check_mean_abundance <- function(alpha, r) {
    solve(alpha, r) %>%
      mean()
  }

  check_stability <- function(alpha, r) {
    if (length(r) > 1) {
      abundance <- solve(alpha, r)
      if (sum(abundance < 0) == 0) {
        jacobian <- -diag(abundance) %*% alpha
        eigenvalues <- eigen(jacobian)$values
        if (sum(Re(eigenvalues) > 0) == 0) {
          return(1)
        } else {
          return(0)
        }
      } else {
        return(0)
      }
    } else {
      return(1)
    }
  }

  check_stability_N <- function(alpha, N) {
    jacobian <- -diag(N) %*% alpha
    eigenvalues <- eigen(jacobian)$values
    if_else(sum(Re(eigenvalues) > 0) == 0, 1, 0)
  }

  check_stability_probability <- function(alpha, steps = 10) {
    feasibility <- c()

    if (nrow(alpha) > 2) {
      points <- sample_points(nrow(alpha), steps)
    } else {
      points_1 <- (1:(steps - 1)) / steps
      points_2 <- 1 - points_1
      points <- tibble(points_1, points_2)
    }

    points %>%
      t() %>%
      as_tibble() %>%
      map_dbl(~ check_stability_N(alpha, .)) %>%
      mean()
  }

  get_result <- function(species_sample) {
    alpha <- alpha_all[species_sample, species_sample]
    r <- r_all[species_sample]
    tibble(
      omega = calculate_Omega(alpha),
      feasibility = check_feasibility(alpha, r),
      stability = check_stability(alpha, r),
      species_sample = paste(as.character(names_all[species_sample]), collapse = ","),
    )
  }

  Nspp <- length(unique(data$species_i))
  names_all <- c("AB", "BD", "BH", "EG", "SP")

  alpha_all <- data %>%
    select(species_i, species_j, alpha_effective) %>%
    pivot_wider(names_from = species_j, values_from = alpha_effective) %>%
    select(AB, BD, BH, EG, SP) %>%
    as.matrix()

  r_all <- data %>%
    pull(r_effective) %>%
    unique()

  bind_rows(
    combn(1:Nspp, 2) %>%
      as_tibble() %>%
      map_dfr(~ get_result(.)),
    combn(1:Nspp, 3) %>%
      as_tibble() %>%
      map_dfr(~ get_result(.)),
    combn(1:Nspp, 4) %>%
      as_tibble() %>%
      map_dfr(~ get_result(.)),
    combn(1:Nspp, 5) %>%
      as_tibble() %>%
      map_dfr(~ get_result(.))
  )
}

get_dynamics_from_simulation_form0 <- function(data, threshold) {
  species <- data$species_i %>%
    unique()
  species_num <- length(species)
  arrival_orders <- gtools::permutations(n = species_num, r = species_num, v = 1:species_num, repeats.allowed = F)
  data_single <- data %>%
    select(species_i, lambda, g, omega, v, type) %>%
    distinct()
  data_pairwise <- data %>%
    select(species_i, species_j, interaction) %>%
    group_split(species_i) %>%
    map(~ arrange(., species_j)) %>%
    map(~ pull(., interaction))

  Nsimu <- 1e3
  initial_abundance <- 10
  total_time <- species_num * Nsimu

  survived <- 1:nrow(arrival_orders) %>%
    map(function(arrival_order_label) {
      arrival_order <- arrival_orders[arrival_order_label, ]
      N <- matrix(NA, nrow = total_time, ncol = species_num)
      N[1, arrival_order[1]] <- initial_abundance
      N[1, arrival_order[2:species_num]] <- 0
      for (time in 1:(total_time - 1)) {
        for (i in 1:species_num) {
          if (N[time, i] < 1e-10) {
            N[time + 1, i] <- 0
          } else {
            N[time + 1, i] <- N[time, i] * data_single$g[i] * data_single$lambda[i] / (1 + sum(data_pairwise[[i]] * data_single$g * N[time, ])) +
              N[time, i] * (1 - data_single$g[i])
          }
        }
        if (time %% Nsimu == 0) {
          N[time + 1, arrival_order[time %/% Nsimu + 1]] <- initial_abundance
        }
      }

      species_survived <- species[which(N[nrow(N), ] > threshold)]
      species_survived
    })

  survived %>%
    unlist() %>%
    table() %>%
    as_tibble() %>%
    rename(species_label = ".") %>%
    right_join(
      tibble(species_label = species),
      by = "species_label"
    ) %>%
    mutate(dynamics = case_when(
      n == nrow(arrival_orders) ~ "Persist",
      n < nrow(arrival_orders) & n > 0 ~ "Contingently excluded",
      is.na(n) ~ "Deteministically excluded"
    ))
}

get_dynamics_from_simulation_form2 <- function(data, threshold) {
  species <- data$species_i %>%
    unique()
  species_num <- length(species)
  species_info <- data %>%
    select(species_i, type) %>%
    distinct() %>%
    mutate(num = if_else(type == "Annual", 1, 2))
  perennial_num <- sum((species_info$type == "Perennial"))
  annual_species <- which(species_info$type == "Annual")
  species_total_num <- sum(species_info$num)
  arrival_orders <- gtools::permutations(n = species_num, r = species_num, v = 1:species_num, repeats.allowed = F)
  data_single <- data %>%
    select(species_i, lambda, g, omega, v, type) %>%
    distinct()
  data_pairwise <- data %>%
    select(species_i, species_j, interaction) %>%
    group_split(species_i) %>%
    map(~ arrange(., species_j)) %>%
    map(~ pull(., interaction))

  Nsimu <- 1e3
  total_time <- species_num * Nsimu
  initial_abundance <- 10

  survived <- 1:nrow(arrival_orders) %>%
    map(function(arrival_order_label) {
      arrival_order <- arrival_orders[arrival_order_label, ]
      N <- matrix(NA, nrow = total_time, ncol = species_total_num)
      N[1, ] <- 0
      if (species_info[arrival_order[1], ]$type == "Annual") {
        N[1, arrival_order[1]] <- initial_abundance
      } else {
        N[1, arrival_order[1]] <- initial_abundance
        N[1, arrival_order[1] + perennial_num] <- initial_abundance
      }
      for (time in 1:(total_time - 1)) {
        for (i in 1:species_num) {
          if (N[time, i] < 1e-10) {
            if (i %in% annual_species) {
              N[time + 1, i] <- 0
            } else {
              N[time + 1, i] <- 0
              N[time + 1, i + perennial_num] <- 0
            }
          } else {
            if (i %in% annual_species) {
              competition_total <- sum(data_pairwise[[i]] * data_single$g * N[time, 1:species_num]) +
                sum(data_pairwise[[i]][(species_num + 1 - perennial_num):(species_total_num - perennial_num)] * N[time, (species_num + 1):species_total_num])
              N[time + 1, i] <- N[time, i] * data_single$g[i] * data_single$lambda[i] / (1 + competition_total) +
                N[time, i] * (1 - data_single$g[i])
            } else {
              competition_total <- sum(data_pairwise[[i]] * data_single$g * N[time, 1:species_num]) +
                sum(data_pairwise[[i]][(species_num + 1 - perennial_num):(species_total_num - perennial_num)] * N[time, (species_num + 1):species_total_num])
              N[time + 1, i] <- N[time, i + perennial_num] * data_single$lambda[i] / (1 + competition_total) +
                N[time, i] * (1 - data_single$g[i])
              N[time + 1, i + perennial_num] <- N[time, i + perennial_num] * data_single$omega[i] +
                N[time, i] * data_single$g[i] * data_single$v[i] / (1 + competition_total)
            }
          }
        }
        if (time %% Nsimu == 0) {
          arrival_species <- arrival_order[time %/% Nsimu + 1]
          if (species_info[arrival_species, ]$type == "Annual") {
            N[time + 1, arrival_species] <- initial_abundance
          } else {
            N[time + 1, arrival_species] <- initial_abundance
            N[time + 1, arrival_species + perennial_num] <- initial_abundance
          }
        }
      }

      species_survived <- species[which(N[nrow(N), 1:species_num] > threshold)]

      species_survived
    })
  
  survived %>%
    unlist() %>%
    table() %>%
    as_tibble() %>%
    rename(species_label = ".") %>%
    right_join(
      tibble(species_label = species),
      by = "species_label"
    ) %>%
    mutate(dynamics = case_when(
      n == nrow(arrival_orders) ~ "Persist",
      n < nrow(arrival_orders) & n > 0 ~ "Contingently excluded",
      is.na(n) ~ "Deteministically excluded"
    ))
}

get_dynamics_from_simulation <- function(data, density_form, threshold = 1e-5) {
  if (density_form == 0) {
    return(get_dynamics_from_simulation_form0(data, threshold))
  } else if (density_form == 2) {
    return(get_dynamics_from_simulation_form2(data, threshold))
  }
}