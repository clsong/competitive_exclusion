source("toolbox.r")

# parallel computing
plan(multisession, workers = 4)

# load original data
df_ori <- bind_rows(
  left_join(get_alpha(0), get_r(0), by = c("sample", "species_i", "density_form")) %>%
    nest(data = -sample) %>%
    mutate(density_form = 0),
  left_join(get_alpha(2), get_r(2), by = c("sample", "species_i", "density_form")) %>%
    nest(data = -sample) %>%
    mutate(density_form = 2)
)

# add dynamics predicted by the structural stability

df_SS <- df_ori %>%
  mutate(dynamics = future_map(data, function(x) {
    x %>%
      select(species_i, species_j, alpha_effective, density_form, r_effective) %>%
      determine_omega_feasibility()
  }, .progress = T)) %>%
  unnest(dynamics) %>%
  mutate(species_names = pro_map(species_sample, ~ unlist(str_split(., ",")))) %>%
  rowwise() %>%
  mutate(species_comb_type = case_when(
    any(str_detect(species_sample, c("AB", "BD", "BH"))) & any(str_detect(species_sample, c("SP", "EG"))) ~ "annual-perennial",
    any(str_detect(species_sample, c("AB", "BD", "BH"))) & !any(str_detect(species_sample, c("SP", "EG"))) ~ "annual-annual",
    TRUE ~ "perennial-perennial"
  )) %>%
  ungroup() %>%
  mutate(num_species = (str_length(species_sample) + 1) / 3)


# add dynamics by the simulation
df_SS_simulation <- df_SS %>%
  filter(species_comb_type == "annual-perennial") %>%
  # filter(density_form == 0) %>%
  # sample_n(1000) %>%
  mutate(data = pro_map2(data, species_names, function(x, y) {
    x %>%
      filter(species_i %in% y & species_j %in% y)
  })) %>%
  mutate(dynamic_pattern = future_map2(data, density_form, ~get_dynamics_from_simulation(.x, .y), .progress = T))