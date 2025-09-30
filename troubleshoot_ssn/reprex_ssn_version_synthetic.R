library(data.table)
library(dplyr)
library(igraph)
library(magrittr)
library(MASS)
library(mvtnorm)
library(OCNet)
library(sf)
library(SSN2)
library(SSNbler)

overwrite <- TRUE
out_dir <- tempdir()  # directory to save LSNs and SSNs
scale_factor <- 500
lsn_path <- file.path(out_dir, "OCNs_lsn")
dir.create(lsn_path, showWarnings = FALSE)

set.seed(123)

# ---------------------- convert_OCN ---------------------- #
OCN2graph <- function(OCN, idcol) {
  graph = graph_from_adjacency_matrix(as.matrix(OCN$RN$W)) %>%
    set_vertex_attr(name ='x', value = OCN$RN$X) %>%
    set_vertex_attr(name ='y', value = OCN$RN$Y) %>%
    set_vertex_attr(name ='DA', value = OCN$RN$A) %>%
    set_vertex_attr(name = idcol, value = 1:gorder(.))
  E(graph)$weight <- OCN$RN$leng[-7]
  graph = set_vertex_attr(graph, "width", value = OCN$RN$A)
  return(graph)
}

convert_OCN <- function(OCN, out_format, out_SSNdir, idcol = 'patch') {
  out_OCN_formatted <- list()
  if ('igraph' %in% out_format) {
    igraph <- OCN2graph(OCN, idcol)
    out_OCN_formatted <- c(out_OCN_formatted, igraph=list(igraph))
  }
  if ('SSN' %in% out_format) {
    SSN <- OCN_to_SSN(OCN = OCN, level = 'RN', obsSites = 1:OCN$RN$nNodes,
                      path = out_SSNdir, randomAllocation = FALSE, importToR = TRUE)
    ssn_df <- SSN2::ssn_get_data(SSN)
    ssn_df[, idcol] <- ssn_df[, 'locID']
    SSN <- SSN::ssn_put_data(ssn_df, x = SSN)
    out_OCN_formatted <- c(out_OCN_formatted, SSN=SSN)
  }
  return(out_OCN_formatted)
}

# ---------------------- generate_OCN_formatted ---------------------- #
generate_OCN_formatted <- function(patches, out_format, out_SSNdir,
                                   cellsize = 0.5, dimX = 25, dimY = 25,
                                   outletPos = 3, expEnergy = 0.1, 
                                   coolingRate = 0.3, slope0 = 0.05) {
  OCN <- create_OCN(dimX, dimY, cellsize = cellsize, outletPos = outletPos,
                    expEnergy = expEnergy, coolingRate = coolingRate)
  OCN <- landscape_OCN(OCN, slope0 = slope0)
  thrA <- as.data.table(OCNet::find_area_threshold_OCN(OCN = OCN, thrValues = seq(0.5,50,0.05))) %>%
    .[which.min(abs(nNodesRN - patches)), thrValues]
  OCN <- aggregate_OCN(OCN, thrA = thrA)
  graph <- convert_OCN(OCN, out_format = 'igraph')
  return(list(ocn = OCN, igraph = graph$igraph))
}

# ---------------------- generate_response (correlated campaigns) ---------------------- #
generate_response <- function(g, range_par = 3, beta = 1.2, nb_size = 40, n_campaigns = 2, rho = 0.5, sigma_campaign = 0.5) {
  n <- vcount(g)
  
  # 1. Network covariance
  D <- distances(g, mode = 'all', weights = E(g)$weight)
  cov_mat <- exp(-D / range_par)
  diag(cov_mat) <- 1
  
  # 2. Base latent spatial effect (shared across campaigns)
  latent_base <- as.numeric(rmvnorm(1, mean = rep(0, n), sigma = cov_mat))
  
  # 3. Campaign-specific deviations
  # Correlated across campaigns? (rho controls correlation)
  corr_campaign <- matrix(rho, nrow = n_campaigns, ncol = n_campaigns)
  diag(corr_campaign) <- 1
  L <- chol(corr_campaign)
  latent_campaign_raw <- matrix(rnorm(n * n_campaigns, 0, sigma_campaign), nrow = n, ncol = n_campaigns)
  latent_campaign <- latent_campaign_raw %*% t(L)   # n x n_campaigns
  
  # 4. Combine base + campaign deviations
  campaigns <- lapply(1:n_campaigns, function(i) {
    log_mu <- beta * log(V(g)$DA + 1) + latent_base + latent_campaign[,i]
    mu <- exp(log_mu)
    richness <- rnbinom(n, mu = mu, size = nb_size)
    data.frame(site = 1:n, basin_area_km2 = V(g)$DA, 
               latent_effect = latent_base + latent_campaign[,i],
               richness = richness, campaign = i)
  })
  
  do.call(rbind, campaigns)
}



# ---------------------- Shift and prune ---------------------- #
shift_ocn <- function(ocn, x_shift = 0, y_shift = 0) {
  g <- ocn$igraph
  V(g)$x <- V(g)$x + x_shift
  V(g)$y <- V(g)$y + y_shift
  ocn$igraph <- g
  return(ocn)
}

prune_complex_confluences <- function(g) {
  repeat {
    in_deg <- degree(g, mode = 'in')
    complex_nodes <- which(in_deg > 2)
    if (length(complex_nodes) == 0) break
    for (v in complex_nodes) {
      in_edges <- incident(g, v, mode = 'in')
      edges_to_remove <- sample(in_edges, length(in_edges) - 2)
      g <- delete_edges(g, edges_to_remove)
    }
  }
  g
}

# ---------------------- Generate OCNs and responses ---------------------- #
ocn1 <- generate_OCN_formatted(patches = 100, out_format = 'igraph')
ocn1$igraph <- prune_complex_confluences(ocn1$igraph)
response_1 <- generate_response(ocn1$igraph, n_campaigns = 2)

ocn2 <- generate_OCN_formatted(patches = 100, out_format = 'igraph') %>% shift_ocn(x_shift = 500000, y_shift = 0)
ocn2$igraph <- prune_complex_confluences(ocn2$igraph)
response_2 <- generate_response(ocn2$igraph, n_campaigns = 2)

# ---------------------- Prepare edges ---------------------- #
ocn_to_edges_sf <- function(ocn, country_name, scale_factor = 1, crs_proj = 3035) {
  g <- ocn$igraph
  edges_df <- data.frame(from = as.integer(ends(g, E(g))[,1]), to = as.integer(ends(g, E(g))[,2]))
  geom_list <- lapply(1:nrow(edges_df), function(i)
    st_linestring(scale_factor * matrix(c(V(g)$x[edges_df$from[i]], V(g)$y[edges_df$from[i]],
                                          V(g)$x[edges_df$to[i]], V(g)$y[edges_df$to[i]]), ncol = 2, byrow = TRUE)))
  st_sf(cat = 1:nrow(edges_df), country = country_name, qsim_avg = 1, geometry = st_sfc(geom_list, crs = crs_proj))
}

edges_ocn1 <- ocn_to_edges_sf(ocn1, 'OCN1', scale_factor)
edges_ocn2 <- ocn_to_edges_sf(ocn2, 'OCN2', scale_factor)
net_eu <- rbind(edges_ocn1, edges_ocn2)

# ---------------------- Prepare sites with campaigns ---------------------- #
ocn_to_sites_sf <- function(ocn, country_name, scale_factor = 1, crs_proj = 3035, response) {
  df <- data.frame(site = response$site,
                   country = country_name,
                   richness = response$richness,
                   basin_area_km2 = response$basin_area_km2,
                   x_coords = V(ocn$igraph)$x[response$site] * scale_factor,
                   y_coords = V(ocn$igraph)$y[response$site] * scale_factor,
                   campaign = response$campaign)
  st_as_sf(df, coords = c('x_coords', 'y_coords'), crs = crs_proj)
}

sites_ocn1 <- ocn_to_sites_sf(ocn1, 'OCN1', scale_factor, response = response_1)
sites_ocn2 <- ocn_to_sites_sf(ocn2, 'OCN2', scale_factor, response = response_2)
sites_eu <- rbind(sites_ocn1, sites_ocn2)

# ---------------------- Build LSN and update distances ---------------------- #
edges_lsn <- SSNbler::lines_to_lsn(
  streams = net_eu, 
  lsn_path = lsn_path, 
  check_topology = TRUE, 
  snap_tolerance = 0.1, 
  topo_tolerance = 20, 
  overwrite = overwrite)

edges_lsn <- SSNbler::updist_edges(edges = edges_lsn, save_local = TRUE, lsn_path = lsn_path, calc_length = TRUE)

sites_lsn <- SSNbler::sites_to_lsn(sites = sites_eu, edges = edges_lsn, lsn_path = lsn_path,
                                   file_name = 'sites', snap_tolerance = 5, save_local = TRUE, overwrite = overwrite)

sites_list_lsn <- list(sites = sites_lsn)
sites_list_lsn <- SSNbler::updist_sites(sites_list_lsn, edges = edges_lsn,
                                        length_col = 'Length', save_local = TRUE, lsn_path = lsn_path)

# ---------------------- Compute AFVs ---------------------- #
edges_lsn$qsim_avg_sqrt <- sqrt(100 * edges_lsn$qsim_avg)
edges_lsn <- SSNbler::afv_edges(edges_lsn, infl_col = 'qsim_avg_sqrt', segpi_col = 'pi_qsqrt', afv_col = 'afv_qsqrt', lsn_path = lsn_path)
sites_list_lsn <- SSNbler::afv_sites(sites_list_lsn, edges_lsn, afv_col = 'afv_qsqrt', save_local = TRUE, lsn_path = lsn_path)

# ---------------------- Assemble SSN ---------------------- #
out_ssn_path <- file.path(out_dir, 'OCNs.ssn')
ssn <- SSNbler::ssn_assemble(edges_lsn, lsn_path, obs_sites = sites_list_lsn$sites, preds_list = NULL, ssn_path = out_ssn_path,
                             import = TRUE, check = TRUE, afv_col = 'afv_qsqrt', overwrite = TRUE)

# ---------------------- Fit SSN model with campaign factor ---------------------- #
SSN2::ssn_create_distmat(ssn)
lm_fit <- ssn_glm(
  formula = sqrt(richness) ~ log10(basin_area_km2), #+ country:log10(basin_area_km2),
  ssn.object = ssn, 
  family = 'gaussian',
  euclid_type = 'exponential', 
  estmethod = 'ml',
  additive = 'afv_qsqrt'
  # ,partition_factor = ~factor(campaign)
)
summary(lm_fit)

