# Server logic for the "Seasonal paths" dashboard page.


# ---------------------------------------------------------------------------
# Seasonal paths: sidebar choices follow the currently loaded dataset
# (same raw_data() reactive used by the other pages)
# ---------------------------------------------------------------------------

observeEvent(raw_data(), {
  x <- raw_data()
  ids4 <- as.character(sort(unique(x$id)))
  updateSelectInput(session, "id4",
    choices  = c("ALL", ids4),
    selected = ids4[1]
  )
  mig_vals4 <- sort(unique(x$migration[!is.na(x$migration)]))
  updateSelectInput(session, "phase4",
    choices  = mig_vals4,
    selected = mig_vals4[1]
  )
  updateSliderInput(session, "year_range4",
    min   = min(x$year),
    max   = max(x$year),
    value = c(min(x$year), max(x$year))
  )
})

# ---------------------------------------------------------------------------
# Filtered subset used for both the map layers and the corridor estimate
# ---------------------------------------------------------------------------

path_data4 <- reactive({
  req(raw_data(), input$id4, input$phase4, input$year_range4)
  x <- raw_data() |>
    filter(
      migration == input$phase4,
      year >= input$year_range4[1], year <= input$year_range4[2]
    )
  if (input$id4 != "ALL") x <- x |> filter(id == input$id4)
  x
})

# --- Line buffer: buffer the ordered track (per id x year) and union -------
# Each id x year segment is buffered on its own and the results are unioned
# together, so this works whether d holds one individual or several.
corridor_buffer4 <- function(d) {
  d |>
    arrange(timestamp) |>
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
    st_transform(32610) |>
    group_by(id, year) |>
    filter(n() >= 2) |>
    summarize(do_union = FALSE, .groups = "drop") |>
    st_cast("LINESTRING") |>
    st_buffer(input$buffer4) |>
    st_union() |>
    st_sf(geometry = _) |>
    st_transform(4326)
}

# Estimate one individual x year chunk's corridor via the Line buffer
# method only (Kriging is handled separately by the workflow functions
# below, since it operates on a shared raster grid rather than per-chunk).
estimate_corridor_one4 <- function(d_one) {
  corridor_buffer4(d_one)
}

# One individual across the selected year range, Line buffer method.
# Returns list(combined = <sf, unioned across years>, by_year = <named
# list of each year's own estimate, or NULL when only one year is
# selected>, years = <integer vector of the year(s) actually used for this
# individual>, by_year_n = <named list of each year's relocation count, or
# NULL when only one year is selected>, n = <total relocation count for
# this individual>).
estimate_corridor_individual4 <- function(d_ind) {
  yrs <- sort(unique(d_ind$year))
  if (length(yrs) <= 1) {
    return(list(combined = estimate_corridor_one4(d_ind), by_year = NULL,
                years = yrs, by_year_n = NULL, n = nrow(d_ind)))
  }

  by_year   <- list()
  by_year_n <- list()
  for (yr in yrs) {
    d_yr   <- d_ind |> filter(year == yr)
    res_yr <- tryCatch(estimate_corridor_one4(d_yr), error = function(e) NULL)
    if (!is.null(res_yr)) {
      by_year[[as.character(yr)]]   <- res_yr
      by_year_n[[as.character(yr)]] <- nrow(d_yr)
    }
  }
  req(length(by_year) > 0)
  list(combined = union_by_level(by_year), by_year = by_year,
       years = as.integer(names(by_year)), by_year_n = by_year_n,
       n = sum(unlist(by_year_n)))
}

# =============================================================================
# Kriging (CTMM) migration-corridor workflow, following corridor_workflow.md
# (Sawyer et al. 2009; Horne et al. 2007; Merkle et al. 2023; McKee et al.
# 2024; Montana FWP pronghorn report; Wyoming Migration Initiative Migration
# Mapper methodology):
#
#  1. Fit a BBMM/CTMM occurrence distribution separately for each individual
#     x year x migration event, all on one shared spatial grid so outputs
#     are directly comparable, excluding gaps of >= 8h between fixes from
#     the bridging (amt::od() doesn't expose this, so occurrence() is
#     called directly here — see bbmm_occurrence4()).
#  2. Normalize each event's raster to sum to 1 (a proper probability
#     surface), regardless of relocation count or migration duration.
#  3. Average an individual's per-year rasters together, then re-normalize
#     — an individual tracked for many years doesn't dominate the signal.
#  4. Extract that individual's migration-route footprint as the 99%
#     isopleth of their averaged raster.
#  5. Stack the individual footprints (as binary rasters) across the
#     population and count, at each pixel, how many individuals' footprints
#     include it.
#  6. Classify the stacked count into use tiers (any use / >=25% / >=50% of
#     sampled individuals) to produce the final corridor polygon(s).
#
# Gap exclusion is a practical parameter, not UI-exposed (see
# MIGRATION_DT_MAX4 below); the footprint isopleth level, grid resolution,
# and tier thresholds ARE UI-exposed, in the sidebar's Kriging/Mixed panel
# (footprint_level4, grid_res4, tier_moderate4, tier_high4).
# =============================================================================

MIGRATION_DT_MAX4 <- 8 * 3600  # exclude >= 8h gaps from bridging (Sawyer et al. 2009)

# Colors for the population-level use-tier layers drawn on the map when the
# Kriging workflow stacks more than one individual's footprint (Step 6).
# Light -> dark = low -> high use. Labels (which include the user-chosen
# tier percentages) are built once per estimate — see tier_labels4() below
# — and carried in the result rather than read live, so a later slider
# tweak can't desync a table/map already drawn from an earlier estimate.
TIER_COLORS4 <- c(Low = "#FEE08B", Moderate = "#FC8D59", High = "#B30000")

tier_labels4 <- function(moderate_pct, high_pct) {
  c(
    Low      = "Corridor: any use",
    Moderate = paste0("Corridor: ≥", moderate_pct, "% of population"),
    High     = paste0("Corridor: ≥", high_pct, "% of population")
  )
}

# Step 1: a single occurrence-distribution raster for one individual x year
# x migration-event chunk, on the shared `trast` grid, excluding gaps of
# MIGRATION_DT_MAX4 or more from the bridging. Mirrors amt:::od.track_xyt()
# internals directly (rather than calling amt::od()) because amt::od()
# does not forward "..." to ctmm::occurrence(), so it cannot exclude gaps.
bbmm_occurrence4 <- function(trk, ctmm_model, trast) {
  model <- fit_ctmm(trk, ctmm_model)
  krige <- ctmm::occurrence(amt::as_telemetry(trk), CTMM = model,
                             res.space = 10, res.time = 10,
                             dt.max = MIGRATION_DT_MAX4)
  r <- 1 - terra::rast(ctmm::raster(krige, DF = "CDF"))
  r <- terra::project(r, y = trast, method = "bilinear", gdal = TRUE, res = terra::res(trast))
  r <- terra::resample(r, trast)
  terra::ifel(is.na(r), 0, r)
}

# Normalize a raster so its values sum to 1 (Steps 2 & 3's "rescale").
normalize_raster4 <- function(r) {
  total <- sum(terra::values(r), na.rm = TRUE)
  if (total <= 0) return(r)
  r / total
}

# Steps 1-3 for one individual: fit + normalize each year's occurrence
# raster, then average across years and re-normalize. Returns list(ud =
# <SpatRaster>, years = <int vector>, by_year_n = <named list>, n = <int>)
# or NULL if no year could be fit.
individual_migration_ud4 <- function(d_ind, ctmm_model, trast) {
  yrs       <- sort(unique(d_ind$year))
  by_year_n <- list()
  rasters   <- list()
  for (yr in yrs) {
    d_yr <- d_ind |> filter(year == yr) |> arrange(timestamp)
    res_yr <- tryCatch({
      trk <- d_yr |>
        make_track(.x = longitude, .y = latitude, .t = timestamp, crs = 4326) |>
        transform_coords(crs_to = 32610)   # UTM zone 10N – adjust for your study area
      normalize_raster4(bbmm_occurrence4(trk, ctmm_model, trast))
    }, error = function(e) NULL)
    if (!is.null(res_yr)) {
      rasters[[as.character(yr)]]   <- res_yr
      by_year_n[[as.character(yr)]] <- nrow(d_yr)
    }
  }
  if (length(rasters) == 0) return(NULL)

  ud <- normalize_raster4(terra::mean(terra::rast(rasters)))
  list(ud = ud, years = as.integer(names(rasters)), by_year_n = by_year_n,
       n = sum(unlist(by_year_n)))
}

# Step 4: this individual's migration-route footprint — the `footprint_level`
# isopleth (sidebar: "Corridor footprint level") of their averaged
# occurrence raster.
migration_footprint4 <- function(ud, footprint_level) {
  hr_isopleths(ud, levels = footprint_level) |>
    filter(what == "estimate") |>
    st_transform(4326)
}

# Steps 5-6: stack a population's binary footprint rasters and classify
# into use tiers. `tier_props` is a named Low/Moderate/High vector of
# population proportions (see estimate_corridor_workflow4()). Returns a
# named list of sf polygons (possibly zero-row if a tier's threshold isn't
# met by any pixel), one per tier_props entry.
stack_footprints4 <- function(footprints_sf, trast, tier_props) {
  n <- length(footprints_sf)
  binary_rasters <- lapply(footprints_sf, function(fp) {
    terra::rasterize(terra::vect(st_transform(fp, terra::crs(trast))),
                      trast, field = 1, background = 0)
  })
  stacked <- Reduce(`+`, binary_rasters)

  counts <- ceiling(tier_props * n)
  counts[1] <- max(counts[1], 1)  # "Low" = at least one individual

  setNames(lapply(counts, function(k) {
    tier_r <- terra::ifel(stacked >= k, 1, NA)
    if (all(is.na(terra::values(tier_r)))) return(st_sf(geometry = st_sfc(crs = 4326)))
    terra::as.polygons(tier_r, dissolve = TRUE) |>
      st_as_sf() |>
      st_set_crs(terra::crs(trast)) |>
      st_transform(4326)
  }), names(tier_props))
}

# Orchestrates the full workflow for the current selection. Returns list(
#   per_id      = <named list of each individual's footprint sf>,
#   per_id_years = <named list of each individual's year(s)>,
#   per_id_by_year_n = <named list (by individual) of named lists (by year)
#                       of relocation counts>,
#   per_id_n    = <named list of each individual's total relocation count>,
#   combined    = <sf, union of all individual footprints — the "any use"
#                  / Low-tier corridor, used when this is the only method
#                  or as the Kriging component of "Mixed approach">,
#   combined_n  = <total relocation count across individuals>,
#   tiers       = <NULL for a single individual, else named list of
#                  Low/Moderate/High tier sf polygons from stacking>,
#   tier_labels = <NULL for a single individual, else the Low/Moderate/High
#                  labels used for "tiers", captured at estimate time so a
#                  later change to the tier sliders can't desync an
#                  already-drawn map/table from this result>
# )
# footprint_level/grid_res/tier_moderate/tier_high are read once here from
# the sidebar's Kriging/Mixed inputs, rather than by the helper functions
# themselves, so every step of one estimate run uses a single consistent
# set of parameters even if the user drags a slider mid-run.
estimate_corridor_workflow4 <- function(d) {
  footprint_level <- input$footprint_level4
  grid_res        <- input$grid_res4
  tier_props      <- c(Low = 0, Moderate = input$tier_moderate4 / 100, High = input$tier_high4 / 100)

  ids <- sort(unique(d$id))
  trk_all <- d |>
    make_track(.x = longitude, .y = latitude, .t = timestamp, crs = 4326) |>
    transform_coords(crs_to = 32610)
  trast <- make_trast(trk_all, res = grid_res)

  per_id           <- list()
  per_id_years     <- list()
  per_id_by_year_n <- list()
  per_id_n         <- list()
  for (i in ids) {
    ud_i <- tryCatch(individual_migration_ud4(d |> filter(id == i), input$ctmm_model4, trast),
                      error = function(e) NULL)
    if (is.null(ud_i)) next
    fp_i <- tryCatch(migration_footprint4(ud_i$ud, footprint_level), error = function(e) NULL)
    if (is.null(fp_i)) next
    per_id[[i]]           <- fp_i
    per_id_years[[i]]     <- ud_i$years
    per_id_by_year_n[[i]] <- ud_i$by_year_n
    per_id_n[[i]]         <- ud_i$n
  }
  req(length(per_id) > 0)

  tiers  <- if (length(per_id) > 1) stack_footprints4(per_id, trast, tier_props) else NULL
  labels <- if (!is.null(tiers)) tier_labels4(input$tier_moderate4, input$tier_high4) else NULL
  combined <- if (!is.null(tiers)) tiers[["Low"]] else per_id[[1]]

  list(
    per_id           = per_id,
    per_id_years     = per_id_years,
    per_id_by_year_n = per_id_by_year_n,
    per_id_n         = per_id_n,
    combined         = combined,
    combined_n       = sum(unlist(per_id_n)),
    tier_labels      = labels,
    tiers            = tiers
  )
}

# ---------------------------------------------------------------------------
# Corridor estimate (recomputed only when the button is pressed, since the
# Kriging method in particular can take a while — "auto" model selection
# and "ALL" individuals especially so)
# ---------------------------------------------------------------------------

# Line buffer method for the current selection: 1) estimate each individual
# (itself a per-year union, see estimate_corridor_individual4()), 2) union
# the per-individual results into a single combined polygon for "ALL".
# Returns list(
#   combined         = <sf, used for the map/table>,
#   per_id           = <named list of per-individual combined sf estimates,
#                       used by the Details tab plots>,
#   per_id_by_year   = <named list (by individual) of named lists (by year)
#                       of sf estimates — each individual's own yearly
#                       contours, NULL when only one year is selected>,
#   per_id_years     = <named list (by individual) of the year(s) actually
#                       used for that individual — used by the estimated-
#                       areas table>,
#   per_id_by_year_n = <named list (by individual) of named lists (by year)
#                       of relocation counts, used by the estimated-areas
#                       table's N column. NULL when only one year>,
#   per_id_n         = <named list (by individual) of total relocation
#                       count used for that individual's own estimate>,
#   combined_n       = <total relocation count behind the "combined"
#                       estimate — the sum of per_id_n>,
#   by_year          = <named list of each year's own estimate, combined
#                       across individuals when "ALL" is selected — the
#                       yearly contours drawn on the map. NULL when only one
#                       year is selected>
# )
compute_line_buffer_result4 <- function(d) {
  if (input$id4 == "ALL") {
    per_id           <- list()
    per_id_by_year   <- list()
    per_id_years     <- list()
    per_id_by_year_n <- list()
    per_id_n         <- list()
    for (i in sort(unique(d$id))) {
      res_i <- tryCatch(estimate_corridor_individual4(d |> filter(id == i)),
                         error = function(e) NULL)
      if (!is.null(res_i)) {
        per_id[[i]]       <- res_i$combined
        per_id_years[[i]] <- res_i$years
        per_id_n[[i]]     <- res_i$n
        if (!is.null(res_i$by_year)) {
          per_id_by_year[[i]]   <- res_i$by_year
          per_id_by_year_n[[i]] <- res_i$by_year_n
        }
      }
    }
    req(length(per_id) > 0)

    # One contour per year: union that year's polygon across whichever
    # individuals have data for it.
    by_year <- NULL
    if (length(per_id_by_year) > 0) {
      yrs <- sort(unique(unlist(lapply(per_id_by_year, names))))
      by_year <- list()
      for (yr in yrs) {
        yr_ests <- lapply(per_id_by_year, `[[`, yr)
        yr_ests <- yr_ests[!vapply(yr_ests, is.null, logical(1))]
        if (length(yr_ests) > 0) by_year[[yr]] <- union_by_level(yr_ests)
      }
    }

    list(
      combined         = union_by_level(per_id),
      per_id           = per_id,
      per_id_by_year   = per_id_by_year,
      per_id_years     = per_id_years,
      per_id_by_year_n = per_id_by_year_n,
      per_id_n         = per_id_n,
      combined_n       = sum(unlist(per_id_n)),
      by_year          = by_year
    )
  } else {
    res <- estimate_corridor_individual4(d)
    list(
      combined         = res$combined,
      per_id           = setNames(list(res$combined), input$id4),
      per_id_by_year   = setNames(list(res$by_year), input$id4),
      per_id_years     = setNames(list(res$years), input$id4),
      per_id_by_year_n = setNames(list(res$by_year_n), input$id4),
      per_id_n         = setNames(list(res$n), input$id4),
      combined_n       = res$n,
      by_year          = res$by_year
    )
  }
}

# Dissolve two sf objects (each possibly multi-row) into a single unioned
# geometry, tolerating either side being NULL/empty (an individual that
# only one of the two pathways could fit).
union_two_sf4 <- function(a, b) {
  a_empty <- is.null(a) || nrow(a) == 0
  b_empty <- is.null(b) || nrow(b) == 0
  if (a_empty && b_empty) return(NULL)
  if (a_empty) return(b)
  if (b_empty) return(a)
  st_sf(geometry = st_union(st_union(st_geometry(a)), st_union(st_geometry(b)))) |>
    st_transform(4326)
}

# Combine the Line buffer and/or Kriging-workflow results into the unified
# shape consumed by the map observer, area table, and population plot below.
# "tiers" is only ever non-NULL for the Kriging workflow's "ALL" case (see
# estimate_corridor_workflow4()); Line buffer never populates it, and Mixed
# approach doesn't propagate it since the mixed "combined"/per-individual
# polygons no longer correspond 1:1 with the stacked-footprint tiers.
combine_corridor_results4 <- function(method, buffer_res, kriging_res) {
  if (method == "Line buffer") {
    return(c(buffer_res, list(tiers = NULL, tier_labels = NULL)))
  }
  if (method == "Kriging (CTMM)") {
    return(list(
      combined         = kriging_res$combined,
      per_id           = kriging_res$per_id,
      per_id_by_year   = NULL,
      per_id_years     = kriging_res$per_id_years,
      per_id_by_year_n = kriging_res$per_id_by_year_n,
      per_id_n         = kriging_res$per_id_n,
      combined_n       = kriging_res$combined_n,
      by_year          = NULL,
      tiers            = kriging_res$tiers,
      tier_labels      = kriging_res$tier_labels
    ))
  }

  # Mixed approach: union each pathway's output, individual by individual.
  ids <- union(names(buffer_res$per_id), names(kriging_res$per_id))
  per_id <- setNames(
    lapply(ids, function(i) union_two_sf4(buffer_res$per_id[[i]], kriging_res$per_id[[i]])),
    ids
  )
  per_id_n <- setNames(
    lapply(ids, function(i) max(buffer_res$per_id_n[[i]] %||% 0, kriging_res$per_id_n[[i]] %||% 0)),
    ids
  )
  per_id_years <- setNames(
    lapply(ids, function(i) sort(unique(c(buffer_res$per_id_years[[i]], kriging_res$per_id_years[[i]])))),
    ids
  )

  list(
    combined         = union_two_sf4(buffer_res$combined, kriging_res$combined),
    per_id           = per_id,
    per_id_by_year   = buffer_res$per_id_by_year,
    per_id_years     = per_id_years,
    per_id_by_year_n = buffer_res$per_id_by_year_n,
    per_id_n         = per_id_n,
    combined_n       = max(buffer_res$combined_n, kriging_res$combined_n),
    by_year          = buffer_res$by_year,
    tiers            = NULL,
    tier_labels      = NULL
  )
}

corridor4 <- eventReactive(input$estimate_path4, {
  d <- path_data4()

  if (input$id4 == "ALL") {
    showNotification(
      "Estimating corridors for all individuals — this could take a while…",
      id = "corridor4_all_busy", duration = NULL, type = "warning"
    )
    on.exit(removeNotification("corridor4_all_busy"), add = TRUE)
  } else if (input$method4 %in% c("Kriging (CTMM)", "Mixed approach")) {
    showNotification("Estimating corridor — Kriging can take a minute…",
                      id = "corridor4_busy", duration = NULL, type = "message")
    on.exit(removeNotification("corridor4_busy"), add = TRUE)
  }

  result <- tryCatch({
    buffer_res  <- if (input$method4 %in% c("Line buffer", "Mixed approach"))
      compute_line_buffer_result4(d) else NULL
    kriging_res <- if (input$method4 %in% c("Kriging (CTMM)", "Mixed approach"))
      estimate_corridor_workflow4(d) else NULL
    combine_corridor_results4(input$method4, buffer_res, kriging_res)
  }, error = function(e) {
    showNotification(
      paste("Could not identify the corridor:", conditionMessage(e)),
      type = "error"
    )
    NULL
  })
  req(result)
  result
})

# ---------------------------------------------------------------------------
# Base map (rendered once)
# ---------------------------------------------------------------------------

output$map4 <- renderLeaflet({
  leaflet(options = leafletOptions(attributionControl = FALSE)) |>
    addProviderTiles("Esri.WorldImagery",    group = "Esri.WorldImagery") |>
    addProviderTiles("Esri.WorldGrayCanvas", group = "Esri.WorldGrayCanvas") |>
    addProviderTiles("Esri.WorldTopoMap",    group = "Esri.WorldTopoMap") |>
    addLayersControl(
      position      = "topright",
      baseGroups    = c("Esri.WorldTopoMap", "Esri.WorldImagery", "Esri.WorldGrayCanvas"),
      overlayGroups = c("Boundary", "Tracks", "GPS", "Corridor"),
      options       = layersControlOptions(collapsed = FALSE)
    )
})
# See the note on output$map above: this page starts hidden, so force the
# render regardless of visibility or leafletProxy() calls made before the
# user first visits this page (e.g. drawing the Boundary as soon as data
# loads) get silently dropped.
outputOptions(output, "map4", suspendWhenHidden = FALSE)

# Population-level KDE boundary, shared with the other pages. Also gives
# the map its initial zoom/extent — without this, the map stays at its
# default world view (so Boundary/GPS/Tracks are all technically drawn but
# imperceptibly small) until "Identify corridor" is pressed.
observeEvent(herd_data(), {
  bnd <- herd_data()$bnd
  leafletProxy("map4") |>
    clearGroup("Boundary") |>
    addPolygons(
      data = bnd,
      fillOpacity = 0, stroke = TRUE,
      color = "black", weight = 2,
      group = "Boundary"
    )

  bb <- st_bbox(bnd)
  leafletProxy("map4") |> fitBounds(bb[["xmin"]], bb[["ymin"]], bb[["xmax"]], bb[["ymax"]])
})

# GPS points / tracks for the current selection. Skipped for "ALL" so only
# the final combined corridor is shown (drawing every individual's points
# and tracks at once would be slow and cluttered).
tracks_hidden_once4 <- reactiveVal(FALSE)

observe({
  d <- path_data4()

  proxy <- leafletProxy("map4") |> clearGroup("GPS") |> clearGroup("Tracks")
  req(nrow(d) > 0, input$id4 != "ALL")

  gps_sf4 <- st_as_sf(d, coords = c("longitude", "latitude"), crs = 4326)
  proxy |>
    addCircles(
      data        = gps_sf4,
      color       = "red", fillColor = "red",
      fillOpacity = 0.5, weight = 1, radius = 1,
      group       = "GPS"
    )

  trk_sf4 <- d |>
    arrange(timestamp) |>
    group_by(id, year, phase) |>
    filter(n() >= 2) |>
    summarise(
      geometry = st_sfc(st_linestring(cbind(longitude, latitude))),
      .groups  = "drop"
    ) |>
    st_as_sf(crs = 4326)

  if (nrow(trk_sf4) > 0) {
    trk_proxy <- leafletProxy("map4") |>
      addPolylines(data = trk_sf4, color = "black", weight = 1, opacity = 0.7, group = "Tracks")
    # Only hide Tracks the first time it's drawn — every redraw after that
    # (e.g. from changing a filter) should respect whatever the user has
    # since toggled it to, not silently re-hide it.
    if (!isTRUE(tracks_hidden_once4())) {
      trk_proxy |> hideGroup("Tracks")
      tracks_hidden_once4(TRUE)
    }
  }
})

# Identified corridor, plus one separately selectable outline-only contour
# layer per year when more than one year was estimated (e.g. "2022",
# "2023", ... in the layer control). For the Kriging workflow's "ALL"
# individuals, the single "Corridor" fill is replaced by one separately
# selectable filled layer per population use tier (see tier_labels4()) —
# light-to-dark from "any individual" to the user-chosen "high use" %.
prev_year_groups4   <- reactiveVal(character(0))
year_groups_hidden4 <- reactiveVal(character(0))
# Tier group names include the user-chosen tier percentages (via
# res$tier_labels), so they can change between estimates — track what was
# actually drawn last time (like prev_year_groups4) rather than a fixed
# set, or a changed slider would leave stale layers uncleared.
prev_tier_groups4   <- reactiveVal(character(0))

observeEvent(corridor4(), {
  res <- corridor4()
  cor <- res$combined
  req(cor, nrow(cor) > 0)

  proxy <- leafletProxy("map4") |> clearGroup("Corridor")
  for (g in prev_tier_groups4()) proxy <- proxy |> clearGroup(g)

  # Remove any year-specific layers left over from a previous estimate.
  for (g in prev_year_groups4()) proxy <- proxy |> clearGroup(g)

  tier_groups <- character(0)
  if (!is.null(res$tiers)) {
    for (tier in names(res$tier_labels)) {
      tpoly <- res$tiers[[tier]]
      if (is.null(tpoly) || nrow(tpoly) == 0) next
      lbl <- res$tier_labels[[tier]]
      tier_groups <- c(tier_groups, lbl)
      proxy <- proxy |>
        addPolygons(
          data        = tpoly,
          color       = TIER_COLORS4[[tier]], weight = 2, fill = TRUE, fillOpacity = 0.5,
          group       = lbl,
          label       = lbl
        )
    }
  } else {
    proxy <- proxy |>
      addPolygons(
        data        = cor,
        color       = "#CC3333", weight = 2, fill = TRUE, fillOpacity = 0.4,
        group       = "Corridor"
      )
  }
  prev_tier_groups4(tier_groups)

  year_groups <- character(0)
  if (!is.null(res$by_year) && length(res$by_year) > 0) {
    yr_pal         <- colorFactor("Dark2", domain = names(res$by_year))
    already_hidden <- year_groups_hidden4()
    newly_hidden   <- character(0)
    for (yr in names(res$by_year)) {
      year_groups <- c(year_groups, yr)
      proxy <- proxy |>
        addPolygons(
          data      = res$by_year[[yr]],
          color     = yr_pal(yr), weight = 3, fill = FALSE, dashArray = "4,4",
          group     = yr,
          label     = yr
        )
      # Only hide a given year layer the first time it ever appears — once
      # the user has toggled it, later re-estimates shouldn't re-hide it.
      if (!(yr %in% already_hidden)) {
        proxy <- proxy |> hideGroup(yr)
        newly_hidden <- c(newly_hidden, yr)
      }
    }
    if (length(newly_hidden) > 0) year_groups_hidden4(c(already_hidden, newly_hidden))
  }
  prev_year_groups4(year_groups)

  proxy |>
    addLayersControl(
      position      = "topright",
      baseGroups    = c("Esri.WorldTopoMap", "Esri.WorldImagery", "Esri.WorldGrayCanvas"),
      overlayGroups = c("Boundary", "Tracks", "GPS",
                         if (is.null(res$tiers)) "Corridor" else tier_groups,
                         year_groups),
      options       = layersControlOptions(collapsed = FALSE)
    )

  bb <- st_bbox(cor)
  leafletProxy("map4") |> fitBounds(bb[["xmin"]], bb[["ymin"]], bb[["xmax"]], bb[["ymax"]])
})

# ---------------------------------------------------------------------------
# Summary stats / area-by-level side panel
# ---------------------------------------------------------------------------

output$summary_stats4 <- renderTable({
  d <- path_data4()
  req(nrow(d) > 0)
  tibble(
    Field = c("Caribou", "Season", "Start date", "End date", "Relocations"),
    Value = c(
      input$id4,
      input$phase4,
      as.character(as.Date(min(d$timestamp))),
      as.character(as.Date(max(d$timestamp))),
      as.character(nrow(d))
    )
  )
}, colnames = FALSE, striped = TRUE, spacing = "xs")

# One row per year for each individual, plus a "Total" row for any
# individual (or the population, when "ALL") that combines more than one
# year — e.g. 4 years for one caribou gives 5 rows: one per year and a
# total. A single-year individual gets just its one row. The Kriging
# workflow doesn't keep separate per-year polygons (years are averaged
# together per corridor_workflow.md Step 3), so a multi-year individual
# there gets one row spanning the year range instead. When the Kriging
# workflow stacked more than one individual's footprint, one extra row per
# population use tier is appended (see tier_labels4()).
output$area_by_level4 <- renderTable({
  res <- corridor4()
  req(res, length(res$per_id) > 0)

  ids <- sort(names(res$per_id))

  rows <- list()
  for (id in ids) {
    yrs_i       <- res$per_id_years[[id]]
    by_year_i   <- res$per_id_by_year[[id]]
    by_year_n_i <- res$per_id_by_year_n[[id]]

    if (!is.null(by_year_i) && length(yrs_i) > 1) {
      for (yr in as.character(sort(yrs_i))) {
        rows[[length(rows) + 1]] <- tibble(
          Caribou = id, Year = yr, N = by_year_n_i[[yr]],
          Area_km2 = area_km2(by_year_i[[yr]])
        )
      }
      rows[[length(rows) + 1]] <- tibble(
        Caribou = id, Year = "Total", N = res$per_id_n[[id]],
        Area_km2 = area_km2(res$per_id[[id]])
      )
    } else if (length(yrs_i) > 1) {
      rows[[length(rows) + 1]] <- tibble(
        Caribou = id, Year = paste0(min(yrs_i), "–", max(yrs_i), " (avg)"),
        N = res$per_id_n[[id]], Area_km2 = area_km2(res$per_id[[id]])
      )
    } else {
      rows[[length(rows) + 1]] <- tibble(
        Caribou = id, Year = as.character(yrs_i[1]), N = res$per_id_n[[id]],
        Area_km2 = area_km2(res$per_id[[id]])
      )
    }
  }

  if (length(ids) > 1) {
    rows[[length(rows) + 1]] <- tibble(
      Caribou = input$id4, Year = "Total", N = res$combined_n,
      Area_km2 = area_km2(res$combined)
    )
  }

  if (!is.null(res$tiers)) {
    for (tier in names(res$tier_labels)) {
      tpoly <- res$tiers[[tier]]
      if (is.null(tpoly) || nrow(tpoly) == 0) next
      rows[[length(rows) + 1]] <- tibble(
        Caribou = "Population", Year = res$tier_labels[[tier]],
        N = res$combined_n, Area_km2 = area_km2(tpoly)
      )
    }
  }

  bind_rows(rows)
}, striped = TRUE, spacing = "xs")

# One panel per individual, shown only when comparing "ALL" individuals
# (a single-individual estimate has nothing to compare, so the tab stays
# empty). facet_wrap()'s default fixed scales give every panel the same
# spatial extent, so corridor size/shape is directly comparable across
# individuals. When more than one year was estimated, each individual's
# own yearly contours (outline only, no fill) are overlaid on their panel.
output$corridor_by_id_plot4 <- renderPlot({
  req(input$id4 == "ALL")
  res    <- corridor4()
  per_id <- res$per_id
  req(length(per_id) > 0)

  plot_df <- bind_rows(lapply(names(per_id), function(nm) {
    p <- per_id[[nm]]
    p$id    <- nm
    p$Level <- if ("level" %in% names(p)) paste0(p$level * 100, "%") else "Corridor"
    p
  }))

  p <- ggplot(plot_df) +
    geom_sf(aes(fill = Level), color = "black", linewidth = 0.3, alpha = 0.6) +
    facet_wrap(~id) +
    scale_fill_brewer(palette = "YlOrRd", name = NULL) +
    theme_minimal() +
    theme(
      axis.text    = element_text(size = 7),
      strip.text   = element_text(size = 14, face = "bold"),
      panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.5)
    ) +
    labs(x = NULL, y = NULL)

  year_df <- per_id_year_to_df(res$per_id_by_year)
  if (!is.null(year_df)) {
    p <- p +
      geom_sf(data = year_df, aes(color = Year), fill = NA,
              linewidth = 1, linetype = "dashed") +
      scale_color_brewer(palette = "Dark2", name = "Year")
  }

  p
})
