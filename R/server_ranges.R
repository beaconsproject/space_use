# Server logic for the "Seasonal ranges" dashboard page.


# ---------------------------------------------------------------------------
# Seasonal ranges: sidebar choices follow the currently loaded dataset
# (same raw_data() reactive used by the Movement Data page)
# ---------------------------------------------------------------------------

observeEvent(raw_data(), {
  x <- raw_data()
  ids3 <- as.character(sort(unique(x$id)))
  updateSelectInput(session, "id3",
    choices  = c("ALL", ids3),
    selected = ids3[1]
  )
  updateSelectInput(session, "phase3",
    choices  = c("Annual", sort(unique(x$season))),
    selected = "Annual"
  )
  updateSliderInput(session, "year_range3",
    min   = min(x$year),
    max   = max(x$year),
    value = c(min(x$year), max(x$year))
  )
})

# ---------------------------------------------------------------------------
# Filtered subset used for both the map layers and the range estimate
# ---------------------------------------------------------------------------

range_data3 <- reactive({
  req(raw_data(), input$id3, input$phase3, input$year_range3)
  x <- raw_data() |>
    filter(year >= input$year_range3[1], year <= input$year_range3[2])
  if (input$id3 != "ALL")    x <- x |> filter(id == input$id3)
  if (input$phase3 != "Annual") x <- x |> filter(season == input$phase3)
  x
})

# ---------------------------------------------------------------------------
# Home range estimate (recomputed only when the button is pressed, since
# aKDE model-fitting can be slow)
# ---------------------------------------------------------------------------

# Estimate isopleths for one individual x year chunk (level/what/area/
# geometry, in 4326) using the selected estimator (MCP/KDE/aKDE).
estimate_hr_fit3 <- function(d_one, levels_sel) {
  trk <- d_one |>
    make_track(.x = longitude, .y = latitude, .t = timestamp, crs = 4326) |>
    transform_coords(crs_to = 32610)   # UTM zone 10N – adjust for your study area

  if (input$hr_method3 == "MCP") {
    fit <- hr_mcp(trk, levels = levels_sel)
  } else if (input$hr_method3 == "KDE") {
    h_ref <- hr_kde_ref(trk)
    h_use <- if (input$h3 > 0) h_ref * input$h3 else h_ref
    fit <- hr_kde(trk, h = h_use, levels = levels_sel)
  } else {
    # aKDE via Kriging of the fitted continuous-time movement model
    # (same amt/ctmm machinery used by the Seasonal paths page).
    model <- fit_ctmm(trk, input$ctmm_model3)
    fit <- hr_akde(trk, model = model, levels = levels_sel)
  }

  # aKDE also returns lower/upper confidence bounds per level; keep only the
  # point estimate so the map shows one polygon per selected isopleth level.
  hr_isopleths(fit) |>
    filter(what == "estimate") |>
    st_transform(4326)
}

# BBMM (Brownian bridge) estimate for migration-period relocations only,
# used to fill in the "Include migration periods?" portion of an Annual
# range — the same amt::od()/ctmm::occurrence() Kriging machinery as the
# Seasonal paths page's Kriging method, with a fixed Brownian-motion model.
estimate_bbmm3 <- function(d_mig, levels_sel) {
  trk <- d_mig |>
    make_track(.x = longitude, .y = latitude, .t = timestamp, crs = 4326) |>
    transform_coords(crs_to = 32610)   # UTM zone 10N – adjust for your study area

  model <- fit_ctmm(trk, "bm")
  od(trk, model = model, trast = make_trast(trk)) |>
    hr_isopleths(levels = levels_sel) |>
    st_transform(4326)
}

# Estimate isopleths for one individual x year chunk. For any season other
# than "Annual", this is just the selected estimator. For "Annual", the
# "Include migration periods?" checkbox controls whether migration-period
# relocations (input$migration is non-NA) are excluded from the estimate,
# or estimated separately via BBMM and unioned into the result — movement
# during migration doesn't reflect stationary-range space use the way
# MCP/KDE/aKDE assume, so it's handled with a movement (bridge) model
# instead of being pooled in with the rest of the year's relocations.
#
# Returns list(combined = <sf, the range for this chunk>, migration = <sf,
# the migration-only BBMM estimate that was unioned into `combined`, or
# NULL if migration periods weren't estimated separately> — `migration` is
# what lets the map highlight which part of an annual range is the
# migration corridor.
estimate_hr_one3 <- function(d_one, levels_sel) {
  if (input$phase3 != "Annual") {
    return(list(combined = estimate_hr_fit3(d_one, levels_sel), migration = NULL))
  }

  is_mig   <- !is.na(d_one$migration)
  d_nonmig <- d_one[!is_mig, , drop = FALSE]
  d_mig    <- d_one[is_mig, , drop = FALSE]

  if (!isTRUE(input$include_migration3) || nrow(d_mig) == 0) {
    return(list(combined = estimate_hr_fit3(d_nonmig, levels_sel), migration = NULL))
  }

  fit_mig <- tryCatch(estimate_bbmm3(d_mig, levels_sel), error = function(e) NULL)

  if (nrow(d_nonmig) == 0) {
    req(fit_mig)
    return(list(combined = fit_mig, migration = fit_mig))
  }

  fit_nonmig <- estimate_hr_fit3(d_nonmig, levels_sel)
  combined   <- if (is.null(fit_mig)) fit_nonmig else union_by_level(list(fit_nonmig, fit_mig))
  list(combined = combined, migration = fit_mig)
}

# Number of relocations actually used for one individual x year chunk,
# consistent with estimate_hr_one3()'s migration-period handling (used for
# the estimated-areas table's N column).
relocation_count3 <- function(d_one) {
  if (input$phase3 == "Annual" && !isTRUE(input$include_migration3)) {
    sum(is.na(d_one$migration))
  } else {
    nrow(d_one)
  }
}

# One individual across the selected year range. Returns list(combined =
# <sf, unioned across years>, by_year = <named list of each year's own
# estimate, or NULL when only one year is selected>, years = <integer
# vector of the year(s) actually used for this individual>, by_year_n =
# <named list of each year's relocation count, or NULL when only one year
# is selected>, n = <total relocation count for this individual>, migration
# = <sf, this individual's migration-only estimate unioned across years, or
# NULL when migration periods weren't estimated separately for any year>).
estimate_hr_individual3 <- function(d_ind, levels_sel) {
  yrs <- sort(unique(d_ind$year))
  if (length(yrs) <= 1) {
    res <- estimate_hr_one3(d_ind, levels_sel)
    return(list(combined = res$combined, by_year = NULL,
                years = yrs, by_year_n = NULL, n = relocation_count3(d_ind),
                migration = res$migration))
  }

  by_year        <- list()
  by_year_n      <- list()
  migration_list <- list()
  for (yr in yrs) {
    d_yr   <- d_ind |> filter(year == yr)
    res_yr <- tryCatch(estimate_hr_one3(d_yr, levels_sel), error = function(e) NULL)
    if (!is.null(res_yr)) {
      by_year[[as.character(yr)]]   <- res_yr$combined
      by_year_n[[as.character(yr)]] <- relocation_count3(d_yr)
      if (!is.null(res_yr$migration)) migration_list[[as.character(yr)]] <- res_yr$migration
    }
  }
  req(length(by_year) > 0)
  migration_combined <- if (length(migration_list) > 0) union_by_level(migration_list) else NULL
  list(combined = union_by_level(by_year), by_year = by_year,
       years = as.integer(names(by_year)), by_year_n = by_year_n,
       n = sum(unlist(by_year_n)), migration = migration_combined)
}

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
#                       year is selected>,
#   migration        = <sf, the migration-only estimate — unioned across
#                       years and, when "ALL", across individuals too —
#                       drawn on the map as the "Migration corridor" layer
#                       so it's visible within the overall annual range.
#                       NULL unless Season is "Annual" and "Include
#                       migration periods?" contributed to the estimate>
# )
# Each individual's own estimate already combines multiple years (see
# estimate_hr_individual3()); for "ALL", those per-individual estimates are
# then unioned again across individuals.
hr3 <- eventReactive(input$estimate_range3, {
  d          <- range_data3()
  levels_sel <- input$levels3
  req(length(levels_sel) > 0)

  uses_bbmm <- input$phase3 == "Annual" && isTRUE(input$include_migration3)

  if (input$id3 == "ALL") {
    showNotification(
      "Estimating ranges for all individuals — this could take a while…",
      id = "hr3_all_busy", duration = NULL, type = "warning"
    )
    on.exit(removeNotification("hr3_all_busy"), add = TRUE)
  } else if (input$hr_method3 == "aKDE" || uses_bbmm) {
    showNotification("Estimating range — aKDE/BBMM Kriging can take a minute…",
                      id = "hr3_busy", duration = NULL, type = "message")
    on.exit(removeNotification("hr3_busy"), add = TRUE)
  }

  result <- tryCatch({
    if (input$id3 == "ALL") {
      # 1) estimate each individual (itself a per-year union, see above),
      # 2) union each level across individuals into a combined polygon.
      per_id           <- list()
      per_id_by_year   <- list()
      per_id_years     <- list()
      per_id_by_year_n <- list()
      per_id_n         <- list()
      per_id_migration <- list()
      for (i in sort(unique(d$id))) {
        res_i <- tryCatch(estimate_hr_individual3(d |> filter(id == i), levels_sel),
                           error = function(e) NULL)
        if (!is.null(res_i)) {
          per_id[[i]]       <- res_i$combined
          per_id_years[[i]] <- res_i$years
          per_id_n[[i]]     <- res_i$n
          if (!is.null(res_i$by_year)) {
            per_id_by_year[[i]]   <- res_i$by_year
            per_id_by_year_n[[i]] <- res_i$by_year_n
          }
          if (!is.null(res_i$migration)) per_id_migration[[i]] <- res_i$migration
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

      migration <- if (length(per_id_migration) > 0) union_by_level(per_id_migration) else NULL

      list(
        combined         = union_by_level(per_id),
        per_id           = per_id,
        per_id_by_year   = per_id_by_year,
        per_id_years     = per_id_years,
        per_id_by_year_n = per_id_by_year_n,
        per_id_n         = per_id_n,
        combined_n       = sum(unlist(per_id_n)),
        by_year          = by_year,
        migration        = migration
      )
    } else {
      res <- estimate_hr_individual3(d, levels_sel)
      list(
        combined         = res$combined,
        per_id           = setNames(list(res$combined), input$id3),
        per_id_by_year   = setNames(list(res$by_year), input$id3),
        per_id_years     = setNames(list(res$years), input$id3),
        per_id_by_year_n = setNames(list(res$by_year_n), input$id3),
        per_id_n         = setNames(list(res$n), input$id3),
        combined_n       = res$n,
        by_year          = res$by_year,
        migration        = res$migration
      )
    }
  }, error = function(e) {
    showNotification(
      paste("Could not estimate the seasonal range:", conditionMessage(e)),
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

output$map3 <- renderLeaflet({
  leaflet(options = leafletOptions(attributionControl = FALSE)) |>
    addProviderTiles("Esri.WorldImagery",    group = "Esri.WorldImagery") |>
    addProviderTiles("Esri.WorldGrayCanvas", group = "Esri.WorldGrayCanvas") |>
    addProviderTiles("Esri.WorldTopoMap",    group = "Esri.WorldTopoMap") |>
    addLayersControl(
      position      = "topright",
      baseGroups    = c("Esri.WorldTopoMap", "Esri.WorldImagery", "Esri.WorldGrayCanvas"),
      overlayGroups = c("Boundary", "Tracks", "GPS", "Seasonal range"),
      options       = layersControlOptions(collapsed = FALSE)
    )
})
# See the note on output$map above: this page starts hidden, so force the
# render regardless of visibility or leafletProxy() calls made before the
# user first visits this page (e.g. drawing the Boundary as soon as data
# loads) get silently dropped.
outputOptions(output, "map3", suspendWhenHidden = FALSE)

# Population-level KDE boundary, shared with the Movement Data page. Also
# gives the map its initial zoom/extent — without this, the map stays at
# its default world view (so Boundary/GPS/Tracks are all technically drawn
# but imperceptibly small) until "Estimate range" is pressed.
observeEvent(herd_data(), {
  bnd <- herd_data()$bnd
  leafletProxy("map3") |>
    clearGroup("Boundary") |>
    addPolygons(
      data = bnd,
      fillOpacity = 0, stroke = TRUE,
      color = "black", weight = 2,
      group = "Boundary"
    )

  bb <- st_bbox(bnd)
  leafletProxy("map3") |> fitBounds(bb[["xmin"]], bb[["ymin"]], bb[["xmax"]], bb[["ymax"]])
})

# GPS points / tracks for the current selection. Skipped for "ALL" so only
# the final combined range is shown (drawing every individual's points and
# tracks at once would be slow and cluttered).
tracks_hidden_once3 <- reactiveVal(FALSE)

observe({
  d <- range_data3()

  proxy <- leafletProxy("map3") |> clearGroup("GPS") |> clearGroup("Tracks")
  req(nrow(d) > 0, input$id3 != "ALL")

  gps_sf3 <- st_as_sf(d, coords = c("longitude", "latitude"), crs = 4326)
  proxy |>
    addCircles(
      data        = gps_sf3,
      color       = "red", fillColor = "red",
      fillOpacity = 0.5, weight = 1, radius = 1,
      group       = "GPS"
    )

  trk_sf3 <- d |>
    arrange(timestamp) |>
    group_by(id, year, phase) |>
    filter(n() >= 2) |>
    summarise(
      geometry = st_sfc(st_linestring(cbind(longitude, latitude))),
      .groups  = "drop"
    ) |>
    st_as_sf(crs = 4326)

  if (nrow(trk_sf3) > 0) {
    trk_proxy <- leafletProxy("map3") |>
      addPolylines(data = trk_sf3, color = "black", weight = 1, opacity = 0.7, group = "Tracks")
    # Only hide Tracks the first time it's drawn — every redraw after that
    # (e.g. from changing a filter) should respect whatever the user has
    # since toggled it to, not silently re-hide it.
    if (!isTRUE(tracks_hidden_once3())) {
      trk_proxy |> hideGroup("Tracks")
      tracks_hidden_once3(TRUE)
    }
  }
})

# Estimated seasonal range for the selected isopleth level (a single fixed
# color, since only one level can be selected at a time), plus:
#  - one separately selectable outline-only contour layer per year when
#    more than one year was estimated (e.g. "2022", "2023", ... in the
#    layer control), and
#  - a "Migration corridor" layer highlighting the migration-period
#    (BBMM) portion of an Annual range, when "Include migration periods?"
#    contributed to the estimate.
RANGE_COLOR3     <- "#2166AC"
MIGRATION_COLOR3 <- "#E67E22"
prev_year_groups3   <- reactiveVal(character(0))
year_groups_hidden3 <- reactiveVal(character(0))

observeEvent(hr3(), {
  res <- hr3()
  rng <- res$combined
  req(rng, nrow(rng) > 0)

  proxy <- leafletProxy("map3") |>
    clearGroup("Seasonal range") |>
    clearGroup("Migration corridor")

  # Remove any year-specific layers left over from a previous estimate.
  for (g in prev_year_groups3()) proxy <- proxy |> clearGroup(g)

  proxy <- proxy |>
    addPolygons(
      data        = rng,
      color       = RANGE_COLOR3, weight = 2, fillOpacity = 0.2,
      group       = "Seasonal range",
      label       = ~paste0(level * 100, "% isopleth")
    )

  has_migration <- !is.null(res$migration) && nrow(res$migration) > 0
  if (has_migration) {
    proxy <- proxy |>
      addPolygons(
        data        = res$migration,
        color       = MIGRATION_COLOR3, weight = 2, fillOpacity = 0.45,
        group       = "Migration corridor",
        label       = "Migration corridor (BBMM)"
      )
  }

  year_groups <- character(0)
  if (!is.null(res$by_year) && length(res$by_year) > 0) {
    yr_pal         <- colorFactor("Dark2", domain = names(res$by_year))
    already_hidden <- year_groups_hidden3()
    newly_hidden   <- character(0)
    for (yr in names(res$by_year)) {
      yr_sf <- res$by_year[[yr]]
      year_groups <- c(year_groups, yr)
      proxy <- proxy |>
        addPolygons(
          data      = yr_sf,
          color     = yr_pal(yr), weight = 3, fill = FALSE, dashArray = "4,4",
          group     = yr,
          label     = paste0(yr, " (", yr_sf$level * 100, "% isopleth)")
        )
      # Only hide a given year layer the first time it ever appears — once
      # the user has toggled it, later re-estimates shouldn't re-hide it.
      if (!(yr %in% already_hidden)) {
        proxy <- proxy |> hideGroup(yr)
        newly_hidden <- c(newly_hidden, yr)
      }
    }
    if (length(newly_hidden) > 0) year_groups_hidden3(c(already_hidden, newly_hidden))
  }
  prev_year_groups3(year_groups)

  proxy |>
    addLayersControl(
      position      = "topright",
      baseGroups    = c("Esri.WorldTopoMap", "Esri.WorldImagery", "Esri.WorldGrayCanvas"),
      overlayGroups = c("Boundary", "Tracks", "GPS", "Seasonal range",
                        if (has_migration) "Migration corridor", year_groups),
      options       = layersControlOptions(collapsed = FALSE)
    )

  bb <- st_bbox(rng)
  leafletProxy("map3") |> fitBounds(bb[["xmin"]], bb[["ymin"]], bb[["xmax"]], bb[["ymax"]])
})

# ---------------------------------------------------------------------------
# Summary stats / area-by-level side panel
# ---------------------------------------------------------------------------

output$summary_stats3 <- renderTable({
  d <- range_data3()
  req(nrow(d) > 0)
  tibble(
    Field = c("Caribou", "Season", "Start date", "End date", "Relocations"),
    Value = c(
      input$id3,
      input$phase3,
      as.character(as.Date(min(d$timestamp))),
      as.character(as.Date(max(d$timestamp))),
      as.character(nrow(d))
    )
  )
}, colnames = FALSE, striped = TRUE, spacing = "xs")

# One row per year for each individual, plus a "Total" row for any
# individual (or the population, when "ALL") that combines more than one
# year — e.g. 4 years for one caribou gives 5 rows: one per year and a
# total. A single-year individual gets just its one row.
output$area_by_level3 <- renderTable({
  res <- hr3()
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
    } else {
      rows[[length(rows) + 1]] <- tibble(
        Caribou = id, Year = as.character(yrs_i[1]), N = res$per_id_n[[id]],
        Area_km2 = area_km2(res$per_id[[id]])
      )
    }
  }

  if (length(ids) > 1) {
    rows[[length(rows) + 1]] <- tibble(
      Caribou = input$id3, Year = "Total", N = res$combined_n,
      Area_km2 = area_km2(res$combined)
    )
  }

  bind_rows(rows)
}, striped = TRUE, spacing = "xs")

# One panel per individual, shown only when comparing "ALL" individuals
# (a single-individual estimate has nothing to compare, so the tab stays
# empty). facet_wrap()'s default fixed scales give every panel the same
# spatial extent, so range size/shape is directly comparable across
# individuals. When more than one year was estimated, each individual's
# own yearly contours (outline only, no fill) are overlaid on their panel.
output$range_by_id_plot3 <- renderPlot({
  req(input$id3 == "ALL")
  res    <- hr3()
  per_id <- res$per_id
  req(length(per_id) > 0)

  plot_df <- bind_rows(lapply(names(per_id), function(nm) {
    per_id[[nm]] |> mutate(id = nm)
  }))

  # A single fixed fill, since only one isopleth level can be selected.
  p <- ggplot(plot_df) +
    geom_sf(fill = RANGE_COLOR3, color = "black", linewidth = 0.3, alpha = 0.5) +
    facet_wrap(~id) +
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
