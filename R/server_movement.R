# Server logic for the "Movement data" dashboard page.


# ---------------------------------------------------------------------------
# Data source: nothing loads until the user picks the demo dataset or
# uploads their own CSV. Either can be changed again at any time.
# ---------------------------------------------------------------------------

data_source_path <- reactiveVal(NULL)

observeEvent(input$use_demo_data, {
  data_source_path(DEMO_DATA_PATH)
})

observeEvent(input$uploaded_data, {
  req(input$uploaded_data)
  data_source_path(input$uploaded_data$datapath)
})

raw_data <- reactive({
  req(data_source_path())
  load_caribou_data(data_source_path())
})

# ---------------------------------------------------------------------------
# Reactive: process the active dataset for the Movement Data map
# ---------------------------------------------------------------------------

herd_data <- reactive({
  raw <- raw_data()

  # --- 1. Build amt track (CRS 4326 → project to UTM for KDE) ---------------
  trk_amt <- raw |>
    make_track(
      .x    = longitude,
      .y    = latitude,
      .t    = timestamp,
      id     = id,
      year   = year,
      season = season,
      crs    = 4326
    ) |>
    transform_coords(crs_to = 32610)   # UTM zone 10N – adjust for your study area

  # --- 2. GPS points as sf (WGS 84) -----------------------------------------
  gps_sf <- raw |>
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
    mutate(
      id        = as.character(id),
      year      = as.integer(year),
      season    = as.character(season),
      timestamp = timestamp
    )

  # --- 3. Tracks (ordered polyline per id × year × season) ------------------
  # Groups with only one point cannot form a linestring, so filter them out first.
  trk_sf <- raw |>
    arrange(id, timestamp) |>
    group_by(id, year, season) |>
    filter(n() >= 2) |>
    summarise(
      geometry = st_sfc(st_linestring(cbind(longitude, latitude))),
      .groups  = "drop"
    ) |>
    st_as_sf(crs = 4326)

  # --- 4. KDE boundary using amt ---------------------------------------------
  # Compute a population-level KDE across all individuals (95 % isopleth).
  # hr_kde() is called directly on the track object.
  bnd_sf <- hr_kde(trk_amt, levels = 0.9999) |>
    hr_isopleths() |>
    st_transform(4326)

  list(
    gps = gps_sf,
    trk = trk_sf,
    bnd = bnd_sf
  )
})

# ---------------------------------------------------------------------------
# Update sidebar widgets when herd data changes
# ---------------------------------------------------------------------------

observeEvent(herd_data(), {
  d   <- herd_data()
  yrs <- sort(unique(d$gps$year))
  updateSelectInput(session, "id",
    choices  = as.character(sort(unique(d$gps$id))),
    selected = as.character(sort(unique(d$gps$id)))[1]
  )
  updateSelectInput(session, "season",
    choices  = c("Annual", sort(unique(d$gps$season))),
    selected = "Annual"
  )
  updateSliderInput(session, "year_range",
    min   = min(yrs),
    max   = max(yrs),
    value = c(min(yrs), max(yrs))
  )
})

# ---------------------------------------------------------------------------
# Filtered GPS reactive
# ---------------------------------------------------------------------------

gps1 <- reactive({
  req(herd_data(), input$id, input$season, input$year_range)
  gps    <- herd_data()$gps
  sel_id <- input$id
  seas   <- input$season
  yr_min <- input$year_range[1]
  yr_max <- input$year_range[2]

  out <- gps |> filter(id == sel_id, year >= yr_min, year <= yr_max)
  if (seas != "Annual") out <- out |> filter(season == seas)
  out
})

# ---------------------------------------------------------------------------
# Base map (rendered once)
# ---------------------------------------------------------------------------

output$map <- renderLeaflet({
  leaflet(options = leafletOptions(attributionControl = FALSE)) |>
    addProviderTiles("Esri.WorldImagery",    group = "Esri.WorldImagery") |>
    addProviderTiles("Esri.WorldGrayCanvas", group = "Esri.WorldGrayCanvas") |>
    addProviderTiles("Esri.WorldTopoMap",    group = "Esri.WorldTopoMap") |>
    addLayersControl(
      position      = "topright",
      baseGroups    = c("Esri.WorldTopoMap", "Esri.WorldImagery",
                        "Esri.WorldGrayCanvas"),
      overlayGroups = c("Boundary", "Tracks", "GPS", "Starting Locations"),
      options       = layersControlOptions(collapsed = FALSE)
    )
})
# This page's maps live on separate dashboard pages that start hidden.
# Shiny suspends renderLeaflet for hidden outputs by default, so the
# widget never registers in the browser until its page is first visited —
# any leafletProxy() calls issued before then (e.g. drawing the Boundary
# as soon as data loads) are silently dropped and never retried. Forcing
# the render to happen regardless of visibility fixes that.
outputOptions(output, "map", suspendWhenHidden = FALSE)

# ---------------------------------------------------------------------------
# Static mode: KDE boundary, shown/hidden and re-flown based on map_mode
# ---------------------------------------------------------------------------

observeEvent(herd_data(), {
  bnd   <- herd_data()$bnd
  proxy <- leafletProxy("map") |>
    clearGroup("Boundary") |>
    addPolygons(
      data = bnd,
      fillOpacity = 0, stroke = TRUE,
      color = "black", weight = 2,
      group = "Boundary"
    )

  if (identical(input$map_mode, "static")) {
    proxy |> showGroup("Boundary")
    bb <- st_bbox(bnd)
    leafletProxy("map") |>
      fitBounds(bb[["xmin"]], bb[["ymin"]], bb[["xmax"]], bb[["ymax"]])
  } else {
    proxy |> hideGroup("Boundary")
  }
})

# Update Tracks / GPS layers
observe({
  req(gps1(), input$id, input$season, input$year_range)
  sel_id <- input$id
  seas   <- input$season
  yr_min <- input$year_range[1]
  yr_max <- input$year_range[2]

  trk <- herd_data()$trk |> filter(id == sel_id, year >= yr_min, year <= yr_max)
  if (seas != "Annual") trk <- trk |> filter(season == seas)

  proxy <- leafletProxy("map") |>
    clearGroup("GPS") |>
    clearGroup("Tracks") |>
    addPolylines(
      data    = trk,
      color   = "black", weight = 1, opacity = 0.7,
      group   = "Tracks"
    ) |>
    addCircles(
      data        = gps1(),
      color       = "red", fillColor = "red",
      fillOpacity = input$alpha,
      weight      = input$size, radius = 1,
      group       = "GPS"
    )

  if (identical(input$map_mode, "static")) {
    proxy |> showGroup("Tracks") |> showGroup("GPS")
  } else {
    proxy |> hideGroup("Tracks") |> hideGroup("GPS")
  }
})

# ---------------------------------------------------------------------------
# Map view toggle: show/hide static vs. animation layers on the shared map
# ---------------------------------------------------------------------------

observeEvent(input$map_mode, {
  proxy <- leafletProxy("map")

  if (identical(input$map_mode, "static")) {
    proxy |>
      showGroup("Boundary") |> showGroup("Tracks") |> showGroup("GPS") |>
      hideGroup("Starting Locations") |> hideGroup("path_ind1") |> hideGroup("marker_ind1")

    bnd <- herd_data()$bnd
    if (!is.null(bnd) && nrow(bnd) > 0) {
      bb <- st_bbox(bnd)
      proxy |> fitBounds(bb[["xmin"]], bb[["ymin"]], bb[["xmax"]], bb[["ymax"]])
    }
  } else {
    proxy |>
      hideGroup("Boundary") |> hideGroup("Tracks") |> hideGroup("GPS") |>
      showGroup("Starting Locations") |> showGroup("path_ind1") |> showGroup("marker_ind1")

    all_pts <- filtered_data()
    if (is.null(all_pts) || nrow(all_pts) == 0) all_pts <- year_first_fixes()
    if (!is.null(all_pts) && nrow(all_pts) > 0) {
      proxy |> fitBounds(min(all_pts$longitude), min(all_pts$latitude),
                          max(all_pts$longitude), max(all_pts$latitude))
    }
  }
}, ignoreNULL = FALSE)

output$summary_stats <- renderUI({
  x <- gps1()
  req(nrow(x) > 0)
  min_date <- as.character(as.Date(min(x$timestamp)))
  max_date <- as.character(as.Date(max(x$timestamp)))
  n        <- nrow(x)

  tagList(
    tags$h6(style = paste0("color:", IND_COLOR, ";margin-bottom:4px;"),
            paste0("Selected - ", input$id)),
    tags$p(style = "font-size:0.85rem;line-height:1.8;margin:0;", HTML(paste0(
      "<b>Season:</b> ",      input$season, "<br>",
      "<b>Years:</b> ",       input$year_range[1], "-", input$year_range[2], "<br>",
      "<b>Start date:</b> ",  min_date, "<br>",
      "<b>End date:</b> ",    max_date, "<br>",
      "<b>Relocations:</b> ", n
    )))
  )
})

output$sampling_duration_plot <- renderPlot({
  req(herd_data())
  herd_data()$gps |>
    st_drop_geometry() |>
    group_by(id) |>
    summarise(start = min(year), end = max(year), .groups = "drop") |>
    mutate(id = fct_rev(factor(id))) |>
    ggplot(aes(y = id)) +
    geom_segment(aes(x = start, xend = end, yend = id)) +
    geom_point(aes(x = start)) +
    geom_point(aes(x = end)) +
    scale_x_continuous(breaks = scales::breaks_width(1)) +
    labs(x = "Year", y = NULL) +
    theme_minimal()
})

output$individuals_by_season_table <- renderTable({
  req(herd_data())
  row_order <- c("Early winter", "Late winter", "Calving", "Summer",
                 "Fall rut", "Spring migration", "Fall migration")

  herd_data()$gps |>
    st_drop_geometry() |>
    group_by(year, phase) |>
    summarise(n_ids = n_distinct(id), .groups = "drop") |>
    pivot_wider(names_from = year, values_from = n_ids) |>
    rename(Phase = phase) |>
    mutate(Phase = factor(Phase, levels = row_order)) |>
    arrange(Phase)
}, striped = TRUE, spacing = "xs", na = "")


# ---------------------------------------------------------------------------
# Animation map (individual movement playback)
# Uses the same Caribou / Years widgets as the Mapview map.
# ---------------------------------------------------------------------------

# First fix per individual within the selected year range (overview overlay)
year_first_fixes <- reactive({
  req(input$year_range)
  yr_min <- input$year_range[1]
  yr_max <- input$year_range[2]
  raw_data() |>
    filter(year >= yr_min, year <= yr_max) |>
    group_by(id) |>
    slice_min(timestamp, n = 1, with_ties = FALSE) |>
    ungroup()
})

# ── Per-individual animation data ─────────────────────────────────────────
filtered_data <- reactive({
  req(input$id, input$year_range)
  yr_min <- input$year_range[1]
  yr_max <- input$year_range[2]
  raw_data() |>
    filter(id == input$id, year >= yr_min, year <= yr_max) |>
    arrange(timestamp)
})

# Calendar timeline for the selected individual
shared_times <- reactive({
  sort(unique(filtered_data()$timestamp))
})

observe({
  n <- length(shared_times())
  updateSliderInput(session, "animation_index",
                    min = 1, max = max(n, 1), value = 1)
})

current_time <- reactive({
  times <- shared_times()
  idx   <- min(input$animation_index, length(times))
  times[idx]
})

output$current_time_label <- renderUI({
  tags$p(style = "font-size:0.85rem;font-weight:600;margin-bottom:2px;",
         paste0("Current time: ", format(current_time(), "%Y-%m-%d %H:%M")))
})

path_up_to <- function(df) {
  ct <- current_time()
  if (is.null(df) || nrow(df) == 0) return(NULL)
  d <- df |> filter(timestamp <= ct)
  if (nrow(d) == 0) return(NULL)
  d
}

animated_data <- reactive({ path_up_to(filtered_data()) })

# ── Starting-location markers (animation mode overview) ───────────────────
observe({
  ff <- year_first_fixes()
  req(nrow(ff) > 0)

  proxy <- leafletProxy("map") |>
    clearGroup("Starting Locations") |>
    addCircleMarkers(
      data        = ff,
      lng         = ~longitude,
      lat         = ~latitude,
      group       = "Starting Locations",
      radius      = 7,
      color       = "#ffffff",
      weight      = 2,
      fillColor   = OVERVIEW_COLOR,
      fillOpacity = 0.85,
      label       = ~as.character(id),
      labelOptions = labelOptions(
        style     = list("font-weight" = "bold", "font-size" = "13px"),
        direction = "top",
        offset    = c(0, -10)
      ),
      popup = ~paste0(
        "<b>ID:</b> ", id, "<br>",
        "<b>First fix (", input$year_range[1], "-", input$year_range[2], "):</b> ",
        format(timestamp, "%Y-%m-%d %H:%M"), "<br>",
        "<b>Lat/Lon:</b> ", round(latitude, 4), " / ", round(longitude, 4)
      )
    )

  if (identical(input$map_mode, "animation")) {
    proxy |> showGroup("Starting Locations")
  } else {
    proxy |> hideGroup("Starting Locations")
  }
})

# ── Animated paths and current-position markers ───────────────────────────
observe({
  pts <- animated_data()

  proxy <- leafletProxy("map") |>
    clearGroup("path_ind1") |> clearGroup("marker_ind1")

  if (is.null(pts) || nrow(pts) == 0) return(invisible(NULL))

  # Assign a segment id that increments whenever migration value changes
  pts <- pts |>
    mutate(
      mig_val  = ifelse(is.na(migration) | !(migration %in% c("Spring", "Fall")),
                        "none", migration),
      seg_id   = cumsum(mig_val != lag(mig_val, default = first(mig_val)))
    )

  seg_ids <- unique(pts$seg_id)
  for (sid in seg_ids) {
    seg      <- pts |> filter(seg_id == sid)
    seg_dash <- switch(seg$mig_val[1],
                       Spring = MIG_DASH[["Spring"]],
                       Fall   = MIG_DASH[["Fall"]],
                       NULL)
    proxy <- proxy |>
      addPolylines(data = seg, lng = ~longitude, lat = ~latitude,
                   group = "path_ind1", color = IND_COLOR,
                   weight = input$line_width, opacity = 0.85,
                   dashArray = seg_dash)
  }

  cur <- pts |> slice_tail(n = 1)
  proxy <- proxy |>
    addCircleMarkers(data = cur, lng = ~longitude, lat = ~latitude,
                     radius = 6, color = "#ffffff", weight = 2,
                     fillColor = IND_COLOR, fillOpacity = 1,
                     group = "marker_ind1",
                     popup = ~paste0("<b>ID:</b> ", id, "<br><b>Time:</b> ",
                                     format(timestamp, "%Y-%m-%d %H:%M")))

  if (identical(input$map_mode, "animation")) {
    proxy |> showGroup("path_ind1") |> showGroup("marker_ind1")
  } else {
    proxy |> hideGroup("path_ind1") |> hideGroup("marker_ind1")
  }
})

# ── Animation details (current location + legend, one combined output ────
# so there's no card-level gap between what used to be four separate
# top-level elements) ───────────────────────────────────────────────────
output$animation_details <- renderUI({
  ct  <- current_time()
  pts <- animated_data()
  id1 <- input$id

  mig_label <- function(val) {
    if (is.null(val) || is.na(val) || trimws(val) == "") "" else val
  }

  make_block <- function(pts, d_all, color, label) {
    if (is.null(pts) || nrow(pts) == 0) {
      start <- if (!is.null(d_all) && nrow(d_all) > 0) min(d_all$timestamp) else NA
      end   <- if (!is.null(d_all) && nrow(d_all) > 0) max(d_all$timestamp) else NA

      status <- if (!is.null(d_all) && nrow(d_all) > 0 && !is.na(start) && ct < start)
        paste0("Not yet active (starts ", format(start, "%b %d"), ")")
      else if (!is.null(d_all) && nrow(d_all) > 0 && !is.na(end) && ct > end)
        paste0("No longer tracked (last seen ", format(end, "%b %d"), ")")
      else
        "No data available"

      return(tagList(
        tags$h6(style = paste0("color:", color, ";margin-bottom:4px;"), label),
        tags$p(style = "color:#888;font-size:0.85rem;font-style:italic;margin:0;", status)
      ))
    }

    cur       <- pts |> slice_tail(n = 1)
    n_shown   <- nrow(pts)
    n_total   <- nrow(d_all)
    phase     <- mig_label(if ("phase" %in% names(cur)) cur$phase else NA)
    phase_disp <- if (nchar(phase) > 0) phase else ""
    phase_html <- paste0("<b>Phase:</b> ", phase_disp, "<br>")

    tagList(
      tags$h6(style = paste0("color:", color, ";margin-bottom:4px;"), label),
      tags$p(style = "font-size:0.85rem;line-height:1.8;margin:0;", HTML(paste0(
        "<b>Points shown:</b> ", n_shown, " / ", n_total, "<br>",
        "<b>Current time:</b> ", format(cur$timestamp, "%Y-%m-%d %H:%M"), "<br>",
        "<b>Latitude:</b> ",  round(cur$latitude,  4), "<br>",
        "<b>Longitude:</b> ", round(cur$longitude, 4), "<br>",
        phase_html
      )))
    )
  }

  location_block <- make_block(pts, filtered_data(), IND_COLOR,
                                paste0("Selected - ", id1))

  swatch <- function(color, shape = "line", label, dash = "solid") {
    icon_style <- if (shape == "circle")
      paste0("display:inline-block;width:12px;height:12px;border-radius:50%;",
             "background:", color, ";margin-right:6px;border:2px solid #fff;",
             "box-shadow:0 0 0 1px #999;")
    else {
      border_style <- switch(dash,
        dashed = "dashed",
        dotted = "dotted",
        "solid"
      )
      paste0("display:inline-block;width:24px;height:0;border-top:4px ",
             border_style, " ", color, ";margin-right:6px;")
    }

    tags$span(
      style = "display:inline-flex;align-items:center;margin-right:14px;font-size:0.82rem;",
      tags$span(style = icon_style, ""),
      label
    )
  }

  legend_items <- tagList(
    swatch(OVERVIEW_COLOR, "circle",
           paste0("All individuals (", input$year_range[1], "-", input$year_range[2], " start)")),
    swatch(IND_COLOR, "line",
           paste0("Selected: ", id1))
  )
  legend_items <- tagList(legend_items,
                   swatch("#666", "line", "Spring migration", dash = "dashed"),
                   swatch("#666", "line", "Fall migration",   dash = "dotted"))

  tagList(
    location_block,
    tags$br(),
    tags$hr(style = "margin:4px 0;"),
    tags$br(),
    tags$p(style = "font-weight:600;margin:0 0 2px;", "Map legend"),
    div(style = "padding:0 10px 3px;display:flex;flex-wrap:wrap;gap:4px;", legend_items)
  )
})
