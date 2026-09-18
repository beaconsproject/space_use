# Shared constants and helper functions used across all three dashboard
# pages (Movement data, Seasonal ranges, Seasonal paths).

DEMO_DATA_PATH <- "data/demo.csv"

# --- Shared data loader ------------------------------------------------------
# Used for both tabs, regardless of whether the active dataset is the demo
# CSV or a file the user uploaded.
load_caribou_data <- function(path) {
  d <- read_csv(path, col_types = cols(
    id        = col_character(),
    longitude = col_double(),
    latitude  = col_double(),
    timestamp = col_character(),
    .default  = col_guess()
  )) |>
    mutate(
      timestamp = as.POSIXct(timestamp, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      year      = year(timestamp)
    ) |>
    arrange(id, timestamp)

  mig_col <- intersect(
    c("migration", "Migration", "season", "Season", "migration_season"),
    names(d)
  )[1]

  if (!is.na(mig_col)) {
    d <- d |>
      rename(migration = all_of(mig_col)) |>
      mutate(migration = case_when(
        tolower(trimws(as.character(migration))) %in%
          c("spring", "spring migration") ~ "Spring",
        tolower(trimws(as.character(migration))) %in%
          c("fall", "autumn", "fall migration") ~ "Fall",
        TRUE ~ NA_character_
      ))
  } else {
    d$migration <- NA_character_
  }
  d
}

IND_COLOR      <- "#E63946"                  # color for the selected individual
OVERVIEW_COLOR <- "#B0B0B0"                  # light grey for overview markers
MIG_DASH       <- c(Spring = "8, 6", Fall = "2, 6")  # dashed / dotted for migration seasons

# Area of an sf object's geometries, in km^2, rounded to 2 decimals.
area_km2 <- function(x) round(as.numeric(units::set_units(st_area(x), "km^2")), 2)

# Union a list of per-unit sf estimates (e.g. one per year, or one per
# individual) into a single result per isopleth level. A placeholder level
# is injected for methods, like the corridor line-buffer, that don't
# produce one, so every estimate can be combined the same way.
union_by_level <- function(estimates) {
  combined_raw <- bind_rows(estimates)
  if (!"level" %in% names(combined_raw)) combined_raw$level <- NA_real_
  combined_raw |>
    group_by(level) |>
    summarise(geometry = st_union(geometry), .groups = "drop")
}

# Flatten a per-individual, per-year list of sf estimates (as produced by
# hr3()/corridor4()'s per_id_by_year) into one data frame with "id" and
# "Year" columns, ready to overlay as yearly contours in a faceted ggplot.
# Returns NULL when there's nothing to draw (e.g. only one year selected).
per_id_year_to_df <- function(per_id_by_year) {
  if (is.null(per_id_by_year) || length(per_id_by_year) == 0) return(NULL)
  df <- bind_rows(lapply(names(per_id_by_year), function(id_nm) {
    by_year <- per_id_by_year[[id_nm]]
    if (is.null(by_year)) return(NULL)
    bind_rows(lapply(names(by_year), function(yr) {
      by_year[[yr]] |> mutate(id = id_nm, Year = yr)
    }))
  }))
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df
}
