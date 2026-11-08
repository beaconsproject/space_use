# Population-level annual range as a time-weighted mixture of seasonal and
# migration-period UDs, with equal weight per animal (companion to gen_ranges.R)
# Isopleths: 95%, 99%, 99.9% and 99.99%
# Ranges are estimated for all years combined and for each year separately,
# and a map (PNG) is saved for each year x period combination.
# Requires tmap >= 4.0 and the maptiles package (for the basemap).

library(sf)
library(amt)
library(terra)
library(tidyverse)
library(tmap)

#-------------------------------------------------------------------------------
# USER INFORMATION

herd <- "little_rancheria"
gps_csv <- paste0("../../data/", herd, "/caribou.csv")
#outdir <- paste0("../../data/", herd, "/ranges/")
outdir <- paste0("output/", herd, "/")
if (!dir.exists(outdir)) dir.create(outdir)
map_dir <- paste0(outdir, "maps/")
if (!dir.exists(map_dir)) dir.create(map_dir)

min_n <- 30      # minimum fixes per animal-period (per year, for yearly ranges)
res_m <- 500     # grid resolution (m) shared by all UDs
levels <- c(0.95, 0.99, 0.999, 0.9999)
make_maps <- TRUE

#-------------------------------------------------------------------------------
# PREPARE DATA

gps_all <- read_csv(gps_csv)
trk_all <- make_track(gps_all, .x=longitude, .y=latitude, .t=timestamp,
                      id=id, year=year, ydat=yday, season=season,
                      migration=migration, crs=4326) |>
  # season is NA during migration, so fill those fixes with the migration
  # label (e.g. "spring migration") and keep them as periods of their own
  filter(!is.na(season) | !is.na(migration)) |>
  mutate(type=if_else(is.na(season), "migration", "season"),
         season=if_else(is.na(season), migration, as.character(season)),
         year=as.character(year)) |>
  transform_coords(3578)

period_type <- trk_all |> as_tibble() |> distinct(season, type) |> deframe()

# One common grid so UDs can be averaged and summed cell by cell
# (extra padding so the outer tails aren't clipped at the grid edge)
trast <- make_trast(trk_all, res=res_m, factor=2)

# Period weights = share of the year in each season/migration period, from the
# distinct days present across all years. The same weights are used for every
# year, so yearly ranges aren't distorted by partial collar coverage.
# Replace with fixed season lengths if you have them, e.g.
# w <- c(winter=151, spring=45, calving=30, summer=61, fall=78) / 365
season_days <- trk_all |>
  as_tibble() |>
  distinct(season, ydat) |>
  count(season, name="days")
w <- setNames(season_days$days, season_days$season)
w <- w / sum(w)
print(round(w, 3))

#-------------------------------------------------------------------------------
# HELPERS

# KDE for one subset, returned as a UD raster that sums to 1
ind_ud <- function(trk) {
  ud <- tryCatch(hr_kde(trk, trast=trast)$ud, error=function(e) NULL)
  if (is.null(ud)) return(NULL)
  ud / global(ud, "sum", na.rm=TRUE)[1, 1]
}

# Isopleth polygon containing `level` of the UD's probability mass
ud_isopleth <- function(ud, level) {
  v <- values(ud, mat=FALSE)
  v[is.na(v)] <- 0
  vs <- sort(v, decreasing=TRUE)
  thr <- vs[which(cumsum(vs) >= level)[1]]
  ifel(ud >= thr, 1, NA) |>
    as.polygons() |>
    st_as_sf() |>
    summarise() |>
    mutate(level=level, area_km2=round(as.numeric(st_area(geometry))/1e6, 1))
}

# Weighted mixture of period UDs, renormalized to sum to 1
mix_ud <- function(uds, w) {
  w <- w / sum(w)
  ud <- Reduce(`+`, lapply(names(w), function(s) uds[[s]] * w[[s]]))
  ud / global(ud, "sum", na.rm=TRUE)[1, 1]
}

# Annual, seasonal and seasonal-union ranges for one subset of the data
# (all years combined, or a single year). Returns polygons and the annual UD.
estimate_ranges <- function(trk, yr) {

  # Period population UDs: mean of individual UDs, equal weight per animal
  period_ud <- list()
  n_ids <- c()
  for (s in unique(trk$season)) {
    trk_s <- filter(trk, season == s)
    uds <- list()
    for (a in unique(trk_s$id)) {
      trk_a <- filter(trk_s, id == a)
      if (nrow(trk_a) >= min_n) uds[[as.character(a)]] <- ind_ud(trk_a)
    }
    uds <- compact(uds)
    if (length(uds) == 0) next
    period_ud[[s]] <- mean(rast(uds))
    n_ids[s] <- length(uds)
  }
  if (length(period_ud) == 0) {
    cat(yr, ": no periods with enough data, skipped\n")
    return(NULL)
  }

  missing <- setdiff(names(w), names(period_ud))
  cat(yr, ": ", paste0(names(n_ids), " (", n_ids, ")", collapse=", "), "\n", sep="")
  if (length(missing) > 0)
    cat("  missing periods (weights renormalized):", paste(missing, collapse=", "), "\n")

  # Annual UD = time-weighted mixture of the periods available
  wy <- w[names(period_ud)]
  wy <- wy / sum(wy)
  annual_ud <- mix_ud(period_ud, wy)

  annual_poly <- map(levels, ~ ud_isopleth(annual_ud, .x)) |>
    bind_rows() |>
    mutate(season="Annual (weighted)", type="annual",
           n_ids=n_distinct(trk$id), weight=1)

  seasonal_poly <- map(names(period_ud), function(s) {
    map(levels, ~ ud_isopleth(period_ud[[s]], .x)) |>
      bind_rows() |>
      mutate(season=s, type=period_type[[s]], n_ids=n_ids[[s]],
             weight=round(wy[[s]], 3))
  }) |> bind_rows()

  # Union of the seasonal and migration ranges at each level
  # (an extent, not a UD)
  extent_poly <- map(levels, function(l) {
    seasonal_poly |>
      filter(level == l) |>
      summarise() |>
      mutate(level=l, season="Seasonal union", type="union",
             area_km2=round(as.numeric(st_area(geometry))/1e6, 1))
  }) |> bind_rows()

  list(poly=bind_rows(annual_poly, seasonal_poly, extent_poly) |> mutate(year=yr),
       annual_ud=annual_ud)
}

#-------------------------------------------------------------------------------
# ESTIMATE RANGES: all years combined, then each year

years <- sort(unique(trk_all$year))
res <- list()
res[["All years"]] <- estimate_ranges(trk_all, "All years")
for (y in years) res[[y]] <- estimate_ranges(filter(trk_all, year == y), y)
res <- compact(res)

#-------------------------------------------------------------------------------
# SAVE

pop_ranges <- map(res, "poly") |>
  bind_rows() |>
  mutate(id="Population", herd=herd) |>
  relocate(herd, id, year, season, type, level, n_ids, weight, area_km2)

# Ranges as GeoJSON (WGS84, as the GeoJSON spec requires)
ranges_file <- paste0(outdir, "pop_ranges.geojson")
pop_ranges |>
  st_transform(4326) |>
  st_write(ranges_file, driver="GeoJSON", delete_dsn=TRUE, quiet=TRUE)

# GPS fixes as CSV with longitude/latitude (WGS84)
gps_file <- paste0(outdir, "gps.csv")
gps_sf <- st_read(paste0("../../data/", herd, "/caribou.gpkg"), "gps", quiet=TRUE)
gps_sf |>
  st_transform(4326) |>
  select(-any_of(c("longitude", "latitude"))) |>
  mutate(longitude=st_coordinates(geometry)[, 1],
         latitude=st_coordinates(geometry)[, 2],
         across(where(is.POSIXct), ~ format(.x, "%Y-%m-%d %H:%M:%S"))) |>
  st_drop_geometry() |>
  write_csv(gps_file)

# Copy to app folder
#app_dir <- paste0("../../apps/04_population_ranges/output/", herd, "/")
#if (!dir.exists(app_dir)) dir.create(app_dir, recursive=TRUE)
#file.copy(c(ranges_file, gps_file), app_dir, overwrite=TRUE)

for (y in names(res)) {
  fname <- paste0(outdir, "annual_ud_", gsub("[^A-Za-z0-9]+", "_", tolower(y)), ".tif")
  writeRaster(res[[y]]$annual_ud, fname, overwrite=TRUE)
}

pop_ranges |>
  st_drop_geometry() |>
  select(year, season, type, level, n_ids, weight, area_km2) #|>
  #print(n=Inf)

#-------------------------------------------------------------------------------
# MAPS: one PNG per year x period, with the four isopleths and GPS locations

if (make_maps) {

  gps_sf <- trk_all |> as_tibble() |> st_as_sf(coords=c("x_", "y_"), crs=3578)

  lvl_lab <- setNames(c("95%", "99%", "99.9%", "99.99%"), as.character(levels))
  lvl_col <- c("95%"="#b2182b", "99%"="#ef8a62", "99.9%"="#2166ac", "99.99%"="#542788")

  combos <- pop_ranges |>
    st_drop_geometry() |>
    filter(season != "Seasonal union") |>
    distinct(year, season)

  for (i in seq_len(nrow(combos))) {
    y <- combos$year[i]
    p <- combos$season[i]

    bnd <- pop_ranges |>
      filter(year == y, season == p) |>
      mutate(isopleth=factor(lvl_lab[as.character(level)], levels=lvl_lab)) |>
      st_boundary()

    pts <- gps_sf
    if (y != "All years") pts <- filter(pts, year == y)
    if (p != "Annual (weighted)") pts <- filter(pts, season == p)

    bb <- st_bbox(c(st_geometry(bnd), st_geometry(pts)))

    m <- tm_basemap("Esri.WorldTopoMap") +
      tm_shape(pts, bbox=bb) +
        tm_dots(fill="black", fill_alpha=0.25, size=0.05) +
      tm_shape(bnd) +
        tm_lines(col="isopleth", lwd=2,
                 col.scale=tm_scale_categorical(values=lvl_col),
                 col.legend=tm_legend(title="Isopleth")) +
      tm_title(paste0(herd, ": ", p, ", ", y,
                      " (", format(nrow(pts), big.mark=","), " fixes)")) +
      tm_scalebar(position=c("left", "bottom"))

    fname <- paste0(map_dir, herd, "_",
                    gsub("[^a-z0-9]+", "_", tolower(y)), "_",
                    gsub("[^a-z0-9]+", "_", tolower(p)), ".png")
    tmap_save(m, fname, width=8, height=8, units="in", dpi=200)
    cat("Saved", fname, "\n")
  }
}