# Migration corridor workflow

This is a well-established pipeline, most directly formalized by Sawyer et al. (2009) and now implemented in tools like the Wyoming Migration Initiative's Migration Mapper. It's a bottom-up, hierarchical process - each individual gets estimated first, then combined - and it deliberately avoids simply pooling all fixes across individuals and years, for the same reason we discussed with population-level aKDE.

**Step 1 - Fit BBMM separately for each individual × year × migration event**

Each spring/fall migration sequence for each individual in each year gets its own BBMM run - not one BBMM per individual pooling multiple years, and not one pooling multiple individuals. Practical details from the applied literature:

- All events are estimated on the same spatial grid/resolution (commonly 50–100 m) covering the full study extent, so outputs are directly comparable/combinable later.
- Gaps of roughly 8 hours or more between successive fixes are typically excluded from bridging (Montana FWP pronghorn report, following Sawyer et al. 2009) - you don't want to build a Brownian bridge across a stopover rest period or a data gap as if the animal moved in a straight diffusive line.
- The motion-variance parameter is estimated per sequence (maximum likelihood, per Horne et al. 2007), and some recent implementations explicitly censor sequences with extreme motion-variance estimates, since these can disproportionately inflate the resulting population-level corridor (McKee et al. 2024, citing Merkle et al. 2023).

**Step 2 - Normalize each event's output to a proper probability surface**

Each individual-event UD is rescaled so its raster sums to 1, keeping every event on equal probabilistic footing regardless of relocation count or migration duration.

**Step 3 - Average within individual, across years, before going to the population**

This is the step that most directly answers "how are they combined," and it mirrors the nested individual-then-population averaging from the aKDE workflow: UDs for a given individual's spring (or fall) migrations are averaged across all years into a single individual-level migration UD, then rescaled to sum to 1 (Montana FWP pronghorn report, following Sawyer et al. 2009). This prevents an individual tracked for 6 years from dominating the population signal relative to one tracked for 2.

**Step 4 - Extract an isopleth footprint per individual**

A migration route footprint is defined for each individual, typically as the 99% isopleth of that individual's averaged UD, converted to a binary raster (inside/outside the footprint) (Montana FWP pronghorn report; Kauffman et al. 2020; Sawyer et al. 2009).

**Step 5 - Stack footprints across individuals**

All individual footprints for the population are stacked, and use tiers (low/moderate/high) are defined by how many individuals' footprints include each pixel (Montana FWP pronghorn report). Combining individual isopleths this way, given a reasonably large sample (often cited as ≥20 individuals), is what delineates the corridor and its intensity of use at the population level (McKee et al. 2024).

**Step 6 - Convert the count raster to the final corridor product**

The stacked count (or proportion-of-population) raster is classified into use tiers and/or contoured at a chosen threshold to produce the final corridor polygon(s) - e.g., "used by ≥25% of the sampled population."

**Optional refinement - separating corridors from stopovers**

Sawyer et al.'s original framework additionally distinguishes route segments functioning as stopover sites from those used primarily as movement corridors, using the time-weighted intensity within the BBMM output (cells where animals linger vs. cells rapidly transited) - useful if you want stopover sites flagged separately, since they typically warrant different protection than a narrow travel corridor.

**A few caveats worth carrying into your analysis**

- Fitting BBMM to sparse location data violates the assumption of conditional random movement between fixes, which can overestimate corridor width (USGS) - worth checking if your Spring/Fall migration phases have coarser fix rates than the seasonal ranges.
- Sample-size guidelines exist specifically for how many tracked individuals are needed for a reliable population-level corridor or distribution map (Beaupre et al. 2025) - worth checking against, given your dataset has 5 individuals in the sample (40 in the real dataset).
- A simpler, line-buffer-based corridor delineation method has been proposed as an alternative to full BBMM isopleth stacking, for cases where transparency/reproducibility matters more than the full probabilistic surface (Merkle et al. 2023).

**References**

- Sawyer, H., Kauffman, M. J., Nielson, R. M., & Horne, J. S. (2009). Identifying and prioritizing ungulate migration routes for landscape-level conservation. *Ecological Applications*, 19(8), 2016–2025.
- Horne, J. S., Garton, E. O., Krone, S. M., & Lewis, J. S. (2007). Analyzing animal movements using Brownian bridges. *Ecology*, 88(9), 2354–2363.
- Merkle, J. A., et al. (2023). Conserving habitat for migratory ungulates: How wide is a migration corridor? *Journal of Applied Ecology*, 60(9), 1743–1754.
- McKee, C., et al. (2024). Estimating ungulate migration corridors from sparse movement data. *Ecosphere*, 15(9), e4983.
- Beaupre, ... et al. (2025). Sample size guidelines for mapping migration corridors and population distributions using tracking data. *Ecological Solutions and Evidence*.
- Kauffman, M., et al. (2020). *Ungulate migrations of the Western United States, Volume 1.* U.S. Geological Survey, Scientific Investigations Report 2020-5101.
- Montana Fish, Wildlife & Parks. Final Report: Grant #F19AP00029, Madison Valley Pronghorn Movements Study.
- U.S. Geological Survey. Estimating ungulate migration corridors from sparse movement data (publication summary).
- Wyoming Migration Initiative. Migration Mapper - methodology and user guide.
