# What makes a CO₂ storage zone "optimal"?

An interactive explorer for the **screening cut-offs** that define optimal CO₂ injection zones in
the Malay Basin. Move a cut-off, and the zones redraw from the underlying property grids in real
time — along with the area, the criteria breakdown and the storage capacity.

Built from the dataset accompanying:

> de Jonge-Anderson, I., Ramachandran, H., Widyanita, A., Busch, A., Doster, F. & Nicholson, U.
> (2025) *Regional screening of saline aquifers in the Malay Basin for CO₂ storage.*
> International Journal of Greenhouse Gas Control **143**, 104347.
> https://doi.org/10.1016/j.ijggc.2025.104347

## Publishing on GitHub Pages

Plain static HTML, no build step, no dependencies.

1. **Settings → Pages → Build and deployment**
2. Source: *Deploy from a branch*; Branch: `main`, folder: **`/docs`**

It also runs by opening `docs/index.html` straight off disk — the rasters are embedded as data
URIs, so there are no `fetch` or CORS constraints.

## What you can do with it

- **Move any cut-off** — porosity, CO₂ density, fault setback, phase, overpressure — separately for
  the optimal and sub-optimal classes, and watch the map reclassify.
- **See what is binding.** The breakdown chart shows how much of each aquifer's footprint each
  criterion rejects on its own. This is where the geology shows up: Group B is 97 % rejected by
  CO₂ phase because it is too shallow for CO₂ to be supercritical, while Group K is 69 % rejected
  by overpressure and 49 % by porosity because it is too deep.
- **Sweep one cut-off** and see area and capacity respond across its whole range — a live,
  generalised version of Fig. 13 in the paper.
- **Watch capacity follow.** Area, porosity and CO₂ density are re-measured inside whatever zone
  you have just defined and fed into Eq. 2, so a cut-off change propagates all the way to gigatonnes.

## Files

| File | Purpose |
|------|---------|
| `index.html` | Page shell, styles, layout |
| `app.js` | Classification, map, charts, projection — all application logic |
| `data/crit_*.js` | Generated. One packed raster per aquifer, base64 PNG, loaded on demand |
| `data/crit_*.png` | The same rasters as standalone files |

## How the data is packed

Each aquifer is one RGB PNG on a 400 m grid:

| Channel | Carries |
|---------|---------|
| R | porosity, 0.2 % per step |
| G | CO₂ density, 4 kg/m³ per step |
| B | bits 0–1 fault setback class (0–3), bit 2 gas phase, bit 3 overpressured, bit 4 in-footprint |

The overpressure flag is derived: the archived pressure grids are cleanly bimodal at 10 MPa/km
(hydrostatic) and 20 MPa/km (overpressured), so the flag is set where the gradient exceeds
15 MPa/km. Files are 110–310 KB each and load one aquifer at a time; reclassifying ~1 M cells takes
about 15 ms, which is what makes the sliders feel live.

## Fidelity

At the published cut-offs the live classification agrees with the archived optimal-zone rasters
cell for cell:

| Group | B | D | E | F | H | I | J | K |
|---|---|---|---|---|---|---|---|---|
| agreement | 99.8 % | 94.3 % | 94.2 % | 95.9 % | 95.5 % | 98.0 % | 97.5 % | 97.9 % |

The residual sits on zone edges, where resampling 200 m → 400 m and quantising porosity and density
move a boundary cell. Areas land within a few per cent of the archived rasters. Note the app does
**not** apply the paper's DBSCAN clustering and 100 km² minimum-area cut-off, so its areas run
slightly higher than Table 2.

## Two things the archived data says that the summary table does not

**1. Overpressure is not excluded from the sub-optimal class.** The summary table marks overpressure
"excluded" for both classes, but `make_optimal_zone_map.py` tests it only on the optimal branch, and
the archived rasters agree with the code — applying it to both drops agreement for Group I from 97 %
to 78 %. The app ships with that checkbox unticked for sub-optimal, reproducing what was run.

**2. The fault rasters are georeferenced to a different extent than they were used at.** Every fault
grid is archived at the full 2000×2050 basin extent, but for Groups D, E, F and H the other property
grids are clipped to a smaller window with a different origin. The published zone maps match a fault
lookup done by *array index* rather than by map position: 94–96 % agreement that way, against 86–92 %
when the faults are placed by their own georeferencing. Effectively the fault control sits about
3.5 km east and 5.3 km north of where its header says, for those four aquifers only. The app follows
the published behaviour so that it reproduces the paper. Groups B, I, J and K are unaffected because
all their grids share one extent.

## A note on the raster format

Despite the `.bil` extension these are **GeoTIFFs** — TIFF magic and full GeoTIFF tags (200 m cells,
EPSG:24548 Kertau 1968 / UTM 48N), with no `.hdr` sidecars needed. Their TIFF strips are not stored
in row order, so a reader must honour the strip offsets rather than assume a contiguous raster.

## Capacity assumptions

The paper gives the means of storage efficiency (E = 2 %) and irreducible water saturation
(S_wirr = 27 %) but not their standard deviations, which the Monte Carlo needs. This app uses
**E = 0.02 ± 0.01** and **S_wirr = 0.27 ± 0.05**; the efficiency spread was chosen because it
reproduces the published P10–P50–P90 most closely across all seven aquifers. Percentiles follow the
paper's exceedance convention, so **P10 is the high estimate**.

## Licence

The underlying paper is open access under CC BY 4.0. Screening estimates, not a substitute for
site-specific evaluation.
