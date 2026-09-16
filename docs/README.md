# Malay Basin CO₂ Storage Capacity Calculator

An interactive web app for calculating probabilistic CO₂ storage capacity in Malay Basin
saline aquifers, built from the dataset accompanying:

> de Jonge-Anderson, I., Ramachandran, H., Widyanita, A., Busch, A., Doster, F. & Nicholson, U.
> (2025) *Regional screening of saline aquifers in the Malay Basin for CO₂ storage.*
> International Journal of Greenhouse Gas Control **143**, 104347.
> https://doi.org/10.1016/j.ijggc.2025.104347

## Publishing on GitHub Pages

The app is plain static HTML with no build step and no dependencies.

1. Push this repository to GitHub.
2. **Settings → Pages → Build and deployment**
   - Source: *Deploy from a branch*
   - Branch: `main`, folder: **`/docs`**
3. The app appears at `https://<user>.github.io/<repo>/` within a minute or so.

It also runs by opening `docs/index.html` directly from disk — the raster data is
embedded as data URIs in `data.js`, so there are no `fetch` or CORS constraints.

## Files

| File | Purpose |
|------|---------|
| `index.html` | Page shell, styles, layout |
| `app.js` | Calculator, map, charts, projection — all application logic |
| `data.js` | Generated. Optimal-zone rasters, base64 PNG + grid georeferencing |
| `data/zones_*.png` | The same rasters as standalone files, if you want them separately |
| `data/grids.json` | Grid extents and per-cluster areas |

## What it does

- **Map.** The eight optimal-zone rasters from `9 Optimal zone grids/` decoded and drawn at
  native 200 m resolution, with a lat/lon graticule computed from the Kertau 1968 / UTM 48N
  projection (EPSG:24548). Pan, zoom, and a live class/coordinate readout.
- **Area selection.** Use the published Table 2 area, click a single connected optimal zone,
  or drag a box to measure optimal ground within any sub-area. Selected area feeds straight
  into the capacity equation.
- **Monte Carlo capacity.** Eq. 2 of the paper,
  `M_CO2 = A · h · NTG · φ · (1 − S_wirr) · E · ρ_CO2`, with every term drawn from a normal
  distribution truncated to physically valid bounds. Every parameter is editable.
- **Sensitivity.** Tornado plot of the swing in capacity from each variable's 10th to 90th
  percentile.
- **Basin summary.** All eight aquifers, with both the column-sum total (as the paper reports)
  and the correct per-trial aggregate. CSV export.

## Data provenance and two things worth knowing

**The `.bil` files are GeoTIFFs.** Despite the extension there are no ESRI `.hdr` sidecars;
the files carry `II*\0` TIFF magic and full GeoTIFF tags (2000×2050 max, 200 m cells,
origin 180000 E / 820000 N, EPSG:24548). They are read here by parsing those tags directly.
Note that the TIFF strips are *not* stored in row order, so strip offsets must be honoured
rather than assuming a contiguous raster.

**Two published values could not be reproduced from the archived grids.** Measuring optimal
zone pixels directly gives D 3,475 · E 9,118 · F 4,467 · H 9,524 · I 24,402 · J 15,465 ·
K 12,597 km², against Table 2's 3,348 · 13,894 · 18,108 · 22,290 · 24,924 · 12,898 ·
10,643 km². D and I agree; E, F and H are markedly smaller in the grids and J and K larger.
Porosity and CO₂ density averaged inside those same grid zones reproduce Table 2 *exactly*
for all seven aquifers, so the decoding is sound and the difference lies in zone extent
alone — most plausibly the halved grid resolution noted in the dataset README, or rasters
regenerated after the table was finalised. **The calculator defaults to the published Table 2
areas.** Grid-derived areas are offered separately and used for any map selection.

**Two standard deviations are assumptions, not published values.** The paper gives the means
of storage efficiency (E = 2 %) and irreducible water saturation (S_wirr = 27 %) but not
their spreads, which a Monte Carlo requires. This app defaults to **E = 0.02 ± 0.01** and
**S_wirr = 0.27 ± 0.05**. The efficiency spread was chosen by fitting — it reproduces the
published P10–P50–P90 most closely across all seven aquifers; S_wirr's spread barely matters.
With these defaults the app returns a basin total of P10 23.6 · P50 7.8 · P90 1.45 Gt
against the published 24.5 · 7.6 · 1.44 Gt. Both values are editable.

## Regenerating `data.js`

`data.js` is generated from the `.bil` grids. If the grids change, re-run the export
described in `docs/data/grids.json`'s provenance — the encoding is one 8-bit greyscale PNG
per aquifer where the pixel value is `0` no data, `1` non-viable, `2` sub-optimal,
`3` optimal with no cluster label, and `4+` optimal with cluster id `value − 3`.

## Licence

The underlying paper is open access under CC BY 4.0. Capacity figures are regional screening
estimates and are not a substitute for site-specific evaluation.
