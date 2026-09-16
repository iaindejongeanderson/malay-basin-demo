# Malay Basin CO2 storage
## Introduction
This repository contains Python scripts, text files and raster files to accompany the journal article: "Regional screening of saline aquifers in the Malay Basin for CO2 storage", authored by Iain de Jonge-Anderson, Hariharan Ramachandran, Ana Widyanita, Andreas Busch, Florian Doster and Uisdean Nicholson
## Loading the raster grids
- Raster grids are presented in .bil format which can be readily imported into GIS software (e.g. ArcGIS). The projection used is Kertau 48N EPSG 24548.
- It was necessary to downsample the grids by half, from the grids presented in the manuscript, to reduce filesize

## Interactive capacity calculator
An interactive web app built from this dataset lives in [`docs/`](docs/). It maps the optimal
injection zones and runs the Monte Carlo storage-capacity calculation (Eq. 2) for any aquifer,
map-selected cluster or user-drawn area. It is plain static HTML with no build step — enable
GitHub Pages on the `/docs` folder to publish it, or open `docs/index.html` directly.
See [`docs/README.md`](docs/README.md) for details and data-provenance notes.

## Note on the raster format
The `.bil` files are GeoTIFFs rather than ESRI BIL: they carry TIFF magic and full GeoTIFF
tags (200 m cells, EPSG:24548, origin 180000 E / 820000 N), and no `.hdr` sidecars are needed.
Most GIS software will open them regardless of the extension.
