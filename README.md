# gmpas demo data

Value tiles for the [gmpas demo](https://nmathewa.github.io/gmpas-lib/): 20
hourly steps of six fields from the MPAS 3.75 km global DYAMOND2 run
(41,943,042 cells), sampled by gmpas from the native mesh.

The page is the gmpas viewer itself, baked in gmpas-lib (`docs/demo/bake.py
--tiles`); its `shim.js` reads `manifest.json` and `tiles/` from this
repository's GitHub Pages site, `https://nmathewa.github.io/gmpas-demo-data/`.

Each tile is a 256-pixel greyscale PNG holding values, not colours: byte 0 is
no value, byte k is `lo + (k - 1) * (hi - lo) / 254`, with `lo` and `hi` the
field's range in `manifest.json`. Every pixel takes the nearest cell to its
centre. The page colours them, so colormaps, ranges and clicks work; values
are rounded to 1/254 of each field's range.

## Rebuilding the tiles

Needs gmpas, about 20 GB of scratch for the raw fields and the mesh file
`x1.41943042.grid.nc` (8.6 GB; its gmpas cache is about 13 GB).

```bash
python render/fetch.py RAW_DIR                        # 6 fields x 20 steps from NCAR GDEX
python render/render.py RAW_DIR x1.41943042.grid.nc   # writes site/tiles and site/manifest.json
```

Fields, colour maps, ranges and times are set in `render/fields.json`.

## Data

Judt, F. (2026). *MPAS Global Storm-Resolving Model Simulations for the
DYAMOND Intercomparison* [Dataset]. NSF National Center for Atmospheric
Research. https://doi.org/10.5065/MWEE-Y657

