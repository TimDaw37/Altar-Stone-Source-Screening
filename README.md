# Altar Stone Orcadian Basin Geochemical Screen

![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)

Open-data geochemical screen for candidate source areas of Stonehenge's Altar Stone within the Orcadian Basin, northeast Scotland. Companion code and data to:

> Daw, T. (2026). *A Multi-Element Geochemical Screen, Verified Against Bedrock Geology, for the Source of the Stonehenge Altar Stone Within the Orcadian Basin.* Preprint. [www.sarsen.org](https://www.sarsen.org)

## Background

The Altar Stone — the six-tonne sandstone megalith at the centre of Stonehenge — was shown in 2024 to originate in the Orcadian Basin of northeast Scotland rather than Wales. The sampled localities on Mainland Orkney were subsequently excluded as its specific source. That leaves most of a roughly 10,000 km² basin, stretching from the Moray Firth to Shetland, effectively unsampled.

This repository screens the whole basin using two free datasets — BGS stream-sediment barium and rubidium grids, and BGS's digital bedrock geology map — to rank locations for field follow-up. The method combines the Altar Stone's published barium signature with a Ba/Rb ratio threshold (more robust to sampling noise than either element alone), then verifies every resulting geochemical anomaly, pixel by pixel, against real mapped bedrock rather than relying on hand-drawn geographic boundaries.

The standout result is a 42.5 km² area of the East Caithness coast near Sarclet — 98.2% confirmed genuine Middle Old Red Sandstone — which independently converges with a peer-reviewed detrital zircon geochronology study (Clarke et al. 2026) that separately identified the same locality as its strongest statistical match to the Altar Stone (p = 0.96), using entirely unrelated data and methods. Full reasoning, caveats, and the negative results this screen also produced are in the paper.

## What's in this repository

| File | Contents |
|---|---|
| `screen_altar_stone_source.py` | The full pipeline: load grids, compute Ba/Rb ratio, threshold, cluster, verify against bedrock, export results. Single documented script. |
| `requirements.txt` | Python dependencies. |
| `outputs/per_cell_results.csv` | Every individual 500 m grid cell passing the composite screen (3,545 rows), with coordinates, Ba, Rb, ratio, and matched bedrock formation/age/lithology. |
| `outputs/per_cluster_results.csv` | One row per cluster (45 rows): area, centroid, mean Ba/Rb/ratio, and percentage of the cluster confirmed as genuine Old Red Sandstone. |
| `data/` | Empty — see below for what to download and place here. |

## Getting the source data

BGS data isn't redistributed here (file size, and so you always get the current release). Download and unzip these three, keeping the folder structure below — all are free, Open Government Licence, no account required:

| Folder | Download page |
|---|---|
| `data/Ba_grid/` | https://www.bgs.ac.uk/download/g-base-for-the-uk-barium_grid/ |
| `data/Rb_grid/` | https://www.bgs.ac.uk/download/g-base-for-the-uk-rubidium_grid/ |
| `data/BGS_Geology_625k_Shapefile/` | https://www.bgs.ac.uk/download/bgs-geology-625k-gis-line-and-polygon-data-shapefile-format/ |

Expected layout after unzipping:

```
data/
  Ba_grid/UK_Kriged_Ba_concentration_in_stream sediments.asc
  Rb_grid/UK_Kriged_Rb_concentration_in_stream sediments.asc
  BGS_Geology_625k_Shapefile/Bedrock/625k_V5_BEDROCK_Geology_Polygons.shp   (+ .dbf/.shx/.prj)
```

## Running it

```bash
pip install -r requirements.txt
python screen_altar_stone_source.py
```

Takes a few minutes — the per-cell bedrock spatial join is the slow step. On completion it prints the P95 Ba/Rb threshold it computed for the study extent (13.761 by default) and the largest clusters found, then writes both CSVs to `outputs/`.

## Changing the study extent

Edit `STUDY_BOX` near the top of the script. The ratio threshold is recalculated from the data within whatever extent you set on every run — it is not a hardcoded constant — so results will shift if the box changes. See Sections 2.6 and 4.3 of the paper for why the default rectangle is a stated pragmatic choice rather than a geological boundary, and for a better alternative (a dissolved outline of real Devonian polygons) not yet implemented here.

## A correction worth knowing about

An earlier draft of the paper cited incorrect cell statistics for the Sarclet cluster, sourced from a separate independent replication (by Grok/xAI, reproduced in the paper's Appendix B) that had matched the wrong cluster to that location — its sandbox lacked the GIS libraries to run the bedrock verification step itself, so it identified a candidate by rough geography rather than a confirmed spatial join. Running this script end-to-end is what caught it. The verified figures for the Sarclet cluster (`cluster_id 18` in `outputs/per_cluster_results.csv`) are: 170 cells, 42.5 km², 98.2% genuine Middle Old Red Sandstone, mean Ba 1453 ppm, mean ratio 18.2. Full account in Appendix B.6 of the paper.

## Citation

If you use this code or data, please cite:

```
Daw, T. (2026). A Multi-Element Geochemical Screen, Verified Against Bedrock Geology,
for the Source of the Stonehenge Altar Stone Within the Orcadian Basin. Preprint.
```

## License

Code and documentation in this repository: **CC BY 4.0** — see `LICENSE`. BGS source data referenced above is separately licensed under the Open Government Licence and is not redistributed in this repository.

## Contact

Tim Daw — tim.daw@gmail.com — [www.sarsen.org](https://www.sarsen.org) — [ORCID: 0000-0002-6377-2177](https://orcid.org/0000-0002-6377-2177)
