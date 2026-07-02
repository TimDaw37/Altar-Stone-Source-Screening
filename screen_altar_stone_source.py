#!/usr/bin/env python3
"""
screen_altar_stone_source.py

Reproduces the Ba/Rb geochemical screen and bedrock verification described in:
Daw, T. (2026). A Multi-Element Geochemical Screen, Verified Against Bedrock
Geology, for the Source of the Stonehenge Altar Stone Within the Orcadian Basin.

WHAT THIS DOES
--------------
1. Loads BGS G-BASE barium and rubidium stream-sediment grids (500 m, whole GB).
2. Computes a Ba/Rb ratio grid.
3. Clips to a study extent (default: Moray Firth to Shetland; see STUDY_BOX below).
4. Screens cells on: Ba >= BA_FLOOR_PPM  AND  (Ba/Rb) >= P95 of the ratio
   within the study extent (P95 is computed fresh from the data each run --
   it is NOT a fixed constant, and will change if you change the study extent).
5. Groups passing cells into 8-connected clusters, discards clusters smaller
   than MIN_CLUSTER_CELLS.
6. Runs a point-in-polygon spatial join of EVERY passing cell (not just cluster
   centroids) against the BGS Geology 625k bedrock polygon layer, classifying
   each cell as genuine Old Red Sandstone or not (see is_genuine_ors()).
7. Writes two CSVs to outputs/:
      per_cell_results.csv    -- every passing grid cell, with its cluster ID,
                                  coordinates, Ba, Rb, ratio, and bedrock match
      per_cluster_results.csv -- one row per cluster: area, centroid, mean
                                  Ba/Rb/ratio, % of cluster on genuine ORS,
                                  and the dominant matching formation name

DATA YOU NEED TO SUPPLY (not included in this repository -- see README.md)
----------------------------------------------------------------------
Place the following, unzipped, under ./data/ before running:

  data/Ba_grid/UK_Kriged_Ba_concentration_in_stream sediments.asc
  data/Rb_grid/UK_Kriged_Rb_concentration_in_stream sediments.asc
  data/BGS_Geology_625k_Shapefile/Bedrock/625k_V5_BEDROCK_Geology_Polygons.shp
      (plus the accompanying .dbf/.shx/.prj files from the same download)

All three are free downloads under the Open Government Licence from:
  https://www.bgs.ac.uk/download/g-base-for-the-uk-barium_grid/
  https://www.bgs.ac.uk/download/g-base-for-the-uk-rubidium_grid/
  https://www.bgs.ac.uk/download/bgs-geology-625k-gis-line-and-polygon-data-shapefile-format/

DEPENDENCIES
------------
numpy, scipy, pyproj, geopandas, shapely  (see requirements.txt)
"""

import json
import numpy as np
import pandas as pd
import geopandas as gpd
from scipy import ndimage
from shapely.geometry import Point
from pyproj import Transformer

# ----------------------------------------------------------------------
# CONFIGURATION -- change these to rerun with a different extent/threshold
# ----------------------------------------------------------------------

BA_PATH = "data/Ba_grid/UK_Kriged_Ba_concentration_in_stream sediments.asc"
RB_PATH = "data/Rb_grid/UK_Kriged_Rb_concentration_in_stream sediments.asc"
BEDROCK_SHP = "data/BGS_Geology_625k_Shapefile/Bedrock/625k_V5_BEDROCK_Geology_Polygons.shp"

# Study extent in OSGB36 easting/northing (metres). Default covers the
# Moray Firth to Shetland. This is a pragmatic rectangle, not a geological
# boundary -- see Section 2.6 / 4.3 of the paper for the edge-effect risk
# this carries, and change it if you have a better basin outline.
STUDY_BOX = dict(e_min=225000, e_max=480000, n_min=790000, n_max=1219700)

BA_FLOOR_PPM = 1025          # Altar Stone barium floor (Bevins et al. 2024)
RATIO_PERCENTILE = 95        # threshold = this percentile of Ba/Rb within STUDY_BOX
MIN_CLUSTER_CELLS = 3        # 3 cells = 0.75 km^2 at 500 m resolution; smaller discarded as noise
CELL_SIZE_M = 500

# Lithology description (RCS_D) substrings that mean "Devonian but NOT
# sedimentary" -- igneous intrusions, lavas, and metamorphic rock of the
# same age, which must be excluded even though they share the Devonian tag.
NON_SEDIMENTARY_TERMS = [
    "IGNEOUS", "LAVA", "TUFF", "SCHIST", "ULTRAMAFIT",
    "PYROCLASTIC", "METABRECCIA", "FELSIC-ROCK", "GNEISS",
]


def read_asc_header(path, n_header_lines=6):
    """Parse an ESRI ASCII grid header into a dict of {KEY: value}."""
    header = {}
    with open(path, "r") as f:
        for _ in range(n_header_lines):
            key, val = f.readline().split()
            header[key.upper()] = float(val)
    return header


def load_grid(path):
    """Load an ESRI ASCII grid, returning (2D numpy array, header dict)."""
    header = read_asc_header(path)
    data = np.loadtxt(path, skiprows=6, dtype=np.float32)
    nrows, ncols = int(header["NROWS"]), int(header["NCOLS"])
    return data.reshape(nrows, ncols), header


def box_to_indices(header, e_min, e_max, n_min, n_max):
    """Convert an OSGB36 bounding box into row/col index bounds for the grid."""
    xll, yll, cs = header["XLLCORNER"], header["YLLCORNER"], header["CELLSIZE"]
    nrows, ncols = int(header["NROWS"]), int(header["NCOLS"])
    c0 = max(0, int((e_min - xll) / cs))
    c1 = min(ncols, int((e_max - xll) / cs) + 1)
    r0 = max(0, int(nrows - 1 - (n_max - yll) / cs))
    r1 = min(nrows, int(nrows - 1 - (n_min - yll) / cs) + 1)
    return r0, r1, c0, c1


def cell_centre_en(header, row_full, col_full):
    """Return (easting, northing) of a grid cell's centre given full-grid row/col."""
    xll, yll, cs = header["XLLCORNER"], header["YLLCORNER"], header["CELLSIZE"]
    nrows = int(header["NROWS"])
    e = xll + col_full * cs + cs / 2
    n = yll + (nrows - 1 - row_full) * cs + cs / 2
    return e, n


def is_genuine_ors(max_period, min_period, rcs_d):
    """True if a bedrock polygon is Devonian-age AND sedimentary (not igneous/metamorphic)."""
    if pd.isna(max_period):
        max_period = ""
    if pd.isna(min_period):
        min_period = ""
    if pd.isna(rcs_d):
        rcs_d = ""
    is_devonian = (max_period == "DEVONIAN") or (min_period == "DEVONIAN")
    is_non_sed = any(term in rcs_d for term in NON_SEDIMENTARY_TERMS)
    return bool(is_devonian and not is_non_sed)


def main():
    print("Loading Ba and Rb grids...")
    ba, ba_header = load_grid(BA_PATH)
    rb, rb_header = load_grid(RB_PATH)
    assert ba.shape == rb.shape, "Ba and Rb grids must have identical dimensions"
    nodata = ba_header["NODATA_VALUE"]

    print("Computing Ba/Rb ratio...")
    ratio = np.full(ba.shape, np.nan, dtype=np.float32)
    valid_all = (ba != nodata) & (rb != nodata) & (rb > 0)
    ratio[valid_all] = ba[valid_all] / rb[valid_all]

    print(f"Clipping to study extent {STUDY_BOX} ...")
    r0, r1, c0, c1 = box_to_indices(ba_header, **STUDY_BOX)
    ba_sub = ba[r0:r1, c0:c1].astype(float)
    rb_sub = rb[r0:r1, c0:c1].astype(float)
    ratio_sub = ratio[r0:r1, c0:c1]
    valid_sub = (ba_sub != nodata) & (rb_sub != nodata) & (~np.isnan(ratio_sub))

    ratio_threshold = float(np.percentile(ratio_sub[valid_sub], RATIO_PERCENTILE))
    print(f"P{RATIO_PERCENTILE} Ba/Rb ratio within study extent: {ratio_threshold:.3f}")

    composite = (ba_sub >= BA_FLOOR_PPM) & (ratio_sub >= ratio_threshold) & valid_sub
    print(f"Cells passing composite screen: {int(composite.sum())}")

    print("Clustering (8-connected)...")
    labeled, n_raw_clusters = ndimage.label(composite, structure=np.ones((3, 3)))

    print("Loading BGS bedrock polygons and running per-cell spatial join "
          "(this is the slow step)...")
    bedrock = gpd.read_file(BEDROCK_SHP)

    rows_idx, cols_idx = np.where(composite)
    points, cluster_ids, cell_e, cell_n = [], [], [], []
    for rr, cc in zip(rows_idx, cols_idx):
        e, n = cell_centre_en(ba_header, r0 + rr, c0 + cc)
        points.append(Point(e, n))
        cluster_ids.append(int(labeled[rr, cc]))
        cell_e.append(e)
        cell_n.append(n)

    pts_gdf = gpd.GeoDataFrame({"cluster_id": cluster_ids}, geometry=points, crs=bedrock.crs)
    joined = gpd.sjoin(
        pts_gdf,
        bedrock[["LEX_D", "RCS_D", "MAX_PERIOD", "MIN_PERIOD", "geometry"]],
        how="left", predicate="within",
    )
    joined = joined[~joined.index.duplicated(keep="first")]  # guard against edge double-hits

    is_ors = [
        is_genuine_ors(row.MAX_PERIOD, row.MIN_PERIOD, row.RCS_D)
        for row in joined.itertuples()
    ]

    inv = Transformer.from_crs(bedrock.crs, "EPSG:4326", always_xy=True)
    lons, lats = inv.transform(cell_e, cell_n)

    per_cell = pd.DataFrame({
        "cluster_id": cluster_ids,
        "easting": cell_e,
        "northing": cell_n,
        "lon": lons,
        "lat": lats,
        "ba_ppm": [ba_sub[rr, cc] for rr, cc in zip(rows_idx, cols_idx)],
        "rb_ppm": [rb_sub[rr, cc] for rr, cc in zip(rows_idx, cols_idx)],
        "ba_rb_ratio": [ratio_sub[rr, cc] for rr, cc in zip(rows_idx, cols_idx)],
        "formation": joined["LEX_D"].values,
        "lithology": joined["RCS_D"].values,
        "max_period": joined["MAX_PERIOD"].values,
        "min_period": joined["MIN_PERIOD"].values,
        "is_genuine_ors": is_ors,
    })

    print("Aggregating per-cluster statistics...")
    cluster_rows = []
    for cid in sorted(per_cell["cluster_id"].unique()):
        sub = per_cell[per_cell["cluster_id"] == cid]
        n_cells = len(sub)
        area_km2 = n_cells * (CELL_SIZE_M / 1000.0) ** 2
        if area_km2 < MIN_CLUSTER_CELLS * (CELL_SIZE_M / 1000.0) ** 2:
            continue
        pct_ors = 100.0 * sub["is_genuine_ors"].mean()
        ors_formations = sub.loc[sub["is_genuine_ors"], "formation"]
        dominant = ors_formations.value_counts().idxmax() if len(ors_formations) else (
            sub["formation"].value_counts().idxmax() if sub["formation"].notna().any() else "NO POLYGON"
        )
        cluster_rows.append({
            "cluster_id": cid,
            "n_cells": n_cells,
            "area_km2": area_km2,
            "centroid_lat": sub["lat"].mean(),
            "centroid_lon": sub["lon"].mean(),
            "mean_ba_ppm": sub["ba_ppm"].mean(),
            "max_ba_ppm": sub["ba_ppm"].max(),
            "mean_rb_ppm": sub["rb_ppm"].mean(),
            "mean_ratio": sub["ba_rb_ratio"].mean(),
            "pct_genuine_ors": pct_ors,
            "dominant_formation": dominant,
        })
    per_cluster = pd.DataFrame(cluster_rows).sort_values("area_km2", ascending=False)

    import os
    os.makedirs("outputs", exist_ok=True)
    per_cell.to_csv("outputs/per_cell_results.csv", index=False)
    per_cluster.to_csv("outputs/per_cluster_results.csv", index=False)

    print(f"\nDone. {len(per_cluster)} clusters >= {MIN_CLUSTER_CELLS} cells written to "
          f"outputs/per_cluster_results.csv")
    print(f"{len(per_cell)} individual cells written to outputs/per_cell_results.csv")
    print("\nTop 5 clusters by area:")
    print(per_cluster.head(5).to_string(index=False))


if __name__ == "__main__":
    main()
