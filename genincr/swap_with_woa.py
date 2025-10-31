#!/usr/bin/env python3
import xarray as xr
import numpy as np
import sys
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description='Replace Temp/Salt data in GDAS increment file with WOA'
    )
    parser.add_argument('src_woa_incr_file', help='Source WOA increments file')
    parser.add_argument('src_jedi_incr_file',
                        help='Source GDAS increments file')
    parser.add_argument('out_file', help='Output file (non-destructive)')
    return parser.parse_args()


# Parse command line arguments
args = parse_args()
src_woa_incr_file = args.src_woa_incr_file
src_jedi_incr_file = args.src_jedi_incr_file
out_file = args.out_file

print(f"📥 Reading source climatology increments: {src_woa_incr_file}")
woa = xr.open_dataset(src_woa_incr_file, decode_times=False)

print(f"📥 Reading target GDAS increment file: {src_jedi_incr_file}")
gdas = xr.open_dataset(src_jedi_incr_file, decode_times=False)

# --- Step 1: Sanity checks ---
for var in ["Temp", "Salt"]:
    if var not in woa:
        sys.exit(f"❌ Variable {var} not found in {src_woa_incr_file}")
    if var not in gdas:
        sys.exit(f"❌ Variable {var} not found in {src_jedi_incr_file}")
    if woa[var].shape != gdas[var].shape:
        print(f"⚠️ Shape mismatch for {var}: "
              f"source {woa[var].shape} vs dest {gdas[var].shape}")
    else:
        print(f"✅ {var} shape OK: {woa[var].shape}")

# --- Step 2: Replace NaNs with zeros before assignment ---


def clean_nan(da):
    arr = da.values
    nan_count = np.isnan(arr).sum()
    if nan_count > 0:
        print(f"🧹 Replacing {nan_count:,} NaNs with 0 in {da.name}")
        arr = np.nan_to_num(arr, nan=0.0)
    return arr


print("🧬 Replacing Temp and Salt in GDAS file...")
gdas["Temp"].data = clean_nan(woa["Temp"])
gdas["Salt"].data = clean_nan(woa["Salt"])

# --- Step 3: Update fill values to 0 ---
for var in ["Temp", "Salt"]:
    gdas[var].attrs["_FillValue"] = 0.0

# --- Step 4: Preserve attributes and record provenance ---
if "history" not in gdas.attrs:
    gdas.attrs["history"] = ""
gdas.attrs["history"] += (
    f"\nReplaced Temp/Salt from {src_woa_incr_file} using Python/xarray; "
    f"NaNs set to 0."
)

# --- Step 5: Write updated file ----------------------------------------------
print(f"💾 Writing updated file → {out_file}")
gdas.to_netcdf(out_file, format="NETCDF4_CLASSIC")

print("✅ Done. The GDAS increment file has been updated cleanly.")
