#!/bin/bash

# WOA to MOM6 Grid Interpolation Script
# This script interpolates WOA18 data to a MOM6 grid using FRE tools
#
# Note: The mosaic grid created by make_hgrid has twice the number of grid points
# along both x and y axes compared to the WOA grid. This is a GFDL convention
# where the grid points represent cell corners/edges rather than cell centers.

# Default values
DEFAULT_INPUT_FILE="woa18_decav_t00_01.nc"
DEFAULT_OUTPUT_FILE="WOA18_on_MOM6.nc"
DEFAULT_NLON=720
DEFAULT_NLAT=360
DEFAULT_XBNDS="-180,180"
DEFAULT_YBNDS="-90,90"
DEFAULT_GRID_NAME="woa100_grid"
DEFAULT_MOSAIC_NAME="woa100_mosaic"
DEFAULT_OUTPUT_MOSAIC="./dst/ocean_mosaic.nc"
DEFAULT_SCALAR_FIELD="t_an"
DEFAULT_INTERP_METHOD="conserve_order1"
DEFAULT_LAYER_FILE="./dst/MOM6_layer_h.nc"
DEFAULT_VERT_OUTPUT="WOA18_on_MOM6_layers.nc"
DEFAULT_DEPTH_VAR="depth"
DEFAULT_LAYER_VAR="h"
# WOA23 download defaults (temperature only)
DEFAULT_WOA23_RES="0.25"       # 0.25 | 1.00 | 5.00
DEFAULT_WOA23_PERIOD="annual"  # annual | monthly
DEFAULT_WOA23_MONTH=""         # 01..12 (required for monthly; defaults to 01)
DEFAULT_WOA23_OUTDIR="./woa23"

# Function to show usage
show_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -i, --input-file FILE       Input WOA file (default: $DEFAULT_INPUT_FILE)"
    echo "  -o, --output-file FILE      Output file (default: $DEFAULT_OUTPUT_FILE)"
    echo "  --nlon N                    Number of longitude points (default: $DEFAULT_NLON)"
    echo "  --nlat N                    Number of latitude points (default: $DEFAULT_NLAT)"
    echo "  --xbnds LON1,LON2          Longitude bounds (default: $DEFAULT_XBNDS)"
    echo "  --ybnds LAT1,LAT2          Latitude bounds (default: $DEFAULT_YBNDS)"
    echo "  --grid-name NAME           Grid name (default: $DEFAULT_GRID_NAME)"
    echo "  --mosaic-name NAME         Mosaic name (default: $DEFAULT_MOSAIC_NAME)"
    echo "  --output-mosaic FILE       Output mosaic file (default: $DEFAULT_OUTPUT_MOSAIC)"
    echo "  --scalar-field FIELD       Variable name(s) for interpolation."
    echo "                             Horizontal uses the first; vertical supports multiple"
    echo "                             when comma/space separated (e.g., t_an,s_an)."
    echo "  --interp-method METHOD     Interpolation method (default: $DEFAULT_INTERP_METHOD)"
    echo "  --vertical                 Perform vertical interpolation after horizontal"
    echo "  --layer-file FILE          MOM6 layer thickness file (default: $DEFAULT_LAYER_FILE)"
    echo "  --vert-output FILE         Vertical interpolation output file (default: $DEFAULT_VERT_OUTPUT)"
    echo "  --depth-var VAR            Depth variable name (default: $DEFAULT_DEPTH_VAR)"
    echo "  --layer-var VAR            Layer thickness variable name (default: $DEFAULT_LAYER_VAR)"
    echo "  -j, --jobs N               Number of parallel processes for vertical interp (default: auto)"
    echo "  --no-extrapolate           Don't use extrapolation in fregrid"
    echo ""
    echo "WOA23 temperature download (regridding input convenience):"
    echo "  --download-woa23           Download WOA23 temperature (decav) file(s)"
    echo "  --woa23-resolution R       0.25 | 1.00 | 5.00 (default: $DEFAULT_WOA23_RES)"
    echo "  --woa23-period P           annual | monthly (default: $DEFAULT_WOA23_PERIOD)"
    echo "  --woa23-month MM           01..12 (used when period=monthly; default: 01)"
    echo "  --woa23-outdir DIR         Download directory (default: $DEFAULT_WOA23_OUTDIR)"
    echo "  -h, --help                 Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0                                    # Use all defaults (horizontal only)"
    echo "  $0 --download-woa23 --woa23-resolution 0.25 --woa23-period annual"
    echo "  $0 --download-woa23 --woa23-resolution 1.00 --woa23-period monthly --woa23-month 07"
    echo "  $0 -i my_woa.nc -o my_output.nc     # Custom input/output files"
    echo "  $0 --nlon 720 --nlat 360            # 1-degree grid (360×180 actual cells)"
    echo "  $0 --vertical                        # Perform both horizontal and vertical interpolation"
    echo "  $0 --vertical --layer-file my_layers.nc  # Vertical interp with custom layer file"
    echo "  $0 --scalar-field t_an,s_an --vertical   # Vertical interpolate multiple variables"
    echo "  $0 --vertical -j 8                   # Use 8 parallel processes for vertical interp"
}

# Parse command line arguments
INPUT_FILE="$DEFAULT_INPUT_FILE"
OUTPUT_FILE="$DEFAULT_OUTPUT_FILE"
NLON="$DEFAULT_NLON"
NLAT="$DEFAULT_NLAT"
XBNDS="$DEFAULT_XBNDS"
YBNDS="$DEFAULT_YBNDS"
GRID_NAME="$DEFAULT_GRID_NAME"
MOSAIC_NAME="$DEFAULT_MOSAIC_NAME"
OUTPUT_MOSAIC="$DEFAULT_OUTPUT_MOSAIC"
SCALAR_FIELD="$DEFAULT_SCALAR_FIELD"
INTERP_METHOD="$DEFAULT_INTERP_METHOD"
DO_VERTICAL=false
LAYER_FILE="$DEFAULT_LAYER_FILE"
VERT_OUTPUT="$DEFAULT_VERT_OUTPUT"
DEPTH_VAR="$DEFAULT_DEPTH_VAR"
LAYER_VAR="$DEFAULT_LAYER_VAR"
NUM_JOBS=""
USE_EXTRAPOLATE=true
# WOA23 download flags
DOWNLOAD_WOA23=false
WOA23_RES="$DEFAULT_WOA23_RES"
WOA23_PERIOD="$DEFAULT_WOA23_PERIOD"
WOA23_MONTH="$DEFAULT_WOA23_MONTH"
WOA23_OUTDIR="$DEFAULT_WOA23_OUTDIR"

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input-file)
            INPUT_FILE="$2"; shift 2 ;;
        -o|--output-file)
            OUTPUT_FILE="$2"; shift 2 ;;
        --nlon)
            NLON="$2"; shift 2 ;;
        --nlat)
            NLAT="$2"; shift 2 ;;
        --xbnds)
            XBNDS="$2"; shift 2 ;;
        --ybnds)
            YBNDS="$2"; shift 2 ;;
        --grid-name)
            GRID_NAME="$2"; shift 2 ;;
        --mosaic-name)
            MOSAIC_NAME="$2"; shift 2 ;;
        --output-mosaic)
            OUTPUT_MOSAIC="$2"; shift 2 ;;
        --scalar-field)
            SCALAR_FIELD="$2"; shift 2 ;;
        --interp-method)
            INTERP_METHOD="$2"; shift 2 ;;
        --vertical)
            DO_VERTICAL=true; shift ;;
        --layer-file)
            LAYER_FILE="$2"; shift 2 ;;
        --vert-output)
            VERT_OUTPUT="$2"; shift 2 ;;
        --depth-var)
            DEPTH_VAR="$2"; shift 2 ;;
        --layer-var)
            LAYER_VAR="$2"; shift 2 ;;
        -j|--jobs)
            NUM_JOBS="$2"; shift 2 ;;
        --no-extrapolate)
            USE_EXTRAPOLATE=false; shift ;;
        # WOA23 download options
        --download-woa23)
            DOWNLOAD_WOA23=true; shift ;;
        --woa23-resolution)
            WOA23_RES="$2"; shift 2 ;;
        --woa23-period)
            WOA23_PERIOD="$2"; shift 2 ;;
        --woa23-month)
            WOA23_MONTH="$2"; shift 2 ;;
        --woa23-outdir)
            WOA23_OUTDIR="$2"; shift 2 ;;
        -h|--help)
            show_usage; exit 0 ;;
        *)
            echo "Unknown option: $1"; show_usage; exit 1 ;;
    esac
done

# If requested, download WOA23 temperature input and set INPUT_FILE accordingly
if [[ "$DOWNLOAD_WOA23" == "true" ]]; then
    BASE_URL="https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/temperature/netcdf/decav"
    case "$WOA23_RES" in
        0.25|0.250)
            RES_DIR="0.25"; RES_CODE="04" ;;
        1.00|1|1.0)
            RES_DIR="1.00"; RES_CODE="01" ;;
        5.00|5|5.0)
            RES_DIR="5.00"; RES_CODE="5d" ;;
        *)
            echo "❌ Invalid --woa23-resolution: $WOA23_RES (use 0.25 | 1.00 | 5.00)"; exit 1 ;;
    esac

    if [[ "$WOA23_PERIOD" == "annual" ]]; then
        FNAME="woa23_decav_t00_${RES_CODE}.nc"
    elif [[ "$WOA23_PERIOD" == "monthly" ]]; then
        MM="$WOA23_MONTH"
        if [[ -z "$MM" ]]; then
            echo "⚠️  --woa23-month not provided for monthly; defaulting to 01"
            MM="01"
        fi
        if ! [[ "$MM" =~ ^(0[1-9]|1[0-2])$ ]]; then
            echo "❌ Invalid --woa23-month: $MM (use 01..12)"; exit 1
        fi
        FNAME="woa23_decav_t${MM}_${RES_CODE}.nc"
    else
        echo "❌ Invalid --woa23-period: $WOA23_PERIOD (use annual | monthly)"; exit 1
    fi

    URL="${BASE_URL}/${RES_DIR}/${FNAME}"
    mkdir -p "$WOA23_OUTDIR"
    TARGET="${WOA23_OUTDIR}/${FNAME}"

    if [[ -f "$TARGET" ]]; then
        echo "ℹ️  Using existing file: $TARGET"
    else
        echo "⬇️  Downloading $URL"
        if ! curl -fSL -o "$TARGET" "$URL"; then
            echo "❌ Download failed: $URL"; exit 1
        fi
    fi

    # Use the downloaded file as interpolation input
    INPUT_FILE="$TARGET"
    echo "✅ WOA23 file ready: $INPUT_FILE"
fi

# Validate input file exists
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: Input file '$INPUT_FILE' not found!"
    exit 1
fi

# Validate vertical interpolation requirements
if [[ "$DO_VERTICAL" == "true" ]]; then
    if [[ ! -f "$LAYER_FILE" ]]; then
        echo "Error: Layer file '$LAYER_FILE' not found! Required for vertical interpolation."
        exit 1
    fi
    if [[ ! -f "vert_interp_3d.py" ]]; then
        echo "Error: vert_interp_3d.py not found in current directory!"
        exit 1
    fi
fi

# Show configuration
echo "🌊 WOA to MOM6 Grid Interpolation"
echo "================================="
SCRIPT_START_TIME=$(date +%s)
echo "Input file:       $INPUT_FILE"
echo "Output file:      $OUTPUT_FILE"
echo "Grid resolution:  ${NLON}×${NLAT}"
echo "Longitude bounds: $XBNDS"
echo "Latitude bounds:  $YBNDS"
echo "Grid name:        $GRID_NAME"
echo "Mosaic name:      $MOSAIC_NAME"
echo "Scalar field:     $SCALAR_FIELD"
echo "Interp method:    $INTERP_METHOD"
if [[ "$DO_VERTICAL" == "true" ]]; then
    echo "Vertical interp:  YES"
    echo "Layer file:       $LAYER_FILE"
    echo "Vert output:      $VERT_OUTPUT"
    echo "Depth variable:   $DEPTH_VAR"
    echo "Layer variable:   $LAYER_VAR"
else
    echo "Vertical interp:  NO"
fi
echo ""

echo "🗂️  Creating grid..."
GRID_START_TIME=$(date +%s)
make_hgrid \
  --grid_type regular_lonlat_grid \
  --nlon "$NLON" \
  --nlat "$NLAT" \
  --nxbnds 2 \
  --nybnds 2 \
  --xbnds "$XBNDS" \
  --ybnds "$YBNDS" \
  --grid_name "$GRID_NAME" \
  --center t_cell
GRID_END_TIME=$(date +%s)
GRID_DURATION=$((GRID_END_TIME - GRID_START_TIME))
echo "⏱️  Grid creation took ${GRID_DURATION}s"

echo "🗂️  Creating mosaic..."
MOSAIC_START_TIME=$(date +%s)
make_solo_mosaic \
  --num_tiles 1 \
  --dir . \
  --tile_file "${GRID_NAME}.nc" \
  --mosaic_name "$MOSAIC_NAME" \
  --periodx 1
MOSAIC_END_TIME=$(date +%s)
MOSAIC_DURATION=$((MOSAIC_END_TIME - MOSAIC_START_TIME))
echo "⏱️  Mosaic creation took ${MOSAIC_DURATION}s"

# Build fregrid command
echo "🔄 Interpolating data..."
FREGRID_START_TIME=$(date +%s)
FREGRID_CMD="fregrid \
  --input_file=$INPUT_FILE \
  --input_mosaic=${MOSAIC_NAME}.nc \
  --remap_file=woa_to_mom6_weights.nc \
  --output_mosaic=$OUTPUT_MOSAIC \
  --output_file=$OUTPUT_FILE \
  --interp_method=$INTERP_METHOD \
  --scalar_field=$SCALAR_FIELD"

if [[ "$USE_EXTRAPOLATE" == "true" ]]; then
    FREGRID_CMD="$FREGRID_CMD --extrapolate"
fi

# Execute fregrid
eval "$FREGRID_CMD"
FREGRID_END_TIME=$(date +%s)
FREGRID_DURATION=$((FREGRID_END_TIME - FREGRID_START_TIME))
echo "⏱️  Horizontal interpolation took ${FREGRID_DURATION}s"
echo "✅ Horizontal interpolation complete! Output file: $OUTPUT_FILE"

# Perform vertical interpolation if requested
if [[ "$DO_VERTICAL" == "true" ]]; then
    echo ""
    echo "🔄 Performing vertical interpolation..."

    # Check if the horizontal output exists
    if [[ ! -f "$OUTPUT_FILE" ]]; then
        echo "❌ Error: Horizontal interpolation output '$OUTPUT_FILE' not found!"
        exit 1
    fi

    # Build vertical interpolation command
    VERT_CMD="./vert_interp_3d.py"
    VERT_CMD="$VERT_CMD --woa-file '$OUTPUT_FILE'"
    VERT_CMD="$VERT_CMD --layer-file '$LAYER_FILE'"
    VERT_CMD="$VERT_CMD --output-file '$VERT_OUTPUT'"
    VERT_CMD="$VERT_CMD --depth-var '$DEPTH_VAR'"
    VERT_CMD="$VERT_CMD --layer-var '$LAYER_VAR'"

    # Support multiple variables (comma or space separated) from SCALAR_FIELD
    VAR_LIST="${SCALAR_FIELD//,/ }"
    for v in $VAR_LIST; do
        VERT_CMD="$VERT_CMD --variable '$v'"
    done

    if [[ -n "$NUM_JOBS" ]]; then
        VERT_CMD="$VERT_CMD --jobs $NUM_JOBS"
    fi

    # Execute vertical interpolation
    VERT_START_TIME=$(date +%s)
    eval "$VERT_CMD"
    VERT_STATUS=$?
    VERT_END_TIME=$(date +%s)
    VERT_DURATION=$((VERT_END_TIME - VERT_START_TIME))

    if [[ $VERT_STATUS -eq 0 ]]; then
        echo "⏱️  Vertical interpolation took ${VERT_DURATION}s"
        echo "✅ Vertical interpolation complete! Final output: $VERT_OUTPUT"
        echo "📋 Summary:"
        echo "   Horizontal:  $OUTPUT_FILE"
        echo "   Vertical:    $VERT_OUTPUT"
    else
        echo "❌ Vertical interpolation failed!"
        exit $VERT_STATUS
    fi
else
    echo "📋 Summary:"
    echo "   Horizontal only:  $OUTPUT_FILE"
fi

SCRIPT_END_TIME=$(date +%s)
TOTAL_DURATION=$((SCRIPT_END_TIME - SCRIPT_START_TIME))
echo ""
echo "⏱️  === TIMING SUMMARY ==="
echo "Grid creation:      ${GRID_DURATION}s"
echo "Mosaic creation:    ${MOSAIC_DURATION}s"
echo "Horizontal interp:  ${FREGRID_DURATION}s"
if [[ "$DO_VERTICAL" == "true" ]]; then
    echo "Vertical interp:    ${VERT_DURATION}s"
fi
echo "Total runtime:      ${TOTAL_DURATION}s"
