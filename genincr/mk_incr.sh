#!/bin/bash

module use ../../sandboxes/global-workflow/sorc/gdas.cd/modulefiles/
module load GDAS/ursa.intel
source /scratch3/NCEPDEV/da/Guillaume.Vernieres/venvs/interpmom6/bin/activate
export PATH="/scratch3/NCEPDEV/da/Guillaume.Vernieres/data/woa/FRE-NCtools/install/bin:$PATH"

yyyymmdd=20251008
hh=12
bkg_file=./gdas.ocean.t${hh}z.inst.f006.nc
layer_file=/scratch3/NCEPDEV/da/Guillaume.Vernieres/runs/gfsv17/mlb1/COMROOT/mlb1/gdas.20250831/00/analysis/ocean/gdas.t00z.ocn.ana.nc
template_incrfile=./gdas.t${hh}z.ocninc.nc

# Extract the climatology for yyyymmdd hh
# ---------------------------------------
./splice_clim.py --date ${yyyymmdd}${hh} \
                 --out ${yyyymmdd}${hh}_woa23_ts.nc \
                 --layer-file $layer_file \
                 --monthly-dir /scratch3/NCEPDEV/da/Guillaume.Vernieres/data/woa \
                 --depth-threshold 1400

# Compute the increment
# ---------------------
# Make a copy of the background file and rename dimensions/variables in it
cp $bkg_file ${yyyymmdd}_bkg.nc
ncrename -d xh,xaxis_1 -d yh,yaxis_1 -d z_l,zaxis_1 -v xh,xaxis_1 -v yh,yaxis_1 -v z_l,zaxis_1 ${yyyymmdd}_bkg.nc
ncrename -d time,Time -v time,Time ${yyyymmdd}_bkg.nc

ncdiff -O -v Temp,Salt ${yyyymmdd}${hh}_woa23_ts.nc ${yyyymmdd}_bkg.nc ${yyyymmdd}${hh}_woa23_incr_ts.nc

# Set fill values to zero
#ncatted -O -a _FillValue,Temp,o,f,0.0 -a _FillValue,Salt,o,f,0.0 ${yyyymmdd}${hh}_incr.nc

# Prepare mom6 iau file
# ---------------------
# make a copy of the template incr file
cp ${template_incrfile} gdas.t${hh}z.ocninc.jedi.nc


# replace Temp and Salt with the computed increment
python swap_with_woa.py ${yyyymmdd}${hh}_woa23_incr_ts.nc gdas.t${hh}z.ocninc.jedi.nc gdas.t${hh}z.ocninc.woa.nc
