#!/bin/bash
#SBATCH --account=fermi:users
#SBATCH --partition=milano
#SBATCH --job-name=LS5039_srcProb
#SBATCH --output=output-%j.txt
#SBATCH --error=output-%j.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=24g
#SBATCH --time=0-04:00:00
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=richard@slac.stanford.edu

set -e

hostname

# Accessing the start and end times passed as arguments
start_time="$1"
end_time="$2"

cd $LSCRATCH
pwd

infile="/sdf/home/r/richard/fermi-user/LS5039/periods/periodicity/LS5039_3deg_gated.fits"
scfile="/sdf/home/r/richard/fermi-user/LS5039/fssc_data/L24090123505104476C3F11_SC00.fits"
srcmdl="/sdf/home/r/richard/fermi-user/LS5039/periods/J1826-1256/fullspan/pickle2_srcmdl_00.xml"
src_list="/sdf/home/r/richard/fermi-user/LS5039/periods/J1826-1256/fullspan/godot/ls5039_list.txt"
outfile="LS5039_3_srcprob.fits"
final_dir="/sdf/home/r/richard/fermi-user/LS5039/periods/J1826-1256/fullspan/godot/diffrsp/"

# Your commands to be executed
echo "Job running with start time: $start_time seconds and end time: $end_time seconds"
echo "Writing to $outfile"

# trim file by time; run gtdiffrsp; run gtsrcprob

echo "running gtselect"
time -p gtselect tmin=$start_time tmax=$end_time emin=100. emax=300000 zmax=90. ra=INDEF dec=INDEF rad=INDEF outfile=LS5039_gts.fits infile=$infile

echo "running gtdiffrsp"
time -p gtdiffrsp evfile=LS5039_gts.fits scfile=$scfile srcmdl=$srcmdl irfs=CALDB

echo "running gtsrcprob"
time -p gtsrcprob evfile=LS5039_gts.fits scfile=$scfile outfile=$outfile srcmdl=$srcmdl srclist=$src_list irfs=CALDB

cp $outfile $final_dir

# clean up $LSCRATCH
echo "ls pwd"
pwd
ls -lh $LSCRATCH
#rm $LSCRATCH/*.par
rm $LSCRATCH/*.fits
