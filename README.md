script.cmd

#PBS -o logfile.log
#PBS -e errorfile_slash.err
#PBS -l walltime=00:60:00
#PBS -l select=1:ncpus=32:mpiprocs=32
#PBS -q rupesh_gpuq

export PMIX_MCA_gds=hash

module load openmpi411

mpicxx -I/lfs/usrhome/oth/rnintern/scratch/MPI_Comparison/boost/install_dir/include /lfs/usrhome/oth/rnintern/scratch/poulomi/pford.cpp -o pffout

mv pffout $PBS_O_WORKDIR

echo "Done Compile #2"

/lfs/sware/openmpi411/bin/mpirun -np 2 $PBS_O_WORKDIR/pffout  /lfs/usrhome/oth/rnintern/scratch/poulomi/input.txt > $PBS_O_WORKDIR/poulomitwo.txt
