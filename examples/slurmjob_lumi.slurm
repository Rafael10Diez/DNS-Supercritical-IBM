#!/bin/bash -l
#SBATCH -J dns
#SBATCH -o %x-j%j.out
#SBATCH --gpus-per-node=8
#SBATCH --ntasks-per-node=8 
#SBATCH -N 3
#SBATCH --account=project_???
#SBATCH --partition=standard-g
#SBATCH --time 2-00:00:00
#SBATCH --exclusive

# https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/lumig-job/
echo "Starting job $SLURM_JOB_ID at `date`"
echo "USER: ${USER}"
echo "PWD: ${PWD}"

export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_OFI_NIC_POLICY=GPU
export ACC_DEVICE_NUM=8

export EXTRA_AMD_FFLAGS=" --offload-arch=${ROCM_GPU} -lhipfft -lrocfft -lhipfort-amdgcn -hnoomp -hacc -hacc_model=auto_async_none:no_fast_addr -I$HIPFORT_HOME/include/hipfort/amdgcn -L$HIPFORT_HOME/lib -lhipfort-amdgcn "

# --------------------------------------------------------------------------------

root_folder_dir=/scratch/project_???/$USER/[work_folder]
repo_name=DNS-Supercritical-IBM
runtime_dir=$root_folder_dir/$repo_name/runtime

# ------------ Please clone this repository ------------
cd $root_folder_dir
git clone https://github.com/Rafael10Diez/DNS-Supercritical-IBM

# copy input folder
cd $root_folder_dir/$repo_name
cp -R ../orig_input input
cp -R src runtime

# copy expressions_rho_mu_cond_buoyancy_cp_div_jacobian.f90
cd $runtime_dir
git clone https://github.com/Rafael10Diez/diezDecomp
mkdir -p thermophysical_properties
cp -f ../input/expressions_rho_mu_cond_buoyancy_cp_div_jacobian.f90 ./thermophysical_properties
mkdir -p metaprogramming
cp -f ../input/prof_1d_variables.txt metaprogramming
cp -f ../input/avg_variables.txt     metaprogramming

python3 ./pycode/quick_includer.py # pre-compiled by Python script
python3 ./pycode/calc_bands.py     # pre-compiled by Python script

cd $runtime_dir
my_objs=""
for my_file in diezDecomp/src/diezdecomp_core \
               diezDecomp/src/diezdecomp_api_generic \
               diezdecomp_api_ibm \
               mod_tictoc \
               mod_utils \
               mod_io \
               poisson_solver_mgpu_dtdma/poisson_solver_multigpu \
               mod_avg_apply \
               mod_get_bulk \
               ibm_setup/mod_cfd_ibm \
               mod_cfd_subroutines \
               main; do
  if [[ "$my_file" == "main" ]] ; then 
    my_flag=""
    my_other=$my_objs
  else
    my_flag="-c"
    my_other=""
  fi 
  # -D_USE_TICTOC_TIMER_ONLY_0_RANK -D_USE_TICTOC_TIMER -D_USE_TICTOC_VERBOSE_TOC
  echo "file: $my_file ; flag: $my_flag ; other: $my_other"
  hipfc  -lmpi -lhipfft -O3 ${EXTRA_AMD_FFLAGS}  "-DHIPFORT_ARCH=\"amd\"" "-s real64" -D_USE_HIPFFT -D_BIN_REDUCTION  ${my_flag} -o ${my_file}.o ${my_other} ${my_file}.f90
  my_objs="$my_objs $my_file.o"
done
chmod u+x main.o
mkdir -p ini

cd $runtime_dir

cat << EOF > select_gpu
#!/bin/bash
export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID
exec \$*
EOF
chmod +x ./select_gpu
export CPU_BIND="map_cpu:49,57,17,25,1,9,33,41"
export my_nproc=$(( 8 * $SLURM_JOB_NUM_NODES ))
export my_divs_z=1
export my_divs_y=$(( 8 * $SLURM_JOB_NUM_NODES ))
echo "my_nproc: $my_nproc ; my_divs_y: $my_divs_y ; my_divs_z: $my_divs_z"
srun -n ${my_nproc} --cpu-bind=${CPU_BIND} --distribution=block:block ./select_gpu ./main.o ${my_divs_z} ${my_divs_y} 8


