export PROCS_PER_SOCKET=1
export OMP_NUM_THREADS=6

export OMP_DISPLAY_ENV=true
# this is required
export OMP_NESTED=false
export OMP_MAX_ACTIVE_LEVELS=1
# this can be tuned
export OMP_PROC_BIND=spread

# this uses all hyper-threading threads
# export OMP_PLACES=threads
# mpirun --map-by ppr:$PROCS_PER_SOCKET:socket:PE=$OMP_NUM_THREADS \
#        --bind-to hwthread --display-map ./aLENS.X > ./outrun.log

# this does not use hyper-threading
export OMP_PLACES=cores
mpirun --map-by ppr:$PROCS_PER_SOCKET:socket:PE=$OMP_NUM_THREADS \
       --bind-to core --display-map ./aLENS.X >./outrun.log
