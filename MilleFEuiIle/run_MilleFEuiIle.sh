# --- Number of cores available ---
ncores=8

# --- Number of cores for MPI ---
n_cores=1

# --- Loop over the presumed age of the bands (Myr) ---
for clathrates in 8 9 10 #0 1 2 3 4 5 6 7 # 10 20
    do
    # --- Loop over the ice shell thickness (km) ---
    for viscosity in 15 #14 # 13
        do
        while true
            do 
                sleep 1        
                if [ $(ps -ef | grep -v grep | grep main_run | wc -l) -lt $ncores ]; then
                    break
                fi
            done
            # --- Copy the main file ---
            main_orig="main.py"
            main_run="main_run_D"$clathrates"km_eta"$viscosity"Pas.py" 

            cp $main_orig $main_run

            param_file1="m_parameters.py"
            param_file2="m_parameters_Titan.py"

            cp $param_file2 $param_file1

            # --- Define the log files ---
            out_file1=$clathrates"km_eta"$viscosity"Pas.out"
            out_file2=$clathrates"km_eta"$viscosity"Pas_e.out"

            # --- Run the code ---
            if [ $n_cores -eq 1 ]; then
                python $main_run $clathrates $viscosity > $out_file1 2> $out_file2&
            else
                mpirun -n $n_cores python $main_run $clathrates $viscosity > $out_file1 2> $out_file2&
            fi

        done
    done