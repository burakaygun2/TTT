# --- Number of cores available ---
ncores=16

# --- Number of cores for MPI ---
n_cores=1

# --- Ice shell thickness ---
for thickness in 40 50 60 70 80 90 100 110 120 130 140 150 160 167
    do
    # --- Loop over the ice shell thickness (km) ---
    for viscosity in 15 14 13
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
            main_run="main_run_D"$thickness"km_eta"$viscosity"Pas.py" 

            cp $main_orig $main_run

            param_file1="m_parameters.py"
            param_file2="m_parameters_Titan.py"

            cp $param_file2 $param_file1

            # --- Define the log files ---
            out_file1=$thickness"km_eta"$viscosity"Pas.out"
            out_file2=$thickness"km_eta"$viscosity"Pas_e.out"

            # --- Run the code ---
            if [ $n_cores -eq 1 ]; then
                python $main_run $thickness $viscosity > $out_file1 2> $out_file2&
            else
                mpirun -n $n_cores python $main_run $thickness $viscosity > $out_file1 2> $out_file2&
            fi

        done
    done