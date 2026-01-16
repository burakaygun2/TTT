module mod_write_h5
    use mod_parameters
    use mod_variables
    use mod_potential
    use hdf5
    implicit none
    contains

    subroutine write_summary_lsingle_h5
        implicit none
        character(300)   :: filename
        real(dp)         :: data_to_write(1)
        integer          :: data_to_write_int(1)
        !! Variables regarding the hdf5
        integer          :: error
        integer(HID_T)   :: ID_FILE, ID_G_summary, ID_DS, ID_G_sum_matp, ID_G_sum_derq, ID_G_sum_numer
        integer(HSIZE_T) :: dims(1)

        filename = trim(dir_output)//"/output.h5"
        ! Open the file
        call h5open_f(error)
            ! Create the root -- When creating the file is always initialized...
            call h5fcreate_f(filename, H5F_ACC_TRUNC_F, ID_FILE, error) !returns the ID_..
                ! Create a group for Summary -- Attached to file via ID_FILE
                call h5gcreate_f(ID_FILE, "Summary", ID_G_summary, error)

                    call h5gcreate_f(ID_G_summary, "Material parameters", ID_G_sum_matp, error)
                        ! Create the dataspace for which the dataset will be written
                        dims = 1
                        call h5screate_simple_f(1, dims, ID_DS, error)
                            ! Create and write the data attach to both group id and dataspace id
                            call write_dataset_double(ID_G_sum_matp, "Radii (km)", lsingle%r_i, dims, ID_DS)
                            data_to_write = lsingle%rho_i(lsingle%n)
                            call write_dataset_double(ID_G_sum_matp, "Density (kgm-3)", data_to_write, dims, ID_DS)
                            ! data_to_write = lsingle%K_i(lsingle%n)
                            ! call write_dataset_double(ID_G_sum_matp, "Bulk modulus (GPa)", data_to_write*Pa_to_GPa, dims, ID_DS) 
                            ! data_to_write = lsingle%mu_i(lsingle%n)
                            ! call write_dataset_double(ID_G_sum_matp, "Shear modulus (GPa)", data_to_write*Pa_to_GPa, dims, ID_DS) 
                            data_to_write = lsingle%eta_i(lsingle%n)
                            call write_dataset_double(ID_G_sum_matp, "Viscosity (Pas)", data_to_write, dims, ID_DS) 

                        call h5sclose_f(ID_DS, error)
                    call h5gclose_f(ID_G_sum_matp, error)

                    call h5gcreate_f(ID_G_summary, "Quantities - Layers, Mass, Orbit", ID_G_sum_derq, error)
                        dims = 1 
                        call h5screate_simple_f(1, dims, ID_DS, error)
                            data_to_write = lsingle%r_i(lsingle%n) - lsingle%r_i(1)
                            call write_dataset_double(ID_G_sum_derq, "Layer thickness", data_to_write, dims, ID_DS) 
                            data_to_write = lsingle%eta_i(1) / (lsingle%rho_i(1)*ang_vel*(lsingle%r_i(lsingle%n))**2)
                            call write_dataset_double(ID_G_sum_derq, "Ekman number", data_to_write, dims, ID_DS) 
                        call h5sclose_f(ID_DS, error)
                        !! Scalar quantities
                        dims=1
                        call h5screate_simple_f(1, dims, ID_DS, error)
                            call write_dataset_double_scalar(ID_G_sum_derq, "Surface gravity (ms-2)", lsingle%g_i(lsingle%n), dims, ID_DS)
                            call write_dataset_double_scalar(ID_G_sum_derq, "Eccentricity", ecc, dims, ID_DS)
                            call write_dataset_double_scalar(ID_G_sum_derq, "Period (s)", rotation_period, dims, ID_DS)
                            call write_dataset_double_scalar(ID_G_sum_derq, "Period (d)", rotation_period/day2sec, dims, ID_DS)
                            call write_dataset_double_scalar(ID_G_sum_derq, "Angular velocity (s-1)", ang_vel, dims, ID_DS)
                        call h5sclose_f(ID_DS, error)
                    call h5gclose_f (ID_G_sum_derq, error)

                    call h5gcreate_f(ID_G_summary, "Numerical Parameters", ID_G_sum_numer, error)
                        dims = 1
                        call h5screate_simple_f(1, dims, ID_DS, error)
                            data_to_write_int = lsingle%n
                            print*, shape(data_to_write_int), 'aa', lsingle%n
                            call write_dataset_int(ID_G_sum_numer, "Number of interface points", data_to_write_int, dims, ID_DS)
                        call h5sclose_f(ID_DS, error)
                        dims = 1
                        call h5screate_simple_f(1, dims, ID_DS, error)
                            call write_dataset_int_scalar(ID_G_sum_numer, "Number of time points", num_time_steps, dims, ID_DS)
                            call write_dataset_double_scalar(ID_G_sum_numer, "Time step (s)", time_step, dims, ID_DS)
                        call h5sclose_f(ID_DS, error)
                    call h5gclose_f (ID_G_sum_numer, error)
                call h5gclose_f(ID_G_summary, error)
            call h5fclose_f (ID_FILE, error)
        call h5close_f(error)

    end subroutine write_summary_lsingle_h5


    subroutine write_radial_profiles_lsingle_h5
        implicit none
        character(300)   :: filename
        integer          :: error, total_points
        integer(HID_T)   :: ID_FILE, ID_G_radprof, ID_G_PROF_INTERFACE, ID_G_PROF_LAYER, ID_DS
        integer(HSIZE_T) :: dims(1)
        integer          :: i

        filename = trim(dir_output)//"/output.h5"

        call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, ID_FILE, error)
                call h5gcreate_f(ID_FILE, "Radial profiles", ID_G_radprof, error)

                    call h5gcreate_f(ID_G_radprof, "Profile (interface)", ID_G_PROF_INTERFACE, error)
                        total_points = lsingle%n
                        dims = total_points
                        call h5screate_simple_f(1, dims, ID_DS, error)
                        call write_dataset_double(ID_G_PROF_INTERFACE, 'Radius (km)'        , [ (lsingle%r_i(1:lsingle%n)*m_to_km) ], dims, ID_DS)
                        call write_dataset_double(ID_G_PROF_INTERFACE, 'Density (kgm-3)'    , [ (lsingle%rho_i(1:lsingle%n)) ], dims, ID_DS)
                        ! call write_dataset_double(ID_G_PROF_INTERFACE, 'Bulk modulus (GPa)' , [ (lsingle%K_i(1:lsingle%n)*Pa_to_GPa) ], dims, ID_DS)
                        ! call write_dataset_double(ID_G_PROF_INTERFACE, 'Shear modulus (GPa)', [ (lsingle%mu_i(1:lsingle%n)*Pa_to_GPa) ], dims, ID_DS)
                        call write_dataset_double(ID_G_PROF_INTERFACE, 'Viscosity (Pa s)'   , [ (lsingle%eta_i(1:lsingle%n)) ], dims, ID_DS)
                        ! call write_dataset_double(ID_G_PROF_INTERFACE, 'Mass (kg)'          , [ (lsingle%mass_i(1:lsingle%n)) ], dims, ID_DS)
                        call write_dataset_double(ID_G_PROF_INTERFACE, 'Gravitational acc. (ms-2)', [ (lsingle%g_i(1:lsingle%n)) ], dims, ID_DS)
                        call h5sclose_f(ID_DS, error)
                    call h5gclose_f (ID_G_PROF_INTERFACE, error)

                    call h5gcreate_f(ID_G_radprof, "Profile (layer)", ID_G_PROF_LAYER, error)
                        total_points = lsingle%n + 1
                        dims = total_points
                        call h5screate_simple_f(1, dims, ID_DS, error)
                        call write_dataset_double(ID_G_PROF_LAYER, 'Radius (km)'        , [ (lsingle%r_l(1:lsingle%n+1)*m_to_km) ], dims, ID_DS)
                        call write_dataset_double(ID_G_PROF_LAYER, 'Density (kgm-3)'    , [ (lsingle%rho_l(1:lsingle%n+1)) ], dims, ID_DS)
                        ! call write_dataset_double(ID_G_PROF_LAYER, 'Mass (kg)'          , [ (lsingle%mass_l(1:lsingle%n+1)) ], dims, ID_DS)
                        ! call write_dataset_double(ID_G_PROF_LAYER, 'Gravitational acc. (ms-2)', [ (lsingle%g_l(1:lsingle%n+1)) ], dims, ID_DS)
                        call h5sclose_f(ID_DS, error)
                    call h5gclose_f (ID_G_PROF_LAYER, error)
                call h5gclose_f(ID_G_radprof, error)
            call h5fclose_f(ID_FILE, error)
        call h5open_f(error)
    end subroutine write_radial_profiles_lsingle_h5



    subroutine write_summary_h5
        implicit none
        character(300)   :: filename
        real(dp)         :: data_to_write(1:numoflayers)
        integer          :: i, data_to_write_int(1:numoflayers)
        !! Variables regarding the hdf5
        integer          :: error
        integer(HID_T)   :: ID_FILE, ID_G_summary, ID_DS, ID_G_sum_matp, ID_G_sum_derq, ID_G_sum_numer
        integer(HSIZE_T) :: dims(1)

        filename = trim(dir_output)//"/output.h5"
        ! Open the file
        call h5open_f(error)
            ! Create the root -- When creating the file is always initialized...
            call h5fcreate_f(filename, H5F_ACC_TRUNC_F, ID_FILE, error) !returns the ID_..
                ! Create a group for Summary -- Attached to file via ID_FILE
                call h5gcreate_f(ID_FILE, "Summary", ID_G_summary, error)

                    call h5gcreate_f(ID_G_summary, "Material parameters", ID_G_sum_matp, error)
                        ! Create the dataspace for which the dataset will be written
                        dims = numoflayers
                        call h5screate_simple_f(1, dims, ID_DS, error)
                            ! Create and write the data attach to both group id and dataspace id
                            call write_dataset_double(ID_G_sum_matp, "Radii (km)", rad*m_to_km, dims, ID_DS)
                            data_to_write = (/ (layer(i)%rho_i(layer(i)%n), i = 1, numoflayers) /)
                            call write_dataset_double(ID_G_sum_matp, "Density (kgm-3)", data_to_write, dims, ID_DS)
                            data_to_write = (/ (layer(i)%K_i(layer(i)%n), i = 1, numoflayers) /)
                            call write_dataset_double(ID_G_sum_matp, "Bulk modulus (GPa)", data_to_write*Pa_to_GPa, dims, ID_DS) 
                            data_to_write = (/ (layer(i)%mu_i(layer(i)%n), i = 1, numoflayers) /)
                            call write_dataset_double(ID_G_sum_matp, "Shear modulus (GPa)", data_to_write*Pa_to_GPa, dims, ID_DS) 
                            data_to_write = (/ (layer(i)%eta_i(layer(i)%n), i = 1, numoflayers) /)
                            call write_dataset_double(ID_G_sum_matp, "Viscosity (Pas)", data_to_write, dims, ID_DS) 

                        call h5sclose_f(ID_DS, error)
                    call h5gclose_f(ID_G_sum_matp, error)

                    call h5gcreate_f(ID_G_summary, "Quantities - Layers, Mass, Orbit", ID_G_sum_derq, error)
                        dims = numoflayers 
                        call h5screate_simple_f(1, dims, ID_DS, error)
                            data_to_write = (/ (rad(i)-rad(i-1), i = 1, numoflayers)/)
                            call write_dataset_double(ID_G_sum_derq, "Layer thickness (km)", data_to_write*m_to_km, dims, ID_DS) 
                            data_to_write = (/ (layer(i)%eta_i(1) / (2.0_dp*layer(i)%rho_i(1)*ang_vel*(rad(i)-rad(i-1))**2), i = 1, numoflayers)/)
                            call write_dataset_double(ID_G_sum_derq, "Ekman number", data_to_write, dims, ID_DS) 
                        call h5sclose_f(ID_DS, error)
                        !! Scalar quantities
                        dims=1
                        call h5screate_simple_f(1, dims, ID_DS, error)
                            call write_dataset_double_scalar(ID_G_sum_derq, "Mass (kg)", layer(numoflayers)%mass_i(layer(numoflayers)%n), dims, ID_DS)
                            call write_dataset_double_scalar(ID_G_sum_derq, "MOI factor", moi, dims, ID_DS)
                            call write_dataset_double_scalar(ID_G_sum_derq, "Surface gravity (ms-2)", layer(numoflayers)%g_i(layer(numoflayers)%n), dims, ID_DS)
                            call write_dataset_double_scalar(ID_G_sum_derq, "Eccentricity", ecc, dims, ID_DS)
                            call write_dataset_double_scalar(ID_G_sum_derq, "Period (s)", rotation_period, dims, ID_DS)
                            call write_dataset_double_scalar(ID_G_sum_derq, "Period (d)", rotation_period/day2sec, dims, ID_DS)
                            call write_dataset_double_scalar(ID_G_sum_derq, "Angular velocity (s-1)", ang_vel, dims, ID_DS)
                            call write_dataset_double_scalar(ID_G_sum_derq, "Mean density (kgm-3)", layer(numoflayers)%mass_i(layer(numoflayers)%n)/(4.0_dp*pi*(rad(numoflayers)**3)/3.0_dp), dims, ID_DS)
                        call h5sclose_f(ID_DS, error)
                    call h5gclose_f (ID_G_sum_derq, error)

                    call h5gcreate_f(ID_G_summary, "Numerical Parameters", ID_G_sum_numer, error)
                        dims = numoflayers
                        call h5screate_simple_f(1, dims, ID_DS, error)
                            data_to_write_int = (/ (layer(i)%n, i = 1, numoflayers)/)
                            call write_dataset_int(ID_G_sum_numer, "Number of interface points", data_to_write_int, dims, ID_DS)
                        call h5sclose_f(ID_DS, error)
                        dims = 1
                        call h5screate_simple_f(1, dims, ID_DS, error)
                            call write_dataset_int_scalar(ID_G_sum_numer, "Number of time points", num_time_steps, dims, ID_DS)
                            call write_dataset_double_scalar(ID_G_sum_numer, "Time step (s)", time_step, dims, ID_DS)
                        call h5sclose_f(ID_DS, error)
                    call h5gclose_f (ID_G_sum_numer, error)
                call h5gclose_f(ID_G_summary, error)
            call h5fclose_f (ID_FILE, error)
        call h5close_f(error)

    end subroutine write_summary_h5


    subroutine write_dataset_double_scalar(group_id, name, data, dims, space_id)
        integer(HID_T), intent(in)   :: group_id, space_id
        character(*), intent(in)     :: name
        real(dp), intent(in)         :: data
        integer(HSIZE_T), intent(in) :: dims(1)
        integer(HID_T)               :: dset_id
        integer                      :: error

        call h5dcreate_f(group_id, name, H5T_NATIVE_DOUBLE, space_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, error)
        call h5dclose_f(dset_id, error)
    end subroutine write_dataset_double_scalar

    subroutine write_dataset_double(group_id, name, data, dims, space_id)
        integer(HID_T), intent(in)   :: group_id, space_id
        character(*), intent(in)     :: name
        real(dp), intent(in)         :: data(*)
        integer(HSIZE_T), intent(in) :: dims(1)
        integer(HID_T)               :: dset_id
        integer                      :: error

        call h5dcreate_f(group_id, name, H5T_NATIVE_DOUBLE, space_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, error)
        call h5dclose_f(dset_id, error)
    end subroutine write_dataset_double

    subroutine write_dataset_int(group_id, name, data, dims, space_id)
        integer(HID_T), intent(in) :: group_id, space_id
        character(*), intent(in)   :: name
        integer, intent(in)        :: data(*)
        integer(HSIZE_T), intent(in) :: dims(1)
        integer(HID_T)             :: dset_id
        integer                    :: error

        call h5dcreate_f(group_id, name, H5T_NATIVE_INTEGER, space_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dims, error)
        call h5dclose_f(dset_id, error)
    end subroutine write_dataset_int

    subroutine write_dataset_int_scalar(group_id, name, data, dims, space_id)
        integer(HID_T), intent(in) :: group_id, space_id
        character(*), intent(in)   :: name
        integer, intent(in)        :: data
        integer(HSIZE_T), intent(in) :: dims(1)
        integer(HID_T)             :: dset_id
        integer                    :: error

        call h5dcreate_f(group_id, name, H5T_NATIVE_INTEGER, space_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dims, error)
        call h5dclose_f(dset_id, error)
    end subroutine write_dataset_int_scalar

    subroutine write_radial_profiles_h5
        implicit none
        character(300)   :: filename
        integer          :: error, total_points
        integer(HID_T)   :: ID_FILE, ID_G_radprof, ID_G_PROF_INTERFACE, ID_G_PROF_LAYER, ID_DS
        integer(HSIZE_T) :: dims(1)
        integer          :: i

        filename = trim(dir_output)//"/output.h5"

        call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, ID_FILE, error)
                call h5gcreate_f(ID_FILE, "Radial profiles", ID_G_radprof, error)

                    call h5gcreate_f(ID_G_radprof, "Profile (interface)", ID_G_PROF_INTERFACE, error)
                        total_points = sum(layer(:)%n)
                        dims = total_points
                        call h5screate_simple_f(1, dims, ID_DS, error)
                        call write_dataset_double(ID_G_PROF_INTERFACE, 'Radius (km)'        , [ (layer(i)%r_i(1:layer(i)%n)*m_to_km, i = 1, numoflayers) ], dims, ID_DS)
                        call write_dataset_double(ID_G_PROF_INTERFACE, 'Density (kgm-3)'    , [ (layer(i)%rho_i(1:layer(i)%n), i = 1, numoflayers) ], dims, ID_DS)
                        call write_dataset_double(ID_G_PROF_INTERFACE, 'Bulk modulus (GPa)' , [ (layer(i)%K_i(1:layer(i)%n)*Pa_to_GPa, i = 1, numoflayers) ], dims, ID_DS)
                        call write_dataset_double(ID_G_PROF_INTERFACE, 'Shear modulus (GPa)', [ (layer(i)%mu_i(1:layer(i)%n)*Pa_to_GPa, i = 1, numoflayers) ], dims, ID_DS)
                        call write_dataset_double(ID_G_PROF_INTERFACE, 'Viscosity (Pa s)'   , [ (layer(i)%eta_i(1:layer(i)%n), i = 1, numoflayers) ], dims, ID_DS)
                        call write_dataset_double(ID_G_PROF_INTERFACE, 'Mass (kg)'          , [ (layer(i)%mass_i(1:layer(i)%n), i = 1, numoflayers) ], dims, ID_DS)
                        call write_dataset_double(ID_G_PROF_INTERFACE, 'Gravitational acc. (ms-2)', [ (layer(i)%g_i(1:layer(i)%n), i = 1, numoflayers) ], dims, ID_DS)
                        call h5sclose_f(ID_DS, error)
                    call h5gclose_f (ID_G_PROF_INTERFACE, error)

                    call h5gcreate_f(ID_G_radprof, "Profile (layer)", ID_G_PROF_LAYER, error)
                        total_points = sum(layer(:)%n) + numoflayers
                        dims = total_points
                        call h5screate_simple_f(1, dims, ID_DS, error)
                        call write_dataset_double(ID_G_PROF_LAYER, 'Radius (km)'        , [ (layer(i)%r_l(1:layer(i)%n+1)*m_to_km, i = 1, numoflayers) ], dims, ID_DS)
                        call write_dataset_double(ID_G_PROF_LAYER, 'Density (kgm-3)'    , [ (layer(i)%rho_l(1:layer(i)%n+1), i = 1, numoflayers) ], dims, ID_DS)
                        call write_dataset_double(ID_G_PROF_LAYER, 'Mass (kg)'          , [ (layer(i)%mass_l(1:layer(i)%n+1), i = 1, numoflayers) ], dims, ID_DS)
                        call write_dataset_double(ID_G_PROF_LAYER, 'Gravitational acc. (ms-2)', [ (layer(i)%g_l(1:layer(i)%n+1), i = 1, numoflayers) ], dims, ID_DS)
                        call h5sclose_f(ID_DS, error)
                    call h5gclose_f (ID_G_PROF_LAYER, error)
                call h5gclose_f(ID_G_radprof, error)
            call h5fclose_f(ID_FILE, error)
        call h5open_f(error)
    end subroutine write_radial_profiles_h5

    subroutine write_spectrum_h5
        implicit none
        character(300)   :: filename, datasetname
        integer          :: error
        integer(HID_T)   :: ID_FILE, ID_G_spectrum, ID_DS, ID_G_stress, ID_G_velocity
        integer(HSIZE_T) :: dims(1)
        integer          :: i
        real(dp)         :: str_spec(1:Max_jmax, 1:numoflayers)
        real(dp)         :: vel_spec(1:Max_jmax, 1:numoflayers)
        logical          :: exists

        filename = trim(dir_output)//"/output.h5"
        do i = 1, numoflayers
            str_spec(:,i) = layer(i)%spectrum
            vel_spec(:,i) = layer(i)%spectrum_v
        end do
        call h5open_f(error)
        call h5fopen_f(filename, H5F_ACC_RDWR_F, ID_FILE, error)
        call h5gopen_f(ID_FILE, "Spectrum", ID_G_spectrum, error)
        call h5gopen_f(ID_G_spectrum, "Stress_spectrum", ID_G_stress, error)
            dims = jmax
            do i = 1, numoflayers
                write(datasetname, '(A,I0)') "stress_spec_layer_", i
                call h5screate_simple_f(1, dims, ID_DS, error)
                call h5lexists_f(ID_G_stress, datasetname, exists, error)
                if (exists) call h5ldelete_f(ID_G_stress, datasetname, error)
                call write_dataset_double(ID_G_stress, datasetname, layer(i)%spectrum, dims, ID_DS)
                call h5sclose_f(ID_DS, error)
            end do
        call h5gclose_f(ID_G_stress, error)
        call h5gopen_f(ID_G_spectrum, "Velocity_spectrum", ID_G_velocity, error)
            do i = 1, numoflayers
                write(datasetname, '(A,I0)') "velocity_spec_layer_", i
                call h5screate_simple_f(1, dims, ID_DS, error)
                call h5lexists_f(ID_G_velocity, datasetname, exists, error)
                if (exists) call h5ldelete_f(ID_G_velocity, datasetname, error)
                call write_dataset_double(ID_G_velocity, datasetname, layer(i)%spectrum_v, dims, ID_DS)
                call h5sclose_f(ID_DS, error)
            end do
        call h5gclose_f(ID_G_velocity, error)
        call h5gclose_f(ID_G_spectrum, error)
        call h5fclose_f(ID_FILE, error)
        call h5close_f(error)
    end subroutine write_spectrum_h5

    subroutine generate_timeseries_space_h5
        implicit none
        character(len=300) :: filename, datasetname
        integer(HID_T)     :: ID_FILE, ID_DS, ID_PLIST, ID_dissip, ID_time, ID_k2, ID_G_timeseries, ID_G_dissipation, ID_G_potential, ID_G_displacement
        integer(HID_T)     :: ID_G_spectrum, ID_G_spec_str, ID_G_spec_vel, ID_G_period
        integer(hsize_t), dimension(1) :: dims, maxdims, chunk_dims
        integer                        :: error, il

        filename=trim(dir_output)//"/output.h5"
        call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, ID_FILE, error)
            ! Create extendable datasets
            dims = 0
            maxdims = H5S_UNLIMITED_F
            chunk_dims = 1000
            call h5screate_simple_f(1, dims, ID_DS, error, maxdims)
            call h5gcreate_f(ID_FILE, "Time_series", ID_G_timeseries, error)
                call h5pcreate_f(H5P_DATASET_CREATE_F, ID_PLIST, error)
                call h5pset_chunk_f(ID_PLIST, 1, chunk_dims, error)

                call h5dcreate_f(ID_G_timeseries, "time", H5T_NATIVE_DOUBLE, ID_DS, ID_time, error, ID_PLIST)
                call h5dclose_f(ID_time, error)
                call h5gcreate_f(ID_G_timeseries, "Dissipation", ID_G_dissipation, error)
                do il = 1, numoflayers
                    write(datasetname, '(A,I0)') 'dissipation_layer_', il
                    call h5dcreate_f(ID_G_dissipation, datasetname, H5T_NATIVE_DOUBLE, ID_DS, ID_dissip, error, ID_PLIST)
                    call h5dclose_f(ID_dissip, error)
                end do
                call h5gclose_f(ID_G_dissipation, error)
                call h5gcreate_f(ID_G_timeseries, "Induced_potential", ID_G_potential, error)
                do il = jm(2,0), jm(2,2)
                    if (jmindx(il, 3) == 0) then
                        write(datasetname, '(A,I0, A,I0)') 'induced_potential_j', jmindx(il, 2), '_m', jmindx(il, 3)
                        call h5dcreate_f(ID_G_potential, datasetname, H5T_NATIVE_DOUBLE, ID_DS, ID_k2, error, ID_PLIST)
                        call h5dclose_f(ID_k2, error)
                    else
                        write(datasetname, '(A,I0, A,I0)') 'induced_potential_RE_j', jmindx(il, 2), '_m', jmindx(il, 3)
                        call h5dcreate_f(ID_G_potential, datasetname, H5T_NATIVE_DOUBLE, ID_DS, ID_k2, error, ID_PLIST)
                        call h5dclose_f(ID_k2, error)
                        write(datasetname, '(A,I0, A,I0)') 'induced_potential_IM_j', jmindx(il, 2), '_m', jmindx(il, 3)
                        call h5dcreate_f(ID_G_potential, datasetname, H5T_NATIVE_DOUBLE, ID_DS, ID_k2, error, ID_PLIST)
                        call h5dclose_f(ID_k2, error)
                    end if
                end do
                ! call h5pclose_f(ID_PLIST, error)
                call h5gclose_f(ID_G_potential, error)

                call h5gcreate_f(ID_G_timeseries, "Top displacement", ID_G_displacement, error)
                do il = jm(2,0), jm(2,2)
                    if (jmindx(il, 3) == 0) then
                        write(datasetname, '(A,I0, A,I0)') 'top_displacement_j', jmindx(il, 2), '_m', jmindx(il, 3)
                        call h5dcreate_f(ID_G_displacement, datasetname, H5T_NATIVE_DOUBLE, ID_DS, ID_k2, error, ID_PLIST)
                        call h5dclose_f(ID_k2, error)
                    else
                        write(datasetname, '(A,I0, A,I0)') 'top_displacement_RE_j', jmindx(il, 2), '_m', jmindx(il, 3)
                        call h5dcreate_f(ID_G_displacement, datasetname, H5T_NATIVE_DOUBLE, ID_DS, ID_k2, error, ID_PLIST)
                        call h5dclose_f(ID_k2, error)
                        write(datasetname, '(A,I0, A,I0)') 'top_displacement_IM_j', jmindx(il, 2), '_m', jmindx(il, 3)
                        call h5dcreate_f(ID_G_displacement, datasetname, H5T_NATIVE_DOUBLE, ID_DS, ID_k2, error, ID_PLIST)
                        call h5dclose_f(ID_k2, error)
                    end if
                end do
                call h5pclose_f(ID_PLIST, error)
                call h5gclose_f(ID_G_displacement, error)
            call h5gclose_f(ID_G_timeseries, error)
            call h5sclose_f(ID_DS, error)
            call h5fclose_f(ID_FILE, error)
        call h5close_f(error)

        call h5open_f(error)
        call h5fopen_f(filename, H5F_ACC_RDWR_F, ID_FILE, error)
            call h5gcreate_f(ID_FILE, "Spectrum", ID_G_spectrum, error)
                call h5gcreate_f(ID_G_spectrum, "Stress_spectrum", ID_G_spec_str, error)
                call h5gclose_f(ID_G_spec_str, error)
                call h5gcreate_f(ID_G_spectrum, "Velocity_spectrum", ID_G_spec_vel, error)
                call h5gclose_f(ID_G_spec_vel, error)
            call h5gclose_f(ID_G_spectrum, error)
        call h5fclose_f(ID_FILE, error)
        call h5close_f(error)

        call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, ID_FILE, error)
            call h5screate_simple_f(1, dims, ID_DS, error, maxdims)
                dims = 0
                maxdims = H5S_UNLIMITED_F
                chunk_dims = 1
                call h5gcreate_f(ID_FILE, "Periodic", ID_G_period, error)
                call h5pcreate_f(H5P_DATASET_CREATE_F, ID_PLIST, error)
                call h5pset_chunk_f(ID_PLIST, 1, chunk_dims, error)
                call h5dcreate_f(ID_G_period, "period", H5T_NATIVE_DOUBLE, ID_DS, ID_time, error, ID_PLIST)
                call h5dclose_f(ID_time, error)
                    call h5gcreate_f(ID_G_period, "Dissipation", ID_G_dissipation, error)
                    do il = 1, numoflayers
                        write(datasetname, '(A,I0)') 'average_dissipation_layer_', il
                        call h5dcreate_f(ID_G_dissipation, datasetname, H5T_NATIVE_DOUBLE, ID_DS, ID_dissip, error, ID_PLIST)
                        call h5dclose_f(ID_dissip, error)
                    end do
                    call h5gclose_f(ID_G_dissipation, error)
                    call h5gcreate_f(ID_G_period, "Love_numbers", ID_G_potential, error)
                    do il = jm(2,0), jm(2,2)
                        if (jmindx(il, 3) == 0) then
                            write(datasetname, '(A,I0, A,I0)') 'k_j', jmindx(il, 2), '_m', jmindx(il, 3)
                            call h5dcreate_f(ID_G_potential, datasetname, H5T_NATIVE_DOUBLE, ID_DS, ID_k2, error, ID_PLIST)
                            call h5dclose_f(ID_k2, error)
                            write(datasetname, '(A,I0, A,I0)') 'h_j', jmindx(il, 2), '_m', jmindx(il, 3)
                            call h5dcreate_f(ID_G_potential, datasetname, H5T_NATIVE_DOUBLE, ID_DS, ID_k2, error, ID_PLIST)
                            call h5dclose_f(ID_k2, error)
                            write(datasetname, '(A,I0, A,I0)') 'l_j', jmindx(il, 2), '_m', jmindx(il, 3)
                            call h5dcreate_f(ID_G_potential, datasetname, H5T_NATIVE_DOUBLE, ID_DS, ID_k2, error, ID_PLIST)
                            call h5dclose_f(ID_k2, error)
                        else
                            write(datasetname, '(A,I0, A,I0)') 'RE_k_j', jmindx(il, 2), '_m', jmindx(il, 3)
                            call h5dcreate_f(ID_G_potential, datasetname, H5T_NATIVE_DOUBLE, ID_DS, ID_k2, error, ID_PLIST)
                            call h5dclose_f(ID_k2, error)
                            write(datasetname, '(A,I0, A,I0)') 'RE_h_j', jmindx(il, 2), '_m', jmindx(il, 3)
                            call h5dcreate_f(ID_G_potential, datasetname, H5T_NATIVE_DOUBLE, ID_DS, ID_k2, error, ID_PLIST)
                            call h5dclose_f(ID_k2, error)
                            write(datasetname, '(A,I0, A,I0)') 'RE_l_j', jmindx(il, 2), '_m', jmindx(il, 3)
                            call h5dcreate_f(ID_G_potential, datasetname, H5T_NATIVE_DOUBLE, ID_DS, ID_k2, error, ID_PLIST)
                            call h5dclose_f(ID_k2, error)
                            write(datasetname, '(A,I0, A,I0)') 'IM_k_j', jmindx(il, 2), '_m', jmindx(il, 3)
                            call h5dcreate_f(ID_G_potential, datasetname, H5T_NATIVE_DOUBLE, ID_DS, ID_k2, error, ID_PLIST)
                            call h5dclose_f(ID_k2, error)
                            write(datasetname, '(A,I0, A,I0)') 'IM_h_j', jmindx(il, 2), '_m', jmindx(il, 3)
                            call h5dcreate_f(ID_G_potential, datasetname, H5T_NATIVE_DOUBLE, ID_DS, ID_k2, error, ID_PLIST)
                            call h5dclose_f(ID_k2, error)
                            write(datasetname, '(A,I0, A,I0)') 'IM_l_j', jmindx(il, 2), '_m', jmindx(il, 3)
                            call h5dcreate_f(ID_G_potential, datasetname, H5T_NATIVE_DOUBLE, ID_DS, ID_k2, error, ID_PLIST)
                            call h5dclose_f(ID_k2, error)
                        end if
                    end do
                    call h5gclose_f(ID_G_potential, error)
                call h5pclose_f(ID_PLIST, error)
                call h5gclose_f(ID_G_period, error)
            call h5sclose_f(ID_DS, error)
            call h5fclose_f(ID_FILE, error)
        call h5close_f(error)
    end subroutine generate_timeseries_space_h5

    subroutine write_time_series_h5(it2)
        implicit none
        character(len=300) :: filename, datasetname
        integer :: nsteps, offset
        integer :: error
        integer(hid_t) :: file_id, dspace_id, memspace_id, ID_G_timeseries, ID_G_dissipation, ID_G_potential, ID_G_displacement
        integer(hid_t) :: dset_time, dset_diss, dset_k2
        integer(hsize_t), dimension(1) :: start, count, new_size
        integer :: i, it2, nn

        nn = layer(numoflayers)%n
        offset = it2 - num_time_steps
        nsteps = num_time_steps
        new_size = offset + nsteps
        filename = trim(dir_output)//'/output.h5'

        call h5open_f(error)
        ! Open or create file
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
        call h5gopen_f(file_id, 'Time_series', ID_G_timeseries, error)
        ! Reopen datasets
        call h5dopen_f(ID_G_timeseries, "time", dset_time, error)
        call h5dset_extent_f(dset_time, new_size, error)
        start = offset
        count = nsteps

        call h5dget_space_f(dset_time, dspace_id, error)
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, start, count, error)
        call h5screate_simple_f(1, count, memspace_id, error)
        call h5dwrite_f(dset_time, H5T_NATIVE_DOUBLE, time, count, error, memspace_id, dspace_id)
        call h5dclose_f(dset_time, error)
        call h5gopen_f(ID_G_timeseries, "Dissipation", ID_G_dissipation, error)
        do i = 1, numoflayers
            write(datasetname, '(A,I0)') 'dissipation_layer_', i
            call h5dopen_f(ID_G_dissipation, datasetname, dset_diss, error)
            call h5dset_extent_f(dset_diss, new_size, error)
            call h5dwrite_f(dset_diss, H5T_NATIVE_DOUBLE, layer(i)%dissipation_time*W_to_GW, count, error, memspace_id, dspace_id)
            call h5dclose_f(dset_diss, error)
        end do
        call h5gclose_f(ID_G_dissipation, error)
        call h5gopen_f(ID_G_timeseries, "Induced_potential", ID_G_potential, error)
        do i = jm(2,0), jm(2,2)

            if (jmindx(i, 3) == 0) then
                write(datasetname, '(A,I0, A,I0)') 'induced_potential_j', jmindx(i, 2), '_m', jmindx(i, 3)
                call h5dopen_f(ID_G_potential, datasetname, dset_k2, error)
                call h5dset_extent_f(dset_k2, new_size, error)
                call h5dwrite_f(dset_k2, H5T_NATIVE_DOUBLE, real(indpot(i, :)), count, error, memspace_id, dspace_id)
                call h5dclose_f(dset_k2, error)
            else
                write(datasetname, '(A,I0, A,I0)') 'induced_potential_RE_j', jmindx(i, 2), '_m', jmindx(i, 3)
                call h5dopen_f(ID_G_potential, datasetname, dset_k2, error)
                call h5dset_extent_f(dset_k2, new_size, error)
                call h5dwrite_f(dset_k2, H5T_NATIVE_DOUBLE, real(indpot(i, :)), count, error, memspace_id, dspace_id)
                call h5dclose_f(dset_k2, error)
                write(datasetname, '(A,I0, A,I0)') 'induced_potential_IM_j', jmindx(i, 2), '_m', jmindx(i, 3)
                call h5dopen_f(ID_G_potential, datasetname, dset_k2, error)
                call h5dset_extent_f(dset_k2, new_size, error)
                call h5dwrite_f(dset_k2, H5T_NATIVE_DOUBLE, imag(indpot(i, :)), count, error, memspace_id, dspace_id)
                call h5dclose_f(dset_k2, error)
            end if
        end do
        call h5gclose_f(ID_G_potential, error)

        call h5gopen_f(ID_G_timeseries, "Top displacement", ID_G_displacement, error)
        do i = jm(2,0), jm(2,2)

            if (jmindx(i, 3) == 0) then
                write(datasetname, '(A,I0, A,I0)') 'top_displacement_j', jmindx(i, 2), '_m', jmindx(i, 3)
                call h5dopen_f(ID_G_displacement, datasetname, dset_k2, error)
                call h5dset_extent_f(dset_k2, new_size, error)
                call h5dwrite_f(dset_k2, H5T_NATIVE_DOUBLE, real(topur(i, :)), count, error, memspace_id, dspace_id)
                call h5dclose_f(dset_k2, error)
            else
                write(datasetname, '(A,I0, A,I0)') 'top_displacement_RE_j', jmindx(i, 2), '_m', jmindx(i, 3)
                call h5dopen_f(ID_G_displacement, datasetname, dset_k2, error)
                call h5dset_extent_f(dset_k2, new_size, error)
                call h5dwrite_f(dset_k2, H5T_NATIVE_DOUBLE, real(topur(i, :)), count, error, memspace_id, dspace_id)
                call h5dclose_f(dset_k2, error)
                write(datasetname, '(A,I0, A,I0)') 'top_displacement_IM_j', jmindx(i, 2), '_m', jmindx(i, 3)
                call h5dopen_f(ID_G_displacement, datasetname, dset_k2, error)
                call h5dset_extent_f(dset_k2, new_size, error)
                call h5dwrite_f(dset_k2, H5T_NATIVE_DOUBLE, imag(topur(i, :)), count, error, memspace_id, dspace_id)
                call h5dclose_f(dset_k2, error)
            end if
        end do
        call h5gclose_f(ID_G_displacement, error)

        call h5sclose_f(memspace_id, error)
        call h5sclose_f(dspace_id, error)
        call h5gclose_f(ID_G_timeseries, error)
        call h5fclose_f(file_id, error)
        call h5close_f (error)
    end subroutine write_time_series_h5

    subroutine write_periodic_h5(ip)
        implicit none
        character(len=300) :: filename, datasetname
        integer :: nsteps, offset
        integer :: error
        integer(hid_t) :: file_id, dspace_id, memspace_id, ID_G_timeseries, ID_G_dissipation, ID_G_potential
        integer(hid_t) :: dset_time, dset_diss, dset_k2
        integer(hsize_t), dimension(1) :: start, count, new_size
        integer :: i, nn, ip

        nn = layer(numoflayers)%n
        offset = ip - 1
        nsteps = 1
        new_size = offset + nsteps
        filename = trim(dir_output)//'/output.h5'

        call h5open_f(error)
        ! Open or create file
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
        call h5gopen_f(file_id, 'Periodic', ID_G_timeseries, error)
        ! Reopen datasets
        call h5dopen_f(ID_G_timeseries, "period", dset_time, error)
        call h5dset_extent_f(dset_time, new_size, error)
        start = offset
        count = nsteps

        call h5dget_space_f(dset_time, dspace_id, error)
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, start, count, error)
        call h5screate_simple_f(1, count, memspace_id, error)
        call h5dwrite_f(dset_time, H5T_NATIVE_DOUBLE, real(ip,dp), count, error, memspace_id, dspace_id)
        call h5dclose_f(dset_time, error)
        call h5gopen_f(ID_G_timeseries, "Dissipation", ID_G_dissipation, error)
        do i = 1, numoflayers
            write(datasetname, '(A,I0)') 'average_dissipation_layer_', i
            call h5dopen_f(ID_G_dissipation, datasetname, dset_diss, error)
            call h5dset_extent_f(dset_diss, new_size, error)
            call h5dwrite_f(dset_diss, H5T_NATIVE_DOUBLE, layer(i)%average_dissip*W_to_GW, count, error, memspace_id, dspace_id)
            call h5dclose_f(dset_diss, error)
        end do
        call h5gclose_f(ID_G_dissipation, error)
        call h5gopen_f(ID_G_timeseries, "Love_numbers", ID_G_potential, error)
        do i = jm(2,0), jm(2,2)
            if (jmindx(i, 3) == 0) then
                write(datasetname, '(A,I0, A,I0)') 'k_j', jmindx(i, 2), '_m', jmindx(i, 3)
                call h5dopen_f(ID_G_potential, datasetname, dset_k2, error)
                call h5dset_extent_f(dset_k2, new_size, error)
                call h5dwrite_f(dset_k2, H5T_NATIVE_DOUBLE, real(love_k(i)), count, error, memspace_id, dspace_id)
                call h5dclose_f(dset_k2, error)
                write(datasetname, '(A,I0, A,I0)') 'h_j', jmindx(i, 2), '_m', jmindx(i, 3)
                call h5dopen_f(ID_G_potential, datasetname, dset_k2, error)
                call h5dset_extent_f(dset_k2, new_size, error)
                call h5dwrite_f(dset_k2, H5T_NATIVE_DOUBLE, real(love_h(i)), count, error, memspace_id, dspace_id)
                call h5dclose_f(dset_k2, error)
                write(datasetname, '(A,I0, A,I0)') 'l_j', jmindx(i, 2), '_m', jmindx(i, 3)
                call h5dopen_f(ID_G_potential, datasetname, dset_k2, error)
                call h5dset_extent_f(dset_k2, new_size, error)
                call h5dwrite_f(dset_k2, H5T_NATIVE_DOUBLE, real(love_l(i)), count, error, memspace_id, dspace_id)
                call h5dclose_f(dset_k2, error)
            else
                write(datasetname, '(A,I0, A,I0)') 'RE_k_j', jmindx(i, 2), '_m', jmindx(i, 3)
                call h5dopen_f(ID_G_potential, datasetname, dset_k2, error)
                call h5dset_extent_f(dset_k2, new_size, error)
                call h5dwrite_f(dset_k2, H5T_NATIVE_DOUBLE, real(love_k(i)), count, error, memspace_id, dspace_id)
                call h5dclose_f(dset_k2, error)
                write(datasetname, '(A,I0, A,I0)') 'RE_h_j', jmindx(i, 2), '_m', jmindx(i, 3)
                call h5dopen_f(ID_G_potential, datasetname, dset_k2, error)
                call h5dset_extent_f(dset_k2, new_size, error)
                call h5dwrite_f(dset_k2, H5T_NATIVE_DOUBLE, real(love_h(i)), count, error, memspace_id, dspace_id)
                call h5dclose_f(dset_k2, error)
                write(datasetname, '(A,I0, A,I0)') 'RE_l_j', jmindx(i, 2), '_m', jmindx(i, 3)
                call h5dopen_f(ID_G_potential, datasetname, dset_k2, error)
                call h5dset_extent_f(dset_k2, new_size, error)
                call h5dwrite_f(dset_k2, H5T_NATIVE_DOUBLE, real(love_l(i)), count, error, memspace_id, dspace_id)
                call h5dclose_f(dset_k2, error)

                write(datasetname, '(A,I0, A,I0)') 'IM_k_j', jmindx(i, 2), '_m', jmindx(i, 3)
                call h5dopen_f(ID_G_potential, datasetname, dset_k2, error)
                call h5dset_extent_f(dset_k2, new_size, error)
                call h5dwrite_f(dset_k2, H5T_NATIVE_DOUBLE, imag(love_k(i)), count, error, memspace_id, dspace_id)
                call h5dclose_f(dset_k2, error)
                write(datasetname, '(A,I0, A,I0)') 'IM_h_j', jmindx(i, 2), '_m', jmindx(i, 3)
                call h5dopen_f(ID_G_potential, datasetname, dset_k2, error)
                call h5dset_extent_f(dset_k2, new_size, error)
                call h5dwrite_f(dset_k2, H5T_NATIVE_DOUBLE, imag(love_h(i)), count, error, memspace_id, dspace_id)
                call h5dclose_f(dset_k2, error)
                write(datasetname, '(A,I0, A,I0)') 'IM_l_j', jmindx(i, 2), '_m', jmindx(i, 3)
                call h5dopen_f(ID_G_potential, datasetname, dset_k2, error)
                call h5dset_extent_f(dset_k2, new_size, error)
                call h5dwrite_f(dset_k2, H5T_NATIVE_DOUBLE, imag(love_l(i)), count, error, memspace_id, dspace_id)
                call h5dclose_f(dset_k2, error)
            end if
        end do
        call h5gclose_f(ID_G_potential, error)

        call h5sclose_f(memspace_id, error)
        call h5sclose_f(dspace_id, error)
        call h5gclose_f(ID_G_timeseries, error)
        call h5fclose_f(file_id, error)
        call h5close_f (error)
    end subroutine write_periodic_h5

    subroutine generate_velocity_shc_h5
        implicit none
        integer            :: error, i, il
        integer(HID_T)     :: ID_file, ID_DS, ID_PLIST, ID_time
        integer(HSIZE_T)   :: dims(1), maxdims(1), chunk_dims(1)
        character(len=400) :: filename, dataname

        write(filename, '(A)') trim(dir_output)//"/sh_coeffs/velocity.h5"
        call h5open_f(error)
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, ID_file, error)
        dims = 0
        maxdims = H5S_UNLIMITED_F
        chunk_dims = 1
        call h5screate_simple_f(1, dims, ID_DS, error, maxdims)
        call h5pcreate_f(H5P_DATASET_CREATE_F, ID_PLIST, error)
        call h5pset_chunk_f(ID_PLIST, 1, chunk_dims, error)
        call h5dcreate_f(ID_file, "time", H5T_NATIVE_INTEGER, ID_DS, ID_time, error, ID_PLIST)
        call h5dclose_f(ID_time, error)
        call h5pclose_f(ID_PLIST, error)
        dims = num_oceans
        call h5screate_simple_f(1, dims, ID_ds, error)
        write(dataname, '(A)') "Ocean_layers"
        call write_dataset_int(ID_file,dataname, layer_ocean(1:num_oceans), dims, ID_ds)
        call h5sclose_f(ID_DS, error)
        do i = 1, num_oceans
            il = layer_ocean(i)
            dims = 1
            call h5screate_simple_f(1, dims, ID_ds, error)
            write(dataname, '(A,I0)') "Npoints_layer_", il
            call write_dataset_int_scalar(ID_file,dataname, layer(il)%n+1, dims, ID_ds)
            call h5sclose_f(ID_DS, error)
        end do
        call h5fclose_f(ID_file, error)
        call h5close_f(error)
    end subroutine generate_velocity_shc_h5


    subroutine write_velocity_shc_h5(it2)
        implicit none
        integer            :: error, it2, i, il
        integer(HID_T)     :: ID_file, ID_G_vel, ID_ds, dset_time, memspace_id
        integer(HSIZE_T)   :: dims(1), dims2(2)
        integer            :: offset, nsteps
        integer(hsize_t), dimension(1) :: start, count, new_size
        character(len=400) :: filename, groupname, dataname
        
        offset = (it2-start_output*num_time_steps)/write_output
        nsteps = 1
        new_size = offset + nsteps
        write(filename, '(A)') trim(dir_output)//"/sh_coeffs/velocity.h5"
        call h5open_f(error)
        call h5fopen_f(filename, H5F_ACC_RDWR_F, ID_file, error)
        call h5dopen_f(ID_file, "time", dset_time, error)
            call h5dset_extent_f(dset_time, new_size, error)
            start = offset
            count = nsteps

            call h5dget_space_f(dset_time, ID_ds, error)
            call h5sselect_hyperslab_f(ID_ds, H5S_SELECT_SET_F, start, count, error)
            call h5screate_simple_f(1, count, memspace_id, error)
            call h5dwrite_f(dset_time, H5T_NATIVE_INTEGER, it2, count, error, memspace_id, ID_ds)
        call h5dclose_f(dset_time, error)
            write(groupname, '(A,I0)') "velocity_", it2
            call h5gcreate_f(ID_file, groupname, ID_G_vel, error)
                dims = 1
                call h5screate_simple_f(1, dims, ID_ds, error)
                call write_dataset_int_scalar(ID_G_vel, "Jmax", jmax, dims, ID_ds)
                call write_dataset_int_scalar(ID_G_vel, "Mmax", mmax, dims, ID_ds)
                call h5sclose_f(ID_DS, error)
                do i = 1, num_oceans
                    il = layer_ocean(i)
                    dims2 = (/compute_vector_harm, layer(il)%n+1/)
                    call h5screate_simple_f(2, dims2, ID_ds, error)
                    write(dataname, '(A,I0,A,I0)') "REvelocity_", il, "_", it2
                    call write_dataset_double(ID_G_vel, dataname, real(layer(il)%vel(1:compute_vector_harm, :)), dims2, ID_ds)
                    write(dataname, '(A,I0,A,I0)') "IMvelocity_", il, "_", it2
                    call write_dataset_double(ID_G_vel, dataname, imag(layer(il)%vel(1:compute_vector_harm, :)), dims2, ID_ds)
                    call h5sclose_f(ID_DS, error)
                end do
            call h5gclose_f(ID_G_vel, error)
        call h5fclose_f(ID_file, error)
        call h5close_f(error)
    end subroutine write_velocity_shc_h5

end module mod_write_h5