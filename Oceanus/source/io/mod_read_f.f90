module mod_read_f
    use mod_parameters
    use mod_variables
    use mod_plib
    implicit none
    
contains

subroutine read_input
    implicit none
    character(200) :: aux
    open(10, file='input')
        read(10, *) aux
        read(10, *) file_radial_profile

        read(10, *) aux
        read(10, *) dir_output

        read(10, *) aux
        read(10, *) body_name

        read(10, *) aux
        read(10, *) SWITCH_one_layer, SWITCH_nondim, SWITCH_bottom_boundary_cond, SWITCH_top_boundary_cond

        read(10, *) aux
        read(10, *) bottom_grav_acc, orb

        read(10, *) aux
        read(10, *) SWITCH_compressible

        read(10, *) aux
        read(10, *) SWITCH_velocity

        read(10, *) aux
        read(10, *) rotation_period, unit_period

        read(10, *) aux
        read(10, *) num_time_steps

        read(10, *) aux
        read(10, *) ecc, obl

        read(10, *) aux
        read(10, *) init_jmax, Max_jmax, mmax

        read(10, *) aux
        read(10, *) SWITCH_inertia

        read(10, *) aux
        read(10, *) SWITCH_coriolis

        read(10, *) aux
        read(10, *) SWITCH_vgradv

        read(10, *) aux
        read(10, *) SWITCH_self_grav

        read(10, *) aux
        read(10, *) final_period

        read(10, *) aux
        read(10, *) SWITCH_ocean_status

        read(10, *) aux
        read(10, *) SWITCH_override_viscosity, ocean_viscosity, SWITCH_hyperviscosity

        read(10, *) aux
        read(10, *) SWITCH_solid_layer

        read(10, *) aux
        read(10, *) SWITCH_hydrostatic_core

        read(10, *) aux
        read(10, *) alpha, zeta

        read(10, *) aux
        read(10, *) checkpoint

        read(10, *) aux
        read(10, *) WRITE_density, WRITE_velocity, WRITE_stress, WRITE_displacement

        read(10, *) aux
        read(10, *) start_output

        read(10, *) aux
        read(10, *) write_output

        read(10, *) aux
        read(10, *) start_simulation

        read(10, *) aux
        read(10, *) SWITCH_mmtides

        read(10, *) aux
        read(10, *) CN_theta
    close(10)
end subroutine read_input


subroutine read_radial_profile
    implicit none
    integer         :: ir, io
    real(dp)        :: radius, density, shear_m, bulk_m, visc
    character(1000) :: header


    ir = 0;
    allocate(temp_radius(0), temp_density(0), temp_shear_m(0), temp_bulk_m(0), temp_visc(0))
    open(10, file=file_radial_profile)
    read(10,'(A)') header
    do
        read(10,*, iostat=io) radius, density, bulk_m, shear_m, visc
        temp_radius   = [temp_radius , radius ]
        temp_density  = [temp_density, density]
        temp_bulk_m   = [temp_bulk_m , bulk_m ]
        temp_shear_m  = [temp_shear_m, shear_m]
        temp_visc     = [temp_visc   , visc   ]
        if(io/=0) exit
        ir = ir + 1
    end do
    close(10)
    N_initial=ir
end subroutine read_radial_profile


subroutine read_initial_cond
    implicit none
    integer :: il, io
    character(300) :: fname
    logical :: fext
    ! This subroutine reads the initial conditions for the simulation.
    ! If the start_simulation flag is not zero, it attempts to read from checkpoint files.
    ! If the checkpoint files do not exist or the specified start_simulation exceeds the last checkpoint,
    ! the subroutine will stop with an error message.
    ! Otherwise, it reads the checkpoint data and initializes the simulation state.
    ! If start_simulation is zero, it sets the initial time step to the value of time_step.
    if(start_simulation /= 0) then
        write(fname, '(A)') trim(dir_output)//'/checkpoint/last_checkpoint_number'
        inquire(file=fname, exist=fext)
        if(.not.fext) stop 'Checkpoint files do not exist...'
        open(10, file=fname, status='old', iostat=io)
        if (io /= 0) then
            stop 'Checkpoint file is empty or does not exist...'
        end if
        read(10, *) start_cp
        close(10)
        call system("rm -f "// trim(dir_output)//"/checkpoint/checkpoint.*")
        print*, '#Reading checkpoint files...'
        if(start_simulation /= -1) then
            if(start_simulation > start_cp) stop 'Exceeded the last checkpoint...'
            write(fname, '(A,I0,A)') trim(dir_output)//"/checkpoint/cp.zip checkpoint.", start_simulation, " -d "// trim(dir_output)//"/checkpoint/"
            call system("unzip  "// trim(fname))
            write(fname, '(A,I0)') trim(dir_output)//"/checkpoint/checkpoint.", start_simulation
            open(10, file=fname, form='unformatted',access='stream',status='old')
            read(10) jmax, t
            do il = 1, numoflayers
                read(10) layer(il)%rhoe
            end do
            do il = 1, numoflayers
                read(10) layer(il)%vel
            end do
            do il = 1, numoflayers
                read(10) layer(il)%pvel
            end do
            do il = 1, numoflayers
                read(10) layer(il)%str
            end do
            do il = 1, numoflayers
                read(10) layer(il)%dis
            end do
            close(10)
        else
            start_simulation = start_cp
            print*, '#starting time step number:', start_cp
            write(fname, '(A,I0,A)') trim(dir_output)//"/checkpoint/cp.zip checkpoint.", start_cp, " -d "// trim(dir_output)//"/checkpoint/"
            call system("unzip  "// trim(fname))
            write(fname, '(A,I0)') trim(dir_output)//"/checkpoint/checkpoint.", start_cp
            open(10, file=fname, form='unformatted',access='stream',status='old')
            read(10) jmax, t
            do il = 1, numoflayers
                read(10) layer(il)%rhoe
            end do
            do il = 1, numoflayers
                read(10) layer(il)%vel
            end do
            do il = 1, numoflayers
                read(10) layer(il)%pvel
            end do
            do il = 1, numoflayers
                read(10) layer(il)%str
            end do
            do il = 1, numoflayers
                read(10) layer(il)%dis
            end do
            close(10)
        end if
    else
        t = 0.0d0
        call system("rm -f "// trim(dir_output)//"/checkpoint/*")
    end if
    call system("rm -f "// trim(dir_output)//"/checkpoint/checkpoint.*")
end subroutine read_initial_cond

!!for going to grid

subroutine read_velocity(file_index, lmax, tt)
    implicit none
    integer                 :: file_index, lmax, il, nn
    character(300)          :: fname
    real(dp)                 :: tt
    intent(in)              :: file_index
    intent(out)             :: lmax, tt

    allocate(layer_grid(numoflayers))
    write(fname, '(A,I0)') trim(dir_output)//"/sh_coeffs/velocity/velocity.", file_index
    open(10, file=fname, form='unformatted',access='stream',status='old')
    read(10) lmax, tt
    do il = 1, numoflayers
        nn = layer(il)%n
        call alloc_complex2(layer_grid(il)%vel, nn, lmax)
        read(10) layer_grid(il)%vel
    end do
    close(10)
    
end subroutine read_velocity

subroutine read_eccentricity
    implicit none
    integer :: io, i
    integer :: numlines
    logical :: fext
    
    if (ecc > 0.0d0)then
        return
    else
        numlines = 0
        inquire(file=trim(dir_output)//'/ecc.dat', exist=fext)
        if (.not. fext) stop 'Error: The file ecc.dat does not exist.'
        do
            read(400, *, iostat=io)
            if(io /= 0) exit
            numlines = numlines + 1
        end do
        close(400)
        rewind(400)
        open(400, file=trim(dir_output)//'/ecc.dat')
        allocate(time_ecc(numlines), ecc_arr(numlines))
        i = 1
        open(400, file=trim(dir_output)//'/ecc.dat')
        do
            read(400, *, iostat=io) time_ecc(i), ecc_arr(i)
            if(io /= 0) exit
            i = i + 1
        end do
        close(400)
    end if
end subroutine read_eccentricity
end module mod_read_f