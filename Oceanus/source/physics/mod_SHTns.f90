module mod_SHTns
    use mod_parameters
    use mod_variables
    use omp_lib, only: omp_get_wtime
    use hdf5
    implicit none
    include './shtns.f03'

    type(shtns_info), pointer :: shtns
    real(dp), pointer :: cosTheta(:), sinTheta(:)
    type(c_ptr) :: shtns_c

    integer  :: lmax, mmax_SHTns, mres, nthreads
    integer  :: nlat, nphi, lm_shtns, l_shtns, m_shtns
    integer  :: nlm, norm, nout, layout
    real(dp) :: eps_polar
    real(dp), allocatable :: curl_j_coef(:)

contains

subroutine SHTns_init_o
    implicit none
    integer :: iph
    lmax        = Max_jmax
    mmax_SHTns  = mmax
    mres     = 1           !! all the orders
    nphi = max(360, 3*mmax_SHTns+1)
    nlat = max(180, 3*lmax/2+5)
    eps_polar = 1.0e-10_dp
    allocate(phi(1:nphi))
    phi = (/(real(iph - 1,dp) * (2.0_dp * pi) / real(Nphi,dp), iph = 1, Nphi)/)
    norm   = SHT_ORTHONORMAL
    layout = SHT_GAUSS + SHT_PHI_CONTIGUOUS
    ! call shtns_verbose(2)
    
    shtns_c = shtns_create(lmax, mmax_SHTns, mres, norm)
    call shtns_set_grid(shtns_c, layout, eps_polar, nlat, nphi)
    
    call c_f_pointer(cptr=shtns_c , fptr=shtns)
    call c_f_pointer(cptr=shtns%ct, fptr=cosTheta, shape=[shtns%nlat])
    call c_f_pointer(cptr=shtns%st, fptr=sinTheta, shape=[shtns%nlat])
    call SHTns_allocation
    allocate(curl_j_coef(1:total_scalar_harm))
    do lm_shtns = 1, total_scalar_harm
        l_shtns = shtns_lm2l(shtns_c, lm_shtns)
        curl_j_coef(lm_shtns) = (real(l_shtns*(l_shtns + 1),dp))
    end do
    print*, 'SHTns has been set up...'

end subroutine SHTns_init_o

subroutine SHTns_allocation
    implicit none
    integer :: il
    do il = 1, numoflayers
        allocate(   layer(il)%Vr(1:shtns%nphi,1:shtns%nlat, 1:layer(il)%n), &
                    layer(il)%Vt(1:shtns%nphi,1:shtns%nlat, 1:layer(il)%n), &
                    layer(il)%Vp(1:shtns%nphi,1:shtns%nlat, 1:layer(il)%n), &
                    layer(il)%VV(1:shtns%nphi,1:shtns%nlat, 1:layer(il)%n))
        allocate(   layer(il)%curlv_Qjm(1:shtns%nlm, 1:layer(il)%n), &
                    layer(il)%curlv_Sjm(1:shtns%nlm, 1:layer(il)%n), &
                    layer(il)%curlv_Tjm(1:shtns%nlm, 1:layer(il)%n))
        allocate(   layer(il)%curlv_Vr(1:shtns%nphi,1:shtns%nlat, 1:layer(il)%n), &
                    layer(il)%curlv_Vt(1:shtns%nphi,1:shtns%nlat, 1:layer(il)%n), &
                    layer(il)%curlv_Vp(1:shtns%nphi,1:shtns%nlat, 1:layer(il)%n))
        allocate(   layer(il)%dSjm(size(layer(il)%curlv_Tjm,1), size(layer(il)%curlv_Tjm,2)))
        allocate(   layer(il)%dTjm(size(layer(il)%curlv_Tjm,1), size(layer(il)%curlv_Tjm,2)))
    end do
end subroutine SHTns_allocation

subroutine convert_to_QST
    implicit none
    integer     :: il, ir, jj, mm, sht_indx, jm1, jp1, jm0, i
    real(dp)    :: j_real, c1Q, c2Q, c1S, c2S, c1T, v_upper, v_lower
    complex(dp) :: avjm1, avjp1, avj

    
    do i = 1, num_oceans
        il = layer_ocean(i)
        
        layer(il)%Qjm = dzero
        layer(il)%Sjm = dzero
        layer(il)%Tjm = dzero
        layer(il)%dSjm = dzero
        layer(il)%dTjm = dzero
        do jj = 1, jmax
            j_real = real(jj,dp)
            c1Q =   j_real          / sqrt(j_real*(2.0 * j_real + 1.0))
            c2Q = -(j_real + 1.0)   / sqrt((j_real+1.0)*(2.0* j_real + 1.0))
            c1S = 1.0 / sqrt(j_real*(2.0 * j_real + 1.0))
            c2S = 1.0 / sqrt((j_real+1.0)*(2.0* j_real + 1.0))
            c1T = 1.0 / sqrt(j_real*(j_real + 1.0))
            do mm = 0, min(jj, mmax)
                jm1 = jml(jj, mm, jj-1)
                jm0 = jml(jj, mm, jj  )
                jp1 = jml(jj, mm, jj+1)
                sht_indx = shtns_lmidx(shtns_c,jj,mm)
                do ir = 1, layer(il)%n
                    v_upper = (layer(il)%r_l(ir+1) - layer(il)%r_i(ir)) / (layer(il)%r_l(ir+1) - layer(il)%r_l(ir))
                    v_lower = (layer(il)%r_i(ir)   - layer(il)%r_l(ir)) / (layer(il)%r_l(ir+1) - layer(il)%r_l(ir))
                    
                    avjm1 = v_upper*layer(il)%vel(jm1, ir)+v_lower*layer(il)%vel(jm1, ir+1)
                    avj   = v_upper*layer(il)%vel(jm0, ir)+v_lower*layer(il)%vel(jm0, ir+1)
                    avjp1 = v_upper*layer(il)%vel(jp1, ir)+v_lower*layer(il)%vel(jp1, ir+1)
                    layer(il)%Qjm(sht_indx,ir) = avjm1*c1Q + avjp1*c2Q

                    layer(il)%Sjm(sht_indx,ir) = avjm1*c1S + avjp1*c2S
                    layer(il)%dSjm(sht_indx, ir) =   c1S * (layer(il)%vel(jm1, ir+1) - layer(il)%vel(jm1, ir)) / (layer(il)%r_l(ir+1) - layer(il)%r_l(ir)) + &
                                                     c2S * (layer(il)%vel(jp1, ir+1) - layer(il)%vel(jp1, ir)) / (layer(il)%r_l(ir+1) - layer(il)%r_l(ir))
                    
                    layer(il)%Tjm(sht_indx,ir)   = im*c1T*avj
                    layer(il)%dTjm(sht_indx, ir) = im * c1T * (layer(il)%vel(jm0, ir+1) - layer(il)%vel(jm0, ir)) / (layer(il)%r_l(ir+1) - layer(il)%r_l(ir))
                end do
            end do
        end do
    end do
end subroutine convert_to_QST


subroutine spectra_to_spatial
    implicit none
    integer :: ir, il, l

    !$OMP PARALLEL
    do il = 1, num_oceans
        !$OMP DO COLLAPSE(2)
        do ir = 1, layer(il)%n
            do l = 1, total_scalar_harm
                layer(il)%curlv_Qjm(l, ir) =   curl_j_coef(l) * layer(il)%Tjm(l, ir) / layer(il)%r_i(ir)
                layer(il)%curlv_Sjm(l, ir) =   (layer(il)%dTjm(l, ir) + layer(il)%Tjm(l, ir) / layer(il)%r_i(ir))
                layer(il)%curlv_Tjm(l, ir) = - layer(il)%Qjm(l, ir) / layer(il)%r_i(ir) + layer(il)%dSjm(l, ir) + layer(il)%Sjm(l, ir) / layer(il)%r_i(ir)
            end do
        end do
        !$OMP END DO
        !$OMP DO
        do ir = 1, layer(il)%n
            call SHqst_to_spat(shtns_c, layer(il)%Qjm(:,ir), layer(il)%Sjm(:,ir), layer(il)%Tjm(:,ir), layer(il)%Vr(:,:,ir), layer(il)%Vt(:,:, ir), layer(il)%Vp(:,:,ir))
            call SHqst_to_spat(shtns_c, layer(il)%curlv_Qjm(:,ir), layer(il)%curlv_Sjm(:,ir), layer(il)%curlv_Tjm(:,ir), layer(il)%curlv_Vr(:,:,ir), layer(il)%curlv_Vt(:,:, ir), layer(il)%curlv_Vp(:,:,ir))
        end do
        !$OMP END DO
    end do
    !$OMP END PARALLEL

end subroutine spectra_to_spatial

subroutine spatial_to_spectra
    implicit none
    integer :: ir, il, ith, iph
    real(dp) :: hup, hlo, denom
    real(dp), allocatable :: cp_curlVr(:,:), cp_curlVt(:,:), cp_curlVp(:,:)
    
    il = 1
    do il = 1, num_oceans
        !$OMP PARALLEL PRIVATE(cp_curlVr, cp_curlVt, cp_curlVp)
        allocate(cp_curlVr(size(layer(il)%curlv_Vr,1), size(layer(il)%curlv_Vr,2)))
        allocate(cp_curlVt(size(layer(il)%curlv_Vt,1), size(layer(il)%curlv_Vt,2)))
        allocate(cp_curlVp(size(layer(il)%curlv_Vp,1), size(layer(il)%curlv_Vp,2)))
        !$OMP DO
        do ir = 1, layer(il)%n
            do ith = 1, Nlat
                do iph = 1, nphi
                    layer(il)%VV(iph, ith, ir) = layer(il)%Vr(iph, ith, ir)*layer(il)%Vr(iph, ith, ir) + &
                                                 layer(il)%Vt(iph, ith, ir)*layer(il)%Vt(iph, ith, ir) + &
                                                 layer(il)%Vp(iph, ith, ir)*layer(il)%Vp(iph, ith, ir)
                    cp_curlVr(iph, ith) = layer(il)%curlv_Vr(iph, ith, ir)
                    cp_curlVt(iph, ith) = layer(il)%curlv_Vt(iph, ith, ir)
                    cp_curlVp(iph, ith) = layer(il)%curlv_Vp(iph, ith, ir)
                    ! curl(v) x v:
                    layer(il)%curlv_Vr(iph, ith, ir) =  (layer(il)%Vp(iph, ith, ir) * cp_curlVt(iph, ith) &
                                                      -  layer(il)%Vt(iph, ith, ir) * cp_curlVp(iph, ith))
                    layer(il)%curlv_Vt(iph, ith, ir) = -(layer(il)%Vp(iph, ith, ir) * cp_curlVr(iph, ith) &
                                                       - layer(il)%Vr(iph, ith, ir) * cp_curlVp(iph, ith))
                    layer(il)%curlv_Vp(iph, ith, ir) =  (layer(il)%Vt(iph, ith, ir) * cp_curlVr(iph, ith) &
                                                      -  layer(il)%Vr(iph, ith, ir) * cp_curlVt(iph, ith))
                end do
            end do
            call spat_to_SH(shtns_c, layer(il)%VV(:,:,ir), layer(il)%Sjm(:,ir))
            call spat_to_SHqst(shtns_c, layer(il)%curlv_Vr(:,:,ir), layer(il)%curlv_Vt(:,:, ir), layer(il)%curlv_Vp(:,:,ir), layer(il)%curlv_Qjm(:,ir), layer(il)%curlv_Sjm(:,ir), layer(il)%curlv_Tjm(:,ir))
        end do
        !$OMP END DO
        deallocate(cp_curlVr, cp_curlVt, cp_curlVp)
        !$OMP END PARALLEL
    end do

    do il = 1, num_oceans
        do ir = 1, layer(il)%n-1
            denom   = layer(il)%r_i(ir+1) - layer(il)%r_i(ir)
            hup     = (layer(il)%r_i(ir+1)   - layer(il)%r_l(ir+1)) / denom
            hlo     = (layer(il)%r_l(ir+1)   - layer(il)%r_i(ir))   / denom

            layer(il)%curlv_Qjm(:, ir) = hup*layer(il)%curlv_Qjm(:, ir) + hlo*layer(il)%curlv_Qjm(:, ir+1)
            layer(il)%curlv_Sjm(:, ir) = hup*layer(il)%curlv_Sjm(:, ir) + hlo*layer(il)%curlv_Sjm(:, ir+1)
            layer(il)%curlv_Tjm(:, ir) = hup*layer(il)%curlv_Tjm(:, ir) + hlo*layer(il)%curlv_Tjm(:, ir+1)
        end do
    end do
end subroutine spatial_to_spectra

subroutine vgradv
    implicit none
    integer  :: j, m, il, ir, idx_jm1, idx_jp1, idx_jm0
    real(dp) :: jj, d1, d2, hup, hlo, denom, r_lip1, ccjm1, ccjp1, ccjm0
    complex(dp) :: grad_minus, grad_plus, avSjm, Sjm_ir, Sjm_irp1
    ! real(dp) :: time_b, time_e
    ! time_b  = omp_get_wtime()
    call convert_to_QST
    ! time_e = omp_get_wtime(); print*, 'Time for convert_to_QST in vgradv: ', time_e - time_b
    ! time_b  = omp_get_wtime()
    call spectra_to_spatial
    ! time_e = omp_get_wtime(); print*, 'Time for spectra_to_spatial in vgradv: ', time_e - time_b
    ! time_b  = omp_get_wtime()
    call spatial_to_spectra
    ! time_e = omp_get_wtime(); print*, 'Time for spatial_to_spectra in vgradv: ', time_e - time_b
    ! time_b  = omp_get_wtime()
    do il = 1, num_oceans
        layer(il)%vgradv = dzero
            do j = 1, jmax
                jj = real(j, dp)
                d1 = sqrt(jj / (2.0*jj + 1.0))
                d2 = sqrt((jj + 1.0) / (2.0*jj + 1.0))
                ccjm1 =  sqrt(jj*(2.0*jj + 1.0)) / (2.0*jj + 1.0)
                ccjp1 = -sqrt((jj + 1.0)*(2.0*jj + 1.0)) / (2.0*jj + 1.0)
                ccjm0 =  sqrt(jj*(jj + 1.0))

                do m = 0, min(j, mmax)
                    lm_shtns = shtns_lmidx(shtns_c, j, m)
                    do ir = 1, layer(il)%n-1
                        denom   = layer(il)%r_i(ir+1)  - layer(il)%r_i(ir)
                        hup     = (layer(il)%r_i(ir+1) - layer(il)%r_l(ir)) / denom
                        hlo     = (layer(il)%r_l(ir)   - layer(il)%r_i(ir)) / denom
                        r_lip1  = layer(il)%r_l(ir+1)

                        ! convert denom and r_lip1 to reciprocals once to avoid repeated divisions
                        denom = 1.0_dp / denom
                        r_lip1 = 1.0_dp / r_lip1
                        ! load spectral values once
                        Sjm_ir   = layer(il)%Sjm(lm_shtns, ir)
                        Sjm_irp1 = layer(il)%Sjm(lm_shtns, ir+1)

                        ! averaged and radial derivative (use reciprocals)
                        avSjm = hup*Sjm_ir + hlo*Sjm_irp1
                        grad_minus = d1 * (Sjm_irp1 - Sjm_ir) * denom
                        grad_minus = grad_minus + d1 * (jj + 1.0_dp) * avSjm * r_lip1
                        idx_jm1 = jml(j, m, j-1)
                        layer(il)%vgradv(ir, idx_jm1) = -0.5_dp * grad_minus

                        grad_plus = -d2 * (Sjm_irp1 - Sjm_ir) * denom
                        grad_plus = grad_plus + d2 * jj * avSjm * r_lip1
                        idx_jp1 = jml(j, m, j+1)
                        layer(il)%vgradv(ir, idx_jp1) = -0.5_dp * grad_plus

                        ! add curl contributions (fetch curl components once)
                        Sjm_ir   = layer(il)%curlv_Qjm(lm_shtns, ir)
                        Sjm_irp1 = layer(il)%curlv_Sjm(lm_shtns, ir)
                        layer(il)%vgradv(ir, idx_jm1) = layer(il)%vgradv(ir, idx_jm1) + ccjm1 * (Sjm_ir + (jj+1.0_dp)*Sjm_irp1)
                        layer(il)%vgradv(ir, idx_jp1) = layer(il)%vgradv(ir, idx_jp1) + ccjp1 * (Sjm_ir - (jj       )*Sjm_irp1)
                        Sjm_ir = layer(il)%curlv_Tjm(lm_shtns, ir)
                        idx_jm0 = jml(j, m, j)
                        layer(il)%vgradv(ir, idx_jm0) = layer(il)%vgradv(ir, idx_jm0) - im * ccjm0 * Sjm_ir
                end do
            end do
        end do
    end do
    ! time_e = omp_get_wtime(); print*, 'Time for vgradv computation: ', time_e - time_b
    ! print*, '================================'
end subroutine vgradv


subroutine write_velocity_grid_h5(it2)
    implicit none
    character(len=300) :: filename
    integer            :: it2, il, ir
    integer            :: error
    integer(HID_T)     :: ID_file, ID_sp, ID_dset, ID_G_coord
    integer(HSIZE_T)   :: dims(3), dims2(1)
    
    il = layer_ocean(1)
    call convert_to_QST
    do ir = 1, layer(il)%n
        call SHqst_to_spat(shtns_c, layer(il)%Qjm(:,ir), &
                                    layer(il)%Sjm(:,ir), &
                                    layer(il)%Tjm(:,ir), &
                                    layer(il)%Vr(:,:,ir), &
                                    layer(il)%Vt(:,:,ir), &
                                    layer(il)%Vp(:,:,ir))
    end do
    dims = shape(layer(il)%Vr)
    write(filename, '(A,I0,A)') trim(dir_output)//'/on_grid/velocity/velocity_', it2, '.h5'
    call h5open_f(error)
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, ID_file, error)

    ! Create dataspace for 3D arrays
    call h5screate_simple_f(3, dims, ID_sp, error)

    ! Create and write each component
    call h5dcreate_f(ID_file, "vr", H5T_NATIVE_DOUBLE, ID_sp, ID_dset, error)
    call h5dwrite_f(ID_dset, H5T_NATIVE_DOUBLE, layer(il)%Vr, dims, error)
    call h5dclose_f(ID_dset, error)

    call h5dcreate_f(ID_file, "vtheta", H5T_NATIVE_DOUBLE, ID_sp, ID_dset, error)
    call h5dwrite_f(ID_dset, H5T_NATIVE_DOUBLE, layer(il)%Vt, dims, error)
    call h5dclose_f(ID_dset, error)

    call h5dcreate_f(ID_file, "vphi", H5T_NATIVE_DOUBLE, ID_sp, ID_dset, error)
    call h5dwrite_f(ID_dset, H5T_NATIVE_DOUBLE, layer(il)%Vt, dims, error)
    call h5dclose_f(ID_dset, error)
    call h5sclose_f(ID_sp, error)

    dims2 = size(layer(il)%r_i(1:layer(il)%n))
    call h5screate_simple_f(1, dims2, ID_sp, error)
    call h5dcreate_f(ID_file, "r", H5T_NATIVE_DOUBLE, ID_sp, ID_dset, error)
    call h5dwrite_f(ID_dset, H5T_NATIVE_DOUBLE, layer(il)%r_l(1:layer(il)%n), dims2, error)
    call h5dclose_f(ID_dset, error)
    call h5sclose_f(ID_sp, error)

    dims2 = size(cosTheta)
    call h5screate_simple_f(1, dims2, ID_sp, error)
    call h5dcreate_f(ID_file, "theta", H5T_NATIVE_DOUBLE, ID_sp, ID_dset, error)
    call h5dwrite_f(ID_dset, H5T_NATIVE_DOUBLE, acos(cosTheta), dims2, error)
    call h5dclose_f(ID_dset, error)
    call h5sclose_f(ID_sp, error)

    dims2 = size(phi)
    call h5screate_simple_f(1, dims2, ID_sp, error)
    call h5dcreate_f(ID_file, "phi", H5T_NATIVE_DOUBLE, ID_sp, ID_dset, error)
    call h5dwrite_f(ID_dset, H5T_NATIVE_DOUBLE, phi, dims2, error)
    call h5dclose_f(ID_dset, error)
    call h5sclose_f(ID_sp, error)

    call h5gcreate_f(ID_file, "Cartesian_coordinates", ID_G_coord, error)

    call h5gclose_f(ID_G_coord, error)
    
    call h5fclose_f(ID_file, error)
    call h5close_f(error)
end subroutine write_velocity_grid_h5

subroutine write_xdmf_file(it2)
    implicit none
    integer :: nr, nth, nph
    integer :: unit, it2, il
    character(len=300)  :: filename

    il  = layer_ocean(1)
    nr  = size(layer(il)%r_i)
    nth = size(cosTheta)
    nph = size(phi)
    write(filename, '(A,I0,A)') trim(dir_output)//'/on_grid/velocity/velocity_', it2, '.xmf'
    open(newunit=unit, file=filename, status='replace', action='write')

    write(unit,'(A)') '<?xml version="1.0" ?>'
    write(unit,'(A)') '<Xdmf Version="3.0">'
    write(unit,'(A)') '  <Domain>'
    write(unit,'(A)') '    <Grid Name="VelocityField" GridType="Uniform">'
    write(unit,'(A)') '      <Topology TopologyType="3DRectMesh" Dimensions="'// &
                      trim(adjustl(itoa(nr)))//' '//trim(adjustl(itoa(nth)))//' '//trim(adjustl(itoa(nph)))//'"/>'

    write(unit,'(A)') '      <Geometry GeometryType="XYZ">'
    write(unit,'(A)') '        <DataItem Dimensions="'//trim(adjustl(itoa(nr)))//' '// &
                      trim(adjustl(itoa(nth)))//' '//trim(adjustl(itoa(nph)))//' 3" Format="HDF">'
    write(unit,'(A)') '          velocity_field.h5:/Coordinates/cartesian'
    write(unit,'(A)') '        </DataItem>'
    write(unit,'(A)') '      </Geometry>'

    write(unit,'(A)') '      <Attribute Name="Velocity" AttributeType="Vector" Center="Node">'
    write(unit,'(A)') '        <DataItem ItemType="Function" Function="JOIN($0, $1, $2)" Dimensions="'// &
                      trim(adjustl(itoa(nr)))//' '//trim(adjustl(itoa(nth)))//' '//trim(adjustl(itoa(nph)))//' 3">'
    write(unit,'(A)') '          <DataItem Dimensions="'//trim(adjustl(itoa(nr)))//' '//trim(adjustl(itoa(nth)))//' '//trim(adjustl(itoa(nph)))//'" Format="HDF">'
    write(unit,'(A)') '            velocity_field.h5:/vr'
    write(unit,'(A)') '          </DataItem>'
    write(unit,'(A)') '          <DataItem Dimensions="'//trim(adjustl(itoa(nr)))//' '//trim(adjustl(itoa(nth)))//' '//trim(adjustl(itoa(nph)))//'" Format="HDF">'
    write(unit,'(A)') '            velocity_field.h5:/vtheta'
    write(unit,'(A)') '          </DataItem>'
    write(unit,'(A)') '          <DataItem Dimensions="'//trim(adjustl(itoa(nr)))//' '//trim(adjustl(itoa(nth)))//' '//trim(adjustl(itoa(nph)))//'" Format="HDF">'
    write(unit,'(A)') '            velocity_field.h5:/vphi'
    write(unit,'(A)') '          </DataItem>'
    write(unit,'(A)') '        </DataItem>'
    write(unit,'(A)') '      </Attribute>'

    write(unit,'(A)') '    </Grid>'
    write(unit,'(A)') '  </Domain>'
    write(unit,'(A)') '</Xdmf>'

    close(unit)
    contains
    function itoa(i) result(str)
        implicit none
        integer, intent(in) :: i
        character(len=20) :: str
        write(str,'(I0)') i
    end function itoa

end subroutine write_xdmf_file



! subroutine analytical_vgradv
!     implicit none
!     integer :: ir, ith, iph, l, j, m, lm_shtns
!     real(dp) :: r, ph, ccjm1, ccjp1, ccjm0, jj
!     real(dp), allocatable :: VgVr(:,:,:), VgVt(:,:,:), VgVp(:,:,:)
!     real(dp), allocatable :: VVr(:,:,:), VVt(:,:,:), VVp(:,:,:)
!     complex(dp), allocatable :: QQjm(:,:), QSjm(:,:), QTjm(:,:)
!     complex(dp), allocatable :: gQjm(:,:), gSjm(:,:), gTjm(:,:)
    

!     allocate(VgVr(1:shtns%nphi,1:shtns%nlat, 1:layer(1)%n+1))
!     allocate(VgVt(1:shtns%nphi,1:shtns%nlat, 1:layer(1)%n+1))
!     allocate(VgVp(1:shtns%nphi,1:shtns%nlat, 1:layer(1)%n+1))
!     allocate(VVr(1:shtns%nphi,1:shtns%nlat, 1:layer(1)%n+1))
!     allocate(VVt(1:shtns%nphi,1:shtns%nlat, 1:layer(1)%n+1))
!     allocate(VVp(1:shtns%nphi,1:shtns%nlat, 1:layer(1)%n+1))
!     allocate(QQjm(size(layer(1)%Qjm,1), size(layer(1)%Qjm,2)))
!     allocate(QSjm(size(layer(1)%Sjm,1), size(layer(1)%Sjm,2)))
!     allocate(QTjm(size(layer(1)%Tjm,1), size(layer(1)%Tjm,2)))
!     allocate(gQjm(size(layer(1)%Qjm,1), layer(1)%n-1))
!     allocate(gSjm(size(layer(1)%Sjm,1), layer(1)%n-1))
!     allocate(gTjm(size(layer(1)%Tjm,1), layer(1)%n-1))
!     QQjm = dzero
!     QSjm = dzero
!     QTjm = dzero
!     gQjm = dzero
!     gSjm = dzero
!     gTjm = dzero
!     VVr = zero
!     VVt = zero
!     VVp = zero
!     VgVr = zero
!     VgVt = zero
!     VgVp = zero
!     print*, rad(1)
!     do ir = 1, layer(1)%n+1
!         r = layer(1)%r_l(ir)
!         do ith = 1, Nlat
!             do iph = 1, Nphi
!                 ph = real(iph - 1,dp) * (2.0_dp * pi) / real(Nphi,dp)
!                 VVr(iph, ith, ir) = zero
!                 VVt(iph, ith, ir) = -2*sin(ph)*1e4
!                 VVp(iph, ith, ir) =-1e4*(-2)*cos(ph)*cosTheta(ith)
!             end do
!         end do
!     end do

!     do ir = 1, layer(1)%n-1
!         r = layer(1)%r_l(ir+1)
!         do ith = 1, Nlat
!             do iph = 1, Nphi
!                 ph = real(iph - 1,dp) * (2.0_dp * pi) / real(Nphi,dp)
!                 VgVr(iph, ith, ir) = 1e8*(-(4*sin(ph)**2)/r - (4*(cos(ph)**2)*cosTheta(ith)**2)/r)
!                 VgVt(iph, ith, ir) = 1e8*(2*cos(ph)*cosTheta(ith)*((2*cos(ph))/(r*sinTheta(ith)) - (2*cos(ph)*cosTheta(ith)**2)/(r*sinTheta(ith))))
!                 VgVp(iph, ith, ir) = 1e8*(-(4*cos(ph)*sin(ph)*sinTheta(ith))/r)

!             end do
!         end do
!     end do
!     do ir = 1, layer(1)%n+1
!         ! call spat_to_SHqst(shtns_c, VVr(:,:,ir), VVt(:,:, ir), VVp(:,:,ir), QQjm(:,ir), QSjm(:,ir), QTjm(:,ir))
!         call spat_to_SH(shtns_c, VVt(:,:,ir), QSjm(:,ir))
!         ! if(ir<layer(1)%n) call spat_to_SHqst(shtns_c, VgVr(:,:,ir), VgVt(:,:, ir), VgVp(:,:,ir), gQjm(:,ir), gSjm(:,ir), gTjm(:,ir))
!     end do
!     do l = 1, total_scalar_harm
!         j = shtns_lm2l(shtns_c, l)
!         m = shtns_lm2m(shtns_c, l)
!         print*,l, j, m, QSjm(l,9)!, QTjm(l,9), &
!                 ! gQjm(l,9), gSjm(l,9), gTjm(l,9)
!     end do

!     do ir = 1, layer(1)%n+1
!         layer(1)%vel(ir, :) = dzero
!         do j = 1, jmax
!             jj = real(j, dp)
!             ccjm1 =  sqrt(jj*(2.0*jj + 1.0)) / (2.0*jj + 1.0)
!             ccjp1 = -sqrt((jj + 1.0)*(2.0*jj + 1.0)) / (2.0*jj + 1.0)
!             ccjm0 =  sqrt(jj*(jj + 1.0))
!             do m = 0, j
!                 lm_shtns = shtns_lmidx(shtns_c, j, m)
!                 layer(1)%vel(ir, jml(j,m,j-1)) = ccjm1 * (QQjm(lm_shtns, ir) + (jj+1.0)*QSjm(lm_shtns, ir))
!                 layer(1)%vel(ir, jml(j,m,j+1)) = ccjp1 * (QQjm(lm_shtns, ir) + (   -jj)*QSjm(lm_shtns, ir))
!                 layer(1)%vel(ir, jml(j,m,j  )) = -im * ccjm0 * QTjm(lm_shtns, ir)
!             end do
!         end do
!     end do
! end subroutine analytical_vgradv

end module mod_SHTns