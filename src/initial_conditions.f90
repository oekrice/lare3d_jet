!  Copyright 2020 University of Warwick

!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at

!      http://www.apache.org/licenses/LICENSE-2.0

!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.

MODULE initial_conditions

  USE shared_data
  USE neutral
  USE boundary
  USE netcdf

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_initial_conditions

CONTAINS

  !****************************************************************************
  ! This function sets up the initial condition for the code
  ! The variables which must be set are:
  !   rho - density
  !   v{x,y,z} - Velocities in x, y, z
  !   b{x,y,z} - Magnetic fields in x, y, z
  !   energy - Specific internal energy
  !
  ! You may also need the neutral fraction. This can be calculated by a
  ! function call to get_neutral(temperature, rho, z). This routine is in
  ! core/neutral.f90 and requires the local temperature and mass density.
  ! For example to set xi_n to the neutral fraction use:
  !   xi_n = get_neutral(temperature, rho, z)
  !****************************************************************************


  SUBROUTINE set_initial_conditions

    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: temperature
    REAL(num) :: xi_v, amp, centre, width
    INTEGER:: i,j,k
    REAL(num):: p0, deltaz, dg
    REAL(num), DIMENSION(:), ALLOCATABLE:: init_pressure !Initial pressure field (1D)

    ! Below are all the variables which must be defined and their sizes
    vx(-2:nx+2, -2:ny+2, -2:nz+2) = 0.0_num
    vy(-2:nx+2, -2:ny+2, -2:nz+2) = 0.0_num
    vz(-2:nx+2, -2:ny+2, -2:nz+2) = 0.0_num

    bx(-2:nx+2, -1:ny+2, -1:nz+2) = 0.0_num
    by(-1:nx+2, -2:ny+2, -1:nz+2) = 0.0_num
    bz(-1:nx+2, -1:ny+2, -2:nz+2) = 0.0_num

    do i = -1,nx+2
    do j = -1,ny+2
    do k = -1,nz+2
    end do
    end do
    end do

    vz(:,:,:) = 0.0_num

    rho(-1:nx+2,-1:ny+2,-1:nz+2) = density_init

    deltaz = zstar/10.0
    do k = -1, nz+2
       energy(:,:,k) = energy_init*(1.0_num + 0.5_num*(chromosphere_temp-1.0_num)*(1.0_num-tanh((zb(k)-zstar)/0.1_num)))
    end do

    energy_reference = energy

    energy = energy_init

    grav(-1:nz+2) = 1e-4

    ALLOCATE(init_pressure(-1:nz+2))

    p0 = energy_init*rho(0,0,nz+2)*(gamma - 1.0_num) !Pressure at upper boundary, then integrate downwards
    init_pressure(nz+2) = p0

    do k = nz+1, -1
    init_pressure(k) = init_pressure(k+1) + 0.5_num*(rho(0,0,k) + rho(0,0,k+1))
    end do

    IF (proc_z_min .ne. MPI_PROC_NULL) THEN
      CALL MPI_RECV(rho(:,:,-1),(nx+4)*(ny+4),mpireal,proc_z_min,tag,comm,status,errcode)
      CALL MPI_RECV(energy(:,:,-1),(nx+4)*(ny+4),mpireal,proc_z_min,tag,comm,status,errcode)
    ENDIF

    !Set up hydrostatic equilibrium (from the flux emergence code).
    !Not entirely sure what this is doing... But it does seem to work.
    DO ix = -1,nx+2
       DO iy = -1,ny+2
          DO iz = 0,nz+2
             dg = 1.0_num/(dzb(iz)+dzb(iz-1))
             rho(ix,iy,iz) = rho(ix,iy,iz-1)*(energy(ix,iy,iz-1)/dzc(iz-1)*(gamma-1.0_num) &
                  - grav(iz-1)*dzb(iz-1)*dg)
             rho(ix,iy,iz) = rho(ix,iy,iz)/(energy(ix,iy,iz)/dzc(iz-1)*(gamma-1.0_num) &
                  + grav(iz-1)*dzb(iz)*dg)
          END DO
       END DO
    END DO

    !Add a constant to the density so the top is the reference
    !rho(:,:,:) = rho(:,:,:) - (rho(nx/2,ny/2,nz) - density_init)
    !print*, maxval(rho), minval(rho)

    IF (proc_z_max .ne. MPI_PROC_NULL) THEN
      CALL MPI_SEND(rho(:,:,nz-1),(nx+4)*(ny+4),mpireal,proc_z_max,tag,comm,errcode)
      CALL MPI_SEND(energy(:,:,nz-1),(nx+4)*(ny+4),mpireal,proc_z_max,tag,comm,errcode)
    ENDIF


    !rho(:,:,:) = density_init
    density_init = rho(0,0,1)  !Switch to bottom reference

    if (rank == 0) then
      print*, 'energy init', energy_init
      print*, 'energy', maxval(energy), minval(energy)
      print*, 'density', maxval(rho), minval(rho)
    end if
    ! set background, non-shock, viscosity
    visc3 = 0.0_num
    IF (use_viscous_damping) THEN
      width = length_z / 10.0_num
      centre = 0.4_num * length_z + width
      amp = 1.e2_num
      DO iz = -1, nz+1
        DO iy = -1, ny + 1
          DO ix = -1, nx + 1
            visc3(ix,iy,iz) = visc3(ix,iy,iz) + amp * (1.0_num + TANH((ABS(zb(iz)) - centre) / width))
          END DO
        END DO
      END DO
    END IF

    ! Import initial magnetic field
    CALL import_bfield

  END SUBROUTINE set_initial_conditions

  SUBROUTINE import_bfield()

    ! Imports the initial condition from the inits file.
    ! Will try to be smart with this, such that the vector potential is only read in in individual processes

    CHARACTER(LEN =64):: init_filename
    INTEGER:: ncid, vid, run_id


    if (run_id < 10) then
      write (init_filename, "(A12, A2, I1, A3)") './inits/init', '00', int(run_id), '.nc'
    else if (run_id < 100) then
      write (init_filename, "(A12, A1, I2, A3)") './inits/init', '0', int(run_id), '.nc'
    else
      write (init_filename, "(A12, I3, A3)") './inits/init', int(run_id), '.nc'
    end if

    if (rank == 0) print*, 'Initial condition filename', trim(init_filename)

    call try(nf90_open(trim(init_filename), nf90_nowrite, ncid))

    call try(nf90_inq_varid(ncid, 'bx', vid))
    call try(nf90_get_var(ncid, vid, bx(0:nx,1:ny,1:nz), &
    start = (/starts(1)+1,starts(2)+1,starts(3)+1/),count = (/nx+1,ny,nz/)))

    call try(nf90_inq_varid(ncid, 'by', vid))
    call try(nf90_get_var(ncid, vid, by(1:nx,0:ny,1:nz), &
    start = (/starts(1)+1,starts(2)+1,starts(3)+1/),count = (/nx,ny+1,nz/)))

    call try(nf90_inq_varid(ncid, 'bz', vid))
    call try(nf90_get_var(ncid, vid, bz(1:nx,1:ny,0:nz), &
    start = (/starts(1)+1,starts(2)+1,starts(3)+1/),count = (/nx,ny,nz+1/)))

    call try(nf90_close(ncid))

    bx = bx*bfield_fact
    by = by*bfield_fact
    bz = bz*bfield_fact

    !print*, 'shape bz', shape(bz), shape(bz_surf_reference)

    CALL bfield_bcs

    allocate(bz_surf_reference(0:nx+1, 0:ny+1))

    bz_surf_reference = bz(0:nx+1, 0:ny+1, 0)

END SUBROUTINE import_bfield

SUBROUTINE try(status)
! Catch error in reading netcdf fild.
INTEGER, INTENT(IN):: status

if (status /= NF90_noerr) THEN
    PRINT*,TRIM(ADJUSTL(NF90_STRERROR(status)))
end if

END SUBROUTINE try


END MODULE initial_conditions
