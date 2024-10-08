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
    REAL(num):: p0, temp_b
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

    rho(-1:nx+2,-1:ny+2,-1:nz+2) = density_init

    temp_b = 1.0_num
    do k = -1, nz+2
       energy(:,:,k) = energy_init*(1.0_num + 0.5_num*(temp_b-1.0_num)*(1.0_num-tanh((zb(k)-1.0_num)/0.1_num)))
    end do

    energy_reference = energy

    energy = energy_init

    grav(-1:nz+2) = 0.0_num

    ALLOCATE(init_pressure(-1:nz+2))

    p0 = energy_init*rho(0,0,nz+2)*(gamma - 1.0_num) !Pressure at upper boundary, then integrate downwards
    init_pressure(nz+2) = p0

    do k = 0, nz+2
        rho(:,:,k) = rho(:,:,k-1) + (dzb(k)/(0.5*(energy(:,:,k)+energy(:,:,k-1))))*&
        (0.25_num*(energy(:,:,k)-energy(:,:,k-1))*(rho(:,:,k) +rho(:,:,k-1))/dzb(k) - &
        0.5_num*(rho(:,:,k) +rho(:,:,k-1))*grav(k)/(gamma-1.0_num))
    end do

    rho(:,:,:) = density_init
    density_init = rho(0,0,1)  !Switch to bottom reference

    do k = 1, nz+2
      vz(:,:,k) = 0.0_num*(zb(k)/zb_global(nz_global))**2
    end do


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

    CALL bfield_bcs


END SUBROUTINE import_bfield

SUBROUTINE try(status)
! Catch error in reading netcdf fild.
INTEGER, INTENT(IN):: status

if (status /= NF90_noerr) THEN
    PRINT*,TRIM(ADJUSTL(NF90_STRERROR(status)))
end if

END SUBROUTINE try


END MODULE initial_conditions
