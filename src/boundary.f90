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

!******************************************************************************
! This module contains the boundary conditions for the entire code
! Any new boundary conditions should be added here
!******************************************************************************

MODULE boundary

  USE shared_data
  USE mpiboundary

  IMPLICIT NONE

CONTAINS

  !****************************************************************************
  ! Set up any necessary variables for the chosen boundary conditions
  !****************************************************************************

  SUBROUTINE set_boundary_conditions

    any_open = .FALSE.
    IF (xbc_min == BC_OPEN .OR. xbc_max == BC_OPEN &
        .OR. ybc_min == BC_OPEN .OR. ybc_max == BC_OPEN &
        .OR. zbc_min == BC_OPEN .OR. zbc_max == BC_OPEN) any_open = .TRUE.

  END SUBROUTINE set_boundary_conditions



  !****************************************************************************
  ! Call all of the boundaries needed by the core Lagrangian solver
  !****************************************************************************

  SUBROUTINE boundary_conditions

    CALL bfield_bcs
    CALL energy_bcs
    CALL density_bcs
    CALL velocity_bcs

  END SUBROUTINE boundary_conditions



  !****************************************************************************
  ! Boundary conditions for magnetic field through plane
  !****************************************************************************

  SUBROUTINE bfield_bcs

    CALL bfield_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      bz(0,0:ny+1,0:nz) = bz(1,0:ny+1,0:nz) - dxb(0)*(bx(0,0:ny+1,1:nz+1) - bx(0, 0:ny+1,0:nz))/dzb(0)
      by(0,0:ny,0:nz+1) = by(1, 0:ny,0:nz+1) - dxb(0)*(bx(0,1:ny+1,0:nz+1) - bx(0,0:ny,0:nz+1))/dyb(0)
      bx(-1,:,:) = 0.0_num
      bx(-2,:,:) = 0.0_num
      by(-1,:,:) = by(2,:,:)
      bz(-1,:,:) = bz(2,:,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      bz(nx+1,0:ny+1,0:nz) = bz(nx,0:ny+1,0:nz) + dxb(0)*(bx(nx,0:ny+1,1:nz+1) - bx(nx, 0:ny+1,0:nz))/dzb(0)
      by(nx+1,0:ny,0:nz+1) = by(nx, 0:ny,0:nz+1) + dxb(0)*(bx(nx,1:ny+1,0:nz+1) - bx(nx,0:ny,0:nz+1))/dyb(0)
      bx(nx+1,:,:) = 0.0_num
      bx(nx+2,:,:) = 0.0_num
      by(nx+2,:,:) = by(nx-1,:,:)
      bz(nx+2,:,:) = bz(nx-1,:,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      bz(0:nx+1,0,0:nz) = bz(0:nx+1, 1,0:nz) - dyb(0)*(by(0:nx+1,0,1:nz+1) - by(0:nx+1,0,0:nz))/dzb(0)
      bx(0:nx,0,0:nz+1) = bx(0:nx,1,0:nz+1) - dyb(0)*(by(1:nx+1,0,0:nz+1) - by(0:nx, 0,0:nz+1))/dxb(0)
      bx(:,-1,:) = bx(:,2,:)
      by(:,-1,:) = 0.0_num
      by(:,-2,:) = 0.0_num
      bz(:,-1,:) = bz(:,2,:)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN

      bz(0:nx+1,ny+1,0:nz) = bz(0:nx+1, ny,0:nz) + dyb(0)*(by(0:nx+1,ny,1:nz+1) - by(0:nx+1,ny,0:nz))/dzb(0)
      bx(0:nx,ny+1,0:nz+1) = bx(0:nx,ny,0:nz+1) + dyb(0)*(by(1:nx+1,ny,0:nz+1) - by(0:nx, ny,0:nz+1))/dxb(0)

      bx(:,ny+2,:) = bx(:,ny-1,:)
      by(:,ny+1,:) = 0.0_num
      by(:,ny+2,:) = 0.0_num
      bz(:,ny+2,:) = bz(:,ny-1,:)
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
      by(0:nx+1,0:ny,0) = by(0:nx+1,0:ny,1) - dzb(0)*(bz(0:nx+1,1:ny+1,0) - bz(0:nx+1, 0:ny,0))/dyb(0)
      bx(0:nx, 0:ny+1,0) = bx(0:nx,0:ny+1,1) - dzb(0)*(bz(1:nx+1,0:ny+1,0) - bz(0:nx,0:ny+1,0))/dxb(0)
      bx(:,:,-1) = bx(:,:,2)
      by(:,:,-1) = by(:,:,2)
      bz(:,:,-1) = bz(:,:,1)
      bz(:,:,-2) = bz(:,:,2)
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
      by(0:nx+1,0:ny,nz+1) = by(0:nx+1,0:ny,nz) + dzb(0)*(bz(0:nx+1,1:ny+1,nz) - bz(0:nx+1, 0:ny,nz))/dyb(0)
      bx(0:nx, 0:ny+1,nz+1) = bx(0:nx,0:ny+1,nz) + dzb(0)*(bz(1:nx+1,0:ny+1,nz) - bz(0:nx,0:ny+1,nz))/dxb(0)
      bx(:,:,nz+2) = bx(:,:,nz-1)
      by(:,:,nz+2) = by(:,:,nz-1)
      bz(:,:,nz+1) = bz(:,:,nz-1)
      bz(:,:,nz+2) = bz(:,:,nz-2)
    END IF

  END SUBROUTINE bfield_bcs



  !****************************************************************************
  ! Boundary conditions for specific internal energy
  !****************************************************************************

  SUBROUTINE energy_bcs

    CALL energy_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      energy( 0,:,:) = energy(1,:,:)
      energy(-1,:,:) = energy(2,:,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      energy(nx+1,:,:) = energy(nx  ,:,:)
      energy(nx+2,:,:) = energy(nx-1,:,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      energy(:, 0,:) = energy(:,1,:)
      energy(:,-1,:) = energy(:,2,:)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      energy(:,ny+1,:) = energy(:,ny  ,:)
      energy(:,ny+2,:) = energy(:,ny-1,:)
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
      energy(:,:, 0) = energy(:,:,1)
      energy(:,:,-1) = energy(:,:,2)
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
      energy(:,:,nz+1) = energy(:,:,nz)
      energy(:,:,nz+2) = energy(:,:,nz-1)
    END IF

  END SUBROUTINE energy_bcs



  !****************************************************************************
  ! Boundary conditions for density
  !****************************************************************************

  SUBROUTINE density_bcs

    CALL density_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      rho( 0,:,:) = rho(1,:,:)
      rho(-1,:,:) = rho(2,:,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      rho(nx+1,:,:) = rho(nx  ,:,:)
      rho(nx+2,:,:) = rho(nx-1,:,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      rho(:, 0,:) = rho(:,1,:)
      rho(:,-1,:) = rho(:,2,:)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      rho(:,ny+1,:) = rho(:,ny  ,:)
      rho(:,ny+2,:) = rho(:,ny-1,:)
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
      rho(:,:, 1) = density_init
      rho(:,:, 0) = rho(:,:,1)
      rho(:,:,-1) = rho(:,:,2)
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
      rho(:,:,nz+1) = rho(:,:,nz)
      rho(:,:,nz+2) = rho(:,:,nz-1)
    END IF

  END SUBROUTINE density_bcs


  !****************************************************************************
  ! Boundary conditions for temperature
  !****************************************************************************

  SUBROUTINE temperature_bcs

    CALL density_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      temperature( 0,:,:) = temperature(1,:,:)
      temperature(-1,:,:) = temperature(2,:,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      temperature(nx+1,:,:) = temperature(nx  ,:,:)
      temperature(nx+2,:,:) = temperature(nx-1,:,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      temperature(:, 0,:) = temperature(:,1,:)
      temperature(:,-1,:) = temperature(:,2,:)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      temperature(:,ny+1,:) = temperature(:,ny  ,:)
      temperature(:,ny+2,:) = temperature(:,ny-1,:)
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
      temperature(:,:, 0) = temperature(:,:,1)
      temperature(:,:,-1) = temperature(:,:,2)
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
      temperature(:,:,nz+1) = temperature(:,:,nz  )
      temperature(:,:,nz+2) = temperature(:,:,nz-1)
    END IF

  END SUBROUTINE temperature_bcs




  !****************************************************************************
  ! Full timestep velocity boundary conditions
  !****************************************************************************

  SUBROUTINE velocity_bcs

    real(num):: br, bl, kb    !Parameters for the boundary motions (taken from Pariat)
    real(num), dimension(:,:):: bzdy(0:nx+1,0:ny), bzdx(0:nx,0:ny+1)
    real(num), dimension(:,:):: bzdy0(0:nx,0:ny), bzdx0(0:nx,0:ny)
    real(num), dimension(:,:):: fact(0:nx+1,0:ny+1), fact0(0:nx,0:ny)
    real(num), dimension(:,:):: vx_surf(0:nx,0:ny), vy_surf(0:nx,0:ny)
    INTEGER:: i,j

    CALL velocity_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      vx(-2:0,:,:) = 0.0_num
      vy(-2:0,:,:) = 0.0_num
      vz(0,:,:) = vz(1,:,:)
      vz(-1,:,:) = vz(0,:,:)
      vz(-2,:,:) = vz(-1,:,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      vx(nx:nx+2,:,:) = 0.0_num
      vy(nx:nx+2,:,:) = 0.0_num
      vz(nx,:,:) = vz(nx-1,:,:)
      vz(nx+1,:,:) = vz(nx,:,:)
      vz(nx+2,:,:) = vz(nx+2,:,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      vx(:,-2:0,:) = 0.0_num
      vy(:,-2:0,:) = 0.0_num
      vz(:,0,:) = vz(:,1,:)
      vz(:,-1,:) = vz(:,0,:)
      vz(:,-2,:) = vz(:,-1,:)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      vx(:,ny:ny+2,:) = 0.0_num
      vy(:,ny:ny+2,:) = 0.0_num
      vz(:,ny,:) = vz(:,ny-1,:)
      vz(:,ny+1,:) = vz(:,ny,:)
      vz(:,ny+2,:) = vz(:,ny+1,:)
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
      vx(:,:,-2:0) = 0.0_num
      vy(:,:,-2:0) = 0.0_num
      vz(:,:,-2:0) = 0.0_num
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
      vx(:,:,nz:nz+2) = 0.0_num
      vy(:,:,nz:nz+2) = 0.0_num

      do i = -2, nx+2
          do j = -2, ny+2
          vz(i, j, nz) = max(0.0,1.25*vz(i,j,nz-1))
          end do
      end do
      vz(:,:,nz+1) = vz(:,:,nz)
      vz(:,:,nz+1) = vz(:,:,nz+1)

    END IF

    if (proc_z_min == MPI_PROC_NULL) then
      br = 13.0_num*bfield_fact; bl = 0.1_num*bfield_fact; kb = 15.0_num

      bzdy = (bz(0:nx+1,1:ny+1,0) - bz(0:nx+1,0:ny,0)) / dyb(1)
      bzdx = (bz(1:nx+1,0:ny+1,0) - bz(0:nx,0:ny+1,0)) / dxb(1)

      bzdy0 = 0.5_num*(bzdy(1:nx+1,0:ny) + bzdy(0:nx,0:ny))
      bzdx0 = 0.5_num*(bzdx(0:nx,1:ny+1) + bzdx(0:nx,0:ny))

      fact = (kb*(br-bl))/(bz(0:nx+1,0:ny+1,0) + 1d-10)*tanh(kb*(bz(0:nx+1,0:ny+1,0)- bl)/(br-bl+1d-10))
      fact0 = 0.25_num*(fact(0:nx,0:ny) + fact(1:nx+1,0:ny) + fact(0:nx,1:ny+1) + fact(1:nx+1, 1:ny+1))

      vx_surf(0:nx, 0:ny) = -fact0*bzdy0
      vy_surf(0:nx, 0:ny) = fact0*bzdx0

      vx(0:nx,0:ny,0) = vx(0:nx,0:ny,0) + shearing_fact*vx_surf
      vy(0:nx,0:ny,0) = vy(0:nx,0:ny,0) + shearing_fact*vy_surf

    end if

  END SUBROUTINE velocity_bcs



  !****************************************************************************
  ! Half timestep velocity boundary conditions
  !****************************************************************************

  SUBROUTINE remap_v_bcs

    real(num):: br, bl, kb    !Parameters for the boundary motions (taken from Pariat)
    real(num), dimension(:,:):: bzdy(0:nx+1,0:ny), bzdx(0:nx,0:ny+1)
    real(num), dimension(:,:):: bzdy0(0:nx,0:ny), bzdx0(0:nx,0:ny)
    real(num), dimension(:,:):: fact(0:nx+1,0:ny+1), fact0(0:nx,0:ny)
    real(num), dimension(:,:):: vx_surf(0:nx,0:ny), vy_surf(0:nx,0:ny)
    INTEGER:: i,j

    CALL remap_v_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      vx1(-2:0,:,:) = 0.0_num
      vy1(-2:0,:,:) = 0.0_num
      vz1(0,:,:) = vz1(1,:,:)
      vz1(-1,:,:) = vz1(0,:,:)
      vz1(-2,:,:) = vz1(-1,:,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      vx1(nx:nx+2,:,:) = 0.0_num
      vy1(nx:nx+2,:,:) = 0.0_num
      vz1(nx,:,:) = vz1(nx-1,:,:)
      vz1(nx+1,:,:) = vz1(nx,:,:)
      vz1(nx+2,:,:) = vz1(nx+2,:,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      vx1(:,-2:0,:) = 0.0_num
      vy1(:,-2:0,:) = 0.0_num
      vz1(:,0,:) = vz1(:,1,:)
      vz1(:,-1,:) = vz1(:,0,:)
      vz1(:,-2,:) = vz1(:,-1,:)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      vx1(:,ny:ny+2,:) = 0.0_num
      vy1(:,ny:ny+2,:) = 0.0_num
      vz1(:,ny,:) = vz1(:,ny-1,:)
      vz1(:,ny+1,:) = vz1(:,ny,:)
      vz1(:,ny+2,:) = vz1(:,ny+1,:)
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
      vx1(:,:,-2:0) = 0.0_num
      vy1(:,:,-2:0) = 0.0_num
      vz1(:,:,-2:0) = 0.0_num
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
      vx1(:,:,nz:nz+2) = 0.0_num
      vy1(:,:,nz:nz+2) = 0.0_num
      do i = -2, nx+2
          do j = -2, ny+2
          vz1(i, j, nz) = max(0.0,1.25*vz1(i,j,nz-1))
          end do
      end do
      vz1(:,:,nz+1) = vz1(:,:,nz)
      vz1(:,:,nz+1) = vz1(:,:,nz+1)
    END IF

    if (proc_z_min == MPI_PROC_NULL) then
        br = 13.0_num*bfield_fact; bl = 0.1_num*bfield_fact; kb = 15.0_num

        bzdy = (bz(0:nx+1,1:ny+1,0) - bz(0:nx+1,0:ny,0)) / dyb(1)
        bzdx = (bz(1:nx+1,0:ny+1,0) - bz(0:nx,0:ny+1,0)) / dxb(1)

        bzdy0 = 0.5_num*(bzdy(1:nx+1,0:ny) + bzdy(0:nx,0:ny))
        bzdx0 = 0.5_num*(bzdx(0:nx,1:ny+1) + bzdx(0:nx,0:ny))

        fact = (kb*(br-bl))/(bz(0:nx+1,0:ny+1,0) + 1d-10)*tanh(kb*(bz(0:nx+1,0:ny+1,0)- bl)/(br-bl+1d-10))
        fact0 = 0.25_num*(fact(0:nx,0:ny) + fact(1:nx+1,0:ny) + fact(0:nx,1:ny+1) + fact(1:nx+1, 1:ny+1))

        vx_surf(0:nx, 0:ny) = -fact0*bzdy0
        vy_surf(0:nx, 0:ny) = fact0*bzdx0

        vx1(0:nx,0:ny,0) = vx1(0:nx,0:ny,0) + shearing_fact*vx_surf
        vy1(0:nx,0:ny,0) = vy1(0:nx,0:ny,0) + shearing_fact*vy_surf
     end if


  END SUBROUTINE remap_v_bcs



END MODULE boundary
