!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Water
! MODULE        : ModuleMP 
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : August, 08th 2024
! REVISION      :  
! DESCRIPTION   : Velocity solution based on inverse problem / iterative approach
!
!------------------------------------------------------------------------------
!
!This program is free software; you can redistribute it and/or
!modify it under the terms of the GNU General Public License 
!version 2, as published by the Free Software Foundation.
!
!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with this program; if not, write to the Free Software
!Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!
!------------------------------------------------------------------------------

    
    module ModuleMP
    implicit none
    
    real :: v0, v1
    real :: MP_Diameter, MP_Density, WaterDensity
    integer :: CORRELATION_TYPE
    character(len=256) :: ReynoldsFilePath
    
    ! Constants
    real, parameter :: gravitational_constant = 9.87
    real, parameter :: water_dynamic_viscosity = 8.9E-4
    
    
contains
    !--------------------------------------------------------------------------
    subroutine Sett_vel(MP_Diameter, MP_Density, WaterDensity, v0, CORRELATION_TYPE, ReynoldsFilePath, v1)
    !--------------------------------------------------------------------------
    !Arguments-------------------------------------------------------------
        real :: MP_Diameter, MP_Density,WaterDensity  
        integer :: CORRELATION_TYPE
        real :: Re, CD, Err, Tolerance, Re0, Re1, water_dynamic_viscosity
        character(len=256) :: ReynoldsFilePath
        real :: v0, v1
    !--------------------------------------------------------------------------
     
        
        
        v0 = vel_ini(gravitational_constant, MP_Diameter, MP_Density, WaterDensity, water_dynamic_viscosity)
        Re0 = reynolds(MP_Diameter, WaterDensity, water_dynamic_viscosity, v0)
        ! Determine drag coefficient based on correlation type
        if     (CORRELATION_TYPE == 1) then
                CD = CalculateMorisson(Re0)  
        elseif (CORRELATION_TYPE == 2) then
                call ReadReynolds(Re0, CD)
        elseif (CORRELATION_TYPE == 3) then
                call cd_sphere(Re0, CD)
        endif
        
        ! Calculate velocity with drag
        v1 = vel_cd(CD, MP_Density, MP_Diameter, gravitational_constant, WaterDensity)
       ! print *, CD, "CD1antes"
        Re1 = reynolds(MP_Diameter, WaterDensity, water_dynamic_viscosity, v1)
      !  print *, v1, "v1antes"
        Err = Error(v0, v1)
      !  print *, Err, "Erro"
        ! Iteratively refine velocity until convergence
        do while (Err > 0.01)
          !  print *, Re0, "Re0"
            Re0 = Re1
         !   print *, Re0, "Re0"
         !   print *, v0, "VO lOOP"
            v0 = v1
           ! print *, v0, "VO lOOP COPY"
           if      (CORRELATION_TYPE == 1) then
                CD = CalculateMorisson(Re0)            
            elseif (CORRELATION_TYPE == 2) then
                call ReadReynolds(Re0, CD)
            elseif (CORRELATION_TYPE == 3) then
                call cd_sphere(Re0, CD)
            endif
        !    print *, CD, "CDDEPOIS"
            v1 = vel_cd(CD, MP_Density, MP_Diameter, gravitational_constant, WaterDensity)
         !   print *, v1, "Last v1 before loop"
            Err = Error(v0,v1)
        end do
        !print *, "Finish loop"

    end subroutine Sett_vel
    !-------------------------------------------------------------
    function vel_ini(gravitational_constant, MP_Diameter, MP_Density, WaterDensity, water_dynamic_viscosity) result(v0)
    !-------------------------------------------------------------
    real :: v0
    real :: MP_Diameter, MP_Density, WaterDensity
    real :: water_dynamic_viscosity
    real :: gravitational_constant
     water_dynamic_viscosity = 8.9E-4
     !-------------------------------------------------------------
     
        v0 = gravitational_constant * MP_Diameter ** 2 * (MP_Density - WaterDensity) / (18 * water_dynamic_viscosity)
        
    end function vel_ini
    !-------------------------------------------------------------
    function reynolds(MP_Diameter, WaterDensity, water_dynamic_viscosity, v0) result(Re)
    !-------------------------------------------------------------
    real :: Re
    real :: v0, MP_Diameter, WaterDensity, water_dynamic_viscosity
    !-------------------------------------------------------------
    
        Re = MP_Diameter * v0 * WaterDensity / water_dynamic_viscosity
        
    end function reynolds
    !-------------------------------------------------------------
    function vel_cd(CD, MP_Density, MP_Diameter, gravitational_constant, WaterDensity) result(v1)
    !-------------------------------------------------------------
    real :: v1
    real :: CD
    real :: MP_Diameter, WaterDensity, MP_Density, gravitational_constant
    !-------------------------------------------------------------
    
        v1 = ((4.0 / 3.0) * (MP_Density - WaterDensity) * MP_Diameter * gravitational_constant / WaterDensity * CD)**0.5
    end function vel_cd
    
    !-------------------------------------------------------------
    function CalculateMorisson(Re) result(CD)
    !-------------------------------------------------------------
    real :: Re
    real :: CD
    !-------------------------------------------------------------
        CD = (24.0 / Re) + ((2.6 * (Re / 5.0)) / (1.0 + (Re / 5.0)**1.52)) &
             + ((0.411 * ((Re / 2.63E5)**(-7.94))) / (1.0 + (Re / (2.63 * 1.0E5))**(-8.0))) &
             + ((0.25 * (Re / 1.0E6)) / (1.0 + (Re / 1.0E6)))
    end function CalculateMorisson
    !-------------------------------------------------------------
    subroutine ReadReynolds(RRey, CDNew)
        integer :: i, numLines
        real :: reynoldsNumber, CDS, RRey, CDNew
        real, dimension(:), allocatable :: Reynoldst, CDt
        
        ! Read Reynolds data from file and interpolate CD
    
        
        open(unit=13, file=ReynoldsFilePath, status='old')
   
        numLines = 0
        do
            read(13, *, iostat=i)
            if (i /= 0) exit
            numLines = numLines + 1
        end do
        
        allocate(Reynoldst(numLines), CDt(numLines))
        rewind(13)
        read(13, *) ! Skip header line
        
        do i = 1, numLines
            read(13, *) Reynoldst(i), CDt(i)
        end do
        
        close(13)
        
        
        deallocate(Reynoldst, CDt) ! Free allocated memory
    end subroutine ReadReynolds
    !-------------------------------------------------------------
    subroutine cd_sphere(Re, CDNew) 
    !-------------------------------------------------------------
        real :: Re, CDNew
        real, dimension(4) :: p
        real :: x1
    !-------------------------------------------------------------
        
        ! Calculate CD for a sphere based on Reynolds number
        
        if (Re <= 0.0) then
            CDNew = 0.0
        else if (Re > 8.0e6) then
            CDNew = 0.2
        else if (Re > 0.0 .and. Re <= 0.5) then
            CDNew = 24.0 / Re
        else if (Re > 0.5 .and. Re <= 100.0) then
            p = [4.22, -14.05, 34.87, 0.658]
            CDNew = polyval(p, 1.0 / Re)
        else if (Re > 100.0 .and. Re <= 1.0e4) then
            p = [-30.41, 43.72, -17.08, 2.41]
            CDNew = polyval(p, 1.0 / log10(Re))
        else if (Re > 1.0e4 .and. Re <= 3.35e5) then
            p = [-0.1584, 2.031, -8.472, 11.932]
            CDNew = polyval(p, log10(Re))
        else if (Re > 3.35e5 .and. Re <= 5.0e5) then
            x1 = log10(Re / 4.5e5)
            CDNew = 91.08 * x1**4 + 0.0764
        else
            p = [-0.06338, 1.1905, -7.332, 14.93]
            CDNew = polyval(p, log10(Re))
        end if
    end subroutine cd_sphere
    
    function Error(v0, v1) result(Err)
    real ::  Err
    real :: v0, v1
        Err = abs(v1 - v0) / v1 * 100
    end function Error
    !-------------------------------------------------------------
    function polyval(coeffs, x) result(resultPoly)
    !-------------------------------------------------------------
        real :: coeffs(:), x, resultPoly
        integer :: i
    !-------------------------------------------------------------
        
        resultPoly = coeffs(1)
        do i = 2, size(coeffs)
            resultPoly = resultPoly * x + coeffs(i)
        end do
    end function polyval

end module ModuleMP
