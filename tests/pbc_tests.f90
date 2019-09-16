subroutine test_pbc_tools_convert
   use iso_fortran_env, wp => real64
   use pbc_tools
   implicit none
   real(wp),parameter :: cellpar(6) = &
      [9.09903131_wp, 9.09903131_wp, 30.46049560_wp, 90.0_wp, 90.0_wp, 120.0_wp]
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[5.9598811567890_wp,      2.1071361905157_wp,      3.6496669404404_wp,    &
      & 0.0000000000000_wp,      6.3214085715472_wp,      3.6496669404404_wp,    &
      & 0.0000000000000_wp,      0.0000000000000_wp,      7.2993338808807_wp],   &
      & shape(lattice))
   real(wp) :: dlat(3,3),rlat(3,3),volume,cpar(6),center(3)

   call cell_to_dlat(cellpar,dlat)
   print'(3g21.14)',dlat
   print*
   call cell_to_rlat(cellpar,rlat)
   print'(3g21.14)',rlat
   print*
   volume = cell_to_dvol(cellpar)
   print'(3g21.14)',volume
   print*

   call dlat_to_cell(lattice,cpar)
   print'(3g21.14)',cpar
   print*
   call dlat_to_rlat(lattice,rlat)
   print'(3g21.14)',rlat
   print*
   volume = dlat_to_dvol(lattice)
   print'(3g21.14)',volume
   print*

   center = get_center_dlat(lattice)
   print'(3g21.14)',center
   print*

end subroutine test_pbc_tools_convert
