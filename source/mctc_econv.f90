!> @brief unit conversion factors
module mctc_econv
   use iso_fortran_env, only : wp => real64
   implicit none
!  convert bohr to Ångström and back
   real(wp),parameter :: autoaa = 0.52917726_wp  !< convert bohr to Ångström
   real(wp),parameter :: aatoau = 1.0_wp/autoaa  !< convert Ångström to bohr
!  convert Hartree to eV and back
   real(wp),parameter :: autoev = 27.21138505_wp !< convert Hartree to eV
   real(wp),parameter :: evtoau = 1.0_wp/autoev  !< convert eV to Hartree
!  convert Hartree to kcal/mol and back
   real(wp),parameter :: autokcal = 627.50947428_wp !< convert Hartree to kcal/mol
   real(wp),parameter :: kcaltoau = 1.0_wp/autokcal !< convert kcal/mol to Hartree
!  convert Hartree to kJ/mol and back (why?)
   real(wp),parameter :: autokj = 2625.49964038_wp  !< convert Hartree to kJ/mol
   real(wp),parameter :: kjtoau = 1.0_wp/autokj     !< convert kJ/mol to Hartree
end module mctc_econv
