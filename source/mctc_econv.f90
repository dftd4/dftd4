module mctc_econv
   use iso_fortran_env, only : wp => real64
   implicit none
!  convert bohr to Ångström and back
   real(wp),parameter :: autoaa = 0.52917726_wp
   real(wp),parameter :: aatoau = 1.0_wp/autoaa
!  convert Hartree to eV and back
   real(wp),parameter :: autoev = 27.21138505_wp
   real(wp),parameter :: evtoau = 1.0_wp/autoev
!  convert Hartree to kcal/mol and back
   real(wp),parameter :: autokcal = 627.50947428_wp
   real(wp),parameter :: kcaltoau = 1.0_wp/autokcal
!  convert Hartree to kJ/mol and back (why?)
   real(wp),parameter :: autokj = 2625.49964038_wp
   real(wp),parameter :: kjtoau = 1.0_wp/autokj
end module mctc_econv
