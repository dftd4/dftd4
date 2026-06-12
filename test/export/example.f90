program tester
  use dftd4_version, only : get_dftd4_version
  implicit none
  character(len=:), allocatable :: version
  call get_dftd4_version(string=version)
  print *, version
end program tester
