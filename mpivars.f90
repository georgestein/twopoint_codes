module mpivars

  use grafic_types

  !-----------------------------------------------------------------------
  ! MPI AND FFTW VARIABLES AND CONSTANTS                                  
  !-----------------------------------------------------------------------

  integer ierr, myid, ntasks
  integer(i8b) plan, iplan
  integer local_z_start,local_nz,local_y_start,local_ny,total_local_sizes
  integer kl
  integer, dimension(MPI_STATUS_SIZE) :: status
end module mpivars
