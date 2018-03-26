program external_field

  !=======================================================================
  !
  ! GENERATES POWERSPECTRUM AND CORRELATION FUNCTION OF HALO CATALOGUE
  !
  !                                                AUTHOR: GEORGE STEIN
  !                                                LAST EDIT:     04.10.15
  !
  !=======================================================================

  !=======================================================================
  ! DECLARATIONS BEGIN  

  !-----------------------------------------------------------------------
  ! INCLUDE NECESSARY MODULES
  !-----------------------------------------------------------------------
  
  use grid
  use memory_management
  use timing_diagnostics
  use normalize
  use mpivars
  use pk_fftw
  use correlation
  use pks2grid
  use pktable
  !-----------------------------------------------------------------------
  ! IMPLICIT STATEMENT
  !-----------------------------------------------------------------------

  implicit none

  real, allocatable :: work(:)
  integer ifield

  !-----------------------------------------------------------------------
  ! EXTERNAL FUNCTIONS
  !-----------------------------------------------------------------------

  real r4arg
  integer i4arg, indlnb, iargc

  ! DATE AND TIME
  character(8)          :: dnt_date
  character(10)         :: dnt_time
  character(5)          :: dnt_zone
  integer, dimension(8) :: dnt_values

  ! DECLARATIONS END
  !=======================================================================

  !=======================================================================
  ! EXCUTABLE BEGIN
  
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call mpi_comm_size(mpi_comm_world,ntasks,ierr)

  !-----------------------------------------------------------------------
  ! GET COMMANDLINE PARAMETERS
  ! ----------------------------------------------------------------------

  if(iargc().lt.5) then
     if(myid==0) write(*,99)
99   format('usage: powerspectra <pksfile1> <outfile> <Lbox in [Mpc]> <nmesh> <fmt=: 0=peaks,1=field, 2=peaks x peaks, 3=field x peaks> <Mmin halo> <Nhalocut> <pksfile2>')
     call mpi_finalize(ierr)  
     stop 
  endif
     
  call getarg(1,mergedfile1)
  call getarg(2,outcode)
  boxsize  = r4arg(3,1.e2)
  n        = i4arg(4,256)
  fmt      = i4arg(5,0)
  Minmass  = r4arg(6,0.)
  Nhalocut = r4arg(7,0.)
  if(fmt>1) call getarg(8,mergedfile2)

  dngrid = 1     !currently not working for new version
  n = n/dngrid

  gridsize = n

  bintype  = 'lin'
  Rmin     = boxsize/n !for calculating sigma(R) 
  Rmax     = boxsize/2
  numRbin  = 10
  Rbinstep = (log10(Rmax)-log10(Rmin))/(numRbin-1)


  ! REPORT DATE AND TIME
  call date_and_time(dnt_date,dnt_time,dnt_zone,dnt_values)
  if(myid==0) write(*,11) dnt_values(3),dnt_values(2),&
       dnt_values(1),dnt_values(5),dnt_values(6),&
       dnt_values(7),dnt_zone

  ! ----------------------------------------------------------------------
  ! ALLOCATE AND DETERMINE LOCAL SIZES FROM PLAN
  ! ----------------------------------------------------------------------

  if(myid==0) call timer_begin

  call rfftw3d_f77_mpi_create_plan(plan,mpi_comm_world,n,n,n,&
       & fftw_real_to_complex,fftw_estimate)
  call rfftwnd_f77_mpi_local_sizes(plan,local_nz,local_z_start, &
       & local_ny, local_y_start, total_local_sizes)

  if(myid==0) then
     timing_diagnostics_code = 'plan creation'
     call timer_end
  endif

  call mpi_barrier(mpi_comm_world,ierr)

  call rfftw3d_f77_mpi_create_plan(iplan,mpi_comm_world,n,n,n,&
       & fftw_complex_to_real,fftw_estimate)
  call rfftwnd_f77_mpi_local_sizes(plan,local_nz,local_z_start, &
       & local_ny, local_y_start, total_local_sizes)

  offset = int(4*n*n,8)*(local_z_start)
  length = int(4*n*n,8)*(local_nz)

  allocate(delta(total_local_sizes))
  allocate(delta2(total_local_sizes))

  n12=n/2
  n21=n12+1
  n2p1=2*(n12+1)

  dk=2.*pi/boxsize
  kmax=n*dk
  d3k=dk**3

  dr = boxsize/n
  d3r = dr**3
  allocate(sigma2(numRbin))
  allocate(sigma2l(numRbin))
  allocate(R_sigma(numRbin))
  do i=1,numRbin
     R_sigma(i) = 10**(log10(Rmin)+Rbinstep*(i-1))
     sigma2(i)=0.
     sigma2l(i)=0.
  enddo

  ! ----------------------------------------------------------------------
  ! Get POWER SPECTRUM AND CORRELATION FUNCTION
  ! ----------------------------------------------------------------------

  if(fmt==0) then 
     call gridpks
     delta2 = delta

  elseif(fmt==1) then
     call readfield
     delta2 = delta

  elseif(fmt>1) then      
     call gridpks              
     delta2 = delta  

     mergedfile1 = mergedfile2
     if(fmt==2) call gridpks   !load mergedfile2 into delta
     if(fmt==3) call readfield 

  endif
  call pk1d
!  call correlate
  call write_pktable

  ! ----------------------------------------------------------------------
  ! DONE
  ! ----------------------------------------------------------------------

  if(myid==0) then
     call date_and_time(dnt_date,dnt_time,dnt_zone,dnt_values)
     write(*,12) dnt_values(3),dnt_values(2),dnt_values(1),&
       dnt_values(5),dnt_values(6),dnt_values(7),&
       dnt_zone
  endif

  call mpi_finalize(ierr)

 11 format(/,3x,61('-'),/,3x,'Power Spectra running on',/,&
         3x,i0.2,'.',i0.2,'.',i4,1x,'at ',&
         i0.2,':',i0.2,':',i0.2,1x,'UTC',a,/,&
         3x,61('-'),/)

 12 format(/,3x,61('-'),/,3x,&
         'Power Spectra finished running on',/,&
         3x,i0.2,'.',i0.2,'.',i4,1x,'at ',&
         i0.2,':',i0.2,':',i0.2,1x,'UTC',a,/,&
         3x,61('-'),/)

end program external_field

