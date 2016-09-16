MODULE correlation

  !======================================================================= 
  !                                                                       
  ! A MODULE TO OBTAIN THE CORRELATION FUNCTION OF A UNIGRID USING THE  FFTW
  ! LIBRARY 
  !                                                AUTHOR: GEORGE STEIN/MARCELO ALVAREZ
  !                                                LAST EDIT:     07.27.16
  !                                                                       
  !=======================================================================

  !-----------------------------------------------------------------------
  ! INCLUDE NECESSARY MODULES                                             
  !-----------------------------------------------------------------------
  
  use timing_diagnostics
  use transform
  use grafic_types
  use mpivars
  use grid

  !-----------------------------------------------------------------------
  ! IMPLICIT STATEMENT                                                   
  !-----------------------------------------------------------------------

  implicit none

  !---------------------------------------------------------------------
  ! DATA ARRAYS                                                       
  !---------------------------------------------------------------------
    
  double precision, allocatable :: corr_1dl(:)
  integer, allocatable :: corr_inbinl(:)
    
  !---------------------------------------------------------------------
  ! BOOKKEEPING AND GLOBAL VARIABLES 
  !---------------------------------------------------------------------

  real xr,yr,zr,xr2,yr2,zr2,rr,corrbin
  double precision corr_1d_val,corr_1d_lval
  double precision sigma,sigmal
  double precision sigma2lint
  integer corr_inbin_val,corr_inbin_lval

  complex ctemp
  
  logical fftback

  ! DECLARATIONS END
  !-----------------------------------------------------------------------
  
CONTAINS

  !=======================================================================
  
  SUBROUTINE correlate
    
    !---------------------------------------------------------------------
    ! pk1d calculates a 1d power spectrum and also sigma(R) for an R array
    ! of arbitrary length
    !---------------------------------------------------------------------
    use timing_diagnostics
   
    implicit none
    integer ii
    ! DECLARATIONS END
    !---------------------------------------------------------------------  

    !=====================================================================
    ! EXCUTABLE BEGIN
   
    ! KEEP COMPLEX FIELD AFTER PK INSTEAD OF C2R and R2C

    !---------------------------------------------------------------------
    ! SQUARE THE FIELD
    !---------------------------------------------------------------------

    if(myid==0) call timer_begin            
    do k=1,local_nz
       do j=1,n
          do i=1,n21
             index = int((k-1)*n+j-1,8)*n2p1+2*i-1 ! 1D INDEX
             delta(index) = delta(index)**2+delta(index+1)**2
             delta(index+1) = 0.0
          enddo
       enddo
    enddo

    if(myid==0) then
       timing_diagnostics_code='Square the Field for Correlation Function'
       call timer_end
    endif

    if(myid==0) call timer_begin            

    call fft_mpi(iplan,delta,total_local_sizes)
    
    if(myid==0) then
       timing_diagnostics_code='pk_fftw complex2real'
       call timer_end
    endif

    sigmal=0.
    do k=1,local_nz
       do j=1,n
          do i=1,n
             index = int((k-1)*n+j-1,8)*n2p1+i ! 1D INDEX 
             sigmal = sigmal + delta(index)**2
          enddo
       enddo
    enddo

    sigmal=sigmal/real(n)**3
    call mpi_allreduce(sigmal,sigma,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    if(myid==0) write(*,102) sigma


    !---------------------------------------------------------------------
    ! SPHERICALLY AVERAGE XI
    !---------------------------------------------------------------------

    rrmin = dr
    rrmax = dr*n/2 
    nrbin = n

    lrrmin = log10(rrmin)
    lrrmax = log10(rrmax)

    drb = (rrmax-rrmin) / nrbin
    dlr = (lrrmax-lrrmin) / nrbin

    if(.not.allocated(corr_1d))     allocate(corr_1d(nrbin))
    if(.not.allocated(corr_1dl))    allocate(corr_1dl(nrbin))
    if(.not.allocated(corr_inbin))  allocate(corr_inbin(nrbin))
    if(.not.allocated(corr_inbinl)) allocate(corr_inbinl(nrbin))

    corr_1dl    = 0.
    corr_inbinl = 0

    sigmal  = 0.
    sigma2l = 0.

    do k=1,local_nz
       zr=(k+local_z_start-1)*dr
       if(k+local_z_start.gt.n12) zr=boxsize-zr
       zr2=zr*zr
       
       do j=1,n
          yr=(j-1)*dr
          if(j.gt.n12) yr=boxsize-yr
          yr2=yr*yr
          
          do i=1,n
             xr=(i-1)*dr
             if(i.gt.n12) xr=boxsize-xr
             xr2=xr*xr

             index = int((k-1)*n+j-1,8)*n2p1+i ! 1D INDEX

             rr = sqrt(xr2+yr2+zr2)             
             if(bintype=='lin') then             
                ibin = floor( (rr - rrmin) / drb) + 1
             elseif(bintype=='log') then             
                ibin = floor( (log10(rr) - lrrmin) / dlr) + 1
             endif

             corrbin = delta(index) 
             sigmal  = sigmal + corrbin
!             corrbin = corrbin / d3r 
                          
             if(ibin>=1.and.ibin<=nrbin) then
                corr_1dl(ibin) = corr_1dl(ibin) + corrbin
                corr_inbinl(ibin) = corr_inbinl(ibin) + 1
             endif
          enddo
       enddo
    enddo
    
    !free up delta
    deallocate(delta)

    call mpi_allreduce(sigmal     ,sigma     ,1    ,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    call mpi_allreduce(corr_1dl   ,corr_1d   ,nrbin,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    call mpi_allreduce(corr_inbinl,corr_inbin,nrbin,mpi_integer         ,mpi_sum,mpi_comm_world,ierr)
    call mpi_barrier(mpi_comm_world,ierr)

    if(myid==0) then
    do i=1,nrbin
       if (corr_inbin(i).ne.0) corr_1d(i) = corr_1d(i)/corr_inbin(i)
    enddo
    endif
   
    if(myid==0) then
       timing_diagnostics_code='Spherically average correlation'
       call timer_end
    endif

    return
101 format(/,10x,"Sigma Pre-Correlation = ",f10.4,/)
102 format(/,10x,"Sigma Xi = ",f10.4,/)
    ! EXECUTABLE END
    !=====================================================================

  END SUBROUTINE correlate

END MODULE correlation

    
