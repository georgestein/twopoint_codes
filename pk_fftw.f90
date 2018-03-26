MODULE pk_fftw

  !======================================================================= 
  !                                                                       
  ! A MODULE TO OBTAIN THE POWER SPECTRUM OF A UNIGRID USING THE  FFTW
  ! LIBRARY 
  !                                                AUTHOR: MARCELO ALVAREZ
  !                                                LAST EDIT:     02.10.10
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
    
  double precision, allocatable :: pk_1dl(:)
  integer, allocatable :: pk_inbinl(:)
    
  !---------------------------------------------------------------------
  ! BOOKKEEPING AND GLOBAL VARIABLES 
  !---------------------------------------------------------------------

  real xk,yk,zk,xk2,yk2,zk2,rk,pkbin,kr
  double precision pk_1d_val,pk_1d_lval
  double precision sigma,sigmal
  double precision avg,avgl
  double precision sigma2lint


  complex ctemp
  
  logical fftback

  ! DECLARATIONS END
  !-----------------------------------------------------------------------
  
CONTAINS

  !=======================================================================
  
  SUBROUTINE pk1d
    
    !---------------------------------------------------------------------
    ! pk1d calculates a 1d power spectrum and also sigma(R) for an R array
    ! of arbitrary length
    !---------------------------------------------------------------------
    use timing_diagnostics
   
    implicit none
    integer ii
    real binround
    ! DECLARATIONS END
    !---------------------------------------------------------------------  

    !=====================================================================
    ! EXCUTABLE BEGIN
   
    !---------------------------------------------------------------------
    ! DO THE FFT
    !---------------------------------------------------------------------
    avgl = 0
    do k=1,local_nz
       do j=1,n
          do i=1,n
             index = int((k-1)*n+j-1,8)*n2p1+i ! 1D INDEX 
             avgl = avgl + delta(index)
          enddo
       enddo
    enddo
    avgl = avgl / real(n)**3
    call mpi_allreduce(avgl,avg,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    sigmal=0.
    do k=1,local_nz
       do j=1,n
          do i=1,n
             index = int((k-1)*n+j-1,8)*n2p1+i ! 1D INDEX 
             sigmal = sigmal + (delta(index)-avg)**2
          enddo
       enddo
    enddo
    sigmal=sigmal/real(n)**3
    call mpi_allreduce(sigmal,sigma,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    if(myid==0) write(*,101) avg,sqrt(sigma)


    if(myid==0) call timer_begin            

    call fft_mpi(plan,delta,total_local_sizes)
    call fft_mpi(plan,delta2,total_local_sizes)

    if(myid==0) then
       timing_diagnostics_code='pk_fftw real2complex'
       call timer_end
    endif

    delta  = delta / real(n)**3
    delta2 = delta2 / real(n)**3

    call mpi_barrier(mpi_comm_world,ierr)  


    if(myid==0) call timer_begin            

    !---------------------------------------------------------------------
    ! SPHERICALLY AVERAGE THE POWER
    !---------------------------------------------------------------------

    kmin = dk
    kmax = dk*n 
    lkmin = log10(kmin)
    lkmax = log10(kmax)

    nk = int(n/2)

    kbinmin = dk
    kbinmax = dk*nk

    lkbinmin = log10(kbinmin)
    lkbinmax = log10(kbinmax)

    dlk = (lkbinmax-lkbinmin) / (nk - 1)
    dkb = (kbinmax-kbinmin) / (nk - 1)   !these will make the bin centers

    if(.not.allocated(pk_1d))     allocate(pk_1d(nk))
    if(.not.allocated(pk_1dl))    allocate(pk_1dl(nk))
    if(.not.allocated(pk_inbin))  allocate(pk_inbin(nk))
    if(.not.allocated(pk_inbinl)) allocate(pk_inbinl(nk))

    pk_1dl=0.
    pk_inbinl=0

    sigmal=0.
    sigma2l=0.

    do k=1,local_nz
       zk=(k+local_z_start-1)*dk
       if(k+local_z_start.gt.n12) zk=zk-kmax
       zk2=zk*zk
       
       do j=1,n
          yk=(j-1)*dk
          if(j.gt.n12) yk=yk-kmax
          yk2=yk*yk
          
          do i=1,n21
             xk=(i-1)*dk
             xk2=xk*xk

             !To deal with wrapping of the fundamental plane. Only use half when i==1
             if( (i.ne.1) .or. (i.eq.1 .and. j.lt.n21 .and. j.gt.1) .or. &
                  (i.eq.1 .and. j.eq.1 .and. (k+local_z_start).lt.n21)) then

             index = int((k-1)*n+j-1,8)*n2p1+2*i-1 ! 1D INDEX
             rk = sqrt(xk2+yk2+zk2)             

             if(bintype=='lin') then
                ibin = nint( (rk - kmin) / dkb - tiny(kmin)) + 1  !nint is nearest int, which puts 
             elseif(bintype=='log') then                          !anything within 1/2dk into the bin
                ibin = nint( (log10(rk) - lkmin) / dlk - tiny(kmin)) + 1
             endif

!             pkbin = delta(index)**2 + delta(index+1)**2 ! variance of mode
             pkbin = delta(index)*delta2(index) + delta(index+1)*delta2(index+1) ! variance of mode

             sigmal = sigmal + pkbin * 2.
             if(rk.ne.0) then
                do ii=1,numRbin
                   kr = rk * R_sigma(ii)
                   sigma2lint = pkbin * 9.0 / kr**6 * (sin(kr)-kr*cos(kr))**2 
                   sigma2l(ii) = sigma2l(ii)+sigma2lint  * 2
                enddo
             endif
             pkbin = pkbin / d3k                               ! p(k) = variance per dk^3 
                          
             if(ibin>=1.and.ibin<=nk) then
                pk_1dl(ibin) = pk_1dl(ibin) + pkbin
                pk_inbinl(ibin) = pk_inbinl(ibin) + 1
             endif
          endif
          enddo
       enddo
    enddo

    call mpi_allreduce(sigmal   ,sigma   ,1      ,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    call mpi_allreduce(sigma2l  ,sigma2  ,numRbin,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    call mpi_allreduce(pk_1dl   ,pk_1d   ,nk     ,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    call mpi_allreduce(pk_inbinl,pk_inbin,nk     ,mpi_integer         ,mpi_sum,mpi_comm_world,ierr)
    call mpi_barrier(mpi_comm_world,ierr)

    if(myid==0) then
    do i=1,nk
       if (pk_inbin(i).ne.0) pk_1d(i) = pk_1d(i)/real(pk_inbin(i))
    enddo
    endif

    !---------------------------------------------------------------------
    ! IMPORTANT!!! IN COSMOLOGICAL CONVENTION, 
    !   variance per dlnk = Delta(k) = k^3*p_cosmo(k)/(2pi^2).
    ! HOWEVER, THE CURRENT DEFINTION OF p(k) is 
    !   p(k) = variance per dk^3 (SEE ABOVE LOOP),
    ! WHICH LEADS TO THE EXPRESSION
    !   variance per dlnk = Delta(k) = 4pi*k^3*p(k).
    ! SO TO GET THE "RIGHT ONE", p_cosmo(k), WE MUST MULTIPLY BY (2*pi)^3,
    !   p_cosmo(k) = (2*pi)^3 * p(k),
    ! TO COMPENSATE.
    !---------------------------------------------------------------------
    
    if(myid==0) pk_1d=pk_1d*(2.*pi)**3

    if(myid==0) then
       timing_diagnostics_code='Spherically average pk and sigma(R)'
       call timer_end
    endif

!    if(myid==0) call timer_begin            

!    call fft_mpi(iplan,delta,total_local_sizes)
    
!    if(myid==0) then
!       timing_diagnostics_code='pk_fftw complex2real'
!       call timer_end
!    endif

!    return
    
101 format(/5x," mean = ",1pe13.6,/5x,"sigma = ",1pe13.6/)

    ! EXECUTABLE END
    !=====================================================================

  END SUBROUTINE pk1d

END MODULE pk_fftw

    
