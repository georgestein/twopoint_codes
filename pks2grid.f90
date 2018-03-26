module pks2grid
  
  use grafic_types
  use grid
  use mpivars
  use timing_diagnostics
    
  implicit none
  integer Non
  real RTHmax,redshiftin
  double precision mass, mean, meanl
contains

  subroutine gridpks

    implicit none
    integer i,j,k
    integer xpki,ypki,zpki
    real xoni,yoni,zoni,RTHi,idum,massi

    delta  = 0.0
    Nhalol = 0
    open(4, file=mergedfile1,access='stream')

    read(4) Non,RTHmax,redshiftin
    
    if(Nhalocut>0) Non = Nhalocut !note this assumes input halos are already ranked in mass

    if(myid==0) write(*,101) Non
    do j=1,Non
       read(4) xoni,yoni,zoni,(idum,i=1,3),&
            RTHi,(idum,i=1,4)

       massi = 4./3 * 3.14159 * RTHi**3 * 2.775e11 * 0.25 * 0.7**2

       xoni = xoni + boxsize/2
       yoni = yoni + boxsize/2
       zoni = zoni + boxsize/2

       !Periodic wrap
       xoni = mod(xoni+boxsize,boxsize)
       yoni = mod(yoni+boxsize,boxsize)
       zoni = mod(zoni+boxsize,boxsize)

       if(massi < Minmass) cycle

       xpki = int(xoni/dr)+1
       ypki = int(yoni/dr)+1
       zpki = int(zoni/dr)+1

       !Add halo to density field
       if((zpki > local_nz*myid).and.(zpki<=local_nz*(myid+1))) then
          zpki = zpki - local_nz*myid
          index = int((zpki-1)*n+ypki-1,8)*n2p1+xpki
          if(index<1.or.index>( ((n-1)*n+n-1)*n2p1+n)) then
             write(*,*) xoni,yoni,zoni,index,xpki,ypki,zpki
          endif
          delta(index) = delta(index) + 1
!          delta(index) = delta(index) + massi
          Nhalol       = Nhalol + 1
       endif
    enddo
    close(4)

    call mpi_allreduce(Nhalol,Nhalo,1,mpi_integer,mpi_sum,mpi_comm_world,ierr)
    if(myid==0) write(*,102) Nhalo

    meanl = 0.0
    do k=1,local_nz
       do j=1,n
          do i=1,n
             index = int((k-1)*n+j-1,8)*n2p1+i
             meanl = meanl + delta(index)
          enddo
       enddo
    enddo
    meanl = meanl/real(n)**3

    call mpi_allreduce(meanl,mean,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

    do k=1,local_nz
       do j=1,n
          do i=1,n
             index = int((k-1)*n+j-1,8)*n2p1+i
             delta(index) = delta(index)/mean - 1
          enddo
       enddo
    enddo

    if(myid==0) then
       timing_diagnostics_code='pks2grid load peaks'
       call timer_end
    endif

    return

101 format(/,10x,i0, ' Halos in file: ',/)
102 format(/,10x,i0, ' Halos after mass cut: ',/)
    end subroutine gridpks

!==================================================================
    subroutine readfield

      implicit none
      integer(kind=mpi_offset_kind) offset_slab, offset_bytes, iorig
      integer idn,jdn,kdn
      real deltain

      !now can downgrid input field
      nin = n*dngrid
      offset_slab = nin*nin*local_z_start*dngrid
      offset_bytes = offset_slab*int(4,8)+1
      
      !read in data
      open(unit=33,file=mergedfile1,access='stream')
      if (dngrid==1) then
         read(33,pos=offset_bytes) (delta(i),i=1,n*n*local_nz)
      else
         do k=1,local_nz*dngrid
            do j=1,n*dngrid
               do i=1,n*dngrid
                  read(33,pos=offset_bytes) deltain
                  idn = int((i+dngrid-1)/dngrid)
                  jdn = int((j+dngrid-1)/dngrid)
                  kdn = int((k+dngrid-1)/dngrid)

                  index = int((kdn-1)*n+jdn-1,8)*n2p1+idn ! 1D INDEX
                  delta(index) = delta(index)+deltain
                  
                  offset_bytes=offset_bytes+4
               enddo
            enddo
         enddo
      endif
      close(33)
      delta = delta / dngrid**3 

      !add back padding for fft's
      iorig=local_nz*n*n+1
      do k=local_nz,1,-1
         do j=n,1,-1
            do i=n,1,-1
               iorig=iorig-1
               index = int((k-1)*n+j-1,8)*n2p1+i ! 1D INDEX
               delta(index)=delta(iorig)
            enddo
         enddo
      enddo



    return

    end subroutine readfield

end module pks2grid
