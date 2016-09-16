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
    real xoni,yoni,zoni,RTHi,idum
    if(myid==0) call timer_begin

    open(4, file=mergedfile,form='binary')
!    read(4) (idum,j=1,2),Non
    read(4) Non
    close(4)

    if(myid==0) write(*,101) Non
    !read in old ppruns
!    open(4, file=mergedfile,form='binary')
!    read(4) Non,(idum,j=1,3),&
!          (xon(i),yon(i),zon(i),(idum,j=1,4),&
!           RTH(i),(idum,j=1,12),i=1,Non)
!    close(4)

    !new ppruns
    open(4, file=mergedfile,form='binary')
    delta = 0.0
    read(4) Non,RTHmax,redshiftin
    do j=1,Non
       read(4) xoni,yoni,zoni,(idum,i=1,3),&
            RTHi,(idum,i=1,16)

       xoni = xoni + boxsize/2
       yoni = yoni + boxsize/2
       zoni = zoni + boxsize/2
       !Periodic wrap
       if(xoni < 0.0)     xoni = boxsize - xoni
       if(xoni > boxsize) xoni = xoni - boxsize
       if(yoni < 0.0)     yoni = boxsize - yoni
       if(yoni > boxsize) yoni = yoni - boxsize
       if(zoni < 0.0)     zoni = boxsize - zoni
       if(zoni > boxsize) zoni = zoni - boxsize
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
          mass = 4*3.14159265359/3*RTHi**3 * 2.775e11 * 0.7**2 * 0.3
          delta(index) = delta(index) + 1 !mass
       endif
    enddo
    close(4)

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

101 format(/,10x,i0, ' peaks read in',/)
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
      open(unit=33,file=mergedfile,access='stream')
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
