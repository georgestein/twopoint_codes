module pktable
  use grid
  use mpivars
  implicit none
  character *128 filepktab
  real pktabmin,pktabmax,pktabminl,pktabmaxl,dpktab,dtab
  real, allocatable:: tsav(:)
  integer npktab,itab
  real minstepk, stepk, minstepr, stepr
contains

! Takes the calculated power spectrum and correlation function from the grid and writes out to file 

  subroutine write_pktable

    implicit none

    double precision, allocatable:: powfileIO(:,:), corrfileIO(:,:),sigma2fileIO(:,:)     
    real kout
    if(myid==0) then
       allocate(powfileIO(nk,3))
       allocate(corrfileIO(nrbin,2))
       allocate(sigma2fileIO(numRbin,2))

       open(unit=2,file='sigma2_'//outcode)
       do i=1,numRbin
          sigma2fileIO(i,2) = sigma2(i)
          sigma2fileIO(i,1) = R_sigma(i)
          write(2,'(E14.7,x,E14.7)') sigma2fileIO(i,1:2)!e10.7
       enddo
       close(2)


       open(unit=1,file='power_'//outcode,status='unknown',form='formatted')
       write(1,*) Nhalo
       do i=1,nk

          !get kbin
          if(bintype=='lin') then
             kout = kmin+dkb*(i-1)          
          elseif(bintype=='log') then
             kout = 10**(lkmin+dlk*(i-1))          
          endif

          write(1,'(2E14.7,I8)') kout, pk_1d(i), pk_inbin(i)

       enddo       
       close(1)


       open(unit=3,file='correlation_'//outcode)
       do i=1,nrbin
          corrfileIO(i,2) = corr_1d(i)
          if(bintype=='lin') then
             corrfileIO(i,1) = rrmin+drb*(i-1)          
          elseif(bintype=='log') then
             corrfileIO(i,1) = 10**(lrrmin+dlr*(i-1))          
          endif
          write(3,'(E14.7,x,E14.7)') corrfileIO(i,1:2)!e10.7
       enddo       
       close(3)

       deallocate(powfileIO)
       deallocate(corrfileIO)
       deallocate(sigma2fileIO)

    endif
 
  end subroutine write_pktable
end module pktable
