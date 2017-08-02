module grid

  use grafic_types

  IMPLICIT NONE

  character *128 code, outcode
  character *128 mergedfile1, mergedfile2, bintype
  real Minmass
  integer Nhalo, Nhalol

  real, allocatable :: delta(:),delta_sub(:),delta_sub_local(:)
  real, allocatable :: delta2(:),delta2_sub(:),delta2_sub_local(:)

  integer(i8b) i,j,k,ii,jj,kk,n,fmt,dngrid,nin
  real dx,boxsize,local_z_origin
  integer(i8b) offset,length,index
  integer total_sub_size    

  integer(i8b) n12,n21,n2p1
  integer nsub12,nsub21,nsub2p1,nsub3
  integer ngrids  
  integer llx,lly,llz
  
  real xo,yo,zo
  integer gridnum 

  real kx,ky,kz,kmax,ak,dk,d3k,dq,kmin,dkb, kbinmin, kbinmax 
  real dr,d3r,rrmin, rrmax, drb, lrrmin,lrrmax, dlr

!from ng_grid
  real, allocatable :: kpk_1d(:), R_sigma(:)
  integer, allocatable :: pk_inbin(:), corr_inbin(:)
  double precision, allocatable :: sigma2(:), sigma2l(:),pk_1d(:),corr_1d(:)
  real Rmin, Rmax, Rbinstep, rho_crit, omegaM_in, Rsmooth
  integer numRbin
  integer nk, ibin, nrbin

 real lkmin, lkmax, dlk, lkbinmin, lkbinmax


end module grid

