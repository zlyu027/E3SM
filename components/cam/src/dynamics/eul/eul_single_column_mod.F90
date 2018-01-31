module eul_single_column_mod
!---------------------------------------------------------
!
! Module for the Eulerian single column model

use prognostics
use scamMod
use constituents, only: cnst_get_ind
use time_manager, only: get_nstep

implicit none

public scm_setinitial
public eul_post_forecast

!===========================================================================
contains
!===========================================================================

subroutine scm_setinitial

  implicit none

  integer i, j, k, thelev
  integer inumliq, inumice, icldliq, icldice

  if (.not. use_camiop) then
    call cnst_get_ind('NUMLIQ', inumliq, abort=.false.)
    call cnst_get_ind('NUMICE', inumice, abort=.false.)
    call cnst_get_ind('CLDLIQ', icldliq)
    call cnst_get_ind('CLDICE', icldice)

    if (have_ps) ps(1,1,n3)=psobs

    thelev=1
    do k=1,PLEV
      if (tobs(k) .ne. 0) then
        thelev=k
        go to 1000
      endif
    enddo

1000 continue
    if (get_nstep() .le. 1) then 
      do k=thelev,PLEV
        if (have_t) t3(1,k,1,n3)=tobs(k)
        if (have_q) q3(1,k,1,1,n3)=qobs(k)
        if (have_u) u3(1,k,1,n3)=uobs(k)
        if (have_v) v3(1,k,1,n3)=vobs(k)
        if (have_numliq) q3(1,k,inumliq,1,n3)=numliqobs(k)
        if (have_cldliq) q3(1,k,icldliq,1,n3)=cldliqobs(k)
        if (have_numice) q3(1,k,inumice,1,n3)=numiceobs(k)
        if (have_cldice) q3(1,k,icldice,1,n3)=cldiceobs(k)
      enddo 

      do k=1,thelev-1
        t3(1,k,1,n3)=t3(1,k,1,1)
        q3(1,k,1,1,n3)=q3(1,k,1,1,1)
        u3(1,k,1,n3)=u3(1,k,1,1)
        v3(1,k,1,n3)=v3(1,k,1,1)
        tobs(k)=t3(1,k,1,1)
        qobs(k)=q3(1,k,1,1,1)
      enddo
    else
      tobs(:)=t3(1,:,1,n3)
      qobs(:)=q3(1,:,1,1,n3)
    endif

  endif

end subroutine scm_setinitial

subroutine eul_post_forecast(lat, psm1, qfcst, cwava, &
                   etamid  ,qminus  ,hw2al   ,hw2bl, &
                   hw3al   ,hw3bl   ,hwxal   ,hwxbl, &
                   nlon)

   use shr_kind_mod,   only: r8 => shr_kind_r8, i8 => shr_kind_i8
   use commap
   use pmgrid
   use scamMod
   use physconst,      only: rga
   use constituents,   only: pcnst, cnst_get_ind
   use eul_control_mod

   implicit none

   real(r8), intent(in) :: psm1(plon)          ! surface pressure (time n)
   real(r8), intent(inout) :: qfcst(plon,plev,pcnst)
   real(r8), intent(in) :: qminus(plon,plev,pcnst)
   real(r8), intent(in) :: cwava              ! normalization factor (1/g*plon)
   real(r8), intent(in) :: etamid(plev)
   real(r8), intent(out) :: hw2al(pcnst)  ! -
   real(r8), intent(out) :: hw2bl(pcnst)  !  | lat contributions to components
   real(r8), intent(out) :: hw3al(pcnst)  !  | of slt global mass integrals 
   real(r8), intent(out) :: hw3bl(pcnst)  ! -
   real(r8), intent(out) :: hwxal(pcnst,4)
   real(r8), intent(out) :: hwxbl(pcnst,4)
   real(r8) sum
   real(r8) dotproda           ! dot product
   real(r8) dotprodb           ! dot product
   real(r8) pdelb(plon,plev)
   real(r8) hcwavaw
   integer lat
   integer nlon
   integer i,k,m

   call pdelb0 (psm1, pdelb, nlon)
!
! Accumulate mass integrals
!
   sum = 0._r8
   do i=1,nlon
      sum = sum + psm1(1)
   end do
   tmass(lat) = w(lat)*rga*sum/nlon

!
! Add spegrd calculations to fix water mass
!
!
! Calculate SLT moisture and constituent integrals
!
   hcwavaw = 0.5_r8*cwava*w(lat)
   do m=1,pcnst
      hw2al(m) = 0._r8
      hw2bl(m) = 0._r8
      hw3al(m) = 0._r8
      hw3bl(m) = 0._r8
      hwxal(m,1) = 0._r8
      hwxal(m,2) = 0._r8
      hwxal(m,3) = 0._r8
      hwxal(m,4) = 0._r8
      hwxbl(m,1) = 0._r8
      hwxbl(m,2) = 0._r8
      hwxbl(m,3) = 0._r8
      hwxbl(m,4) = 0._r8
      do k=1,plev
         dotproda = 0._r8
         dotprodb = 0._r8
         do i=1,nlon
            dotproda = dotproda + qfcst(i,k,m)*pdela(i,k)
            dotprodb = dotprodb + qfcst(i,k,m)*pdelb(i,k)
         end do
         hw2al(m) = hw2al(m) + hcwavaw*dotproda
         hw2bl(m) = hw2bl(m) + hcwavaw*dotprodb
      end do
   end do

   call qmassd (cwava, etamid, w(lat), qminus, qfcst, &
                pdela, hw3al, nlon)

   call qmassd (cwava, etamid, w(lat), qminus, qfcst, &
                pdelb, hw3bl, nlon)

   if (pcnst.gt.1) then
      call xqmass (cwava, etamid, w(lat), qminus, qfcst, &
                   qminus, qfcst, pdela, pdelb, hwxal, &
                   hwxbl, nlon)
   end if

   return

end subroutine eul_post_forecast

end module eul_single_column_mod
