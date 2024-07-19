!=======================================================================

! Microplastics tracer within sea ice
! Microplastics (MP) are tracked in kg in the code.
! That means sources (ocean and atmosphere) are in kg/kg, and fluxes are in kg/m2/s,
! using the density of Neutrally buoyant MPs (1023 kg/m3) for the conversion
! from particles/m3 to kg/kg
! Neutral buoyant density if based on Mountford, A. S., & Morales Maqueda, M. A. (2021).
!
! Currently written with a max_mp of 6. Those could be different sources (Atlantic, Pacific, Arctic, Southern Ocean, Atmospheric deposition (by default in Artcic as neglected everywhere else)), in which case the scavenging factors set below should all be the same, as they would be behaving the same way physically. Or they could be used to try out different scavenging factors within the same model run, by setting the scavenging factors to different values
!
! Code authors: Alxandra Jahn (CU Boulder), David Bailey (NCAR), and Lingwei Li (CU Boulder)

      module icepack_microplastics

      use icepack_kinds
      use icepack_parameters, only: c0, c1, c2, puny, rhoi, rhos, hs_min
      use icepack_parameters, only: hi_ssl, hs_ssl
      use icepack_tracers, only: max_mp
      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      use icepack_zbgc_shared, only: kscavz

      implicit none

      private
      public :: update_microplastics, update_snow_bgc_mp

!---------------------------------------------------------------------
! Uptake and release coefficients for microplastic tracers 
! for melting and ice formation processes
!---------------------------------------------------------------------

      real (kind=dbl_kind), parameter, dimension(max_mp), public :: &
    !     krelm, krelsi,      & ! scavenging by melt water (krelm) and snow-ice formation (krelsi)
    !     kupfb, kupff          ! uptaking by ice formation (basal (kupfb) and frazil (kupff))

         krelm =   (/ c0, c0, c0, c0, c0, c0 /),      &       ! melt water release factor   
         krelsi =  (/ c0, c0, c0, c0, c0, c0 /),      &       ! snow-ice formation release factor
         kupfb =   (/ c1, c1, c1, c1, c1, c1 /),      &       ! basal ice formation uptake factor
         kupff =   (/ c1, c1, c1, c1, c1, c1 /)               ! frazil ice formation uptake factor
                                                        ! kupff used in icepack_therm_itd.F90 but set here
!=======================================================================

      contains

!=======================================================================

!  Increase microplastics in ice or snow surface due to surface deposition
!  and ocean uptake during freezing
!  and vertical cycling

      subroutine update_microplastics(dt, nilyr, nslyr,  n_mp,   &
                                meltt,    melts,      &
                                meltb,    congel,     &
                                snoice,   fsnow,      &
                                mpsno,    mpice,      &
                                aice_old,             &
                                vice_old, vsno_old,   &
                                vicen, vsnon, aicen,  &
                                fmp_atm, fmp_ocn, mp_ocn)


      integer (kind=int_kind), intent(in) :: &
         nilyr, nslyr, n_mp

      real (kind=dbl_kind), intent(in) :: &
         dt,       & ! time step
         meltt,      & ! top melt rate (m/step)
         melts,      & ! snow melt rate (m/step)
         meltb,      & ! bottom melt rate (m/step)
         congel,     & ! congelation ice growth rate   (m/step)
         snoice,     & ! ice thickness increase rate from snow ice formation  (m/step)
         fsnow,      & ! snowfall rate (kg/m^2/s of water)
         vicen,      & ! ice volume (m)
         vsnon,      & ! snow volume (m)
         aicen,      & ! ice area fraction
         aice_old,   & ! values prior to thermodynamic changes
         vice_old,   &
         vsno_old

      real (kind=dbl_kind), dimension(:), &
         intent(inout) :: &
         mp_ocn        ! ocean tracer concentrations of microplastics tracer (kg MP/kg sea water)


      real (kind=dbl_kind), dimension(:), &
         intent(in) :: &
         fmp_atm   ! microplastics atmospheric deposition rate (kg/m^2/s)

      real (kind=dbl_kind), dimension(:), &
         intent(inout) :: &
         fmp_ocn   ! microplastics flux to ocean (kg/m^2 s)

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         mpsno,  mpice    ! mass of microplastics in snow and ice (kg)

      !  local variables
      integer (kind=int_kind) :: k

      real (kind=dbl_kind) :: &
         dzssl,  dzssl_new,      & ! snow ssl thickness
         dzint,                  & ! snow interior thickness
         dzssli, dzssli_new,     & ! ice ssl thickness
         dzinti,                 & ! ice interior thickness
         dznew,                  & ! tracks thickness changes
         hs, hi,                 & ! snow/ice thickness (m)
         dhs_evap, dhi_evap,     & ! snow/ice thickness change due to evap
         dhs_melts, dhi_meltt,   & ! ... due to surface melt
         dhs_snoice, dhi_snoice, & ! ... due to snow-ice formation
         dhi_congel, dhi_meltb,  & ! ... due to bottom growth, melt
         hslyr, hilyr,           & ! snow, ice layer thickness (m)
         hslyr_old, hilyr_old,   & ! old snow, ice layer thickness (m)
         hs_old, hi_old,         & ! old snow, ice thickness (m)
         sloss1, sloss2,         & ! microplastics mass loss (kg/m^2)
         ar                        ! 1/aicen(i,j)

!      real (kind=dbl_kind), parameter, dimension(max_mp), public :: &
!         krelm, krelsi,      & ! scavenging by melt water (krelm) and snow-ice formation (krelsi)
!         kupfb, kupff          ! uptaking by ice formation (basal (kupfb) and frazil (kupff))
      
      real (kind=dbl_kind), dimension(n_mp) :: &
         mptot, mptot0,     & ! for conservation check
         focn_old             ! for conservation check

      ! The below assumes max_mp=6, need to expand table for more MPs
      !krelm =   (/ 0.0_dbl_kind, 0.0_dbl_kind, 0.0_dbl_kind, &               !melt water scavenging
      !             0.0_dbl_kind, 0.0_dbl_kind, 0.0_dbl_kind /)
      !krelsi = (/ 0.0_dbl_kind, 0.0_dbl_kind, 0.0_dbl_kind, &               !snow-ice formation scavenging
      !             0.0_dbl_kind, 0.0_dbl_kind, 0.0_dbl_kind /)
      !kupfb =   (/ 0.0_dbl_kind, 0.0_dbl_kind, 0.0_dbl_kind, &               !LLW: basal ice formation uptake factor
      !             0.0_dbl_kind, 0.0_dbl_kind, 0.0_dbl_kind /)
      !kupff =   (/ c1, c1, c1, &               !frazil ice formation uptake factor
       !            c1, c1, c1 /)               !LLW: uptake factor updated in therm_itd.F90

     character(len=*),parameter :: subname='(update_microplastics)'

    !-------------------------------------------------------------------
    ! initialize
    !-------------------------------------------------------------------

      focn_old(:) = fmp_ocn(:)

      hs_old    = vsno_old/aice_old
      hi_old    = vice_old/aice_old
      hslyr_old = hs_old/real(nslyr,kind=dbl_kind)
      hilyr_old = hi_old/real(nilyr,kind=dbl_kind)

      dzssl  = min(hslyr_old/c2, hs_ssl)
      dzssli = min(hilyr_old/c2, hi_ssl)
      dzint  = hs_old - dzssl
      dzinti = hi_old - dzssli

      if (aicen > c0) then
         ar = c1/aicen
      else ! ice disappeared during time step
         ar = c1/aice_old
      endif

      hs = vsnon*ar
      hi = vicen*ar
      dhs_melts  = -melts
      dhi_snoice = snoice
      dhs_snoice = dhi_snoice*rhoi/rhos
      dhi_meltt  = -meltt
      dhi_meltb  = -meltb
      dhi_congel = congel

      dhs_evap = hs - (hs_old + dhs_melts - dhs_snoice &
                              + fsnow/rhos*dt)
      dhi_evap = hi - (hi_old + dhi_meltt + dhi_meltb &
                              + dhi_congel + dhi_snoice)

      do k = 1, n_mp
         mptot0(k) = mpsno(k,2) + mpsno(k,1) &
                     + mpice(k,2) + mpice(k,1)
      enddo

    !-------------------------------------------------------------------
    ! evaporation
    !-------------------------------------------------------------------
      dzint  = dzint  + min(dzssl  + dhs_evap, c0)
      dzinti = dzinti + min(dzssli + dhi_evap, c0)
      dzssl  = max(dzssl  + dhs_evap, c0)
      dzssli = max(dzssli + dhi_evap, c0)

    !-------------------------------------------------------------------
    ! basal ice growth
    !-------------------------------------------------------------------

      if (dhi_congel > c0) then
         do k = 1,n_mp
            !LLW: uptake factor for basal ice formation, just for internal layer now
            !sloss1 = c0
            
            sloss2 = c0
            sloss2 = kupfb(k)*mp_ocn(k)*rhoi*dhi_congel*aicen  !

            mpice(k,2) = mpice(k,2) + sloss2
            fmp_ocn(k) = fmp_ocn(k) - sloss2/dt

            !if (dzint > puny)  &     ! Internal layer
            !     sloss2 = min(dhi_congel, dzint)  &
            !            *mp_ocn(k)/dzint*rhoi*aicen
            !if (dzssl > puny)  & ! Surface scattering layer
            !     sloss1 = max(dhi_congel-dzint, c0)  &
            !            *mp_ocn(k)/dzint*rhoi*aicen
            !mpice(k,1) = mpice(k,1) &
            !             + (c1-kupfb(k))*(sloss2+sloss1)
            !fmp_ocn(k) = fmp_ocn(k) &
            !             - kupfb(k)*(sloss2+sloss1)/dt
         enddo

            dzinti = dzinti + dhi_congel
      endif


      dzssli = dzssli + min(dzinti+dhi_meltb, c0)
      dzinti = max(dzinti+dhi_meltb, c0)
    !-------------------------------------------------------------------
    ! surface snow melt
    !-------------------------------------------------------------------
      if (-dhs_melts > puny) then
         do k = 1, n_mp
            sloss1 = c0
            sloss2 = c0
            if (dzssl > puny)  &
                 sloss1 = krelm(k)*mpsno(k,1)  &
                                  *min(-dhs_melts,dzssl)/dzssl
            mpsno(k,1) = mpsno(k,1) - sloss1
            if (dzint > puny)  &
                 sloss2 = krelm(k)*mpsno(k,2) &
                                  *max(-dhs_melts-dzssl,c0)/dzint
            mpsno(k,2) = mpsno(k,2) - sloss2
            fmp_ocn(k) = fmp_ocn(k) + (sloss1+sloss2)/dt
         enddo  ! n_mp

         ! update snow thickness
         dzint=dzint+min(dzssl+dhs_melts, c0)
         dzssl=max(dzssl+dhs_melts, c0)

         if ( dzssl <= puny ) then ! ssl melts away
            mpsno(:,2) = mpsno(:,1) + mpsno(:,2)
            mpsno(:,1) = c0
            dzssl = max(dzssl, c0)
         endif
         if (dzint <= puny ) then  ! all snow melts away
            mpice(:,1) = mpice(:,1) &
                         + mpsno(:,1) + mpsno(:,2)
            mpsno(:,:) = c0
            dzint = max(dzint, c0)
         endif
      endif

    !-------------------------------------------------------------------
    ! surface ice melt
    !-------------------------------------------------------------------
      if (-dhi_meltt > puny) then
         do k = 1, n_mp
            sloss1 = c0
            sloss2 = c0
            if (dzssli > puny)  &
                 sloss1 = krelm(k)*mpice(k,1)  &
                                  *min(-dhi_meltt,dzssli)/dzssli
            mpice(k,1) = mpice(k,1) - sloss1
            if (dzinti > puny)  &
                 sloss2 = krelm(k)*mpice(k,2)  &
                                  *max(-dhi_meltt-dzssli,c0)/dzinti
            mpice(k,2) = mpice(k,2) - sloss2
            fmp_ocn(k) = fmp_ocn(k) + (sloss1+sloss2)/dt
         enddo

         dzinti = dzinti + min(dzssli+dhi_meltt, c0)
         dzssli = max(dzssli+dhi_meltt, c0)
         if (dzssli <= puny) then   ! ssl ice melts away
            do k = 1, n_mp
               mpice(k,2) = mpice(k,1) + mpice(k,2)
               mpice(k,1) = c0
            enddo
            dzssli = max(dzssli, c0)
         endif
         if (dzinti <= puny) then   ! all ice melts away
            do k = 1, n_mp
               fmp_ocn(k) = fmp_ocn(k)  &
                            + (mpice(k,1)+mpice(k,2))/dt
               mpice(k,:)=c0
            enddo
            dzinti = max(dzinti, c0)
         endif
      endif

    !-------------------------------------------------------------------
    ! basal ice melt.  Assume all microplastics lost in basal melt
    !-------------------------------------------------------------------
      if (-dhi_meltb > puny) then
         do k=1,n_mp
            sloss1=c0
            sloss2=c0
            if (dzssli > puny)  &
                 sloss1 = max(-dhi_meltb-dzinti, c0)  &
                        *mpice(k,1)/dzssli
            mpice(k,1) = mpice(k,1) - sloss1
            if (dzinti > puny)  &
                 sloss2 = min(-dhi_meltb, dzinti)  &
                        *mpice(k,2)/dzinti
            mpice(k,2) = mpice(k,2) - sloss2
            fmp_ocn(k) = fmp_ocn(k) + (sloss1+sloss2)/dt
         enddo

         dzssli = dzssli + min(dzinti+dhi_meltb, c0)
         dzinti = max(dzinti+dhi_meltb, c0)
      endif

    !-------------------------------------------------------------------
    ! snowfall
    !-------------------------------------------------------------------
      if (fsnow > c0) dzssl = dzssl + fsnow/rhos*dt

    !-------------------------------------------------------------------
    ! snow-ice formation
    !-------------------------------------------------------------------
      if (dhs_snoice > puny) then
         do k = 1, n_mp
            sloss1 = c0
            sloss2 = c0
            if (dzint > puny)  &
                 sloss2 = min(dhs_snoice, dzint)  &
                        *mpsno(k,2)/dzint
            mpsno(k,2) = mpsno(k,2) - sloss2
            if (dzssl > puny)  &
                 sloss1 = max(dhs_snoice-dzint, c0)  &
                        *mpsno(k,1)/dzssl
            mpsno(k,1) = mpsno(k,1) - sloss1
            mpice(k,1) = mpice(k,1) &
                         + (c1-krelsi(k))*(sloss2+sloss1)
            fmp_ocn(k) = fmp_ocn(k) &
                         + krelsi(k)*(sloss2+sloss1)/dt
         enddo
         dzssl  = dzssl - max(dhs_snoice-dzint, c0)
         dzint  = max(dzint-dhs_snoice, c0)
         dzssli = dzssli + dhi_snoice
      endif

    !-------------------------------------------------------------------
    ! microplastic deposition
    !-------------------------------------------------------------------
      if (aicen > c0) then
         hs = vsnon * ar
      else
         hs = c0
      endif
      if (hs > hs_min) then    ! AJ: should this really be hs_min or 0?
         ! should use same hs_min value as in radiation
         do k=1,n_mp
            mpsno(k,1) = mpsno(k,1) &
                         + fmp_atm(k)*dt*aicen
         enddo
      else
         do k=1,n_mp
            mpice(k,1) = mpice(k,1) &
                         + fmp_atm(k)*dt*aicen
         enddo
      endif

    !-------------------------------------------------------------------
    ! redistribute microplastics within vertical layers - from aero
    !-------------------------------------------------------------------
      if (aicen > c0) then
         hs = vsnon  * ar     ! new snow thickness
         hi = vicen  * ar     ! new ice thickness
      else
         hs = c0
         hi = c0
      endif
      if (dzssl <= puny) then   ! nothing in SSL
         do k=1,n_mp
            mpsno(k,2) = mpsno(k,2) + mpsno(k,1)
            mpsno(k,1) = c0
         enddo
      endif
      if (dzint <= puny) then   ! nothing in Snow Int
         do k = 1, n_mp
            mpice(k,1) = mpice(k,1) + mpsno(k,2)
            mpsno(k,2) = c0
         enddo
      endif
      if (dzssli <= puny) then  ! nothing in Ice SSL
         do k = 1, n_mp
            mpice(k,2) = mpice(k,2) + mpice(k,1)
            mpice(k,1) = c0
         enddo
      endif

      if (dzinti <= puny) then  ! nothing in Ice INT
         do k = 1, n_mp
            fmp_ocn(k) = fmp_ocn(k) &
                         + (mpice(k,1)+mpice(k,2))/dt
            mpice(k,:)=c0
         enddo
      endif

      hslyr      = hs/real(nslyr,kind=dbl_kind)
      hilyr      = hi/real(nilyr,kind=dbl_kind)
      dzssl_new  = min(hslyr/c2, hs_ssl)
      dzssli_new = min(hilyr/c2, hi_ssl)

      if (hs > hs_min) then
         do k = 1, n_mp
            dznew = min(dzssl_new-dzssl, c0)
            sloss1 = c0
            if (dzssl > puny) &
                 sloss1 = dznew*mpsno(k,1)/dzssl ! not neccesarily a loss
            dznew = max(dzssl_new-dzssl, c0)
            if (dzint > puny) &
                 sloss1 = sloss1 + mpsno(k,2)*dznew/dzint
            mpsno(k,1) = mpsno(k,1) + sloss1
            mpsno(k,2) = mpsno(k,2) - sloss1
         enddo
      else
         mpice(:,1) = mpice(:,1)  &
                      + mpsno(:,1) + mpsno(:,2)
         mpsno(:,:) = c0
      endif

      if (vicen > puny) then ! AJ: may want a limit on hi instead?
         do k = 1, n_mp
            sloss2 = c0
            dznew = min(dzssli_new-dzssli, c0)
            if (dzssli > puny) &
                 sloss2 = dznew*mpice(k,1)/dzssli
            dznew = max(dzssli_new-dzssli, c0)
            if (dzinti > puny) &
                 sloss2 = sloss2 + mpice(k,2)*dznew/dzinti
            mpice(k,1) = mpice(k,1) + sloss2
            mpice(k,2) = mpice(k,2) - sloss2
         enddo
      else
         fmp_ocn(:) = fmp_ocn(:) + (mpice(:,1)+mpice(:,2))/dt
         mpice(:,:) = c0
      endif

! AJ: DO WE NEED TO ADD THIS, TO GET THE DISGNIOSTICS TO WORK. THIS IS FROM ISO CODE, NEEDS CHANGES IN merge_fluxes
! scale fmp_ocnn. It will be re-scaled by aicen later in merge_fluxes
!   if (aicen > puny) then
!       fmp_ocnn(:) = fmp_ocnn(:)/aicen
!   endif
    !-------------------------------------------------------------------
    ! check conservation
    !-------------------------------------------------------------------
      do k = 1, n_mp
         mptot(k) = mpsno(k,2) + mpsno(k,1) &
                    + mpice(k,2) + mpice(k,1)
         if ((mptot(k)-mptot0(k)) &
              - (   fmp_atm(k)*aicen &
              - (fmp_ocn(k)-focn_old(k)) )*dt  > puny) then

            write(warnstr,*) subname, 'microplastics tracer:  ',k
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'mptot-mptot0 ',mptot(k)-mptot0(k)
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'fmp_atm-fmp_ocn      ', &
                 (fmp_atm(k)*aicen-(fmp_ocn(k)-focn_old(k)))*dt
            call icepack_warnings_add(warnstr)
         endif
      enddo

    !-------------------------------------------------------------------
    ! check for negative values
    !-------------------------------------------------------------------

!echmod:  note that this does not test or fix all microplastics tracers
      if (mpice(1,1) < -puny .or. &
          mpice(1,2) < -puny .or. &
          mpsno(1,1) < -puny .or. &
          mpsno(1,2) < -puny) then

         write(warnstr,*) subname, 'microplastics negative in microplastics code'
         call icepack_warnings_add(warnstr)

         mpice(1,1) = max(mpice(1,1), c0)
         mpice(1,2) = max(mpice(1,2), c0)
         mpsno(1,1) = max(mpsno(1,1), c0)
         mpsno(1,2) = max(mpsno(1,2), c0)

      endif

      end subroutine update_microplastics

!=======================================================================
! UNCLEAR IF THIS IS NEEDED. Seems to be NEEDED FOR ZMP? CAN'T BE SAME NAME AS FOR AEROSOLS, CAUSES BUILD ERROR
! AJ. SO RENAMED TO _MP. But I get the same results for zmp as for zaero for the same forcing with calling this and without calling it for mp, so does it do somehting?
!  Increase microplastics in snow surface due to deposition
!  and vertical cycling : after update_microplastics

     subroutine update_snow_bgc_mp (dt,     nblyr,       &
                                nslyr,                 &
                                meltt,    melts,       &
                                meltb,    congel,      &
                                snoice,   nbtrcr,      &
                                fsnow,    ntrcr,       &
                                trcrn,    bio_index,   &
                                aice_old, zbgc_snow,   &
                                vice_old, vsno_old,    &
                                vicen,    vsnon,       &
                                aicen,    flux_bio_atm,&
                                zbgc_atm, flux_bio)

      integer (kind=int_kind), intent(in) :: &
         nbtrcr,             & ! number of distinct snow tracers
         nblyr,              & ! number of bio layers
         nslyr,              & ! number of snow layers
         ntrcr                 ! number of tracers

      integer (kind=int_kind), dimension (nbtrcr), intent(in) :: &
         bio_index

      real (kind=dbl_kind), intent(in) :: &
         dt                    ! time step

      real (kind=dbl_kind), intent(in) :: &
         meltt,    & ! top melt rate (m/step)
         melts,    & ! snow melt rate (m/step)
         meltb,    & ! bottom melt rate (m/step)
         congel,   & ! congelation ice growth rate (m/step)
         snoice,   & !snow-ice formation (m/step)
         fsnow,    &
         vicen,    & ! ice volume (m)
         vsnon,    & ! snow volume (m)
         aicen,    & ! ice area fraction
         aice_old, & ! values prior to thermodynamic changes
         vice_old, &
         vsno_old

      real (kind=dbl_kind),dimension(nbtrcr), intent(inout) :: &
         zbgc_snow, & ! microplastics contribution from snow to ice
         zbgc_atm,  & ! and atm to ice concentration * volume (kg m^3*m)
         flux_bio     ! total ocean tracer flux (mmol/m^2/s)

      real (kind=dbl_kind), dimension(nbtrcr), &
         intent(in) :: &
         flux_bio_atm   ! microplastics deposition rate (kg/m^2 s)

      real (kind=dbl_kind), dimension(ntrcr), &
         intent(inout) :: &
         trcrn       ! ice/snow tracer array

      !  local variables

      integer (kind=int_kind) ::  k, n

      real (kind=dbl_kind) :: &
         dzssl,  dzssl_new,      & ! snow ssl thickness
         dzint,  dzint_new,      & ! snow interior thickness
         hs,                     & ! snow thickness (m)
         dhs_evap,               & ! snow thickness change due to evap
         dhs_melts,              & ! ... due to surface melt
         dhs_snoice,             & ! ... due to snow-ice formation
         hslyr,                  & ! snow layer thickness (m)
         hslyr_old,              & ! old snow layer thickness (m)
         hs_old,                 & ! old snow thickness (m)
         dznew,                  & ! change in the snow sl (m)
         sloss1, sloss2,         & ! microplastics mass loss (kg/m^2)
         ar                        ! 1/aicen(i,j)

      real (kind=dbl_kind), dimension(nbtrcr) :: &
         mptot, mptot0,          & ! for conservation check (mmol/m^3)
         mp_cons      ,          & ! for conservation check (mmol/m^2)
         flux_bio_o                ! initial ocean tracer flux (mmol/m^2/s)

      real (kind=dbl_kind), dimension(nbtrcr,2) :: &
         mpsno,  & ! mass of microplastics in snow (kg)
         mpsno0    ! mass of microplastics in (kg)
                   ! both for diagnostic prints

      character(len=*),parameter :: subname='(update_snow_bgc)'

    !-------------------------------------------------------------------
    ! initialize
    !-------------------------------------------------------------------
         mpsno (:,:) = c0
         mpsno0(:,:) = c0
         mp_cons(:) = c0
         zbgc_snow(:) = c0
         zbgc_atm(:) = c0

         hs_old    = vsno_old/aice_old
         hslyr_old = hs_old/real(nslyr,kind=dbl_kind)

         dzssl  = min(hslyr_old/c2, hs_ssl)
         dzint  = hs_old - dzssl

         if (aicen > c0) then
            ar = c1/aicen
            hs = vsnon*ar
            dhs_melts  = -melts
            dhs_snoice = snoice*rhoi/rhos
         else ! ice disappeared during time step
            ar = c1
            hs = vsnon/aice_old
            dhs_melts  = -melts
            dhs_snoice = snoice*rhoi/rhos
         endif

         dhs_evap = hs - (hs_old + dhs_melts - dhs_snoice &
                                 + fsnow/rhos*dt)

         ! trcrn() has units kg/m^3

      if ((vsno_old .le. puny) .or. (vsnon .le. puny)) then

         do k=1,nbtrcr
            flux_bio(k) = flux_bio(k) +  &
                         (trcrn(bio_index(k)+ nblyr+1)*dzssl+ &
                          trcrn(bio_index(k)+ nblyr+2)*dzint)/dt
            trcrn(bio_index(k) + nblyr+1) = c0
            trcrn(bio_index(k) + nblyr+2) = c0
            zbgc_atm(k) = zbgc_atm(k) &
                                + flux_bio_atm(k)*dt
         enddo

      else

         do k=1,nbtrcr
            flux_bio_o(k) = flux_bio(k)
            mpsno (k,1) = trcrn(bio_index(k)+ nblyr+1) * dzssl
            mpsno (k,2) = trcrn(bio_index(k)+ nblyr+2) * dzint
            mpsno0(k,:) = mpsno(k,:)
            mptot0(k)   = mpsno(k,2) + mpsno(k,1)
         enddo

    !-------------------------------------------------------------------
    ! evaporation
    !-------------------------------------------------------------------
         dzint  = dzint  + min(dzssl  + dhs_evap, c0)
         dzssl  = max(dzssl  + dhs_evap, c0)

    !-------------------------------------------------------------------
    ! surface snow melt
    !-------------------------------------------------------------------
         if (-dhs_melts > puny) then
            do k = 1, nbtrcr
               sloss1 = c0
               sloss2 = c0
               if (dzssl > puny)  &
                  sloss1 = kscavz(k)*mpsno(k,1)  &
                                   *min(-dhs_melts,dzssl)/dzssl
               mpsno(k,1) = mpsno(k,1) - sloss1
               if (dzint > puny)  &
                   sloss2 = kscavz(k)*mpsno(k,2) &
                                    *max(-dhs_melts-dzssl,c0)/dzint
               mpsno(k,2) = mpsno(k,2) - sloss2
               zbgc_snow(k) = zbgc_snow(k) + (sloss1+sloss2)
            enddo  !

            ! update snow thickness
            dzint=dzint+min(dzssl+dhs_melts, c0)
            dzssl=max(dzssl+dhs_melts, c0)

            if ( dzssl <= puny ) then ! ssl melts away
               mpsno(:,2) = mpsno(:,1) + mpsno(:,2)
               mpsno(:,1) = c0
               dzssl = max(dzssl, c0)
            endif
            if (dzint <= puny ) then  ! all snow melts away
               zbgc_snow(:) = zbgc_snow(:) &
                                + max(c0,mpsno(:,1) + mpsno(:,2))
               mpsno(:,:) = c0
               dzint = max(dzint, c0)
            endif
         endif

    !-------------------------------------------------------------------
    ! snowfall
    !-------------------------------------------------------------------
         if (fsnow > c0) dzssl = dzssl + fsnow/rhos*dt

    !-------------------------------------------------------------------
    ! snow-ice formation
    !-------------------------------------------------------------------
         if (dhs_snoice > puny) then
            do k = 1, nbtrcr
               sloss1 = c0
               sloss2 = c0
               if (dzint > puny)  &
                  sloss2 = min(dhs_snoice, dzint)  &
                           *mpsno(k,2)/dzint
               mpsno(k,2) = mpsno(k,2) - sloss2
               if (dzssl > puny)  &
                  sloss1 = max(dhs_snoice-dzint, c0)  &
                           *mpsno(k,1)/dzssl
               mpsno(k,1) = mpsno(k,1) - sloss1
               zbgc_snow(k) = zbgc_snow(k) &
                               + (sloss2+sloss1)
            enddo
            dzssl  = dzssl - max(dhs_snoice-dzint, c0)
            dzint  = max(dzint-dhs_snoice, c0)
         endif

    !-------------------------------------------------------------------
    ! microplastics deposition
    !-------------------------------------------------------------------
         if (aicen > c0) then
            hs = vsnon * ar
         else
            hs = c0
         endif
         if (hs >= hs_min)  then !should this really be hs_min or 0?
                                  ! should use same hs_min value as in radiation
            do k=1,nbtrcr
               mpsno(k,1) = mpsno(k,1) &
                                + flux_bio_atm(k)*dt
            enddo
         else
            do k=1,nbtrcr
               zbgc_atm(k) = zbgc_atm(k) &
                                + flux_bio_atm(k)*dt
            enddo
         endif

    !-------------------------------------------------------------------
    ! redistribute microplastics within vertical layers
    !-------------------------------------------------------------------
         if (aicen > c0) then
            hs = vsnon  * ar     ! new snow thickness
         else
            hs = c0
         endif
         if (dzssl <= puny) then   ! nothing in SSL
            do k=1,nbtrcr
               mpsno(k,2) = mpsno(k,2) + mpsno(k,1)
               mpsno(k,1) = c0
            enddo
         endif
         if (dzint <= puny) then   ! nothing in Snow Int
            do k = 1, nbtrcr
               zbgc_snow(k) = zbgc_snow(k) + max(c0,mpsno(k,2))
               mpsno(k,2) = c0
            enddo
         endif

         hslyr      = hs/real(nslyr,kind=dbl_kind)
         dzssl_new  = min(hslyr/c2, hs_ssl)
         dzint_new  = hs - dzssl_new

         if (hs > hs_min) then !should this really be hs_min or 0?
            do k = 1, nbtrcr
               dznew = min(dzssl_new-dzssl, c0)
               sloss1 = c0
               if (dzssl > puny) &
                  sloss1 = dznew*mpsno(k,1)/dzssl ! not neccesarily a loss
                  dznew = max(dzssl_new-dzssl, c0)
               if (dzint > puny) &
                  sloss1 = sloss1 + mpsno(k,2)*dznew/dzint
               mpsno(k,1) = mpsno(k,1) + sloss1
               mpsno(k,2) = mpsno(k,2) - sloss1
            enddo
         else
            zbgc_snow(:) = zbgc_snow(:)  &
                             + max(c0,mpsno(:,1) + mpsno(:,2))
            mpsno(:,:) = c0
         endif

    !-------------------------------------------------------------------
    ! check conservation
    !-------------------------------------------------------------------
         do k = 1, nbtrcr
            mptot(k) = mpsno(k,2) + mpsno(k,1) &
                       + zbgc_snow(k) + zbgc_atm(k)
            mp_cons(k) = mptot(k)-mptot0(k) &
                          - (    flux_bio_atm(k) &
                          - (flux_bio(k)-flux_bio_o(k))) * dt
            if (mp_cons(k)  > puny .or. zbgc_snow(k) + zbgc_atm(k) < c0) then
               write(warnstr,*) subname, 'Conservation failure: microplasticss in snow'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'test microplastics 1'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'microplastics tracer:  ',k
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'mp_cons(k),puny:', mp_cons(k),puny
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'mptot,mptot0 ',mptot(k),mptot0(k)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, ' mpsno(k,2),mpsno(k,1) ', mpsno(k,2),mpsno(k,1)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'flux_bio_atm(k)*aicen*dt', &
                    flux_bio_atm(k)*aicen*dt
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'zbgc_snow(k)', &
                    zbgc_snow(k)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'zbgc_atm(k)', &
                    zbgc_atm(k)
               call icepack_warnings_add(warnstr)
            endif
         enddo

    !-------------------------------------------------------------------
    ! reload tracers
    !-------------------------------------------------------------------
         if (vsnon > puny) then
         do k = 1,nbtrcr
             trcrn(bio_index(k)+nblyr+1)=mpsno(k,1)/dzssl_new
             trcrn(bio_index(k)+nblyr+2)=mpsno(k,2)/dzint_new
         enddo
         else
         do k = 1,nbtrcr
            zbgc_snow(k) = (zbgc_snow(k) + mpsno(k,1) + mpsno(k,2))
            trcrn(bio_index(k)+nblyr+1)= c0
            trcrn(bio_index(k)+nblyr+2)= c0
         enddo
         endif
    !-------------------------------------------------------------------
    ! check for negative values
    !-------------------------------------------------------------------
         if (minval(mpsno(:,1)) < -puny  .or. &
            minval(mpsno(:,2)) < -puny) then

            write(warnstr,*) subname, 'Snow microplastics negative in update_snow_bgc'
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'aicen= '      ,aicen
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'vicen= '      ,vicen
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'vsnon= '      ,vsnon
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'viceold= '    ,vice_old
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'vsnoold= '    ,vsno_old
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'melts= '      ,melts
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'meltt= '      ,meltt
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'meltb= '      ,meltb
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'congel= '     ,congel
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'snoice= '     ,snoice
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'mp evap snow= '  ,dhs_evap
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'fsnow= '      ,fsnow
            call icepack_warnings_add(warnstr)
            do k = 1, nbtrcr
              write(warnstr,*) subname, 'NBTRCR value k = ', k
              call icepack_warnings_add(warnstr)
              write(warnstr,*) subname, 'mp snowssl (k)= '    ,mpsno0(k,1)
              call icepack_warnings_add(warnstr)
              write(warnstr,*) subname, 'mp new snowssl (k)= ',mpsno(k,1)
              call icepack_warnings_add(warnstr)
              write(warnstr,*) subname, 'mp snowint (k)= '    ,mpsno0(k,2)
              call icepack_warnings_add(warnstr)
              write(warnstr,*) subname, 'mp new snowint(k)= ',mpsno(k,2)
              call icepack_warnings_add(warnstr)
              write(warnstr,*) subname, 'flux_bio_atm(k)= ' , flux_bio_atm(k)
              call icepack_warnings_add(warnstr)
              write(warnstr,*) subname, 'zbgc_snow(k)= '  ,zbgc_snow(k)
              call icepack_warnings_add(warnstr)
              write(warnstr,*) subname, 'zbgc_atm(k)= '  ,zbgc_atm(k)
              call icepack_warnings_add(warnstr)

              do n = 1,2
                trcrn(bio_index(k)+nblyr+n)=max(trcrn(bio_index(k)+nblyr+n), c0)
              enddo
            enddo
         endif
        endif

      end subroutine update_snow_bgc_mp

!=======================================================================

      end module icepack_microplastics

!=======================================================================
