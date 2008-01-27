MODULE from_wrf_to_grads

CONTAINS

!--------------------------------------------------------

  subroutine interp_to_z( data_in , nx_in , ny_in , nz_in , &
                              data_out, nx_out, ny_out, nz_out, &
                              z_in, z_out, missing_value, &
                              vertical_type, debug   )
  implicit none
  integer, intent(in)                                  :: nx_in , ny_in , nz_in 
  integer, intent(in)                                  :: nx_out, ny_out, nz_out
  real, intent(in)                                     :: missing_value
  real, dimension(nx_in , ny_in , nz_in ), intent(in ) :: data_in, z_in
  real, dimension(nx_out, ny_out, nz_out), intent(out) :: data_out
  real, dimension(nz_out), intent(in)                  :: z_out
  logical, intent(in)                                  :: debug
  character (len=1)                                    :: vertical_type

  real, dimension(nz_in)                               :: data_in_z, zz_in
  real, dimension(nz_out)                              :: data_out_z

  integer :: i,j,k

    do i=1,nx_in
    do j=1,ny_in

      do k=1,nz_in
        data_in_z(k) = data_in(i,j,k)
        zz_in(k) = z_in(i,j,k)
      enddo
      do k=1,nz_out
        data_out_z(k) = data_out(i,j,k)
      enddo

      call interp_1d( data_in_z, zz_in, nz_in, &
                      data_out_z, z_out, nz_out, &
                      vertical_type, missing_value )

      do k=1,nz_out
        data_out(i,j,k) = data_out_z(k)
      enddo


    enddo
    enddo

  end subroutine interp_to_z

!----------------------------------------------

  subroutine interp_1d( a, xa, na, &
                        b, xb, nb, vertical_type, missing_value )
  implicit none
  integer, intent(in)              ::  na, nb
  real, intent(in), dimension(na)  :: a, xa
  real, intent(in), dimension(nb)  :: xb
  real, intent(out), dimension(nb) :: b
  real, intent(in)                 :: missing_value

  integer                          :: n_in, n_out
  logical                          :: interp
  real                             :: w1, w2
  character (len=1)                :: vertical_type


  if ( vertical_type == 'p' ) then

  do n_out = 1, nb

    b(n_out) = missing_value
    interp = .false.
    n_in = 1

    do while ( (.not.interp) .and. (n_in < na) )

      if( (xa(n_in)   >= xb(n_out)) .and. &
          (xa(n_in+1) <= xb(n_out))        ) then
        interp = .true.
        w1 = (xa(n_in+1)-xb(n_out))/(xa(n_in+1)-xa(n_in))
        w2 = 1. - w1
        b(n_out) = w1*a(n_in) + w2*a(n_in+1)
      end if
      n_in = n_in +1

    enddo

  enddo
  
  else

  do n_out = 1, nb

    b(n_out) = missing_value
    interp = .false.
    n_in = 1

    do while ( (.not.interp) .and. (n_in < na) )

      if( (xa(n_in)   <= xb(n_out)) .and. &
          (xa(n_in+1) >= xb(n_out))        ) then
        interp = .true.
        w1 = (xa(n_in+1)-xb(n_out))/(xa(n_in+1)-xa(n_in))
        w2 = 1. - w1
        b(n_out) = w1*a(n_in) + w2*a(n_in+1)
      end if
      n_in = n_in +1

    enddo

  enddo
  
  endif

  end subroutine interp_1d

!-------------------------------------------------------------------------
!
! This routines has been taken "as is" from wrf_user_fortran_util_0.f
!
! This routine assumes
!    index order is (i,j,k)
!    wrf staggering
!    units: pressure (Pa), temperature(K), height (m), mixing ratio (kg kg{-1})
!    availability of 3d p, t, and qv; 2d terrain; 1d half-level zeta string
!    output units of SLP are Pa, but you should divide that by 100 for the
!          weather weenies.
!    virtual effects are included
!
! Dave

      subroutine compute_seaprs ( nx , ny , nz  ,         &
                                  z, t , p , q ,          &
                                  sea_level_pressure,debug)
!     &                            t_sea_level, t_surf, level )
      IMPLICIT NONE
!     Estimate sea level pressure.
      INTEGER nx , ny , nz
!f2py integer, optional, depend(z), check(shape(z,0)==nx) :: nx=shape(z,0)
!f2py integer, optional, depend(z), check(shape(z,1)==ny) :: ny=shape(z,1)
!f2py integer, optional, depend(z), check(shape(z,2)==nz) :: nz=shape(z,2)
      REAL    z(nx,ny,nz), t(nx,ny,nz) , p(nx,ny,nz) , q(nx,ny,nz)
!f2py real, intent(in) :: z
!f2py real, intent(in), depend(nx, ny, nz) :: t
!f2py real, intent(in), depend(nx, ny, nz) :: p
!f2py real, intent(in), depend(nx, ny, nz) :: q
!     The output is the 2d sea level pressure.
      REAL    sea_level_pressure(nx,ny)
!f2py real, intent(out), depend(nx,ny) :: sea_level_pressure
      INTEGER level(nx,ny)
      REAL t_surf(nx,ny) , t_sea_level(nx,ny)
      LOGICAL debug
!f2py logical, intent(in) :: debug

!     Some required physical constants:

      REAL R, G, GAMMA
      PARAMETER (R=287.04, G=9.81, GAMMA=0.0065)

!     Specific constants for assumptions made in this routine:

      REAL    TC, PCONST
      PARAMETER (TC=273.16+17.5, PCONST = 10000)
      LOGICAL ridiculous_mm5_test
      PARAMETER (ridiculous_mm5_test = .TRUE.)
!      PARAMETER (ridiculous_mm5_test = .false.)

!     Local variables:

      INTEGER i , j , k
      INTEGER klo , khi


      REAL plo , phi , tlo, thi , zlo , zhi
      REAL p_at_pconst , t_at_pconst , z_at_pconst
      REAL z_half_lowest

      REAL    , PARAMETER :: cp           = 7.*R/2.
      REAL    , PARAMETER :: rcp          = R/cp
      REAL    , PARAMETER :: p1000mb      = 100000.

      LOGICAL  l1 , l2 , l3, found

!     Find least zeta level that is PCONST Pa above the surface.  We later use this
!     level to extrapolate a surface pressure and temperature, which is supposed
!     to reduce the effect of the diurnal heating cycle in the pressure field.

      t(:,:,:) = (t(:,:,:)+300.)*(p(:,:,:)/p1000mb)**rcp

      DO j = 1 , ny
         DO i = 1 , nx
            level(i,j) = -1

            k = 1
            found = .false.
            do while( (.not. found) .and. (k.le.nz))
               IF ( p(i,j,k) .LT. p(i,j,1)-PCONST ) THEN
                  level(i,j) = k
                  found = .true.
               END IF
               k = k+1
            END DO

            IF ( level(i,j) .EQ. -1 ) THEN
            PRINT '(A,I4,A)','Troubles finding level ',   &
                        NINT(PCONST)/100,' above ground.'
            PRINT '(A,I4,A,I4,A)',                        &
                  'Problems first occur at (',i,',',j,')'
            PRINT '(A,F6.1,A)',                           &
                  'Surface pressure = ',p(i,j,1)/100,' hPa.'
            STOP 'Error_in_finding_100_hPa_up'
         END IF


         END DO
      END DO

!     Get temperature PCONST Pa above surface.  Use this to extrapolate
!     the temperature at the surface and down to sea level.

      DO j = 1 , ny
         DO i = 1 , nx

            klo = MAX ( level(i,j) - 1 , 1      )
            khi = MIN ( klo + 1        , nz - 1 )

            IF ( klo .EQ. khi ) THEN
               PRINT '(A)','Trapping levels are weird.'
               PRINT '(A,I3,A,I3,A)','klo = ',klo,', khi = ',khi, &
                            ': and they should not be equal.'
               STOP 'Error_trapping_levels'
            END IF

         plo = p(i,j,klo)
         phi = p(i,j,khi)
         tlo = t(i,j,klo)*(1. + 0.608 * q(i,j,klo) )
         thi = t(i,j,khi)*(1. + 0.608 * q(i,j,khi) )
!         zlo = zetahalf(klo)/ztop*(ztop-terrain(i,j))+terrain(i,j)
!         zhi = zetahalf(khi)/ztop*(ztop-terrain(i,j))+terrain(i,j)
         zlo = z(i,j,klo)
         zhi = z(i,j,khi)

         p_at_pconst = p(i,j,1) - pconst
         t_at_pconst = thi-(thi-tlo)*LOG(p_at_pconst/phi)*LOG(plo/phi)
         z_at_pconst = zhi-(zhi-zlo)*LOG(p_at_pconst/phi)*LOG(plo/phi)

         t_surf(i,j) = t_at_pconst*(p(i,j,1)/p_at_pconst)**(gamma*R/g)
         t_sea_level(i,j) = t_at_pconst+gamma*z_at_pconst

         END DO
      END DO

!     If we follow a traditional computation, there is a correction to the sea level
!     temperature if both the surface and sea level temnperatures are *too* hot.

      IF ( ridiculous_mm5_test ) THEN
         DO j = 1 , ny
            DO i = 1 , nx
               l1 = t_sea_level(i,j) .LT. TC
               l2 = t_surf     (i,j) .LE. TC
               l3 = .NOT. l1
               IF ( l2 .AND. l3 ) THEN
                  t_sea_level(i,j) = TC
               ELSE
                  t_sea_level(i,j) = TC - 0.005*(t_surf(i,j)-TC)**2
               END IF
            END DO
         END DO
      END IF

!     The grand finale: ta da!

      DO j = 1 , ny
      DO i = 1 , nx
!         z_half_lowest=zetahalf(1)/ztop*(ztop-terrain(i,j))+terrain(i,j)
         z_half_lowest=z(i,j,1)
         sea_level_pressure(i,j) = p(i,j,1) *              &
                               EXP((2.*g*z_half_lowest)/   &
                                   (R*(t_sea_level(i,j)+t_surf(i,j))))

         sea_level_pressure(i,j) = sea_level_pressure(i,j)*0.01

      END DO
      END DO

!    if (debug) then
!      print *,'sea pres input at weird location i=20,j=1,k=1'
!      print *,'t=',t(20,1,1),t(20,2,1),t(20,3,1)
!      print *,'z=',z(20,1,1),z(20,2,1),z(20,3,1)
!      print *,'p=',p(20,1,1),p(20,2,1),p(20,3,1)
!      print *,'slp=',sea_level_pressure(20,1),     &
!               sea_level_pressure(20,2),sea_level_pressure(20,3)
!    endif
!      print *,'t=',t(10:15,10:15,1),t(10:15,2,1),t(10:15,3,1)
!      print *,'z=',z(10:15,1,1),z(10:15,2,1),z(10:15,3,1)
!      print *,'p=',p(10:15,1,1),p(10:15,2,1),p(10:15,3,1)
!      print *,'slp=',sea_level_pressure(10:15,10:15),     &
!         sea_level_pressure(10:15,10:15),sea_level_pressure(20,10:15)

      end subroutine compute_seaprs

!------------------------------------------------------------------


  subroutine rotate_wind (u,v,d1,d2,d3,var,                 &
                          map_proj,cen_lon,xlat,xlon,       &
                          truelat1,truelat2,data_out)

  implicit none

  integer, intent(in)            ::  d1, d2, d3

  real, dimension(d1,d2,d3)      :: data_out
  integer                        :: map_proj,i,j,k
  real                           :: cen_lon, truelat1, truelat2, cone
  real, dimension(d1,d2,d3)      :: u,v 
  real, dimension(d1,d2)         :: xlat, xlon, diff, alpha

  character (len=10), intent(in) :: var

   REAL    , PARAMETER           :: pii = 3.14159265
   REAL    , PARAMETER           :: radians_per_degree = pii/180.




       cone = 1.                                          !  PS
       if( map_proj .eq. 1) then                          !  Lambert Conformal mapping
         IF (ABS(truelat1-truelat2) .GT. 0.1) THEN
            cone=(ALOG(COS(truelat1*radians_per_degree))-            &
                  ALOG(COS(truelat2*radians_per_degree))) /          &
            (ALOG(TAN((90.-ABS(truelat1))*radians_per_degree*0.5 ))- &
             ALOG(TAN((90.-ABS(truelat2))*radians_per_degree*0.5 )) )
         ELSE
            cone = SIN(ABS(truelat1)*radians_per_degree )  
         ENDIF
       end if


       diff = xlon - cen_lon

       do i = 1, d1
       do j = 1, d2
         if(diff(i,j) .gt. 180.) then
           diff(i,j) = diff(i,j) - 360.
         end if
         if(diff(i,j) .lt. -180.) then
           diff(i,j) = diff(i,j) + 360.
         end if
       end do
       end do

       do i = 1, d1
       do j = 1, d2
          if(xlat(i,j) .lt. 0.) then
            alpha(i,j) = - diff(i,j) * cone * radians_per_degree
          else
            alpha(i,j) = diff(i,j) * cone * radians_per_degree
          end if
       end do
       end do


       if(var(1:1) .eq. "U") then
         do k=1,d3
           data_out(:,:,k) = v(:,:,k)*sin(alpha) + u(:,:,k)*cos(alpha)
         end do
       else if(var(1:1) .eq. "V") then    
         do k=1,d3
           data_out(:,:,k) = v(:,:,k)*cos(alpha) - u(:,:,k)*sin(alpha)
         end do
       end if


  end subroutine rotate_wind

END MODULE from_wrf_to_grads
