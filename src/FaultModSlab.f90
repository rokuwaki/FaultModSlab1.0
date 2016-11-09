program FaultModSlab
  use geodesic
  use, intrinsic:: iso_fortran_env, only: ERROR_UNIT
  implicit none
  integer                   :: k, nk, l, nl, k0, l0
  integer, parameter        :: omask = 0
  real                      :: num_dip, num_stk, x, y, xk, yk, stk, stk1, lonp, latp, w0
  real                      :: depth, strike, dip, rake, dk,dl
  real, parameter           :: pi=3.1415927, pideg=pi/180.0
  real(kind=8), allocatable :: lons(:), lats(:), depths(:), strikes(:), dips(:)
  real(kind=8), parameter   :: a = 6378137d0, f = 1/298.257223563d0
  real(kind=8)              :: azi1, azi2, s12a12, dummy1, dummy2, dummy3, dummy4, dummy5
  real(kind=8)              :: lat0, lon0, lat, lon
  logical                   :: isDepthInterpolated, isStrikeInterpolated, isDipInterpolated
  character(len=2**10)      :: depth_xyz
  character(len=2**10)      :: strike_xyz
  character(len=2**10)      :: dip_xyz
  
  open(11, file = 'work/input.dat',status = 'old', action = 'read')
  read(11, *) lat0, lon0  
  read(11, *) latp, lonp, w0
  read(11, *) dl, dk, nl, nk, l0, k0, stk
  stk = stk * pideg
  read(11, *) depth_xyz
  read(11, *) strike_xyz
  read(11, *) dip_xyz
  open(14, file = 'work/knot_value.dat', status = 'replace', action = 'write')
  open(15, file = 'work/i_greenf', status = 'replace', action = 'write')
  open(16, file = 'work/epicenter.dat', status = 'replace', action = 'write')

  call load_slab1_file(trim(depth_xyz), lons, lats, depths)
  write(6, *) 'Loaded slab1.0 depth data'
  call load_slab1_file(trim(strike_xyz), lons, lats, strikes)
  write(6, *) 'Loaded slab1.0 strike data'
  call load_slab1_file(trim(dip_xyz), lons, lats, dips)
  write(6, *) 'Loaded slab1.0 dip data'

  do k = 1, nk !! along-dip
     num_dip =  k - k0
     xk = 0.0 - dk * 1000 * num_dip * cos(stk)
     yk = 0.0 + dk * 1000 * num_dip * sin(stk)
     do l = 1, nl
        num_stk = l - l0 !! along-strike
        x = xk + dl * 1000 * num_stk * sin(stk)
        y = yk + dl * 1000 * num_stk * cos(stk)
        s12a12 = sqrt(x**2 + y**2)
        call CalAzi(x, y, azi1)
        call direct(a, f, lat0, lon0, azi1, s12a12, .False., &
             lat, lon, azi2, omask, dummy1, dummy2, dummy3, dummy4 , dummy5)
        call interpolate_from_lons_lats_values(lon, lat, lons, lats, depths, depth, isDepthInterpolated)
        call interpolate_from_lons_lats_values(lon, lat, lons, lats, strikes, strike, isStrikeInterpolated)
        call interpolate_from_lons_lats_values(lon, lat, lons, lats, dips, dip, isDipInterpolated)
        call cal_rake(lat, lon, latp, lonp, w0, strike, rake)
        if(isDepthInterpolated .and. isStrikeInterpolated .and. isDipInterpolated)then
           write(14, '(2i4, 3f10.3, 3f7.1, l2)') l, k, lat, lon, abs(depth), strike, abs(dip), rake, .true.
        else
           write(14, '(2i4, 3f10.3, 3f7.1, l2)') l, k, lat, lon, abs(depth), strike, abs(dip), rake, .false.
        end if
        if ((l .eq. l0) .and. (k .eq. k0)) then
           write(15, '(a4, a5, 3f6.1, f8.3, a4)') '2048', '0.05', strike, abs(dip), rake, abs(depth), '2.0'
           write(15, '(2f8.1, 4i5, 2a3)') dl, dk, nl, nk, l0, k0, '1', '45'
           write(15, '(a10)') 'ak135 P'
           write(16, '(6f10.3)') lat, lon, abs(depth), strike, abs(dip), rake
        end if
     end do
  end do
contains

  function is_nan(x) result(answer)
    Real, intent(in):: x
    Logical:: answer

    answer = (x /= x)
  end function is_nan

  subroutine load_slab1_file(path, lons, lats, values)
    Character(len=*), intent(in):: path
    Real(kind=8), allocatable, intent(inout):: lons(:), lats(:), values(:)
    Integer:: io
    Integer:: ios
    Integer:: n_lines
    Real:: lon, lat, value
    Integer:: i_line
    Real:: rTrash, nanTest

    open(newunit=io, file=trim(path), status='old', action='read')
    n_lines = 0
    do
       read(io, *, iostat=ios) rTrash, rTrash, nanTest
       if(ios == 0 .and. (.not. is_nan(nanTest)))then
          n_lines = n_lines + 1
       else if(ios < 0)then
          exit
       end if
    end do

    if(allocated(lons)) deallocate(lons)
    if(allocated(lats)) deallocate(lats)
    if(allocated(values)) deallocate(values)
    allocate(lons(1:n_lines))
    allocate(lats(1:n_lines))
    allocate(values(1:n_lines))
    
    rewind(io)
    i_line = 0
    do
       if(i_line >= n_lines)exit

       read(io, *, iostat = ios) lon, lat, value
       if(ios == 0 .and. (.not. is_nan(value)))then
          i_line = i_line + 1
          if(i_line > n_lines)then
             write(ERROR_UNIT, *) 'MUST NOT HAPPEN: i_line > n_line'
             stop 1
          end if
          lons(i_line) = lon
          lats(i_line) = lat
          values(i_line) = value
       else if(ios < 0)then
          exit
       end if
    end do

    close(io)
  end subroutine load_slab1_file


  subroutine  interpolate_from_lons_lats_values(lon, lat, lons, lats, values, value, isValueInterpolated)
    Real(kind=8), intent(in):: lon, lat, lons(:), lats(:), values(:)
    Real, intent(out):: value
    Logical, intent(out):: isValueInterpolated
    Real:: latA, lonA, valueA
    Real:: latB, lonB, valueB
    Real:: latC, lonC, valueC
    Real:: latD, lonD, valueD
    Logical:: isAExist, isBExist, isCExist, isDExist
    Integer:: i

    if(.not.(size(lons) == size(lats) .and. size(lats) == size(values)))then
       write(ERROR_UNIT, *) 'lons, lats and values should have same sizes'
       stop 1
    end if

    isAExist = .false.
    isBExist = .false.
    isCExist = .false.
    isDExist = .false.
    do i = 1, size(values)
       if(isAExist .and. isBExist .and. isCExist .and. isDExist) exit

       !A
       if (((lon - 0.02 < lons(i)) .and. (lons(i) <= lon)) &
            .and. ((lat - 0.02 .lt. lats(i)) .and. (lats(i) .le. lat))) then 
          latA = lats(i)
          lonA = lons(i)
          valueA = values(i)
          isAExist = .true.
          !B
       else if  (((lon < lons(i)) .and. (lons(i) <= lon + 0.02)) &
            .and. ((lat - 0.02 < lats(i)) .and. (lats(i) <= lat))) then
          latB = lats(i)
          lonB = lons(i)
          valueB = values(i)
          isBExist = .true.

          !C
       else if  (((lon - 0.02 < lons(i)) .and. (lons(i) <= lon)) &
            .and. ((lat < lats(i)) .and. (lats(i) <= lat + 0.02))) then
          latC = lats(i)
          lonC = lons(i)
          valueC = values(i)
          isCExist = .true.

          !D
       else if  (((lon < lons(i)) .and. (lons(i) <= lon + 0.02)) &
            .and. ((lat < lats(i)) .and. (lats(i) <= lat + 0.02))) then
          latD = lats(i)
          lonD = lons(i)
          valueD = values(i)
          isDExist = .true.
       end if
    end do

    isValueInterpolated = .false.
    if(isAExist .and. isBExist .and. isCExist .and. isDExist)then
       value = linear_interpolate(lat, lon, latA, lonA, valueA, latB, &
            lonB, valueB, latC, lonC, valueC, latD, lonD, valueD)
       isValueInterpolated = .true.
    end if
  end subroutine interpolate_from_lons_lats_values


  subroutine cal_rake(lat, lon, latp, lonp, w0, strike, rake)
    implicit none
    real(kind=8):: lat, lon
    real w, w0, wx, wy, wz, rlat, rlon, latp, lonp, rlatp, rlonp, px, py, pz, vx, vy, vz, vn, ve, vd, v, azi, azi0
    real, parameter::r=6370.8e6 ! radious in millimeters
    real, intent(in)::strike
    real, intent(out)::rake
    rlat = lat * pideg
    rlon = lon * pideg
    rlatp = latp * pideg
    rlonp = lonp * pideg
    w = w0 * 1e-06 * pideg ! degree/Myr ---> degree/year
    
    
    wx = w * cos(rlatp) * cos(rlonp)
    wy = w * cos(rlatp) * sin(rlonp)
    wz = w * sin(rlatp)
    
    px = cos(rlat) * cos(rlon)
    py = cos(rlat) * sin(rlon)
    pz = sin(rlat)
    
    vx = r * ( wy * pz - py * wz )
    vy = r * ( wz * px - pz * wx )
    vz = r * ( wx * py - px * wy )
    
    vn = vx * (-sin(rlat) * cos(rlon)) + vy * (-sin(rlat) * sin(rlon)) + vz * cos(rlat)
    ve = vx * (-sin(rlon)) + vy * cos(rlon)
!    vd = vx * (-cos(rlat) * cos(rlon)) + vy * (-cos(rlat) * sin(rlon)) - vz * sin(rlat)
    
!    v = sqrt( vn**2 + ve**2 + vd**2 )
    
    azi0 = atan( ve / vn ) * 1/pideg
    if (azi0 .lt. 0) then 
       azi = 360 - abs(azi0)
       rake = 180 - ( azi - strike)
       if (rake > 180) then
          rake = rake - 360
       end if
    else
       rake = 180 - ( azi0 - strike )
       if (rake > 180) then
          rake = rake - 360
       end if
    end if
  end subroutine cal_rake
  
  
  function linear_interpolate(lat, lon, latA, lonA, valueA, latB, lonB, &
       valueB, latC, lonC, valueC, latD, lonD, valueD) result(value)
    Real(kind=8), intent(in):: lat, lon
    Real, intent(in):: latA, lonA, valueA
    Real, intent(in):: latB, lonB, valueB
    Real, intent(in):: latC, lonC, valueC
    Real, intent(in):: latD, lonD, valueD
    Real:: p, q
    Real:: value

    p = (lon - lonA) / (lonB - lonA)
    q = (lat - latA) / (latC - latA)

    value = (1 - p)*(1 - q)*valueA + p*(1 - q)*valueB + q*(1 - p)*valueC + p*q*valueD
  end function linear_interpolate

  subroutine CalAzi(x, y, azi1)
    implicit none
    real, intent(in)         :: x, y
    real(kind=8), intent(out):: azi1
    real, parameter          :: pi = 3.141519265, pideg = pi/180.0
    
    if (x >= 0 .and. y >= 0) then
       azi1 = pi/2.0 - abs(atan(y/x))
    else if (x < 0 .and. y >= 0) then
       azi1 = 2.0*pi - (pi/2.0 - abs(atan(y/x)))
    else if (x < 0 .and. y < 0) then
       azi1 = 2.0*pi - (pi/2.0 + abs(atan(y/x)))
    else if (x >= 0 .and. y < 0) then
       azi1 = pi/2.0 + abs(atan(y/x))
    end if
    azi1 = azi1 / pideg
    if (x == 0 .and. y == 0) azi1 = 0
    return
  end subroutine CalAzi

end program FaultModSlab
