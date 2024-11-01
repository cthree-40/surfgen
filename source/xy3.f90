! This module contains routine necessary to evaluate an XY3 system, treated
! with the G4 (two equivalent atoms) CNPI group in surfgen.
module xy3
  implicit none
  
  ! Six distances are evaluated: Y12, Y13, and Y23; and XY1, XY2, and XY3
  integer, dimension(3,3)          :: ypairs
  double precision, dimension(3)   :: ydists
  double precision, dimension(3)   :: xdists
  
  ! Tolerance to use second distance criterion (XYa distances)
  double precision, parameter :: xdtol = 1d-1
  
  ! Ordered geometry
  double precision, dimension(3,4) :: canon_geom

  ! Mapping
  integer, dimension(4)            :: ymap

  ! For trajectory-use: check which "lobe" we are in
  logical :: checklobe
  integer :: lastlobe
  
contains
  !*
  ! compute_yatomdists: compute Y atom distances
  !*
  subroutine compute_xyatomdists(cgeom, natom)
    implicit none
    integer, intent(in) :: natom
    double precision, dimension(3, natom), intent(in) :: cgeom
    integer :: i
    ydists = 0d0
    do i = 1, 3
            ydists(i) = distance_3d(cgeom(:,ypairs(1,i)), cgeom(:,ypairs(2,i)))
            xdists(i) = distance_3d(cgeom(:,1), cgeom(:,ypairs(3,i)))
    end do
    return
  end subroutine compute_xyatomdists
  !*
  ! distance_3d: compute distance between two points in 3D
  !*
  function distance_3d(pt1, pt2) result(dist)
    implicit none
    double precision, dimension(3), intent(in) :: pt1, pt2
    double precision, dimension(3) :: d
    double precision :: dist
    double precision, external :: dnrm2
    d = pt2 - pt1
    dist = dnrm2(3,d,1)
    return
  end function distance_3d
  !*
  ! find_min_array: find minimum of array
  !*
  function find_min_array(array, len) result(index)
    implicit none
    integer, intent(in) :: len
    double precision, dimension(len), intent(in) :: array
    integer :: i, index
    index = 1
    do i = 2, len
      if (array(i) .lt. array(index)) then
        index = i
      end if
    end do
    return
  end function find_min_array
  !*
  ! init_xy3: initialize xy3 module values
  !*
  subroutine init_xy3()
    implicit none
    ypairs(:,1) = (/ 2, 3, 4 /)
    ypairs(:,2) = (/ 2, 4, 3 /)
    ypairs(:,3) = (/ 4, 3, 2 /)

    checklobe = .true.
    lastlobe  = 1
    
    return
  end subroutine init_xy3

  !*
  ! make_ymap: generate mapping for Y atoms
  !*
  subroutine make_ymap (cgeom, na)
    implicit none
    integer, intent(in) :: na
    double precision, dimension(3,na), intent(in) :: cgeom
    integer :: index
    call compute_xyatomdists(cgeom, na)
    index = find_min_array(ydists, 3)
    write(*,"('YDISTS: ',3f10.5)") ydists
    write(*,"('XDISTS: ',3f10.5)") xdists
    ! Check if all distances are within XDTOL
    if (    abs(ydists(1)-ydists(2)) .le. xdtol .and. &
            abs(ydists(1)-ydists(3)) .le. xdtol .and. &
            abs(ydists(2)-ydists(3)) .le. xdtol ) then
            index = find_min_array(xdists, 3)
    end if
    ymap(1) = 1
    ymap(2) = ypairs(1,index)
    ymap(3) = ypairs(2,index)
    ymap(4) = ypairs(3,index)
    if (checklobe) then
       if (lastlobe .ne. index) then
          print "(A)", "WARNING: Lobe switch!"
          lastlobe = index
       end if
    end if
    return
  end subroutine make_ymap

end module xy3
