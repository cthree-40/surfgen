! This module contains routine necessary to evaluate an XY3 system, treated
! with the G4 (two equivalent atoms) CNPI group in surfgen.
module xy3
  implicit none

  ! Three distances are evaluated: Y12, Y13, and Y23
  integer, dimension(3,3)          :: ypairs
  double precision, dimension(3)   :: ydists
  ! Ordered geometry
  double precision, dimension(3,4) :: canon_geom
  ! Mapping
  integer, dimension(4)            :: ymap

contains
  !*
  ! compute_yatomdists: compute Y atom distances
  !*
  subroutine compute_yatomdists(cgeom, natom)
    implicit none
    integer, intent(in) :: natom
    double precision, dimension(3, natom), intent(in) :: cgeom
    integer :: i
    ydists = 0d0
    do i = 1, 3
      ydists(i) = distance_3d(cgeom(:,ypairs(1,i)), cgeom(:,ypairs(2,i)))
    end do
    return
  end subroutine compute_yatomdists
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
    call compute_yatomdists(cgeom, na)
    index = find_min_array(ydists, 3)
    ymap(1) = 1
    ymap(2) = ypairs(1,index)
    ymap(3) = ypairs(2,index)
    ymap(4) = ypairs(3,index)
    return
  end subroutine make_ymap

end module xy3
