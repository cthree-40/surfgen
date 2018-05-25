program getneighbor
  ! Get closest point and distance to point for a geometry.
  use potdata
  use hddata, only: ncoord
  implicit none
  integer :: num_atoms, num_states
  double precision, dimension(:), allocatable :: input_geometry
  
  character(255) :: ingeom_file
  integer :: ios
  integer :: neighbor
  double precision, dimension(:), allocatable :: eerr
  double precision :: dist

  double precision, dimension(:,:),   allocatable :: bmat
  double precision, dimension(:),     allocatable :: igeom
  double precision, dimension(:),     allocatable :: energy
  double precision, dimension(:,:,:), allocatable :: cgrads, dcgrads
  double precision, dimension(:,:),   allocatable :: hmat

  integer :: nc, j
  double precision, external :: dnrm2

  call initpotential()
  call getinfo(num_atoms, num_states)
  allocate(input_geometry(num_atoms*3))
  allocate(eerr(num_states))
  allocate(bmat(ncoord, 3*num_atoms))
  allocate(igeom(ncoord))
  allocate(hmat(num_states, num_states))
  allocate(cgrads(num_atoms*3,num_states, num_states))
  allocate(dcgrads(num_atoms*3,num_states,num_states))
  allocate(energy(num_states))
  call get_command_argument(1, value = ingeom_file, status = ios)
  if (ios .ne. 0) stop "Error reading geometry file name."
  print *, "File: ", trim(adjustl(ingeom_file))

  call read_geometry_from_file(ingeom_file, num_atoms, input_geometry)
  print *, "Input geometry: "
  print "(3f14.8)", input_geometry
  call EvaluateSurfgen(input_geometry, energy, cgrads, hmat, dcgrads)
  call buildWBMat(input_geometry, igeom, bmat)
  call getmind(igeom, dist, neighbor, eerr)
  if (cvanish <= 0) then
          nc = ndcoord
  else
          nc=0
          do j=1, ndcoord
                  if (dnrm2(3*num_atoms,bmat(dcoordls(j),1),ncoord)>cvanish)nc=nc+1
          end do
  end if
  print *, dist
  dist = dist / sqrt(dble(nc))
  print "(f14.2)", energy*219474.63
  print *, "Neighbor: ", neighbor
  print *, "Distance: ", dist
                  
contains
  subroutine read_geometry_from_file(gmfl, na, gm)
    implicit none
    character(255), intent(in) :: gmfl
    integer, intent(in) :: na
    double precision, dimension(3, na), intent(inout) :: gm
    integer :: i, ios
    
    open(file = trim(adjustl(gmfl)), unit = 15, status = "old", &
            action = "read", iostat = ios)
    if (ios .ne. 0) stop "Error reading geometry from file."
    do i = 1, na
            read(unit = 15, fmt = "(10x,f14.8,f14.8,f14.8,14x)") &
                    gm(1,i), gm(2,i), gm(3,i)
    end do
    close(unit = 15)
    return
  end subroutine read_geometry_from_file
end program getneighbor
