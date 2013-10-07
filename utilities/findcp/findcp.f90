module opttools 
 implicit none
 contains


!---print out geometry info ---
subroutine analysegeom(natoms,geom)
  implicit none
  integer, intent(in)   :: natoms
  double precision,intent(in)  ::  geom(3,natoms)
  double precision, parameter  ::  bohr2ang=0.529177249d0
  integer   ::  i,j,k,l
  double precision  ::  distmat(natoms,natoms), TLen = 2.2D0, d,d1(3),d2(3),d3(3),cpd(3)
  double precision,external  ::   dnrm2
  print *,"Cartesian Geometries in Atomic Units"
  print "(2x,3F15.10)",geom
  distmat = 0d0
  print "(/,A)","   Atom1   Atom2   Bond Length(Ang)"
  do i=1,natoms
    do j=i+1,natoms
      d = dnrm2(3,geom(:,i)-geom(:,j),1)*bohr2ang
      distmat(i,j) = d
      distmat(j,i) = d
      if(d<TLen)print "(2x,I5,3x,I5,3x,F12.8)",i,j,d
    end do
  end do

  print "(/,A)","   Atom1   Atom2   Atom3   Bond Angle (Degrees)"
  do i=1,natoms
    do j=1,natoms
     if(j==i .or. distmat(i,j)>TLen)cycle
     d1 = geom(:,j)-geom(:,i)
     d1 = d1/dnrm2(3,d1,1)
     do k=j+1,natoms
       if(k==i .or. distmat(i,k)>TLen)cycle
       d2 = geom(:,k)-geom(:,i)
       d2 = d2/dnrm2(3,d2,1)
       print "(2x,3(I5,3x),F12.4)",J,I,K, 90/Acos(0d0)* &
           ACOS(dot_product(d1,d2))
     end do!k
    end do !j
  end do   !i    

  print "(/,A)","   Atom1   Atom2   Atom3   Atom4   Out-of-Plane Angle (Degrees)"
  do i=1,natoms
    do j=1,natoms
     if(j==i .or. distmat(i,j)>TLen)cycle
     d1 = geom(:,j)-geom(:,i)
     d1 = d1/dnrm2(3,d1,1)
     do k=j+1,natoms
       if(k==i .or. distmat(i,k)>TLen)cycle
       d2 = geom(:,k)-geom(:,i)
       d2 = d2/dnrm2(3,d2,1)
       do l=k+1,natoms
         if(l==i .or. distmat(i,l)>TLen)cycle
         d3 = geom(:,l)-geom(:,i)
         d3 = d3/dnrm2(3,d3,1)
         cpd(1) = d1(3)*(d2(2)-d3(2))+d2(3)*d3(2)-d2(2)*d3(3)+d1(2)*(d3(3)-d2(3))
         cpd(2) =-d2(3)*d3(1)+d1(3)*(d3(1)-d2(1))+d1(1)*(d2(3)-d3(3)) +d2(1)*d3(3)
         cpd(3) = d1(2)*(d2(1)-d3(1))+d2(2)*d3(1)-d2(1)*d3(2)+d1(1)*(d3(2)-d2(2))
         print "(2x,4(I5,3x),F12.4)",I,J,K,L, 90/Acos(0d0)* &
           asin((-d1(3)*d2(2)*d3(1)+d1(2)*d2(3)*d3(1)+d1(3)*d2(1)*d3(2)       &
                -d1(1)*d2(3)*d3(2)-d1(2)*d2(1)*d3(3)+d1(1)*d2(2)*d3(3))/      &
               dnrm2(3,cpd,1))   
       end do! l
     end do!k
    end do !j
  end do   !i    
end subroutine analysegeom


!---calculate harmonic frequencies from hessian matrix
subroutine getFreq(natoms,masses,hess,w)
  implicit none
  integer,intent(in)          :: natoms
  double precision,intent(in) :: masses(natoms),hess(3*natoms,3*natoms)
  double precision,intent(out):: w(3*natoms)
  double precision,  parameter  :: amu2au=1.822888484514D3,au2cm1=219474.6305d0
  integer  :: i,j
  double precision  :: sqrm,  hmw(3*natoms,3*natoms), tmp(1)
  double precision,dimension(:),allocatable :: WORK
  integer,dimension(:),allocatable :: IWORK
  integer           :: LIWORK, LWORK, itmp(1),INFO
  ! convert hessian into mass weighed coordinates
  hmw = hess/amu2au
  do i=1,natoms
    sqrm = sqrt(masses(i))
    do j=i*3-2,i*3
      hmw(j,:)=hmw(j,:)/sqrm
      hmw(:,j)=hmw(:,j)/sqrm
    end do
  end do 
  !calculate eigenvalues of hmw
  call DSYEVD('N','U',3*natoms,hmw,3*natoms,w,tmp,-1,itmp,-1,INFO)
  if(info/=0)print *,"DSYEVD allocation investigation failed.  info=",info
  LWORK = int(tmp(1))
  LIWORK= itmp(1)
  allocate(WORK(LWORK))
  allocate(IWORK(LIWORK))
  
  call DSYEVD('N','U',3*natoms,hmw,3*natoms,w,WORK,LWORK,IWORK,LIWORK,INFO)
  if(info/=0)print *,"DSYEVD failed.  info=",info
  do i=1,3*natoms
    if(w(i)<0)then
      w(i) = -sqrt(-w(i))*au2cm1
    else 
      w(i) = sqrt(w(i))*au2cm1
    end if
  end do
end subroutine getFreq


!---calculate hessian at a certain geometry
subroutine calcHess(natoms,cgeom,nstate,istate,stepsize,hessian,skip)
  implicit none
  integer, intent(in)          :: natoms, nstate,istate
  logical, intent(in),optional :: skip(natoms)
  double precision,intent(in)  :: stepsize
  double precision,intent(in)  :: cgeom(3*natoms)
  double precision,intent(out) :: hessian(3*natoms,3*natoms)
  double precision   ::  mdif
  
  integer   ::   i,  j
  logical   ::   skipdisp(natoms*3)
  double precision  ::  dispgeom(3*natoms), dgrd(3*natoms)
  real*8    ::  h(nstate,nstate),cg(3*natoms,nstate,nstate),dcg(3*natoms,nstate,nstate),e(nstate)
  skipdisp=.false.
  if(present(skip))then
    do i=1,natoms
      if(skip(i))skipdisp(i*3-2:i*3)=.true.
    end do
  end if
  hessian = 0d0
  do i=1,3*natoms
    if(skipdisp(i))cycle
    dispgeom=cgeom
    dispgeom(i)=dispgeom(i) - stepsize
    call EvaluateSurfgen(dispgeom,e,cg,h,dcg)
    dgrd  =-cg(:,istate,istate)
    dispgeom=cgeom
    dispgeom(i)=dispgeom(i) + stepsize
    call EvaluateSurfgen(dispgeom,e,cg,h,dcg)
    dgrd = dgrd+cg(:,istate,istate)
    hessian(i,:)= dgrd/2/stepsize
  end do!o=1,3*natoms
  do i=1,3*natoms
    if(skipdisp(i))hessian(:,i)=0d0
  end do
  mdif = maxval(abs(hessian-transpose(hessian)))
  if(mdif>1d-5)print *,"maximum hermiticity breaking : ",mdif
  hessian = (hessian+transpose(hessian))/2
end subroutine calcHess


!----search for minimum on adiabatic surfaces 
subroutine findmin(natoms,nstate,cgeom,isurf,maxiter,shift,Etol,Stol)
  implicit none
  integer, intent(in)                                 ::  natoms,isurf,maxiter,nstate
  double precision,dimension(3*natoms),intent(inout)  ::  cgeom
  double precision,intent(in)                         ::  shift,Etol,Stol

  real*8    ::  h(nstate,nstate),cg(3*natoms,nstate,nstate),dcg(3*natoms,nstate,nstate),e(nstate)
  double precision,dimension(3*natoms)  ::  grad, b1, b2, w
  double precision,dimension(3*natoms,3*natoms)  :: hess,hinv
  double precision,dimension(:),allocatable :: WORK
  integer,dimension(:),allocatable :: IWORK
  integer           :: LIWORK, LWORK, itmp(1),INFO  
  integer  ::  iter  , i
  double precision            :: nrmG, nrmD, tmp(1)
  double precision, external  :: dnrm2
  double precision,  parameter  :: amu2au=1.822888484514D3,au2cm1=219474.6305d0
  double precision, parameter   :: MAXD = 1D-1

  ! initialize work spaces
  call DSYEVD('V','U',3*natoms,hess,3*natoms,w,tmp,-1,itmp,-1,INFO)
  if(info/=0)print *,"DSYEVD allocation failed.  info=",info
  LWORK = int(tmp(1))
  LIWORK= itmp(1)
  allocate(WORK(LWORK))
  allocate(IWORK(LIWORK))

  print "(A,I4,A,I4,A)","Searching for minimum on surface ",isurf," in ",maxiter," iterations."
  print "(A)","  Convergence Tolerances"
  print "(A,E10.2,A,E10.2)","  Energy Gradient: ",Etol,"   Displacement:",Stol
  do iter=1,maxiter
     call EvaluateSurfgen(cgeom,e,cg,h,dcg)
     grad = cg(:,isurf,isurf)
     nrmG=dnrm2(3*natoms,grad,1)
     call calcHess(natoms,cgeom,nstate,isurf,1D-4,hess)
     hinv = hess
     ! invert the hessian
     call DSYEVD('V','U',3*natoms,hess,3*natoms,w,WORK,LWORK,IWORK,LIWORK,INFO)
     if(info/=0)print *,"DSYEVD failed.  info=",info
     ! hess. w. hess^T = hess_old
     ! so x= hess_old^-1.b = hess. w^-1. hess^T. b
     ! b1= hess^T.g       =>
     call DGEMV('T',3*natoms,3*natoms,1d0,hess,3*natoms,grad,1,0d0,b1,1)
     ! b1' = w^-1*b1
     do i=1,3*natoms
      if(abs(w(i))<shift)then
         b1(i)=dble(0)
      else!
         b1(i)=b1(i)/w(i)
      end if!
     end do!
     ! b2 = hess.b1'
     call DGEMV('N',3*natoms,3*natoms,1d0,hess,3*natoms,b1,1,0d0,b2,1)
     ! cgeom' = cgeom - H^-1.g
     nrmD=dnrm2(3*natoms,b2,1)
     if(nrmD>maxD)then
       b2=b2/nrmD*maxD
       nrmD=maxD
     end if
     print 1000,iter,e(isurf)*au2cm1,nrmG,nrmD
     cgeom = cgeom - b2
     if(nrmG<Etol.and.nrmD<Stol)then
       print *,"Optimization converged"
       return
     end if
  end do 
1000 format("Iteration ",I4,": E=",F20.4,", |Grd|=",E12.5,", |Disp|=",E12.5)
end subroutine findmin
end module opttools 

! problems:  wrong bound for array r() in subroutine 
! rx1 and rx2 have single precision default implicit type
! variable in subroutine set as intent(INOUT) as default and r(6) got switched every run
! inappropriate bound for grid

program findcp

 use opttools 
 implicit none

  integer :: natm , nst, ngeoms
  integer :: i,j, isurf, ios
  double precision,allocatable  :: anum(:),masses(:),  hess(:,:), w(:)
  character*3,allocatable       :: aname(:)
  double precision,allocatable  :: cgeom(:),hvec(:)
  logical,allocatable           :: skip(:)
  double precision,allocatable  :: h(:,:),cg(:,:,:),dcg(:,:,:),e(:),evec(:,:)
  double precision,external     :: dnrm2
  character*300                 ::  geomfile,str


  print *," ***************************************** "
  print *," *    findcp.x                           * "
  print *," ***************************************** "
  print *," Critical point search on Surfgen surface" 
  call initPotential()
  call getInfo(natm,nst)

! allocate arrays

  print "(A,I6)","Number of Atoms:  ",natm
  print "(A,I6)","Number of States: ",nst 

  allocate(masses(natm))
  allocate(anum(natm))
  allocate(aname(natm))
  allocate(hess(natm*3,natm*3))
  allocate(w(3*natm))
  allocate(cgeom(3*natm))
  allocate(hvec(3*natm))
  allocate(skip(natm))
  skip = .false.
  allocate(h(nst,nst))
  allocate(cg(3*natm,nst,nst))
  allocate(dcg(3*natm,nst,nst))
  allocate(e(nst))
  allocate(evec(nst,nst)) 

! process arguments
! synposis:    findcp.x geom type states
! Default values:
! geom        geom
! type        min
! states      1
  call get_command_argument(number=1,value=geomfile,status=ios)
  if(ios<-1)  &     ! input argument larger than container
      stop "Filename truncated.  Please use less than 300 characters"
  if(ios>0) then
    print *,"No filename suppied.  Using default."
    write(geomfile,"(A)"), "geom"
    isurf = 1
  else
    call get_command_argument(number=2,value=str,status=ios)
    if(ios/=0)then
      print *,"Cannot get surface number from command line options.  Using default."
      isurf = 1
    else
      read(str,*,iostat=ios)isurf
      if(ios/=0)then
        print *,"Cannot get surface number from command line options.  Using default."
        isurf = 1
      end if
    end if
  end if

  print *,"Reading input from input file "//trim(geomfile)
  call readColGeom(geomfile,1,natm,aname,anum,cgeom,masses)

  print "(/,A)","-------------- Initial Geometry ------------"
  call analysegeom(natm,cgeom)
  print "(/,A)","----------- Geometry Optimizations ---------"
  call findmin(natm,nst,cgeom,isurf,100,1d-3,1d-7 ,1d-5)
  print "(/,A)","--------------- Final Geometry -------------"
  call analysegeom(natm,cgeom)
  print "(/,A)","------------ Harmonic Frequencies ----------"
  call calcHess(natm,cgeom,nst,isurf,1D-3,hess,skip)
  call getFreq(natm,masses,hess,w)
  do i=1,3*natm
    print "(I5,F12.2)",i,w(i)
  end do

! deallocate arrays
  deallocate(masses)
  deallocate(anum)
  deallocate(aname)
  deallocate(hess)
  deallocate(w)
  deallocate(cgeom)
  deallocate(hvec)
  deallocate(skip)
  deallocate(h)
  deallocate(cg)
  deallocate(dcg)
  deallocate(e)
  deallocate(evec) 

end