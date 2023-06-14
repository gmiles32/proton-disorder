PROGRAM iceIh_disorder
  ! generates a cubic box of iceIh with correct proton disorder 
  ! and vanishing dipole moment
  IMPLICIT NONE

  TYPE half_link1
     INTEGER :: neighbour
     LOGICAL :: bond
  END TYPE half_link1

  TYPE oxy
     REAL(8), DIMENSION (3) :: r ! Coordinate
     TYPE (half_link1), DIMENSION (4) :: part ! The bonds, so each is a link object
     INTEGER :: nbonds ! The number of protons that it is associated with
  END TYPE oxy

  TYPE half_link2
     INTEGER, POINTER :: nbonds
     LOGICAL, POINTER :: bond
  END TYPE half_link2

  TYPE lin
     TYPE (half_link2), DIMENSION (2) :: part
  END TYPE lin

  TYPE at
     REAL(8), DIMENSION (3) :: r
     INTEGER :: nummol
     CHARACTER (len=2) :: alabel
     CHARACTER (len=5) :: alabel_gro
  END TYPE at

  TYPE cell_type
     REAL ( 8 ), DIMENSION ( 3, 3 ) :: hmat
     REAL ( 8 ), DIMENSION ( 3, 3 ) :: h_inv
     REAL ( 8 ) :: deth
  END TYPE cell_type

  ! global variables
  TYPE(oxy), ALLOCATABLE, TARGET, DIMENSION (:) :: oxygen
  TYPE(oxy), ALLOCATABLE, TARGET, DIMENSION (:) :: oxygen_old
  TYPE(lin), ALLOCATABLE, DIMENSION (:) :: link
  TYPE(at), ALLOCATABLE, DIMENSION (:) :: atom
  !
  REAL(8) :: dipole, dipole_old
  TYPE ( cell_type ) :: box
  REAL(8), PARAMETER :: pi = 3.14159265358979d0
  REAL(8), PARAMETER  :: dipole_small=1.d-6
  INTEGER :: nmols, nshake
!!$  INTEGER :: seed(47)

!!$  OPEN(16,file='seed')
!!$  READ(16,*)seed
!!$  CLOSE(16)

  !**************************START******************************************
  ! initialize rng
  CALL RANDOM_SEED !(PUT = seed(:))        
  ! generate a hexagonal diamond structure
  CALL lattice
  ! find the nearest neighbours
  CALL neighbour_list
  ! initialize the hydrogens randomly
  CALL init_hydrogens
  ! generate first correct proton disorder by an MC procedure
  CALL mc
  CALL get_dipole(dipole)
  WRITE(33,*)'dipole moment = ', dipole
  nshake = nmols/10
  ! minimize the dipole moment
  DO
     oxygen_old = oxygen
     dipole_old = dipole
     ! propose a trial new configuration
     CALL shake_bonds(nshake)     
     CALL mc
     ! now the sample has again a correct proton disorder
     CALL get_dipole(dipole)
     WRITE(33,*)'dipole moment = ', dipole
     IF(dipole < dipole_old) THEN 
        WRITE(33,*)'new configuration accepted'
     ELSE
        oxygen = oxygen_old
        dipole = dipole_old 
     END IF
     IF(dipole < dipole_small) EXIT
  END DO

  ! construct the positions of hydrogens for given intramolecular geometry
  CALL put_hydrogens
  ! write the atomic coordinates
  CALL output
  STOP
  !**************************END******************************************
CONTAINS 
  !*************************************************************************
  SUBROUTINE lattice
    IMPLICIT NONE 

    ! lattice size
    INTEGER, PARAMETER :: imax=5,jmax=6,kmax=3,ibnum=4  ! input params (cells replicas)
    ! lattice parameters 
    REAL(8), PARAMETER :: a=4.5d0,c=7.32d0,u=(3.d0/8.d0)

    REAL(8), DIMENSION(3), PARAMETER :: &
         &x=(/1.d0,0.d0,0.d0/),y=(/0.d0,1.d0,0.d0/),z=(/0.d0,0.d0,1.d0/)
    REAL(8), DIMENSION(3) :: rc,r,aa,bb,cc,aa0
    INTEGER :: i,j,k,ib,iatom
    REAL(8), DIMENSION(ibnum,3) :: base

    nmols=ibnum*imax*jmax*kmax
    ALLOCATE (oxygen(nmols), oxygen_old(nmols))

    aa0 = (1.d0/2.d0) *  a * X - (1.d0/2.d0) * SQRT(3.d0) * a * Y
    ! primitive cell vectors of the hexagonal diamond structure
    aa = a * X 
    bb = (1.d0/2.d0) *  a * X + (1.d0/2.d0) * SQRT(3.d0) * a * Y
    cc = c * z
    ! basis vectors 1,...,ibnum
    base(1,:)  =  (1.d0/3.d0) * aa0 + (2.d0/3.d0) * bb 
    base(2,:)  =  (2.d0/3.d0) * aa0 + (1.d0/3.d0) * bb + (1.d0/2.d0) * cc  
    base(3,:)  =  (1.d0/3.d0) * aa0 + (2.d0/3.d0) * bb + u * cc  
    base(4,:)  =  (2.d0/3.d0) * aa0 + (1.d0/3.d0) * bb + ((1.d0/2.d0) + u) * cc 

    ! generate lattice
    iatom = 0
!!$    DO i = 0,imax-1
!!$       DO j = 0,jmax-1
    DO j = 0,jmax-1
       DO i = -INT(j/2),-INT(j/2) + imax-1
          DO k = 0,kmax-1
             rc = i*aa + j*bb + k*cc
             DO ib=1,ibnum
                iatom = iatom + 1
                oxygen(iatom) % r = rc + base(ib,:)
             END DO
          END DO
       END DO
    END DO

    box % hmat(:,1) = DBLE(imax*aa)
!!$    box % hmat(:,2) = DBLE(jmax*bb)
    box % hmat(:,2) = DBLE(jmax*a*SQRT(3.d0)/2.d0)*y
    box % hmat(:,3) = DBLE(kmax*cc)

    ! debug
    WRITE(*,*) 'oxygen lattice generated',nmols,iatom

  END SUBROUTINE lattice
  !*************************************************************************
  SUBROUTINE neighbour_list
    ! creates nearest neighbour_list
    IMPLICIT NONE 

    REAL(8), PARAMETER :: rcut=3.d0, rcut2=rcut**2 ! Declared a real valued constant
    INTEGER :: near(nmols), i, j, ilink
    REAL(8) :: rin(3), rout(3), s(3)

    CALL get_hinv ( box )
    ALLOCATE (link(2*nmols))

    ! calculate all distances between the atoms
    near(:) = 0 ! matrix counting the number of nearby neighors
    ilink = 0
    DO i=1,nmols-1
       DO j=i+1,nmols 
          rin(:) = oxygen(i) % r - oxygen(j) % r
          s = MATMUL ( box % h_inv, rin )
          s = s - NINT ( s )
          rout = MATMUL ( box % hmat, s )
          IF (SUM(rout*rout) <= rcut2) THEN 
             near(i) = near(i) + 1 ! update the number of neighbors in list
             near(j) = near(j) + 1 
             IF(near(i) <= 4 .AND. near(j) <= 4) THEN ! if not 4 neighbors, then update the oxygen neighbors
                oxygen(i) % part(near(i)) % neighbour = j
                oxygen(j) % part(near(j)) % neighbour = i
                ilink = ilink + 1
                !debug
                WRITE(*,*)i,j,near(i),near(j),ilink
                link(ilink) % part(1) % bond => oxygen(i) %part(near(i)) % bond
                link(ilink) % part(1) % nbonds => oxygen(i) % nbonds
                link(ilink) % part(2) % bond => oxygen(j) %part(near(j)) % bond
                link(ilink) % part(2) % nbonds => oxygen(j) % nbonds
             ELSE
                !debug
                STOP '# of neighbours > 4, wrong lattice coordination'
             END IF
          END IF
       END DO
    END DO

    ! debug
    WRITE(*,*) 'neighbour_list created'

  END SUBROUTINE neighbour_list
  !*************************************************************************
  SUBROUTINE init_hydrogens
    IMPLICIT NONE 
    REAL :: z
    INTEGER :: i

    ! put 1 hydrogen on each link, and 'give' it to one of the oxygens
    DO i = 1, 2*nmols
       CALL RANDOM_NUMBER (z)
       IF(z < 0.5) THEN
          link(i) % part(1) % bond  = .TRUE.
       ELSE
          link(i) % part(1) % bond  = .FALSE.
       END IF
       link(i) % part(2) % bond  = .NOT. link(i) % part(1) % bond
    END DO

    DO i = 1, nmols
       oxygen(i) % nbonds = COUNT(oxygen(i) % part(:) % bond)
       ! debug
       WRITE(*,*)i, oxygen(i) % nbonds
    END DO
    ! debug
    WRITE(*,*) 'hydrogens initialized randomly'

  END SUBROUTINE init_hydrogens
  !*************************************************************************
  SUBROUTINE shake_bonds(nshake)
    IMPLICIT NONE 

    INTEGER, INTENT(in) :: nshake
    INTEGER :: i, nbonds(2), icount
    REAL :: z
!!$    LOGICAL :: accept

    DO icount = 1, nshake
       DO
          CALL RANDOM_NUMBER (z)
          i = INT(2*nmols*z) + 1
          IF(1 <= i .AND. i <= 2*nmols) EXIT    
       END DO

       nbonds(1) = link(i) % part(1) % nbonds 
       nbonds(2) = link(i) % part(2) % nbonds 
!!$       diff_old = ABS(nbonds(1) - nbonds(2))
       IF( link(i) % part(1) % bond ) THEN
          nbonds(1) = nbonds(1) - 1
          nbonds(2) = nbonds(2) + 1
       ELSE
          nbonds(1) = nbonds(1) + 1
          nbonds(2) = nbonds(2) - 1
       END IF
!!$       diff_new = ABS(nbonds(1) - nbonds(2))
!!$       accept = .FALSE.
!!$       IF (diff_new <= diff_old) THEN
!!$          accept = .TRUE.
!!$          ! accept the move
       link(i) % part(1) % bond   = .NOT. link(i) % part(1) % bond 
       link(i) % part(2) % bond   = .NOT. link(i) % part(2) % bond 
       link(i) % part(1) % nbonds = nbonds(1)
       link(i) % part(2) % nbonds = nbonds(2)
!!$       END IF
!!$
!!$       icount = icount + 1
       ! stop when all oxygens have 2 bonds
!!$       IF(ALL(oxygen(:) % nbonds == 2)) EXIT
       ! debug
!!$       WRITE(*,*)icount, ALL(oxygen(:) % nbonds == 2), accept
    END DO
    !debug
!!$    WRITE(*,*) icount, ' MCS performed'
!!$    DO i = 1,nmols
!!$       WRITE(*,*)i, COUNT(oxygen(i) % part(:) % bond)
!!$    END DO

  END SUBROUTINE shake_bonds
  !*************************************************************************
  SUBROUTINE mc
    IMPLICIT NONE 

    INTEGER :: i, nbonds(2), diff_old, diff_new, icount
    REAL :: z
    LOGICAL :: accept

    icount = 0
    DO 
       DO
          CALL RANDOM_NUMBER (z)
          i = INT(2*nmols*z) + 1
          IF(1 <= i .AND. i <= 2*nmols) EXIT    
       END DO

       nbonds(1) = link(i) % part(1) % nbonds ! dot operator for fortran
       nbonds(2) = link(i) % part(2) % nbonds 
       diff_old = ABS(nbonds(1) - nbonds(2))
       IF( link(i) % part(1) % bond ) THEN
          nbonds(1) = nbonds(1) - 1
          nbonds(2) = nbonds(2) + 1
       ELSE
          nbonds(1) = nbonds(1) + 1
          nbonds(2) = nbonds(2) - 1
       END IF
       diff_new = ABS(nbonds(1) - nbonds(2))
       accept = .FALSE.
       IF (diff_new <= diff_old) THEN
          accept = .TRUE.
          ! accept the move
          link(i) % part(1) % bond   = .NOT. link(i) % part(1) % bond 
          link(i) % part(2) % bond   = .NOT. link(i) % part(2) % bond 
          link(i) % part(1) % nbonds = nbonds(1)
          link(i) % part(2) % nbonds = nbonds(2)
       END IF

       icount = icount + 1
       ! stop when all oxygens have 2 bonds
       IF(ALL(oxygen(:) % nbonds == 2)) EXIT
!!$       ! debug
!!$       WRITE(*,*)icount, ALL(oxygen(:) % nbonds == 2), accept
    END DO
    !debug
    WRITE(*,*) icount, ' MCS performed'
    DO i = 1,nmols
       WRITE(*,*)i, COUNT(oxygen(i) % part(:) % bond)
    END DO

  END SUBROUTINE mc
  !*************************************************************************
  SUBROUTINE get_dipole(dipole)
    IMPLICIT NONE 

    REAL(8), INTENT (out) :: dipole
    INTEGER :: i, j, k, in
    REAL(8), DIMENSION(3) :: rin, s, u, dip
    REAL(8), DIMENSION(2,3) :: ro
    REAL(8) :: cosal2, sinal2 

    dip = 0.d0
    DO i = 1,nmols
       ! now we construct the hydrogen positions
       k = 0
       DO j = 1,4
          IF(oxygen(i) % part(j) % bond) THEN
             in = oxygen(i) % part(j) % neighbour
             rin = oxygen(in) % r - oxygen(i) % r
             s = MATMUL ( box % h_inv, rin )
             s = s - NINT ( s )
             rin = MATMUL ( box % hmat, s )
             k = k + 1
             ro(k,:) = rin 
          END IF
       END DO
       u = ro(1,:) + ro(2,:)
       u = u/SQRT(SUM(u*u))
!!$       v = ro(1,:) - ro(2,:)
!!$       v = v/SQRT(SUM(v*v))
!!$       rh(1,:) = (cosal2 * u + sinal2 * v)
!!$       rh(2,:) = (cosal2 * u - sinal2 * v)
       dip = dip + u
    END DO
    dipole = SQRT(SUM(dip*dip))

  END SUBROUTINE get_dipole
  !*************************************************************************  
  SUBROUTINE put_hydrogens
    IMPLICIT NONE 
    ! TIP4P parameters
    REAL(8), PARAMETER :: r_oh = 0.9572d0 ! Angstroem
    REAL(8), PARAMETER :: ang_hoh=104.52d0 * (pi/180.d0)
    REAL(8), PARAMETER :: r_om = 0.15d0   ! Angstroem
    REAL(8), DIMENSION(3) :: r, rin, s, u, v
    REAL(8), DIMENSION(2,3) :: ro, rh
    REAL(8) :: cosal2, sinal2 
    INTEGER :: i, j, k, in, icount

    cosal2 = COS(ang_hoh/2.d0)
    sinal2 = SIN(ang_hoh/2.d0)
    ALLOCATE (atom(4*nmols))

    icount = 0
    DO i = 1,nmols
       ! now we construct the hydrogen positions
       k = 0
       DO j = 1,4
          IF(oxygen(i) % part(j) % bond) THEN
             in = oxygen(i) % part(j) % neighbour
             rin = oxygen(in) % r - oxygen(i) % r
             s = MATMUL ( box % h_inv, rin )
             s = s - NINT ( s )
             rin = MATMUL ( box % hmat, s )
             k = k + 1
             ro(k,:) = rin 
          END IF
       END DO
       u = ro(1,:) + ro(2,:)
       u = u/SQRT(SUM(u*u))
       v = ro(1,:) - ro(2,:)
       v = v/SQRT(SUM(v*v))
       rh(1,:) = oxygen(i) % r + r_oh * (cosal2 * u + sinal2 * v)
       rh(2,:) = oxygen(i) % r + r_oh * (cosal2 * u - sinal2 * v)
       icount = icount + 1
       atom(icount) % nummol = i
       atom(icount) % r = oxygen(i) % r + r_om * u
       atom(icount) % alabel = 'M '
       atom(icount) % alabel_gro = '   DW'
       icount = icount + 1
       atom(icount) % nummol = i
       atom(icount) % r = rh(1,:)
       atom(icount) % alabel = 'H '
       atom(icount) % alabel_gro = '  HW2'
       icount = icount + 1
       atom(icount) % nummol = i
       atom(icount) % r = rh(2,:)
       atom(icount) % alabel = 'H '
       atom(icount) % alabel_gro = '  HW3'
       icount = icount + 1
       atom(icount) % nummol = i
       atom(icount) % r = oxygen(i) % r
       atom(icount) % alabel = 'O '
       atom(icount) % alabel_gro = '  OW1'
    END DO
    ! debug
    WRITE(*,*) 'hydrogens placed'

  END SUBROUTINE put_hydrogens
  !************************************************************************* 
  SUBROUTINE output
    IMPLICIT NONE 
    INTEGER :: i
    REAL(8), PARAMETER :: zero = 0.d0
    CHARACTER(len=5), PARAMETER :: mol_label='WATER'
    CHARACTER(len=5) :: c_nmols
    CHARACTER(len=11) :: filename

    ! output of proton disordered structure
    WRITE(c_nmols,'(i5)') nmols
    filename = c_nmols//'_tip4p'
    OPEN(10,file=TRIM(filename)//'.dat')
    OPEN(12,file=TRIM(filename)//'.xyz')
    OPEN(14,file=TRIM(filename)//'.gro') ! for gromacs
    WRITE(10,*)4*nmols
    WRITE(12,*)3*nmols
    WRITE(12,*)'ice_Ih'
    WRITE(14,'(a6,3x,i5,3x,a7)')'MD of ', nmols, ' waters'
    WRITE(14,'(i5)')4*nmols

    DO i = 1,4*nmols
       ! we have to write the atoms molecule by molecule in the order D,H,H,O
       IF(MOD(i,4) /= 1) THEN
          WRITE(12,'(a2,4x,3(f12.6,4x))')atom(i) % alabel, atom(i) % r 
       END IF
       ! output for cp2k is in A, for gromacs in nm
       WRITE(10,'(3(f12.6,4x))')atom(i) % r 
       WRITE(14,'(i5,2a5,i5,3f8.3,3f8.4)')atom(i) % nummol, mol_label, &
            &atom(i) % alabel_gro, i, 0.1d0*atom(i) % r, zero, zero, zero 
    END DO

    DO i = 1, 3
       WRITE(10,'(3(f12.6,4x))') box % hmat(i,:)
    END DO
    WRITE(14,'(9f9.4)')&
         &0.1d0*box % hmat(1,1), 0.1d0*box % hmat(2,2), 0.1d0*box % hmat(3,3), &
         &0.1d0*box % hmat(2,1), 0.1d0*box % hmat(3,1), 0.1d0*box % hmat(1,2), &
         &0.1d0*box % hmat(3,2), 0.1d0*box % hmat(1,3), 0.1d0*box % hmat(2,3)

    CLOSE(10)
    CLOSE(12)
    CLOSE(14)

    WRITE(*,*)'# of atoms = ',3*nmols

  END SUBROUTINE output
  !************************************************************************* 
  SUBROUTINE get_hinv ( box )

    IMPLICIT NONE

    ! Arguments
    TYPE ( cell_type ), INTENT ( INOUT ) :: box

    ! Locals
    REAL ( 8 ), DIMENSION ( 3, 3 ) :: hmat, hmati
    REAL ( 8 ) :: odet

    !-----------------------------------------------------------------------

    hmat = box % hmat
    box % deth = &
         hmat(1,1) * ( hmat(2,2)*hmat(3,3)-hmat(2,3)*hmat(3,2) ) + &
         hmat(1,2) * ( hmat(2,3)*hmat(3,1)-hmat(2,1)*hmat(3,3) ) + &
         hmat(1,3) * ( hmat(2,1)*hmat(3,2)-hmat(2,2)*hmat(3,1) )
    IF ( box % deth < 1.0d-10 ) STOP 'get_hinv, box determinant too small'
    odet = 1.0d0 / box % deth
    hmati(1,1) = (hmat(2,2)*hmat(3,3)-hmat(2,3)*hmat(3,2))*odet
    hmati(2,2) = (hmat(1,1)*hmat(3,3)-hmat(1,3)*hmat(3,1))*odet
    hmati(3,3) = (hmat(1,1)*hmat(2,2)-hmat(1,2)*hmat(2,1))*odet
    hmati(1,2) = (hmat(1,3)*hmat(3,2)-hmat(1,2)*hmat(3,3))*odet
    hmati(2,1) = (hmat(3,1)*hmat(2,3)-hmat(2,1)*hmat(3,3))*odet
    hmati(1,3) = (hmat(1,2)*hmat(2,3)-hmat(1,3)*hmat(2,2))*odet
    hmati(3,1) = (hmat(2,1)*hmat(3,2)-hmat(3,1)*hmat(2,2))*odet
    hmati(2,3) = (hmat(1,3)*hmat(2,1)-hmat(2,3)*hmat(1,1))*odet
    hmati(3,2) = (hmat(3,1)*hmat(1,2)-hmat(3,2)*hmat(1,1))*odet

    box % h_inv = hmati

  END SUBROUTINE get_hinv
  !************************************************************************* 
END PROGRAM iceIh_disorder