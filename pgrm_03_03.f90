      program pgrm_03_03
!
!     This program computes components of the Hartree-Fock energy and the number
!     of electrons using contraction of various matrices loaded from
!     user-provided files (that presumably come from Gaussian calculations).
!
!     At run time, the program expects 5 command line arguments:
!       (1) the number of electrons;
!       (2) the number of atomic orbital basis functions;
!       (3) an input file containing core-Hamiltonian matrix elements (in
!           symmetric upper/column storage form);
!       (4) an input file containing the overlap matrix elements (in
!           symmetric upper/column storage form); and
!       (5) an input file containing Fock matrix elements (in symmetric
!           upper/column storage form).
!
!     At run time, the program outputs 6 items:
!       (1) the MO energies;
!       (2) the MO coefficients;
!       (3) the density matrix;
!       (4) the one-electron contribution to the total electronic energy;
!       (5) the two-electron contribution to the total electronic energy; and
!       (6) the tr(PS), which should equal the number of electrons.
!
!
!
      implicit none
      integer,parameter::unitIn=10
      integer::i,iError,nElectrons,nOcc,nBasis,lenSym,j,k
      real::oneElectronEnergy,twoElectronEnergy,tracePS
      real,dimension(:),allocatable::symFock,symCoreHamiltonian,symOverlap,  &
        moEnergies,tempSymMatrix
      real,dimension(:,:),allocatable::sqFock,fockTilde,sqCoreHamiltonian,  &
        InvSqrtOverlap,moCoefficients,densityMatrix,tempSqMatrix,  &
        SSPEV_Scratch, temp_core
      character(len=256)::cmdlineArg
!
!
!     Begin by reading the number of basis functions, allocating array memory,
!     and loading the symmetric Fock and overlap matrices from input files
!     provided on the command line.
!
      call Get_Command_Argument(1,cmdlineArg)
      read(cmdlineArg,'(I)') nElectrons
      nOcc = nElectrons/2
!
      call Get_Command_Argument(2,cmdlineArg)
      read(cmdlineArg,'(I)') nBasis
      lenSym = (nBasis*(nBasis+1))/2
      allocate(symFock(lenSym),symCoreHamiltonian(lenSym),  &
        symOverlap(lenSym),tempSymMatrix(lenSym))
      allocate(moEnergies(nBasis))
      allocate(sqFock(nBasis,nBasis),fockTilde(nBasis,nBasis),  &
        sqCoreHamiltonian(nBasis,nBasis),invSqrtOverlap(nBasis,nBasis),  &
        moCoefficients(nBasis,nBasis),densityMatrix(nBasis,nBasis),  &
        tempSqMatrix(nBasis,nBasis),SSPEV_Scratch(nBasis,3),temp_core(nBasis,nBasis))

      call Get_Command_Argument(3,cmdlineArg)
      open(Unit=unitIn,File=TRIM(cmdlineArg),Status='OLD',IOStat=iError)
      if(iError.ne.0) then
        write(*,*)' Error opening input file.'
        STOP
      endIf
      do i = 1,lenSym
        read(unitIn,*) symCoreHamiltonian(i)
      endDo
      close(Unit=unitIn)

      call Get_Command_Argument(4,cmdlineArg)
      open(Unit=unitIn,File=TRIM(cmdlineArg),Status='OLD',IOStat=iError)
      if(iError.ne.0) then
        write(*,*)' Error opening input file.'
        STOP
      endIf
      do i = 1,lenSym
        read(unitIn,*) symOverlap(i)
      endDo
      close(Unit=unitIn)

      call Get_Command_Argument(5,cmdlineArg)
      open(Unit=unitIn,File=TRIM(cmdlineArg),Status='OLD',IOStat=iError)
      if(iError.ne.0) then
        write(*,*)' Error opening input file.'
        STOP
      endIf
      do i = 1,lenSym
        read(unitIn,*) symFock(i)
      endDo
      close(Unit=unitIn)

!
!


!
!     Form the square-root of the overlap matrix.
!
      call InvSQRT_SymMatrix(nBasis,symOverlap,invSqrtOverlap)
!      write(*,*)' Inverse SQRT Matrix:'
!      call Print_Matrix_Full_Real(invSqrtOverlap,nBasis,nBasis)
!     Form fTilde and solve for the MO energies and coefficients.
!
      call SymmetricPacked2Matrix_UpperPacked(nBasis,symFock,sqFock)
      tempSqMatrix = MatMul(invSqrtOverlap,sqFock)
      fockTilde = MatMul(tempSqMatrix,invSqrtOverlap)
      call Sq2SymMatrix(nBasis,fockTilde,tempSymMatrix)
      call SSPEV('V','U',nBasis,tempSymMatrix,moEnergies,  &
        moCoefficients,nBasis,SSPEV_Scratch,iError)
      If(iError.ne.0) then
        write(*,*)' Failure in SSPEV.'
        STOP
      endIf
      write(*,*)' MO Energies:'
      call Print_Matrix_Full_Real(Reshape(moEnergies,[nBasis,1]),nBasis,1)
      tempSqMatrix = moCoefficients
      moCoefficients = MatMul(invSqrtOverlap,tempSqMatrix)
      write(*,*)' MO Coefficients:'
      call Print_Matrix_Full_Real(moCoefficients,nBasis,nBasis)
!
!
      densityMatrix = 0
      do i = 1, nBasis
        do j = 1, nBasis
           do k = 1, nOcc
            densityMatrix(i,j) = densityMatrix(i,j)+2*moCoefficients(i,k)*moCoefficients(j,k)
           endDo
        endDo
      endDo
      write(*,*) 'Density Matrix:'
      call Print_Matrix_Full_Real(densityMatrix,nBasis,nBasis)
      
      call SymmetricPacked2Matrix_UpperPacked(nBasis,symCoreHamiltonian,sqCoreHamiltonian)
!      call Print_Matrix_Full_Real(sqCoreHamiltonian,nBasis,nBasis)
      temp_core = MatMul(densityMatrix,sqCoreHamiltonian)
      oneElectronEnergy = 0
      do i = 1, nBasis
        oneElectronEnergy = oneElectronEnergy + temp_core(i,i)
      enddo

      write(*,*), 'One Electron Energy contribution = ',oneElectronEnergy     
!
      
      temp_core = 0.5*(MatMul(densityMatrix, (sqFock - sqCoreHamiltonian)))
      twoElectronEnergy = 0
      do i = 1, nBasis
        twoElectronEnergy = twoElectronEnergy + temp_core(i,i)
      endDo
      write(*,*),'Two-electron Energy contribution = ', twoElectronEnergy

      call SymmetricPacked2Matrix_UpperPacked(nBasis,symOverlap,temp_core)
      temp_core = MatMul(densityMatrix,temp_core)
      tracePS = 0
      do i = 1, nBasis
        tracePS = tracePS + temp_core(i,i)
      endDo
      write(*,*),'tr(PS) = <PS>', tracePS
!

      End program




      Subroutine InvSQRT_SymMatrix(nDim,inputSymMatrix,invSqrtInputMatrix)
      
!     Totally Done By Ajay Khanna, Yeah !!!
      implicit none
      integer:: nDim,i
      real:: IError
      real,dimension(nDim)::Evals
      real,dimension((nDim*(nDim+1)/2))::inputSymMatrix, temp_mat
      real,dimension(nDim,nDim):: invSqrtInputMatrix, Evecs,last_mat,mat_evals
      real,dimension(nDim*3)::temp_Vector

      temp_mat = inputSymMatrix
      mat_evals = 0

      call SSPEV('V','U',nDim, temp_mat,Evals,Evecs,ndim, &
         temp_Vector,IError)
!      call Print_Matrix_Full_Real(mat_evals,nDim,nDim)
      do i = 1, nDim
        mat_evals(i,i) = 1/(sqrt(Evals(i)))
      enddo
!      call Print_Matrix_Full_Real(mat_evals,nDim,nDim)

      invSqrtInputMatrix = MatMul(MatMul(Evecs,mat_evals),transpose(Evecs))
      return
      End Subroutine InvSQRT_SymMatrix

      Subroutine SymmetricPacked2Matrix_UpperPacked(N,ArrayIn,AMatOut)
!
!     This subroutine accepts an array, ArrayIn, that is (N*(N+1))/2
!     long.
!     It then converts that form to the N-by-N matrix AMatOut taking
!     ArrayIn to be in upper-packed storage form. Note: The storage mode
!     also assumes the upper-packed storage is packed by columns.
!
      Implicit None
      Integer,Intent(In)::N
      Real,Dimension((N*(N+1))/2),Intent(In)::ArrayIn
      Real,Dimension(N,N),Intent(Out)::AMatOut
!
      Integer::i,j,k
!
!     Loop through the elements of AMatOut and fill them appropriately
!     from
!     Array_Input.
!
!
! *************************************************************************
      k = 1
      do i = 1,N
       do j = 1,i
        AMatOut(j,i) = ArrayIn(k)
        AMatout(i,j) = AMatOut(j,i)
        k = k + 1
       enddo
      enddo
! *************************************************************************
!
!
      Return
      End Subroutine SymmetricPacked2Matrix_UpperPacked

      Subroutine Sq2SymMatrix(nDim,sqMatrix,symMatrix)
      implicit none
      integer::nDim, i, j, k = 1
      real,dimension(nDim,nDim)::sqMatrix
      real,dimension((nDim*(nDim+1))/2)::symMatrix

      Do i = 1, nDim
        Do j = 1, i
                symMatrix(k) = sqMatrix(i,j)
                k = k + 1
        EndDo
      EndDo
!      print *, 'lookhere',symMatrix
      return
      End Subroutine Sq2SymMatrix


      Subroutine Print_Matrix_Full_Real(AMat,M,N)
!
!     This subroutine prints a real matrix that is fully dimension -
!     i.e.,
!     not stored in packed form. AMat is the matrix, which is
!     dimensioned
!     (M,N).
!
!     The output of this routine is sent to unit number 6 (set by the
!     local
!     parameter integer IOut).
!
!
!     Variable Declarations
!
      implicit none
      integer,intent(in)::M,N
      real,dimension(M,N),intent(in)::AMat
!
!     Local variables
      integer,parameter::IOut=6,NColumns=5
      integer::i,j,IFirst,ILast
!
 1000 Format(1x,A)
 2000 Format(5x,5(7x,I7))
 2010 Format(1x,I7,5F14.6)
!
      Do IFirst = 1,N,NColumns
        ILast = Min(IFirst+NColumns-1,N)
        write(IOut,2000) (i,i=IFirst,ILast)
        Do i = 1,M
          write(IOut,2010) i,(AMat(i,j),j=IFirst,ILast)
        endDo
      endDo
!
      Return
      End Subroutine Print_Matrix_Full_Real

