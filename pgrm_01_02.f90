      Program prgm_01_02
!
!     This program reads two 3x3 matrix from  user-provided input files.
!     After the
!     file is opened and read, it is closed and then printed.
!
!     Ajay Khanna, 2019.
!
      implicit none
      integer,parameter::inFileUnitA=10, inFileUnitB=11
      integer::errorFlag,i
      real,dimension(3,3)::matrixInA, matrixInB
      character(len=128)::fileNameA, fileNameB
!
!
!     Start by asking the user for the name of the data file.
!
      write(*,*)' What is the name of the input data files?'
      read(*,*) fileNameA, fileNameB
!
!     Open the data files and read matrixInA & matrixInB from the files.
!
      open(unit=inFileUnitA,file=TRIM(fileNameA),status='old',  &
        iostat=errorFlag)
      open(unit=inFileUnitB,file=TRIM(fileNameB),status='old',  &
        iostat=errorFlag)

      if(errorFlag.ne.0) then
        write(*,*)' There was a problem opening the input file(s).'
        goto 999
      endIf
      do i = 1,3
        read(inFileUnitA,*) matrixInA(1,i),matrixInA(2,i),matrixInA(3,i)
        read(inFileUnitB,*) matrixInB(1,i),matrixInB(2,i),matrixInB(3,i)
      endDo
      close(inFileUnitA)
      close(inFileUnitB)
!
!     Call the subroutine PrintMatrix to print matrixInA & matrixInB.
!
      call PrintMatrix3x3(matrixInA, matrixInB)
!
  999 continue
      End Program prgm_01_02


      Subroutine PrintMatrix3x3(matrixA, matrixB)
!
!     This subroutine prints two 3x3 real matrix. The output is written to
!     StdOut.
!
      implicit none
      real,dimension(3,3),intent(in)::matrixA, matrixB
      integer::i
!
!     Format statements.
!
 1000 format(3(2x,f5.1))
!
!     Do the printing job.
!
      print*,'Programmed by Ajay Khanna'
      print*
      write(*,*)' Printing Matrices'
!
      print*, 'Matrix A'
      do i = 1,3
        write(*,1000),matrixA(i,1), matrixA(i,2), matrixA(i,3)
      enddo
!      
      Print *

      print *, 'Matrix B'
      do i = 1,3
        write(*,1000),matrixB(i,1), matrixB(i,2), matrixB(i,3)
      enddo

      return
      End Subroutine PrintMatrix3x3
