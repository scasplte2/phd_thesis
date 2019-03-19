#include "fintrf.h"

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C-----------------------------------------------------------------------
C
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs, INFO
      mwPointer mxcreatedoublematrix
      mwPointer mxgetpr
      mwPointer A, E, U
      mwSignedIndex mxgetm, mxgetn
      mwSignedIndex n, LDA, numel, LWORK, k
      double precision ar 
      character JOBZ, UPLO
      real*8, allocatable :: WORK(:)

      JOBZ = 'V'
      UPLO = 'L'
      LWORK = -1

      A = mxgetpr(prhs(1))
      LDA = mxgetm(prhs(1))
      N = mxgetn(prhs(1))
      k = 1

      allocate(WORK(1))
      WORK(1)=0

      plhs(1) = mxcreatedoublematrix(LDA, k, 0.0)
      E = mxgetpr(plhs(1))
      numel = N*N
      plhs(2) = mxcreatedoublematrix(N,N,0.0)
      U = mxgetpr(plhs(2))
      call mxcopyptrtoreal8(A,%val(U),numel)
C      call mxcopyptrtoreal8(A, ar, numel)

      call dsyev(JOBZ, UPLO, N, %val(U), LDA, %val(E), WORK, LWORK,
     + INFO)

      LWORK=WORK( 1 )
      deallocate(WORK)
      allocate(WORK(LWORK))

      call dsyev(JOBZ, UPLO, N, %val(U), LDA, %val(E), WORK, LWORK,
     + INFO)

      deallocate(WORK)

      return
      end
