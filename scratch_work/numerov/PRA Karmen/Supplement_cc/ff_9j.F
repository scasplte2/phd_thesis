#include "fintrf.h"


      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C-----------------------------------------------------------------------
C
      mwpointer plhs(*), prhs(*)
      mwpointer mxGetPr, mxCreateDoubleMatrix
      mwpointer x_pr, y_pr
C-----------------------------------------------------------------------
C

      mwsize nlhs, nrhs
      integer mxIsNumeric
      mwsize mxGetM, mxGetN
      mwsize m, n, size, sizey, my, ny
      real*8  y
      real*8  x(100)

C     Check for proper number of arguments. 
      if(nrhs .ne. 1) then
         call mexErrMsgTxt('One input required.')
      elseif(nlhs .gt. 1) then
         call mexErrMsgTxt('Too many output arguments.')
      endif

C     Get the size of the input array.
      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      size = m*n
      sizey = 1

C     Check to insure the input is a number.
      if(mxIsNumeric(prhs(1)) .eq. 0) then
         call mexErrMsgTxt('Input must be a number.')
      endif

C     Create matrix for the return argument.
      my = 1
      ny = 1
      plhs(1) = mxCreateDoubleMatrix(my,ny,0)
      x_pr = mxGetPr(prhs(1))
      y_pr = mxGetPr(plhs(1))
      call mxCopyPtrToReal8(x_pr,x, size)
c      call mxCopyPtrToInteger4(x_pr,x,size)

C     Call the computational subroutine.
      call ff_ninej(x, y)

C     Load the data into y_pr, which is the output to MATLAB
      call mxCopyReal8ToPtr(y,y_pr,sizey)

      return
      end



c
c**********************************************************************
c
      subroutine ff_threej(qnarray, tj)
      implicit none
      real*8 qnarray(2, 3), tj
c
      real*8 j1, j2, j3, m1, m2, m3, sign, binom
      integer i1, i2, i3, k1, k2, k3, l1, l2, l3, n1, n2, n3,
     +  imin1, imax1, iteller
c
      character*11 head/'FF_THREEJ: '/
c
      j1 = qnarray(1,1)
      j2 = qnarray(1,2)
      j3 = qnarray(1,3)
      m1 = qnarray(2,1)
      m2 = qnarray(2,2)
      m3 = qnarray(2,3)
C hack voor afronden
      j1=0.5d0*int(2.00001d0*j1)
      m1=0.5d0*int(2.00001d0*m1)
      j2=0.5d0*int(2.00001d0*j2)
      m2=0.5d0*int(2.00001d0*m2)
      j3=0.5d0*int(2.00001d0*j3)
      m3=0.5d0*int(2.00001d0*m3)

      tj = 0d0
      if (2*j1.ne.int(2*j1)) then
         call mexErrMsgTxt('j1 must be (half) integer')
         return
      endif
      if (2*j2.ne.int(2*j2)) then
         call mexErrMsgTxt('j2 must be (half) integer')
         return
      endif
      if (abs(m1).gt.j1) then
c         call mexErrMsgTxt('|m1| > j1')
         return
      endif
      if (abs(m2).gt.j2) then
c         call mexErrMsgTxt('|m2| > j2')
         return
      endif
      if (abs(m3).gt.j3) then
c         call mexErrMsgTxt('|m3| > j3')
         return
      endif
      if (j1+m1.ne.int(j1+m1)) then
         call mexErrMsgTxt('j1+m1 is not integer')
         return
      endif
      if (j2+m2.ne.int(j2+m2)) then
         call mexErrMsgTxt('j2+m2 is not integer')
         return
      endif
      if (j3+m3.ne.int(j3+m3)) then
         call mexErrMsgTxt('j3+m3 is not integer')
         return
      endif

      if (m1+m2+m3.ne.0d0) return
      i1 = -j1+j2+j3
      if (i1.lt.0) return
      i2 = j1-j2+j3
      if (i2.lt.0) return
      i3 = j1+j2-j3
      if (i3.lt.0) return
      k1 = j1+m1
      if (k1.lt.0) return
      k2 = j2+m2
      if (k2.lt.0) return
      k3 = j3+m3
      if (k3.lt.0) return
      l1 = j1-m1
      if (l1.lt.0) return
      l2 = j2-m2
      if (l2.lt.0) return
      l3 = j3-m3
      if (l3.lt.0) return
      n1 = -j1-m2+j3
      n2 =  m1-j2+j3
      n3 =  j1-j2-m3

      imin1 = max(-n1,-n2,0)
      imax1 = min(l1,k2,i3)
      if (imin1.gt.imax1) return
      sign = -1d0
      do iteller = imin1,imax1
        sign = -sign
        tj   = tj + sign * binom(i1,n1+iteller) * 
     .                     binom(i2,n2+iteller) *
     .                     binom(i3,iteller)
      end do
      
      tj = tj * sqrt(binom(int(j2+j2),i3)*binom(int(j1+j1),i2)/
     .              (binom(int(j1+j2+j3+1),i3)*(j3+j3+1)*
     .               binom(int(j1+j1),l1)*binom(int(j2+j2),l2)*
     .               binom(int(j3+j3),l3)))
      if (2*int((n3+imin1)/2).ne.n3+imin1) then
        tj = -tj
      end if
      end
c
c**********************************************************************
c
      subroutine ff_sixj(qnarray, sj)
      implicit none
      real*8 qnarray(2,3), sj
c
c  QNARRAY=[J1 J2 J3; J4 J5 J6]
c  Wigner 6j-symbol, integer and half integer quantum numbers.
c
c  GCG, 18-Mar-2003

      real*8  j1, j2, j3, j4, j5, j6, ak1, al1, fact, phase, am1, an3,
     +        binom
      integer k1, k2, k3, l1, l2, l3, m1, m2, m3, n1, n2, n3, i1, i2,
     +        imin1, imax1, iteller
      character*(9) head/'SF_SIXJ: '/

      j1 = qnarray(1,1)
      j2 = qnarray(1,2)
      j3 = qnarray(1,3)
      j4 = qnarray(2,1)
      j5 = qnarray(2,2)
      j6 = qnarray(2,3)
C hack voor afronden
      j1=0.5d0*int(2.00001d0*j1)
      j2=0.5d0*int(2.00001d0*j2)
      j3=0.5d0*int(2.00001d0*j3)
      j4=0.5d0*int(2.00001d0*j4)
      j5=0.5d0*int(2.00001d0*j5)
      j6=0.5d0*int(2.00001d0*j6)

      if (2*j1.ne.int(2*j1)) then
        call mexErrMsgTxt('j1 must be (half) integer')
        return
      endif
      if (2*j2.ne.int(2*j2)) then
        call mexErrMsgTxt('j2 must be (half) integer')
        return
      endif
      if (2*j6.ne.int(2*j6)) then
        call mexErrMsgTxt('j6 must be (half) integer')
        return
      endif
      sj  =  0d0
      ak1 = -j1+j2+j3
      k1  =  ak1
      if (ak1.ne.k1) then
        call mexErrMsgTxt('j1+j2+j3 must be integer')
        return
      endif
      k2  =  j1-j2+j3
      k3  =  j1+j2-j3
      al1 = -j1+j5+j6
      l1  =  al1
      if (l1.ne.al1) then
        call mexErrMsgTxt('j1+j5+j6 must be integer')
        return
      endif
      l2  =  j1-j5+j6
      l3  =  j1+j5-j6
      am1 = -j4+j2+j6
      m1  =  am1
      if (m1.ne.am1) then
        call mexErrMsgTxt('j2+j4+j6 must be integer')
        return
      endif
      m2  =  j4-j2+j6
      m3  =  j4+j2-j6
      n1  = -j4+j5+j3
      n2  =  j4-j5+j3
      an3 =  j4+j5-j3
      n3  =  an3
      if (n3.ne.an3) then
        call mexErrMsgTxt('j3+j4+j5 must be integer')
        return
      endif

c  Triangular conditions OK?

      if (min(k1,k2,k3,l1,l2,l3,m1,m2,m3,n1,n2,n3).lt.0) return
c     j1+j2+j3 is integer
c     j4+j5+j3 is integer
c     j1+j2+j4+j5 is integer
      i1   = j1+j2+j4+j5+1
      imin1 = max(0,k3-m1,k3-l2)
      imax1 = min(k3,l3,m3,n3,i1)
      if (imin1.gt.imax1) return
      i2   = j2+j4+j6+1
      fact = dsqrt( binom(i2,m1) * binom(int(j4+j4),n2) / 
     .      ((k1+1) * (l1+1) * binom(int(j4+j4),m3) * 
     .       binom(int(j3+j4+j5)+1,n1) * binom(int(j1+j2+j3)+1,k3) *
     .       binom(int(j3+j3)+1,k2) * binom(int(j1+j5+j6)+1,l3) *
     .       binom(int(j6+j6)+1,l2)))
      phase = 1d0
      if (2*((imin1+i1)/2).ne.(imin1+i1)) phase=-phase
      do iteller=imin1,imax1
        phase = -phase
        sj = sj + phase * binom(n3,iteller)*binom(m1,k3-iteller)*
     .       binom(n2,m3-iteller) * binom(i1-iteller,i2)
      enddo
      sj = sj * fact
      end
c
c**********************************************************************
c
      real*8 function binom(n, k0)
      implicit none
      integer n, k0
      integer i, k
      real*8 dnum, den, round
c
      k = k0
      if (k.lt.0.or.k.gt.n) then
         call mexErrMsgTxt('BINOM: k<0 or k>n')
         return
      endif
      if (2*k.gt.n) k = n-k
      if (k.le.1) then
        if (k.eq.0) then
          binom = 1
        else
          binom = n
        endif
      else
        binom  = 1d0
        dnum = n-k+1d0
        den  = 1d0
        do i=1,k
          binom  = binom*(dnum/den)
          dnum = dnum+1d0
          den  = den+1d0
        end do
      endif
      binom = round(binom)
      end
c
c**********************************************************************
c
      subroutine ff_ninej(qnarray, a9j)
      implicit none
      real*8 qnarray(3,3), a9j
c
c     SF_NINEJ(QNARRAY)  Wigner's 9-j symbol.
c     QNARRAY = [ J11 J12 J13
c                 J21 J22 J23
c                 J31 J32 J33]
c     returns the 9-j symbol for both integer and half-integer
c     quantum number.
c      
c     From: Quantum theory of Angular Momentum,
c           D. A. Varshalovich, A. N. Moskalev, V. K. Kheersonskii,
c           10.2.4, eq (20), p 340
c      
c     GCG, 18-Mar-2003
c
      real*8 j11, j12, j13, j21, j22, j23, j31, j32, j33, dim, zmin,
     +       zmax, az, sxj1, sxj2, sxj3, q1(2,3), q2(2,3), q3(2,3)
      integer i, j, n
c
      j11 = qnarray(1,1)
      j12 = qnarray(1,2)
      j13 = qnarray(1,3)
      j21 = qnarray(2,1)
      j22 = qnarray(2,2)
      j23 = qnarray(2,3)
      j31 = qnarray(3,1)
      j32 = qnarray(3,2)
      j33 = qnarray(3,3)

c	hack voor afronden
      j11=0.5d0*int(2.00001d0*j11)
      j12=0.5d0*int(2.00001d0*j12)
      j13=0.5d0*int(2.00001d0*j13)
      j21=0.5d0*int(2.00001d0*j21)
      j22=0.5d0*int(2.00001d0*j22)
      j23=0.5d0*int(2.00001d0*j23)
      j31=0.5d0*int(2.00001d0*j31)
      j32=0.5d0*int(2.00001d0*j32)
      j33=0.5d0*int(2.00001d0*j33)

      do j=1,3
        q1(1,j) = qnarray(1,j)
        q2(1,j) = qnarray(2,j)
        q3(1,j) = qnarray(3,j)
      enddo
      q1(2,1) = j23
      q1(2,2) = j33
      q2(2,1) = j12
      q2(2,3) = j32
      q3(2,2) = j11
      q3(2,3) = j21
       
      a9j  = 0d0
      zmin = max(abs(j11-j33),abs(j32-j21),abs(j12-j23))
      zmax = min(j11+j33, j32+j21, j12+j23)
      n    = nint(zmax-zmin)+1
      
      dim   = 2d0*zmin + 1d0
      az    = zmin
      do i = 1, n
        q1(2,3) = az
        q2(2,2) = az
        q3(2,1) = az
        call ff_sixj(q1, sxj1)
        call ff_sixj(q2, sxj2)
        call ff_sixj(q3, sxj3)
        a9j     = a9j + dim * sxj1 * sxj2 * sxj3
        dim     = dim + 2d0
        az      = az + 1d0
      enddo
      a9j    = a9j*(-1)**(nint(2d0*zmin))
      end
c
c**********************************************************************
c
      real*8 function round(x1)
c
c     Copyright INRIA
      real*8 x1,x,y,z,e,h
      data h/1.0d+9/
c
      x=x1
      if (x.eq.0d0) then 
         round=x
         return
      endif
      if ((2.0d0*x).eq.dble(int(2.d0*x))) then
c     changing the signs gives round(0.5)=0 round(-0.5)=0
         if (x.gt.0.d0) x=x+1.d-10
         if (x.lt.0.d0) x=x-1.d-10
      endif
      z = abs(x)
c     -----testing Nans 
      if (isnan(x)) then 
         round=x
         return
      endif
      y = z + 1.0d+0
      if (y .eq. z) then 
         round=x
         return
      endif
      y = 0.0d+0
      e = h
 10   if (e .ge. z) go to 20
      e = 2.0d+0*e
      go to 10
 20   if (e .le. h) go to 30
      if (e .le. z) y = y + e
      if (e .le. z) z = z - e
      e = e/2.0d+0
      go to 20
 30   z = int(z + 0.50d+0)
      y = y + z
      if (x .lt. 0.0d+0) y = -y
      round = y
      return
      end
