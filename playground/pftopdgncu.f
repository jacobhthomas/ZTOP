c======================================================================
      subroutine pftopdgn(x,q,dxpdf)
c======================================================================
C   Before using the PDF, it is necessary to do the initialization by
C       Call SetCtq(Idist,Iset)
C     Where `Idist' is 4,5, or 6 and `Iset' is from the tables.
C
C   This routine requires STRUCTM.  It forwards the call to STRUCTM, and
C    therefore assumes (valid through CTEQ61) s=s_,c=c_,b=b_.
C   This is intended soley as a light-weight interface for use with an
C    optimized STRUCTM.
C
C   Written by Z. Sullivan, 5/20/04
C
c     *****************************************************************
c     * This is an interface for the CTEQ4/5/6 tables functions:      *
c     *  CtqPdf (Iparton, X, Q)                                       *
c     * This is called like pftopdg:                                  *
c     *   call pftopdg(x,q,dxpdf(-6:6))                               *
c     *****************************************************************
      implicit none
c
c     Local
c
      integer sets, xpdfmax
      parameter (sets = 59, xpdfmax = 13*(sets-1)+6)
      double precision x, q, dxpdf(-6:xpdfmax)

c      double precision x,q,dxpdf(-6:6)
      integer i
c
c     External
c
      double precision CtqPdf
      external CtqPdf

      integer MXX, MXQ, MXF, MaxVal, MXPQX, MXPQXTOT
      PARAMETER (MXX = 201, MXQ = 40, MXF = 6, MaxVal=4)
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)
      parameter (MXPQXTOT = MXPQX*59)
      integer icallset
      double precision updmulti(MXPQXTOT)
      common / ctmulti / updmulti, icallset

      integer iSetch, iParton, nx, nt, npts, NfMx, MxVal
      double precision XX, QQ, qB, ret
      cudaTextureObject_t xTex, tTex, pTex
      external partonx12_wrapper_

      iSetch = 1  ! Example value, set appropriately
      nx = MXX
      nt = MXQ
      npts = MXPQX
      NfMx = MXF
      MxVal = MaxVal
      qB = 1.0  ! Example value, set appropriately

c     ----------
c     Begin Code
c     ----------
c t~ b~ c~ s~ u~ d~ 0 d u s c b t
c      call structm(x,q,dxpdf(2),dxpdf(1),dxpdf(-2),dxpdf(-1),
c     &     dxpdf(-3),dxpdf(-4),dxpdf(-5),dxpdf(-6),dxpdf(0))

      do icallset = 0, 58
         dxpdf(-6 + 13*icallset) = 0d0
         iParton = 5
         call partonx12_wrapper_(iSetch, iParton, nx, nt, npts, NfMx, MxVal, x, q, qB, xTex, tTex, pTex, ret)
         dxpdf(-5 + 13*icallset) = x * ret
         iParton = 4
         call partonx12_wrapper_(iSetch, iParton, nx, nt, npts, NfMx, MxVal, x, q, qB, xTex, tTex, pTex, ret)
         dxpdf(-4 + 13*icallset) = x * ret
         iParton = 3
         call partonx12_wrapper_(iSetch, iParton, nx, nt, npts, NfMx, MxVal, x, q, qB, xTex, tTex, pTex, ret)
         dxpdf(-3 + 13*icallset) = x * ret
         iParton = -1
         call partonx12_wrapper_(iSetch, iParton, nx, nt, npts, NfMx, MxVal, x, q, qB, xTex, tTex, pTex, ret)
         dxpdf(-2 + 13*icallset) = x * ret
         iParton = -2
         call partonx12_wrapper_(iSetch, iParton, nx, nt, npts, NfMx, MxVal, x, q, qB, xTex, tTex, pTex, ret)
         dxpdf(-1 + 13*icallset) = x * ret
         iParton = 0
         call partonx12_wrapper_(iSetch, iParton, nx, nt, npts, NfMx, MxVal, x, q, qB, xTex, tTex, pTex, ret)
         dxpdf( 0 + 13*icallset) = x * ret
         iParton = 2
         call partonx12_wrapper_(iSetch, iParton, nx, nt, npts, NfMx, MxVal, x, q, qB, xTex, tTex, pTex, ret)
         dxpdf( 1 + 13*icallset) = x * ret
         iParton = 1
         call partonx12_wrapper_(iSetch, iParton, nx, nt, npts, NfMx, MxVal, x, q, qB, xTex, tTex, pTex, ret)
         dxpdf( 2 + 13*icallset) = x * ret
         dxpdf( 3 + 13*icallset) = dxpdf(-3 + 13*icallset)
         dxpdf( 4 + 13*icallset) = dxpdf(-4 + 13*icallset)
         dxpdf( 5 + 13*icallset) = dxpdf(-5 + 13*icallset)
         dxpdf( 6 + 13*icallset) = 0d0
      enddo
c
      return
      end