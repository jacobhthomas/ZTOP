
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
      double precision CtqPdf, partonx12
      external CtqPdf, partonx12

      integer MXX, MXQ, MXF, MaxVal, MXPQX, MXPQXTOT
      PARAMETER (MXX = 201, MXQ = 40, MXF = 6, MaxVal=4)
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)
      parameter (MXPQXTOT = MXPQX*59)
      integer icallset
      double precision updmulti(MXPQXTOT)
      common / ctmulti / updmulti, icallset

c     ----------
c     Begin Code
c     ----------
c t~ b~ c~ s~ u~ d~ 0 d u s c b t
c      call structm(x,q,dxpdf(2),dxpdf(1),dxpdf(-2),dxpdf(-1),
c     &     dxpdf(-3),dxpdf(-4),dxpdf(-5),dxpdf(-6),dxpdf(0))

      do icallset = 0, 58
         dxpdf(-6 + 13*icallset) = 0d0
         dxpdf(-5 + 13*icallset) = x*partonx12(5,x,q)
         dxpdf(-4 + 13*icallset) = x*partonx12(4,x,q)
         dxpdf(-3 + 13*icallset) = x*partonx12(3,x,q)
         dxpdf(-2 + 13*icallset) = x*partonx12(-1,x,q)
         dxpdf(-1 + 13*icallset) = x*partonx12(-2,x,q)
         dxpdf( 0 + 13*icallset) = x*partonx12(0,x,q)
         dxpdf( 1 + 13*icallset) = x*partonx12(2,x,q)
         dxpdf( 2 + 13*icallset) = x*partonx12(1,x,q)
         dxpdf( 3 + 13*icallset) = dxpdf(-3 + 13*icallset)
         dxpdf( 4 + 13*icallset) = dxpdf(-4 + 13*icallset)
         dxpdf( 5 + 13*icallset) = dxpdf(-5 + 13*icallset)
         dxpdf( 6 + 13*icallset) = 0d0
      enddo
c
      return
      end

