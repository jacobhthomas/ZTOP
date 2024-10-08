
c======================================================================
      subroutine pftopdg(x,q,dxpdf)
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
      double precision x,q,dxpdf(-6:6)
      integer i
c
c     External
c
      double precision CtqPdf
      external CtqPdf
c     ----------
c     Begin Code
c     ----------
c t~ b~ c~ s~ u~ d~ 0 d u s c b t
      call structm(x,q,dxpdf(2),dxpdf(1),dxpdf(-2),dxpdf(-1),
     &     dxpdf(-3),dxpdf(-4),dxpdf(-5),dxpdf(-6),dxpdf(0))
c
      do i=3,6
         dxpdf(i) = dxpdf(-i)
      enddo
      dxpdf(2) = dxpdf(2) + dxpdf(-2)
      dxpdf(1) = dxpdf(1) + dxpdf(-1)
c
      return
      end

