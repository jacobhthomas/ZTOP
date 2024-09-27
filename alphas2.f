C======================================================================
      DOUBLE PRECISION FUNCTION ALPHAS2(Q)
C======================================================================
c     This is like the function alphas2 from CERNLIB:
c        a = alphas2(q)
c
      implicit none
      double precision q, alphaspdf
      external alphaspdf

      alphas2 = alphasPDF(q)

      RETURN
      END
