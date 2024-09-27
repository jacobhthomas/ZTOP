C======================================================================
      DOUBLE PRECISION FUNCTION ALPHAS2(Q)
C======================================================================
c     This is like the function alphas2 from CERNLIB:
c        a = alphas2(q)
c
      implicit none
      double precision q, CT18Alphas
      external CT18Alphas

      alphas2 = CT18Alphas(q)

      RETURN
      END


      double precision Function CT18Alphas (QQ)

      Implicit Double Precision (A-H,O-Z)

      PARAMETER (MXX = 201, MXQ = 40, MXF = 6, MaxVal=4)
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)
      double precision Alsout
      
      Common
     > / CtqPar1 / qBase,XV(0:MXX),TV(0:MXQ),UPD(MXPQX), AlsCTEQ(0:mxq)
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / XQrange / Qini, Qmax, Xmin
     > / Masstbl / Amass(6)
     > / QCDtbl /  AlfaQ, Qalfa, Ipk, Iorder, Nfl !for external use
     > /Setchange/ Isetch, ipdsset, ipdsformat

      Data Q, JQ /-1D0, 0/
      save

      if (ipdsset.ne.1) 
     >  STOP 'CT18Alphas: the PDF table was not initialized'

      
      if (ipdsformat.lt.11) then
        print *
        print *, 
     >    'STOP in CT18alphas: the PDF table file has an older format'
        print *,
     >    'and does not include the table of QCD coupling values.'
        print *, 
     >    'You can still compute the PDFs, but do not call'
        print *,
     >    'the CT18alphas function for the interpolation of alpha_s.'
        stop
      endif

      Q = QQ
      tt = dlog(dlog(Q/qBase))

c         --------------   Find lower end of interval containing Q, i.e.,
c                          get jq such that qv(jq) .le. q .le. qv(jq+1)...
      JLq = -1
      JU = NT+1
 13   If (JU-JLq .GT. 1) Then
        JM = (JU+JLq) / 2
        If (tt .GE. TV(JM)) Then
            JLq = JM
          Else
            JU = JM
          Endif
          Goto 13
       Endif

      If     (JLq .LE. 0) Then
         Jq = 0
      Elseif (JLq .LE. Nt-2) Then
C                                  keep q in the middle, as shown above
         Jq = JLq - 1
      Else
C                         JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq.
        Jq = Nt - 3

      Endif
C                                 This is the interpolation variable in Q
      Call Polint4F (TV(jq), AlsCTEQ(jq), tt, Alsout)
      
      CT18Alphas = Alsout

      Return
C                                       ********************
      End
