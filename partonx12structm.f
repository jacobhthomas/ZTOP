c     QUARKS FOR IPARTON
c     1  2  3  4  5  6
c     u  d  s  c  b  t

c     This program makes use of a few cuda routines
      Subroutine SetCT18(Tablefile)    
      Implicit Double Precision (A-H,O-Z)
      Character Tablefile*40
      Common /Setchange/ Isetch, ipdsset, ipdsformat
      data ipdsset, ipdsformat/0,0/
      save

      IU= NextUn()
      Open(IU, File=Tablefile, Status='OLD', Err=100)
      Call Readpds0 (IU)
      Close (IU)
      Isetch=1; ipdsset=1
      Return

 100  Print *, ' Data file ', Tablefile, ' cannot be opened '
     >  //'in SetCT18!!'
      Stop
C                             ********************
      End

      Function NextUn()
C                                 Returns an unallocated FORTRAN i/o unit.
      Logical EX
C
      Do 10 N = 10, 300
         INQUIRE (UNIT=N, OPENED=EX)
         If (.NOT. EX) then
            NextUn = N
            Return
         Endif
 10   Continue
      Stop ' There is no available I/O unit. '
C               *************************
      End


      Subroutine Readpds0 (Nu)
      Implicit Double Precision (A-H,O-Z)
      Character Line*80
      integer ipdsformat
      PARAMETER (MXX = 201, MXQ = 40, MXF = 6, MaxVal=4)
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)
      double precision qv(0:mxq)

      Common
     > / CtqPar1 / qBase,XV(0:MXX),TV(0:MXQ),UPD(MXPQX), AlsCTEQ(0:mxq)
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / XQrange / Qini, Qmax, Xmin
     > / Masstbl / Amass(6)
     > / QCDtbl /  AlfaQ, Qalfa, Ipk, Iorder, Nfl !for external use
     > /Setchange/ Isetch, ipdsset, ipdsformat

      Read  (Nu, '(A)') Line 
      Read  (Nu, '(A)') Line

      if (Line(1:11) .eq. '  ipk, Ordr') then !post-CT10 .pds format;
c Set alphas(MZ) at scale Zm, quark masses, and evolution type
        ipdsformat = 10           !Post-CT10 .pds format
        Read (Nu, *) ipk, Dr, Qalfa, AlfaQ, (amass(i),i=1,6) 
        Iorder = Nint(Dr)        
        read (Nu, '(A)') Line
        if (Line(1:7) .eq. '  IMASS' ) then
          ipdsformat = 11         !CT12 .pds format
          read (Nu, *) aimass, fswitch, N0, N0, N0, Nfmx, MxVal
          Nfl=Nfmx
        else                      !Pre-CT12 format
          Read  (Nu, *) N0, N0, N0, NfMx, MxVal
        endif                     !Line(1:7)
        
      else                        !old .pds format;      
        ipdsformat = 6            !CTEQ6.6 .pds format; alpha_s  is not specified        
        Read (Nu, *) Dr, fl, Alambda, (amass(i),i=1,6)  !set Lambda_QCD
        Iorder = Nint(Dr); Nfl = Nint(fl)

        Read  (Nu, '(A)') Line
        Read  (Nu, *) dummy,dummy,dummy, NfMx, MxVal, N0
      endif                       !Line(1:11...
      
      Read  (Nu, '(A)') Line
      Read  (Nu, *) NX,  NT, N0, NG, N0
      
      if (ng.gt.0) Read  (Nu, '(A)') (Line, i=1,ng+1)

      Read  (Nu, '(A)') Line
      if (ipdsformat.ge.11) then  !CT12 format with alpha_s values
         Read  (Nu, *) QINI, QMAX, (QV(I),TV(I), AlsCTEQ(I), I =0, NT)
      else                        !pre-CT12 format
         Read  (Nu, *) QINI, QMAX, (qv(i),TV(I), I =0, NT)
      endif                       !ipdsformat.ge.11

c check that qBase is consistent with the definition of Tv(0:nQ) for 2 values of Qv
      qbase1 = Qv(1)/Exp(Exp(Tv(1)))
      qbase2 = Qv(nT)/Exp(Exp(Tv(NT)))
      if (abs(qbase1-qbase2).gt.1e-5) then
        print *, 'Readpds0: something wrong with qbase'
        print *,'qbase1, qbase2=',qbase1,qbase2
        stop
      else
        qbase=(qbase1+qbase2)/2.0d0
      endif                     !abs(qbase1...

      Read  (Nu, '(A)') Line
      Read  (Nu, *) XMIN, aa, (XV(I), I =1, NX)
      XV(0)=0D0
      
      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+1+MxVal) ! n flavors max = nfmax, mxval = max valence
      Read  (Nu, '(A)') Line
      Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

      Return
C                        ****************************
      End

c     helper function in return pdf value
c     find where you are in x and Q
c     maybe try not caching at first
      FUNCTION PartonX12 (IPARTON, XX, QQ)

c  Given the parton distribution function in the array U in
c  COMMON / PEVLDT / , this routine interpolates to find
c  the parton distribution at an arbitray point in x and q.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (MXX = 201, MXQ = 40, MXF = 6, MaxVal=4)
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)

c New declarations      

      INTEGER*8 :: texX, texT, texPDF
      REAL*8 :: XVA(0:MXX),TVA(0:MXQ),UPDA(MXPQX)     

      Common
     > / CtqPar1 / qBase, XV(0:MXX), TV(0:MXQ),UPD(MXPQX),AlsCTEQ(0:mxq)
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / XQrange / Qini, Qmax, Xmin
     > / Setchange / Isetch, ipdsset, ipdsformat
     > / TEXTURES  / texX, texT, texPDF
      
      Dimension fvec(4), fij(4)
      Dimension xvpow(0:mxx)
      Data OneP / 1.00001 /
      Data xpow / 0.3d0 /       !**** choice of interpolation variable
      Data nqvec / 4 /
      Data ientry / 0 /
      Data X, Q, JX, JQ /-1D0, -1D0, 0, 0/
      Save xvpow
      Save X, Q, JX, JQ, JLX, JLQ
      Save ss, const1, const2, const3, const4, const5, const6
      Save sy2, sy3, s23, tt, t12, t13, t23, t24, t34, ty2, ty3
      Save tmp1, tmp2, tdet

c     get data into arrays from CUDA
c$$$      XVA  = yieldarr(1,NX+1)
c$$$      TVA  = yieldarr(2,NT+1)
c$$$      UPDA = yieldarr(3,MXPQX)

c$$$      print *, "XVA:", XVA
      
      
c      COULD BE SETUP
c      store the powers used for interpolation on first call...
      if(Isetch .eq. 1) then
         Isetch = 0

         xvpow(0) = 0D0
         do i = 1, nx
            xvpow(i) = xv(i)**xpow
         enddo
      elseIf((XX.eq.X).and.(QQ.eq.Q)) then
      	goto 99
      endif

      X = XX
      Q = QQ
      tt = log(log(Q/qBase))
c$$$      print *, tt

c      -------------    find lower end of interval containing x, i.e.,
c                       get jx such that xv(jx) .le. x .le. xv(jx+1)...
      JLx = -1
      JU = Nx+1
 
c     BIN SEARCH
      
 11   If (JU-JLx .GT. 1) Then
         JM = (JU+JLx) / 2
         If (X .Ge. XV(JM)) Then
            JLx = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif
C                     Ix    0   1   2      Jx  JLx         Nx-2     Nx
C                           |---|---|---|...|---|-x-|---|...|---|---|
C                     x     0  Xmin               x                 1


      If     (JLx .LE. -1) Then
        Print '(A,1pE12.4)','Severe error: x <= 0 in PartonX12! x = ',x
        Stop
      ElseIf (JLx .Eq. 0) Then
         Jx = 0
      Elseif (JLx .LE. Nx-2) Then

c               For interrior points, keep x in the middle, as shown above
         Jx = JLx - 1
      Elseif (JLx.Eq.Nx-1 .or. x.LT.OneP) Then

C                  We tolerate a slight over-shoot of one (OneP=1.00001),
C              perhaps due to roundoff or whatever, but not more than that.
C                                      Keep at least 4 points >= Jx
         Jx = JLx - 2
      Else
        Print '(A,1pE12.4)','Severe error: x > 1 in PartonX12! x = ',x
        Stop
      Endif
C          ---------- Note: JLx uniquely identifies the x-bin; Jx does not.

C$$$    CONFIRMED: Fort and CUDA return same x val in bin search
C     This is the variable to be interpolated in
      
      
      ss = x**xpow

      If (JLx.Ge.2 .and. JLx.Le.Nx-2) Then
         
c     initiation work for "interior bins": store the lattice points in s...
      svec1 = xvpow(jx)
      svec2 = xvpow(jx+1)
      svec3 = xvpow(jx+2)
      svec4 = xvpow(jx+3)

      s12 = svec1 - svec2
      s13 = svec1 - svec3
      s23 = svec2 - svec3
      s24 = svec2 - svec4
      s34 = svec3 - svec4

      sy2 = ss - svec2
      sy3 = ss - svec3

c constants needed for interpolating in s at fixed t lattice points...
      const1 = s13/s23
      const2 = s12/s23
      const3 = s34/s23
      const4 = s24/s23
      s1213 = s12 + s13
      s2434 = s24 + s34
      sdet = s12*s34 - s1213*s2434
      tmp = sy2*sy3/sdet
      const5 = (s34*sy2-s2434*sy3)*tmp/s12
      const6 = (s1213*sy2-s12*sy3)*tmp/s34

      
      EndIf
c
c         --------------Now find lower end of interval containing Q, i.e.,
c                         get jq such that qv(jq) .le. q .le. qv(jq+1)...
      JLq = -1
      JU = NT+1
 12   If (JU-JLq .GT. 1) Then
         JM = (JU+JLq) / 2
         If (tt .GE. TV(JM)) Then
            JLq = JM
         Else
            JU = JM
         Endif
         Goto 12
       Endif

      If     (JLq .LE. 0) Then
c$$$         print *, "Q case 1"
         Jq = 0
      Elseif (JLq .LE. Nt-2) Then
C                                  keep q in the middle, as shown above
c$$$         print *, "Q case 2"
         Jq = JLq - 1
      Else
C                         JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq.
c$$$         print *, "Q case 3"
         Jq = Nt - 3

      Endif
c$$$  CONFIRMED: Fort and CUDA Return same q bin search values
C                                   This is the interpolation variable in Q

      If (JLq.GE.1 .and. JLq.LE.Nt-2) Then
c                                        store the lattice points in t...
         
         tvec1 = Tv(jq)
         tvec2 = Tv(jq+1)
         tvec3 = Tv(jq+2)
         tvec4 = Tv(jq+3)
         
         t12 = tvec1 - tvec2
         t13 = tvec1 - tvec3
         t23 = tvec2 - tvec3
         t24 = tvec2 - tvec4
         t34 = tvec3 - tvec4

         ty2 = tt - tvec2
         ty3 = tt - tvec3

         tmp1 = t12 + t13
         tmp2 = t24 + t34

         tdet = t12*t34 - tmp1*tmp2

      EndIf


c get the pdf function values at the lattice points...

 99   If (IPARTON .Gt. MxVal) Then
         Ip = - IPARTON
      Else
         Ip = IPARTON
      EndIf
      jtmp = ((Ip + NfMx)*(NT+1)+(jq-1))*(NX+1)+jx+1

      Do it = 1, nqvec
         J1  = jtmp + it*(NX+1)

         If (Jx .Eq. 0) Then
C                      For the first 4 x points, interpolate x^2*f(x,Q)
C                      This applies to the two lowest bins JLx = 0, 1
C            We can not put the JLx.eq.1 bin into the "interrior" section
C                           (as we do for q), since Upd(J1) is undefined.
            Fij(1) = 0
            Fij(2) = Upd(J1+1) * XV(1)**2
            Fij(3) = Upd(J1+2) * XV(2)**2
            Fij(4) = Upd(J1+3) * XV(3)**2

C
C                 Use Polint which allows x to be anywhere w.r.t. the grid
c$$$            print *, "hello world 1"
         call polint_wrapper (XVpow(0), Fij(1), ss, Fx)
c$$$         call polint4f (XVpow(0), Fij(1), ss, Fx)
c$$$         print *, Fx
         
         If (x .GT. 0D0) Fvec(it) =  Fx / x**2
C                                              Pdf is undefined for x.eq.0
       ElseIf  (JLx .Eq. Nx-1) Then
C                                                This is the highest x bin:
        call polint_wrapper (XVpow(Nx-3), Upd(J1), ss, Fx)
c$$$        call polint4f (XVpow(Nx-3), Upd(J1), ss, Fx)
c$$$        print *, Fx
        
        Fvec(it) = Fx

       Else
C                       for all interior points, use Jon's in-line function
C                              This applied to (JLx.Ge.2 .and. JLx.Le.Nx-2)
          
         sf2 = Upd(J1+1)
         sf3 = Upd(J1+2)

         g1 =  sf2*const1 - sf3*const2
         g4 = -sf2*const3 + sf3*const4


         Fvec(it) = (const5*(Upd(J1)-g1)
     &               + const6*(Upd(J1+3)-g4)
     &               + sf2*sy3 - sf3*sy2) / s23

      Endif

      enddo

      

C                                   We now have the four values Fvec(1:4)
c     interpolate in t...

      If (JLq .LE. 0) Then
         
C     1st Q-bin, as well as extrapolation to lower Q
         call polint4f (TV(0), Fvec(1), tt, ff)

      ElseIf (JLq .GE. Nt-1) Then
c     Last Q-bin, as well as extrapolation to higher Q
         call polint4f (TV(Nt-3), Fvec(1), tt, ff)
      Else
C                         Interrior bins : (JLq.GE.1 .and. JLq.LE.Nt-2)
C       which include JLq.Eq.1 and JLq.Eq.Nt-2, since Upd is defined for
C                         the full range QV(0:Nt)  (in contrast to XV)
         
        tf2 = fvec(2)
        tf3 = fvec(3)

        g1 = ( tf2*t13 - tf3*t12) / t23
        g4 = (-tf2*t34 + tf3*t24) / t23

        h00 = ((t34*ty2-tmp2*ty3)*(fvec(1)-g1)/t12
     &    +  (tmp1*ty2-t12*ty3)*(fvec(4)-g4)/t34)

        ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23
      EndIf

      PartonX12 = ff

      Return
C                                       ********************
      END

      SUBROUTINE POLINT4F (XA,YA,X,Y)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C  The POLINT4 routine is based on the POLINT routine from "Numerical Recipes",
C  but assuming N=4, and ignoring the error estimation.
C  suggested by Z. Sullivan.
      DIMENSION XA(*),YA(*)

      H1=XA(1)-X
      H2=XA(2)-X
      H3=XA(3)-X
      H4=XA(4)-X

      W=YA(2)-YA(1)
      DEN=W/(H1-H2)
      D1=H2*DEN
      C1=H1*DEN

      W=YA(3)-YA(2)
      DEN=W/(H2-H3)
      D2=H3*DEN
      C2=H2*DEN

      W=YA(4)-YA(3)
      DEN=W/(H3-H4)
      D3=H4*DEN
      C3=H3*DEN

      W=C2-D1
      DEN=W/(H1-H3)
      CD1=H3*DEN
      CC1=H1*DEN

      W=C3-D2
      DEN=W/(H2-H4)
      CD2=H4*DEN
      CC2=H2*DEN

      W=CC2-CD1
      DEN=W/(H1-H4)
      DD1=H4*DEN
      DC1=H1*DEN

      If((H3+H4).lt.0D0) Then
         Y=YA(4)+D3+CD2+DD1
      Elseif((H2+H3).lt.0D0) Then
         Y=YA(3)+D2+CD1+DC1
      Elseif((H1+H2).lt.0D0) Then
         Y=YA(2)+C2+CD1+DC1
      ELSE
         Y=YA(1)+C1+CC1+DC1
      ENDIF

      RETURN
C               *************************
      END

c$$$      PROGRAM TEST
c$$$      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
c$$$      Character Tablefile*40
c$$$
c$$$      PARAMETER (MXX = 201, MXQ = 40, MXF = 6, MaxVal=4)
c$$$      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)
c$$$      double precision Alsout
c$$$      INTEGER*8 :: i, setup, texX, texT, texPDF
c$$$      INTEGER*4 :: NPTS
c$$$
c$$$      Common
c$$$     > / CtqPar1 / qBase,XV(0:MXX), TV(0:MXQ), UPD(MXPQX),AlsCTEQ(0:mxq)
c$$$     > / CtqPar2 / Nx, Nt, NfMx, MxVal
c$$$     > / XQrange / Qini, Qmax, Xmin
c$$$     > /Setchange/ Isetch, ipdsset, ipdsformat
c$$$     > / TEXTURES/ texX, texT, texPDF
c$$$
c$$$      REAL*8 :: fval, cval, XVAL, QVAL
c$$$      
c$$$      double precision QAry(5), xAry(9) 
c$$$      data QAry/ 1.5d0, 4.5d0, 10d0, 91.187d0, 200d0/,
c$$$     >  xAry /1d-5,1d-4,1d-3,1d-2,1d-1,0.3d0,0.5d0,0.7d0,0.9d0/
c$$$      
c$$$      data Tablefile /'pdf.pds'/
c$$$
c$$$c     sets up and reads data file
c$$$      call SetCT18(Tablefile)
c$$$
c$$$c     Here is the sizes associated with everything
c$$$c     | QUANTITY | XV | QV | TV | UPD   |
c$$$c     | SIZE     | NX | NT | NT | MXPQX |
c$$$c     when moving to CUDA, we need to add one to each of the sizes to maintain
c$$$      NPTS = (NX+1) * (NT+1) * (NfMx+1+MxVal)
c$$$
c$$$c     for now we are gonna just do VX QV and UPD
c$$$      texX = setup(XV,NX+1)
c$$$      texT = setup(TV,NT+1)
c$$$      texPDF = setup(UPD,NPTS)
c$$$
c$$$      
c$$$c     do for iparton = -2, -1, 0, 1, 2, 3, 4, 5
c$$$c     x = 10^-6, 0.001, 0.01, 0.1, 0.3, 0.95
c$$$c     Q = 0.4, 1.2, 1.5, 5, 10, 100, 1000
c$$$
c$$$      xval = 0.3
c$$$      qval = 1000
c$$$      
c$$$      print *, '***************************'
c$$$      print *, "x = ", xval
c$$$      print *, "q = ", qval
c$$$      print *, '***************************'
c$$$      print *, ""
c$$$
c$$$      DO jpart=-2, 5
c$$$         
c$$$         print *, "iParton = ", jpart
c$$$         print *, "************************************"
c$$$         fval =  partonx12(jpart,xval,qval)
c$$$         print *, "FORT RETURNED: ", fval
c$$$         print *, ""                  
c$$$         iSetch = 1
c$$$         call partonx12_wrapper(iSetch,jpart,nx,nt,npts,nfmx,mxval,
c$$$     &        xval,qval,qBase,texX,texT,texPDF,cval)
c$$$         
c$$$         print *, "CUDA RETURNED: ", cval
c$$$         print *, ""
c$$$         print *, "DIFFERENCE: ", fval-cval
c$$$         print *, "************************************"
c$$$      END DO
c$$$c$$$      END DO
c$$$c$$$      call free_cuda_memory()
c$$$
c$$$      END PROGRAM
      
c======================================================================
      subroutine structm(x,q,uv,dv,us,ds,s,c,b,t,g)
c======================================================================
c     *****************************************************************
c     * This is an interface for the CTEQ4/5/6 tables functions:      *
c     *  structm4, structm5, structm6
c     * This is called like structm:                                  *
c     *   call structm(x,q,u_v,d_v,u_s,d_s,s,c,b,t,g)                 *
c     *****************************************************************
C   Before using the PDF, it is necessary to do the initialization by
C       Call SetCtq(Idist,Iset)
C     Where `Idist' is 4,5, or 6 and `Iset' is from the tables.
C
C   Written by Z. Sullivan, June 2004
C
      implicit none
c
c     Local
c
      double precision x,q,uv,dv,us,ds,s,c,b,t,g
      double precision parton(8), partonx12
      integer i, jpart
      
      DO i=1, 8
         jpart = 6 - i
         parton(i)  =  partonx12(jpart,x,q)

      ENDDO
      
      t = 0d0
      b = x*parton(1)
      c = x*parton(2)
      s = x*parton(3)
      dv = x*(parton(4)-parton(8))
      ds = x*parton(8)
      uv = x*(parton(5)-parton(7))
      us = x*parton(7)
      g = x*parton(6)
c
      return
      end

