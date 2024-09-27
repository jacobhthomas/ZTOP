CZ Released by Zack Sullivan, 6/14/04
CZ
CZ HISTORY:
CZ
CZ 8/26/04 Added checks to structm5 and structm6 that set pieces to 0d0 if
CZ  if they are < 0d0.  This is to exactly match what is done on the base
CZ  CTEQ5/6, but only arises if x ~ 0.99, so is not generally an issue.  Note
CZ  that this means structm may be returning negative values for some pieces,
CZ  e.g. d_v = -d_s, but d_v+d_s >= 0.  So unless different rounding schemes
CZ  are used in different places, this should be OK.  This is not a bug.
CZ
CZ
CZ This file includes: structm5.f, structm6.f, structmn.f, spftopdg.f,
CZ     alphas2.f, cteqbase.f, cteqpdfset.f, zpolint3.f, polint4.f
CZ
CZ Note: You may have to comment out checks of Q^2 range to work with HERWIG.
CZ
CZ This has been updated by Z. Sullivan on 2/1/02 to include CTEQ6.
CZ   Updated on 2/14/02 to include fixed CTEQ6 PartonX6 routine.
CZ   Updated on 9/20/02 to include cteq6L1, with proper LO alphas for LO PDF.
CZ   Updated on 5/14/03 to include cteq61.
CZ
CZ This has been modified by Z. Sullivan on 3/27/00 to work for either
CZ  CTEQ4 or CTEQ5.
CZ
CZ These COMMON blocks were added.
CZ      common / SETCTAL / setctcl
CZ      common / CTEQDIST / idistz
CZ      common / CtqParz1 / Al6, XV6, TV, UPD6   ! For CTEQ6
CZ
CZ
CZ Change the path on line beginning with CZPATH to where you put the tables
CZ e.g.,      Tablefile='/home/zack/lib/pdf/tables/' // Flnm
CZ
CZ First set up the PDFs with SetCtq(n,I); n=4/5/6 ; I=iset
CZ Then call PDFs with CtqPDF(Iparton, X, Q)
CZ   (Adaptors exist for SetCtq4, SetCtq5, SetCtq6 and
CZ                       Ctq4Pdf, Ctq5Pdf, Ctq6Pdf.)
CZ
CZ The CTEQ4/5/6 headers follow:
CZ
C============================================================================
C                CTEQ Parton Distribution Functions: Version 4.6
C                             June 21, 1996
C                   Modified: 10/17/96, 1/7/97, 1/15/97
C                             2/17/97, 2/21/97
C                   Last Modified on April 2, 1997
C
C   Ref[1]: "IMPROVED PARTON DISTRIBUTIONS FROM GLOBAL ANALYSIS OF RECENT DEEP
C         INELASTIC SCATTERING AND INCLUSIVE JET DATA"
C   By: H.L. Lai, J. Huston, S. Kuhlmann, F. Olness, J. Owens, D. Soper
C       W.K. Tung, H. Weerts
C       Phys. Rev. D55, 1280 (1997)
C
C   Ref[2]: "CHARM PRODUCTION AND PARTON DISTRIBUTIONS"
C   By: H.L. Lai and W.K. Tung
C       MSU-HEP-61222, CTEQ-622, e-Print Archive: hep-ph/9701256
C       to appear in Z. Phys.
C
C   This package contains 13 sets of CTEQ4 PDF's. Details are:
C ---------------------------------------------------------------------------
C  Iset   PDF        Description       Alpha_s(Mz)  Lam4  Lam5   Table_File
C ---------------------------------------------------------------------------
C Ref[1]
C   1    CTEQ4M   Standard MSbar scheme   0.116     298   202    cteq4m.tbl
C   2    CTEQ4D   Standard DIS scheme     0.116     298   202    cteq4d.tbl
C   3    CTEQ4L   Leading Order           0.132     236   181    cteq4l.tbl
C   4    CTEQ4A1  Alpha_s series          0.110     215   140    cteq4a1.tbl
C   5    CTEQ4A2  Alpha_s series          0.113     254   169    cteq4a2.tbl
C   6    CTEQ4A3            ( same as CTEQ4M )
C   7    CTEQ4A4  Alpha_s series          0.119     346   239    cteq4a4.tbl
C   8    CTEQ4A5  Alpha_s series          0.122     401   282    cteq4a5.tbl
C   9    CTEQ4HJ  High Jet                0.116     303   206    cteq4hj.tbl
C   10   CTEQ4LQ  Low Q0                  0.114     261   174    cteq4lq.tbl
C ---------------------------------------------------------------------------
C Ref[2]
C   11   CTEQ4HQ  Heavy Quark             0.116     298   202    cteq4hq.tbl
C   12   CTEQ4HQ1 Heavy Quark:Q0=1,Mc=1.3 0.116     298   202    cteq4hq1.tbl
C        (Improved version of CTEQ4HQ, recommended)
C   13   CTEQ4F3  Nf=3 FixedFlavorNumber  0.106     (Lam3=385)   cteq4f3.tbl
C   14   CTEQ4F4  Nf=4 FixedFlavorNumber  0.111     292   XXX    cteq4f4.tbl
C ---------------------------------------------------------------------------
C   
C   The available applied range is 10^-5 < x < 1 and 1.6 < Q < 10,000 (GeV) 
C   except CTEQ4LQ(4HQ1) for which Q starts at a lower value of 0.7(1.0) GeV.  
C   Lam5 (Lam4, Lam3) represents Lambda value (in MeV) for 5 (4,3) flavors. 
C   The matching alpha_s between 4 and 5 flavors takes place at Q=5.0 GeV,  
C   which is defined as the bottom quark mass, whenever it can be applied.
C
C   The Table_Files are assumed to be in the working directory.
C   
C   Before using the PDF, it is necessary to do the initialization by
CZ       Call SetCtq(4,Iset)  or
C        Call SetCtq4(Iset)
C   where Iset is the desired PDF specified in the above table.
C   
CZ   The function CtqPdf (Iparton, X, Q)
C   The function Ctq4Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton] 
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
C      whereas CTEQ4F3 has, by definition, only 3 flavors and gluon;
C              CTEQ4F4 has only 4 flavors and gluon.
C   
C   For detailed information on the parameters used, e.q. quark masses, 
C   QCD Lambda, ... etc.,  see info lines at the beginning of the 
C   Table_Files.
C
C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single 
C   precision.
C   
C   If you have detailed questions concerning these CTEQ4 distributions, 
C   or if you find problems/bugs using this package, direct inquires to 
C   Hung-Liang Lai(Lai_H@pa.msu.edu) or Wu-Ki Tung(Tung@pa.msu.edu).
C   
C===========================================================================
C============================================================================
C                CTEQ Parton Distribution Functions: Version 5.0
C                             Nov. 1, 1999
C
C   Ref: "GLOBAL QCD ANALYSIS OF PARTON STRUCTURE OF THE NUCLEON:
C         CTEQ5 PPARTON DISTRIBUTIONS"
C
C  hep-ph/9903282; to be published in Eur. Phys. J. C 1999.
C
C  These PDF's use quadratic interpolation of attached tables. A parametrized 
C  version of the same PDF's without external tables is under construction.  
C  They will become available later.
C
C   This package contains 7 sets of CTEQ5 PDF's; plus two updated ones.
C   The undated CTEQ5M1 and CTEQHQ1 use an improved evolution code.
C   Both the original and the updated ones fit current data with comparable
C   accuracy.  The CTEQHQ1 set also involve a different choice of scale,
C   hence differs from CTEQHQ slightly more.  It is preferred over CTEQ5HQ.
C
C   Details are:
C ---------------------------------------------------------------------------
C  Iset   PDF        Description       Alpha_s(Mz)  Lam4  Lam5   Table_File
C ---------------------------------------------------------------------------
C   1    CTEQ5M   Standard MSbar scheme   0.118     326   226    cteq5m.tbl
C   2    CTEQ5D   Standard DIS scheme     0.118     326   226    cteq5d.tbl
C   3    CTEQ5L   Leading Order           0.127     192   146    cteq5l.tbl
C   4    CTEQ5HJ  Large-x gluon enhanced  0.118     326   226    cteq5hj.tbl
C   5    CTEQ5HQ  Heavy Quark             0.118     326   226    cteq5hq.tbl
C   6    CTEQ5F3  Nf=3 FixedFlavorNumber  0.106     (Lam3=395)   cteq5f3.tbl
C   7    CTEQ5F4  Nf=4 FixedFlavorNumber  0.112     309   XXX    cteq5f4.tbl
C         --------------------------------------------------------
C   8    CTEQ5M1  Improved CTEQ5M         0.118     326   226    cteq5m1.tbl
C   9    CTEQ5HQ1 Improved CTEQ5HQ        0.118     326   226    ctq5hq1.tbl
C ---------------------------------------------------------------------------
C   
C  The available applied range is 10^-5 << x << 1 and 1.0 << Q << 10,000 (GeV).
C   Lam5 (Lam4, Lam3) represents Lambda value (in MeV) for 5 (4,3) flavors. 
C   The matching alpha_s between 4 and 5 flavors takes place at Q=4.5 GeV,  
C   which is defined as the bottom quark mass, whenever it can be applied.
C
C   The Table_Files are assumed to be in the working directory.
C   
C   Before using the PDF, it is necessary to do the initialization by
CZ       Call SetCtq(5,Iset)  or
C        Call SetCtq5(Iset)
C   where Iset is the desired PDF specified in the above table.
C   
CZ   The function CtqPdf (Iparton, X, Q)
C   The function Ctq5Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton] 
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
C      whereas CTEQ5F3 has, by definition, only 3 flavors and gluon;
C              CTEQ5F4 has only 4 flavors and gluon.
C   
C   For detailed information on the parameters used, e.q. quark masses, 
C   QCD Lambda, ... etc.,  see info lines at the beginning of the 
C   Table_Files.
C
C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single 
C   precision.
C   
C   If you have detailed questions concerning these CTEQ5 distributions, 
C   or if you find problems/bugs using this package, direct inquires to 
C   Hung-Liang Lai(lai@phys.nthu.edu.tw) or Wu-Ki Tung(Tung@pa.msu.edu).
C   
C===========================================================================
C============================================================================
C                             April 10, 2002, v6.01
C                             February 23, 2003, v6.1
C
C   Ref[1]: "New Generation of Parton Distributions with Uncertainties from
C            Global QCD Analysis"
C       By: J. Pumplin, D.R. Stump, J.Huston, H.L. Lai, P. Nadolsky, W.K. Tung
C       JHEP 0207:012(2002), hep-ph/0201195
C
C   Ref[2]: "Inclusive Jet Production, Parton Distributions, and the Search for
C            New Physics"
C       By : D. Stump, J. Huston, J. Pumplin, W.K. Tung, H.L. Lai, S. Kuhlmann,
C             J. Owens, hep-ph/0303013
C
C   This package contains
C   (1) 4 standard sets of CTEQ6 PDF's (CTEQ6M, CTEQ6D, CTEQ6L, CTEQ6L1) ;
C   (2) 40 up/down sets (with respect to CTEQ6M) for uncertainty studies from
C        Ref[1];
C   (3) updated version of the above: CTEQ6.1M and its 40 up/down eigenvector
C        sets from Ref[2].
C
C  The CTEQ6.1M set provides a global fit that is almost equivalent in every
C  respect to the published CTEQ6M, Ref[1], although some parton distributions
C  (e.g., the gluon) may deviate from CTEQ6M in some kinematic ranges by
C  amounts that are well within the specified uncertainties.
C  The more significant improvements of the new version are associated with
C  some of the 40 eigenvector sets, which are made more symmetrical and
C  reliable in (3), compared to (2).
C
C  Details about calling convention are:
C ---------------------------------------------------------------------------
C  Iset   PDF-set     Description       Alpha_s(Mz)**Lam4  Lam5   Table_File
C ===========================================================================
C Standard, "best-fit", sets:
C --------------------------
C   1    CTEQ6M   Standard MSbar scheme   0.118     326   226    cteq6m.tbl
C   2    CTEQ6D   Standard DIS scheme     0.118     326   226    cteq6d.tbl
C   3    CTEQ6L   Leading Order           0.118**   326** 226    cteq6l.tbl
C   4    CTEQ6L1  Leading Order           0.130**   215** 165    cteq6l1.tbl
C ============================================================================
C For uncertainty calculations using eigenvectors of the Hessian:
C ---------------------------------------------------------------
C     central + 40 up/down sets along 20 eigenvector directions
C                             -----------------------------
C                Original version, Ref[1]:  central fit: CTEQ6M (=CTEQ6M.00)
C                             -----------------------
C  1xx  CTEQ6M.xx  +/- sets               0.118     326   226    cteq6m1xx.tbl
C        where xx = 01-40: 01/02 corresponds to +/- for the 1st eigenvector,
C         ... etc.
C        e.g. 100      is CTEQ6M.00 (=CTEQ6M),
C             101/102 are CTEQ6M.01/02, +/- sets of 1st eigenvector, ... etc.
C                              -----------------------
C                Updated version, Ref[2]:  central fit: CTEQ6.1M (=CTEQ61.00)
C                              -----------------------
C  2xx  CTEQ61.xx  +/- sets               0.118     326   226    ctq61.xx.tbl
C        where xx = 01-40: 01/02 corresponds to +/- for the 1st eigenvector,
C         ... etc.
C        e.g. 200      is CTEQ61.00 (=CTEQ6.1M),
C             201/202 are CTEQ61.01/02, +/- sets of 1st eigenvector, ... etc.
C ===========================================================================
C   ** ALL fits are obtained by using the same coupling strength
C   \alpha_s(Mz)=0.118 and the NLO running \alpha_s formula, except CTEQ6L1
C   which uses the LO running \alpha_s and its value determined from the fit.
C   For the LO fits, the evolution of the PDF and the hard cross sections are
C   calculated at LO.  More detailed discussions are given in the references.
C
C   The table grids are generated for 10^-6 < x < 1 and 1.3 < Q < 10,000 (GeV).
C   PDF values outside of the above range are returned using extrapolation.
C   Lam5 (Lam4) represents Lambda value (in MeV) for 5 (4) flavors.
C   The matching alpha_s between 4 and 5 flavors takes place at Q=4.5 GeV,
C   which is defined as the bottom quark mass, whenever it can be applied.
C
C   The Table_Files are assumed to be in the working directory.
C
C   Before using the PDF, it is necessary to do the initialization by
CZ       Call SetCtq(6,Iset)  or
C        Call SetCtq6(Iset) 
C   where Iset is the desired PDF specified in the above table.
C
CZ   The function CtqPdf (Iparton, X, Q)
C   The function Ctq6Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton]
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
C
C   For detailed information on the parameters used, e.q. quark masses,
C   QCD Lambda, ... etc.,  see info lines at the beginning of the
C   Table_Files.
C
C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single
C   precision.
C
C   If you have detailed questions concerning these CTEQ6 distributions,
C   or if you find problems/bugs using this package, direct inquires to
C   Pumplin@pa.msu.edu or Tung@pa.msu.edu.
C
C===========================================================================

      double precision function Ctq4Pdf (Iparton, X, Q)
      integer          iparton
      double precision x, q, ctqpdf
      ctq4pdf=ctqpdf(Iparton, x, q)
      return
      end

      double precision function Ctq5Pdf (Iparton, X, Q)
      integer          iparton
      double precision x, q, ctqpdf
      ctq5pdf=ctqpdf(Iparton, x, q)
      return
      end

      double precision function Ctq6Pdf (Iparton, X, Q)
      integer          iparton
      double precision x, q, ctqpdf
      ctq6pdf=ctqpdf(Iparton, x, q)
      return
      end

      Function CtqPdf (Iparton, X, Q)
      Implicit Double Precision (A-H,O-Z)
      Logical Warn
      Common
     > / CtqPar2 / Nx, Nt, NfMx
     > / QCDtable /  Alambda, Nfl, Iorder
      integer idist
      common / CTEQDIST / idist

      Data Warn /.true./
      save Warn

      If (X .lt. 0D0 .or. X .gt. 1D0) Then
	Print *, 'X out of range in CtqPdf: ', X
	Stop
      Endif
      If (Q .lt. Alambda) Then
	Print *, 'Q out of range in CtqPdf: ', Q
	Stop
      Endif
      If ((Iparton .lt. -NfMx .or. Iparton .gt. NfMx)) Then
         If (Warn) Then
C        put a warning for calling extra flavor.
	     Warn = .false.
	     Print *, 'Warning: Iparton out of range in CtqPdf: '
     >              , Iparton
         Endif
         CtqPdf = 0D0
         Return
      Endif

      if(idist.le.5) then
         CtqPdf = PartonX5 (Iparton, X, Q)
      else                                       ! idist=6
         CtqPdf = PartonX6 (Iparton, X, Q)
      endif
      if(CtqPdf.lt.0.D0)  CtqPdf = 0.D0

      Return

C                             ********************
      End

      Subroutine SetCtq (IDist,Iset)
      Implicit Double Precision (A-H,O-Z)
CZS Logical added to initialize alphas.
      logical setctcl
      common / SETCTAL / setctcl
      integer idistz
      common / CTEQDIST / idistz
CZS
      Parameter (Iset4max=14)
      Parameter (Iset5max=9)
      Parameter (Iset6max0=5)
cz      Parameter (Iset6max0=4)
CZS      Character Flnm(Isetmax)*12, Tablefile*40
      Character Flnm4(Iset4max)*12, Tablefile*50
      Character Flnm5(Iset5max)*12,Flnm*14
      Character Flnm6(Iset6max0)*6, nn*3
      Data (Flnm4(I), I=1,Iset4max)
     > / 'cteq4m.tbl', 'cteq4d.tbl', 'cteq4l.tbl'
     > , 'cteq4a1.tbl', 'cteq4a2.tbl', 'cteq4m.tbl', 'cteq4a4.tbl'
     > , 'cteq4a5.tbl', 'cteq4hj.tbl', 'cteq4lq.tbl'
     > , 'cteq4hq.tbl', 'cteq4hq1.tbl', 'cteq4f3.tbl', 'cteq4f4.tbl' /
      Data (Flnm5(I), I=1,Iset5max)
     > / 'cteq5m.tbl', 'cteq5d.tbl', 'cteq5l.tbl', 'cteq5hj.tbl'
     > , 'cteq5hq.tbl', 'cteq5f3.tbl', 'cteq5f4.tbl'
     > , 'cteq5m1.tbl', 'ctq5hq1.tbl'  /
      Data (Flnm6(I), I=1,Iset6max0)
     > / 'cteq6m', 'cteq6d', 'cteq6l', 'cteq6l', 'ctq61.' /
      Data Tablefile / 'test.tbl' /
      Data Isetold, Isetmin, Isettest / -987, 1, 911 /
      Data Iset6min0, Iset6min1, Iset6max1 /1,100,140/
      Data Iset6min2,Iset6max2 /200,240/
      save

C             If data file not initialized, do so.
CZ      If(Iset.ne.Isetold) then
      If(((IDist-4)*15+Iset).ne.Isetold) then
	 IU= NextUn5()
         If (Iset .eq. Isettest) then
            Print* ,'Opening ', Tablefile
 21         Open(IU, File=Tablefile, Status='OLD', Err=101)
            GoTo 22
 101        Print*, Tablefile, ' cannot be opened '
            Print*, 'Please input the .tbl file:'
            Read (*,'(A)') Tablefile
            Goto 21
 22         Continue
CZ         ElseIf (Iset.lt.Isetmin .or. Iset.gt.Isetmax) Then
         ElseIf (Iset.lt.Isetmin .or. ((Iset.gt.Iset4max).and.
     >           (Idist.eq.4)) .or. ((Iset.gt.Iset5max).and.
     >           (Idist.eq.5)) .or. Idist.lt.4 .or. Idist.gt.6
     >           .or. ((Idist.eq.6).and.((Iset.lt.Iset6min0).or.
     >           ((Iset.gt.Iset6max0).and.(Iset.lt.Iset6min1)).or.
     >           ((Iset.gt.Iset6max1).and.(Iset.lt.Iset6min2)).or.
     >           (Iset.gt.Iset6max2)))) Then
	    Print *, 'Invalid Iset number in SetCtq :', Iset
	    Stop
         Else
C
CZS Local directory where tables reside. Modified 3/23/00, ZS.
C
            if (Idist.eq.4) Flnm=Flnm4(Iset)
            if (Idist.eq.5) Flnm=Flnm5(Iset)
            If (Idist.eq.6) then
               if (Iset.ge.Iset6min0 .and. Iset.le.3) Then
                  Flnm=Flnm6(Iset)//'.tbl'
               elseif (Iset.eq.4) Then
                  Flnm=Flnm6(Iset)//'1.tbl'
               elseif (Iset.eq.100) Then
                  Flnm=Flnm6(1)//'.tbl'
               Elseif (Iset.ge.Iset6min1 .and. Iset.le.Iset6max1) Then
                  write(nn,'(I3)') Iset
                  Flnm=Flnm6(1)//nn//'.tbl'
               Elseif (Iset.ge.Iset6min2 .and. Iset.le.Iset6max2) Then
                  write(nn,'(I3)') Iset
                  Flnm=Flnm6(5)//nn(2:3)//'.tbl'
               endif
            endif
            Tablefile=Flnm
CZPATH            Tablefile='/home/zack/lib/pdf/tables/' // Flnm
            Open(IU, File=Tablefile, Status='OLD', Err=100)
	 Endif
         idistz=idist
         Call ReadTbl5 (IU)
         Close (IU)
CZ	 Isetold=Iset
	 Isetold=(IDist-4)*15+Iset
         write (*,'(a10,i1,a9,a14)') 'Using CTEQ',Idist,', table: ',
     >        Flnm
      Endif
CZ Initialize alphas2 for the new set.
      setctcl = .true.
      aldummy = alphas2(91.187d0)
      setctcl = .false.
      Return

 100  Print *, ' Data file ', Tablefile, ' cannot be opened '
     >//'in SetCtq!!'
      Stop
C                             ********************
      End

      Subroutine SetCtq4 (Iset)
      integer iset
      call SetCtq(4,Iset)
      return
      end

      Subroutine SetCtq5 (Iset)
      integer iset
      call SetCtq(5,Iset)
      return
      end

      Subroutine SetCtq6 (Iset)
      integer iset
      call SetCtq(6,Iset)
      return
      end

      Subroutine ReadTbl5 (Nu)
      Implicit Double Precision (A-H,O-Z)
      Character Line*80
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      PARAMETER (MXX6 = 96, MXQ6 = 20, MXF6 = 5)
      PARAMETER (MXPQX6 = (MXF6 + 3) * MXQ6 * MXX6)
      Common 
     > / CtqPar1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / CtqParz1 / Al6, XV6(0:MXX6), TV(0:MXQ6), UPD6(MXPQX6)
     > / CtqPar2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin
     > / QCDtable /  Alambda, Nfl, Iorder
     > / Masstbl / Amass(6)
      integer idist
      common / CTEQDIST / idist

      Read  (Nu, '(A)') Line     
      Read  (Nu, '(A)') Line
      Read  (Nu, *) Dr, Fl, Al, (Amass(I),I=1,6)
      al6=al
      Iorder = Nint(Dr)
      Nfl = Nint(Fl)
      Alambda = Al

      Read  (Nu, '(A)') Line 
      Read  (Nu, *) NX,  NT, NfMx

      Read  (Nu, '(A)') Line
cz      Read  (Nu, *) QINI, QMAX, (QL(I), I =0, NT)
      if(idist.le.5) Read  (Nu, *) QINI, QMAX, (QL(I), I =0, NT)
      if(idist.eq.6) Read  (Nu, *) QINI, QMAX, (TV(I), I =0, NT)

      Read  (Nu, '(A)') Line
cz      Read  (Nu, *) XMIN, (XV(I), I =0, NX)
      if(idist.le.5) Read  (Nu, *) XMIN, (XV(I), I =0, NX)
      if(idist.eq.6) Read  (Nu, *) XMIN, (XV6(I), I =0, NX)

      if(idist.le.5)then
      Do 11 Iq = 0, NT
         QL(Iq) = DLog (QL(Iq) /Al)
cz         QL(Iq) = Log (QL(Iq) /Al)
   11 Continue
      else                     ! idist = 6
      Do 16 Iq = 0, NT
         TV(Iq) = DLog(DLog (TV(Iq) /Al6))
   16 Continue
      endif
C
C                  Since quark = anti-quark for nfl>2 at this stage, 
C                  we Read  out only the non-redundent data points
C     No of flavors = NfMx (sea) + 1 (gluon) + 2 (valence) 

      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+3)
      Read  (Nu, '(A)') Line
      if(idist.le.5) Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)
      if(idist.eq.6) Read  (Nu, *, IOSTAT=IRET) (UPD6(I), I=1,Npts)

      Return
C                        ****************************
      End

      Function NextUn5()
C                                 Returns an unallocated FORTRAN i/o unit.
      Logical EX
C
      Do 10 N = 10, 300
         INQUIRE (UNIT=N, OPENED=EX)
         If (.NOT. EX) then
            NextUn5 = N
            Return
         Endif
 10   Continue
      Stop ' There is no available I/O unit. '
C               *************************
      End
C


C======================================================================
      DOUBLE PRECISION FUNCTION ALPHAS2(Q)
C======================================================================
c     Written by Z. Sullivan, March 30, 2000
c     Based on Collins, Tung NPB 278 (86) 934, eq. 5.
c     This is an interface for the 1 or 2-loop alpha_s to go with
c      the CTEQ 4/5/6 PDf's.
c     NOTE: This is reinitialized when SetCtq is called.
c     This is like the function alphas2 from CERNLIB:
c        a = alphas2(q)
c
      implicit none
      double precision q

      double precision pi
      parameter (pi = 3.1415926535897931d0)
      double precision Qini, Qmax, Xmin, Alambda, Amass(6)
      integer Nfl,Iorder
      common / XQrange / Qini, Qmax, Xmin
      common / QCDtable /  Alambda, Nfl, Iorder
      common / Masstbl / Amass
      logical setctcl
      common / SETCTAL / setctcl
      double precision almu,allamin
      integer alor,alnf,alwhch
      common /ZALCTFNC/ almu,allamin,alor,alnf,alwhch
      DOUBLE PRECISION QCDL4,QCDL5
      COMMON/W50512/ QCDL4,QCDL5
      double precision lam3,lam4,lam5,lam6
      double precision beta1,beta2,nf,xlam,lql
      double precision alrtbis
      external alrtbis

      LOGICAL INIT,FIX

      SAVE INIT,FIX,lam3,lam4,lam5,lam6
      DATA INIT /.FALSE./
      DATA FIX  /.FALSE./
      IF ((.NOT.INIT).or.(setctcl)) THEN
c     ** lamqcd(4) = 303 MeV for mZ=91.187, als=0.117 **
c     CTEQ PDF's fix lambda_5 in the table, unless it is a FFS.
         if (Nfl .ne. 5) then
            fix = .true.
            if (Nfl .eq. 3) then
               lam3 = Alambda
               lam4 = 1d8
            else
               lam4 = Alambda
               lam3 = 1d-8
            endif
            lam5=1d8
            lam6=1d8
            goto 10
         endif
         lam5 = Alambda     ! VFN scheme.
c     Compute lam6,lam4,lam3
         alor = Iorder
         almu = Amass(6)
         alnf = 5
         allamin = lam5
         alwhch = 1      ! n+1 from n
         lam6 = alrtbis(1d-3,9d-1,1d-5)
c
         almu = Amass(5)
         alnf = 5
         allamin = lam5
         alwhch = -1      ! n-1 from n
         lam4 = alrtbis(1d-3,9d-1,1d-5)
c
         almu = Amass(4)
         alnf = 4
         allamin = lam4
         alwhch = -1      ! n-1 from n
         lam3 = alrtbis(1d-3,9d-1,1d-5)
c
 10      INIT=.TRUE.

         QCDL4=lam4
         QCDL5=lam5
c         print *,lam3,lam4,lam5,lam6,Iorder,Nfl
      ENDIF

      IF (q .gt. amass(6)) then
         nf = 6d0
         xlam = lam6
      elseif (q .gt. amass(5)) then
         nf = 5d0
         xlam = lam5
      elseif (q .gt. amass(4)) then
         nf = 4d0
         xlam = lam4
      elseif (q .gt. amass(3)) then
         nf = 3d0
         xlam = lam3
      endif
      if (fix) then
         nf = Nfl
         xlam = Alambda
      endif

      if (Iorder .eq. 1) then
         alphas2 = 12d0*pi/(33d0-2d0*nf)/log(q*q/xlam/xlam)
         RETURN
      endif            ! Otherwise 2-loop

      beta1 = 12d0*pi/(33d0-2d0*nf)
      beta2 = 24d0*pi*pi/(153d0-19d0*nf)
      lql = dlog(q*q/xlam/xlam)
      alphas2 = beta1/lql*(1d0 - beta1*beta1/beta2*dlog(lql)/lql)

      RETURN
      END

C======================================================================
C The following are used to calculate lambda_n/n+1 for alphas2.
C======================================================================
      DOUBLE PRECISION FUNCTION alrtbis(x1,x2,xacc)
C
C Modified from Numerical Recipies rtbis to be double precision.
C
      IMPLICIT NONE
      INTEGER JMAX
      DOUBLE PRECISION x1,x2,xacc,alrtfunc
      PARAMETER (JMAX=40)
      INTEGER j
      DOUBLE PRECISION dx,f,fmid,xmid
      EXTERNAL alrtfunc
      fmid=alrtfunc(x2)
      f=alrtfunc(x1)
      if(f*fmid.ge.0.) pause 'root must be bracketed in rtbis'
      if(f.lt.0.)then
        alrtbis=x1
        dx=x2-x1
      else
        alrtbis=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*.5
        xmid=alrtbis+dx
        fmid=alrtfunc(xmid)
        if(fmid.le.0.)alrtbis=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.) return
11    continue
      pause 'too many bisections in rtbis'
      END

      double precision function alrtfunc(allamout)
      implicit none
      double precision pi
      parameter (pi = 3.1415926535897931d0)
      double precision allamout
      double precision almu,allamin
      integer alor,alnf,alwhch
      common /ZALCTFNC/ almu,allamin,alor,alnf,alwhch
      double precision beta1i,beta1o,beta2i,beta2o,lqli,lqlo

      beta1i = 12d0*pi/(33d0-2d0*dble(alnf))
      lqli = log(almu*almu/allamin/allamin)
      beta1o = 12d0*pi/(33d0-2d0*dble(alnf+alwhch))
      lqlo = log(almu*almu/allamout/allamout)

      if (alor .eq. 1) then
         alrtfunc = beta1i/lqli - beta1o/lqlo
         RETURN
      endif            ! Otherwise 2-loop

      beta2i = 24d0*pi*pi/(153d0-19d0*dble(alnf))
      beta2o = 24d0*pi*pi/(153d0-19d0*dble(alnf+alwhch))

      alrtfunc = beta1i/lqli*(1d0-beta1i*beta1i/beta2i*log(lqli)/lqli)
     &     - beta1o/lqlo*(1d0-beta1o*beta1o/beta2o*log(lqlo)/lqlo)

      RETURN
      end


cz Modified partonx5 to save variables between runs.  This removes common
cz   code between iterations in structm/etc. - i.e. x1=x0, Q1=Q0, u->d etc.
cz   There is no reason to recalculate many parameters if x, Q do not change.
cz Also specified a call to POLINT3 rather than POLINT.
cz
cz Rewritten by Zack Sullivan, May, 2004
cz
      FUNCTION PartonX5 (IPRTN, X, Q)
C
C   Given the parton distribution function in the array Upd in
C   COMMON / CtqPar1 / , this routine fetches u(fl, x, q) at any value of
C   x and q using Mth-order polynomial interpolation for x and Ln(Q/Lambda).
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      PARAMETER (M= 2, M1 = M + 1)
C
      Logical First
      Common 
     > / CtqPar1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin
C
      Dimension Fq(M1), Df(M1)

      Data First /.true./
      save First

cz Simple speedup for structm, etc. 5/2/04
      double precision xlast, qlast, qg
      integer jx,jq
      data xlast, qlast / -1d0, -1d0 /
      data jx, jq / 0, 0 /
      save xlast,qlast,qg,jx,jq

      if ((x.eq.xlast).and.(q.eq.qlast)) goto 99
      xlast=x
      qlast=q

C                                                 Work with Log (Q)
      QG  = LOG (Q/AL)

C                           Find lower end of interval containing X
      JL = -1
      JU = Nx+1
 11   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (X .GT. XV(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif

      Jx = JL - (M-1)/2
      If (X .lt. Xmin .and. First ) Then
         First = .false.
         Print '(A, 2(1pE12.4))', 
     >     ' WARNING: X << Xmin, extrapolation used; X, Xmin =', X, Xmin
         If (Jx .LT. 0) Jx = 0
      Elseif (Jx .GT. Nx-M) Then
         Jx = Nx - M
      Endif
C                                    Find the interval where Q lies
      JL = -1
      JU = NT+1
 12   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (QG .GT. QL(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 12
      Endif

      Jq = JL - (M-1)/2
      If (Jq .LT. 0) Then
         Jq = 0
         If (Q .lt. Qini)  Print '(A, 2(1pE12.4))', 
     >     ' WARNING: Q << Qini, extrapolation used; Q, Qini =', Q, Qini
      Elseif (Jq .GT. Nt-M) Then
         Jq = Nt - M
         If (Q .gt. Qmax)  Print '(A, 2(1pE12.4))', 
     >     ' WARNING: Q > Qmax, extrapolation used; Q, Qmax =', Q, Qmax
      Endif

 99   If (Iprtn .GE. 3) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
C                             Find the off-set in the linear array Upd
      JFL = Ip + NfMx
      J0  = (JFL * (NT+1) + Jq) * (NX+1) + Jx
C
C                                           Now interpolate in x for M1 Q's
      Do 21 Iq = 1, M1
         J1 = J0 + (Nx+1)*(Iq-1) + 1
         Call Polint3 (XV(Jx), Upd(J1), M1, X, Fq(Iq), Df(Iq))
 21   Continue
C                                          Finish off by interpolating in Q
      Call Polint3 (QL(Jq), Fq(1), M1, QG, Ftmp, Ddf)

      PartonX5 = Ftmp
C
      RETURN
C                        ****************************
      END


cz Modified partonx6 to save variables between runs.  This removes common
cz   code between iterations in structm/etc. - i.e. x1=x0, Q1=Q0, u->d etc.
cz   There is no reason to recalculate many parameters if x, Q do not change.
cz Also specified a call to POLINT4 rather than POLINT.
cz
cz Rewritten by Zack Sullivan, May, 2004
cz
      Function PartonX6 (IPRTN, XX, QQ)

c  Given the parton distribution function in the array U in
c  COMMON / PEVLDT / , this routine interpolates to find 
c  the parton distribution at an arbitray point in x and q.
c
      Implicit Double Precision (A-H,O-Z)
 
      Parameter (MXX = 96, MXQ = 20, MXF = 5)
      Parameter (MXQX= MXQ * MXX,   MXPQX = MXQX * (MXF+3))
 
      Common
     > / CtqParz1 / Al, XV(0:MXX), TV(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin
 
      Dimension fvec(4), fij(4)
      Dimension xvpow(0:mxx)
      Data OneP / 1.00001 /
      Data xpow / 0.3d0 /	!**** choice of interpolation variable
      Data nqvec / 4 /
      Data ientry / 0 /
      Save ientry,xvpow

cz Simple speedup for structm, etc. 5/2/04
      double precision x, q
      integer jx, jq
      data x, q / -1d0, -1d0 /
      data jx, jq / 0, 0 /
      save x, q, jx, jq, jlx, jlq
      save ss, const1, const2, const3, const4, const5, const6
      save sy2, sy3, s23, tt, t12, t13, t23, t24, t34, ty2, ty3
      save tmp1, tmp2, tdet

      if ((xx.eq.x).and.(qq.eq.q)) goto 99

c store the powers used for interpolation on first call...
      if(ientry .eq. 0) then
         ientry = 1

         xvpow(0) = 0D0
         do i = 1, nx
            xvpow(i) = xv(i)**xpow
         enddo
      endif

      X = XX
      Q = QQ
      tt = log(log(Q/Al))

c      -------------    find lower end of interval containing x, i.e.,
c                       get jx such that xv(jx) .le. x .le. xv(jx+1)...
      JLx = -1
      JU = Nx+1
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
C
      If     (JLx .LE. -1) Then
        Print '(A,1pE12.4)', 'Severe error: x <= 0 in PartonX6! x = ', x 
        Stop 
      ElseIf (JLx .Eq. 0) Then
         Jx = 0
      Elseif (JLx .LE. Nx-2) Then

C                For interrior points, keep x in the middle, as shown above
         Jx = JLx - 1
      Elseif (JLx.Eq.Nx-1 .or. x.LT.OneP) Then

C                  We tolerate a slight over-shoot of one (OneP=1.00001),
C              perhaps due to roundoff or whatever, but not more than that.
C                                      Keep at least 4 points >= Jx
         Jx = JLx - 2
      Else
        Print '(A,1pE12.4)', 'Severe error: x > 1 in PartonX6! x = ', x 
        Stop 
      Endif
C          ---------- Note: JLx uniquely identifies the x-bin; Jx does not.

C                       This is the variable to be interpolated in
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

c         --------------Now find lower end of interval containing Q, i.e.,
c                          get jq such that qv(jq) .le. q .le. qv(jq+1)...
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
         Jq = 0
      Elseif (JLq .LE. Nt-2) Then
C                                  keep q in the middle, as shown above
         Jq = JLq - 1
      Else
C                         JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq.
        Jq = Nt - 3

      Endif
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

 99   If (Iprtn .GE. 3) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
      jtmp = ((Ip + NfMx)*(NT+1)+(jq-1))*(NX+1)+jx+1

      Do it = 1, nqvec

         J1  = jtmp + it*(NX+1) 

       If (Jx .Eq. 0) Then
C                          For the first 4 x points, interpolate x^2*f(x,Q)
C                           This applies to the two lowest bins JLx = 0, 1
C            We can not put the JLx.eq.1 bin into the "interrior" section
C                           (as we do for q), since Upd(J1) is undefined.
         fij(1) = 0
         fij(2) = Upd(J1+1) * XV(1)**2
         fij(3) = Upd(J1+2) * XV(2)**2
         fij(4) = Upd(J1+3) * XV(3)**2
C 
C                 Use Polint which allows x to be anywhere w.r.t. the grid

         Call Polint4 (XVpow(0), Fij(1), 4, ss, Fx, Dfx) 
         
         If (x .GT. 0D0)  Fvec(it) =  Fx / x**2 
C                                              Pdf is undefined for x.eq.0
       ElseIf  (JLx .Eq. Nx-1) Then
C                                                This is the highest x bin:

        Call Polint4 (XVpow(Nx-3), Upd(J1), 4, ss, Fx, Dfx)

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
C                         1st Q-bin, as well as extrapolation to lower Q
        Call Polint4 (TV(0), Fvec(1), 4, tt, ff, Dfq)

      ElseIf (JLq .GE. Nt-1) Then
C                         Last Q-bin, as well as extrapolation to higher Q
        Call Polint4 (TV(Nt-3), Fvec(1), 4, tt, ff, Dfq)
      Else
C                         Interrior bins : (JLq.GE.1 .and. JLq.LE.Nt-2)
C       which include JLq.Eq.1 and JLq.Eq.Nt-2, since Upd is defined for
C                         the full range QV(0:Nt)  (in contrast to XV)
        tf2 = fvec(2)
        tf3 = fvec(3)

        g1 = ( tf2*t13 - tf3*t12) / t23
        g4 = (-tf2*t34 + tf3*t24) / t23

        h00 = ((t34*ty2-tmp2*ty3)*(fvec(1)-g1)/t12 
     &	  +  (tmp1*ty2-t12*ty3)*(fvec(4)-g4)/t34)

        ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23
      EndIf

      PartonX6 = ff

      Return
C					********************
      End


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
c
c     Global
c
      integer idist
      common / CTEQDIST / idist
c     ----------
c     Begin Code
c     ----------
      if(idist.eq.6) call structm6(x,q,uv,dv,us,ds,s,c,b,t,g)
      if(idist.le.5) call structm5(x,q,uv,dv,us,ds,s,c,b,t,g)  ! 4,5 are same
c      if(idist.eq.5) call structm5(x,q,uv,dv,us,ds,s,c,b,t,g)
c      if(idist.eq.4) call structm4(x,q,uv,dv,us,ds,s,c,b,t,g)
c
      return
      end


C   Written by Z. Sullivan, 3/26/04
C
C   This code is a PDFSET interface for use with the CTEQ 4-6 codes, and
C    cteqbase.f.  It will fill the common blocks required by HERWIG and
C    PYTHIA.
C
C   PDFSET(PARM,VAL)
C      PARM other than 'MODE' or 'NSET' are ignored.
C      VAL is parsed to determine which set is desired: SET*1000+PDFSET
C         ex) PARM='MODE', VAL=5003 is CTEQ5L, while 6202 is CTEQ6_202
C       If VAL < 1000, then standard mappings are used (to the extent known).
C         HERWIG ex) autpdf(2)='MODE'
C                    modpdf(2)=6004
C                    autpdf(1)='MODE'
C                    modpdf(1)=6004
C         PYTHIA ex) mstp(51)=4604 ! CTEQ6L1, PYTHIA will strip off leading 4
C                    mstp(52)=2    ! PDFLIB
C   All common blocks as of PDFLIB 8.0

      SUBROUTINE PDFSET(PARM,VAL)
      implicit none

      CHARACTER*20 PARM(20),STRING
      DOUBLE PRECISION VAL(20)
C Local
      integer i,idist,iset

C CTEQ common blocks
      double precision Qinict, Qmaxct, Xminct, Alambdact, Amassct(6)
      integer Nflct,Iorderct
      common / XQrange / Qinict, Qmaxct, Xminct
      common / QCDtable /  Alambdact, Nflct, Iorderct
      common / Masstbl / Amassct

C PDFLIB common blocks
      INTEGER        NPTYPE,NGROUP,NSET,MODE,NFL,LO
      DOUBLE PRECISION TMAS
      COMMON/W50511/ NPTYPE,NGROUP,NSET,MODE,NFL,LO,TMAS
      DOUBLE PRECISION QCDL4,QCDL5
      COMMON/W50512/ QCDL4,QCDL5
      DOUBLE PRECISION XMIN,XMAX,Q2MIN,Q2MAX
      COMMON/W50513/ XMIN,XMAX,Q2MIN,Q2MAX
      INTEGER        IFLPRTP,NPTYPP,NGROPP,NSETP,NFLP,LOP
      DOUBLE PRECISION TMASP,QCDL4P,QCDL5P,XMINP,XMAXP,Q2MINP,Q2MAXP
      COMMON/W50518/ IFLPRTP,
     +               NPTYPP,NGROPP,NSETP,NFLP,LOP,TMASP,
     +               QCDL4P,QCDL5P,
     +               XMINP,XMAXP,Q2MINP,Q2MAXP
      LOGICAL NEWVER
      COMMON/W50519/ NEWVER
      CHARACTER*10 PDFVER(3)
      COMMON/W505190/ PDFVER

      INTEGER ISTART, INSET
      DATA ISTART/0/
      DATA INSET/0/
      SAVE ISTART, INSET

C------------------------------------
C Begin code
C------------------------------------
      do i=1,20
         string=parm(i)
         if((string(1:4).eq.'MODE').or.(string(1:4).eq.'NSET'))
     &        nset=INT(val(i))
      enddo

      if((istart.ne.0).and.(nset.eq.inset)) then
         nptype=nptypp
         ngroup=ngropp
         nset=nsetp
         nfl=nflp
         lo=lop
         tmas=tmasp
         qcdl4=qcdl4p
         qcdl5=qcdl5p
         xmin=xminp
         xmax=xmaxp
         q2min=q2minp
         q2max=q2maxp
      else
         istart=1
         inset=nset

         nptype = 1     ! Fill W50511
         ngroup = 4
         mode = nset    ! Not standard, but allowing both versions.
         if (inset>999) then
            idist = inset/1000
            iset = inset-1000*idist
         else
            if(inset.lt.46) then       ! CTEQ4
               idist=4
               if(inset.eq.32) iset=3
               if(inset.eq.33) iset=2
               if(inset.eq.34) iset=1
               if(inset.ge.35) iset=inset-31
            elseif(inset.lt.55) then   ! CTEQ5
               idist=5
               if(inset.eq.46) iset=3
               if(inset.eq.47) iset=2
               if(inset.eq.48) iset=1
               if(inset.ge.49) iset=inset-45
            else                       ! CTEQ6
cz Added 5/3/04  CTEQ6: 600+61-64, 600+1xx, 600+2xx
               idist=6
               iset=inset-600
            endif
         endif

         call setctq(idist,iset)

         nfl=nflct         ! Fill W50511
         lo=iorderct
         tmas=amassct(6)
c W50512 filled in my alphas2(q)
         xmin=xminct       ! Fill W50513
         xmax=0.9999999d0
         q2min=qinict*qinict
         q2max=qmaxct*qmaxct
         newver=.true.     ! Fill W50519 (fake)
         pdfver(1)='8.04(cteq)'  ! Fill W505190
         pdfver(2)='2004-03-26'
         pdfver(3)='          '

         iflprtp=0         ! Fill W50518 (for safety)
         nptypp=nptype
         ngropp=ngroup
         nsetp=nset
         nflp=nfl
         lop=lo
         tmasp=tmas
         qcdl4p=qcdl4
         qcdl5p=qcdl5
         xminp=xmin
         xmaxp=xmax
         q2minp=q2min
         q2maxp=q2max
      endif
      return
      end

C The following common blocks are not used:
C      INTEGER        IFLPRT
C      COMMON/W50510/ IFLPRT
C      INTEGER         NATYPE,NAGROUP,NASET
C      COMMON/W50511A/ NATYPE,NAGROUP,NASET
C      DOUBLE PRECISION WXMIN,WXMAX,WQ2MIN,WQ2MAX
C      DOUBLE PRECISION WTXMIN,WTXMAX,WTQ2MIN,WTQ2MAX
C      COMMON/W50514/ WXMIN,WXMAX,WQ2MIN,WQ2MAX,
C     +               WTXMIN,WTXMAX,WTQ2MIN,WTQ2MAX
C      DOUBLE PRECISION PDFWGT
C      COMMON/W50514W/ PDFWGT
C      INTEGER        IFLSET,IFLSTA
C      COMMON/W50515/ IFLSET,IFLSTA
C      LOGICAL FIRST
C      COMMON/W50516/ FIRST
C      INTEGER        N6
C      COMMON/W50517/ N6
C      INTEGER         IFLPRTA,NATYPP,NAGROPP,NASETP
C      COMMON/W50518A/ IFLPRTA,
C     +               NATYPP,NAGROPP,NASETP
C      CHARACTER*8 SFNAME(NPTYMX,NGRMAX,NSETMX)
C      COMMON/W505110/ SFNAME
C      INTEGER         NPGSMX(NPTYMX,NGRMAX),NSETFL(NPTYMX,NGRMAX,NSETMX)
C      COMMON/W505120/ NPGSMX(NPTYMX,NGRMAX),NSETFL(NPTYMX,NGRMAX,NSETMX)
C      INTEGER         NPTYCR(MODEMX),NGROCR(MODEMX),NSETCR(MODEMX)
C      COMMON/W505121/ NPTYCR(MODEMX),NGROCR(MODEMX),NSETCR(MODEMX)
C      INTEGER         MODECR(NPTYMX,NGRMAX,NSETMX)
C      COMMON/W505122/ MODECR(NPTYMX,NGRMAX,NSETMX)


c======================================================================
      subroutine structm5(x,q,uv,dv,us,ds,s,c,b,t,g)
c======================================================================
c     *****************************************************************
c     * This is an interface for the CTEQ4/5 tables functions:        *
c     *  CtqPdf (Iparton, X, Q)                                       *
c     * This is called like structm:                                  *
c     *   call structm(x,q,u_v,d_v,u_s,d_s,s,c,b,t,g)                 *
c     *****************************************************************
C   Before using the PDF, it is necessary to do the initialization by
C       Call SetCtq(Idist,Iset)
C     Where `Idist' is 4 or 5, and `Iset' is from the tables.
C   *** The version ONLY works with CTEQ4/5 PDFs. ***
C
C   Written by Z. Sullivan, May 2004
C     Absorbs and optimizes PartonX5 from the base CTEQ4/5 code.
C
      implicit none

      integer mxx, mxq, mxf, mxpqx, m, m1
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      PARAMETER (M= 2, M1 = M + 1)
c
c     Local
c
      double precision x,q,uv,dv,us,ds,s,c,b,t,g
      double precision parton(8)
      integer i
c     ----------
c     Begin Code
c     ----------
cz Begin New code

C   Given the parton distribution function in the array Upd in
C   COMMON / CtqPar1 / , this routine fetches u(fl, x, q) at any value of
C   x and q using Mth-order polynomial interpolation for x and Ln(Q/Lambda).
C
c      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
C
      Logical First
      double precision Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
      integer Nx, Nt, NfMx
      double precision Qini, Qmax, Xmin
      Common 
     > / CtqPar1 / Al, XV, QL, UPD
     > / CtqPar2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin
C
      double precision Fq(M1), Df(M1)
      double precision qg,ftmp,ddf
      integer jl,ju,jm,jz,jq,iprtn,ip,jfl,j0,j1,jx,iq

      DOUBLE PRECISION C1,HO,HP,HP2,W,D1,D2,DEN,QHO,QHP,QHP2
      DOUBLE PRECISION XA(3),YA(3)
      logical branch2,branch3

      Data First /.true./
      save First
C                                                 Work with Log (Q)
      If (X .lt. 0D0 .or. X .gt. 1D0) Then
	Print *, 'X out of range in structm5: ', X
	Stop
      Endif

      QG  = DLOG (Q/AL)

C                           Find lower end of interval containing X
      JL = -1
      JU = Nx+1
 11   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (X .GT. XV(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif

c      Jx = JL - (M-1)/2
      Jx = JL
      If (X .lt. Xmin .and. First ) Then
         First = .false.
         Print '(A, 2(1pE12.4))', 
     >     ' WARNING: X << Xmin, extrapolation used; X, Xmin =', X, Xmin
         If (Jx .LT. 0) Jx = 0
      Elseif (Jx .GT. Nx-M) Then
         Jx = Nx - M
      Endif
C                                    Find the interval where Q lies
      JL = -1
      JU = NT+1
 12   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (QG .GT. QL(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 12
      Endif

c      Jq = JL - (M-1)/2
      Jq = JL
      If (Jq .LT. 0) Then
         Jq = 0
         If (Q .lt. Qini)  Print '(A, 2(1pE12.4))', 
     >     ' WARNING: Q << Qini, extrapolation used; Q, Qini =', Q, Qini
      Elseif (Jq .GT. Nt-M) Then
         Jq = Nt - M
         If (Q .gt. Qmax)  Print '(A, 2(1pE12.4))', 
     >     ' WARNING: Q > Qmax, extrapolation used; Q, Qmax =', Q, Qmax
      Endif

c Can jump to here, all of x dependence comes down to these branches
cz Reduce redundant calls to HO, HP, HP2
      HO=xv(jx)-X
      HP=xv(jx+1)-X
      HP2=xv(jx+2)-X
      branch2=.false.
      branch3=.false.
      if((x+x-xv(jx)-xv(jx+1)).gt.0d0) branch2=.true.
      if((x+x-xv(jx+1)-xv(jx+2)).gt.0d0) branch3=.true.

      do i=1,8
cz         parton(i) = CtqPdf(6-i, X, Q)
         iprtn=6-i
         If (Iprtn .GE. 3) Then
            Ip = - Iprtn
         Else
            Ip = Iprtn
         EndIf
C                             Find the off-set in the linear array Upd
         JFL = Ip + NfMx
         J0  = (JFL * (NT+1) + Jq) * (NX+1) + Jx
C
C                                           Now interpolate in x for M1 Q's
c Iqold = 1
      w=upd(j0+2)-upd(j0+1)
      DEN=HO-HP
      DEN=W/DEN
      D1=HP*DEN
      C1=HO*DEN

      w=upd(j0+3)-upd(j0+2)
      DEN=HP-HP2
      DEN=W/DEN
      D2=HP2*DEN

      W=HP*DEN-D1
      DEN=HO-HP2

      if(branch3) then
         fq(1)=upd(j0+3)+d2+hp2*w/den
      elseif(branch2) then
         fq(1)=upd(j0+2)+d1+ho*w/den
      else
         fq(1)=upd(j0+1)+c1+ho*w/den
      endif

c Iqold = 2
      w=upd(J0 + Nx+3)-upd(J0 + Nx+2)
      DEN=HO-HP
      DEN=W/DEN
      D1=HP*DEN
      C1=HO*DEN

      w=upd(J0 + Nx+4)-upd(J0 + Nx+3)
      DEN=HP-HP2
      DEN=W/DEN
      D2=HP2*DEN

      W=HP*DEN-D1
      DEN=HO-HP2

      if(branch3) then
         fq(2)=upd(J0 + Nx+4)+d2+hp2*w/den
      elseif(branch2) then
         fq(2)=upd(J0 + Nx+3)+d1+ho*w/den
      else
         fq(2)=upd(J0 + Nx+2)+c1+ho*w/den
      endif

c Iqold = 3
      j1=j0+(nx+1)*2   ! My temp ya(i)=upd(j1+i)
      w=upd(j1+2)-upd(j1+1)
      DEN=HO-HP
      DEN=W/DEN
      D1=HP*DEN
      C1=HO*DEN

      w=upd(j1+3)-upd(j1+2)
      DEN=HP-HP2
      DEN=W/DEN
      D2=HP2*DEN

      W=HP*DEN-D1
      DEN=HO-HP2

      if(branch3) then
         fq(3)=upd(j1+3)+d2+hp2*w/den
      elseif(branch2) then
         fq(3)=upd(j1+2)+d1+ho*w/den
      else
         fq(3)=upd(j1+1)+c1+ho*w/den
      endif

C                                          Finish off by interpolating in Q
      QHO=QL(JQ)-qg
      QHP=QL(JQ+1)-qg
      w=fq(2)-fq(1)
      DEN=QHO-QHP
      DEN=W/DEN
      D1=QHP*DEN
      C1=QHO*DEN

      QHP2=QL(JQ+2)-qg
      w=fq(3)-fq(2)
      DEN=QHP-QHP2
      DEN=W/DEN
      D2=QHP2*DEN

      W=QHP*DEN-D1
      DEN=QHO-QHP2

      if((qg+qg-ql(jq+1)-ql(jq+2)).gt.0d0) then
         ftmp=fq(3)+d2+qhp2*w/den
      elseif((qg+qg-ql(jq)-ql(jq+1)).gt.0d0) then
         ftmp=fq(2)+d1+qho*w/den
      else
         ftmp=fq(1)+c1+qho*w/den
      endif

         parton(i)=ftmp
c Check that makes this equivalent to CTEQ5, though questionable. 8/26/04
c         if (parton(i).lt.0d0) print *,x,parton(i)
         if (parton(i).lt.0d0) parton(i)=0d0
      enddo
C

cz End New code
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


c======================================================================
      subroutine structm6(x,q,uv,dv,us,ds,s,c,b,t,g)
c======================================================================
c     *****************************************************************
c     * This is an interface for the CTEQ6 table functions:           *
c     *  CtqPdf (Iparton, X, Q)                                       *
c     * This is called like structm:                                  *
c     *   call structm(x,q,u_v,d_v,u_s,d_s,s,c,b,t,g)                 *
c     *****************************************************************
C   Before using the PDF, it is necessary to do the initialization by
C       Call SetCtq(6,Iset)
C     Where `Iset' is from the tables.
C   *** The version ONLY works with CTEQ6 PDFs. ***
C
C   Written by Z. Sullivan, May 2004
C     Absorbs and optimizes PartonX6 from the base CTEQ6 code.
C     Requires common block / CtqParz1 / to be filled.
C
      Implicit Double Precision (A-H,O-Z)

      integer mxx, mxq, mxf, mxpqx
      Parameter (MXX = 96, MXQ = 20, MXF = 5)
      Parameter (MXQX= MXQ * MXX,   MXPQX = MXQX * (MXF+3))
c
c     Local
c
      double precision x,q,uv,dv,us,ds,s,c,b,t,g
      double precision parton(8)
      integer i
c     ----------
c     Begin Code
c     ----------
cz Begin New code

c  Given the parton distribution function in the array U in
c  COMMON / PEVLDT / , this routine interpolates to find 
c  the parton distribution at an arbitray point in x and q.
c
      double precision tt

      double precision Al, XV(0:MXX), TV(0:MXQ), UPD(MXPQX)
      integer Nx, Nt, NfMx
      double precision Qini, Qmax, Xmin
 
      Common
     > / CtqParz1 / Al, XV, TV, UPD
     > / CtqPar2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin

      double precision fvec(4), fij(4), xvpow(0:mxx)
      double precision onep, xpow
      integer nqvec, ientry
      parameter (onep = 1.00001d0, xpow = 0.3d0)
      parameter (nqvec = 4)
cz      Data OneP / 1.00001 /
cz      Data xpow / 0.3d0 /	!**** choice of interpolation variable
cz      Data nqvec / 4 /
      Data ientry / 0 /
      Save ientry,xvpow

c store the powers used for interpolation on first call...
      if(ientry .eq. 0) then
         ientry = 1

         xvpow(0) = 0D0
         do i = 1, nx
            xvpow(i) = xv(i)**xpow
         enddo
      endif

cz      X = XX
cz      Q = QQ
      tt = log(log(Q/Al))

c      -------------    find lower end of interval containing x, i.e.,
c                       get jx such that xv(jx) .le. x .le. xv(jx+1)...
      JLx = -1
      JU = Nx+1
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
C
      If     (JLx .LE. -1) Then
        Print '(A,1pE12.4)', 'Severe error: x <= 0 in structm6! x = ', x 
        Stop 
      ElseIf (JLx .Eq. 0) Then
         Jx = 0
      Elseif (JLx .LE. Nx-2) Then

C                For interrior points, keep x in the middle, as shown above
         Jx = JLx - 1
      Elseif (JLx.Eq.Nx-1 .or. x.LT.OneP) Then

C                  We tolerate a slight over-shoot of one (OneP=1.00001),
C              perhaps due to roundoff or whatever, but not more than that.
C                                      Keep at least 4 points >= Jx
         Jx = JLx - 2
      Else
        Print '(A,1pE12.4)', 'Severe error: x > 1 in structm6! x = ', x 
        Stop 
      Endif
C          ---------- Note: JLx uniquely identifies the x-bin; Jx does not.

C                       This is the variable to be interpolated in
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

c         --------------Now find lower end of interval containing Q, i.e.,
c                          get jq such that qv(jq) .le. q .le. qv(jq+1)...
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
         Jq = 0
      Elseif (JLq .LE. Nt-2) Then
C                                  keep q in the middle, as shown above
         Jq = JLq - 1
      Else
C                         JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq.
        Jq = Nt - 3

      Endif
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


      do i=1,8
cz         parton(i) = CtqPdf(6-i, X, Q)
         iprtn=6-i

c get the pdf function values at the lattice points...

      If (Iprtn .GE. 3) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
      jtmp = ((Ip + NfMx)*(NT+1)+(jq-1))*(NX+1)+jx+1

      Do it = 1, nqvec

         J1  = jtmp + it*(NX+1) 

       If (Jx .Eq. 0) Then
C                          For the first 4 x points, interpolate x^2*f(x,Q)
C                           This applies to the two lowest bins JLx = 0, 1
C            We can not put the JLx.eq.1 bin into the "interrior" section
C                           (as we do for q), since Upd(J1) is undefined.
         fij(1) = 0
         fij(2) = Upd(J1+1) * XV(1)**2
         fij(3) = Upd(J1+2) * XV(2)**2
         fij(4) = Upd(J1+3) * XV(3)**2
C 
C                 Use Polint which allows x to be anywhere w.r.t. the grid

         Call Polint4 (XVpow(0), Fij(1), 4, ss, Fx, Dfx) 
         
         If (x .GT. 0D0)  Fvec(it) =  Fx / x**2 
C                                              Pdf is undefined for x.eq.0
       ElseIf  (JLx .Eq. Nx-1) Then
C                                                This is the highest x bin:

        Call Polint4 (XVpow(Nx-3), Upd(J1), 4, ss, Fx, Dfx)

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
C                         1st Q-bin, as well as extrapolation to lower Q
        Call Polint4 (TV(0), Fvec(1), 4, tt, ff, Dfq)

      ElseIf (JLq .GE. Nt-1) Then
C                         Last Q-bin, as well as extrapolation to higher Q
        Call Polint4 (TV(Nt-3), Fvec(1), 4, tt, ff, Dfq)
      Else
C                         Interrior bins : (JLq.GE.1 .and. JLq.LE.Nt-2)
C       which include JLq.Eq.1 and JLq.Eq.Nt-2, since Upd is defined for
C                         the full range QV(0:Nt)  (in contrast to XV)
        tf2 = fvec(2)
        tf3 = fvec(3)

        g1 = ( tf2*t13 - tf3*t12) / t23
        g4 = (-tf2*t34 + tf3*t24) / t23

        h00 = ((t34*ty2-tmp2*ty3)*(fvec(1)-g1)/t12 
     &	  +  (tmp1*ty2-t12*ty3)*(fvec(4)-g4)/t34)

        ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23
      EndIf

         parton(i)=ff
c Check that makes this equivalent to CTEQ6, though questionable. 8/26/04
c         if (parton(i).lt.0d0) print *,x,parton(i)
         if (parton(i).lt.0d0) parton(i)=0d0
      enddo
C

cz End New code
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


C  This is a specialized recoding of Neville's algorithm based on the
C   POLINT routine from "Numerical Recipes", but assuming N=3, and
C   ignoring the error estimation.
C  Written by Z. Sullivan, May 2004
C  This file uses a minimal number of instructions to do 3-point fitting.
C     SUBROUTINE POLINT  (XA,YA,3,X,Y,IGNORED)
      SUBROUTINE POLINT3 (XA,YA,N,X,Y,DY)
      IMPLICIT NONE
      DOUBLE PRECISION XA(3),YA(3),X,Y,DY
      DOUBLE PRECISION C1,HO,HP,HP2,W,D1,D2,DEN
      INTEGER N

      HO=XA(1)-X
      HP=XA(2)-X
      W=YA(2)-YA(1)
      DEN=HO-HP
      DEN=W/DEN
      D1=HP*DEN
      C1=HO*DEN

      HP2=XA(3)-X
      W=YA(3)-YA(2)
      DEN=HP-HP2
      DEN=W/DEN
      D2=HP2*DEN

      W=HP*DEN-D1
      DEN=HO-HP2

      IF((X+X-XA(2)-XA(3)).GT.0D0) THEN
         Y=YA(3)+D2+HP2*W/DEN
      ELSEIF((X+X-XA(1)-XA(2)).GT.0D0) THEN
         Y=YA(2)+D1+HO*W/DEN
      ELSE
         Y=YA(1)+C1+HO*W/DEN
      ENDIF

      RETURN
      END

cz This is a specialized recoding that assumes N=4.
cz Modified from CTEQ source by Z. Sullivan, March 2004
cz      SUBROUTINE POLINT (XA,YA,4,X,Y,DY)
      SUBROUTINE POLINT4 (XA,YA,N,X,Y,DY)
C     Adapted from "Numerical Recipes" 
      IMPLICIT NONE
      DOUBLE PRECISION XA(4),YA(4),X,Y,DY
      DOUBLE PRECISION C(4),DIF,DIFT,HO,HP,W,D(4),DEN
      INTEGER NS,I,M,N

      NS=1
      DIF=DABS(X-XA(1))
      DO 11 I=1,4
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,3
        DO 12 I=1,4-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
cz          IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.4-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
