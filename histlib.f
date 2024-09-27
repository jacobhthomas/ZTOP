C.--------------------------------------------------------------------.C
C. POOR MAN'S HISTOGRAMMING PACKAGE VERSION 1.0                       .C
C  BY B.W.HARRIS 8/96                                                 .C
C.                                                                    .C
C. VERY SIMPLE HISTOGRAMMER DESIGNED TO LOOK LIKE HBOOK BUT HAVE LESS .C
C. OVERHEAD AND PERFORM ALL NUMERICAL OPERATIONS IN DOUBLE PRECISION  .C
C.                                                                    .C
C. ROUTINES:                                                          .C
C. =========                                                          .C
C. IHIST  - INTIALIZE HISTOGRAMS                                      .C
C. BBOOK  - BOOK HISTOGRAMS                                           .C
C. BFILL  - FILL HISTOGRAMS                                           .C
C. BWRITE - WRITE OUT HISTOGRAMS TO FILE                              .C
C. BINTE  - INTEGRATE UP A HISTOGRAM                                  .C
C.                                                                    .C
cz modified 4/12/01 to be 30 histograms
C. MAXIMUM NUMBER OF HISTOGRAMS = 30                                  .C
cz modified 2/23/01 to be 1000 points
C. MAXIMUM NUMBER OF BINS PER HISTOGRAM = 1000                        .C
C.                                                                    .C
C. COMMON/BHIST/                                                      .C
C.                                                                    .C
C. HMAX(N)    MAX X (RIGHT SIDE OF LAST BIN)                          .C
C. HMIN(N)    MIN X (LEFT SIDE OF FIRST BIN)                          .C
C. DX(N)      BIN WIDTH = (XMAX-XMIN)/NBINS                           .C
C. HIST(N,K)  ACCUMULATED VALUE IN KTH BIN                            .C
C. HIST2(N,K) ACCUMULATED SQUARES IN KTH BIN                          .C
C. XHIS(N,K)  CENTRAL X VALUE OF KTH BIN                              .C
C.                                                                    .C
C. NBINS(N)   NUMBER OF BINS                                          .C
C. NOVER(N)   NUMBER OF OVERFLOWS  (X>XMAX)                           .C
C. NUNDER(N)  NUMBER OF UNDERFLOWS (X<XMIN)                           .C
C. NENTRY(N)  NUMBER OF ENTRIES (XMIN<X<XMAX)                         .C
C. IHIS(N,K)  NUMBER OF ENTRIES IN KTH BIN                            .C
C.                                                                    .C
C. OUTFILE(N) OUTFILE                                                 .C
C.--------------------------------------------------------------------.C
      SUBROUTINE IHIST
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*40 OUTFILE
      PARAMETER (MH=30)
cz      PARAMETER (MH=25)
cz      PARAMETER (MB=100)
      PARAMETER (MB=1000)
      COMMON/BHIST/HMAX(MH),HMIN(MH),DX(MH),HIST(MH,MB),HIST2(MH,MB),
     &XHIS(MH,MB),NBINS(MH),NOVER(MH),NUNDER(MH),NENTRY(MH),IHIS(MH,MB),
     &OUTFILE(MH)
      DO 1 I=1,MH
         HMAX(I)=0D0
         HMIN(I)=0D0
         DX(I)=0D0
         NBINS(I)=0
         NOVER(I)=0
         NUNDER(I)=0
         NENTRY(I)=0
 1    CONTINUE
      DO 2 J=1,MB
         DO 2 I=1,MH
            HIST(I,J)=0D0
            XHIS(I,J)=0D0
            IHIS(I,J)=0
 2    CONTINUE
      RETURN
      END

      SUBROUTINE BBOOK(N,NAME,NB,XMAX,XMIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*40 OUTFILE,NAME
      PARAMETER (MH=30)
cz      PARAMETER (MH=25)
cz      PARAMETER (MB=100)
      PARAMETER (MB=1000)
      COMMON/BHIST/HMAX(MH),HMIN(MH),DX(MH),HIST(MH,MB),HIST2(MH,MB),
     &XHIS(MH,MB),NBINS(MH),NOVER(MH),NUNDER(MH),NENTRY(MH),IHIS(MH,MB),
     &OUTFILE(MH)
      IF(N.GT.MH)THEN
         WRITE(*,*) 'TOO MANY HISTOGRAMS BOOKED: ',N,MH
         STOP
      ENDIF
      OUTFILE(N)=NAME
      HMAX(N)=XMAX
      HMIN(N)=XMIN
      NBINS(N)=NB
      DX(N)=(XMAX-XMIN)/NB
      DO 1 I=1,NB
        XHIS(N,I)=XMIN+DX(N)*(I-0.5D0)
 1    CONTINUE
      RETURN
      END

      SUBROUTINE BFILL(N,X,WGT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*40 OUTFILE
      PARAMETER (MH=30)
cz      PARAMETER (MH=25)
cz      PARAMETER (MB=100)
      PARAMETER (MB=1000)
      COMMON/BHIST/HMAX(MH),HMIN(MH),DX(MH),HIST(MH,MB),HIST2(MH,MB),
     &XHIS(MH,MB),NBINS(MH),NOVER(MH),NUNDER(MH),NENTRY(MH),IHIS(MH,MB),
     &OUTFILE(MH)
      F=(X-HMIN(N))/DX(N)
      I=INT(F)+1
      IF(X.LT.HMIN(N))I=0
      IF(X.GT.HMAX(N))I=NBINS(N)+1
      IF(I.LT.1)THEN
         NUNDER(N)=NUNDER(N)+1
      ELSEIF(I.GT.NBINS(N))THEN
         NOVER(N)=NOVER(N)+1
      ELSE
         NENTRY(N)=NENTRY(N)+1
         IHIS(N,I)=IHIS(N,I)+1
         FT=WGT/DX(N)
         HIST(N,I)=HIST(N,I)+FT
         HIST2(N,I)=HIST2(N,I)+FT*FT
      ENDIF
      RETURN
      END

      SUBROUTINE BWRITE(N,FAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*40 OUTFILE
      PARAMETER (MH=30)
cz      PARAMETER (MH=25)
cz      PARAMETER (MB=100)
      PARAMETER (MB=1000)
      COMMON/BHIST/HMAX(MH),HMIN(MH),DX(MH),HIST(MH,MB),HIST2(MH,MB),
     &XHIS(MH,MB),NBINS(MH),NOVER(MH),NUNDER(MH),NENTRY(MH),IHIS(MH,MB),
     &OUTFILE(MH)
      NUNIT=13
      OPEN(UNIT=NUNIT,FILE=OUTFILE(N),STATUS='UNKNOWN')
      DO 1 I=1,NBINS(N)
         WRITE(NUNIT,*) XHIS(N,I),HIST(N,I)/FAC,DSQRT(HIST2(N,I))/FAC
 1    CONTINUE
      CLOSE(NUNIT)
      RETURN
      END

      SUBROUTINE BINTE(N,SUM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*40 OUTFILE
      PARAMETER (MH=30)
cz      PARAMETER (MH=25)
cz      PARAMETER (MB=100)
      PARAMETER (MB=1000)
      COMMON/BHIST/HMAX(MH),HMIN(MH),DX(MH),HIST(MH,MB),HIST2(MH,MB),
     &XHIS(MH,MB),NBINS(MH),NOVER(MH),NUNDER(MH),NENTRY(MH),IHIS(MH,MB),
     &OUTFILE(MH)
      SUM=0D0
      DO 1 I=1,NBINS(N)
        SUM=SUM+HIST(N,I)*DX(N)
 1    CONTINUE
      RETURN
      END

      SUBROUTINE BGET(N,I,X,H)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*40 OUTFILE
      PARAMETER (MH=30)
cz      PARAMETER (MH=25)
cz      PARAMETER (MB=100)
      PARAMETER (MB=1000)
      COMMON/BHIST/HMAX(MH),HMIN(MH),DX(MH),HIST(MH,MB),HIST2(MH,MB),
     &XHIS(MH,MB),NBINS(MH),NOVER(MH),NUNDER(MH),NENTRY(MH),IHIS(MH,MB),
     &OUTFILE(MH)
      X=XHIS(N,I)
      H=HIST(N,I)
      RETURN
      END
C.--------------------------------------------------------------------.C
C.                                                                    .C
C.                      END OF FILE                                   .C
C.                                                                    .C
C.--------------------------------------------------------------------.C
