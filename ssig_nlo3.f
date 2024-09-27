c.
c. 2 -> 3 nlo cross section, corrected 6/12/02
c.
      subroutine ssig_nlo3(p,tot)
      implicit none
C  
C ARGUMENTS 
C  
      double precision p(0:3,5),tot
C  
C CONSTANTS
C  
      double precision pi
      parameter       (pi=3.14159265358979d0)
      double precision zeta2
      parameter       (zeta2=1.64493406684823d0)
C
C LOCAL VARIABLES 
C
      DOUBLE PRECISION XPQ1(-6:6),XPQ2(-6:6)
      double precision f1(-6:6),f2(-6:6)
      double precision xl(4),ckm1,ckm2,ckm3
      double precision tmp1,tmp2,tmp3,tmp4
      double precision e3,e4,e5,soft,coll,rs12
      double precision s12,t13,t14,t15,t23,t24,t25,s34,s35,s45
      double precision t14p,t24p,s34p,s45p
      integer i13,i14,i15,i23,i24,i25,i34,i35,i45,is3,is4,is5
      integer i,itest
C
C EXTERNAL FUNCTIONS
C
      DOUBLE PRECISION dot,psi1,psi3,xlog,xlog2,dilog
C
C GLOBAL VARIABLES
C
      double precision deltas,deltac
      common/slice/deltas,deltac
      double precision hbarc2,gweak2,xm2,wm2
      common/parm/hbarc2,gweak2,xm2,wm2
      double precision xa,xb,scale
      common/toevent/  xa,xb,scale
      integer              icollide,ismear,igroup,iset,iq2
      common /to_collider2/icollide,ismear,igroup,iset,iq2
      double precision vud2,vus2,vub2,vcd2,vcs2,vcb2,vtd2,vts2,vtb2
      common /ckm/     vud2,vus2,vub2,vcd2,vcs2,vcb2,vtd2,vts2,vtb2
      integer      itop,iatop
      common/ttbar/itop,iatop
c
c Begin code
c
      s12=2d0*dot(p(0,1),p(0,2))
      t13=-2d0*dot(p(0,1),p(0,3))
      t14p=-2d0*dot(p(0,1),p(0,4))
      t23=-2d0*dot(p(0,2),p(0,3))
      t24p=-2d0*dot(p(0,2),p(0,4))
      t24=t24p+xm2

      s35=s12+t14p+t24
      s45=s12+t13+t23
      s45p=s45-xm2
      s34=s12-s35-s45p
      t15=t23+t24p+s34
      t25=t13+t14p+s34
c
      t14=t14p+xm2
      s34p=s34-xm2
c
      rs12=dsqrt(s12)
      soft=0.5d0*deltas*rs12
      coll=deltac*s12

      i13=0
      i14=0
      i15=0
      i23=0
      i24=0
      i25=0
      i34=0
      i35=0
      i45=0

      is3=0
      is4=0
      is5=0

c collinear
      if (dabs(t13).le.coll) i13=1
      if (dabs(t14).le.coll) i14=1
      if (dabs(t15).le.coll) i15=1
      if (dabs(t23).le.coll) i23=1
      if (dabs(t24).le.coll) i24=1
      if (dabs(t25).le.coll) i25=1
      if (dabs(s34).le.coll) i34=1
      if (dabs(s35).le.coll) i35=1
      if (dabs(s45).le.coll) i45=1
c soft
      e3=0.5d0*(s12-s45)/rs12
      e4=0.5d0*(s12-s35)/rs12
      e5=0.5d0*(s12-s34)/rs12
      if (e3.le.soft) is3=1
      if (e4.le.soft) is4=1
      if (e5.le.soft) is5=1

      tot=0d0

      CALL PFTOPDG(xa,SCALE,XPQ1)
      CALL PFTOPDG(xb,SCALE,XPQ2)
      DO I=-6,6
         XPQ1(i)=XPQ1(i)/xa
         XPQ2(i)=XPQ2(i)/xb
      ENDDO
      if (icollide .le. 2) then ! f1 is proton, f2 is antiproton
         do i=-6,6
            f1(i)=xpq1(i)
            f2(i)=xpq2(-i)
         enddo
      elseif (icollide .eq. 3) then !PP
         do i=-6,6
            f1(i)=xpq1(i)
            f2(i)=xpq2(i)
         enddo
      endif
c
ct      write (*,*) 's12,t13,t14,t15,t23,t24,t25,s34,s35,s45:'
ct      write (*,*) s12,t13,t14,t15,t23,t24,t25,s34,s35,s45
ct      write (*,*) 't14p,t24p,s34p,s45p:'
ct      write (*,*) t14p,t24p,s34p,s45p

      itest=is5+i15+i25+i35
      if (itest.eq.0) then
         call setlum(f1,f2,xl)
c    u d~ -> b~ t g
         tmp1=psi1(t13,s12,t14,t15,t23,s34,s35,t24,t25,s45,xm2,wm2)
c    d~ u -> b~ t g
         tmp2=psi1(t23,s12,t24,t25,t13,s34,s35,t14,t15,s45,xm2,wm2)
         if(itop.eq.1) tot=tot+xl(1)*tmp1+xl(2)*tmp2
c    u~ d -> b t~ g
c         tot=tot+xl(3)*tmp1
c    d u~ -> b t~ g
c         tot=tot+xl(4)*tmp2
         if(iatop.eq.1) tot=tot+xl(3)*tmp1+xl(4)*tmp2
      endif

c These were corrected 6/12/02.  Exchanged 3 <-> 5.
      itest=i15+i25
      if (itest.eq.0) then
c    v_ud^2 + v_us^2 + v_ub^2 == 1 therefore no CKM, but put in anyway
         ckm1 = vud2+vus2+vub2
         ckm2 = vcd2+vcs2+vcb2
c    g (u,c) -> b~ t (d+s+b)
         tmp1=psi3(t13,t15,t14,s12,s35,s34,t23,s45,t25,t24,xm2,wm2)
c    u g -> b~ t d
         tmp2=psi3(t23,t25,t24,s12,s35,s34,t13,s45,t15,t14,xm2,wm2)
         if(itop.eq.1) tot=tot+tmp1*f1(0)*(f2(2)*ckm1+f2(4)*ckm2)
     &        +tmp2*f2(0)*(f1(2)*ckm1+f1(4)*ckm2)
c    g u~ -> b t~ d~
c         tot=tot+tmp1*f1(0)*(f2(-2)*ckm1+f2(-4)*ckm2)
c    u~ g -> b t~ d~
c         tot=tot+tmp2*f2(0)*(f1(-2)*ckm1+f1(-4)*ckm2)
         if(iatop.eq.1) tot=tot+tmp1*f1(0)*(f2(-2)*ckm1+f2(-4)*ckm2)
     &        +tmp2*f2(0)*(f1(-2)*ckm1+f1(-4)*ckm2)
      endif

c These were corrected 6/12/02.  Exchanged 3 <-> 5.
      itest=i15+i25
      if (itest.eq.0) then
         ckm1 = vud2+vcd2
         ckm2 = vus2+vcs2
         ckm3 = vub2+vcb2
c    g (d~,s~,b~) -> b~ t (u~+c~)
         tmp1=psi3(t13,s12,t14,t15,t23,s34,s35,t24,t25,s45,xm2,wm2)
c    d~ g -> b~ t u~
         tmp2=psi3(t23,s12,t24,t25,t13,s34,s35,t14,t15,s45,xm2,wm2)
         if(itop.eq.1)tot=tot+tmp1*f1(0)*(f2(-1)*ckm1+f2(-3)*ckm2
     &        +f2(-5)*ckm3) + tmp2*f2(0)*(f1(-1)*ckm1+f1(-3)*ckm2
     &        +f1(-5)*ckm3)
c    g d -> b t~ u
c         tot=tot+tmp1*f1(0)*(f2(1)*ckm1+f2(3)*ckm2+f2(5)*ckm3)
c    d g -> b t~ u
c         tot=tot+tmp2*f2(0)*(f1(1)*ckm1+f1(3)*ckm2+f1(5)*ckm3)
         if(iatop.eq.1)tot=tot+tmp1*f1(0)*(f2(1)*ckm1+f2(3)*ckm2
     &        +f2(5)*ckm3) + tmp2*f2(0)*(f1(1)*ckm1+f1(3)*ckm2
     &        +f1(5)*ckm3)
      endif
c
      return
      end

      double precision function psi1(s12,t13,t14,t15,t23,t24,t25,
     &     s34,s35,s45,xm2,wm2)
c
c. matrix element
c. u b --> d t g
c
      double precision s12,t13,t14,t15,t23,t24,t25,s34,s35,s45
      double precision xm2,wm2
      double precision cf
      parameter       (cf=4d0/3d0)
      double precision t14p,t24p,s34p,s45p
      double precision fq1,fq2

      t14p = t14 - xm2
      t24p = t24 - xm2
      s34p = s34 - xm2
      s45p = s45 - xm2
      fq1 = 1d0 / ( t24 - wm2 ) / ( t24 - wm2 )
      fq2 = 1d0 / ( t13 - wm2 ) / ( t13 - wm2 )
      psi1 = fq1 * ( - 2d0*s12*t13/t15*s34p/s35 - s12*t13/t15*s45p/s35 
     . + s12*t14p/t15 - s12/t15*s34p - s12*s34p/s35 - s12*s45p/s35 
     . - t13/t15*s34p/s35*t25 - s34p*t25/t15 + t23*s34p/s35 )
     . + fq2 * ( 2d0*xm2*s12*s34p/s45p/s45p + 2d0*xm2*s12/s45p/s45p*s35 
     . + s12*t23/t25 - 2d0*s12*t24p*s34p/s45p/t25 
     . - s12*t24p/s45p*s35/t25 - s12*s34p/s45p - s12*s34p/t25 
     . - s12/s45p*s35 + t14p*s34p/s45p - t15*t24p*s34p/s45p/t25 
     . - t15*s34p/t25 )
      psi1=psi1*2d0*cf
      return
      end

      double precision function psi2(s12,t13,t14,t15,t23,t24,t25,
     &     s34,s35,s45,xm2,wm2)
c
c. matrix element
c. u g --> d t b~
c
      double precision s12,t13,t14,t15,t23,t24,t25,s34,s35,s45
      double precision xm2,wm2
      double precision cf
      parameter       (cf=4d0/3d0)
      double precision t14p,t24p,s34p,s45p
      double precision fq1,fq2

      t14p = t14 - xm2
      t24p = t24 - xm2
      s34p = s34 - xm2
      s45p = s45 - xm2
      psi2 = -2d0*xm2*t15/t24p/t24p * ( t23 + s34p ) 
     .     + s12/t24p/t25*s34p*s45p + s12/t25*s34p 
     .     + t15/t24p/t25*s45p* ( t23 + 2d0*s34p )
     .     + t15/t24p * ( t23 + s34p ) + t15/t25 * ( s34p - s35 )
     .     - t14p*s34p/t24p
      psi2=psi2/(t13-wm2)/(t13-wm2)
      return
      end


      double precision function psi3(s12,t13,t14,t15,t23,t24,t25,
     &     s34,s35,s45,xm2,wm2)
c
c. matrix element
c. g b --> d t u~
c
      double precision s12,t13,t14,t15,t23,t24,t25,s34,s35,s45
      double precision xm2,wm2
      double precision cf
      parameter       (cf=4d0/3d0)
      double precision t14p,t24p,s34p,s45p
      double precision fq1,fq2

      t14p = t14 - xm2
      s34p = s34 - xm2
      s45p = s45 - xm2
      psi3 = s12/t15 * ( s34p*s35/t13 + s34p ) 
     .     + t25/t15 * ( t14p*s35/t13 
     .     + 2d0*s34p*s35/t13 + s34p - s45p ) 
     .     + t25/t13 * ( t14p + s34p ) 
     .     - t23*s34p/t13
      psi3=psi3/(t24-wm2)/(t24-wm2)
      return
      end
