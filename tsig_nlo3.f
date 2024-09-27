c.
c. 2 -> 3 nlo cross section
c.
c Modified 9/17/02 to allow arbitrary scaling via isclf/isclfb, and renorm/fac
c  scales treated separately.
c
      subroutine tsig_nlo3(p,tot)
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
      DOUBLE PRECISION XPQ1(-6:6),XPQ2(-6:6),XPQ3(-6:6),XPQ4(-6:6)
      DOUBLE PRECISION XPQB1(-6:6),XPQB2(-6:6),XPQB3(-6:6),XPQB4(-6:6)
      double precision f1(-6:6),f2(-6:6),f3(-6:6),f4(-6:6)
      double precision fb1(-6:6),fb2(-6:6),fb3(-6:6),fb4(-6:6)
      double precision xl1(8),xl2(8),ckm1,ckm2,ckm3
      double precision tmp1,tmp2,tmp3,tmp4,btmp
      double precision e3,e4,e5,soft,coll,rs12
      double precision s12,t13,t14,t15,t23,t24,t25,s34,s35,s45
      double precision t14p,t24p,s34p,s45p
      integer i13,i14,i15,i23,i24,i25,i34,i35,i45,is3,is4,is5
      integer i,itest

      double precision scale1,scale2,scale3,scale4
      double precision scaleb1,scaleb2,scaleb3,scaleb4
      double precision rscale1,rscale2,rscale3,rscale4
      double precision rscaleb1,rscaleb2,rscaleb3,rscaleb4
      double precision asl1,asl2,asl3,asl4
      double precision ash1,ash2,ash3,ash4
      double precision ptmp(0:3)
C
C EXTERNAL FUNCTIONS
C
      DOUBLE PRECISION tpsi1,tpsi2,tpsi3,xlog,xlog2,dilog
      DOUBLE PRECISION dot,pt,alphas2
C
C GLOBAL VARIABLES
C
      double precision deltas,deltac
      common/slice/deltas,deltac
      double precision hbarc2,gweak2,xm2,wm2
      common/parm/hbarc2,gweak2,xm2,wm2
      integer      isclf,ialphas,iscale
      common/flags/isclf,ialphas,iscale
      integer       isclfb,iscaleb
      common/flagsb/isclfb,iscaleb
      integer       irsclf,irscale,irsclfb,irscaleb
      common/rflags/irsclf,irscale,irsclfb,irscaleb
      double precision xa,xb,scale
      common/toevent/  xa,xb,scale
      integer              icollide,ismear,igroup,iset,iq2
      common /to_collider2/icollide,ismear,igroup,iset,iq2
      double precision vud2,vus2,vub2,vcd2,vcs2,vcb2,vtd2,vts2,vtb2
      common /ckm/     vud2,vus2,vub2,vcd2,vcs2,vcb2,vtd2,vts2,vtb2
      integer      itop,iatop
      common/ttbar/itop,iatop
      double precision bwgt
      common/btagging/ bwgt

      double precision alph1,alph2
      common /sctopsi/ alph1,alph2
c
c Begin code
c
      bwgt=0d0  ! fraction of this part that has b~ in it.
      btmp=0d0

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
c
c Set scale, pdf's
c
      if(iscale.gt.7) then
      if(iscale.eq.8) then
         scale1=dsqrt(2d0*dot(p(0,5),p(0,3)))  ! M_jj^2 (light cor)
         scale2=pt(p(0,3))                     ! for 2 -> 2, M_jj=0 (heavy cor)
         scale3=scale1
         scale4=scale2
      elseif (iscale.eq.9) then                    ! Q^2+M_jj^2
         scale1=dsqrt(-t24+2d0*dot(p(0,5),p(0,3))) ! -t = Q^2 out of pbar (l)
         scale2=dsqrt(-t13)                        ! -t = Q^2 out of p    (h)
         scale3=dsqrt(-t14+2d0*dot(p(0,5),p(0,3))) ! -u = Q^2 out of p    (l)
         scale4=dsqrt(-t23)                        ! -u = Q^2 out of pbar (h)
      elseif (iscale.eq.10) then
         do i=0,3
            ptmp(i)=p(i,3)+p(i,5)
         enddo
         scale1=pt(ptmp)                 ! MT_jj^2                (light cor)
         scale2=pt(p(0,3))               ! for 2 -> 2, MT_jj=pt_d (heavy cor)
         scale3=scale1
         scale4=scale2
      elseif (iscale.eq.11) then
         do i=0,3
            ptmp(i)=p(i,3)+p(i,5)
         enddo                              ! Q^2 + MT_jj^2
         scale1=dsqrt(-t24+(pt(ptmp))**2)   ! -t = Q^2 out of pbar (l)
         scale2=dsqrt(-t13+(pt(p(0,3)))**2) ! -t = Q^2 out of p (h), MT_jj=pt_d
         scale3=dsqrt(-t14+(pt(ptmp))**2)   ! -u = Q^2 out of p    (l)
         scale4=dsqrt(-t23+(pt(p(0,3)))**2) ! -u = Q^2 out of pbar (h)
      elseif (iscale.eq.12) then
         scale1=pt(p(0,3))+pt(p(0,4))+pt(p(0,5))                 ! HT = sum PT
         scale2=scale1
         scale3=scale1
         scale4=scale1
      elseif (iscale.eq.13) then
         scale1=pt(p(0,3))+dsqrt((pt(p(0,4)))**2+xm2)+pt(p(0,5)) ! HT~ = sum MT
         scale2=scale1
         scale3=scale1
         scale4=scale1
      endif
      else            ! normal scale choices
      if(iscale.eq.4) then
         scale1=dsqrt(xm2)               ! mt
         scale2=scale1
         scale3=scale1
         scale4=scale1
      elseif (iscale.eq.5) then
         scale1=rs12                     ! shat
         scale2=scale1
         scale3=scale1
         scale4=scale1
      elseif (iscale.eq.6) then
         scale1=dsqrt(-t24)              ! -t = Q^2 out of pbar (light cor)
         scale2=dsqrt(-t13)              ! -t = Q^2 out of p    (heavy cor)
         scale3=dsqrt(-t14)              ! -u = Q^2 out of p    (light cor)
         scale4=dsqrt(-t23)              ! -u = Q^2 out of pbar (heavy cor)
      elseif (iscale.eq.7) then
         scale1=dsqrt(-t24+xm2)          ! -t = Q^2 out of pbar (light cor)
         scale2=dsqrt(-t13+xm2)          ! -t = Q^2 out of p    (heavy cor)
         scale3=dsqrt(-t14+xm2)          ! -u = Q^2 out of p    (light cor)
         scale4=dsqrt(-t23+xm2)          ! -u = Q^2 out of pbar (heavy cor)
      elseif (iscale.eq.0) then
         scale1=pt(p(0,4))
         scale2=scale1
         scale3=scale1
         scale4=scale1
      elseif (iscale.eq.1) then
         scale1=dsqrt((pt(p(0,4)))**2+xm2)
         scale2=scale1
         scale3=scale1
         scale4=scale1
      elseif (iscale.eq.2) then
         if(pt(p(0,3)).ge.pt(p(0,5))) then
            scale1=pt(p(0,3))            ! Highest PT jet is 3
         else
            scale1=pt(p(0,5))            ! Highest PT jet is 5
         endif
         scale2=scale1
         scale3=scale1
         scale4=scale1
      elseif (iscale.eq.3) then
         if(pt(p(0,3)).ge.pt(p(0,5))) then
            scale1=dsqrt(xm2+2d0*dot(p(0,4),p(0,3)))  ! 3 is highest PT jet
         else
            scale1=dsqrt(xm2+2d0*dot(p(0,4),p(0,5)))  ! 5 is highest PT jet
         endif
         scale2=scale1
         scale3=scale1
         scale4=scale1
      endif
      endif    ! iscale > 7
      if(isclf.ne.0) then
         scale1=scale1*(dble(isclf)/1d2)
         scale2=scale2*(dble(isclf)/1d2)
         scale3=scale3*(dble(isclf)/1d2)
         scale4=scale4*(dble(isclf)/1d2)
      endif

c
      if(iscaleb.gt.7) then
      if(iscaleb.eq.8) then
         scaleb1=dsqrt(xm2)                        ! for 2 -> 2, M_tj=mt (l)
         scaleb2=dsqrt(xm2+2d0*dot(p(0,4),p(0,5))) ! M_tj^2              (h)
         scaleb3=scaleb1
         scaleb4=scaleb2
      elseif (iscaleb.eq.9) then                        ! Q^2+M_tj^2
         scaleb1=dsqrt(-t24+xm2)                        !-t=Q^2 out of pbar(l)
         scaleb2=dsqrt(-t13+2d0*dot(p(0,4),p(0,5))+xm2) ! -t= Q^2 out of p (h)
         scaleb3=dsqrt(-t14+xm2)                        ! -u= Q^2 out of p (l)
         scaleb4=dsqrt(-t23+2d0*dot(p(0,4),p(0,5))+xm2) !-u=Q^2 out of pbar(h)
      elseif (iscaleb.eq.10) then
         do i=0,3
            ptmp(i)=p(i,4)+p(i,5)
         enddo
         scaleb1=dsqrt(xm2+pt(p(0,4))**2)          ! MT_tj^2  (light cor)
         scaleb2=dsqrt(xm2+pt(ptmp)**2)            !          (heavy cor)
         scaleb3=scaleb1
         scaleb4=scaleb2
      elseif (iscaleb.eq.11) then
         do i=0,3
            ptmp(i)=p(i,4)+p(i,5)
         enddo                                   ! Q^2 + MT_tj^2
         scaleb1=dsqrt(-t24+xm2+(pt(p(0,4)))**2) ! -t = Q^2 out of pbar (l)
         scaleb2=dsqrt(-t13+xm2+(pt(ptmp))**2)   ! -t = Q^2 out of p    (h)
         scaleb3=dsqrt(-t14+xm2+(pt(p(0,4)))**2) ! -u = Q^2 out of p    (l)
         scaleb4=dsqrt(-t23+xm2+(pt(ptmp))**2)   ! -u = Q^2 out of pbar (h)
      elseif (iscaleb.eq.12) then
         scaleb1=pt(p(0,3))+pt(p(0,4))+pt(p(0,5))                 ! HT= sum PT
         scaleb2=scaleb1
         scaleb3=scaleb1
         scaleb4=scaleb1
      elseif (iscaleb.eq.13) then
         scaleb1=pt(p(0,3))+dsqrt((pt(p(0,4)))**2+xm2)+pt(p(0,5)) ! HT~= sum MT
         scaleb2=scaleb1
         scaleb3=scaleb1
         scaleb4=scaleb1
      endif
      else            ! normal scale choices
      if(iscaleb.eq.4) then
         scaleb1=dsqrt(xm2)               ! mt
         scaleb2=scaleb1
         scaleb3=scaleb1
         scaleb4=scaleb1
      elseif (iscaleb.eq.5) then
         scaleb1=rs12                     ! shat
         scaleb2=scaleb1
         scaleb3=scaleb1
         scaleb4=scaleb1
      elseif (iscaleb.eq.6) then
         scaleb1=dsqrt(-t24)              ! -t = Q^2 out of pbar (light cor)
         scaleb2=dsqrt(-t13)              ! -t = Q^2 out of p    (heavy cor)
         scaleb3=dsqrt(-t14)              ! -u = Q^2 out of p    (light cor)
         scaleb4=dsqrt(-t23)              ! -u = Q^2 out of pbar (heavy cor)
      elseif (iscaleb.eq.7) then
         scaleb1=dsqrt(-t24+xm2)          ! -t = Q^2 out of pbar (light cor)
         scaleb2=dsqrt(-t13+xm2)          ! -t = Q^2 out of p    (heavy cor)
         scaleb3=dsqrt(-t14+xm2)          ! -u = Q^2 out of p    (light cor)
         scaleb4=dsqrt(-t23+xm2)          ! -u = Q^2 out of pbar (heavy cor)
      elseif (iscaleb.eq.0) then
         scaleb1=pt(p(0,4))
         scaleb2=scaleb1
         scaleb3=scaleb1
         scaleb4=scaleb1
      elseif (iscaleb.eq.1) then
         scaleb1=dsqrt((pt(p(0,4)))**2+xm2)
         scaleb2=scaleb1
         scaleb3=scaleb1
         scaleb4=scaleb1
      elseif (iscaleb.eq.2) then
         if(pt(p(0,3)).ge.pt(p(0,5))) then
            scaleb1=pt(p(0,3))            ! Highest PT jet is 3
         else
            scaleb1=pt(p(0,5))            ! Highest PT jet is 5
         endif
         scaleb2=scaleb1
         scaleb3=scaleb1
         scaleb4=scaleb1
      elseif (iscaleb.eq.3) then
         if(pt(p(0,3)).ge.pt(p(0,5))) then
            scaleb1=dsqrt(xm2+2d0*dot(p(0,4),p(0,3)))  ! 3 is highest PT jet
         else
            scaleb1=dsqrt(xm2+2d0*dot(p(0,4),p(0,5)))  ! 5 is highest PT jet
         endif
         scaleb2=scaleb1
         scaleb3=scaleb1
         scaleb4=scaleb1
      endif
      endif
      if(isclfb.ne.0) then
         scaleb1=scaleb1*(dble(isclfb)/1d2)
         scaleb2=scaleb2*(dble(isclfb)/1d2)
         scaleb3=scaleb3*(dble(isclfb)/1d2)
         scaleb4=scaleb4*(dble(isclfb)/1d2)
      endif

cz Scales for renormalization
      if(irscale.gt.7) then
      if(irscale.eq.8) then
         rscale1=dsqrt(2d0*dot(p(0,5),p(0,3)))  ! M_jj^2 (light cor)
         rscale2=pt(p(0,3))                     ! for 2 -> 2,M_jj=0 (heavy cor)
         rscale3=rscale1
         rscale4=rscale2
      elseif (irscale.eq.9) then                    ! Q^2+M_jj^2
         rscale1=dsqrt(-t24+2d0*dot(p(0,5),p(0,3))) ! -t = Q^2 out of pbar (l)
         rscale2=dsqrt(-t13)                        ! -t = Q^2 out of p    (h)
         rscale3=dsqrt(-t14+2d0*dot(p(0,5),p(0,3))) ! -u = Q^2 out of p    (l)
         rscale4=dsqrt(-t23)                        ! -u = Q^2 out of pbar (h)
      elseif (irscale.eq.10) then
         do i=0,3
            ptmp(i)=p(i,3)+p(i,5)
         enddo
         rscale1=pt(ptmp)                 ! MT_jj^2                (light cor)
         rscale2=pt(p(0,3))               ! for 2 -> 2, MT_jj=pt_d (heavy cor)
         rscale3=rscale1
         rscale4=rscale2
      elseif (irscale.eq.11) then
         do i=0,3
            ptmp(i)=p(i,3)+p(i,5)
         enddo                            ! Q^2 + MT_jj^2
         rscale1=dsqrt(-t24+pt(ptmp)**2)   ! -t = Q^2 out of pbar (l)
         rscale2=dsqrt(-t13+pt(p(0,3))**2) ! -t = Q^2 out of p (h),MT_jj=pt_d
         rscale3=dsqrt(-t14+pt(ptmp)**2)   ! -u = Q^2 out of p    (l)
         rscale4=dsqrt(-t23+pt(p(0,3))**2) ! -u = Q^2 out of pbar (h)
      elseif (irscale.eq.12) then
         rscale1=pt(p(0,3))+pt(p(0,4))+pt(p(0,5))              ! HT = sum PT
         rscale2=rscale1
         rscale3=rscale1
         rscale4=rscale1
      elseif (irscale.eq.13) then                              ! HT~ = sum MT
         rscale1=pt(p(0,3))+dsqrt((pt(p(0,4)))**2+xm2)+pt(p(0,5))
         rscale2=rscale1
         rscale3=rscale1
         rscale4=rscale1
      endif
      else            ! normal rscale choices
      if(irscale.eq.4) then
         rscale1=dsqrt(xm2)               ! mt
         rscale2=rscale1
         rscale3=rscale1
         rscale4=rscale1
      elseif (irscale.eq.5) then
         rscale1=rs12                     ! shat
         rscale2=rscale1
         rscale3=rscale1
         rscale4=rscale1
      elseif (irscale.eq.6) then
         rscale1=dsqrt(-t24)              ! -t = Q^2 out of pbar (light cor)
         rscale2=dsqrt(-t13)              ! -t = Q^2 out of p    (heavy cor)
         rscale3=dsqrt(-t14)              ! -u = Q^2 out of p    (light cor)
         rscale4=dsqrt(-t23)              ! -u = Q^2 out of pbar (heavy cor)
      elseif (irscale.eq.7) then
         rscale1=dsqrt(-t24+xm2)          ! -t = Q^2 out of pbar (light cor)
         rscale2=dsqrt(-t13+xm2)          ! -t = Q^2 out of p    (heavy cor)
         rscale3=dsqrt(-t14+xm2)          ! -u = Q^2 out of p    (light cor)
         rscale4=dsqrt(-t23+xm2)          ! -u = Q^2 out of pbar (heavy cor)
      elseif (irscale.eq.0) then
         rscale1=pt(p(0,4))
         rscale2=rscale1
         rscale3=rscale1
         rscale4=rscale1
      elseif (irscale.eq.1) then
         rscale1=dsqrt((pt(p(0,4)))**2+xm2)
         rscale2=rscale1
         rscale3=rscale1
         rscale4=rscale1
      elseif (irscale.eq.2) then
         if(pt(p(0,3)).ge.pt(p(0,5))) then
            rscale1=pt(p(0,3))            ! Highest PT jet is 3
         else
            rscale1=pt(p(0,5))            ! Highest PT jet is 5
         endif
         rscale2=rscale1
         rscale3=rscale1
         rscale4=rscale1
      elseif (irscale.eq.3) then
         if(pt(p(0,3)).ge.pt(p(0,5))) then
            rscale1=dsqrt(xm2+2d0*dot(p(0,4),p(0,3)))  ! 3 is highest PT jet
         else
            rscale1=dsqrt(xm2+2d0*dot(p(0,4),p(0,5)))  ! 5 is highest PT jet
         endif
         rscale2=rscale1
         rscale3=rscale1
         rscale4=rscale1
      endif
      endif    ! irscale > 7
      if(irsclf.ne.0) then
         rscale1=rscale1*(dble(irsclf)/1d2)
         rscale2=rscale2*(dble(irsclf)/1d2)
         rscale3=rscale3*(dble(irsclf)/1d2)
         rscale4=rscale4*(dble(irsclf)/1d2)
      endif
c
      if(irscaleb.gt.7) then
      if(irscaleb.eq.8) then
         rscaleb1=dsqrt(xm2)                        ! for 2 -> 2, M_tj=mt (l)
         rscaleb2=dsqrt(xm2+2d0*dot(p(0,4),p(0,5))) ! M_tj^2              (h)
         rscaleb3=rscaleb1
         rscaleb4=rscaleb2
      elseif (irscaleb.eq.9) then                        ! Q^2+M_tj^2
         rscaleb1=dsqrt(-t24+xm2)                        !-t=Q^2 out of pbar(l)
         rscaleb2=dsqrt(-t13+2d0*dot(p(0,4),p(0,5))+xm2) ! -t= Q^2 out of p (h)
         rscaleb3=dsqrt(-t14+xm2)                        ! -u= Q^2 out of p (l)
         rscaleb4=dsqrt(-t23+2d0*dot(p(0,4),p(0,5))+xm2) !-u=Q^2 out of pbar(h)
      elseif (irscaleb.eq.10) then
         do i=0,3
            ptmp(i)=p(i,4)+p(i,5)
         enddo
         rscaleb1=dsqrt(xm2+pt(p(0,4))**2)          ! MT_tj^2  (light cor)
         rscaleb2=dsqrt(xm2+pt(ptmp)**2)            !          (heavy cor)
         rscaleb3=rscaleb1
         rscaleb4=rscaleb2
      elseif (irscaleb.eq.11) then
         do i=0,3
            ptmp(i)=p(i,4)+p(i,5)
         enddo                                 ! Q^2 + MT_tj^2
         rscaleb1=dsqrt(-t24+xm2+pt(p(0,4))**2) ! -t = Q^2 out of pbar (l)
         rscaleb2=dsqrt(-t13+xm2+pt(ptmp)**2)   ! -t = Q^2 out of p    (h)
         rscaleb3=dsqrt(-t14+xm2+pt(p(0,4))**2) ! -u = Q^2 out of p    (l)
         rscaleb4=dsqrt(-t23+xm2+pt(ptmp)**2)   ! -u = Q^2 out of pbar (h)
      elseif (irscaleb.eq.12) then
         rscaleb1=pt(p(0,3))+pt(p(0,4))+pt(p(0,5))              ! HT= sum PT
         rscaleb2=rscaleb1
         rscaleb3=rscaleb1
         rscaleb4=rscaleb1
      elseif (irscaleb.eq.13) then                              ! HT~= sum MT
         rscaleb1=pt(p(0,3))+dsqrt((pt(p(0,4)))**2+xm2)+pt(p(0,5))
         rscaleb2=rscaleb1
         rscaleb3=rscaleb1
         rscaleb4=rscaleb1
      endif
      else            ! normal rscale choices
      if(irscaleb.eq.4) then
         rscaleb1=dsqrt(xm2)               ! mt
         rscaleb2=rscaleb1
         rscaleb3=rscaleb1
         rscaleb4=rscaleb1
      elseif (irscaleb.eq.5) then
         rscaleb1=rs12                     ! shat
         rscaleb2=rscaleb1
         rscaleb3=rscaleb1
         rscaleb4=rscaleb1
      elseif (irscaleb.eq.6) then
         rscaleb1=dsqrt(-t24)              ! -t = Q^2 out of pbar (light cor)
         rscaleb2=dsqrt(-t13)              ! -t = Q^2 out of p    (heavy cor)
         rscaleb3=dsqrt(-t14)              ! -u = Q^2 out of p    (light cor)
         rscaleb4=dsqrt(-t23)              ! -u = Q^2 out of pbar (heavy cor)
      elseif (irscaleb.eq.7) then
         rscaleb1=dsqrt(-t24+xm2)          ! -t = Q^2 out of pbar (light cor)
         rscaleb2=dsqrt(-t13+xm2)          ! -t = Q^2 out of p    (heavy cor)
         rscaleb3=dsqrt(-t14+xm2)          ! -u = Q^2 out of p    (light cor)
         rscaleb4=dsqrt(-t23+xm2)          ! -u = Q^2 out of pbar (heavy cor)
      elseif (irscaleb.eq.0) then
         rscaleb1=pt(p(0,4))
         rscaleb2=rscaleb1
         rscaleb3=rscaleb1
         rscaleb4=rscaleb1
      elseif (irscaleb.eq.1) then
         rscaleb1=dsqrt((pt(p(0,4)))**2+xm2)
         rscaleb2=rscaleb1
         rscaleb3=rscaleb1
         rscaleb4=rscaleb1
      elseif (irscaleb.eq.2) then
         if(pt(p(0,3)).ge.pt(p(0,5))) then
            rscaleb1=pt(p(0,3))            ! Highest PT jet is 3
         else
            rscaleb1=pt(p(0,5))            ! Highest PT jet is 5
         endif
         rscaleb2=rscaleb1
         rscaleb3=rscaleb1
         rscaleb4=rscaleb1
      elseif (irscaleb.eq.3) then
         if(pt(p(0,3)).ge.pt(p(0,5))) then
            rscaleb1=dsqrt(xm2+2d0*dot(p(0,4),p(0,3)))  ! 3 is highest PT jet
         else
            rscaleb1=dsqrt(xm2+2d0*dot(p(0,4),p(0,5)))  ! 5 is highest PT jet
         endif
         rscaleb2=rscaleb1
         rscaleb3=rscaleb1
         rscaleb4=rscaleb1
      endif
      endif
      if(irsclfb.ne.0) then
         rscaleb1=rscaleb1*(dble(irsclfb)/1d2)
         rscaleb2=rscaleb2*(dble(irsclfb)/1d2)
         rscaleb3=rscaleb3*(dble(irsclfb)/1d2)
         rscaleb4=rscaleb4*(dble(irsclfb)/1d2)
      endif

c
c test
c
c      print *,'scales:',scale1,scale2,scale3,scale4
c      print *,'scaleb:',scaleb1,scaleb2,scaleb3,scaleb4
c      print *,'rscale:',rscale1,rscale2,rscale3,rscale4
c      print *,'rsclb :',rscaleb1,rscaleb2,rscaleb3,rscaleb4
c      print *,'s12,t13,t14,t15,t23,t24,t25,s34,s35,s45:',s12,t13,t14,
c     &     t15,t23,t24,t25,s34,s35,s45
c      print *,'t14p,t24p,s34p,s45p:',t14p,t24p,s34p,s45p
c      print *,'i13,i14,i15,i23,i24,i25,i34,i35,i45,is3,is4,is5:',
c     &     i13,i14,i15,i23,i24,i25,i34,i35,i45,is3,is4,is5

c
c Get PDFs, Alpha_s
c
cz Need both factorization and renormalization scales > 1 GeV to work.
      if ((scale1.gt.1d0).and.(rscale1.gt.1d0)) then
         CALL PFTOPDG(xa,SCALE1,XPQ1)
         asl1=alphas2(rscale1)
      else
         call fillxpq0(XPQ1)
         asl1=0d0           ! 10 would fine as well
      endif
      if (scale2.gt.1d0) then
         CALL PFTOPDG(xa,SCALE2,XPQ2)
c        asl2=alphas2(rscale2)  ! not used, so comment to speed up
      else
         call fillxpq0(XPQ2)
         asl2=0d0           ! 10 would fine as well
      endif
      if ((scale3.gt.1d0).and.(rscale3.gt.1d0)) then
         CALL PFTOPDG(xb,SCALE3,XPQ3)
         asl3=alphas2(rscale3)
      else
         call fillxpq0(XPQ3)
         asl3=0d0           ! 10 would fine as well
      endif
      if (scale4.gt.1d0) then
         CALL PFTOPDG(xb,SCALE4,XPQ4)
c         asl4=alphas2(rscale4)
      else
         call fillxpq0(XPQ4)
         asl4=0d0           ! 10 would fine as well
      endif
c
      if (scaleb1.gt.1d0) then
         CALL PFTOPDG(xb,SCALEB1,XPQB1)
c         ash1=alphas2(rscaleb1)
      else
         call fillxpq0(XPQB1)
         ash1=0d0           ! 10 would fine as well
      endif
      if ((scaleb2.gt.1d0).and.(rscaleb2.gt.1d0)) then
         CALL PFTOPDG(xb,SCALEB2,XPQB2)
         ash2=alphas2(rscaleb2)
      else
         call fillxpq0(XPQB2)
         ash2=0d0           ! 10 would fine as well
      endif
      if (scaleb3.gt.1d0) then
         CALL PFTOPDG(xa,SCALEB3,XPQB3)
c         ash3=alphas2(rscaleb3)
      else
         call fillxpq0(XPQB3)
         ash3=0d0           ! 10 would fine as well
      endif
      if ((scaleb4.gt.1d0).and.(rscaleb4.gt.1d0)) then
         CALL PFTOPDG(xa,SCALEB4,XPQB4)
         ash4=alphas2(rscaleb4)
      else
         call fillxpq0(XPQB4)
         ash4=0d0           ! 10 would fine as well
      endif
c
      DO I=-6,6
         XPQ1(i)=XPQ1(i)/xa
         XPQ2(i)=XPQ2(i)/xa
         XPQ3(i)=XPQ3(i)/xb
         XPQ4(i)=XPQ4(i)/xb
         XPQB1(i)=XPQB1(i)/xb
         XPQB2(i)=XPQB2(i)/xb
         XPQB3(i)=XPQB3(i)/xa
         XPQB4(i)=XPQB4(i)/xa
      ENDDO
c
      if (icollide .le. 2) then ! f1 is proton, f2 is antiproton
         do i=-6,6
            f1(i)=xpq1(i)
            f2(i)=xpq2(i)
            f3(i)=xpq3(-i)
            f4(i)=xpq4(-i)
            fb1(i)=xpqb1(-i)
            fb2(i)=xpqb2(-i)
            fb3(i)=xpqb3(i)
            fb4(i)=xpqb4(i)
         enddo
      elseif (icollide .eq. 3) then !PP
         do i=-6,6
            f1(i)=xpq1(i)
            f2(i)=xpq2(i)
            f3(i)=xpq3(i)
            f4(i)=xpq4(i)
            fb1(i)=xpqb1(i)
            fb2(i)=xpqb2(i)
            fb3(i)=xpqb3(i)
            fb4(i)=xpqb4(i)
         enddo
      endif
c
c Matrix Element Calls
c
c
c test
c
c      print *,'f:',f1,f2,f3,f4
c      print *,'fb:',fb1,fb2,fb3,fb4
c      print *,asl1,asl3,ash2,ash4

c quark-initiated
      itest=is5+i15+i25+i35
      if (itest.eq.0) then
c *** change luminosities... ***
         call setlumt(f1,f3,fb1,fb3,xl1)
         call setlumt(f2,f4,fb2,fb4,xl2)
c  (u,c) (d,s,b) -> W+ -> d t g
         if(itop.eq.1)then
            alph1=asl1*xl1(1)
            alph2=ash2*xl2(1)
            tmp1=tpsi1(s12,t13,t14,t15,t23,t24,t25,s34,s35,s45,xm2,wm2)
c  (d~,s~,b~) (d,s,b) -> W+, s <-> u
            alph1=asl1*xl1(2)
            alph2=ash2*xl2(2)
            tmp2=tpsi1(t23,t13,s34,s35,s12,t24,t25,t14,t15,s45,xm2,wm2)
c  b  u  -> W+, t <-> u
            alph1=asl3*xl1(3)
            alph2=ash4*xl2(3)
            tmp3=tpsi1(s12,t23,t24,t25,t13,t14,t15,s34,s35,s45,xm2,wm2)
c  b  d~ -> W+, (3) s <-> t
            alph1=asl3*xl1(4)
            alph2=ash4*xl2(4)
            tmp4=tpsi1(t13,t23,s34,s35,s12,t14,t15,t24,t25,s45,xm2,wm2)
c
            tot=tot+tmp1+tmp2+tmp3+tmp4
         endif
c  u~ b~ -> W- = tmp1
         if(iatop.eq.1)then
            alph1=asl1*xl1(5)
            alph2=ash2*xl2(5)
            tmp1=tpsi1(s12,t13,t14,t15,t23,t24,t25,s34,s35,s45,xm2,wm2)
c  d  b~ -> W- = tmp2
            alph1=asl1*xl1(6)
            alph2=ash2*xl2(6)
            tmp2=tpsi1(t23,t13,s34,s35,s12,t24,t25,t14,t15,s45,xm2,wm2)
c  b~ u~ -> W- = tmp3
            alph1=asl3*xl1(7)
            alph2=ash4*xl2(7)
            tmp3=tpsi1(s12,t23,t24,t25,t13,t14,t15,s34,s35,s45,xm2,wm2)
c  b~ d  -> W- = tmp4
            alph1=asl3*xl1(8)
            alph2=ash4*xl2(8)
            tmp4=tpsi1(t13,t23,s34,s35,s12,t14,t15,t24,t25,s45,xm2,wm2)
c
            tot=tot+tmp1+tmp2+tmp3+tmp4
         endif
      endif

c gluon-initiated / light line
       itest=i15+i25
      if (itest.eq.0) then
c  (u,c) g -> d t b~
c         alph1=asl1           ! dummy
         alph2=ash2
         tmp1=tpsi2(s12,t13,t14,t15,t23,t24,t25,s34,s35,s45,xm2,wm2)
c  (d~,s~,b~) g -> u~ t b~
         tmp2=tpsi2(t23,t13,s34,s35,s12,t24,t25,t14,t15,s45,xm2,wm2)
c  g u  -> d t b~
c         alph1=asl3           ! dummy
         alph2=ash4
         tmp3=tpsi2(s12,t23,t24,t25,t13,t14,t15,s34,s35,s45,xm2,wm2)
c  g d~ -> u~ t b~
         tmp4=tpsi2(t13,t23,s34,s35,s12,t14,t15,t24,t25,s45,xm2,wm2)
c
c  g=0,d=1,u=2,s=3,c=4,b=5,anti = -n
c         ckm1 = vtd2+vts2+vtb2                ! (d,s,b) + W -> t or reverse
         ckm1=1d0                              ! more accurate
         if(itop.eq.1) tot=tot+ckm1*((f2(2)+f2(4))*fb2(0)*tmp1 +
     &        (f2(-1)*(vud2+vcd2)+ f2(-3)*(vus2+vcs2)+f2(-5)*(vub2+
     &        vcb2))*fb2(0)*tmp2 + (f4(2)+f4(4))*fb4(0)*tmp3 + (f4(-1)
     &        *(vud2+vcd2)+f4(-3)*(vus2+vcs2)+f4(-5)*(vub2+vcb2))*
     &        fb4(0)*tmp4 )                    ! more accurate
c         print *,'vud2+vus2+vub2',vud2+vus2+vub2
c         print *,'vcd2+vcs2+vcb2',vcd2+vcs2+vcb2
c         print *,'vud2+vcd2',vud2+vcd2
c         print *,'vus2+vcs2',vus2+vcs2
c         print *,'vub2+vcb2',vub2+vcb2
c         stop
         if(iatop.eq.1)tot=tot+ckm1*( (f2(-2)+f2(-4))*fb2(0)*tmp1 +
     &        (f2(1)*(vud2+vcd2)+ f2(3)*(vus2+vcs2)+f2(5)*(vub2+vcb2))
     &        *fb2(0)*tmp2 + (f4(-2)+f4(-4))*fb4(0)*tmp3 + (f4(1)*
     &        (vud2+vcd2)+f4(3)*(vus2+vcs2)+f4(5)*(vub2+vcb2))*
     &        fb4(0)*tmp4 )                    ! more accurate
cz Added to pull out b fraction 8/21/02
         if(itop.eq.1) btmp=btmp+ckm1*((f2(2)+f2(4))*fb2(0)*tmp1 +
     &        (f2(-1)*(vud2+vcd2)+ f2(-3)*(vus2+vcs2)+f2(-5)*(vub2+
     &        vcb2))*fb2(0)*tmp2 + (f4(2)+f4(4))*fb4(0)*tmp3 + (f4(-1)
     &        *(vud2+vcd2)+f4(-3)*(vus2+vcs2)+f4(-5)*(vub2+vcb2))*
     &        fb4(0)*tmp4 )                    ! more accurate
         if(iatop.eq.1)btmp=btmp+ckm1*( (f2(-2)+f2(-4))*fb2(0)*tmp1 +
     &        (f2(1)*(vud2+vcd2)+ f2(3)*(vus2+vcs2)+f2(5)*(vub2+vcb2))
     &        *fb2(0)*tmp2 + (f4(-2)+f4(-4))*fb4(0)*tmp3 + (f4(1)*
     &        (vud2+vcd2)+f4(3)*(vus2+vcs2)+f4(5)*(vub2+vcb2))*
     &        fb4(0)*tmp4 )                    ! more accurate
      endif

c gluon-initiated / heavy line
      itest=i13+i15+i23+i25
      if (itest.eq.0) then
c  g b(s,d) -> (d,s,b) t (u~,c~,NOT t~ here)
         alph1=asl1
c         alph2=ash1           ! dummy
         tmp1=tpsi3(s12,t13,t14,t15,t23,t24,t25,s34,s35,s45,xm2,wm2)
c  b g -> (d,s,b) t (u~,c~,NOT t~ here)
         alph1=asl3
c         alph2=ash2           ! dummy
         tmp2=tpsi3(s12,t23,t24,t25,t13,t14,t15,s34,s35,s45,xm2,wm2)
c
c         ckm1 = (vud2+vcd2) + (vus2+vcs2) + (vub2+vcb2) ! (d,s,b) or (u,c)
         ckm1 = 2d0                                      ! more accurate
         if(itop.eq.1) tot=tot+f1(0)*(vtd2*fb1(1)+vts2*fb1(3)+
     &        vtb2*fb1(5))*ckm1*tmp1 + f3(0)*(vtd2*fb3(1)+vts2*fb3(3)+
     &        vtb2*fb3(5))*ckm1*tmp2
         if(iatop.eq.1) tot=tot+f1(0)*(vtd2*fb1(-1)+vts2*fb1(-3)+
     &        vtb2*fb1(-5))*ckm1*tmp1 + f3(0)*(vtd2*fb3(-1)+vts2*
     &        fb3(-3)+vtb2*fb3(-5))*ckm1*tmp2
      endif
c
c Add safety valve for bad Monte Carlo points
c
c      if (dabs(tot).gt.1d3) then
c         print *,'OOPS!'
c         print *,dsqrt(s12),dsqrt(-t13),dsqrt(-t14),dsqrt(-t15),
c     &        dsqrt(-t23),dsqrt(-t24),dsqrt(-t25),dsqrt(s34),
c     &        dsqrt(s35),dsqrt(s45)
c         tot=0d0
c      endif
c
cz Calculate final b fraction here: (Note this probably fails at pt~0, since
cz   subtraction is not accounted for.)
      if(dabs(tot).gt.1d-25) bwgt=btmp/tot
cz      print *,btmp,bwgt,tot

      return
      end

      double precision function tpsi1(s12,t13,t14,t15,t23,t24,t25,
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
      double precision alph1,alph2
      common /sctopsi/ alph1,alph2

      t14p = t14 - xm2
      t24p = t24 - xm2
      s34p = s34 - xm2
      s45p = s45 - xm2

      fq1 = 1d0 / ( t24 - wm2 ) / ( t24 - wm2 )*alph1
      fq2 = 1d0 / ( t13 - wm2 ) / ( t13 - wm2 )*alph2
      tpsi1 = fq1 * ( - 2d0*s12*t13/t15*s34p/s35 - s12*t13/t15*s45p/s35 
     . + s12*t14p/t15 - s12/t15*s34p - s12*s34p/s35 - s12*s45p/s35 
     . - t13/t15*s34p/s35*t25 - s34p*t25/t15 + t23*s34p/s35 )
     . + fq2 * ( 2d0*xm2*s12*s34p/s45p/s45p + 2d0*xm2*s12/s45p/s45p*s35 
     . + s12*t23/t25 - 2d0*s12*t24p*s34p/s45p/t25 
     . - s12*t24p/s45p*s35/t25 - s12*s34p/s45p - s12*s34p/t25 
     . - s12/s45p*s35 + t14p*s34p/s45p - t15*t24p*s34p/s45p/t25 
     . - t15*s34p/t25 )
      tpsi1=tpsi1*2d0*cf
      return
      end

      double precision function tpsi2(s12,t13,t14,t15,t23,t24,t25,
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
      double precision alph1,alph2
      common /sctopsi/ alph1,alph2

      t14p = t14 - xm2
      t24p = t24 - xm2
      s34p = s34 - xm2
      s45p = s45 - xm2

      tpsi2 = -2d0*xm2*t15/t24p/t24p * ( t23 + s34p ) 
     .     + s12/t24p/t25*s34p*s45p + s12/t25*s34p 
     .     + t15/t24p/t25*s45p* ( t23 + 2d0*s34p )
     .     + t15/t24p * ( t23 + s34p ) + t15/t25 * ( s34p - s35 )
     .     - t14p*s34p/t24p
      tpsi2=tpsi2/(t13-wm2)/(t13-wm2)*alph2
      return
      end


      double precision function tpsi3(s12,t13,t14,t15,t23,t24,t25,
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
      double precision alph1,alph2
      common /sctopsi/ alph1,alph2

      t14p = t14 - xm2
      s34p = s34 - xm2
      s45p = s45 - xm2

      tpsi3 = s12/t15 * ( s34p*s35/t13 + s34p ) 
     .     + t25/t15 * ( t14p*s35/t13 
     .     + 2d0*s34p*s35/t13 + s34p - s45p ) 
     .     + t25/t13 * ( t14p + s34p ) 
     .     - t23*s34p/t13
      tpsi3=tpsi3/(t24-wm2)/(t24-wm2)*alph1
      return
      end

      subroutine fillxpq0(xpq)
      implicit none
      double precision xpq(-6:6)
      integer i
      do i=-6,6
         xpq(i)=0d0
      enddo
      return
      end
