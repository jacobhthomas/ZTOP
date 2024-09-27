       double precision function dilog(x)
       implicit double precision (a-z)
       parameter (pi6=1.644934066848226d+00)
       parameter (een=1.d+00)
       parameter (vier=0.25d+00)
       parameter (b2=+2.7777777777777778D-02)
       parameter (b3=-2.7777777777777778D-04)
       parameter (b4=+4.7241118669690098D-06)
       parameter (b5=-9.1857730746619641D-08)
       parameter (b6=+1.8978869988971001D-09)
       parameter (b7=-4.0647616451442256D-11)
       parameter (b8=+8.9216910204564523D-13)
1      if(x.lt.0.d0)go to 3
       if(x.gt.0.5d0)go to 4
       z=-dlog(1.d0-x)
7      z2=z*z
       dilog=z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b8+b7)+b6)
     1 +b5)+b4)+b3)+b2)+een)-z2*vier
       if(x.gt.een)dilog=-dilog-.5d0*u*u+2.d0*pi6
       return
3      if(x.gt.-een)go to 5
       y=een/(een-x)
       z=-dlog(een-y)
       z2=z*z
       u=dlog(y)
       dilog=z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b8+b7)+b6)
     1 +b5)+b4)+b3)+b2)+een)-z2*vier-u*(z+.5d0*u)-pi6
       return
4      if(x.ge.een)go to 10
       y=een-x
       z=-dlog(x)
6      u=dlog(y)
       z2=z*z
       dilog=-z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b8+b7)+b6)
     1 +b5)+b4)+b3)+b2)+een-u)+z2*vier+pi6
       if(x.gt.een)dilog=-dilog-.5d0*z*z+pi6*2.d0
       return
5      y=een/(een-x)
       z=-dlog(y)
       z2=z*z
       dilog=-z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b8+b7)+b6)
     1 +b5)+b4)+b3)+b2)+een)-z2*vier
       return
10     if(x.eq.een)go to 20
       xx=1.d0/x
       if(x.gt.2.d0)go to 11
       z=dlog(x)
       y=1.d0-xx
       go to 6
11     u=dlog(x)
       z=-dlog(1.d0-xx)
       go to 7
20     dilog=pi6
       return
       end
