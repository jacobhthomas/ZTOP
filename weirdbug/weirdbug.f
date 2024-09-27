      program weird
      implicit double precision (A-H, O-Z)

      common / block1 / var1, var2(0:9)

      call setblock

      write (*,*) var1
      write (*,*) (var2(i), i=0, 9)

      end
      
