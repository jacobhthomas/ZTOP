      program simpletest
      real*8 v1, vout(8)
      integer i1,i

      i1 = 2
      v1 = 3.4d0

      vout(1) = 5.1d0
      vout(2) = 5.2d0

      call simple(i1, v1, vout)

      write(*,*) (vout(i), i=1, 8)
      
      end
      
