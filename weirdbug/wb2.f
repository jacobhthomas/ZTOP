      subroutine setblock
      implicit double precision (A-H, O-Z)

      common / block1 / var1, var2(0:9)

      var1 = .24d0
      do i = 0, 9
         var2(i) = .1d0*i + .26d0
      enddo

      return
      end
      
