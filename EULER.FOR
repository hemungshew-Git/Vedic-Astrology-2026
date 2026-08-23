      subroutine euler(x0,y0,h,n,abcissae,ordinatae)
c
c     INPUT:
c     (x0,y0) = starting point
c     h       = net spacing
c     n       = number of points (i=0,1,...n)
c
c     OUTPUT:
c     (abcissae,ordinatae) = solution of difequ
c
      implicit double precision (a-h,o-z)
      double precision abcissae(0:n),ordinatae(0:n)
      x=x0
      y=y0
      abcissae(0)=x0
      ordinatae(0)=y0
      do, i=1,n
          y=y+h*x*y**(1.d0/3.d0)
          x=x0+h*dble(i)
          abcissae(i)=x
          ordinatae(i)=y
      enddo
      return
      end

