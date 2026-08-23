c  page- 287, 292-Newton's Interpolation method better then lagrange c                 interpolation method

      Dimension t(20),q(20),a(20),d(20,20),q1(200),tpp(100)
      real a,b

      open(unit= 4, file='qt.dat', status='old')
      open(unit= 6, file='qt1.dat', status='old')
      open(unit= 5, file='qt.out', status='new')

      write(*,*)'input size of the system, data at i=1 to i=11
     + total=(11-1)+1= at 11 data points '

      read(*,*) n                        

      do 1 i=1,n

      read(4,*) q(i)
      read(6,*) t(i)
      
1     continue
c-------------------------------------------------------------------
c     construction of the difference table d(i,j) for j=1

c     the first column (j=1) of the d(20,20) matrix are nothing but the given
c     values of the Q (discharge) in the qt.dat file for different time
c     intervals as given in the qt1.dat file

      do 2 i=1,n
      d(i,1)= q(i)
2     continue

c     construction of the difference table d(i,j) for j=2 to n

      do 3 j=2,n
      do 4 i=1,n-(j-1)

      d(i,j)=( d(i+1,j-1)-d(i,j-1) ) /(  t(i+  j-1) - t(i) )

4     continue
3     continue
c------------------------------------------------------------------------
c     The a(j) matric are nothing but the first row of the D(i,j) matrix

      do 5 j=1,n
      a(j)= d(1,j)
5     continue
c-------------------------------------------------------------------------
c      write(*,*)'pt where interpolation is required'
c      read(*,*) tp
           
c------------------------------------------
      write(*,*)'Number of points where interpolation value required'
      read(*,*) m
      write(*,*)'write the values of time tp in the above points, where
     + interpolated values of discharge Q (i.e f(x) ) is required'

      do 111 k=1,m
      read(*,*) tpp(k)
      tp=tpp(k)

c------------------------------------------
      sum=a(1)
      do 6 i=2,n
c----------------------
      p=1.0
      do 7 j=1,i-1

      p = p * ( tp-t(j) )

7     continue
c-----------------------
      sum=sum+a(i)*p
6     continue

      write(5,*) 'interpolated values for these time intervals=' , sum

             
111   continue

      
c-------------------------------------------------------
      stop

      end
