!**********************************************************************
! Sample FORTRAN Program from "First Course in Fourier Analysis"      !
!                               by D.W. Kammler                       !
!**********************************************************************
! FORTRAN code for a radix 2 FFT
      subroutine fft(fr,fi,mp,isw)
! This FORTRAN subroutine computes the discrete Fourier transform
!
!      ft(k+1)=scale*sum f(j+1)*exp(isn*2*pi*i*j*k/n), k=0,1,...n-1
!
! of the complex array
!
!      f(j+1)=fr(j+1)+i*fi(j+1), j=0,1,...n-1 where
!
!      n    =2**mp
!      i    =sqrt(-1)
!      isn  =+1 is isw.gt.0
!           =-1 if isw.lt.0
!      scale=1         if isw=-1 or +1
!           =1/n       if isw=-2 or +2
!           =1/sqrt(n) if isw=-3 or +3
!
! Computations are done in place, and at the time of return the real
! arrays fr,fi have been overwritten with the real and imaginary parts
! of the desired Fourier transform
!
! The code is based on the presentation of the FFT given in the
! text 'A First Course in Fourier Analysis' by David W. Kammler.
!
      implicit real*8 (a-h,o-z)
      parameter(zero=0.d0,one=1.d0,two=2.d0)
!
! To change the code to single precision, delete the above implicit
! real*8 statement, delete the d0's in the above parameter statement
! that defines zero,one,two, and replace dsqrt by sqrt in the equations
! defining rrootn, rroot2, and in the first equation of the do 6 loop.
!
      parameter(maxmp=10,maxds=257,maxdir=32)
      dimension fr(*),fi(*),s(maxds),ir(maxdir)
!
! To accomodate a larger vector f, insert a larger value of maxmp and
! corresponding values of maxds=2**(maxmp-2)+1 and
! maxdir=2**[(maxmp+1)/2] in the above parameter statement.
!
      data lastmp,maxs/0,0/
      isws=isw**2
      if((mp.lt.0).or.(mp.gt.maxmp).or.(isws.eq.0).or.(isws.gt.9)) then
          write(*,2)
2         format(' Improper argument in FFT Subroutine')
          stop
      endif
      if(mp.eq.0) return
      if(mp.eq.lastmp) then
          if(mp.eq.1) go to 20
          if(mp.gt.1) go to 14
      endif
!
! At the time of the first call of the subroutine with a give value of
! mp.gt.0, initialize various constants and arrays that depend on mp
! but not fr,fi. On subsequent calls with the same mp, bypass this
! initialization process.
!
      lastmp=mp
      n     =2**mp
      temp  =n
      rn    =one/temp
      rrootn=dsqrt(rn)
      nh    =n/2
      if(mp.eq.1) go to 20
      if(mp.eq.2) go to 8
      if(mp.lt.maxs) then
          nsn=2**(maxs-mp)
      else
          nsn=1
      endif
!
! When mp.le.maxs, the spacing parameter nsn is used to retreive sine
! values from a previously computed table.
!
      if(mp.le.maxs) go to 8
! When n=8,16,32,... and mp.gt.maxs, precompute
! s(j+1)=sin((pi/2)*(j*nq)), j=0,1,2,...nq, for use in
! the subsequent calculations.
!
      nq     =n/4
      ne     =n/8
      rroot2 =one/dsqrt(two)
      maxs   =mp
      s(1)   =zero
!            =sin(0*pi/4)
      s(ne+1)=rroot2
!            =sin(1*pi/4)
      s(nq+1)=one
!            =sin(2*pi/4)
      if(mp.eq.3) go to 8
        h=rroot2
!        =.5*sec(pi/4)
!
! Pass from the course grid to a finer one by using the trig identity
!
!            sin(a)=(.5*sec(b))*(sin(a-b)+sin(a+b)).
!
! (This clever idea is due to O.Buneman, cf.   SIAM J. SCI. STAT.
! COMPT. 7 (1986), pp. 624-638.)
!
      k=ne
      do 6 i=4,mp
         h  =one/dsqrt(two+one/h)
!           =half secant of half the previous angle
         kt2=k
         k  =k/2
!           =n/2**i
         do 4 j=k,nq,kt2
4           s(j+1)=h*(s(j-k+1)+s(j+k+1))
6        continue
8     continue
!
! Prepare a short table of bit reversed integers to use in the
! subsequent bit reversal permutation.
!
      muplus=(mp+1)/2
      m     =1
      ir(1) =0
      do 12 nu=1,mplus
          do 10 k=1,m
              it      =2*ir(k)
              ir(k)  =it
10        ir(k+m)=it+1
12    m=m+m
!
! If mp is odd, then m=m/2.
!
      itemp=2*(mp/2)
      if(itemp.ne.mp) then
          m=m/2
      endif
!
! The parameters n,nh,nq,ne,m,lastmp,maxs,nsn,rn,rrootn,rroot2 and the
! arrays ir(*), s(*) are now suitably initialized.
!
14    continue
! Apply the bit reversal permutation to the array f using the
! Bracewell-Buneman scheme.
!
      do 18 iq=1,m-1
          npr =iq-m
          irpp=ir(iq+1)*m
          do 16 ip=0,ir(iq+1)-1
              npr     =npr+m
              npr1    =npr+1
              irp1    =irpp+ir(ip+1)+1
              tempr   =fr(npr1)
              tempi   =fr(npr1)
              fr(npr1)=fr(irp1)
              fi(npr1)=fi(irp1)
              fr(irp1)=tempr
16            fi(irp1)=tempi
18    continue
20    continue
!
! Apply the mp Q-matrices
!
! Carry out stage 1 of the FFT using blocks of size 2x2 and apply
! the desired scale factor.
!
      if(isws.eq.1) then
          scale=one
      elseif(isws.eq.4) then
          scale=rn
      else
          scale=rrootn
      endif
      do 22 k=0,nh-1
          k1    =2*k+1
          k2    =k1+1
          tempr =(fr(k1)-fr(k2))*scale
          tempi =(fi(k1)-fi(k2))*scale
          fr(k1)=(fr(k1)+fr(k2))*scale
          fi(k1)=(fi(k1)+fi(k2))*scale
          fr(k2)=tempr
22        fi(k2)=tempi
      if(mp.eq.1) return
!
! Carry out stages 2,3,...,mp of the FFt using blocks
! of size mxm = 4x4,8x8,...nxn.
!
      mcap=1
      kcap=n/4
      do 32 mu=2,mp
!
! At this point mcap=2**(mu-2) and kcap=2**(mp-mu).
!
! Deal first with the quadruplet of components where sin=0 or cos=0.
!
      do 24 k=0,kcap-1
          k0=k*4*mcap+1
          k1=k0+mcap
          k2=k0+2*mcap
          k3=k0+3*mcap
          tempr =fr(k0)-fr(k2)
          tempi =fi(k0)-fi(k2)
          fr(k0)=fr(k0)+fr(k2)
          fi(k0)=fi(k0)+fi(k2)
          fr(k2)=tempr
          fi(k2)=tempi
          if(isw.lt.0) then
              fr(k3)=-fr(k3)
              fi(k3)=-fi(k3)
          endif
          temp1 =fr(k1)+fi(k3)
          temp2 =fi(k1)-fr(k3)
          fr(k1)=fr(k1)-fi(k3)
          fi(k1)=fi(k1)+fr(k3)
          fr(k3)=temp1
          fi(k3)=temp2
24    continue
      if(mcap.eq.1) go to 30
!
! Now deal with the remaining mcap-1 quadruplets of components where sin,
! cos are both nonzero.
!
      do 28 lamda=1,mcap-1
          indx =nsn*lamda*kcap+1
          sn   =s(indx)
!              =sin((2*pi)*(lamda/mcap))
          indx =nq-indx+2
          cs   =s(indx)
!              =cos((2*pi)*(lamda/mcap))
          if(isw.lt.0) then
              sn=-sn
          endif
          do 26 k=0,kcap-1
              k4m=k*4*mcap+1
              k0 =k4m+lamda
              k1 =k4m+2*mcap-lamda
              k2 =k4m+2*mcap+lamda
              k3 =k4m+4*mcap-lamda
              r1    =cs*fr(k2)-sn*fi(k2)
              r2    =cs*fi(k2)+sn*fr(k2)
              temp1 =fr(k0)-r1
              temp2 =fi(k0)-r2
              fr(k0)=fr(k0)+r1
              fi(k0)=fi(k0)+r2
              r1    =cs*fr(k3)+sn*fi(k3)
              r2    =cs*fi(k3)-sn*fr(k3)
              temp1 =fr(k1)+r1
              temp2 =fi(k1)+r2
              fr(k1)=fr(k1)-r1
              fi(k1)=fi(k1)-r2
              fr(k3)=temp1
26            fi(k3)=temp2
28    continue
30    mcap=mcap*2
      kcap=kcap/2
32    continue
      return
      end