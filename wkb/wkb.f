ccccccc
      program main
        use parameter
        use pot
        use system
        use potential
        use variable
        use mesh
        implicit none
        integer::i,j,k
        real*8::s,t
        real*8::t1,t2
        !4He+209Bi
        call get_info()
        call cpu_time(t1)
        n=1000000
        mass_1=4d0   
        mass_2=207d0
        z_1=2d0
        z_2=82d0
        L=5d0
        v0=162.3d0 !MeV
        a=0.4 !fm
        P=0.03d0
        Q=7.599d0!9.815d0 !MeV
        hcm=0.0001d0 !fm
        r1=0.8d0 !fm
        r0=7.642d0
ccccccc
        z12=z_1*z_2
        mu=amu*mass_1*mass_2/(mass_1+mass_2)
ccccccc
        write(*,*) 'Q value:',Q,'MeV'
ccccccc
        allocate(rr(n))
        allocate(fr(n))
        allocate(r(10))
ccccccc
        G=21d0
ccccccc
        do i=1,n
          rr(i)=r1+hcm*i
        end do
ccccccc
        do i=1,n        !!k^2=2mu/h^2(Q-V(r))=k^2,|k|=sqrt(abs(k^2))
          fr(i)=kr(Q,mu,v0,a,r0,z12,l,rr(i))
          write(33,*) rr(i),fr(i)
        end do
ccccccc
        k=1
        do i=1,n-1      !!where Q=V(r),r(i) -> rr(r(i)) -> r_i
          if(fr(i)*fr(i+1)<0) then 
            r(k)=i
            k=k+1
          end if
        end do

ccccccc
        s=0d0
        do i=r(1),r(2)    !!\int_{r_1}^{r_2}dr\sqrt{2\mu/\hbar^2(Q-V(r))}=\int |k(r)|dr
          s=s+hcm*sqrt(abs(fr(i)))
        end do

ccccccc
         write(*,*) 'r0=',r0
         write(*,*) rr(r(1)),rr(r(2)),rr(r(3))
         write(*,*) 'the number of zeroes where Q=V(r):',k-1

ccccccc
          F=0d0 
          do i=r(1),r(2)
            s=0d0 
            do j=r(1),i           !!int_r1^r dr'
              s=s+hcm*sqrt(abs(fr(j)))
            end do                                    
            F=F+hcm*(cos(s-pi/4d0))**2/sqrt(abs(fr(i)))     !!int_r1^r2 dr
          end do
        !   do i=r(1),r(2)
        !     F=F+hcm/2d0/sqrt(abs(fr(i)))
        !   end do
          !!F goes reciprocally to get normalization factor F
          F=1d0/F
          write(*,*) 'normalization factor F:',F
ccccccc
          t=0d0
          do i=r(2),r(3)
            t=t+hcm*sqrt(abs(fr(i)))
          end do                  !!width gamma
ccccccc
          write(*,*) 't:',t
          write(*,*) 'exp(-2d0*t):',exp(-2d0*t)
          write(*,*) 'P*F*hbarc**2/4d0/mu:',P*F*hbarc**2/4d0/mu
ccccccc
          gamma=P*F*hbarc**2/4d0/mu*exp(-2d0*t)
          t_half=hbarc*log(2d0)/gamma                   !!fm
          t_half=t_half/3d0/1e23                            !!s
ccccccc
          write(*,*) 'gamma:',gamma,'MeV'
          write(*,*) 'half life T:',t_half,'s'
          write(*,*) 'half life T:',t_half/3600d0/24d0/365d0,'yr'
ccccccc
        deallocate(rr,r,fr)
ccccccc
        call cpu_time(t2)
        write(*,*) 'total cputime:',t2-t1
ccccccc
        contains
ccccccc
            subroutine get_info()
            ! 使用预处理器检查宏是否被定义
#ifdef BASE
        print *, 'Base directory: ', BASE
#endif

#ifdef VERDATE
        print *, 'Version date: ', VERDATE
#endif

#ifdef VERREV
        print *, 'Version revision: ', VERREV
#endif

#ifdef COMPDATE
        print *, 'Compilation date: ', COMPDATE
#endif
            end subroutine
      end program main