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
ccccccc
        namelist /systems/ mass_1,z_1,mass_2,z_2,L,Q
        namelist /meshs/ n,hcm,r1
        namelist /Pfactor/ P
ccccccc
        call get_info()
        call cpu_time(t1)
ccccccc
        read(5,nml=systems)
        read(5,nml=meshs)
        read(5,nml=Pfactor)
ccccccc

        v0=162.3d0 !MeV
        a=0.4 !fm
    !    P=0.35d0
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
        if (mass_1+mass_2-z_1-z_2>126d0) then
          G=24d0
        else if(mass_1+mass_2-z_1-z_2<=126d0.and.mass_1+mass_2-z_1-z_2>82) then
          G=22d0
        else 
          G=20d0
        end if
ccccccc
        G=23d0
ccccccc
        write(*,*) 'global quantum num G:',G
ccccccc
        do i=1,n
          rr(i)=r1+hcm*i
        end do
ccccccc
         s=0d0
         r0=mass_2**(1d0/3d0)*1.2d0!7.64d0
         do while(abs(s-(G-L+1d0)*pi/2d0)>1e-1)
ccccccc
        do i=1,n        !!k^2=2mu/h^2(Q-V(r))=k^2,|k|=sqrt(abs(k^2))
          fr(i)=kr(Q,mu,v0,a,r0,z12,l,rr(i))
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
         r0=r0+0.001d0    !!do while loop body
         end do
ccccccc
         write(*,*) 'r0=',r0
         write(*,*) rr(r(1)),rr(r(2)),rr(r(3))
         write(*,*) 'the number of zeroes where Q=V(r):',k-1
        write(*,*) 'the integration of the first well:',s
        write(*,*) '(G-L+1)*pi/2=',(G-L+1)*pi/2
ccccccc
          F=0d0 
        !   do i=r(1),r(2)
        !     s=0d0 
        !     do j=r(1),i           !!int_r1^r dr'
        !       s=s+hcm*sqrt(abs(fr(j)))
        !     end do                                    
        !     F=F+hcm*(cos(s-pi/4d0))**2/sqrt(abs(fr(i)))     !!int_r1^r2 dr
        !   end do
          do i=r(1),r(2)
            F=F+hcm/2d0/sqrt(abs(fr(i)))
          end do
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
          t_half=t_half*ratio!/3/1e23                            !!s
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