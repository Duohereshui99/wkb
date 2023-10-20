ccccccc
      module pot     
        use parameter
        implicit none
        contains
!
!!Nuclear potential: depth v0, diffuseness a, radius parameter r0
!
        real*8 function vpotnn(v0,a,r0,r)
        real*8::v0,a,r0,r
        vpotnn=-v0*(1d0+cosh(r0/a))/(cosh(r/a)+cosh(r0/a))
        end function vpotnn
!
!!Coulomb potential: z12=z1*z2, radius parameter r0
!
        real*8 function vpotcoul(z12,r0,r)
        real*8::z12,r0,r
        if (r>=r0) then
            vpotcoul=z12*e2/r
        else
            vpotcoul=z12*e2/2/r0*(3d0-r**2/r0**2)
        end if 
        end function vpotcoul
!
!!Centrifugal barrier potential: orbital angular momentum l, reduced mass mu
!
        real*8 function vpotcent(mu,l,r)
        real*8::mu,l,r
        vpotcent=hbarc**2/2d0/mu*(l+0.5d0)**2/r**2
        end function vpotcent
!
!!alpha core potential in the cluster model
!
        real*8 function vpot(v0,a,r0,z12,mu,l,r)
        real*8::v0,a,r0,z12,mu,l,r
        vpot=vpotnn(v0,a,r0,r)+vpotcoul(z12,r0,r)
     &   +vpotcent(mu,l,r)
        end function vpot
!
!!kr=2mu/h^2(Q-V(r))=k^2,|k|=sqrt(abs(k^2))
!
        real*8 function kr(Q,mu,v0,a,r0,z12,l,r)
        real*8::Q,mu,v0,a,r0,z12,l,r
        kr=2d0*mu/hbarc**2*(Q-vpot(v0,a,r0,z12,mu,l,r))
        end function kr

      end module