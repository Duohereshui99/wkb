ccccccc      
      module parameter    !全局参数,也就是各个（物理）常数
            implicit none
            real*8,parameter :: hbarc=197.3269718d0  !hbar      ! NIST Ref 02.12.2014   ! MeV.fm           
            real*8,parameter :: finec=137.03599d0
            real*8,parameter :: amu=931.49432d0      !MeV
            real*8,parameter :: e2=1.43997d0         !MeV.fm
            real*8,parameter :: PI=acos(-1.0)        !圆周率pi
      end module