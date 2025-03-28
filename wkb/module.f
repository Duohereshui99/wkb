ccccccc
      module system   !information about parent and daughter nucleus
            implicit none          
            real*8::z_1             !charge
            real*8::z_2
            real*8::z12             !z1*z2
            real*8::m_1          !mass number,not MeV
            real*8::m_2
            real*8::mass_1
            real*8::mass_2          !mass=mass number*amu+mass excess
            real*8::mass_excess_1   !mass excess (MeV)
            real*8::mass_excess_2
            real*8::mu              !reduced mass num
            real*8::L               !angular momentum num
            real*8::Q               !Q value
            real*8::P               !preformation factor
            real*8::G               !global quantum num
            real*8::E               !decay kinetic energy Eα
            real*8::gamma           !decay width
            real*8::t_half          !half life
      end module
ccccccc
      module potential
            implicit none
            real*8::v0          !!depth
            real*8::r0          !!
            real*8::a           !!
      end module
ccccccc
      module variable
            implicit none
            real*8,allocatable::rr(:)   !radial coordinate
            integer,allocatable::r(:)    !index of r_1,r_2,r_3, ...
            real*8,allocatable::fr(:)   !!f=k(r)=2mu/h^2(Q-V(r))=k^2
            real*8::F                 !!normalization factor
      end module
ccccccc
      module mesh 
            implicit none
            real*8::hcm
            real*8::r1
            integer::n
      end module