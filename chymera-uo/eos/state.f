      subroutine State() ! this routine calculates the pressure
        use eos
        implicit none
     
#include "hydroparam.h"
#include "globals.h"
#include "units.h"

        integer JSTART,J,K,L
        real*8::mu

        if (jmin.gt.2) then
           jstart=jmin2
        else   
           jstart=1
        endif

C$OMP DO SCHEDULE(STATIC) PRIVATE(mu)
      DO L=1,LMAX
         DO K=1,KMAX2
            DO J=jstart,jmax2

c              call get_gamma(eps(J,K,L),rho(J,K,L),tempk(J,K,L),
c     &             mu,gamma)
c              p(J,K,L) = bkmpcgs*rho(J,K,L)*tempk(J,K,L)*rhoconv
c     &            / (mu*pconv)
               p(j,k,l) = eps(j,k,l)*(gamma - one)
            end do
         ENDDO
      ENDDO
C$OMP END DO
c        do l =1,8
c           do k = 1,4
c              do j = 1,7
c        write(97,*) p(j,k,l), (gamma-1)*eps(j,k,l)
c              end do
c           end do
c        end do

      RETURN
      END
 

