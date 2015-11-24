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

cpoly              call get_gamma2(eps(J,K,L),rho(J,K,L),tempk(J,K,L),
cpoly     &             mu,gamma1(J,K,L))
cpoly              p(J,K,L) = bkmpcgs*rho(J,K,L)*tempk(J,K,L)*rhoconv
cpoly     &            / (mu*pconv)
              p(J,K,L) = (gamma-1.0)*eps(j,k,l)

            end do
         ENDDO
      ENDDO
C$OMP END DO

      RETURN
      END
 

