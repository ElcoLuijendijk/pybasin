!aft.f90

subroutine reduced_ln(timestepduration, temperature, &
                      rmr0, kappa, alpha, C0, C1, C2, C3, &
                      rcf, nsteps)

!     fission track annealing algorithm, using Ketcham 2005, 2007 eqs.
!     calculate reduced track lengths
   
    implicit none

    !integer, intent(in) :: nsteps
    integer :: nsteps
    
    real*8, intent(in) :: timestepduration(nsteps)
    real*8, intent(in) :: temperature(nsteps)
    real*8, intent(in) :: rmr0
    real*8, intent(in) :: kappa
    
    real*8, intent(in) :: alpha
    real*8, intent(in) :: C0
    real*8, intent(in) :: C1
    real*8, intent(in) :: C2
    real*8, intent(in) :: C3
    
    real*8, intent(out) :: rcf(nsteps)
    
    integer segments, timestepcount
    
    real*8 gf(nsteps), teq(nsteps), rcmod(nsteps)
    real*8 tt, temp1
    real*8 f1, f2

    do segments=1, nsteps
        gf(:) = 0
        teq(:) = 0
        do timestepcount=segments, nsteps
            tt = timestepduration(timestepcount)
            temp1 = temperature(timestepcount)
            
            if (timestepcount.gt.segments) then
                teq(timestepcount) = exp( ((gf(timestepcount-1) - C0)/C1 & 
                * (log(1.0/temp1)-C3)) + C2 )
            else
                teq(timestepcount) = 0
            endif
            
            gf(timestepcount) = ( C0 + C1 * ((log(tt+teq(timestepcount))-C2) &
            / (log(1./temp1)-C3)) )

        end do  
            
    !   f1 = math.pow(g[-1], (1./alpha))
        f1 = (gf(nsteps)**(1.0/alpha))
        f2 = f1+1
        rcmod(segments) = 1.0 / f2
        
        rcf(segments) =  ((rcmod(segments) -rmr0)/ (1.0-rmr0)) ** kappa
        
        if (rcmod(segments).lt.rmr0) then
            rcf(segments) = 0
        endif

    end do

end subroutine reduced_ln

!C END FILE AFT.F

            
        
        
            
         
            
       
        
  
