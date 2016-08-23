!aft.f90

subroutine reduced_ln(timestepduration, temperature, &
                      rmr0, kappa, annealing_eq_f90, alpha, C0, C1, C2, C3, &
                      rcf, nsteps)

!   fission track annealing algorithm
!   uses eqs. by Ketcham 2005, 2007 to calculate reduced track lengths
!   use annealing_eq_f90=1 for fanning the Arrhenius model of Lasslett (1987)
!   use annealing_eq_f90=2 for fanning curvelinear model (as used by HeFTy)
!
!   compile this by running this command:
!   f2py -c calculate_reduced_AFT_lengths.f90 -m calculate_reduced_AFT_lengths
!   in the same folder as this module
!   f2py requires numpy to be installed, you can find this at www.numpy.org

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

    integer :: annealing_eq_f90

    real*8, intent(out) :: rcf(nsteps)
    
    integer segments, timestepcount
    
    real*8 gf(nsteps), teq(nsteps), rcmod(nsteps)
    real*8 tt, temp1, temp_term
    real*8 f1, f2

    do segments=1, nsteps
        gf(:) = 0
        teq(:) = 0
        do timestepcount=segments, nsteps
            tt = timestepduration(timestepcount)
            temp1 = temperature(timestepcount)

            if (annealing_eq_f90.eq.2) then
             temp_term = log(1.0/temp1)
            else
             temp_term = 1.0/temp1
            endif

            if (timestepcount.gt.segments) then
              teq(timestepcount) = exp( ((gf(timestepcount-1) - C0)/C1 &
                   * (temp_term-C3)) + C2 )
            else
                teq(timestepcount) = 0
            endif

            gf(timestepcount) = ( C0 + C1 * ((log(tt+teq(timestepcount))-C2) &
             / (temp_term-C3)) )

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

            
        
        
            
         
            
       
        
  
