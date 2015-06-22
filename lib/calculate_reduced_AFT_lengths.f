
C FILE: AFT.F

      SUBROUTINE REDUCED_LN(rmf, rcf, timestepDuration,temperature,
     $   rmr0, kappa, Nsteps)
C
C     fission track annealing algorithm, using Ketcham 2005, 2007 eqs.
C     calculate reduced track lengths
C   
      IMPLICIT NONE

      INTEGER Nsteps, i
      INTEGER segments, timestepcount
      
      REAL*8 temperature(Nsteps)
      REAL*8 alfa, C0, C1, C2, C3
      REAL*8 rmr0, kappa
      REAL*8 tt, temp1
      REAL*8 f1, f2
      REAL*8 rc_a, rc_b, rc_c
      
      REAL*8 gf(Nsteps), teq(Nsteps), rcmod(Nsteps)
      REAL*8 rcf(Nsteps), rmf(Nsteps)
      REAL*8 timestepDuration(Nsteps)


Cf2py intent(in) time, temperature, rmr0, kappa, Nsteps

Cf2py intent(out) rmf, rcf


      !! Nsteps = Ntime - 1
      

      
       !! g function parameters:
       alfa =  0.04672
       !! fanning curvelinear function
       C0 =  0.39528
       C1 =  0.01073
       C2 =  -65.12969
       C3 =  -7.91715
    
      
      DO 101 segments=1, Nsteps
       gf(:) = 0
       teq(:) = 0
       DO 102 timestepcount=segments, Nsteps
        tt = timestepDuration(timestepcount)
        temp1 = temperature(timestepcount)
        
        if (timestepcount.gt.segments) then
         teq(timestepcount) = exp( ((gf(timestepcount-1) - C0)/C1 
     $    * (log(1.0/temp1)-C3)) + C2 )
        else
         teq(timestepcount) = 0
        endif
        
        gf(timestepcount) = ( C0 + C1 * 
     $   ((log(tt+teq(timestepcount))-C2) / 
     $   (log(1./temp1)-C3)) )

  102 continue    
        
C       f1 = math.pow(g[-1], (1./alfa))
        f1 = (gf(Nsteps)**(1.0/alfa))
        f2 = f1+1
        rcmod(segments) = 1.0 / f2
        
        rcf(segments) =  ((rcmod(segments) -rmr0)/
     $   (1.0-rmr0)) ** kappa
        
        if (rcmod(segments).lt.rmr0) then
            rcf(segments) = 0
        endif

  101 continue
        
        
C project c-axis corrected lenghts to normal lengths:            
        rc_a = -1.499
        rc_b = 4.150
        rc_c = -1.656
        rmf = rc_a * (rcf**2) + 
     $   rc_b * rcf + rc_c
    
        DO 103, i=1, Nsteps 
         if (rmf(i).lt.0.0) then
         rmf(i) = 0.0
         endif
  103 continue
        !write (*,*) rmf
        
        
       end

C END FILE AFT.F

            
        
        
            
         
            
       
        
  
