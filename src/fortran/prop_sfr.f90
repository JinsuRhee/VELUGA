!234567 
      Subroutine prop_sfr(larr, darr, xx, yy, zz, age, mass, &
                xc, yc, zc)

      USE omp_lib
      IMPLICIT NONE

      INTEGER(KIND=4) larr(20)
      REAL(KIND=8) darr(20)
        
      REAL(KIND=8) xc, yc, zc
      REAL(KIND=8) xx(larr(2)), yy(larr(2)), zz(larr(2)), age(larr(2)), mass(larr(2))

      !!!!!!
      INTEGER(KIND=4) i, j, k, l
      INTEGER(KIND=4) n_gal, n_part, n_thread
      REAL(KIND=8) dum, dsize, dtime, sfr

      n_part    = larr(2)
      n_thread  = larr(3)
      dsize     = darr(1)
      dtime     = darr(11)
      call omp_set_num_threads(n_thread)

      sfr = 0.

      !$OMP PARALLEL DO default(shared) private(dum, j) schedule(static) reduction(+:sfr)
      DO i=1, n_part
        IF(age(i) .GT. dtime .OR. age(i) .LT. 0) CYCLE
        IF(mass(i) .LT. 0) CYCLE
        
        IF(dsize .LT. 0) THEN
               sfr = sfr + mass(i)
        ELSE
                dum = (xx(i) - xc)**2 + (yy(i) - yc)**2 + (zz(i) - zc)**2
                dum = sqrt(dum)
                IF(dum .LT. dsize) sfr = sfr + mass(i)
        ENDIF
      ENDDO

      sfr = sfr / dtime / 1.0d9
      darr(20) = sfr
      Return
      End subroutine
