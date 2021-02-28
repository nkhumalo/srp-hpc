c-----------------------------------------------------------------------
C>
C> \brief Compute b-coordinate from damping factor
C>
C>
      double precision function wfn1_func_b(d)
      implicit none
C>    Damping factor
      double precision, intent(in) :: d
c
      wfn1_func_b=sin(0.5d0*atan(2.0d0*d))**2
c
      return
      end
c
c-----------------------------------------------------------------------
C>
C> \brief Compute damping factor from b-coordinate
C>
C>
      double precision function wfn1_func_d(b)
      implicit none
C>    Damping factor
      double precision, intent(in) :: b
c
      wfn1_func_d=0.5d0*tan(2.0d0*asin(sqrt(b)))
c
      return
      end
c
c-----------------------------------------------------------------------
C>
C> \brief Compute gradient scale factor
C>
C>
      double precision function wfn1_func_dd(d)
      implicit none
C>    Damping factor
      double precision, intent(in) :: d
      double precision :: t
c
      t=0.5d0*atan(d)
      wfn1_func_dd=-(d*d+1.0d0)/(cos(t)*sin(t))
c
      return
      end
c
c-----------------------------------------------------------------------
