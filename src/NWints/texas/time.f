      SUBROUTINE GVBMP3
c $Id: time.f,v 1.3 1996-02-09 01:15:47 d3g681 Exp $
C     EMPTY
      END
C
      SUBROUTINE txs_SECOND(T)
      double precision T
c     REAL*8 T
      double precision util_cpusec
      external util_cpusec
C      DIMENSION ITIME(4) 
C     I1=TIMES(ITIME)
C     WRITE(*,*) ITIME
C     T=DBLE(ITIME(1)+ITIME(2)+ITIME(3)+ITIME(4))/100.0D0         
C      T=MCLOCK()/100.0D0
      T = util_cpusec()
C  
      RETURN
      END
C
*rak:      SUBROUTINE FDATE(X)
*rak:      CHARACTER*24 BLANK,X
*rak:      DATA BLANK /'                        '/
*rak:      X=BLANK
*rak:      END
