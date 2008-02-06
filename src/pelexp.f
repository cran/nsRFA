      SUBROUTINE PELEXP(XMOM,PARA)
C***********************************************************************
C*                                                                     *
C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
C*                                                                     *
C*  J. R. M. HOSKING                                                   *
C*  IBM RESEARCH DIVISION                                              *
C*  T. J. WATSON RESEARCH CENTER                                       *
C*  YORKTOWN HEIGHTS                                                   *
C*  NEW YORK 10598, U.S.A.                                             *
C*                                                                     *
C*  VERSION 3     AUGUST 1996                                          *
C*                                                                     *
C***********************************************************************
C
C  PARAMETER ESTIMATION VIA L-MOMENTS FOR THE EXPONENTIAL DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  XMOM   * INPUT* ARRAY OF LENGTH 2. CONTAINS THE L-MOMENTS LAMBDA-1,
C                  LAMBDA-2.
C  PARA   *OUTPUT* ARRAY OF LENGTH 2. ON EXIT, CONTAINS THE PARAMETERS
C                  IN THE ORDER XI, ALPHA (LOCATION, SCALE).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XMOM(2),PARA(2)
      DATA ZERO/0D0/,TWO/2D0/
C
      IF(XMOM(2).LE.ZERO)GOTO 1000
      PARA(2)=TWO*XMOM(2)
      PARA(1)=XMOM(1)-PARA(2)
      RETURN
C
 1000 WRITE(6,7000)
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE PELEXP : L-MOMENTS INVALID')
      END
