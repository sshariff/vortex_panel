       PROGRAM MAIN
C
C  THIS MAIN PROGRAM BELONGS IN THE INITIALIZATION PART OF THE PANEL
C  METHOD CODE.
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  BNODEX(I,J)  =    NEW FISH NODE INDEX AND POSITIONS  (CM)
C  CM(I)        =    INSTANTANEOUS FISH CENTER OF MASS POSITION  (CM)
C  DT           =    INCREMENTAL TIME STEP  (S)
C  F(I)         =    INVISCID FLUID DYNAMIC FORCE ON FISH  (KG-CM/S^2)
C  IMAX         =    MAXIMUM NUMBER OF FISH NODES
C  M1           =    FIRST NODE INDEX TO BE OUTPUT
C  M2           =    LAST NODE INDEX TO BE OUTPUT
C  N            =    NUMBER OF TIME STEPPING LOOPS PERFORMED
C  OLDBNDX(I,J) =    OLD FISH NODE INDEX AND POSITIONS  (CM)
C  PER          =    TAIL BEAT PERIOD  (S)
C  THICK(I,J)   =    UNDEFORMED FISH NODE POSITIONS  (CM)
C  TIME         =    TIME OF THE PREVIOUS OR SOLVED TIME STEP  (S)
C  VEL(I)       =    INSTANTANEOUS FISH SWIMMING VELOCITIES  (CM/S)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
       IMPLICIT REAL *8 (A-H, O-Z)
C
       PARAMETER (IMAX=1400)
C
       COMMON   /MOTION/ TIME, DT, PER
C
C  ARRAYS FOR THE NODAL POSITION CALCULATIONS.
C
       DIMENSION F(3), CM(3), VEL(3), BNODEX(4,IMAX), THICK(3,IMAX),
     .           OLDBNDX(4,IMAX), COEF(7,14)
C
       DATA ZERO/0.0D+00/
C
C  OPEN THE DATA FILES USED TO STORE INPUTS AND OUTPUTS.
C
       OPEN (UNIT=13, STATUS='UNKNOWN', FILE='SUPER.DAT')
       OPEN (UNIT=14, STATUS='UNKNOWN', FILE='SHAPE.DAT')
       OPEN (UNIT=15, STATUS='UNKNOWN', FILE='MOTION.DAT')
C
C  INITIALIZE THE FLUID DYNAMIC FORCES.
C
       F(1) = ZERO
       F(2) = ZERO
       F(3) = ZERO
C
C  INITIALIZE THE CENTER OF MASS POSITION AT THE ORIGIN.  THE X-AXIS
C  POINTS TOWARDS THE TAIL, THE Y-AXIS POINTS UP, AND THE Z-AXIS POINTS
C  TO THE LEFT OF THE SWIMMING FISH.
C
       CM(1) = ZERO
       CM(2) = ZERO
       CM(3) = ZERO
C
C  INITIALIZE THE INSTANTANEOUS FORWARD SWIMMING VELOCITY IN CM/S AND
C  THE TAILBEAT PERIOD IN S.  IT'S BETTER TO START WITH THE KNOWN
C  SWIMMING VELOCITY AND SEE WHERE WE END UP THEN TO ASK THE FISH TO
C  SHAKE AROUND TRYING TO ATTAIN TERMINAL VELOCITY.
C
       VEL(1) = -24.93D+00
       VEL(2) = ZERO
       VEL(3) = ZERO
       PER = 0.233D+00
C
C  READ IN THE FISH THICKNESS FILE THAT IS ALSO USED FOR BUILDING THE
C  CONNECTIVITY MATRIX.  THE CURRENT FILE SUPER.DAT DOES NOT HAVE A
C  KUTTA INDEX, CAUDAL INDEX, STRIP INDEX, OR NUMBER OF NODES PER STRIP.
C  THE ARRAY THICK(I,J) SHOULD BE SET UP FROM SUPER.DAT AT THE SAME TIME
C  AND IN THE SAME ORDER AS THE NODES IN BNODEX(I,J).  THE METHOD USED
C  HERE TO ENTER THE VALUES OF THICK(I,J) FOR J>697 IS BOGUS.
C  THICK(I,J) CONTAINS (X,Y,Z) NODE POSITIONS FOR A STRAIGHT SPINE.
C
       DO I=1,697
          READ(13,*) THICK(1,I), THICK(2,I), THICK(3,I)
       END DO
       DO I=1,697
          THICK(1,I+697) = THICK(1,I)
          THICK(2,I+697) = THICK(2,I)
          THICK(3,I+697) = - THICK(3,I)
       END DO
C
C  READ IN THE POLYNOMIAL MOTION COEFFICIENTS.  THE INDEX I=1 IS
C  THE NONDIMENSIONAL TIME AT WHICH THE SHAPE IS VALID.  THE NEXT SIX
C  VALUES I=2-7 ARE THE POLYNOMIAL COEFFICIENTS OF FISH DEFLECTION
C  IN THE Z-DIRECTION Z=Z(X), STARTING WITH THE COEFFICIENT OF X
C  WHEN I=2 AND ENDING WITH X^6 WHEN I=7.  THE SECOND INDEX J=1-14
C  CORRESPONDS TO THE NUMBER OF THE CURVE FIT.  BOTH X AND Z ARE IN
C  CM.
C
       DO J=1,14
          READ(15,*) COEF(1,J), COEF(2,J), COEF(3,J), COEF(4,J),
     .               COEF(5,J), COEF(6,J), COEF(7,J)
       END DO
C
C  INITIALIZE THE ARRAYS BNODEX AND OLDBNDX.  THE ZEROES MAY CAUSE
C  PROBLEMS WHEN CALCULATNG VELOCITIES DURING THE FIRST TIME STEP.
C  THE NODE NUMBER IS GIVEN TO AVOID HAVING TO COPY IT LATER ON.  ONE
C  OPTION IS TO MAKE VELOCITIES ZERO IF THE VALUES OF OLDBNDX ARE ALSO
C  ZERO.
C
       DO I=1,IMAX
          BNODEX(1,I) = DFLOAT(I)
          BNODEX(2,I) = ZERO
          BNODEX(3,I) = ZERO
          BNODEX(4,I) = ZERO
          OLDBNDX(1,I) = DFLOAT(I)
          OLDBNDX(2,I) = ZERO
          OLDBNDX(3,I) = ZERO
          OLDBNDX(4,I) = ZERO
       END DO
C
C  PROVIDE A TIME, TIME STEP, NUMBER OF LOOPS, AND OUTPUT NODES.  THIS
C  WILL BE ELIMINATED IN YOUR OWN PROGRAM.
C
       WRITE (*,*) 'ENTER THE INITIAL TIME (0.0)'
       READ (*,*) TIME
       WRITE (*,*)
       WRITE (*,*) 'ENTER THE TIME STEP (0.008)'
       READ (*,*) DT
       WRITE (*,*)
       WRITE (*,*) 'ENTER THE NUMBER OF LOOPS'
       READ (*,*) N
       WRITE (*,*)
       WRITE (*,*) 'ENTER THE FIRST NODE TO OUTPUT'
       READ (*,*) M1
       WRITE (*,*)
       WRITE (*,*) 'ENTER THE LAST NODE TO OUTPUT'
       READ (*,*) M2
C
C  IMPOSE AN IMPORTANT CONDITION ON THE MAXIMUM TIME STEP TO ENSURE THAT
C  THE INTERPOLATION ALGORITHM IN SUBROUTINE BMOTION WORKS.  THIS LIMITS
C  THE COURANT NUMBER BUT IS A NECESSARY EVIL TO RESOLVE SWIMMING
C  MOTIONS.  KEEP THIS!
C
       IF (DT.GT.(PER / 1.4D+01)) THEN
          DT = PER / 3.0D+01
       END IF
C
C  DO A NUMBER OF LOOPS TO SEE THE FISH SWIM.  ELIMINATE THIS.
C
       DO I=1,N
C
C  CALL THE SUBROUTINE BMOTION TO BEGIN MOVEMENT CALCULATIONS.
C
          CALL BMOTION(THICK,BNODEX,OLDBNDX,F,CM,VEL,COEF)
C
C  UPDATE THE TIME.
C
          TIME = TIME + DT
C
C  WRITE THE FISH SHAPE INTO THE FILE SHAPE.DAT.  ELIMINATE THIS.
C
          DO J=M1,M2
             WRITE (14,*) BNODEX(2,J), BNODEX(3,J), BNODEX(4,J)
          END DO
       END DO
C
C  CLOSE THE DATA FILES BEFORE STOPPING THE PROGRAM.
C
       CLOSE (UNIT=13, STATUS='KEEP')
       CLOSE (UNIT=14, STATUS='KEEP')
       CLOSE (UNIT=15, STATUS='KEEP')
C
       STOP
C
       END
C
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
       SUBROUTINE BMOTION(THICK,BNODEX,OLDBNDX,F,CM,VEL,COEF)
C
C
C  THIS SUBROUTINE UPDATES THE NODE POSITIONS AND REPOSITIONS THE FISH
C  CENTER OF MASS.
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  A            =    INTERPOLATION QUADRATIC COEFFICIENT  (CM)
C  AREA         =    AREA OF ONE FISH FLANK  (CM^2)
C  B            =    INTERPOLATION QUADRATIC COEFFICIENT  (CM/S)
C  BNODEX(I,J)  =    NEW FISH NODE INDEX AND POSITIONS  (CM)
C  C            =    INTERPOLATION QUADRATIC COEFFICIENT  (CM/S^2)
C  CM(I)        =    INSTANTANEOUS FISH CENTER OF MASS POSITION  (CM)
C  COEF(I,J)    =    TIMES AND COEFFICIENTS (CM) FOR 14 CURVE FITS
C  DER          =    INTERPOLATED FINAL SPINE SLOPE FOR TRANSFORMATION
C  DER1         =    FINAL SPINE SLOPE FOR CURVE FIT N1
C  DER2         =    FINAL SPINE SLOPE FOR CURVE FIT N2
C  DER3         =    FINAL SPINE SLOPE FOR CURVE FIT N3
C  DRAG         =    ESTIMATE OF TOTAL FRICTIONAL FORCE  (KG-CM/S^2)
C  DT           =    INCREMENTAL TIME STEP  (S)
C  DT1          =    DIFFERENCE BETWEEN TIME AT N2 AND N1  (S)
C  DT3          =    DIFFERENCE BETWEEN TIME AT N3 AND N2  (S)
C  DX           =    ARCLENGTH INTEGRATION STEP  (CM)
C  DZDX         =    LOCAL SPINE SLOPE DURING ARCLENGTH INTEGRATION
C  F(I)         =    INVISCID FLUID DYNAMIC FORCE ON FISH  (KG-CM/S^2)
C  FMAS         =    FISH MASS  (KG)
C  FRAC1        =    NONDIMENSIONAL FRACTION OF PERIOD AT TIME
C  FRAC2        =    NONDIMENSIONAL FRACTION OF PERIOD AT TIME+DT
C  IMAX         =    MAXIMUM NUMBER OF FISH NODES
C  N            =    NUMBER OF TIME STEPPING LOOPS PERFORMED
C  N1           =    FIRST CURVE FIT INDEX
C  N2           =    SECOND CURVE FIT INDEX
C  N3           =    THIRD CURVE FIT INDEX
C  OLDBNDX(I,J) =    OLD FISH NODE INDEX AND POSITIONS  (CM)
C  PER          =    TAIL BEAT PERIOD  (S)
C  PI           =    THAT TRANSCENDENTAL NUMBER 3.14592 ...
C  QUCK         =    SPEEDS UP RUN TIME BY CALLING DSQRT() ONLY ONCE
C  RHO          =    DENSITY OF FRESH WATER  (KG/CM^3)
C  S1           =    ARCLENGTH ALONG THE CURVE FIT N1  (CM)
C  S2           =    ARCLENGTH ALONG THE CURVE FIT N2  (CM)
C  S3           =    ARCLENGTH ALONG THE CURVE FIT N3  (CM)
C  TEST         =    ARCLENGTH TO TEST IS ANOTHER STEP IS NEEDED (CM)
C  THICK(I,J)   =    UNDEFORMED FISH NODE POSITIONS  (CM)
C  TIME         =    TIME OF THE PREVIOUS OR SOLVED TIME STEP  (S)
C  UMAX         =    EXPECTED FISH SWIMMING VELOCITY  (CM/S)
C  VEL(I)       =    INSTANTANEOUS FISH SWIMMING VELOCITIES  (CM/S)
C  X            =    CURRENT ARCLEGNTH INTEGRATION POSITION  (CM)
C  X1           =    X POSITION OF FINAL POINT ALONG CURVE FIT N1  (CM)
C  X2           =    X POSITION OF FINAL POINT ALONG CURVE FIT N2  (CM)
C  X3           =    X POSITION OF FINAL POINT ALONG CURVE FIT N3  (CM)
C  XCM          =    INTERPOLATED CENTER OF MASS X POSITION  (CM)
C  XCOM         =    DISTANCE FROM NOSE TO CENTER OF MASS  (CM)
C  XNOT         =    X VALUE OF THE CURRENT FISH SLICE  (CM)
C  XO           =    INTERPOLATED X POSITION FOR TRANSFORMATION  (CM)
C  XTOT         =    DISTANCE FROM NOSE TO TAIL  (CM)
C  Z1           =    Z POSITION OF FINAL POINT ALONG CURVE FIT N1  (CM)
C  Z2           =    Z POSITION OF FINAL POINT ALONG CURVE FIT N2  (CM)
C  Z3           =    Z POSITION OF FINAL POINT ALONG CURVE FIT N3  (CM)
C  ZCM          =    INTERPOLATED CENTER OF MASS Z POSITION  (CM)
C  ZO           =    INTERPOLATED Z POSITION FOR TRANSFORMATION  (CM)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
       IMPLICIT REAL *8 (A-H, O-Z)
C
       PARAMETER (IMAX=1400)
C
       COMMON   /MOTION/ TIME, DT, PER
C
C  ARRAYS FOR THE NODAL POSITION CALCULATIONS.
C
       DIMENSION F(3), CM(3), VEL(3), BNODEX(4,IMAX), THICK(3,IMAX),
     .           OLDBNDX(4,IMAX), COEF(7,14)
C
       DATA HALF/0.5D+00/, TWO/2.0D+00/, ZERO/0.0D+00/, TRE/3.0D+00/,
     .     FOR/4.0D+00/, FIV/5.0D+00/, SIX/6.0D+00/, ONE/1.0D+00/
C
C  CALCULATE THE VALUE OF PI.
C
       PI = FOR * DATAN(ONE)
C
C  COPY THE CURRENT VALUES OF BNODEX INTO OLDBNDX SO THAT THEY BECOME
C  THE PREVIOUS VALUES OF BNODEX.  THERE IS NO APPARENT NEED TO COPY
C  THE NODE INDEX AS WELL, AT LEAST FOR THE FISH BODY.  RECALL THAT
C  THE FLUID DYNAMICS SOLUTION IS FOR TIME + DT.
C
       DO I=2,4
          DO J=1,IMAX
             OLDBNDX(I,J) = BNODEX(I,J)
          END DO
       END DO
C
C  STORE THE FISH TOTAL LENGTH IN CM AND SWIMMING SPEED IN CM/S.  THE
C  FISH MASS IS IN KG.  THE AREA OF ONE SIDE IS IN CM^2.  THE DENSITY
C  OF WATER IS IN KG/CM^3.  THE DISTANCE OF THE CENTER OF MASS FROM THE
C  FISH NOSE IS XCOM GIVEN IN CM.  NOTE THAT THE FISH "CHORD" IS NO
C  LONGER UNITY.  THE OLD VERSION OF THE PROGRAM ASSUMED C=1 IN PLACES.
C
       XTOT = 12.02D+00
       UMAX = 24.93D+00
       FMAS = 0.0298D+00
       AREA = 42.3D+00
       XCOM = 4.476D+00
       RHO = 1.0D-03
C
C  ESTIMATE THE DRAG FORCE ACTING ON THE FISH BY ASSUMING THAT THE
C  FISH IS ON A FLAT PLATE -- BAD PUN.  THE FACTOR OF TWO IS NEEDED
C  BECAUSE THE FISH HAS TWO SIDES.  THE FRICTION COEFFICIENT 0.005 IS
C  PROBABLY WRONG BECASUE 1) IT IS A GUESS AND 2) THESE ARE UNSTEADY
C  BOUNDARY LAYERS.
C
       DRAG = TWO * 0.005 * AREA * RHO * VEL(1)**2 / TWO
C
C  GIVE SOME FAKE FISH BODY FORCES FROM THE FLUID DYNAMICS JUST TO
C  MAKE THIS SECTION COMPLETE.  THESE FORCES WILL HAVE TO HAVE UNITS
C  OF KG-CM/S^2 IN ORDER TO MATCH THE CURRENT FISH UNITS.  THESE WILL
C  BE ELIMINATED IN THE FUTURE AND REPLACED WITH REAL FORCES.
C
       F(1) = -TWO * 0.005 * AREA * RHO * UMAX*2 / TWO
       F(2) = ZERO
       F(3) = ZERO
C
C  CALCULATE CENTER OF MASS MOTION BASED ON A FORCE BALANCE TAKEN OVER
C  THE ENTIRE FISH BODY.  THIS WILL EVENTUALLY BE ELIMINATED SINCE ROB
C  WILL SEND US POSITION DATA WITH THE CENTER OF MASS ALREADY MOVED.
C  A SECOND ORDER TAYLOR SERIES EXPANSION OF POSITION IS USED HERE.
C
       CM(1) = CM(1) + VEL(1) * DT + HALF * DT**2 * (F(1) + DRAG) /
     .         FMAS
       CM(2) = CM(2) + VEL(2) * DT + HALF * DT**2 * F(2) / FMAS
       CM(3) = CM(3) + VEL(3) * DT + HALF * DT**2 * F(3) / FMAS
C
C  UPDATE THE SWIMMING VELOCITY IN CM/S.  THIS MUST BE KEPT IN THE
C  FUTURE TO CALCULATE DRAG FORCE.  WHAT DO WE DO IF THERE IS A NET
C  VERTICAL FORCE ACTING ON THE FISH?
C
       VEL(1) = VEL(1) + DT * (F(1) + DRAG) / FMAS
       VEL(2) = VEL(2) + DT * F(2) / FMAS
       VEL(3) = VEL(3) + DT * F(3) / FMAS
C
C  NOW FIGURE OUT WHERE THE CURRENT TIME AND NEXT TIME ARE LOCATED
C  WITHIN THE 14 VALUES OF COEF(1,J).  THE FRACTION OF THE TAILBEAT
C  PERIOD BETWEEN ZERO AND ONE IS NEEDED HERE.
C
       FRAC1 = DMOD((TIME / PER),ONE)
       FRAC2 = DMOD(((TIME + DT) / PER),ONE)
C
C  SOME CARE IS NEEDED WHEN ASSIGNING THE THREE INTEGER VALUES OF THE
C  SPINE SHAPES THAT WILL BE USED FOR INTERPOLATION SINCE THE FISH TAIL
C  DOES NOT UNDERGO A SMOOTH MOTION.  THE INTEGERS N1 AND N2 INDICATE
C  THE VALUES OF COEF(1,J) THAT PRECEED FRAC1 AND FRAC2, RESPECIVELY.
C
       N1 = 0
       N2 = 0
       DO I=1,13
          IF ((FRAC1.GE.COEF(1,I)).AND.(FRAC1.LT.COEF(1,I+1))) THEN
             N1 = I
          ELSE IF (FRAC1.GE.COEF(1,14)) THEN
             N1 = 14
          END IF
          IF ((FRAC2.GE.COEF(1,I)).AND.(FRAC2.LT.COEF(1,I+1))) THEN
             N2 = I
          ELSE IF (FRAC2.GE.COEF(1,14)) THEN
             N2 = 14
          END IF
       END DO
       IF ((N1.EQ.0).OR.(N2.EQ.0)) THEN
          WRITE (*,*)
          WRITE (*,*) 'N1 AND/OR N2 IS ZERO'
          WRITE (*,*) 'THE PROGRAM HAS QUIT, PRESS RETURN'
          PAUSE
          STOP
       END IF
C
C  CHECK IF N1 OR N2 ARE NEAR THE TRANSITION 14-1-2 OR 7-8-9 OF TAIL
C  MOTION FROM ONE DIRECTION TO THE OTHER.  INTEGERS 1 AND 8 REPRESENT
C  EXTREME LATERAL MOTION.  HENCE, N2 CANNOT EQUAL 1 OR 8 SINCE THAT
C  WOULD FORCE INTERPOLATION TO STRADDLE A SHARP TRANSITION.  MAKE SURE
C  THAT N3 IS NOT GREATER THAN 14.
C
       IF ((N1.EQ.7).OR.(N2.EQ.7)) THEN
          IF ((N1.EQ.7).AND.(N2.EQ.8)) THEN
             N1 = 8
             N2 = 9
             N3 = 10
          ELSE
             N1 = 6
             N2 = 7
             N3 = 8
          END IF
       ELSE IF ((N1.EQ.14).OR.(N2.EQ.14)) THEN
          IF ((N1.EQ.14).AND.(N2.EQ.1)) THEN
             N1 = 1
             N2 = 2
             N3 = 3
          ELSE
             N1 = 13
             N2 = 14
             N3 = 1
          END IF
       ELSE
          N2 = N1 + 1
          N3 = N1 + 2
       END IF
       IF (N3.GT.14) THEN
          N3 = N3 - 14
       END IF
C
C  THE VALUES OF COEF(I,J) ARE USED IN A WRAP AROUND FASHION.  HENCE,
C  THE VALUE I=1 NEEDS TO BE BOTH ZERO AND ONE, DEPENDING ON THE VALUE
C  OF N3.  IF N2=14 AND N3=1, THEN COEF(1,1) IS AT THE END OF A TAILBEAT
C  AND NEEDS TO INDICATE A NONDIMENSIONAL TIME OF ONE.   OTHERWISE, ITS
C  VALUE IS ZERO.
C
       IF (N3.EQ.1) THEN
          COEF(1,1) = ONE
       END IF
C
C  BEGIN THE LOOP THAT TRANSLATES AND ROTATES THE NODE POSITIONS FROM
C  THE ARRAY THICK INTO THE ARRAY BNODEX FOR THE NEXT TIME STEP AT
C  TIME+DT.  THE VARIABLE XNOT KEEPS THE PROGRAM FROM RECREATING THE
C  COORDINATE TRANSFORMATION DATA IF NODES SHARE THE SAME VALUE OF
C  THICK(1,I).  THIS WILL SAVE A LOT OF TIME.
C
       XNOT = -1.0D10
       DO I=1,1394
          IF (THICK(1,I).NE.XNOT) THEN
             XNOT = THICK(1,I)
C
C  CALCULATE THE POSITION (X,Z) WHERE THE ARCLENGTH S=XNOT ALONG THE
C  THREE CURVE FITS IDENTIFIED BY N1, N2, N3.  THE THREE ARCLENGTHS S1,
C  S2, S3 CORRESPOND TO N1, N2, N3, RESPECTIVELY.  THE ARCLENGTHS ARE
C  INTEGRATED UNTIL S=XNOT, WHICH HAS A DISTINCT POSITION (X,Z) NOT
C  EQUAL TO S OR XNOT AND THAT MUST BE FOUND FROM THE CURVE FIT.
C  ACTUAL FISH SHAPE WILL BE FOUND BY SECOND ORDER INTERPOLATION
C  BETWEEN THREE PAIRS OF COORDINATES (X1,Z1), (X2,Z2), (X3,Z3).  THE
C  DERIVATIVE DZDX IS NEEDED TO FIND THE ORIENTATION OF THE VECTOR
C  NORMAL TO THE SPINE WHEN ROTATING FISH BODY NODES.  SECOND ORDER
C  INTEGRATION IS PERFORMED AND THE LAST STEP DX IS ADJUSTED TO FALL
C  PRECISELY ON XNOT.  ACCURACY IS IMPORTANT SINCE THESE POSITIONS WILL
C  DETERMINE NODE VELOCITIES.  CALLING THE DSQRT() ALGORITHM ONLY ONCE
C  WILL SAVE COMPUTATIONAL TIME.
C
             S1 = ZERO
             S2 = ZERO
             S3 = ZERO
             DX = DABS(XTOT / 5.0D+03)
             X = -DX / TWO
100          CONTINUE
             X = X + DX
             DZDX = COEF(2,N1) + TWO * COEF(3,N1) * X + TRE *
     .              COEF(4,N1) * X**2 + FOR * COEF(5,N1) * X**3 +
     .              FIV * COEF(6,N1) * X**4 + SIX * COEF(7,N1) * X**5
             QUCK = DSQRT(ONE + DZDX**2)
             S1 = S1 + QUCK * DX
             TEST = S1 + QUCK * DX
             IF (TEST.LT.XNOT) GO TO 100
             X = X + DX / TWO + (XNOT - S1) / QUCK
             X1 = X
             Z1 = COEF(2,N1) * X + COEF(3,N1) * X**2 + COEF(4,N1) *
     .            X**3 + COEF(5,N1) * X**4 + COEF(6,N1) * X**5 +
     .            COEF(7,N1) * X**6
             DER1 = COEF(2,N1) + TWO * COEF(3,N1) * X + TRE *
     .              COEF(4,N1) * X**2 + FOR * COEF(5,N1) * X**3 +
     .              FIV * COEF(6,N1) * X**4 + SIX * COEF(7,N1) * X**5
             X = -DX / TWO
200          CONTINUE
             X = X + DX
             DZDX = COEF(2,N2) + TWO * COEF(3,N2) * X + TRE *
     .              COEF(4,N2) * X**2 + FOR * COEF(5,N2) * X**3 +
     .              FIV * COEF(6,N2) * X**4 + SIX * COEF(7,N2) * X**5
             QUCK = DSQRT(ONE + DZDX**2)
             S2 = S2 + QUCK * DX
             TEST = S2 + QUCK * DX
             IF (TEST.LT.XNOT) GO TO 200
             X = X + DX / TWO + (XNOT - S2) / QUCK
             X2 = X
             Z2 = COEF(2,N2) * X + COEF(3,N2) * X**2 + COEF(4,N2) *
     .            X**3 + COEF(5,N2) * X**4 + COEF(6,N2) * X**5 +
     .            COEF(7,N2) * X**6
             DER2 = COEF(2,N2) + TWO * COEF(3,N2) * X + TRE *
     .              COEF(4,N2) * X**2 + FOR * COEF(5,N2) * X**3 +
     .              FIV * COEF(6,N2) * X**4 + SIX * COEF(7,N2) * X**5
             X = -DX / TWO
300          CONTINUE
             X = X + DX
             DZDX = COEF(2,N3) + TWO * COEF(3,N3) * X + TRE *
     .              COEF(4,N3) * X**2 + FOR * COEF(5,N3) * X**3 +
     .              FIV * COEF(6,N3) * X**4 + SIX * COEF(7,N3) * X**5
             QUCK = DSQRT(ONE + DZDX**2)
             S3 = S3 + QUCK * DX
             TEST = S3 + QUCK * DX
             IF (TEST.LT.XNOT) GO TO 300
             X = X + DX / TWO + (XNOT - S3) / QUCK
             X3 = X
             Z3 = COEF(2,N3) * X + COEF(3,N3) * X**2 + COEF(4,N3) *
     .            X**3 + COEF(5,N3) * X**4 + COEF(6,N3) * X**5 +
     .            COEF(7,N3) * X**6
             DER3 = COEF(2,N3) + TWO * COEF(3,N3) * X + TRE *
     .              COEF(4,N3) * X**2 + FOR * COEF(5,N3) * X**3 +
     .              FIV * COEF(6,N3) * X**4 + SIX * COEF(7,N3) * X**5
C
C  INTERPOLATE TO FIND SECOND ORDER ACCURATE X, Z, AND DZDX VALUES FOR
C  THIS SUNFISH SLICE.  THE COEFFICIENTS OF A QUADRATIC CURVE THROUGH
C  THE THREE KNOWN VALUES ARE CALCULATED.  RECALL THAT THE TIMES IN
C  THE ARRAY COEF(1,I) ARE MADE NONDIMENSIONAL BY THE PERIOD PER AND
C  ARE BETWEEN ZERO AND ONE.  HENCE, THE CORRECT TIME TO USE IN THE
C  INTERPOLATION IS THE NONDIMENSIONAL TIME FRAC2.
C
C
             DT1 = COEF(1,N2) - COEF(1,N1)
	     DT3 = COEF(1,N3) - COEF(1,N2)
	     C = ((X1 - X2) / DT1 + (X3 - X2) / DT3) / (DT1 + DT3)
	     B = C * DT1 - (X1 - X2) / DT1 - TWO * C * COEF(1,N2)
	     A = X2 - B * COEF(1,N2) - C * COEF(1,N2)**2
	     XO = A + B * FRAC2 + C * FRAC2**2
	     C = ((Z1 - Z2) / DT1 + (Z3 - Z2) / DT3) / (DT1 + DT3)
	     B = C * DT1 - (Z1 - Z2) / DT1 - TWO * C * COEF(1,N2)
	     A = Z2 - B * COEF(1,N2) - C * COEF(1,N2)**2
	     ZO = A + B * FRAC2 + C * FRAC2**2
	     C = ((DER1 - DER2) / DT1 + (DER3 - DER2) / DT3) /
     .            (DT1 + DT3)
             B = C * DT1 - (DER1 - DER2) / DT1 - TWO * C * COEF(1,N2)
             A = DER2 - B * COEF(1,N2) - C * COEF(1,N2)**2
             DER = A + B * FRAC2 + C * FRAC2**2
C
C  THIS COMPLETES THE IF STATEMENT THAT PROVIDES VALUES OF XO, ZO,
C  AND DER FOR COORDINATE TRANSFORMATIONS.
C
          END IF
C
C  APPLY A TRANSLATION AND ROTATION TO PLACE THE NODES IN THE RIGHT
C  PLACES.  THICK(1,I) HAS ALREADY BEEN USED TO PROVIDE (XO,ZO).  THE
C  QUANTITY DER / DSQRT(ONE + DER**2) REPRESENTS SINE OF THE SPINE
C  ANGLE WHILE 1 / DSQRT(ONE + DER**2) IS ITS COSINE.  THIS ANGLE
C  ROTATES THE CONTRIBUTION OF THICK(3,I) TO A NODE POSITION.  THE
C  NEGATIVE SIGN IS A RESULT OF AXES DEFINITIONS.  THE VERTICAL NODE
C  POSITION IS ONLY MODIFIED BY CENTER OF MASS MOTION.
C
          BNODEX(2,I) = XO - THICK(3,I) * DER / DSQRT(ONE + DER**2)
          BNODEX(3,I) = THICK(2,I)
          BNODEX(4,I) = ZO + THICK(3,I) / DSQRT(ONE + DER**2)
C
C  THIS COMPLETES THE DO LOOP OVER ALL NODE VALUES IN THICK(I,J).
C
       END DO
C
C  FIND THE FISH CENTER OF MASS POSITION ACCORDING TO THE CURVE FITS.
C  THESE THEN NEED TO BE SUBTRACTED OFF TO PLACE THE CENTER OF MASS
C  WHERE IT ACTUALLY BELONGS, THAT IS AT CM(I).
C
       XNOT = XCOM
       S1 = ZERO
       S2 = ZERO
       S3 = ZERO
       X = -DX / TWO
400    CONTINUE
       X = X + DX
       DZDX = COEF(2,N1) + TWO * COEF(3,N1) * X + TRE *
     .        COEF(4,N1) * X**2 + FOR * COEF(5,N1) * X**3 +
     .        FIV * COEF(6,N1) * X**4 + SIX * COEF(7,N1) * X**5
       QUCK = DSQRT(ONE + DZDX**2)
       S1 = S1 + QUCK * DX
       TEST = S1 + QUCK * DX
       IF (TEST.LT.XNOT) GO TO 400
       X = X + DX / TWO + (XNOT - S1) / QUCK
       X1 = X
       Z1 = COEF(2,N1) * X + COEF(3,N1) * X**2 + COEF(4,N1) *
     .      X**3 + COEF(5,N1) * X**4 + COEF(6,N1) * X**5 +
     .      COEF(7,N1) * X**6
       X = -DX / TWO
500    CONTINUE
       X = X + DX
       DZDX = COEF(2,N2) + TWO * COEF(3,N2) * X + TRE *
     .        COEF(4,N2) * X**2 + FOR * COEF(5,N2) * X**3 +
     .        FIV * COEF(6,N2) * X**4 + SIX * COEF(7,N2) * X**5
       QUCK = DSQRT(ONE + DZDX**2)
       S2 = S2 + QUCK * DX
       TEST = S2 + QUCK * DX
       IF (TEST.LT.XNOT) GO TO 500
       X = X + DX / TWO + (XNOT - S2) / QUCK
       X2 = X
       Z2 = COEF(2,N2) * X + COEF(3,N2) * X**2 + COEF(4,N2) *
     .      X**3 + COEF(5,N2) * X**4 + COEF(6,N2) * X**5 +
     .      COEF(7,N2) * X**6
       X = -DX / TWO
600    CONTINUE
       X = X + DX
       DZDX = COEF(2,N3) + TWO * COEF(3,N3) * X + TRE *
     .        COEF(4,N3) * X**2 + FOR * COEF(5,N3) * X**3 +
     .        FIV * COEF(6,N3) * X**4 + SIX * COEF(7,N3) * X**5
       QUCK = DSQRT(ONE + DZDX**2)
       S3 = S3 + QUCK * DX
       TEST = S3 + QUCK * DX
       IF (TEST.LT.XNOT) GO TO 600
       X = X + DX / TWO + (XNOT - S3) / QUCK
       X3 = X
       Z3 = COEF(2,N3) * X + COEF(3,N3) * X**2 + COEF(4,N3) *
     .      X**3 + COEF(5,N3) * X**4 + COEF(6,N3) * X**5 +
     .      COEF(7,N3) * X**6
C
C  INTERPOLATE TO FIND AN ACCURATE CENTER OF MASS POSITION.  THE VALUE
C  OF YCM IS ASSUMED TO BE ZERO.
C
       DT1 = COEF(1,N2) - COEF(1,N1)
       DT3 = COEF(1,N3) - COEF(1,N2)
       C = ((X1 - X2) / DT1 + (X3 - X2) / DT3) / (DT1 + DT3)
       B = C * DT1 - (X1 - X2) / DT1 - TWO * C * COEF(1,N2)
       A = X2 - B * COEF(1,N2) - C * COEF(1,N2)**2
       XCM = A + B * FRAC2 + C * FRAC2**2
       C = ((Z1 - Z2) / DT1 + (Z3 - Z2) / DT3) / (DT1 + DT3)
       B = C * DT1 - (Z1 - Z2) / DT1 - TWO * C * COEF(1,N2)
       A = Z2 - B * COEF(1,N2) - C * COEF(1,N2)**2
       ZCM = A + B * FRAC2 + C * FRAC2**2
C
C  MODIFY THE NODE POSITIONS SO THAT THE CENTER OF MASS STARTS AT THE
C  ORIGIN AT TIME=0 AND THEN MOVES WITH CM(I).
C
       DO I=1,1394
          BNODEX(2,I) = BNODEX(2,I) - XCOM - XCM + CM(1)
          BNODEX(3,I) = BNODEX(3,I) + CM(2)
          BNODEX(4,I) = BNODEX(4,I) - ZCM + CM(3)
       END DO
C
C  ENSURE THAT THE ORIGINAL VALUE IS RESTORED UNDER ALL CONDITIONS.
C
       COEF(1,1) = ZERO
C
       RETURN
C
       END
C
C

