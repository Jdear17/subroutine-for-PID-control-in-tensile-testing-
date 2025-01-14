        SUBROUTINE UAMP(
     *     ampName, time, ampValueOld, dt, nProps, props, nSvars, 
     *     svars, lFlagsInfo,
     *     nSensor, sensorValues, sensorNames, jSensorLookUpTable, 
     *     AmpValueNew, 
     *     lFlagsDefine,
     *     AmpDerivative, AmpSecDerivative, AmpIncIntegral,
     *     AmpDoubleIntegral)
C
      INCLUDE 'ABA_PARAM.INC'

C     time indices
      parameter (iStepTime        = 1,
     *           iTotalTime       = 2,
     *           nTime            = 2)
C     flags passed in for information
      parameter (iInitialization   = 1,
     *           iRegularInc       = 2,
     *           iCuts             = 3,
     *           ikStep            = 4,
     *           nFlagsInfo        = 4)
C     optional flags to be defined
      parameter (iComputeDeriv       = 1,
     *           iComputeSecDeriv    = 2,
     *           iComputeInteg       = 3,
     *           iComputeDoubleInteg = 4,
     *           iStopAnalysis       = 5,
     *           iConcludeStep       = 6,
     *           nFlagsDefine        = 6)
      dimension time(nTime), lFlagsInfo(nFlagsInfo),
     *          lFlagsDefine(nFlagsDefine)
      dimension jSensorLookUpTable(*)
      dimension sensorValues(nSensor), svars(nSvars), props(nProps)
      character*80 sensorNames(nSensor)
      character*80 ampName
      real*8 :: Tc
      real*8 :: error_old
	  real*8 :: Kp, Ki, Kd
	  real*8 :: integral, derivative
      real*8 :: tstar

c     Get the sensor value
      iR_TEMP  = IGetSensorID('TEMP-SENSOR',jSensorLookUpTable)
      temp_sensor = sensorValues(iR_TEMP)
      ampValueNew = 1.0
      Tdef=450.0
	  tstar=(535.0-Tdef)/(-1.9)
	  Kp = 0.12
	  Ki = 0.12
	  Kd = 0.002
	        ! it is best to start with Kd=0.
			! then keep increasing Kp until best choice is found 
			! them increase Kd & ki from zero if needed.
			! note Ki can be quite big if a small stepsize is used (e.g. ki>kp for small time step.)
c     
	  if(time(iTotalTime)==0)then
	     svars(1) = zero
		 svars(2) = zero
	  end if 
      if(time(iTotalTime).le.16.0)then 
          Texp=(30.0*time(iTotalTime))+20.0
      else if(time(iTotalTime).gt.16.0 .AND.  
     &   time(iTotalTime).le.23.0)then
          Texp=(((time(iTotalTime))-(16.0))*5.0)+500.0
      else if(time(iTotalTime).gt.23.0 .AND. 
     &   time(iTotalTime).le.83.0)then
          Texp=535.0
      else if(time(iTotalTime).gt.83.0 .AND. 
     &    time(iTotalTime).le.84.9)then
          Texp=(((time(iTotalTime))-(83.0))*(tstar))+535.0
      else if(time(iTotalTime).gt.84.9)then
          Texp=Tdef
      end if
	  error_old=svars(1)
	  integral=svars(2)
      error=Texp-temp_sensor
	  derivative=(error-error_old)/dt
	  integral=integral+dt*(error+error_old)/2.0
      if(abs(time(iTotalTime) - 16.0) .lt. 1e-4)then
          integral = 0
          derivative = 0
		  Print*, "SWITCH"
      end if
      if(abs(time(iTotalTime) - 23.0) .lt. 1e-4)then
          integral = 0
          derivative = 0
		  Print*, "SWITCH"
      end if
      if(abs(time(iTotalTime) - 83.0) .lt. 1e-4)then
          integral = 0
          derivative = 0
		  Print*, "SWITCH"
      end if
      if(abs(time(iTotalTime) - 84.9) .lt. 1e-4)then
          integral = 0
          derivative = 0
		  Print*, "SWITCH"
      end if
	  derivative=MIN(1.0,derivative)
      ampValueNew = EXP(1*(Kp*error+Ki*integral+Kd*derivative))

      Print*, "Time:", time(iTotalTime)
	  Print*, "Texp:", Texp, "Temp Sensor:", temp_sensor
	  Print*, "DT_old:", error_old, "DT:", error
	  Print*, "AmpValueNew:", ampValueNew, "dt:", dt
	  Print*, "Derivative:", derivative, "Integral:", integral
	  Print*, "-------------"
	  error_old=error
	  svars(1)=error_old
	  svars(2)=integral
      return 
      end
      SUBROUTINE UHARD(SYIELD, HARD, EQPLAS, EQPLASRT, TIME, DTIME,
     & TEMP, DTEMP, NOEL, NPT, LAYER, KSPT, KSTEP, KINC, CMNAME,
     & NSTATV, STATEV, NUMFIELDV, PREDEF, DPREDEF, NUMPROPS, PROPS)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION HARD(3), STATEV(NSTATV), TIME(*), PREDEF(NUMFIELDV),
     & DPREDEF(NUMFIELDV), PROPS(NUMPROPS)
C
      REAL*8 SYIELD, EQPLAS, EQPLASRT, DTIME, TEMP, DTEMP
      INTEGER NOEL, NPT, LAYER, KSPT, KSTEP, KINC, NSTATV, NUMFIELDV,
     & NUMPROPS
C
      REAL*8 A, n, B, m, T0, Tmax, C, SR_0, theta, Norm_SR
C
      A = PROPS(1)
	  B = PROPS(2)
      n = PROPS(3)
      m = PROPS(4)
      Tmax = PROPS(5)
	  T0 = PROPS(6)
	  C = PROPS(7)
	  SR_0 = PROPS(8)
C
      theta=(TEMP-T0)/(Tmax-T0)
	  Norm_SR=EQPLASRT/SR_0
      IF (TEMP.LE.(T0)) THEN
          theta=1.e-5
      else IF(TEMP.GT.(Tmax)) THEN
          theta=1.0
      END IF
      IF (EQPLASRT.LE.(SR_0)) THEN
          Norm_SR=1+1.e-5
      END IF
      SYIELD = (A + B*(EQPLAS)**n)*(1 - (theta)**m)*(1+C*LOG(Norm_SR))
      HARD(1) = B * n * (EQPLAS)**(n-1)*(1 - (theta)**m)*(1+C*LOG(Norm_SR))
      HARD(2) = (C/(Norm_SR*SR_0))*(A + B*(EQPLAS)**n)*(1 - (theta)**m)
      HARD(3) = -(m/(Tmax-T0))*(theta)**(m-1)*(A + B*(EQPLAS)**n)*(1+C*LOG(Norm_SR))
C
      IF (SYIELD.LE.(B/100.)) THEN
          SYIELD = B/100.
          HARD(1) = 0.D0
		  HARD(2) = 0.D0
          HARD(3) = 0.D0
      END IF
C
      RETURN
      END
