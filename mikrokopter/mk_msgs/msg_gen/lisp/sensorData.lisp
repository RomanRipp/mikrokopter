; Auto-generated. Do not edit!


(cl:in-package mk_msgs-msg)


;//! \htmlinclude sensorData.msg.html

(cl:defclass <sensorData> (roslisp-msg-protocol:ros-message)
  ((header
    :reader header
    :initarg :header
    :type std_msgs-msg:Header
    :initform (cl:make-instance 'std_msgs-msg:Header))
   (Longitude
    :reader Longitude
    :initarg :Longitude
    :type cl:float
    :initform 0.0)
   (Latitude
    :reader Latitude
    :initarg :Latitude
    :type cl:float
    :initform 0.0)
   (gpsAltitude
    :reader gpsAltitude
    :initarg :gpsAltitude
    :type cl:float
    :initform 0.0)
   (SatsInUse
    :reader SatsInUse
    :initarg :SatsInUse
    :type cl:integer
    :initform 0)
   (NumberOfWPs
    :reader NumberOfWPs
    :initarg :NumberOfWPs
    :type cl:integer
    :initform 0)
   (CurrentWP
    :reader CurrentWP
    :initarg :CurrentWP
    :type cl:integer
    :initform 0)
   (ErrorCode
    :reader ErrorCode
    :initarg :ErrorCode
    :type cl:integer
    :initform 0)
   (Free
    :reader Free
    :initarg :Free
    :type cl:integer
    :initform 0)
   (PH
    :reader PH
    :initarg :PH
    :type cl:integer
    :initform 0)
   (CH
    :reader CH
    :initarg :CH
    :type cl:integer
    :initform 0)
   (RangeLimit
    :reader RangeLimit
    :initarg :RangeLimit
    :type cl:integer
    :initform 0)
   (NoSerialLink
    :reader NoSerialLink
    :initarg :NoSerialLink
    :type cl:integer
    :initform 0)
   (TargetReached
    :reader TargetReached
    :initarg :TargetReached
    :type cl:integer
    :initform 0)
   (Manual
    :reader Manual
    :initarg :Manual
    :type cl:integer
    :initform 0)
   (GPSOK
    :reader GPSOK
    :initarg :GPSOK
    :type cl:integer
    :initform 0)
   (MotorsOn
    :reader MotorsOn
    :initarg :MotorsOn
    :type cl:integer
    :initform 0)
   (Flying
    :reader Flying
    :initarg :Flying
    :type cl:integer
    :initform 0)
   (LowBat
    :reader LowBat
    :initarg :LowBat
    :type cl:integer
    :initform 0)
   (CareFree
    :reader CareFree
    :initarg :CareFree
    :type cl:integer
    :initform 0)
   (AltHld
    :reader AltHld
    :initarg :AltHld
    :type cl:integer
    :initform 0)
   (Failsafe
    :reader Failsafe
    :initarg :Failsafe
    :type cl:integer
    :initform 0)
   (Altitude
    :reader Altitude
    :initarg :Altitude
    :type cl:float
    :initform 0.0)
   (FlyingTime
    :reader FlyingTime
    :initarg :FlyingTime
    :type cl:float
    :initform 0.0)
   (Battery
    :reader Battery
    :initarg :Battery
    :type cl:float
    :initform 0.0)
   (GroundSpeed
    :reader GroundSpeed
    :initarg :GroundSpeed
    :type cl:float
    :initform 0.0)
   (Heading
    :reader Heading
    :initarg :Heading
    :type cl:float
    :initform 0.0)
   (CompassHeading
    :reader CompassHeading
    :initarg :CompassHeading
    :type cl:float
    :initform 0.0)
   (Nick
    :reader Nick
    :initarg :Nick
    :type cl:float
    :initform 0.0)
   (Roll
    :reader Roll
    :initarg :Roll
    :type cl:float
    :initform 0.0)
   (RCQuality
    :reader RCQuality
    :initarg :RCQuality
    :type cl:float
    :initform 0.0)
   (zSpeed
    :reader zSpeed
    :initarg :zSpeed
    :type cl:float
    :initform 0.0)
   (TargetHoldTime
    :reader TargetHoldTime
    :initarg :TargetHoldTime
    :type cl:float
    :initform 0.0)
   (Gas
    :reader Gas
    :initarg :Gas
    :type cl:float
    :initform 0.0)
   (Current
    :reader Current
    :initarg :Current
    :type cl:float
    :initform 0.0)
   (UsedCapacity
    :reader UsedCapacity
    :initarg :UsedCapacity
    :type cl:float
    :initform 0.0))
)

(cl:defclass sensorData (<sensorData>)
  ())

(cl:defmethod cl:initialize-instance :after ((m <sensorData>) cl:&rest args)
  (cl:declare (cl:ignorable args))
  (cl:unless (cl:typep m 'sensorData)
    (roslisp-msg-protocol:msg-deprecation-warning "using old message class name mk_msgs-msg:<sensorData> is deprecated: use mk_msgs-msg:sensorData instead.")))

(cl:ensure-generic-function 'header-val :lambda-list '(m))
(cl:defmethod header-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:header-val is deprecated.  Use mk_msgs-msg:header instead.")
  (header m))

(cl:ensure-generic-function 'Longitude-val :lambda-list '(m))
(cl:defmethod Longitude-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:Longitude-val is deprecated.  Use mk_msgs-msg:Longitude instead.")
  (Longitude m))

(cl:ensure-generic-function 'Latitude-val :lambda-list '(m))
(cl:defmethod Latitude-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:Latitude-val is deprecated.  Use mk_msgs-msg:Latitude instead.")
  (Latitude m))

(cl:ensure-generic-function 'gpsAltitude-val :lambda-list '(m))
(cl:defmethod gpsAltitude-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:gpsAltitude-val is deprecated.  Use mk_msgs-msg:gpsAltitude instead.")
  (gpsAltitude m))

(cl:ensure-generic-function 'SatsInUse-val :lambda-list '(m))
(cl:defmethod SatsInUse-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:SatsInUse-val is deprecated.  Use mk_msgs-msg:SatsInUse instead.")
  (SatsInUse m))

(cl:ensure-generic-function 'NumberOfWPs-val :lambda-list '(m))
(cl:defmethod NumberOfWPs-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:NumberOfWPs-val is deprecated.  Use mk_msgs-msg:NumberOfWPs instead.")
  (NumberOfWPs m))

(cl:ensure-generic-function 'CurrentWP-val :lambda-list '(m))
(cl:defmethod CurrentWP-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:CurrentWP-val is deprecated.  Use mk_msgs-msg:CurrentWP instead.")
  (CurrentWP m))

(cl:ensure-generic-function 'ErrorCode-val :lambda-list '(m))
(cl:defmethod ErrorCode-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:ErrorCode-val is deprecated.  Use mk_msgs-msg:ErrorCode instead.")
  (ErrorCode m))

(cl:ensure-generic-function 'Free-val :lambda-list '(m))
(cl:defmethod Free-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:Free-val is deprecated.  Use mk_msgs-msg:Free instead.")
  (Free m))

(cl:ensure-generic-function 'PH-val :lambda-list '(m))
(cl:defmethod PH-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:PH-val is deprecated.  Use mk_msgs-msg:PH instead.")
  (PH m))

(cl:ensure-generic-function 'CH-val :lambda-list '(m))
(cl:defmethod CH-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:CH-val is deprecated.  Use mk_msgs-msg:CH instead.")
  (CH m))

(cl:ensure-generic-function 'RangeLimit-val :lambda-list '(m))
(cl:defmethod RangeLimit-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:RangeLimit-val is deprecated.  Use mk_msgs-msg:RangeLimit instead.")
  (RangeLimit m))

(cl:ensure-generic-function 'NoSerialLink-val :lambda-list '(m))
(cl:defmethod NoSerialLink-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:NoSerialLink-val is deprecated.  Use mk_msgs-msg:NoSerialLink instead.")
  (NoSerialLink m))

(cl:ensure-generic-function 'TargetReached-val :lambda-list '(m))
(cl:defmethod TargetReached-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:TargetReached-val is deprecated.  Use mk_msgs-msg:TargetReached instead.")
  (TargetReached m))

(cl:ensure-generic-function 'Manual-val :lambda-list '(m))
(cl:defmethod Manual-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:Manual-val is deprecated.  Use mk_msgs-msg:Manual instead.")
  (Manual m))

(cl:ensure-generic-function 'GPSOK-val :lambda-list '(m))
(cl:defmethod GPSOK-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:GPSOK-val is deprecated.  Use mk_msgs-msg:GPSOK instead.")
  (GPSOK m))

(cl:ensure-generic-function 'MotorsOn-val :lambda-list '(m))
(cl:defmethod MotorsOn-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:MotorsOn-val is deprecated.  Use mk_msgs-msg:MotorsOn instead.")
  (MotorsOn m))

(cl:ensure-generic-function 'Flying-val :lambda-list '(m))
(cl:defmethod Flying-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:Flying-val is deprecated.  Use mk_msgs-msg:Flying instead.")
  (Flying m))

(cl:ensure-generic-function 'LowBat-val :lambda-list '(m))
(cl:defmethod LowBat-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:LowBat-val is deprecated.  Use mk_msgs-msg:LowBat instead.")
  (LowBat m))

(cl:ensure-generic-function 'CareFree-val :lambda-list '(m))
(cl:defmethod CareFree-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:CareFree-val is deprecated.  Use mk_msgs-msg:CareFree instead.")
  (CareFree m))

(cl:ensure-generic-function 'AltHld-val :lambda-list '(m))
(cl:defmethod AltHld-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:AltHld-val is deprecated.  Use mk_msgs-msg:AltHld instead.")
  (AltHld m))

(cl:ensure-generic-function 'Failsafe-val :lambda-list '(m))
(cl:defmethod Failsafe-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:Failsafe-val is deprecated.  Use mk_msgs-msg:Failsafe instead.")
  (Failsafe m))

(cl:ensure-generic-function 'Altitude-val :lambda-list '(m))
(cl:defmethod Altitude-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:Altitude-val is deprecated.  Use mk_msgs-msg:Altitude instead.")
  (Altitude m))

(cl:ensure-generic-function 'FlyingTime-val :lambda-list '(m))
(cl:defmethod FlyingTime-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:FlyingTime-val is deprecated.  Use mk_msgs-msg:FlyingTime instead.")
  (FlyingTime m))

(cl:ensure-generic-function 'Battery-val :lambda-list '(m))
(cl:defmethod Battery-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:Battery-val is deprecated.  Use mk_msgs-msg:Battery instead.")
  (Battery m))

(cl:ensure-generic-function 'GroundSpeed-val :lambda-list '(m))
(cl:defmethod GroundSpeed-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:GroundSpeed-val is deprecated.  Use mk_msgs-msg:GroundSpeed instead.")
  (GroundSpeed m))

(cl:ensure-generic-function 'Heading-val :lambda-list '(m))
(cl:defmethod Heading-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:Heading-val is deprecated.  Use mk_msgs-msg:Heading instead.")
  (Heading m))

(cl:ensure-generic-function 'CompassHeading-val :lambda-list '(m))
(cl:defmethod CompassHeading-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:CompassHeading-val is deprecated.  Use mk_msgs-msg:CompassHeading instead.")
  (CompassHeading m))

(cl:ensure-generic-function 'Nick-val :lambda-list '(m))
(cl:defmethod Nick-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:Nick-val is deprecated.  Use mk_msgs-msg:Nick instead.")
  (Nick m))

(cl:ensure-generic-function 'Roll-val :lambda-list '(m))
(cl:defmethod Roll-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:Roll-val is deprecated.  Use mk_msgs-msg:Roll instead.")
  (Roll m))

(cl:ensure-generic-function 'RCQuality-val :lambda-list '(m))
(cl:defmethod RCQuality-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:RCQuality-val is deprecated.  Use mk_msgs-msg:RCQuality instead.")
  (RCQuality m))

(cl:ensure-generic-function 'zSpeed-val :lambda-list '(m))
(cl:defmethod zSpeed-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:zSpeed-val is deprecated.  Use mk_msgs-msg:zSpeed instead.")
  (zSpeed m))

(cl:ensure-generic-function 'TargetHoldTime-val :lambda-list '(m))
(cl:defmethod TargetHoldTime-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:TargetHoldTime-val is deprecated.  Use mk_msgs-msg:TargetHoldTime instead.")
  (TargetHoldTime m))

(cl:ensure-generic-function 'Gas-val :lambda-list '(m))
(cl:defmethod Gas-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:Gas-val is deprecated.  Use mk_msgs-msg:Gas instead.")
  (Gas m))

(cl:ensure-generic-function 'Current-val :lambda-list '(m))
(cl:defmethod Current-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:Current-val is deprecated.  Use mk_msgs-msg:Current instead.")
  (Current m))

(cl:ensure-generic-function 'UsedCapacity-val :lambda-list '(m))
(cl:defmethod UsedCapacity-val ((m <sensorData>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-msg:UsedCapacity-val is deprecated.  Use mk_msgs-msg:UsedCapacity instead.")
  (UsedCapacity m))
(cl:defmethod roslisp-msg-protocol:serialize ((msg <sensorData>) ostream)
  "Serializes a message object of type '<sensorData>"
  (roslisp-msg-protocol:serialize (cl:slot-value msg 'header) ostream)
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'Longitude))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'Latitude))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'gpsAltitude))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let* ((signed (cl:slot-value msg 'SatsInUse)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'NumberOfWPs)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'CurrentWP)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'ErrorCode)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'Free)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'PH)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'CH)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'RangeLimit)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'NoSerialLink)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'TargetReached)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'Manual)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'GPSOK)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'MotorsOn)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'Flying)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'LowBat)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'CareFree)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'AltHld)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'Failsafe)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'Altitude))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'FlyingTime))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'Battery))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'GroundSpeed))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'Heading))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'CompassHeading))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'Nick))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'Roll))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'RCQuality))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'zSpeed))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'TargetHoldTime))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'Gas))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'Current))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'UsedCapacity))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
)
(cl:defmethod roslisp-msg-protocol:deserialize ((msg <sensorData>) istream)
  "Deserializes a message object of type '<sensorData>"
  (roslisp-msg-protocol:deserialize (cl:slot-value msg 'header) istream)
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'Longitude) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'Latitude) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'gpsAltitude) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'SatsInUse) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'NumberOfWPs) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'CurrentWP) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'ErrorCode) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'Free) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'PH) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'CH) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'RangeLimit) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'NoSerialLink) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'TargetReached) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'Manual) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'GPSOK) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'MotorsOn) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'Flying) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'LowBat) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'CareFree) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'AltHld) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'Failsafe) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'Altitude) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'FlyingTime) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'Battery) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'GroundSpeed) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'Heading) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'CompassHeading) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'Nick) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'Roll) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'RCQuality) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'zSpeed) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'TargetHoldTime) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'Gas) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'Current) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'UsedCapacity) (roslisp-utils:decode-single-float-bits bits)))
  msg
)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql '<sensorData>)))
  "Returns string type for a message object of type '<sensorData>"
  "mk_msgs/sensorData")
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'sensorData)))
  "Returns string type for a message object of type 'sensorData"
  "mk_msgs/sensorData")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql '<sensorData>)))
  "Returns md5sum for a message object of type '<sensorData>"
  "7e8cef3e80536e777b04973b40e71d35")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql 'sensorData)))
  "Returns md5sum for a message object of type 'sensorData"
  "7e8cef3e80536e777b04973b40e71d35")
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql '<sensorData>)))
  "Returns full string definition for message of type '<sensorData>"
  (cl:format cl:nil "Header header~%~%float32 Longitude~%float32 Latitude~%float32 gpsAltitude~%~%int32 SatsInUse~%int32 NumberOfWPs~%int32 CurrentWP~%int32 ErrorCode~%#int32 NCMode~%#int32 FCMode1~%#int32 FCMode2~%~%int32 Free~%int32 PH~%int32 CH~%int32 RangeLimit~%int32 NoSerialLink~%int32 TargetReached~%int32 Manual~%int32 GPSOK~%~%int32 MotorsOn~%int32 Flying~%#int32 Calibrate~%#int32 Start~%#int32 EmergencyLanding~%int32 LowBat~%#int32 VarioTimUp~%#int32 VarioTrimDown~%~%int32 CareFree~%int32 AltHld~%int32 Failsafe ~%~%float32 Altitude~%float32 FlyingTime~%float32 Battery~%float32 GroundSpeed~%float32 Heading~%float32 CompassHeading~%float32 Nick~%float32 Roll~%float32 RCQuality~%float32 zSpeed~%float32 TargetHoldTime~%float32 Gas~%float32 Current~%float32 UsedCapacity~%~%~%================================================================================~%MSG: std_msgs/Header~%# Standard metadata for higher-level stamped data types.~%# This is generally used to communicate timestamped data ~%# in a particular coordinate frame.~%# ~%# sequence ID: consecutively increasing ID ~%uint32 seq~%#Two-integer timestamp that is expressed as:~%# * stamp.secs: seconds (stamp_secs) since epoch~%# * stamp.nsecs: nanoseconds since stamp_secs~%# time-handling sugar is provided by the client library~%time stamp~%#Frame this data is associated with~%# 0: no frame~%# 1: global frame~%string frame_id~%~%~%"))
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql 'sensorData)))
  "Returns full string definition for message of type 'sensorData"
  (cl:format cl:nil "Header header~%~%float32 Longitude~%float32 Latitude~%float32 gpsAltitude~%~%int32 SatsInUse~%int32 NumberOfWPs~%int32 CurrentWP~%int32 ErrorCode~%#int32 NCMode~%#int32 FCMode1~%#int32 FCMode2~%~%int32 Free~%int32 PH~%int32 CH~%int32 RangeLimit~%int32 NoSerialLink~%int32 TargetReached~%int32 Manual~%int32 GPSOK~%~%int32 MotorsOn~%int32 Flying~%#int32 Calibrate~%#int32 Start~%#int32 EmergencyLanding~%int32 LowBat~%#int32 VarioTimUp~%#int32 VarioTrimDown~%~%int32 CareFree~%int32 AltHld~%int32 Failsafe ~%~%float32 Altitude~%float32 FlyingTime~%float32 Battery~%float32 GroundSpeed~%float32 Heading~%float32 CompassHeading~%float32 Nick~%float32 Roll~%float32 RCQuality~%float32 zSpeed~%float32 TargetHoldTime~%float32 Gas~%float32 Current~%float32 UsedCapacity~%~%~%================================================================================~%MSG: std_msgs/Header~%# Standard metadata for higher-level stamped data types.~%# This is generally used to communicate timestamped data ~%# in a particular coordinate frame.~%# ~%# sequence ID: consecutively increasing ID ~%uint32 seq~%#Two-integer timestamp that is expressed as:~%# * stamp.secs: seconds (stamp_secs) since epoch~%# * stamp.nsecs: nanoseconds since stamp_secs~%# time-handling sugar is provided by the client library~%time stamp~%#Frame this data is associated with~%# 0: no frame~%# 1: global frame~%string frame_id~%~%~%"))
(cl:defmethod roslisp-msg-protocol:serialization-length ((msg <sensorData>))
  (cl:+ 0
     (roslisp-msg-protocol:serialization-length (cl:slot-value msg 'header))
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
     4
))
(cl:defmethod roslisp-msg-protocol:ros-message-to-list ((msg <sensorData>))
  "Converts a ROS message object to a list"
  (cl:list 'sensorData
    (cl:cons ':header (header msg))
    (cl:cons ':Longitude (Longitude msg))
    (cl:cons ':Latitude (Latitude msg))
    (cl:cons ':gpsAltitude (gpsAltitude msg))
    (cl:cons ':SatsInUse (SatsInUse msg))
    (cl:cons ':NumberOfWPs (NumberOfWPs msg))
    (cl:cons ':CurrentWP (CurrentWP msg))
    (cl:cons ':ErrorCode (ErrorCode msg))
    (cl:cons ':Free (Free msg))
    (cl:cons ':PH (PH msg))
    (cl:cons ':CH (CH msg))
    (cl:cons ':RangeLimit (RangeLimit msg))
    (cl:cons ':NoSerialLink (NoSerialLink msg))
    (cl:cons ':TargetReached (TargetReached msg))
    (cl:cons ':Manual (Manual msg))
    (cl:cons ':GPSOK (GPSOK msg))
    (cl:cons ':MotorsOn (MotorsOn msg))
    (cl:cons ':Flying (Flying msg))
    (cl:cons ':LowBat (LowBat msg))
    (cl:cons ':CareFree (CareFree msg))
    (cl:cons ':AltHld (AltHld msg))
    (cl:cons ':Failsafe (Failsafe msg))
    (cl:cons ':Altitude (Altitude msg))
    (cl:cons ':FlyingTime (FlyingTime msg))
    (cl:cons ':Battery (Battery msg))
    (cl:cons ':GroundSpeed (GroundSpeed msg))
    (cl:cons ':Heading (Heading msg))
    (cl:cons ':CompassHeading (CompassHeading msg))
    (cl:cons ':Nick (Nick msg))
    (cl:cons ':Roll (Roll msg))
    (cl:cons ':RCQuality (RCQuality msg))
    (cl:cons ':zSpeed (zSpeed msg))
    (cl:cons ':TargetHoldTime (TargetHoldTime msg))
    (cl:cons ':Gas (Gas msg))
    (cl:cons ':Current (Current msg))
    (cl:cons ':UsedCapacity (UsedCapacity msg))
))
