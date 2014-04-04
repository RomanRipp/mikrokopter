; Auto-generated. Do not edit!


(cl:in-package mk_msgs-srv)


;//! \htmlinclude SetGPSWaypoint-request.msg.html

(cl:defclass <SetGPSWaypoint-request> (roslisp-msg-protocol:ros-message)
  ((index
    :reader index
    :initarg :index
    :type cl:integer
    :initform 0)
   (latitude
    :reader latitude
    :initarg :latitude
    :type cl:float
    :initform 0.0)
   (longitude
    :reader longitude
    :initarg :longitude
    :type cl:float
    :initform 0.0)
   (altitude
    :reader altitude
    :initarg :altitude
    :type cl:float
    :initform 0.0)
   (status
    :reader status
    :initarg :status
    :type cl:integer
    :initform 0)
   (holdTime
    :reader holdTime
    :initarg :holdTime
    :type cl:integer
    :initform 0)
   (altitudeRate
    :reader altitudeRate
    :initarg :altitudeRate
    :type cl:float
    :initform 0.0)
   (speed
    :reader speed
    :initarg :speed
    :type cl:float
    :initform 0.0)
   (heading
    :reader heading
    :initarg :heading
    :type cl:float
    :initform 0.0)
   (toleranceRadius
    :reader toleranceRadius
    :initarg :toleranceRadius
    :type cl:integer
    :initform 0)
   (eventChannel
    :reader eventChannel
    :initarg :eventChannel
    :type cl:integer
    :initform 0)
   (eventFlag
    :reader eventFlag
    :initarg :eventFlag
    :type cl:integer
    :initform 0)
   (type
    :reader type
    :initarg :type
    :type cl:integer
    :initform 0)
   (camAngle
    :reader camAngle
    :initarg :camAngle
    :type cl:integer
    :initform 0))
)

(cl:defclass SetGPSWaypoint-request (<SetGPSWaypoint-request>)
  ())

(cl:defmethod cl:initialize-instance :after ((m <SetGPSWaypoint-request>) cl:&rest args)
  (cl:declare (cl:ignorable args))
  (cl:unless (cl:typep m 'SetGPSWaypoint-request)
    (roslisp-msg-protocol:msg-deprecation-warning "using old message class name mk_msgs-srv:<SetGPSWaypoint-request> is deprecated: use mk_msgs-srv:SetGPSWaypoint-request instead.")))

(cl:ensure-generic-function 'index-val :lambda-list '(m))
(cl:defmethod index-val ((m <SetGPSWaypoint-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:index-val is deprecated.  Use mk_msgs-srv:index instead.")
  (index m))

(cl:ensure-generic-function 'latitude-val :lambda-list '(m))
(cl:defmethod latitude-val ((m <SetGPSWaypoint-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:latitude-val is deprecated.  Use mk_msgs-srv:latitude instead.")
  (latitude m))

(cl:ensure-generic-function 'longitude-val :lambda-list '(m))
(cl:defmethod longitude-val ((m <SetGPSWaypoint-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:longitude-val is deprecated.  Use mk_msgs-srv:longitude instead.")
  (longitude m))

(cl:ensure-generic-function 'altitude-val :lambda-list '(m))
(cl:defmethod altitude-val ((m <SetGPSWaypoint-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:altitude-val is deprecated.  Use mk_msgs-srv:altitude instead.")
  (altitude m))

(cl:ensure-generic-function 'status-val :lambda-list '(m))
(cl:defmethod status-val ((m <SetGPSWaypoint-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:status-val is deprecated.  Use mk_msgs-srv:status instead.")
  (status m))

(cl:ensure-generic-function 'holdTime-val :lambda-list '(m))
(cl:defmethod holdTime-val ((m <SetGPSWaypoint-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:holdTime-val is deprecated.  Use mk_msgs-srv:holdTime instead.")
  (holdTime m))

(cl:ensure-generic-function 'altitudeRate-val :lambda-list '(m))
(cl:defmethod altitudeRate-val ((m <SetGPSWaypoint-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:altitudeRate-val is deprecated.  Use mk_msgs-srv:altitudeRate instead.")
  (altitudeRate m))

(cl:ensure-generic-function 'speed-val :lambda-list '(m))
(cl:defmethod speed-val ((m <SetGPSWaypoint-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:speed-val is deprecated.  Use mk_msgs-srv:speed instead.")
  (speed m))

(cl:ensure-generic-function 'heading-val :lambda-list '(m))
(cl:defmethod heading-val ((m <SetGPSWaypoint-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:heading-val is deprecated.  Use mk_msgs-srv:heading instead.")
  (heading m))

(cl:ensure-generic-function 'toleranceRadius-val :lambda-list '(m))
(cl:defmethod toleranceRadius-val ((m <SetGPSWaypoint-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:toleranceRadius-val is deprecated.  Use mk_msgs-srv:toleranceRadius instead.")
  (toleranceRadius m))

(cl:ensure-generic-function 'eventChannel-val :lambda-list '(m))
(cl:defmethod eventChannel-val ((m <SetGPSWaypoint-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:eventChannel-val is deprecated.  Use mk_msgs-srv:eventChannel instead.")
  (eventChannel m))

(cl:ensure-generic-function 'eventFlag-val :lambda-list '(m))
(cl:defmethod eventFlag-val ((m <SetGPSWaypoint-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:eventFlag-val is deprecated.  Use mk_msgs-srv:eventFlag instead.")
  (eventFlag m))

(cl:ensure-generic-function 'type-val :lambda-list '(m))
(cl:defmethod type-val ((m <SetGPSWaypoint-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:type-val is deprecated.  Use mk_msgs-srv:type instead.")
  (type m))

(cl:ensure-generic-function 'camAngle-val :lambda-list '(m))
(cl:defmethod camAngle-val ((m <SetGPSWaypoint-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:camAngle-val is deprecated.  Use mk_msgs-srv:camAngle instead.")
  (camAngle m))
(cl:defmethod roslisp-msg-protocol:serialize ((msg <SetGPSWaypoint-request>) ostream)
  "Serializes a message object of type '<SetGPSWaypoint-request>"
  (cl:let* ((signed (cl:slot-value msg 'index)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'latitude))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'longitude))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'altitude))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let* ((signed (cl:slot-value msg 'status)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'holdTime)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'altitudeRate))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'speed))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'heading))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let* ((signed (cl:slot-value msg 'toleranceRadius)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'eventChannel)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'eventFlag)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'type)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'camAngle)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
)
(cl:defmethod roslisp-msg-protocol:deserialize ((msg <SetGPSWaypoint-request>) istream)
  "Deserializes a message object of type '<SetGPSWaypoint-request>"
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'index) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'latitude) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'longitude) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'altitude) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'status) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'holdTime) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'altitudeRate) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'speed) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'heading) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'toleranceRadius) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'eventChannel) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'eventFlag) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'type) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'camAngle) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
  msg
)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql '<SetGPSWaypoint-request>)))
  "Returns string type for a service object of type '<SetGPSWaypoint-request>"
  "mk_msgs/SetGPSWaypointRequest")
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'SetGPSWaypoint-request)))
  "Returns string type for a service object of type 'SetGPSWaypoint-request"
  "mk_msgs/SetGPSWaypointRequest")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql '<SetGPSWaypoint-request>)))
  "Returns md5sum for a message object of type '<SetGPSWaypoint-request>"
  "167ab490b44feaf311e87f25a5b4ea94")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql 'SetGPSWaypoint-request)))
  "Returns md5sum for a message object of type 'SetGPSWaypoint-request"
  "167ab490b44feaf311e87f25a5b4ea94")
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql '<SetGPSWaypoint-request>)))
  "Returns full string definition for message of type '<SetGPSWaypoint-request>"
  (cl:format cl:nil "~%~%~%int32 index~%float32 latitude~%float32 longitude~%float32 altitude~%int32 status~%int32 holdTime~%float32 altitudeRate~%float32 speed~%~%float32 heading~%int32 toleranceRadius~%int32 eventChannel~%int32 eventFlag~%int32 type~%int32 camAngle~%~%~%"))
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql 'SetGPSWaypoint-request)))
  "Returns full string definition for message of type 'SetGPSWaypoint-request"
  (cl:format cl:nil "~%~%~%int32 index~%float32 latitude~%float32 longitude~%float32 altitude~%int32 status~%int32 holdTime~%float32 altitudeRate~%float32 speed~%~%float32 heading~%int32 toleranceRadius~%int32 eventChannel~%int32 eventFlag~%int32 type~%int32 camAngle~%~%~%"))
(cl:defmethod roslisp-msg-protocol:serialization-length ((msg <SetGPSWaypoint-request>))
  (cl:+ 0
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
(cl:defmethod roslisp-msg-protocol:ros-message-to-list ((msg <SetGPSWaypoint-request>))
  "Converts a ROS message object to a list"
  (cl:list 'SetGPSWaypoint-request
    (cl:cons ':index (index msg))
    (cl:cons ':latitude (latitude msg))
    (cl:cons ':longitude (longitude msg))
    (cl:cons ':altitude (altitude msg))
    (cl:cons ':status (status msg))
    (cl:cons ':holdTime (holdTime msg))
    (cl:cons ':altitudeRate (altitudeRate msg))
    (cl:cons ':speed (speed msg))
    (cl:cons ':heading (heading msg))
    (cl:cons ':toleranceRadius (toleranceRadius msg))
    (cl:cons ':eventChannel (eventChannel msg))
    (cl:cons ':eventFlag (eventFlag msg))
    (cl:cons ':type (type msg))
    (cl:cons ':camAngle (camAngle msg))
))
;//! \htmlinclude SetGPSWaypoint-response.msg.html

(cl:defclass <SetGPSWaypoint-response> (roslisp-msg-protocol:ros-message)
  ((ack
    :reader ack
    :initarg :ack
    :type cl:fixnum
    :initform 0))
)

(cl:defclass SetGPSWaypoint-response (<SetGPSWaypoint-response>)
  ())

(cl:defmethod cl:initialize-instance :after ((m <SetGPSWaypoint-response>) cl:&rest args)
  (cl:declare (cl:ignorable args))
  (cl:unless (cl:typep m 'SetGPSWaypoint-response)
    (roslisp-msg-protocol:msg-deprecation-warning "using old message class name mk_msgs-srv:<SetGPSWaypoint-response> is deprecated: use mk_msgs-srv:SetGPSWaypoint-response instead.")))

(cl:ensure-generic-function 'ack-val :lambda-list '(m))
(cl:defmethod ack-val ((m <SetGPSWaypoint-response>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:ack-val is deprecated.  Use mk_msgs-srv:ack instead.")
  (ack m))
(cl:defmethod roslisp-msg-protocol:serialize ((msg <SetGPSWaypoint-response>) ostream)
  "Serializes a message object of type '<SetGPSWaypoint-response>"
  (cl:let* ((signed (cl:slot-value msg 'ack)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 65536) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    )
)
(cl:defmethod roslisp-msg-protocol:deserialize ((msg <SetGPSWaypoint-response>) istream)
  "Deserializes a message object of type '<SetGPSWaypoint-response>"
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'ack) (cl:if (cl:< unsigned 32768) unsigned (cl:- unsigned 65536))))
  msg
)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql '<SetGPSWaypoint-response>)))
  "Returns string type for a service object of type '<SetGPSWaypoint-response>"
  "mk_msgs/SetGPSWaypointResponse")
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'SetGPSWaypoint-response)))
  "Returns string type for a service object of type 'SetGPSWaypoint-response"
  "mk_msgs/SetGPSWaypointResponse")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql '<SetGPSWaypoint-response>)))
  "Returns md5sum for a message object of type '<SetGPSWaypoint-response>"
  "167ab490b44feaf311e87f25a5b4ea94")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql 'SetGPSWaypoint-response)))
  "Returns md5sum for a message object of type 'SetGPSWaypoint-response"
  "167ab490b44feaf311e87f25a5b4ea94")
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql '<SetGPSWaypoint-response>)))
  "Returns full string definition for message of type '<SetGPSWaypoint-response>"
  (cl:format cl:nil "~%~%~%int16 ack~%~%~%~%"))
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql 'SetGPSWaypoint-response)))
  "Returns full string definition for message of type 'SetGPSWaypoint-response"
  (cl:format cl:nil "~%~%~%int16 ack~%~%~%~%"))
(cl:defmethod roslisp-msg-protocol:serialization-length ((msg <SetGPSWaypoint-response>))
  (cl:+ 0
     2
))
(cl:defmethod roslisp-msg-protocol:ros-message-to-list ((msg <SetGPSWaypoint-response>))
  "Converts a ROS message object to a list"
  (cl:list 'SetGPSWaypoint-response
    (cl:cons ':ack (ack msg))
))
(cl:defmethod roslisp-msg-protocol:service-request-type ((msg (cl:eql 'SetGPSWaypoint)))
  'SetGPSWaypoint-request)
(cl:defmethod roslisp-msg-protocol:service-response-type ((msg (cl:eql 'SetGPSWaypoint)))
  'SetGPSWaypoint-response)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'SetGPSWaypoint)))
  "Returns string type for a service object of type '<SetGPSWaypoint>"
  "mk_msgs/SetGPSWaypoint")