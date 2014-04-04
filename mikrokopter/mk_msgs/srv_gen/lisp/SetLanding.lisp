; Auto-generated. Do not edit!


(cl:in-package mk_msgs-srv)


;//! \htmlinclude SetLanding-request.msg.html

(cl:defclass <SetLanding-request> (roslisp-msg-protocol:ros-message)
  ((height
    :reader height
    :initarg :height
    :type cl:float
    :initform 0.0)
   (lat
    :reader lat
    :initarg :lat
    :type cl:float
    :initform 0.0)
   (lon
    :reader lon
    :initarg :lon
    :type cl:float
    :initform 0.0))
)

(cl:defclass SetLanding-request (<SetLanding-request>)
  ())

(cl:defmethod cl:initialize-instance :after ((m <SetLanding-request>) cl:&rest args)
  (cl:declare (cl:ignorable args))
  (cl:unless (cl:typep m 'SetLanding-request)
    (roslisp-msg-protocol:msg-deprecation-warning "using old message class name mk_msgs-srv:<SetLanding-request> is deprecated: use mk_msgs-srv:SetLanding-request instead.")))

(cl:ensure-generic-function 'height-val :lambda-list '(m))
(cl:defmethod height-val ((m <SetLanding-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:height-val is deprecated.  Use mk_msgs-srv:height instead.")
  (height m))

(cl:ensure-generic-function 'lat-val :lambda-list '(m))
(cl:defmethod lat-val ((m <SetLanding-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:lat-val is deprecated.  Use mk_msgs-srv:lat instead.")
  (lat m))

(cl:ensure-generic-function 'lon-val :lambda-list '(m))
(cl:defmethod lon-val ((m <SetLanding-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:lon-val is deprecated.  Use mk_msgs-srv:lon instead.")
  (lon m))
(cl:defmethod roslisp-msg-protocol:serialize ((msg <SetLanding-request>) ostream)
  "Serializes a message object of type '<SetLanding-request>"
  (cl:let ((bits (roslisp-utils:encode-double-float-bits (cl:slot-value msg 'height))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 32) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 40) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 48) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 56) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-double-float-bits (cl:slot-value msg 'lat))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 32) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 40) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 48) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 56) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-double-float-bits (cl:slot-value msg 'lon))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 32) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 40) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 48) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 56) bits) ostream))
)
(cl:defmethod roslisp-msg-protocol:deserialize ((msg <SetLanding-request>) istream)
  "Deserializes a message object of type '<SetLanding-request>"
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 32) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 40) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 48) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 56) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'height) (roslisp-utils:decode-double-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 32) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 40) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 48) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 56) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'lat) (roslisp-utils:decode-double-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 32) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 40) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 48) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 56) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'lon) (roslisp-utils:decode-double-float-bits bits)))
  msg
)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql '<SetLanding-request>)))
  "Returns string type for a service object of type '<SetLanding-request>"
  "mk_msgs/SetLandingRequest")
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'SetLanding-request)))
  "Returns string type for a service object of type 'SetLanding-request"
  "mk_msgs/SetLandingRequest")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql '<SetLanding-request>)))
  "Returns md5sum for a message object of type '<SetLanding-request>"
  "f85bad350c9de80d1ec2e364a6cdc138")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql 'SetLanding-request)))
  "Returns md5sum for a message object of type 'SetLanding-request"
  "f85bad350c9de80d1ec2e364a6cdc138")
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql '<SetLanding-request>)))
  "Returns full string definition for message of type '<SetLanding-request>"
  (cl:format cl:nil "~%~%~%float64 height~%float64 lat~%float64 lon~%~%~%"))
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql 'SetLanding-request)))
  "Returns full string definition for message of type 'SetLanding-request"
  (cl:format cl:nil "~%~%~%float64 height~%float64 lat~%float64 lon~%~%~%"))
(cl:defmethod roslisp-msg-protocol:serialization-length ((msg <SetLanding-request>))
  (cl:+ 0
     8
     8
     8
))
(cl:defmethod roslisp-msg-protocol:ros-message-to-list ((msg <SetLanding-request>))
  "Converts a ROS message object to a list"
  (cl:list 'SetLanding-request
    (cl:cons ':height (height msg))
    (cl:cons ':lat (lat msg))
    (cl:cons ':lon (lon msg))
))
;//! \htmlinclude SetLanding-response.msg.html

(cl:defclass <SetLanding-response> (roslisp-msg-protocol:ros-message)
  ((ack
    :reader ack
    :initarg :ack
    :type cl:fixnum
    :initform 0))
)

(cl:defclass SetLanding-response (<SetLanding-response>)
  ())

(cl:defmethod cl:initialize-instance :after ((m <SetLanding-response>) cl:&rest args)
  (cl:declare (cl:ignorable args))
  (cl:unless (cl:typep m 'SetLanding-response)
    (roslisp-msg-protocol:msg-deprecation-warning "using old message class name mk_msgs-srv:<SetLanding-response> is deprecated: use mk_msgs-srv:SetLanding-response instead.")))

(cl:ensure-generic-function 'ack-val :lambda-list '(m))
(cl:defmethod ack-val ((m <SetLanding-response>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:ack-val is deprecated.  Use mk_msgs-srv:ack instead.")
  (ack m))
(cl:defmethod roslisp-msg-protocol:serialize ((msg <SetLanding-response>) ostream)
  "Serializes a message object of type '<SetLanding-response>"
  (cl:let* ((signed (cl:slot-value msg 'ack)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 65536) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    )
)
(cl:defmethod roslisp-msg-protocol:deserialize ((msg <SetLanding-response>) istream)
  "Deserializes a message object of type '<SetLanding-response>"
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'ack) (cl:if (cl:< unsigned 32768) unsigned (cl:- unsigned 65536))))
  msg
)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql '<SetLanding-response>)))
  "Returns string type for a service object of type '<SetLanding-response>"
  "mk_msgs/SetLandingResponse")
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'SetLanding-response)))
  "Returns string type for a service object of type 'SetLanding-response"
  "mk_msgs/SetLandingResponse")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql '<SetLanding-response>)))
  "Returns md5sum for a message object of type '<SetLanding-response>"
  "f85bad350c9de80d1ec2e364a6cdc138")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql 'SetLanding-response)))
  "Returns md5sum for a message object of type 'SetLanding-response"
  "f85bad350c9de80d1ec2e364a6cdc138")
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql '<SetLanding-response>)))
  "Returns full string definition for message of type '<SetLanding-response>"
  (cl:format cl:nil "~%~%~%int16 ack~%~%~%~%"))
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql 'SetLanding-response)))
  "Returns full string definition for message of type 'SetLanding-response"
  (cl:format cl:nil "~%~%~%int16 ack~%~%~%~%"))
(cl:defmethod roslisp-msg-protocol:serialization-length ((msg <SetLanding-response>))
  (cl:+ 0
     2
))
(cl:defmethod roslisp-msg-protocol:ros-message-to-list ((msg <SetLanding-response>))
  "Converts a ROS message object to a list"
  (cl:list 'SetLanding-response
    (cl:cons ':ack (ack msg))
))
(cl:defmethod roslisp-msg-protocol:service-request-type ((msg (cl:eql 'SetLanding)))
  'SetLanding-request)
(cl:defmethod roslisp-msg-protocol:service-response-type ((msg (cl:eql 'SetLanding)))
  'SetLanding-response)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'SetLanding)))
  "Returns string type for a service object of type '<SetLanding>"
  "mk_msgs/SetLanding")