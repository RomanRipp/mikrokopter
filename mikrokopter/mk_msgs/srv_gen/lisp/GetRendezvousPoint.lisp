; Auto-generated. Do not edit!


(cl:in-package mk_msgs-srv)


;//! \htmlinclude GetRendezvousPoint-request.msg.html

(cl:defclass <GetRendezvousPoint-request> (roslisp-msg-protocol:ros-message)
  ((currLatitude
    :reader currLatitude
    :initarg :currLatitude
    :type cl:float
    :initform 0.0)
   (currLongitude
    :reader currLongitude
    :initarg :currLongitude
    :type cl:float
    :initform 0.0))
)

(cl:defclass GetRendezvousPoint-request (<GetRendezvousPoint-request>)
  ())

(cl:defmethod cl:initialize-instance :after ((m <GetRendezvousPoint-request>) cl:&rest args)
  (cl:declare (cl:ignorable args))
  (cl:unless (cl:typep m 'GetRendezvousPoint-request)
    (roslisp-msg-protocol:msg-deprecation-warning "using old message class name mk_msgs-srv:<GetRendezvousPoint-request> is deprecated: use mk_msgs-srv:GetRendezvousPoint-request instead.")))

(cl:ensure-generic-function 'currLatitude-val :lambda-list '(m))
(cl:defmethod currLatitude-val ((m <GetRendezvousPoint-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:currLatitude-val is deprecated.  Use mk_msgs-srv:currLatitude instead.")
  (currLatitude m))

(cl:ensure-generic-function 'currLongitude-val :lambda-list '(m))
(cl:defmethod currLongitude-val ((m <GetRendezvousPoint-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:currLongitude-val is deprecated.  Use mk_msgs-srv:currLongitude instead.")
  (currLongitude m))
(cl:defmethod roslisp-msg-protocol:serialize ((msg <GetRendezvousPoint-request>) ostream)
  "Serializes a message object of type '<GetRendezvousPoint-request>"
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'currLatitude))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'currLongitude))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
)
(cl:defmethod roslisp-msg-protocol:deserialize ((msg <GetRendezvousPoint-request>) istream)
  "Deserializes a message object of type '<GetRendezvousPoint-request>"
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'currLatitude) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'currLongitude) (roslisp-utils:decode-single-float-bits bits)))
  msg
)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql '<GetRendezvousPoint-request>)))
  "Returns string type for a service object of type '<GetRendezvousPoint-request>"
  "mk_msgs/GetRendezvousPointRequest")
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'GetRendezvousPoint-request)))
  "Returns string type for a service object of type 'GetRendezvousPoint-request"
  "mk_msgs/GetRendezvousPointRequest")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql '<GetRendezvousPoint-request>)))
  "Returns md5sum for a message object of type '<GetRendezvousPoint-request>"
  "60fb96562d024b790f9c37a37ae940f9")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql 'GetRendezvousPoint-request)))
  "Returns md5sum for a message object of type 'GetRendezvousPoint-request"
  "60fb96562d024b790f9c37a37ae940f9")
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql '<GetRendezvousPoint-request>)))
  "Returns full string definition for message of type '<GetRendezvousPoint-request>"
  (cl:format cl:nil "~%~%~%float32 currLatitude~%float32 currLongitude~%~%~%"))
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql 'GetRendezvousPoint-request)))
  "Returns full string definition for message of type 'GetRendezvousPoint-request"
  (cl:format cl:nil "~%~%~%float32 currLatitude~%float32 currLongitude~%~%~%"))
(cl:defmethod roslisp-msg-protocol:serialization-length ((msg <GetRendezvousPoint-request>))
  (cl:+ 0
     4
     4
))
(cl:defmethod roslisp-msg-protocol:ros-message-to-list ((msg <GetRendezvousPoint-request>))
  "Converts a ROS message object to a list"
  (cl:list 'GetRendezvousPoint-request
    (cl:cons ':currLatitude (currLatitude msg))
    (cl:cons ':currLongitude (currLongitude msg))
))
;//! \htmlinclude GetRendezvousPoint-response.msg.html

(cl:defclass <GetRendezvousPoint-response> (roslisp-msg-protocol:ros-message)
  ((ack
    :reader ack
    :initarg :ack
    :type cl:fixnum
    :initform 0)
   (rendezLatitude
    :reader rendezLatitude
    :initarg :rendezLatitude
    :type cl:float
    :initform 0.0)
   (rendezLongitude
    :reader rendezLongitude
    :initarg :rendezLongitude
    :type cl:float
    :initform 0.0))
)

(cl:defclass GetRendezvousPoint-response (<GetRendezvousPoint-response>)
  ())

(cl:defmethod cl:initialize-instance :after ((m <GetRendezvousPoint-response>) cl:&rest args)
  (cl:declare (cl:ignorable args))
  (cl:unless (cl:typep m 'GetRendezvousPoint-response)
    (roslisp-msg-protocol:msg-deprecation-warning "using old message class name mk_msgs-srv:<GetRendezvousPoint-response> is deprecated: use mk_msgs-srv:GetRendezvousPoint-response instead.")))

(cl:ensure-generic-function 'ack-val :lambda-list '(m))
(cl:defmethod ack-val ((m <GetRendezvousPoint-response>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:ack-val is deprecated.  Use mk_msgs-srv:ack instead.")
  (ack m))

(cl:ensure-generic-function 'rendezLatitude-val :lambda-list '(m))
(cl:defmethod rendezLatitude-val ((m <GetRendezvousPoint-response>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:rendezLatitude-val is deprecated.  Use mk_msgs-srv:rendezLatitude instead.")
  (rendezLatitude m))

(cl:ensure-generic-function 'rendezLongitude-val :lambda-list '(m))
(cl:defmethod rendezLongitude-val ((m <GetRendezvousPoint-response>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:rendezLongitude-val is deprecated.  Use mk_msgs-srv:rendezLongitude instead.")
  (rendezLongitude m))
(cl:defmethod roslisp-msg-protocol:serialize ((msg <GetRendezvousPoint-response>) ostream)
  "Serializes a message object of type '<GetRendezvousPoint-response>"
  (cl:let* ((signed (cl:slot-value msg 'ack)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 65536) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    )
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'rendezLatitude))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
  (cl:let ((bits (roslisp-utils:encode-single-float-bits (cl:slot-value msg 'rendezLongitude))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream))
)
(cl:defmethod roslisp-msg-protocol:deserialize ((msg <GetRendezvousPoint-response>) istream)
  "Deserializes a message object of type '<GetRendezvousPoint-response>"
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'ack) (cl:if (cl:< unsigned 32768) unsigned (cl:- unsigned 65536))))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'rendezLatitude) (roslisp-utils:decode-single-float-bits bits)))
    (cl:let ((bits 0))
      (cl:setf (cl:ldb (cl:byte 8 0) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) bits) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) bits) (cl:read-byte istream))
    (cl:setf (cl:slot-value msg 'rendezLongitude) (roslisp-utils:decode-single-float-bits bits)))
  msg
)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql '<GetRendezvousPoint-response>)))
  "Returns string type for a service object of type '<GetRendezvousPoint-response>"
  "mk_msgs/GetRendezvousPointResponse")
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'GetRendezvousPoint-response)))
  "Returns string type for a service object of type 'GetRendezvousPoint-response"
  "mk_msgs/GetRendezvousPointResponse")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql '<GetRendezvousPoint-response>)))
  "Returns md5sum for a message object of type '<GetRendezvousPoint-response>"
  "60fb96562d024b790f9c37a37ae940f9")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql 'GetRendezvousPoint-response)))
  "Returns md5sum for a message object of type 'GetRendezvousPoint-response"
  "60fb96562d024b790f9c37a37ae940f9")
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql '<GetRendezvousPoint-response>)))
  "Returns full string definition for message of type '<GetRendezvousPoint-response>"
  (cl:format cl:nil "~%~%~%int16 ack~%float32 rendezLatitude~%float32 rendezLongitude~%~%~%~%"))
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql 'GetRendezvousPoint-response)))
  "Returns full string definition for message of type 'GetRendezvousPoint-response"
  (cl:format cl:nil "~%~%~%int16 ack~%float32 rendezLatitude~%float32 rendezLongitude~%~%~%~%"))
(cl:defmethod roslisp-msg-protocol:serialization-length ((msg <GetRendezvousPoint-response>))
  (cl:+ 0
     2
     4
     4
))
(cl:defmethod roslisp-msg-protocol:ros-message-to-list ((msg <GetRendezvousPoint-response>))
  "Converts a ROS message object to a list"
  (cl:list 'GetRendezvousPoint-response
    (cl:cons ':ack (ack msg))
    (cl:cons ':rendezLatitude (rendezLatitude msg))
    (cl:cons ':rendezLongitude (rendezLongitude msg))
))
(cl:defmethod roslisp-msg-protocol:service-request-type ((msg (cl:eql 'GetRendezvousPoint)))
  'GetRendezvousPoint-request)
(cl:defmethod roslisp-msg-protocol:service-response-type ((msg (cl:eql 'GetRendezvousPoint)))
  'GetRendezvousPoint-response)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'GetRendezvousPoint)))
  "Returns string type for a service object of type '<GetRendezvousPoint>"
  "mk_msgs/GetRendezvousPoint")