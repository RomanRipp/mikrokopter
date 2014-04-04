; Auto-generated. Do not edit!


(cl:in-package mk_msgs-srv)


;//! \htmlinclude SetNavigation-request.msg.html

(cl:defclass <SetNavigation-request> (roslisp-msg-protocol:ros-message)
  ((numofwps
    :reader numofwps
    :initarg :numofwps
    :type cl:integer
    :initform 0))
)

(cl:defclass SetNavigation-request (<SetNavigation-request>)
  ())

(cl:defmethod cl:initialize-instance :after ((m <SetNavigation-request>) cl:&rest args)
  (cl:declare (cl:ignorable args))
  (cl:unless (cl:typep m 'SetNavigation-request)
    (roslisp-msg-protocol:msg-deprecation-warning "using old message class name mk_msgs-srv:<SetNavigation-request> is deprecated: use mk_msgs-srv:SetNavigation-request instead.")))

(cl:ensure-generic-function 'numofwps-val :lambda-list '(m))
(cl:defmethod numofwps-val ((m <SetNavigation-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:numofwps-val is deprecated.  Use mk_msgs-srv:numofwps instead.")
  (numofwps m))
(cl:defmethod roslisp-msg-protocol:serialize ((msg <SetNavigation-request>) ostream)
  "Serializes a message object of type '<SetNavigation-request>"
  (cl:let* ((signed (cl:slot-value msg 'numofwps)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
)
(cl:defmethod roslisp-msg-protocol:deserialize ((msg <SetNavigation-request>) istream)
  "Deserializes a message object of type '<SetNavigation-request>"
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'numofwps) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
  msg
)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql '<SetNavigation-request>)))
  "Returns string type for a service object of type '<SetNavigation-request>"
  "mk_msgs/SetNavigationRequest")
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'SetNavigation-request)))
  "Returns string type for a service object of type 'SetNavigation-request"
  "mk_msgs/SetNavigationRequest")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql '<SetNavigation-request>)))
  "Returns md5sum for a message object of type '<SetNavigation-request>"
  "bebf409cbd73a6674c483f81f26e37d8")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql 'SetNavigation-request)))
  "Returns md5sum for a message object of type 'SetNavigation-request"
  "bebf409cbd73a6674c483f81f26e37d8")
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql '<SetNavigation-request>)))
  "Returns full string definition for message of type '<SetNavigation-request>"
  (cl:format cl:nil "~%~%int32 numofwps~%~%~%"))
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql 'SetNavigation-request)))
  "Returns full string definition for message of type 'SetNavigation-request"
  (cl:format cl:nil "~%~%int32 numofwps~%~%~%"))
(cl:defmethod roslisp-msg-protocol:serialization-length ((msg <SetNavigation-request>))
  (cl:+ 0
     4
))
(cl:defmethod roslisp-msg-protocol:ros-message-to-list ((msg <SetNavigation-request>))
  "Converts a ROS message object to a list"
  (cl:list 'SetNavigation-request
    (cl:cons ':numofwps (numofwps msg))
))
;//! \htmlinclude SetNavigation-response.msg.html

(cl:defclass <SetNavigation-response> (roslisp-msg-protocol:ros-message)
  ((ack
    :reader ack
    :initarg :ack
    :type cl:fixnum
    :initform 0))
)

(cl:defclass SetNavigation-response (<SetNavigation-response>)
  ())

(cl:defmethod cl:initialize-instance :after ((m <SetNavigation-response>) cl:&rest args)
  (cl:declare (cl:ignorable args))
  (cl:unless (cl:typep m 'SetNavigation-response)
    (roslisp-msg-protocol:msg-deprecation-warning "using old message class name mk_msgs-srv:<SetNavigation-response> is deprecated: use mk_msgs-srv:SetNavigation-response instead.")))

(cl:ensure-generic-function 'ack-val :lambda-list '(m))
(cl:defmethod ack-val ((m <SetNavigation-response>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:ack-val is deprecated.  Use mk_msgs-srv:ack instead.")
  (ack m))
(cl:defmethod roslisp-msg-protocol:serialize ((msg <SetNavigation-response>) ostream)
  "Serializes a message object of type '<SetNavigation-response>"
  (cl:let* ((signed (cl:slot-value msg 'ack)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 65536) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    )
)
(cl:defmethod roslisp-msg-protocol:deserialize ((msg <SetNavigation-response>) istream)
  "Deserializes a message object of type '<SetNavigation-response>"
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'ack) (cl:if (cl:< unsigned 32768) unsigned (cl:- unsigned 65536))))
  msg
)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql '<SetNavigation-response>)))
  "Returns string type for a service object of type '<SetNavigation-response>"
  "mk_msgs/SetNavigationResponse")
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'SetNavigation-response)))
  "Returns string type for a service object of type 'SetNavigation-response"
  "mk_msgs/SetNavigationResponse")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql '<SetNavigation-response>)))
  "Returns md5sum for a message object of type '<SetNavigation-response>"
  "bebf409cbd73a6674c483f81f26e37d8")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql 'SetNavigation-response)))
  "Returns md5sum for a message object of type 'SetNavigation-response"
  "bebf409cbd73a6674c483f81f26e37d8")
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql '<SetNavigation-response>)))
  "Returns full string definition for message of type '<SetNavigation-response>"
  (cl:format cl:nil "~%~%~%int16 ack~%~%~%~%"))
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql 'SetNavigation-response)))
  "Returns full string definition for message of type 'SetNavigation-response"
  (cl:format cl:nil "~%~%~%int16 ack~%~%~%~%"))
(cl:defmethod roslisp-msg-protocol:serialization-length ((msg <SetNavigation-response>))
  (cl:+ 0
     2
))
(cl:defmethod roslisp-msg-protocol:ros-message-to-list ((msg <SetNavigation-response>))
  "Converts a ROS message object to a list"
  (cl:list 'SetNavigation-response
    (cl:cons ':ack (ack msg))
))
(cl:defmethod roslisp-msg-protocol:service-request-type ((msg (cl:eql 'SetNavigation)))
  'SetNavigation-request)
(cl:defmethod roslisp-msg-protocol:service-response-type ((msg (cl:eql 'SetNavigation)))
  'SetNavigation-response)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'SetNavigation)))
  "Returns string type for a service object of type '<SetNavigation>"
  "mk_msgs/SetNavigation")