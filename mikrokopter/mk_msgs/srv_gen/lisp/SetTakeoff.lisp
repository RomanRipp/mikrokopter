; Auto-generated. Do not edit!


(cl:in-package mk_msgs-srv)


;//! \htmlinclude SetTakeoff-request.msg.html

(cl:defclass <SetTakeoff-request> (roslisp-msg-protocol:ros-message)
  ((height
    :reader height
    :initarg :height
    :type cl:float
    :initform 0.0))
)

(cl:defclass SetTakeoff-request (<SetTakeoff-request>)
  ())

(cl:defmethod cl:initialize-instance :after ((m <SetTakeoff-request>) cl:&rest args)
  (cl:declare (cl:ignorable args))
  (cl:unless (cl:typep m 'SetTakeoff-request)
    (roslisp-msg-protocol:msg-deprecation-warning "using old message class name mk_msgs-srv:<SetTakeoff-request> is deprecated: use mk_msgs-srv:SetTakeoff-request instead.")))

(cl:ensure-generic-function 'height-val :lambda-list '(m))
(cl:defmethod height-val ((m <SetTakeoff-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:height-val is deprecated.  Use mk_msgs-srv:height instead.")
  (height m))
(cl:defmethod roslisp-msg-protocol:serialize ((msg <SetTakeoff-request>) ostream)
  "Serializes a message object of type '<SetTakeoff-request>"
  (cl:let ((bits (roslisp-utils:encode-double-float-bits (cl:slot-value msg 'height))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 32) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 40) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 48) bits) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 56) bits) ostream))
)
(cl:defmethod roslisp-msg-protocol:deserialize ((msg <SetTakeoff-request>) istream)
  "Deserializes a message object of type '<SetTakeoff-request>"
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
  msg
)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql '<SetTakeoff-request>)))
  "Returns string type for a service object of type '<SetTakeoff-request>"
  "mk_msgs/SetTakeoffRequest")
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'SetTakeoff-request)))
  "Returns string type for a service object of type 'SetTakeoff-request"
  "mk_msgs/SetTakeoffRequest")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql '<SetTakeoff-request>)))
  "Returns md5sum for a message object of type '<SetTakeoff-request>"
  "b563bf2328afa304f456e28b10cefd41")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql 'SetTakeoff-request)))
  "Returns md5sum for a message object of type 'SetTakeoff-request"
  "b563bf2328afa304f456e28b10cefd41")
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql '<SetTakeoff-request>)))
  "Returns full string definition for message of type '<SetTakeoff-request>"
  (cl:format cl:nil "~%~%float64 height~%~%~%"))
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql 'SetTakeoff-request)))
  "Returns full string definition for message of type 'SetTakeoff-request"
  (cl:format cl:nil "~%~%float64 height~%~%~%"))
(cl:defmethod roslisp-msg-protocol:serialization-length ((msg <SetTakeoff-request>))
  (cl:+ 0
     8
))
(cl:defmethod roslisp-msg-protocol:ros-message-to-list ((msg <SetTakeoff-request>))
  "Converts a ROS message object to a list"
  (cl:list 'SetTakeoff-request
    (cl:cons ':height (height msg))
))
;//! \htmlinclude SetTakeoff-response.msg.html

(cl:defclass <SetTakeoff-response> (roslisp-msg-protocol:ros-message)
  ((ack
    :reader ack
    :initarg :ack
    :type cl:fixnum
    :initform 0))
)

(cl:defclass SetTakeoff-response (<SetTakeoff-response>)
  ())

(cl:defmethod cl:initialize-instance :after ((m <SetTakeoff-response>) cl:&rest args)
  (cl:declare (cl:ignorable args))
  (cl:unless (cl:typep m 'SetTakeoff-response)
    (roslisp-msg-protocol:msg-deprecation-warning "using old message class name mk_msgs-srv:<SetTakeoff-response> is deprecated: use mk_msgs-srv:SetTakeoff-response instead.")))

(cl:ensure-generic-function 'ack-val :lambda-list '(m))
(cl:defmethod ack-val ((m <SetTakeoff-response>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:ack-val is deprecated.  Use mk_msgs-srv:ack instead.")
  (ack m))
(cl:defmethod roslisp-msg-protocol:serialize ((msg <SetTakeoff-response>) ostream)
  "Serializes a message object of type '<SetTakeoff-response>"
  (cl:let* ((signed (cl:slot-value msg 'ack)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 65536) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    )
)
(cl:defmethod roslisp-msg-protocol:deserialize ((msg <SetTakeoff-response>) istream)
  "Deserializes a message object of type '<SetTakeoff-response>"
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'ack) (cl:if (cl:< unsigned 32768) unsigned (cl:- unsigned 65536))))
  msg
)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql '<SetTakeoff-response>)))
  "Returns string type for a service object of type '<SetTakeoff-response>"
  "mk_msgs/SetTakeoffResponse")
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'SetTakeoff-response)))
  "Returns string type for a service object of type 'SetTakeoff-response"
  "mk_msgs/SetTakeoffResponse")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql '<SetTakeoff-response>)))
  "Returns md5sum for a message object of type '<SetTakeoff-response>"
  "b563bf2328afa304f456e28b10cefd41")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql 'SetTakeoff-response)))
  "Returns md5sum for a message object of type 'SetTakeoff-response"
  "b563bf2328afa304f456e28b10cefd41")
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql '<SetTakeoff-response>)))
  "Returns full string definition for message of type '<SetTakeoff-response>"
  (cl:format cl:nil "~%~%~%int16 ack~%~%~%~%"))
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql 'SetTakeoff-response)))
  "Returns full string definition for message of type 'SetTakeoff-response"
  (cl:format cl:nil "~%~%~%int16 ack~%~%~%~%"))
(cl:defmethod roslisp-msg-protocol:serialization-length ((msg <SetTakeoff-response>))
  (cl:+ 0
     2
))
(cl:defmethod roslisp-msg-protocol:ros-message-to-list ((msg <SetTakeoff-response>))
  "Converts a ROS message object to a list"
  (cl:list 'SetTakeoff-response
    (cl:cons ':ack (ack msg))
))
(cl:defmethod roslisp-msg-protocol:service-request-type ((msg (cl:eql 'SetTakeoff)))
  'SetTakeoff-request)
(cl:defmethod roslisp-msg-protocol:service-response-type ((msg (cl:eql 'SetTakeoff)))
  'SetTakeoff-response)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'SetTakeoff)))
  "Returns string type for a service object of type '<SetTakeoff>"
  "mk_msgs/SetTakeoff")