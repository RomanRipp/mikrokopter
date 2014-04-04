; Auto-generated. Do not edit!


(cl:in-package mk_msgs-srv)


;//! \htmlinclude listwp-request.msg.html

(cl:defclass <listwp-request> (roslisp-msg-protocol:ros-message)
  ()
)

(cl:defclass listwp-request (<listwp-request>)
  ())

(cl:defmethod cl:initialize-instance :after ((m <listwp-request>) cl:&rest args)
  (cl:declare (cl:ignorable args))
  (cl:unless (cl:typep m 'listwp-request)
    (roslisp-msg-protocol:msg-deprecation-warning "using old message class name mk_msgs-srv:<listwp-request> is deprecated: use mk_msgs-srv:listwp-request instead.")))
(cl:defmethod roslisp-msg-protocol:serialize ((msg <listwp-request>) ostream)
  "Serializes a message object of type '<listwp-request>"
)
(cl:defmethod roslisp-msg-protocol:deserialize ((msg <listwp-request>) istream)
  "Deserializes a message object of type '<listwp-request>"
  msg
)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql '<listwp-request>)))
  "Returns string type for a service object of type '<listwp-request>"
  "mk_msgs/listwpRequest")
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'listwp-request)))
  "Returns string type for a service object of type 'listwp-request"
  "mk_msgs/listwpRequest")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql '<listwp-request>)))
  "Returns md5sum for a message object of type '<listwp-request>"
  "efd27fb04e2a6ce3c9ff1f47eb32e7bb")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql 'listwp-request)))
  "Returns md5sum for a message object of type 'listwp-request"
  "efd27fb04e2a6ce3c9ff1f47eb32e7bb")
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql '<listwp-request>)))
  "Returns full string definition for message of type '<listwp-request>"
  (cl:format cl:nil "~%~%~%~%~%"))
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql 'listwp-request)))
  "Returns full string definition for message of type 'listwp-request"
  (cl:format cl:nil "~%~%~%~%~%"))
(cl:defmethod roslisp-msg-protocol:serialization-length ((msg <listwp-request>))
  (cl:+ 0
))
(cl:defmethod roslisp-msg-protocol:ros-message-to-list ((msg <listwp-request>))
  "Converts a ROS message object to a list"
  (cl:list 'listwp-request
))
;//! \htmlinclude listwp-response.msg.html

(cl:defclass <listwp-response> (roslisp-msg-protocol:ros-message)
  ((ack
    :reader ack
    :initarg :ack
    :type cl:fixnum
    :initform 0))
)

(cl:defclass listwp-response (<listwp-response>)
  ())

(cl:defmethod cl:initialize-instance :after ((m <listwp-response>) cl:&rest args)
  (cl:declare (cl:ignorable args))
  (cl:unless (cl:typep m 'listwp-response)
    (roslisp-msg-protocol:msg-deprecation-warning "using old message class name mk_msgs-srv:<listwp-response> is deprecated: use mk_msgs-srv:listwp-response instead.")))

(cl:ensure-generic-function 'ack-val :lambda-list '(m))
(cl:defmethod ack-val ((m <listwp-response>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:ack-val is deprecated.  Use mk_msgs-srv:ack instead.")
  (ack m))
(cl:defmethod roslisp-msg-protocol:serialize ((msg <listwp-response>) ostream)
  "Serializes a message object of type '<listwp-response>"
  (cl:let* ((signed (cl:slot-value msg 'ack)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 65536) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    )
)
(cl:defmethod roslisp-msg-protocol:deserialize ((msg <listwp-response>) istream)
  "Deserializes a message object of type '<listwp-response>"
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'ack) (cl:if (cl:< unsigned 32768) unsigned (cl:- unsigned 65536))))
  msg
)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql '<listwp-response>)))
  "Returns string type for a service object of type '<listwp-response>"
  "mk_msgs/listwpResponse")
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'listwp-response)))
  "Returns string type for a service object of type 'listwp-response"
  "mk_msgs/listwpResponse")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql '<listwp-response>)))
  "Returns md5sum for a message object of type '<listwp-response>"
  "efd27fb04e2a6ce3c9ff1f47eb32e7bb")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql 'listwp-response)))
  "Returns md5sum for a message object of type 'listwp-response"
  "efd27fb04e2a6ce3c9ff1f47eb32e7bb")
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql '<listwp-response>)))
  "Returns full string definition for message of type '<listwp-response>"
  (cl:format cl:nil "~%~%~%int16 ack~%~%~%~%"))
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql 'listwp-response)))
  "Returns full string definition for message of type 'listwp-response"
  (cl:format cl:nil "~%~%~%int16 ack~%~%~%~%"))
(cl:defmethod roslisp-msg-protocol:serialization-length ((msg <listwp-response>))
  (cl:+ 0
     2
))
(cl:defmethod roslisp-msg-protocol:ros-message-to-list ((msg <listwp-response>))
  "Converts a ROS message object to a list"
  (cl:list 'listwp-response
    (cl:cons ':ack (ack msg))
))
(cl:defmethod roslisp-msg-protocol:service-request-type ((msg (cl:eql 'listwp)))
  'listwp-request)
(cl:defmethod roslisp-msg-protocol:service-response-type ((msg (cl:eql 'listwp)))
  'listwp-response)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'listwp)))
  "Returns string type for a service object of type '<listwp>"
  "mk_msgs/listwp")