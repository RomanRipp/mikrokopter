; Auto-generated. Do not edit!


(cl:in-package mk_msgs-srv)


;//! \htmlinclude SetSerialPoti-request.msg.html

(cl:defclass <SetSerialPoti-request> (roslisp-msg-protocol:ros-message)
  ((poti1
    :reader poti1
    :initarg :poti1
    :type cl:integer
    :initform 0)
   (poti2
    :reader poti2
    :initarg :poti2
    :type cl:integer
    :initform 0)
   (poti3
    :reader poti3
    :initarg :poti3
    :type cl:integer
    :initform 0)
   (poti4
    :reader poti4
    :initarg :poti4
    :type cl:integer
    :initform 0))
)

(cl:defclass SetSerialPoti-request (<SetSerialPoti-request>)
  ())

(cl:defmethod cl:initialize-instance :after ((m <SetSerialPoti-request>) cl:&rest args)
  (cl:declare (cl:ignorable args))
  (cl:unless (cl:typep m 'SetSerialPoti-request)
    (roslisp-msg-protocol:msg-deprecation-warning "using old message class name mk_msgs-srv:<SetSerialPoti-request> is deprecated: use mk_msgs-srv:SetSerialPoti-request instead.")))

(cl:ensure-generic-function 'poti1-val :lambda-list '(m))
(cl:defmethod poti1-val ((m <SetSerialPoti-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:poti1-val is deprecated.  Use mk_msgs-srv:poti1 instead.")
  (poti1 m))

(cl:ensure-generic-function 'poti2-val :lambda-list '(m))
(cl:defmethod poti2-val ((m <SetSerialPoti-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:poti2-val is deprecated.  Use mk_msgs-srv:poti2 instead.")
  (poti2 m))

(cl:ensure-generic-function 'poti3-val :lambda-list '(m))
(cl:defmethod poti3-val ((m <SetSerialPoti-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:poti3-val is deprecated.  Use mk_msgs-srv:poti3 instead.")
  (poti3 m))

(cl:ensure-generic-function 'poti4-val :lambda-list '(m))
(cl:defmethod poti4-val ((m <SetSerialPoti-request>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:poti4-val is deprecated.  Use mk_msgs-srv:poti4 instead.")
  (poti4 m))
(cl:defmethod roslisp-msg-protocol:serialize ((msg <SetSerialPoti-request>) ostream)
  "Serializes a message object of type '<SetSerialPoti-request>"
  (cl:let* ((signed (cl:slot-value msg 'poti1)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'poti2)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'poti3)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
  (cl:let* ((signed (cl:slot-value msg 'poti4)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 4294967296) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) unsigned) ostream)
    )
)
(cl:defmethod roslisp-msg-protocol:deserialize ((msg <SetSerialPoti-request>) istream)
  "Deserializes a message object of type '<SetSerialPoti-request>"
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'poti1) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'poti2) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'poti3) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 16) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 24) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'poti4) (cl:if (cl:< unsigned 2147483648) unsigned (cl:- unsigned 4294967296))))
  msg
)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql '<SetSerialPoti-request>)))
  "Returns string type for a service object of type '<SetSerialPoti-request>"
  "mk_msgs/SetSerialPotiRequest")
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'SetSerialPoti-request)))
  "Returns string type for a service object of type 'SetSerialPoti-request"
  "mk_msgs/SetSerialPotiRequest")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql '<SetSerialPoti-request>)))
  "Returns md5sum for a message object of type '<SetSerialPoti-request>"
  "d1224392e458021cfd856dc94ab7891a")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql 'SetSerialPoti-request)))
  "Returns md5sum for a message object of type 'SetSerialPoti-request"
  "d1224392e458021cfd856dc94ab7891a")
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql '<SetSerialPoti-request>)))
  "Returns full string definition for message of type '<SetSerialPoti-request>"
  (cl:format cl:nil "~%~%~%int32 poti1~%int32 poti2~%int32 poti3~%int32 poti4~%~%~%~%"))
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql 'SetSerialPoti-request)))
  "Returns full string definition for message of type 'SetSerialPoti-request"
  (cl:format cl:nil "~%~%~%int32 poti1~%int32 poti2~%int32 poti3~%int32 poti4~%~%~%~%"))
(cl:defmethod roslisp-msg-protocol:serialization-length ((msg <SetSerialPoti-request>))
  (cl:+ 0
     4
     4
     4
     4
))
(cl:defmethod roslisp-msg-protocol:ros-message-to-list ((msg <SetSerialPoti-request>))
  "Converts a ROS message object to a list"
  (cl:list 'SetSerialPoti-request
    (cl:cons ':poti1 (poti1 msg))
    (cl:cons ':poti2 (poti2 msg))
    (cl:cons ':poti3 (poti3 msg))
    (cl:cons ':poti4 (poti4 msg))
))
;//! \htmlinclude SetSerialPoti-response.msg.html

(cl:defclass <SetSerialPoti-response> (roslisp-msg-protocol:ros-message)
  ((ack
    :reader ack
    :initarg :ack
    :type cl:fixnum
    :initform 0))
)

(cl:defclass SetSerialPoti-response (<SetSerialPoti-response>)
  ())

(cl:defmethod cl:initialize-instance :after ((m <SetSerialPoti-response>) cl:&rest args)
  (cl:declare (cl:ignorable args))
  (cl:unless (cl:typep m 'SetSerialPoti-response)
    (roslisp-msg-protocol:msg-deprecation-warning "using old message class name mk_msgs-srv:<SetSerialPoti-response> is deprecated: use mk_msgs-srv:SetSerialPoti-response instead.")))

(cl:ensure-generic-function 'ack-val :lambda-list '(m))
(cl:defmethod ack-val ((m <SetSerialPoti-response>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_msgs-srv:ack-val is deprecated.  Use mk_msgs-srv:ack instead.")
  (ack m))
(cl:defmethod roslisp-msg-protocol:serialize ((msg <SetSerialPoti-response>) ostream)
  "Serializes a message object of type '<SetSerialPoti-response>"
  (cl:let* ((signed (cl:slot-value msg 'ack)) (unsigned (cl:if (cl:< signed 0) (cl:+ signed 65536) signed)))
    (cl:write-byte (cl:ldb (cl:byte 8 0) unsigned) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) unsigned) ostream)
    )
)
(cl:defmethod roslisp-msg-protocol:deserialize ((msg <SetSerialPoti-response>) istream)
  "Deserializes a message object of type '<SetSerialPoti-response>"
    (cl:let ((unsigned 0))
      (cl:setf (cl:ldb (cl:byte 8 0) unsigned) (cl:read-byte istream))
      (cl:setf (cl:ldb (cl:byte 8 8) unsigned) (cl:read-byte istream))
      (cl:setf (cl:slot-value msg 'ack) (cl:if (cl:< unsigned 32768) unsigned (cl:- unsigned 65536))))
  msg
)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql '<SetSerialPoti-response>)))
  "Returns string type for a service object of type '<SetSerialPoti-response>"
  "mk_msgs/SetSerialPotiResponse")
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'SetSerialPoti-response)))
  "Returns string type for a service object of type 'SetSerialPoti-response"
  "mk_msgs/SetSerialPotiResponse")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql '<SetSerialPoti-response>)))
  "Returns md5sum for a message object of type '<SetSerialPoti-response>"
  "d1224392e458021cfd856dc94ab7891a")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql 'SetSerialPoti-response)))
  "Returns md5sum for a message object of type 'SetSerialPoti-response"
  "d1224392e458021cfd856dc94ab7891a")
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql '<SetSerialPoti-response>)))
  "Returns full string definition for message of type '<SetSerialPoti-response>"
  (cl:format cl:nil "~%~%~%int16 ack~%~%~%~%"))
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql 'SetSerialPoti-response)))
  "Returns full string definition for message of type 'SetSerialPoti-response"
  (cl:format cl:nil "~%~%~%int16 ack~%~%~%~%"))
(cl:defmethod roslisp-msg-protocol:serialization-length ((msg <SetSerialPoti-response>))
  (cl:+ 0
     2
))
(cl:defmethod roslisp-msg-protocol:ros-message-to-list ((msg <SetSerialPoti-response>))
  "Converts a ROS message object to a list"
  (cl:list 'SetSerialPoti-response
    (cl:cons ':ack (ack msg))
))
(cl:defmethod roslisp-msg-protocol:service-request-type ((msg (cl:eql 'SetSerialPoti)))
  'SetSerialPoti-request)
(cl:defmethod roslisp-msg-protocol:service-response-type ((msg (cl:eql 'SetSerialPoti)))
  'SetSerialPoti-response)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'SetSerialPoti)))
  "Returns string type for a service object of type '<SetSerialPoti>"
  "mk_msgs/SetSerialPoti")