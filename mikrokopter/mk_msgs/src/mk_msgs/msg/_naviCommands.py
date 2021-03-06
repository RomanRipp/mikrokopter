"""autogenerated by genpy from mk_msgs/naviCommands.msg. Do not edit."""
import sys
python3 = True if sys.hexversion > 0x03000000 else False
import genpy
import struct

import std_msgs.msg

class naviCommands(genpy.Message):
  _md5sum = "86916e842c23439ce3ae34b26ef90d4f"
  _type = "mk_msgs/naviCommands"
  _has_header = True #flag to mark the presence of a Header object
  _full_text = """Header header

# between -pi to +pi
float32 longitude
float32 latitude
float32 altitude
float32 heading
int32 toleranceRadius
int32 holdTime
int32 index
int32 eventChannel
int32 altitudeRate
float32 speed
int32 status
int32 eventFlag
int32 type
int32 camAngle


================================================================================
MSG: std_msgs/Header
# Standard metadata for higher-level stamped data types.
# This is generally used to communicate timestamped data 
# in a particular coordinate frame.
# 
# sequence ID: consecutively increasing ID 
uint32 seq
#Two-integer timestamp that is expressed as:
# * stamp.secs: seconds (stamp_secs) since epoch
# * stamp.nsecs: nanoseconds since stamp_secs
# time-handling sugar is provided by the client library
time stamp
#Frame this data is associated with
# 0: no frame
# 1: global frame
string frame_id

"""
  __slots__ = ['header','longitude','latitude','altitude','heading','toleranceRadius','holdTime','index','eventChannel','altitudeRate','speed','status','eventFlag','type','camAngle']
  _slot_types = ['std_msgs/Header','float32','float32','float32','float32','int32','int32','int32','int32','int32','float32','int32','int32','int32','int32']

  def __init__(self, *args, **kwds):
    """
    Constructor. Any message fields that are implicitly/explicitly
    set to None will be assigned a default value. The recommend
    use is keyword arguments as this is more robust to future message
    changes.  You cannot mix in-order arguments and keyword arguments.

    The available fields are:
       header,longitude,latitude,altitude,heading,toleranceRadius,holdTime,index,eventChannel,altitudeRate,speed,status,eventFlag,type,camAngle

    :param args: complete set of field values, in .msg order
    :param kwds: use keyword arguments corresponding to message field names
    to set specific fields.
    """
    if args or kwds:
      super(naviCommands, self).__init__(*args, **kwds)
      #message fields cannot be None, assign default values for those that are
      if self.header is None:
        self.header = std_msgs.msg.Header()
      if self.longitude is None:
        self.longitude = 0.
      if self.latitude is None:
        self.latitude = 0.
      if self.altitude is None:
        self.altitude = 0.
      if self.heading is None:
        self.heading = 0.
      if self.toleranceRadius is None:
        self.toleranceRadius = 0
      if self.holdTime is None:
        self.holdTime = 0
      if self.index is None:
        self.index = 0
      if self.eventChannel is None:
        self.eventChannel = 0
      if self.altitudeRate is None:
        self.altitudeRate = 0
      if self.speed is None:
        self.speed = 0.
      if self.status is None:
        self.status = 0
      if self.eventFlag is None:
        self.eventFlag = 0
      if self.type is None:
        self.type = 0
      if self.camAngle is None:
        self.camAngle = 0
    else:
      self.header = std_msgs.msg.Header()
      self.longitude = 0.
      self.latitude = 0.
      self.altitude = 0.
      self.heading = 0.
      self.toleranceRadius = 0
      self.holdTime = 0
      self.index = 0
      self.eventChannel = 0
      self.altitudeRate = 0
      self.speed = 0.
      self.status = 0
      self.eventFlag = 0
      self.type = 0
      self.camAngle = 0

  def _get_types(self):
    """
    internal API method
    """
    return self._slot_types

  def serialize(self, buff):
    """
    serialize message into buffer
    :param buff: buffer, ``StringIO``
    """
    try:
      _x = self
      buff.write(_struct_3I.pack(_x.header.seq, _x.header.stamp.secs, _x.header.stamp.nsecs))
      _x = self.header.frame_id
      length = len(_x)
      if python3 or type(_x) == unicode:
        _x = _x.encode('utf-8')
        length = len(_x)
      buff.write(struct.pack('<I%ss'%length, length, _x))
      _x = self
      buff.write(_struct_4f5if4i.pack(_x.longitude, _x.latitude, _x.altitude, _x.heading, _x.toleranceRadius, _x.holdTime, _x.index, _x.eventChannel, _x.altitudeRate, _x.speed, _x.status, _x.eventFlag, _x.type, _x.camAngle))
    except struct.error as se: self._check_types(se)
    except TypeError as te: self._check_types(te)

  def deserialize(self, str):
    """
    unpack serialized message in str into this message instance
    :param str: byte array of serialized message, ``str``
    """
    try:
      if self.header is None:
        self.header = std_msgs.msg.Header()
      end = 0
      _x = self
      start = end
      end += 12
      (_x.header.seq, _x.header.stamp.secs, _x.header.stamp.nsecs,) = _struct_3I.unpack(str[start:end])
      start = end
      end += 4
      (length,) = _struct_I.unpack(str[start:end])
      start = end
      end += length
      if python3:
        self.header.frame_id = str[start:end].decode('utf-8')
      else:
        self.header.frame_id = str[start:end]
      _x = self
      start = end
      end += 56
      (_x.longitude, _x.latitude, _x.altitude, _x.heading, _x.toleranceRadius, _x.holdTime, _x.index, _x.eventChannel, _x.altitudeRate, _x.speed, _x.status, _x.eventFlag, _x.type, _x.camAngle,) = _struct_4f5if4i.unpack(str[start:end])
      return self
    except struct.error as e:
      raise genpy.DeserializationError(e) #most likely buffer underfill


  def serialize_numpy(self, buff, numpy):
    """
    serialize message with numpy array types into buffer
    :param buff: buffer, ``StringIO``
    :param numpy: numpy python module
    """
    try:
      _x = self
      buff.write(_struct_3I.pack(_x.header.seq, _x.header.stamp.secs, _x.header.stamp.nsecs))
      _x = self.header.frame_id
      length = len(_x)
      if python3 or type(_x) == unicode:
        _x = _x.encode('utf-8')
        length = len(_x)
      buff.write(struct.pack('<I%ss'%length, length, _x))
      _x = self
      buff.write(_struct_4f5if4i.pack(_x.longitude, _x.latitude, _x.altitude, _x.heading, _x.toleranceRadius, _x.holdTime, _x.index, _x.eventChannel, _x.altitudeRate, _x.speed, _x.status, _x.eventFlag, _x.type, _x.camAngle))
    except struct.error as se: self._check_types(se)
    except TypeError as te: self._check_types(te)

  def deserialize_numpy(self, str, numpy):
    """
    unpack serialized message in str into this message instance using numpy for array types
    :param str: byte array of serialized message, ``str``
    :param numpy: numpy python module
    """
    try:
      if self.header is None:
        self.header = std_msgs.msg.Header()
      end = 0
      _x = self
      start = end
      end += 12
      (_x.header.seq, _x.header.stamp.secs, _x.header.stamp.nsecs,) = _struct_3I.unpack(str[start:end])
      start = end
      end += 4
      (length,) = _struct_I.unpack(str[start:end])
      start = end
      end += length
      if python3:
        self.header.frame_id = str[start:end].decode('utf-8')
      else:
        self.header.frame_id = str[start:end]
      _x = self
      start = end
      end += 56
      (_x.longitude, _x.latitude, _x.altitude, _x.heading, _x.toleranceRadius, _x.holdTime, _x.index, _x.eventChannel, _x.altitudeRate, _x.speed, _x.status, _x.eventFlag, _x.type, _x.camAngle,) = _struct_4f5if4i.unpack(str[start:end])
      return self
    except struct.error as e:
      raise genpy.DeserializationError(e) #most likely buffer underfill

_struct_I = genpy.struct_I
_struct_3I = struct.Struct("<3I")
_struct_4f5if4i = struct.Struct("<4f5if4i")
