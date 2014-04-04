/* Auto-generated by genmsg_cpp for file /home/loki/Documents/ros_workspace/mikrokopter/mk_msgs/srv/SetGPSWaypoint.srv */
#ifndef MK_MSGS_SERVICE_SETGPSWAYPOINT_H
#define MK_MSGS_SERVICE_SETGPSWAYPOINT_H
#include <string>
#include <vector>
#include <map>
#include <ostream>
#include "ros/serialization.h"
#include "ros/builtin_message_traits.h"
#include "ros/message_operations.h"
#include "ros/time.h"

#include "ros/macros.h"

#include "ros/assert.h"

#include "ros/service_traits.h"




namespace mk_msgs
{
template <class ContainerAllocator>
struct SetGPSWaypointRequest_ {
  typedef SetGPSWaypointRequest_<ContainerAllocator> Type;

  SetGPSWaypointRequest_()
  : index(0)
  , latitude(0.0)
  , longitude(0.0)
  , altitude(0.0)
  , status(0)
  , holdTime(0)
  , altitudeRate(0.0)
  , speed(0.0)
  , heading(0.0)
  , toleranceRadius(0)
  , eventChannel(0)
  , eventFlag(0)
  , type(0)
  , camAngle(0)
  {
  }

  SetGPSWaypointRequest_(const ContainerAllocator& _alloc)
  : index(0)
  , latitude(0.0)
  , longitude(0.0)
  , altitude(0.0)
  , status(0)
  , holdTime(0)
  , altitudeRate(0.0)
  , speed(0.0)
  , heading(0.0)
  , toleranceRadius(0)
  , eventChannel(0)
  , eventFlag(0)
  , type(0)
  , camAngle(0)
  {
  }

  typedef int32_t _index_type;
  int32_t index;

  typedef float _latitude_type;
  float latitude;

  typedef float _longitude_type;
  float longitude;

  typedef float _altitude_type;
  float altitude;

  typedef int32_t _status_type;
  int32_t status;

  typedef int32_t _holdTime_type;
  int32_t holdTime;

  typedef float _altitudeRate_type;
  float altitudeRate;

  typedef float _speed_type;
  float speed;

  typedef float _heading_type;
  float heading;

  typedef int32_t _toleranceRadius_type;
  int32_t toleranceRadius;

  typedef int32_t _eventChannel_type;
  int32_t eventChannel;

  typedef int32_t _eventFlag_type;
  int32_t eventFlag;

  typedef int32_t _type_type;
  int32_t type;

  typedef int32_t _camAngle_type;
  int32_t camAngle;


  typedef boost::shared_ptr< ::mk_msgs::SetGPSWaypointRequest_<ContainerAllocator> > Ptr;
  typedef boost::shared_ptr< ::mk_msgs::SetGPSWaypointRequest_<ContainerAllocator>  const> ConstPtr;
  boost::shared_ptr<std::map<std::string, std::string> > __connection_header;
}; // struct SetGPSWaypointRequest
typedef  ::mk_msgs::SetGPSWaypointRequest_<std::allocator<void> > SetGPSWaypointRequest;

typedef boost::shared_ptr< ::mk_msgs::SetGPSWaypointRequest> SetGPSWaypointRequestPtr;
typedef boost::shared_ptr< ::mk_msgs::SetGPSWaypointRequest const> SetGPSWaypointRequestConstPtr;


template <class ContainerAllocator>
struct SetGPSWaypointResponse_ {
  typedef SetGPSWaypointResponse_<ContainerAllocator> Type;

  SetGPSWaypointResponse_()
  : ack(0)
  {
  }

  SetGPSWaypointResponse_(const ContainerAllocator& _alloc)
  : ack(0)
  {
  }

  typedef int16_t _ack_type;
  int16_t ack;


  typedef boost::shared_ptr< ::mk_msgs::SetGPSWaypointResponse_<ContainerAllocator> > Ptr;
  typedef boost::shared_ptr< ::mk_msgs::SetGPSWaypointResponse_<ContainerAllocator>  const> ConstPtr;
  boost::shared_ptr<std::map<std::string, std::string> > __connection_header;
}; // struct SetGPSWaypointResponse
typedef  ::mk_msgs::SetGPSWaypointResponse_<std::allocator<void> > SetGPSWaypointResponse;

typedef boost::shared_ptr< ::mk_msgs::SetGPSWaypointResponse> SetGPSWaypointResponsePtr;
typedef boost::shared_ptr< ::mk_msgs::SetGPSWaypointResponse const> SetGPSWaypointResponseConstPtr;

struct SetGPSWaypoint
{

typedef SetGPSWaypointRequest Request;
typedef SetGPSWaypointResponse Response;
Request request;
Response response;

typedef Request RequestType;
typedef Response ResponseType;
}; // struct SetGPSWaypoint
} // namespace mk_msgs

namespace ros
{
namespace message_traits
{
template<class ContainerAllocator> struct IsMessage< ::mk_msgs::SetGPSWaypointRequest_<ContainerAllocator> > : public TrueType {};
template<class ContainerAllocator> struct IsMessage< ::mk_msgs::SetGPSWaypointRequest_<ContainerAllocator>  const> : public TrueType {};
template<class ContainerAllocator>
struct MD5Sum< ::mk_msgs::SetGPSWaypointRequest_<ContainerAllocator> > {
  static const char* value() 
  {
    return "c1ce5d208dc30cef244a48ead733472b";
  }

  static const char* value(const  ::mk_msgs::SetGPSWaypointRequest_<ContainerAllocator> &) { return value(); } 
  static const uint64_t static_value1 = 0xc1ce5d208dc30cefULL;
  static const uint64_t static_value2 = 0x244a48ead733472bULL;
};

template<class ContainerAllocator>
struct DataType< ::mk_msgs::SetGPSWaypointRequest_<ContainerAllocator> > {
  static const char* value() 
  {
    return "mk_msgs/SetGPSWaypointRequest";
  }

  static const char* value(const  ::mk_msgs::SetGPSWaypointRequest_<ContainerAllocator> &) { return value(); } 
};

template<class ContainerAllocator>
struct Definition< ::mk_msgs::SetGPSWaypointRequest_<ContainerAllocator> > {
  static const char* value() 
  {
    return "\n\
\n\
\n\
int32 index\n\
float32 latitude\n\
float32 longitude\n\
float32 altitude\n\
int32 status\n\
int32 holdTime\n\
float32 altitudeRate\n\
float32 speed\n\
\n\
float32 heading\n\
int32 toleranceRadius\n\
int32 eventChannel\n\
int32 eventFlag\n\
int32 type\n\
int32 camAngle\n\
\n\
";
  }

  static const char* value(const  ::mk_msgs::SetGPSWaypointRequest_<ContainerAllocator> &) { return value(); } 
};

template<class ContainerAllocator> struct IsFixedSize< ::mk_msgs::SetGPSWaypointRequest_<ContainerAllocator> > : public TrueType {};
} // namespace message_traits
} // namespace ros


namespace ros
{
namespace message_traits
{
template<class ContainerAllocator> struct IsMessage< ::mk_msgs::SetGPSWaypointResponse_<ContainerAllocator> > : public TrueType {};
template<class ContainerAllocator> struct IsMessage< ::mk_msgs::SetGPSWaypointResponse_<ContainerAllocator>  const> : public TrueType {};
template<class ContainerAllocator>
struct MD5Sum< ::mk_msgs::SetGPSWaypointResponse_<ContainerAllocator> > {
  static const char* value() 
  {
    return "efd27fb04e2a6ce3c9ff1f47eb32e7bb";
  }

  static const char* value(const  ::mk_msgs::SetGPSWaypointResponse_<ContainerAllocator> &) { return value(); } 
  static const uint64_t static_value1 = 0xefd27fb04e2a6ce3ULL;
  static const uint64_t static_value2 = 0xc9ff1f47eb32e7bbULL;
};

template<class ContainerAllocator>
struct DataType< ::mk_msgs::SetGPSWaypointResponse_<ContainerAllocator> > {
  static const char* value() 
  {
    return "mk_msgs/SetGPSWaypointResponse";
  }

  static const char* value(const  ::mk_msgs::SetGPSWaypointResponse_<ContainerAllocator> &) { return value(); } 
};

template<class ContainerAllocator>
struct Definition< ::mk_msgs::SetGPSWaypointResponse_<ContainerAllocator> > {
  static const char* value() 
  {
    return "\n\
\n\
\n\
int16 ack\n\
\n\
\n\
";
  }

  static const char* value(const  ::mk_msgs::SetGPSWaypointResponse_<ContainerAllocator> &) { return value(); } 
};

template<class ContainerAllocator> struct IsFixedSize< ::mk_msgs::SetGPSWaypointResponse_<ContainerAllocator> > : public TrueType {};
} // namespace message_traits
} // namespace ros

namespace ros
{
namespace serialization
{

template<class ContainerAllocator> struct Serializer< ::mk_msgs::SetGPSWaypointRequest_<ContainerAllocator> >
{
  template<typename Stream, typename T> inline static void allInOne(Stream& stream, T m)
  {
    stream.next(m.index);
    stream.next(m.latitude);
    stream.next(m.longitude);
    stream.next(m.altitude);
    stream.next(m.status);
    stream.next(m.holdTime);
    stream.next(m.altitudeRate);
    stream.next(m.speed);
    stream.next(m.heading);
    stream.next(m.toleranceRadius);
    stream.next(m.eventChannel);
    stream.next(m.eventFlag);
    stream.next(m.type);
    stream.next(m.camAngle);
  }

  ROS_DECLARE_ALLINONE_SERIALIZER;
}; // struct SetGPSWaypointRequest_
} // namespace serialization
} // namespace ros


namespace ros
{
namespace serialization
{

template<class ContainerAllocator> struct Serializer< ::mk_msgs::SetGPSWaypointResponse_<ContainerAllocator> >
{
  template<typename Stream, typename T> inline static void allInOne(Stream& stream, T m)
  {
    stream.next(m.ack);
  }

  ROS_DECLARE_ALLINONE_SERIALIZER;
}; // struct SetGPSWaypointResponse_
} // namespace serialization
} // namespace ros

namespace ros
{
namespace service_traits
{
template<>
struct MD5Sum<mk_msgs::SetGPSWaypoint> {
  static const char* value() 
  {
    return "167ab490b44feaf311e87f25a5b4ea94";
  }

  static const char* value(const mk_msgs::SetGPSWaypoint&) { return value(); } 
};

template<>
struct DataType<mk_msgs::SetGPSWaypoint> {
  static const char* value() 
  {
    return "mk_msgs/SetGPSWaypoint";
  }

  static const char* value(const mk_msgs::SetGPSWaypoint&) { return value(); } 
};

template<class ContainerAllocator>
struct MD5Sum<mk_msgs::SetGPSWaypointRequest_<ContainerAllocator> > {
  static const char* value() 
  {
    return "167ab490b44feaf311e87f25a5b4ea94";
  }

  static const char* value(const mk_msgs::SetGPSWaypointRequest_<ContainerAllocator> &) { return value(); } 
};

template<class ContainerAllocator>
struct DataType<mk_msgs::SetGPSWaypointRequest_<ContainerAllocator> > {
  static const char* value() 
  {
    return "mk_msgs/SetGPSWaypoint";
  }

  static const char* value(const mk_msgs::SetGPSWaypointRequest_<ContainerAllocator> &) { return value(); } 
};

template<class ContainerAllocator>
struct MD5Sum<mk_msgs::SetGPSWaypointResponse_<ContainerAllocator> > {
  static const char* value() 
  {
    return "167ab490b44feaf311e87f25a5b4ea94";
  }

  static const char* value(const mk_msgs::SetGPSWaypointResponse_<ContainerAllocator> &) { return value(); } 
};

template<class ContainerAllocator>
struct DataType<mk_msgs::SetGPSWaypointResponse_<ContainerAllocator> > {
  static const char* value() 
  {
    return "mk_msgs/SetGPSWaypoint";
  }

  static const char* value(const mk_msgs::SetGPSWaypointResponse_<ContainerAllocator> &) { return value(); } 
};

} // namespace service_traits
} // namespace ros

#endif // MK_MSGS_SERVICE_SETGPSWAYPOINT_H
