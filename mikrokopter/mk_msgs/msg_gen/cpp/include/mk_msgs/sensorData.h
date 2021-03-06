/* Auto-generated by genmsg_cpp for file /home/loki/Documents/ros_workspace/mikrokopter/mk_msgs/msg/sensorData.msg */
#ifndef MK_MSGS_MESSAGE_SENSORDATA_H
#define MK_MSGS_MESSAGE_SENSORDATA_H
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

#include "std_msgs/Header.h"

namespace mk_msgs
{
template <class ContainerAllocator>
struct sensorData_ {
  typedef sensorData_<ContainerAllocator> Type;

  sensorData_()
  : header()
  , Longitude(0.0)
  , Latitude(0.0)
  , gpsAltitude(0.0)
  , SatsInUse(0)
  , NumberOfWPs(0)
  , CurrentWP(0)
  , ErrorCode(0)
  , Free(0)
  , PH(0)
  , CH(0)
  , RangeLimit(0)
  , NoSerialLink(0)
  , TargetReached(0)
  , Manual(0)
  , GPSOK(0)
  , MotorsOn(0)
  , Flying(0)
  , LowBat(0)
  , CareFree(0)
  , AltHld(0)
  , Failsafe(0)
  , Altitude(0.0)
  , FlyingTime(0.0)
  , Battery(0.0)
  , GroundSpeed(0.0)
  , Heading(0.0)
  , CompassHeading(0.0)
  , Nick(0.0)
  , Roll(0.0)
  , RCQuality(0.0)
  , zSpeed(0.0)
  , TargetHoldTime(0.0)
  , Gas(0.0)
  , Current(0.0)
  , UsedCapacity(0.0)
  {
  }

  sensorData_(const ContainerAllocator& _alloc)
  : header(_alloc)
  , Longitude(0.0)
  , Latitude(0.0)
  , gpsAltitude(0.0)
  , SatsInUse(0)
  , NumberOfWPs(0)
  , CurrentWP(0)
  , ErrorCode(0)
  , Free(0)
  , PH(0)
  , CH(0)
  , RangeLimit(0)
  , NoSerialLink(0)
  , TargetReached(0)
  , Manual(0)
  , GPSOK(0)
  , MotorsOn(0)
  , Flying(0)
  , LowBat(0)
  , CareFree(0)
  , AltHld(0)
  , Failsafe(0)
  , Altitude(0.0)
  , FlyingTime(0.0)
  , Battery(0.0)
  , GroundSpeed(0.0)
  , Heading(0.0)
  , CompassHeading(0.0)
  , Nick(0.0)
  , Roll(0.0)
  , RCQuality(0.0)
  , zSpeed(0.0)
  , TargetHoldTime(0.0)
  , Gas(0.0)
  , Current(0.0)
  , UsedCapacity(0.0)
  {
  }

  typedef  ::std_msgs::Header_<ContainerAllocator>  _header_type;
   ::std_msgs::Header_<ContainerAllocator>  header;

  typedef float _Longitude_type;
  float Longitude;

  typedef float _Latitude_type;
  float Latitude;

  typedef float _gpsAltitude_type;
  float gpsAltitude;

  typedef int32_t _SatsInUse_type;
  int32_t SatsInUse;

  typedef int32_t _NumberOfWPs_type;
  int32_t NumberOfWPs;

  typedef int32_t _CurrentWP_type;
  int32_t CurrentWP;

  typedef int32_t _ErrorCode_type;
  int32_t ErrorCode;

  typedef int32_t _Free_type;
  int32_t Free;

  typedef int32_t _PH_type;
  int32_t PH;

  typedef int32_t _CH_type;
  int32_t CH;

  typedef int32_t _RangeLimit_type;
  int32_t RangeLimit;

  typedef int32_t _NoSerialLink_type;
  int32_t NoSerialLink;

  typedef int32_t _TargetReached_type;
  int32_t TargetReached;

  typedef int32_t _Manual_type;
  int32_t Manual;

  typedef int32_t _GPSOK_type;
  int32_t GPSOK;

  typedef int32_t _MotorsOn_type;
  int32_t MotorsOn;

  typedef int32_t _Flying_type;
  int32_t Flying;

  typedef int32_t _LowBat_type;
  int32_t LowBat;

  typedef int32_t _CareFree_type;
  int32_t CareFree;

  typedef int32_t _AltHld_type;
  int32_t AltHld;

  typedef int32_t _Failsafe_type;
  int32_t Failsafe;

  typedef float _Altitude_type;
  float Altitude;

  typedef float _FlyingTime_type;
  float FlyingTime;

  typedef float _Battery_type;
  float Battery;

  typedef float _GroundSpeed_type;
  float GroundSpeed;

  typedef float _Heading_type;
  float Heading;

  typedef float _CompassHeading_type;
  float CompassHeading;

  typedef float _Nick_type;
  float Nick;

  typedef float _Roll_type;
  float Roll;

  typedef float _RCQuality_type;
  float RCQuality;

  typedef float _zSpeed_type;
  float zSpeed;

  typedef float _TargetHoldTime_type;
  float TargetHoldTime;

  typedef float _Gas_type;
  float Gas;

  typedef float _Current_type;
  float Current;

  typedef float _UsedCapacity_type;
  float UsedCapacity;


  typedef boost::shared_ptr< ::mk_msgs::sensorData_<ContainerAllocator> > Ptr;
  typedef boost::shared_ptr< ::mk_msgs::sensorData_<ContainerAllocator>  const> ConstPtr;
  boost::shared_ptr<std::map<std::string, std::string> > __connection_header;
}; // struct sensorData
typedef  ::mk_msgs::sensorData_<std::allocator<void> > sensorData;

typedef boost::shared_ptr< ::mk_msgs::sensorData> sensorDataPtr;
typedef boost::shared_ptr< ::mk_msgs::sensorData const> sensorDataConstPtr;


template<typename ContainerAllocator>
std::ostream& operator<<(std::ostream& s, const  ::mk_msgs::sensorData_<ContainerAllocator> & v)
{
  ros::message_operations::Printer< ::mk_msgs::sensorData_<ContainerAllocator> >::stream(s, "", v);
  return s;}

} // namespace mk_msgs

namespace ros
{
namespace message_traits
{
template<class ContainerAllocator> struct IsMessage< ::mk_msgs::sensorData_<ContainerAllocator> > : public TrueType {};
template<class ContainerAllocator> struct IsMessage< ::mk_msgs::sensorData_<ContainerAllocator>  const> : public TrueType {};
template<class ContainerAllocator>
struct MD5Sum< ::mk_msgs::sensorData_<ContainerAllocator> > {
  static const char* value() 
  {
    return "7e8cef3e80536e777b04973b40e71d35";
  }

  static const char* value(const  ::mk_msgs::sensorData_<ContainerAllocator> &) { return value(); } 
  static const uint64_t static_value1 = 0x7e8cef3e80536e77ULL;
  static const uint64_t static_value2 = 0x7b04973b40e71d35ULL;
};

template<class ContainerAllocator>
struct DataType< ::mk_msgs::sensorData_<ContainerAllocator> > {
  static const char* value() 
  {
    return "mk_msgs/sensorData";
  }

  static const char* value(const  ::mk_msgs::sensorData_<ContainerAllocator> &) { return value(); } 
};

template<class ContainerAllocator>
struct Definition< ::mk_msgs::sensorData_<ContainerAllocator> > {
  static const char* value() 
  {
    return "Header header\n\
\n\
float32 Longitude\n\
float32 Latitude\n\
float32 gpsAltitude\n\
\n\
int32 SatsInUse\n\
int32 NumberOfWPs\n\
int32 CurrentWP\n\
int32 ErrorCode\n\
#int32 NCMode\n\
#int32 FCMode1\n\
#int32 FCMode2\n\
\n\
int32 Free\n\
int32 PH\n\
int32 CH\n\
int32 RangeLimit\n\
int32 NoSerialLink\n\
int32 TargetReached\n\
int32 Manual\n\
int32 GPSOK\n\
\n\
int32 MotorsOn\n\
int32 Flying\n\
#int32 Calibrate\n\
#int32 Start\n\
#int32 EmergencyLanding\n\
int32 LowBat\n\
#int32 VarioTimUp\n\
#int32 VarioTrimDown\n\
\n\
int32 CareFree\n\
int32 AltHld\n\
int32 Failsafe \n\
\n\
float32 Altitude\n\
float32 FlyingTime\n\
float32 Battery\n\
float32 GroundSpeed\n\
float32 Heading\n\
float32 CompassHeading\n\
float32 Nick\n\
float32 Roll\n\
float32 RCQuality\n\
float32 zSpeed\n\
float32 TargetHoldTime\n\
float32 Gas\n\
float32 Current\n\
float32 UsedCapacity\n\
\n\
\n\
================================================================================\n\
MSG: std_msgs/Header\n\
# Standard metadata for higher-level stamped data types.\n\
# This is generally used to communicate timestamped data \n\
# in a particular coordinate frame.\n\
# \n\
# sequence ID: consecutively increasing ID \n\
uint32 seq\n\
#Two-integer timestamp that is expressed as:\n\
# * stamp.secs: seconds (stamp_secs) since epoch\n\
# * stamp.nsecs: nanoseconds since stamp_secs\n\
# time-handling sugar is provided by the client library\n\
time stamp\n\
#Frame this data is associated with\n\
# 0: no frame\n\
# 1: global frame\n\
string frame_id\n\
\n\
";
  }

  static const char* value(const  ::mk_msgs::sensorData_<ContainerAllocator> &) { return value(); } 
};

template<class ContainerAllocator> struct HasHeader< ::mk_msgs::sensorData_<ContainerAllocator> > : public TrueType {};
template<class ContainerAllocator> struct HasHeader< const ::mk_msgs::sensorData_<ContainerAllocator> > : public TrueType {};
} // namespace message_traits
} // namespace ros

namespace ros
{
namespace serialization
{

template<class ContainerAllocator> struct Serializer< ::mk_msgs::sensorData_<ContainerAllocator> >
{
  template<typename Stream, typename T> inline static void allInOne(Stream& stream, T m)
  {
    stream.next(m.header);
    stream.next(m.Longitude);
    stream.next(m.Latitude);
    stream.next(m.gpsAltitude);
    stream.next(m.SatsInUse);
    stream.next(m.NumberOfWPs);
    stream.next(m.CurrentWP);
    stream.next(m.ErrorCode);
    stream.next(m.Free);
    stream.next(m.PH);
    stream.next(m.CH);
    stream.next(m.RangeLimit);
    stream.next(m.NoSerialLink);
    stream.next(m.TargetReached);
    stream.next(m.Manual);
    stream.next(m.GPSOK);
    stream.next(m.MotorsOn);
    stream.next(m.Flying);
    stream.next(m.LowBat);
    stream.next(m.CareFree);
    stream.next(m.AltHld);
    stream.next(m.Failsafe);
    stream.next(m.Altitude);
    stream.next(m.FlyingTime);
    stream.next(m.Battery);
    stream.next(m.GroundSpeed);
    stream.next(m.Heading);
    stream.next(m.CompassHeading);
    stream.next(m.Nick);
    stream.next(m.Roll);
    stream.next(m.RCQuality);
    stream.next(m.zSpeed);
    stream.next(m.TargetHoldTime);
    stream.next(m.Gas);
    stream.next(m.Current);
    stream.next(m.UsedCapacity);
  }

  ROS_DECLARE_ALLINONE_SERIALIZER;
}; // struct sensorData_
} // namespace serialization
} // namespace ros

namespace ros
{
namespace message_operations
{

template<class ContainerAllocator>
struct Printer< ::mk_msgs::sensorData_<ContainerAllocator> >
{
  template<typename Stream> static void stream(Stream& s, const std::string& indent, const  ::mk_msgs::sensorData_<ContainerAllocator> & v) 
  {
    s << indent << "header: ";
s << std::endl;
    Printer< ::std_msgs::Header_<ContainerAllocator> >::stream(s, indent + "  ", v.header);
    s << indent << "Longitude: ";
    Printer<float>::stream(s, indent + "  ", v.Longitude);
    s << indent << "Latitude: ";
    Printer<float>::stream(s, indent + "  ", v.Latitude);
    s << indent << "gpsAltitude: ";
    Printer<float>::stream(s, indent + "  ", v.gpsAltitude);
    s << indent << "SatsInUse: ";
    Printer<int32_t>::stream(s, indent + "  ", v.SatsInUse);
    s << indent << "NumberOfWPs: ";
    Printer<int32_t>::stream(s, indent + "  ", v.NumberOfWPs);
    s << indent << "CurrentWP: ";
    Printer<int32_t>::stream(s, indent + "  ", v.CurrentWP);
    s << indent << "ErrorCode: ";
    Printer<int32_t>::stream(s, indent + "  ", v.ErrorCode);
    s << indent << "Free: ";
    Printer<int32_t>::stream(s, indent + "  ", v.Free);
    s << indent << "PH: ";
    Printer<int32_t>::stream(s, indent + "  ", v.PH);
    s << indent << "CH: ";
    Printer<int32_t>::stream(s, indent + "  ", v.CH);
    s << indent << "RangeLimit: ";
    Printer<int32_t>::stream(s, indent + "  ", v.RangeLimit);
    s << indent << "NoSerialLink: ";
    Printer<int32_t>::stream(s, indent + "  ", v.NoSerialLink);
    s << indent << "TargetReached: ";
    Printer<int32_t>::stream(s, indent + "  ", v.TargetReached);
    s << indent << "Manual: ";
    Printer<int32_t>::stream(s, indent + "  ", v.Manual);
    s << indent << "GPSOK: ";
    Printer<int32_t>::stream(s, indent + "  ", v.GPSOK);
    s << indent << "MotorsOn: ";
    Printer<int32_t>::stream(s, indent + "  ", v.MotorsOn);
    s << indent << "Flying: ";
    Printer<int32_t>::stream(s, indent + "  ", v.Flying);
    s << indent << "LowBat: ";
    Printer<int32_t>::stream(s, indent + "  ", v.LowBat);
    s << indent << "CareFree: ";
    Printer<int32_t>::stream(s, indent + "  ", v.CareFree);
    s << indent << "AltHld: ";
    Printer<int32_t>::stream(s, indent + "  ", v.AltHld);
    s << indent << "Failsafe: ";
    Printer<int32_t>::stream(s, indent + "  ", v.Failsafe);
    s << indent << "Altitude: ";
    Printer<float>::stream(s, indent + "  ", v.Altitude);
    s << indent << "FlyingTime: ";
    Printer<float>::stream(s, indent + "  ", v.FlyingTime);
    s << indent << "Battery: ";
    Printer<float>::stream(s, indent + "  ", v.Battery);
    s << indent << "GroundSpeed: ";
    Printer<float>::stream(s, indent + "  ", v.GroundSpeed);
    s << indent << "Heading: ";
    Printer<float>::stream(s, indent + "  ", v.Heading);
    s << indent << "CompassHeading: ";
    Printer<float>::stream(s, indent + "  ", v.CompassHeading);
    s << indent << "Nick: ";
    Printer<float>::stream(s, indent + "  ", v.Nick);
    s << indent << "Roll: ";
    Printer<float>::stream(s, indent + "  ", v.Roll);
    s << indent << "RCQuality: ";
    Printer<float>::stream(s, indent + "  ", v.RCQuality);
    s << indent << "zSpeed: ";
    Printer<float>::stream(s, indent + "  ", v.zSpeed);
    s << indent << "TargetHoldTime: ";
    Printer<float>::stream(s, indent + "  ", v.TargetHoldTime);
    s << indent << "Gas: ";
    Printer<float>::stream(s, indent + "  ", v.Gas);
    s << indent << "Current: ";
    Printer<float>::stream(s, indent + "  ", v.Current);
    s << indent << "UsedCapacity: ";
    Printer<float>::stream(s, indent + "  ", v.UsedCapacity);
  }
};


} // namespace message_operations
} // namespace ros

#endif // MK_MSGS_MESSAGE_SENSORDATA_H

