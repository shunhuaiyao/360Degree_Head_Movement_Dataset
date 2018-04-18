/**
 Author: Xavier Corbillon
 IMT Atlantique

 This class manage a head position Log
*/
#pragma once

//library includes
#include <Eigen/Geometry>
#include <osvr/Util/EigenInterop.h>

namespace IMT {
  typedef Eigen::Quaterniond Quaternion;

  inline Quaternion ToQuaternion(const OSVR_Quaternion& quat) {return osvr::util::fromQuat(quat);}
  //return Quaternion(osvrQuatGetW(&quat), -osvrQuatGetZ(&quat), -osvrQuatGetX(&quat), osvrQuatGetY(&quat));}

  inline bool operator==(const Quaternion& quat1, const Quaternion& quat2)
  {
    return quat1.w() == quat2.w() && quat1.x() == quat2.x() &&
           quat1.y() == quat2.y() && quat1.z() == quat2.z();
  }
  inline bool operator!=(const Quaternion& quat1, const Quaternion& quat2)
  {
    return !(quat1 == quat2);
  }

  inline std::ostream& operator<< (std::ostream& stream, const Quaternion& quat)
  {
    stream << quat.w() << " " << quat.x() << " " << quat.y() << " " << quat.z();
	return stream;
  }
}
