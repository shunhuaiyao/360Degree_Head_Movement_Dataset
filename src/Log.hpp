/**
 Author: Xavier Corbillon
 IMT Atlantique

 The class represent one Log line (Head position at a specific time for a specific video frame)
*/
#pragma once

//internal includes
#include "Quaternion.hpp"
#include "Timestamp.hpp"

//standard includes
#include <iostream>
#include <vector>

namespace IMT {

class Log
{
public:
  Log(Timestamp t, Timestamp pts, Quaternion q, size_t frameId, double yaw, double pitch): m_t(t), m_pts(pts), m_q(q), m_frameId(frameId), m_yaw(yaw), m_pitch(pitch) {};

  const Timestamp& GetTimestamp(void) const {return m_t;};
  const Quaternion& GetQuat(void) const {return m_q; };
  const size_t& GetFrameId(void) const {return m_frameId;};
  const double& GetYaw(void) const { return m_yaw; };
  const double& GetPitch(void) const { return m_pitch; };
  friend std::ostream& operator<< (std::ostream& stream, const Log& log)
  {
    stream << log.m_t << " " << log.m_frameId << " " /*<< log.m_pts << " "*/ << log.m_q << " " << log.m_yaw << " " << log.m_pitch;
	return stream;
  }
  Log operator-(const Timestamp& t) const
  {
    Log out(*this);
    out.m_t -= t;
    return out;
  }
  const auto& GetQuaternion(void) const {return m_q;}
private:
  Timestamp m_t;
  Timestamp m_pts;
  Quaternion m_q;
  size_t m_frameId;
  double m_yaw;
  double m_pitch;
};

}
