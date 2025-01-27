/* 
This file is part of FAST-LIVO2: Fast, Direct LiDAR-Inertial-Visual Odometry.

Developer: Chunran Zheng <zhengcr@connect.hku.hk>

For commercial use, please contact me at <zhengcr@connect.hku.hk> or
Prof. Fu Zhang at <fuzhang@hku.hk>.

This file is subject to the terms and conditions outlined in the 'LICENSE' file,
which is included as part of this source code package.
*/

#ifndef LIVO_POINT_H_
#define LIVO_POINT_H_

#include <boost/noncopyable.hpp>
#include "common_lib.h"
#include "frame.h"

class Feature;

/// A visual map point on the surface of the scene.
class VisualPoint : boost::noncopyable
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /// \brief 3d pos of the point in the world coordinate frame.
  Vector3d pos_;
  /// \brief Surface normal at point.
  Vector3d normal_;             
  /// \brief Inverse covariance matrix of normal estimation.
  Matrix3d normal_information_; 
  /// \brief Last updated normal vector.
  Vector3d previous_normal_;    
  /// \brief Reference patches which observe the point.
  list<Feature *> obs_;         
  /// \brief Covariance of the point.
  Eigen::Matrix3d covariance_; 
  /// \brief True if the point is converged.
  bool is_converged_;           
  /// \brief True if the normal is initialized.
  bool is_normal_initialized_;  
  /// \brief True if the point has a reference patch.
  bool has_ref_patch_;          
  /// \brief Reference patch of the point.
  Feature *ref_patch;          

  VisualPoint(const Vector3d &pos);
  ~VisualPoint();
  void findMinScoreFeature(const Vector3d &framepos, Feature *&ftr) const;
  void deleteNonRefPatchFeatures();
  void deleteFeatureRef(Feature *ftr);
  void addFrameRef(Feature *ftr);
  
  /// \brief Check if pos is closed-view to this::pos_
  bool getCloseViewObs(const Vector3d &pos, Feature *&obs) const;
};

#endif // LIVO_POINT_H_
