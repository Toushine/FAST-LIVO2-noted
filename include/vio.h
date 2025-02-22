/* 
This file is part of FAST-LIVO2: Fast, Direct LiDAR-Inertial-Visual Odometry.

Developer: Chunran Zheng <zhengcr@connect.hku.hk>

For commercial use, please contact me at <zhengcr@connect.hku.hk> or
Prof. Fu Zhang at <fuzhang@hku.hk>.

This file is subject to the terms and conditions outlined in the 'LICENSE' file,
which is included as part of this source code package.
*/

#ifndef VIO_H_
#define VIO_H_

#include "voxel_map.h"
#include "feature.h"
#include <opencv2/imgproc/imgproc_c.h>
#include <pcl/filters/voxel_grid.h>
#include <set>
#include <vikit/math_utils.h>
#include <vikit/robust_cost.h>
#include <vikit/vision.h>
#include <vikit/pinhole_camera.h>

struct ScaledPixel {
  int scale;

  float u_ref;
  float v_ref;

  int u_ref_i;
  int v_ref_i;

  float subpix_u_ref;
  float subpix_v_ref;

  float w_ref_tl;
  float w_ref_tr;
  float w_ref_bl;
  float w_ref_br;

  ScaledPixel(const V2D pc, const int scale) {
    this->scale = scale;
    u_ref = pc[0];
    v_ref = pc[1];

    u_ref_i = floorf(u_ref / scale) * scale;
    v_ref_i = floorf(v_ref / scale) * scale;

    subpix_u_ref = (u_ref - u_ref_i) / scale;
    subpix_v_ref = (v_ref - v_ref_i) / scale;

    w_ref_tl = (1.0 - subpix_u_ref) * (1.0 - subpix_v_ref);
    w_ref_tr = subpix_u_ref * (1.0 - subpix_v_ref);
    w_ref_bl = (1.0 - subpix_u_ref) * subpix_v_ref;
    w_ref_br = subpix_u_ref * subpix_v_ref;
  }

  std::array<float, 2> GetDuDv(const uint8_t *img_ptr, const int width) {
    std::array<float, 2> du_dv = {0,0};
    du_dv[0] =   0.5f *
              ((w_ref_tl * img_ptr[scale] + w_ref_tr * img_ptr[scale * 2] + w_ref_bl * img_ptr[scale * width + scale] +
                w_ref_br * img_ptr[scale * width + scale * 2]) -
               (w_ref_tl * img_ptr[-scale] + w_ref_tr * img_ptr[0] + w_ref_bl * img_ptr[scale * width - scale] + w_ref_br * img_ptr[scale * width]));

    du_dv[1] =  0.5f *
              ((w_ref_tl * img_ptr[scale * width] + w_ref_tr * img_ptr[scale + scale * width] + w_ref_bl * img_ptr[width * scale * 2] +
                w_ref_br * img_ptr[width * scale * 2 + scale]) -
               (w_ref_tl * img_ptr[-scale * width] + w_ref_tr * img_ptr[-scale * width + scale] + w_ref_bl * img_ptr[0] + w_ref_br * img_ptr[scale]));

    return du_dv;
  }

  V3F GetInterpolatedPixel(const cv::Mat img, const int width) {
    V3F pixel;
    uint8_t *img_ptr = (uint8_t *)img.data + ((v_ref_i)*width + (u_ref_i)) * 3;
    for (int i = 0; i < 3; i++) {
      pixel[i] = w_ref_tl * img_ptr[i] + w_ref_tr * img_ptr[i + 3] +
                 w_ref_bl * img_ptr[i + width * 3] +
                 w_ref_br * img_ptr[width * 3 + i + 3];
    }

    return pixel;
  }
};

struct SubSparseMap
{
  vector<float> propa_errors;
  vector<float> errors;
  vector<vector<float>> warp_patch;
  vector<int> search_levels;
  vector<VisualPoint *> voxel_points;
  vector<double> inv_expo_list;
  vector<pointWithVar> add_from_voxel_map;

  SubSparseMap()
  {
    propa_errors.reserve(SIZE_LARGE);
    errors.reserve(SIZE_LARGE);
    warp_patch.reserve(SIZE_LARGE);
    search_levels.reserve(SIZE_LARGE);
    voxel_points.reserve(SIZE_LARGE);
    inv_expo_list.reserve(SIZE_LARGE);
    add_from_voxel_map.reserve(SIZE_SMALL);
  };

  void reset()
  {
    ROS_DEBUG("Reset SubSparseMap.");
    propa_errors.clear();
    errors.clear();
    warp_patch.clear();
    search_levels.clear();
    voxel_points.clear();
    inv_expo_list.clear();
    add_from_voxel_map.clear();
  }
};

class Warp
{
public:
  Matrix2d A_cur_ref;
  int search_level;
  Warp(int level, Matrix2d warp_matrix) : search_level(level), A_cur_ref(warp_matrix) {}
  ~Warp() {}
};

class VOXEL_POINTS
{
public:
  /// \brief The visual-points in the voxel
  std::vector<VisualPoint *> voxel_points;
  /// \brief The number of points in the voxel
  int count;

public:
  /// \brief Constructor with the number of points in the voxel
  VOXEL_POINTS(int num) : count(num) {}

  /// \brief Destructor to release the memory
  ~VOXEL_POINTS() 
  { 
    for (VisualPoint* vp : voxel_points) 
    {
      if (vp != nullptr) { delete vp; vp = nullptr; }
    }
  }
};

/// \brief VIO manager
class VIOManager
{
public:
  /// \brief Divide the image into grids
  int grid_size;
  /// \brief abstract camera model
  vk::AbstractCamera *cam;
  /// \brief pinhole camera model
  vk::PinholeCamera *pinhole_cam;

  /// \brief State of vio system
  StatesGroup *state;
  /// \brief State of vio system after propagation
  StatesGroup *state_propagat;

  /// \brief jacobian of the rotatin vector w.r.t rotation matrix
  /// \note rotation matrix itself
  M3D Jdphi_dR;
  /// \brief jacobian of the position w.r.t translation
  M3D Jdp_dt;
  /// \brief jacobian of the position w.r.t rotation
  M3D Jdp_dR;
  
  /// \brief rotation of imu frame in lidar frame
  M3D Rli;
  /// \brief translation of imu frame in lidar frame
  V3D Pli;
  
  /// \brief rotation of imu frame in camera frame
  M3D Rci;
  /// \brief translation of imu frmae in camera frame
  V3D Pci;

  /// \brief rotation of lidar frame in camera frame
  M3D Rcl;
  /// \brief translation of lidar frame in camera frame
  V3D Pcl;

  /// \brief rotation of camera frame in world frame
  M3D Rcw;
  /// \brief translation of world frame in camera frame
  V3D Pcw;

  /* all resize to grid_n_width*gtid_n_height */

  /// \brief save CellType
  vector<int> grid_num;
  // vector<int> map_index; // useless
  /// \brief whether the grid is in the border of the image
  vector<int> border_flag;
  /// \brief whether the grid is updated
  vector<int> update_flag;
  /// \brief save minimum distance
  vector<float> map_dist;
  vector<float> scan_value;
  vector<float> patch_buffer;

  /// \brief whether to enable Normal Refine V.E
  bool normal_en;
  /// \brief whether to enable Inverse Composional formulation in formula(23)
  bool inverse_composition_en;
  /// \brief whether to enable extimate the image exposure time
  bool exposure_estimate_en;
  /// \brief whether to enable on-demand voxel raycasting VII.A.2
  bool raycast_en;
  /// \brief whether has reference patch cache
  bool has_ref_patch_cache;
  /// \brief whether to use NCC(Normalized Cross Correlation) to compute the
  /// similarity between two patches
  bool ncc_en = false;
  /// \brief whether to output colmap
  bool colmap_output_en = false;

  /// \brief Image width(pixel)
  int width;
  /// \brief Image height(pixel)
  int height;
  /// \brief resize factor of Image
  double image_resize_factor;

  /// \brief pixels num in grid width
  int grid_n_width;
  /// \brief pixels num in grid height
  int grid_n_height;
  /// \brief pixels num in grid 
  int length;

  /// \brief focal length in x direction
  double fx;
  /// \brief focal length in y direction
  double fy;
  /// \brief Coordinates of the principal point in x direction
  double cx;
  /// \brief Coordinates of the principal point in y direction
  double cy;

  int patch_pyrimid_level, patch_size, patch_size_total, patch_size_half, border, warp_len;
  int max_iterations, total_points;

  double img_point_cov, outlier_threshold, ncc_thre;

  SubSparseMap *visual_submap;

  /// \brief sample points in image grids(length)
  std::vector<std::vector<V3D>> rays_with_sample_points;

  double compute_jacobian_time, update_ekf_time;
  double ave_total = 0;
  // double ave_build_residual_time = 0;
  // double ave_ekf_time = 0;

  int frame_count = 0;
  bool plot_flag;

  /* Related to ESIKF */

  /// \brief K*H
  Matrix<double, DIM_STATE, DIM_STATE> G;
  /// \brief H^T*H
  Matrix<double, DIM_STATE, DIM_STATE> H_T_H;
  
  MatrixXd K, H_sub_inv;

  ofstream fout_camera, fout_colmap;
  /// \brief Voxel Map(Local Mapping)
  unordered_map<VOXEL_LOCATION, VOXEL_POINTS *> feat_map;
  unordered_map<VOXEL_LOCATION, int> sub_feat_map; 
  unordered_map<int, Warp *> warp_map;
  vector<VisualPoint *> retrieve_voxel_points;
  vector<pointWithVar> append_voxel_points;
  FramePtr new_frame_;
  cv::Mat img_cp, img_rgb, img_test;

  enum CellType
  {
    TYPE_MAP = 1,
    TYPE_POINTCLOUD,
    TYPE_UNKNOWN
  };

  VIOManager();
  ~VIOManager();
  void updateStateInverse(cv::Mat img, int level);
  void updateState(cv::Mat img, int level);
  void processFrame(cv::Mat &img, vector<pointWithVar> &pg, const unordered_map<VOXEL_LOCATION, VoxelOctoTree *> &feat_map, double img_time);
  void retrieveFromVisualSparseMap(cv::Mat img, vector<pointWithVar> &pg, const unordered_map<VOXEL_LOCATION, VoxelOctoTree *> &plane_map);
  void generateVisualMapPoints(cv::Mat img, vector<pointWithVar> &pg);
  void setImuToLidarExtrinsic(const V3D &transl, const M3D &rot);
  void setLidarToCameraExtrinsic(vector<double> &R, vector<double> &P);
  void initializeVIO();
  void getImagePatch(cv::Mat img, V2D pc, float *patch_tmp, int level);
  void computeProjectionJacobian(V3D p, MD(2, 3) & J);
  void computeJacobianAndUpdateEKF(cv::Mat img);
  void resetGrid();
  void updateVisualMapPoints(cv::Mat img);
  void getWarpMatrixAffine(const vk::AbstractCamera &cam, const Vector2d &px_ref, const Vector3d &f_ref, const double depth_ref, const SE3 &T_cur_ref,
                           const int level_ref, 
                           const int pyramid_level, const int halfpatch_size, Matrix2d &A_cur_ref);
  void getWarpMatrixAffineHomography(const vk::AbstractCamera &cam, const V2D &px_ref,
                                     const V3D &xyz_ref, const V3D &normal_ref, const SE3 &T_cur_ref, const int level_ref, Matrix2d &A_cur_ref);
  void warpAffine(const Matrix2d &A_cur_ref, const cv::Mat &img_ref, const Vector2d &px_ref, const int level_ref, const int search_level,
                  const int pyramid_level, const int halfpatch_size, float *patch);
  void insertPointIntoVoxelMap(VisualPoint *pt_new);
  void plotTrackedPoints();
  void updateFrameState(StatesGroup state);
  void projectPatchFromRefToCur(
      const unordered_map<VOXEL_LOCATION, VoxelOctoTree *> &plane_map);

  /// \brief Perform V.D.Reference Patch Update
  /// \note choose one reference patch for image alignment in the visual update
  void updateReferencePatch(
      const unordered_map<VOXEL_LOCATION, VoxelOctoTree *> &plane_map);
  
  void precomputeReferencePatches(int level);
  void dumpDataForColmap();
  
  /// \brief Peform Normalized Cross-Correlation(NCC) to measure the similarity
  /// between two patch
  /// \note Refer to V.D.Reference Patch Update,formula(12)
  double calculateNCC(float *ref_patch, float *cur_patch, int patch_size);
  int getBestSearchLevel(const Matrix2d &A_cur_ref, const int max_level);
  V3F getInterpolatedPixel(cv::Mat img, V2D pc);

  /// \brief Peform Ray-Casting
  void PeformRayCasting(
      int loc_xyz[3],
      const unordered_map<VOXEL_LOCATION, VoxelOctoTree *> &plane_map);

  /// \brief Update Visual Map
  void UpdateVisualMap(cv::Mat img, float *image_data);

  /// \brief Get most similar patch as reference patch
  /// \note: similarity score : NCC + cosine
  void GetReferencePatch(VisualPoint *pt);
  
  // void resetRvizDisplay();
  // deque<VisualPoint *> map_cur_frame;
  // deque<VisualPoint *> sub_map_ray;
  // deque<VisualPoint *> sub_map_ray_fov;
  // deque<VisualPoint *> visual_sub_map_cur;
  // deque<VisualPoint *> visual_converged_point;
  // std::vector<std::vector<V3D>> sample_points;

  // PointCloudXYZI::Ptr pg_down;
  // pcl::VoxelGrid<PointType> downSizeFilter;
};
typedef std::shared_ptr<VIOManager> VIOManagerPtr;

#endif // VIO_H_