//
// Created by xiang on 2022/7/18.
//

#ifndef SLAM_IN_AUTO_DRIVING_DIRECT_ICP_LO_H
#define SLAM_IN_AUTO_DRIVING_DIRECT_ICP_LO_H

#include <deque>

#include "ch7/icp_3d.h"
#include "common/eigen_types.h"
#include "common/point_types.h"
#include "tools/pcl_map_viewer.h"

namespace sad
{

/**
 * 使用直接ICP方法进行递推的Lidar Odometry
 * 使用历史几个关键帧作为local map，进行ICP定位
 */
class DirectICPLO
{
  public:
    struct Options
    {
        Options()
        {
        }
        double kf_distance_ = 0.5;           // 关键帧距离
        double kf_angle_deg_ = 30;           // 旋转角度
        int num_kfs_in_local_map_ = 30;      // 局部地图含有多少个关键帧
        bool display_realtime_cloud_ = true; // 是否显示实时点云

        Icp3d::Options icp3d_options_; // ICP3D 的配置
    };

    DirectICPLO(Options options = Options()) : options_(options)
    {
        if (options_.display_realtime_cloud_)
        {
            viewer_ = std::make_shared<PCLMapViewer>(0.5);
        }

        icp_ = Icp3d(options_.icp3d_options_);
    }

    /**
     * 往LO中增加一个点云
     * @param scan  当前帧点云
     * @param pose 估计pose
     */
    void AddCloud(CloudPtr scan, SE3 &pose);

    /// 存储地图(viewer里）
    void SaveMap(const std::string &map_path);

  private:
    /// 与local map进行配准
    SE3 AlignWithLocalMap(CloudPtr scan);

    /// 判定是否为关键帧
    bool IsKeyframe(const SE3 &current_pose);

  private:
    Options options_;
    CloudPtr local_map_ = nullptr;
    std::deque<CloudPtr> scans_in_local_map_;
    std::vector<SE3> estimated_poses_; // 所有估计出来的pose，用于记录轨迹和预测下一个帧
    SE3 last_kf_pose_;                 // 上一关键帧的位姿

    Icp3d icp_;

    std::shared_ptr<PCLMapViewer> viewer_ = nullptr;
};

} // namespace sad

#endif // SLAM_IN_AUTO_DRIVING_DIRECT_NDT_LO_H
