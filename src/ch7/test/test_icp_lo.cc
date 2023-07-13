//
// Created by xiang on 2022/7/18.
//

#include <gflags/gflags.h>
#include <glog/logging.h>

#include "ch7/direct_icp_lo.h"
#include "ch7/icp_3d.h"
#include "common/dataset_type.h"
#include "common/io_utils.h"
#include "common/timer/timer.h"

/// 本程序以ULHK数据集为例
/// 测试以点面ICP为主的Lidar Odometry
DEFINE_string(bag_path, "../dataset/ulhk/test2.bag", "path to rosbag");
DEFINE_string(dataset_type, "ULHK", "NCLT/ULHK/KITTI/WXB_3D"); // 数据集类型
DEFINE_bool(display_map, true, "display map?");

int main(int argc, char **argv)
{
    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = google::INFO;
    FLAGS_colorlogtostderr = true;
    google::ParseCommandLineFlags(&argc, &argv, true);

    sad::RosbagIO rosbag_io(fLS::FLAGS_bag_path, sad::Str2DatasetType(FLAGS_dataset_type));

    sad::DirectICPLO::Options options;
    options.display_realtime_cloud_ = FLAGS_display_map;
    sad::DirectICPLO icp_lo(options);

    rosbag_io
        .AddAutoPointCloudHandle([&icp_lo](sensor_msgs::PointCloud2::Ptr msg) -> bool {
            sad::common::Timer::Evaluate(
                [&]() {
                    SE3 pose;
                    icp_lo.AddCloud(sad::VoxelCloud(sad::PointCloud2ToCloudPtr(msg), 0.5), pose);
                },
                "ICP registration");
            return true;
        })
        .Go();

    if (FLAGS_display_map)
    {
        // 把地图存下来
        icp_lo.SaveMap("./data/ch7/map.pcd");
    }

    sad::common::Timer::PrintAll();
    LOG(INFO) << "done.";

    return 0;
}
