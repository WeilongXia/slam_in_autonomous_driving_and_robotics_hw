//
// Created by xiang on 2022/7/7.
//

#include "icp_3d.h"
#include "common/math_utils.h"

#include <execution>

namespace sad
{

void Icp3d::AddCloud(CloudPtr cloud_world)
{
    if (local_map_ == nullptr)
    {
        // 第一个帧，直接加入local map
        local_map_.reset(new PointCloudType);
        // operator += 用来拼接点云
        *local_map_ += *cloud_world;
        return;
    }
    else
    {
        // 重建local map
        scans_in_local_map_.emplace_back(cloud_world);
        if (scans_in_local_map_.size() > 30)
        {
            scans_in_local_map_.pop_front();
        }

        local_map_.reset(new PointCloudType);
        for (auto &scan : scans_in_local_map_)
        {
            *local_map_ += *scan;
        }
    }
}

bool Icp3d::AlignP2P(SE3 &init_pose)
{
    LOG(INFO) << "aligning with point to point";
    assert(target_ != nullptr && source_ != nullptr);

    SE3 pose = init_pose;
    if (!options_.use_initial_translation_)
    {
        pose.translation() = target_center_ - source_center_; // 设置平移初始值
    }

    // 对点的索引，预先生成
    std::vector<int> index(source_->points.size());
    for (int i = 0; i < index.size(); ++i)
    {
        index[i] = i;
    }

    // 我们来写一些并发代码
    std::vector<bool> effect_pts(index.size(), false);
    std::vector<Eigen::Matrix<double, 3, 6>> jacobians(index.size());
    std::vector<Vec3d> errors(index.size());

    for (int iter = 0; iter < options_.max_iteration_; ++iter)
    {
        // gauss-newton 迭代
        // 最近邻，可以并发
        std::for_each(std::execution::par_unseq, index.begin(), index.end(), [&](int idx) {
            auto q = ToVec3d(source_->points[idx]);
            Vec3d qs = pose * q; // 转换之后的q
            std::vector<int> nn;
            kdtree_->GetClosestPoint(ToPointType(qs), nn, 1);

            if (!nn.empty())
            {
                Vec3d p = ToVec3d(target_->points[nn[0]]);
                double dis2 = (p - qs).squaredNorm();
                if (dis2 > options_.max_nn_distance_)
                {
                    // 点离的太远了不要
                    effect_pts[idx] = false;
                    return;
                }

                effect_pts[idx] = true;

                // build residual
                Vec3d e = p - qs;
                double weight = CauchyLoss(e.norm());

                Eigen::Matrix<double, 3, 6> J;
                J.block<3, 3>(0, 0) = weight * pose.so3().matrix() * SO3::hat(q);
                J.block<3, 3>(0, 3) = -weight * Mat3d::Identity();

                jacobians[idx] = J;
                errors[idx] = weight * e;
            }
            else
            {
                effect_pts[idx] = false;
            }
        });

        // 累加Hessian和error,计算dx
        // 原则上可以用reduce并发，写起来比较麻烦，这里写成accumulate
        double total_res = 0;
        int effective_num = 0;
        auto H_and_err = std::accumulate(
            index.begin(), index.end(), std::pair<Mat6d, Vec6d>(Mat6d::Zero(), Vec6d::Zero()),
            [&jacobians, &errors, &effect_pts, &total_res, &effective_num](const std::pair<Mat6d, Vec6d> &pre,
                                                                           int idx) -> std::pair<Mat6d, Vec6d> {
                if (!effect_pts[idx])
                {
                    return pre;
                }
                else
                {
                    total_res += errors[idx].dot(errors[idx]);
                    effective_num++;
                    return std::pair<Mat6d, Vec6d>(pre.first + jacobians[idx].transpose() * jacobians[idx],
                                                   pre.second - jacobians[idx].transpose() * errors[idx]);
                }
            });

        if (effective_num < options_.min_effective_pts_)
        {
            LOG(WARNING) << "effective num too small: " << effective_num;
            return false;
        }

        Mat6d H = H_and_err.first;
        Vec6d err = H_and_err.second;

        Vec6d dx = H.inverse() * err;
        pose.so3() = pose.so3() * SO3::exp(dx.head<3>());
        pose.translation() += dx.tail<3>();

        // 更新
        LOG(INFO) << "iter " << iter << " total res: " << total_res << ", eff: " << effective_num
                  << ", mean res: " << total_res / effective_num << ", dxn: " << dx.norm();

        if (gt_set_)
        {
            double pose_error = (gt_pose_.inverse() * pose).log().norm();
            LOG(INFO) << "iter " << iter << " pose error: " << pose_error;
        }

        if (dx.norm() < options_.eps_)
        {
            LOG(INFO) << "converged, dx = " << dx.transpose();
            break;
        }
    }

    init_pose = pose;
    return true;
}

bool Icp3d::AlignP2Plane(SE3 &init_pose)
{
    LOG(INFO) << "aligning with point to plane";
    assert(target_ != nullptr && source_ != nullptr);
    // 整体流程与p2p一致，读者请关注变化部分

    SE3 pose = init_pose;
    if (!options_.use_initial_translation_)
    {
        pose.translation() = target_center_ - source_center_; // 设置平移初始值
    }

    std::vector<int> index(source_->points.size());
    for (int i = 0; i < index.size(); ++i)
    {
        index[i] = i;
    }

    std::vector<bool> effect_pts(index.size(), false);
    std::vector<Eigen::Matrix<double, 1, 6>> jacobians(index.size());
    std::vector<double> errors(index.size());

    for (int iter = 0; iter < options_.max_iteration_; ++iter)
    {
        // gauss-newton 迭代
        // 最近邻，可以并发
        std::for_each(std::execution::par_unseq, index.begin(), index.end(), [&](int idx) {
            auto q = ToVec3d(source_->points[idx]);
            Vec3d qs = pose * q; // 转换之后的q
            std::vector<int> nn;
            kdtree_->GetClosestPoint(ToPointType(qs), nn, 5); // 这里取5个最近邻
            if (nn.size() > 3)
            {
                // convert to eigen
                std::vector<Vec3d> nn_eigen;
                for (int i = 0; i < nn.size(); ++i)
                {
                    nn_eigen.emplace_back(ToVec3d(target_->points[nn[i]]));
                }

                Vec4d n;
                if (!math::FitPlane(nn_eigen, n))
                {
                    // 失败的不要
                    effect_pts[idx] = false;
                    return;
                }

                double dis = n.head<3>().dot(qs) + n[3];
                if (fabs(dis) > options_.max_plane_distance_)
                {
                    // 点离的太远了不要
                    effect_pts[idx] = false;
                    return;
                }

                // double weight = CauchyLoss(dis);
                double weight = 1.0;

                effect_pts[idx] = true;

                // build residual
                Eigen::Matrix<double, 1, 6> J;
                J.block<1, 3>(0, 0) = -weight * n.head<3>().transpose() * pose.so3().matrix() * SO3::hat(q);
                J.block<1, 3>(0, 3) = weight * n.head<3>().transpose();

                jacobians[idx] = J;
                errors[idx] = weight * dis;
            }
            else
            {
                effect_pts[idx] = false;
            }
        });

        // 累加Hessian和error,计算dx
        // 原则上可以用reduce并发，写起来比较麻烦，这里写成accumulate
        double total_res = 0;
        int effective_num = 0;
        auto H_and_err = std::accumulate(
            index.begin(), index.end(), std::pair<Mat6d, Vec6d>(Mat6d::Zero(), Vec6d::Zero()),
            [&jacobians, &errors, &effect_pts, &total_res, &effective_num](const std::pair<Mat6d, Vec6d> &pre,
                                                                           int idx) -> std::pair<Mat6d, Vec6d> {
                if (!effect_pts[idx])
                {
                    return pre;
                }
                else
                {
                    total_res += errors[idx] * errors[idx];
                    effective_num++;
                    return std::pair<Mat6d, Vec6d>(pre.first + jacobians[idx].transpose() * jacobians[idx],
                                                   pre.second - jacobians[idx].transpose() * errors[idx]);
                }
            });

        if (effective_num < options_.min_effective_pts_)
        {
            LOG(WARNING) << "effective num too small: " << effective_num;
            return false;
        }

        Mat6d H = H_and_err.first;
        Vec6d err = H_and_err.second;

        Vec6d dx = H.inverse() * err;
        pose.so3() = pose.so3() * SO3::exp(dx.head<3>());
        pose.translation() += dx.tail<3>();

        // 更新
        LOG(INFO) << "iter " << iter << " total res: " << total_res << ", eff: " << effective_num
                  << ", mean res: " << total_res / effective_num << ", dxn: " << dx.norm();

        if (gt_set_)
        {
            double pose_error = (gt_pose_.inverse() * pose).log().norm();
            LOG(INFO) << "iter " << iter << " pose error: " << pose_error;
        }

        if (dx.norm() < options_.eps_)
        {
            LOG(INFO) << "converged, dx = " << dx.transpose();
            break;
        }
    }

    init_pose = pose;
    return true;
}

void Icp3d::ComputeResidualAndJacobians(const SE3 &input_pose, Mat18d &HTVH, Vec18d &HTVr)
{
    // LOG(INFO) << "compute residual and jacobians";
    assert(local_map_ != nullptr && source_ != nullptr);

    // 整体流程与AlignP2Plane一致，不需要迭代更新位姿，只需要计算HTVH和HTVr

    SE3 pose = input_pose;
    if (!options_.use_initial_translation_)
    {
        pose.translation() = target_center_ - source_center_; // 设置平移初始值
    }

    std::vector<int> index(source_->points.size());
    for (int i = 0; i < index.size(); ++i)
    {
        index[i] = i;
    }

    std::vector<bool> effect_pts(index.size(), false);
    std::vector<Eigen::Matrix<double, 1, 18>> jacobians(index.size());
    std::vector<double> errors(index.size());

    // 最近邻，可以并发
    std::for_each(std::execution::par_unseq, index.begin(), index.end(), [&](int idx) {
        auto q = ToVec3d(source_->points[idx]);
        Vec3d qs = pose * q; // 转换之后的q
        std::vector<int> nn;
        localmap_kdtree_->GetClosestPoint(ToPointType(qs), nn, 5); // 这里取5个最近邻
        if (nn.size() > 3)
        {
            // convert to eigen
            std::vector<Vec3d> nn_eigen;
            for (int i = 0; i < nn.size(); ++i)
            {
                nn_eigen.emplace_back(ToVec3d(local_map_->points[nn[i]]));
            }

            Vec4d n;
            if (!math::FitPlane(nn_eigen, n))
            {
                // 失败的不要
                effect_pts[idx] = false;
                return;
            }

            double dis = n.head<3>().dot(qs) + n[3];
            if (fabs(dis) > options_.max_plane_distance_)
            {
                // 点离的太远了不要
                effect_pts[idx] = false;
                return;
            }

            // double weight = CauchyLoss(dis);
            double weight = 1.0;

            effect_pts[idx] = true;

            // calculate jacobians
            Eigen::Matrix<double, 1, 18> J = Eigen::Matrix<double, 1, 18>::Zero();
            J.block<1, 3>(0, 0) = weight * n.head<3>().transpose();
            J.block<1, 3>(0, 6) = -weight * n.head<3>().transpose() * pose.so3().matrix() * SO3::hat(q);

            jacobians[idx] = J;
            errors[idx] = weight * dis;
            effect_pts[idx] = true;
        }
        else
        {
            effect_pts[idx] = false;
        }
    });

    // 累加Hessian和error，计算dx
    double total_res = 0;
    int effective_num = 0;

    HTVH.setZero();
    HTVr.setZero();

    const double info_ratio = 10; // 每个点的反馈因子

    for (int idx = 0; idx < effect_pts.size(); ++idx)
    {
        if (!effect_pts[idx])
        {
            continue;
        }

        total_res += errors[idx] * errors[idx];
        effective_num++;

        HTVH += jacobians[idx].transpose() * jacobians[idx] * info_ratio;
        HTVr += -jacobians[idx].transpose() * errors[idx] * info_ratio;
    }

    LOG(INFO) << "effective: " << effective_num << std::endl;
}

void Icp3d::BuildTargetKdTree()
{
    kdtree_ = std::make_shared<KdTree>();
    kdtree_->BuildTree(target_);
    kdtree_->SetEnableANN();
}

void Icp3d::BuildLocalMapKdTree()
{
    localmap_kdtree_ = std::make_shared<KdTree>();
    localmap_kdtree_->BuildTree(local_map_);
    localmap_kdtree_->SetEnableANN();
}

double Icp3d::CauchyLoss(double residual, double c)
{
    return c * c * std::log(1.0 + (residual / c) * (residual / c));
}

bool Icp3d::AlignP2Line(SE3 &init_pose)
{
    LOG(INFO) << "aligning with point to line";
    assert(target_ != nullptr && source_ != nullptr);
    // 点线与点面基本是完全一样的

    SE3 pose = init_pose;
    pose.translation() = target_center_ - source_center_; // 设置平移初始值
    LOG(INFO) << "init trans set to " << pose.translation().transpose();

    std::vector<int> index(source_->points.size());
    for (int i = 0; i < index.size(); ++i)
    {
        index[i] = i;
    }

    std::vector<bool> effect_pts(index.size(), false);
    std::vector<Eigen::Matrix<double, 3, 6>> jacobians(index.size());
    std::vector<Vec3d> errors(index.size());

    for (int iter = 0; iter < options_.max_iteration_; ++iter)
    {
        // gauss-newton 迭代
        // 最近邻，可以并发
        std::for_each(std::execution::par_unseq, index.begin(), index.end(), [&](int idx) {
            auto q = ToVec3d(source_->points[idx]);
            Vec3d qs = pose * q; // 转换之后的q
            std::vector<int> nn;
            kdtree_->GetClosestPoint(ToPointType(qs), nn, 5); // 这里取5个最近邻
            if (nn.size() == 5)
            {
                // convert to eigen
                std::vector<Vec3d> nn_eigen;
                for (int i = 0; i < 5; ++i)
                {
                    nn_eigen.emplace_back(ToVec3d(target_->points[nn[i]]));
                }

                Vec3d d, p0;
                if (!math::FitLine(nn_eigen, p0, d, options_.max_line_distance_))
                {
                    // 失败的不要
                    effect_pts[idx] = false;
                    return;
                }

                Vec3d err = SO3::hat(d) * (qs - p0);

                if (err.norm() > options_.max_line_distance_)
                {
                    // 点离的太远了不要
                    effect_pts[idx] = false;
                    return;
                }

                double weight = CauchyLoss(err.norm());

                effect_pts[idx] = true;

                // build residual
                Eigen::Matrix<double, 3, 6> J;
                J.block<3, 3>(0, 0) = -weight * SO3::hat(d) * pose.so3().matrix() * SO3::hat(q);
                J.block<3, 3>(0, 3) = weight * SO3::hat(d);

                jacobians[idx] = J;
                errors[idx] = weight * err;
            }
            else
            {
                effect_pts[idx] = false;
            }
        });

        // 累加Hessian和error,计算dx
        // 原则上可以用reduce并发，写起来比较麻烦，这里写成accumulate
        double total_res = 0;
        int effective_num = 0;
        auto H_and_err = std::accumulate(
            index.begin(), index.end(), std::pair<Mat6d, Vec6d>(Mat6d::Zero(), Vec6d::Zero()),
            [&jacobians, &errors, &effect_pts, &total_res, &effective_num](const std::pair<Mat6d, Vec6d> &pre,
                                                                           int idx) -> std::pair<Mat6d, Vec6d> {
                if (!effect_pts[idx])
                {
                    return pre;
                }
                else
                {
                    total_res += errors[idx].dot(errors[idx]);
                    effective_num++;
                    return std::pair<Mat6d, Vec6d>(pre.first + jacobians[idx].transpose() * jacobians[idx],
                                                   pre.second - jacobians[idx].transpose() * errors[idx]);
                }
            });

        if (effective_num < options_.min_effective_pts_)
        {
            LOG(WARNING) << "effective num too small: " << effective_num;
            return false;
        }

        Mat6d H = H_and_err.first;
        Vec6d err = H_and_err.second;

        Vec6d dx = H.inverse() * err;
        pose.so3() = pose.so3() * SO3::exp(dx.head<3>());
        pose.translation() += dx.tail<3>();

        if (gt_set_)
        {
            double pose_error = (gt_pose_.inverse() * pose).log().norm();
            LOG(INFO) << "iter " << iter << " pose error: " << pose_error;
        }

        // 更新
        LOG(INFO) << "iter " << iter << " total res: " << total_res << ", eff: " << effective_num
                  << ", mean res: " << total_res / effective_num << ", dxn: " << dx.norm();

        if (dx.norm() < options_.eps_)
        {
            LOG(INFO) << "converged, dx = " << dx.transpose();
            break;
        }
    }

    init_pose = pose;
    return true;
}

} // namespace sad