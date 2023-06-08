//
// Created by xiang on 23-1-19.
//

#include "ch4/g2o_types.h"
#include "common/g2o_types.h"

namespace sad
{

Vec3d r_Rij(Mat3d dR_ij, Mat3d R_i, Mat3d R_j)
{
    // dR_ij为预积分观测量
    // R_i和R_j分别为i和j时刻旋转估计值
    return SO3((dR_ij.inverse() * R_i.inverse() * R_j)).log();
}

Vec3d r_vij(Mat3d R_i, Vec3d v_j, Vec3d v_i, Vec3d g, double dt_ij, Vec3d dv_ij)
{
    // R_i, v_j, v_i为估计值
    // dv_ij为预积分观测值
    return R_i.inverse() * (v_j - v_i - g * dt_ij) - dv_ij;
}

Vec3d r_pij(Mat3d R_i, Vec3d p_j, Vec3d p_i, Vec3d v_i, double dt_ij, Vec3d g, Vec3d dp_ij)
{
    // R_i, p_j, p_i, v_i为估计值
    // dp_ij为预积分观测值
    return R_i.inverse() * (p_j - p_i - v_i * dt_ij - 0.5 * g * dt_ij * dt_ij) - dp_ij;
}

Mat3d ch4_skew(const Vec3d &v)
{
    Mat3d m;
    m << 0.0, -v[2], v[1], v[2], 0.0, -v[0], -v[1], v[0], 0.0;
    return m;
}

Mat3d ch4_Exp(const Vec3d &ang)
{
    double ang_norm = ang.norm();
    Mat3d Eye3 = Mat3d::Identity();
    if (ang_norm > 0.0000001)
    {
        Vec3d r_axis = ang / ang_norm;
        Mat3d K;
        K = ch4_skew(r_axis);
        /// Roderigous Tranformation
        return Eye3 + std::sin(ang_norm) * K + (1.0 - std::cos(ang_norm)) * K * K;
    }
    else
    {
        return Eye3;
    }
}

EdgeInertial::EdgeInertial(std::shared_ptr<IMUPreintegration> preinteg, const Vec3d &gravity, double weight)
    : preint_(preinteg), dt_(preinteg->dt_)
{
    resize(6); // 6个关联顶点
    grav_ = gravity;
    setInformation(preinteg->cov_.inverse() * weight);
}

void EdgeInertial::computeError()
{
    auto *p1 = dynamic_cast<const VertexPose *>(_vertices[0]);
    auto *v1 = dynamic_cast<const VertexVelocity *>(_vertices[1]);
    auto *bg1 = dynamic_cast<const VertexGyroBias *>(_vertices[2]);
    auto *ba1 = dynamic_cast<const VertexAccBias *>(_vertices[3]);
    auto *p2 = dynamic_cast<const VertexPose *>(_vertices[4]);
    auto *v2 = dynamic_cast<const VertexVelocity *>(_vertices[5]);

    Vec3d bg = bg1->estimate();
    Vec3d ba = ba1->estimate();

    const SO3 dR = preint_->GetDeltaRotation(bg);
    const Vec3d dv = preint_->GetDeltaVelocity(bg, ba);
    const Vec3d dp = preint_->GetDeltaPosition(bg, ba);

    /// 预积分误差项（4.41）
    const Vec3d er = (dR.inverse() * p1->estimate().so3().inverse() * p2->estimate().so3()).log();
    Mat3d RiT = p1->estimate().so3().inverse().matrix();
    const Vec3d ev = RiT * (v2->estimate() - v1->estimate() - grav_ * dt_) - dv;
    const Vec3d ep = RiT * (p2->estimate().translation() - p1->estimate().translation() - v1->estimate() * dt_ -
                            grav_ * dt_ * dt_ / 2) -
                     dp;
    _error << er, ev, ep;
}

void EdgeInertial::linearizeOplus()
{
    auto *p1 = dynamic_cast<const VertexPose *>(_vertices[0]);
    auto *v1 = dynamic_cast<const VertexVelocity *>(_vertices[1]);
    auto *bg1 = dynamic_cast<const VertexGyroBias *>(_vertices[2]);
    auto *ba1 = dynamic_cast<const VertexAccBias *>(_vertices[3]);
    auto *p2 = dynamic_cast<const VertexPose *>(_vertices[4]);
    auto *v2 = dynamic_cast<const VertexVelocity *>(_vertices[5]);

    Vec3d bg = bg1->estimate();
    Vec3d ba = ba1->estimate();
    Vec3d dbg = bg - preint_->bg_;

    // 一些中间符号
    const SO3 R1 = p1->estimate().so3();
    const SO3 R1T = R1.inverse();
    const SO3 R2 = p2->estimate().so3();

    auto dR_dbg = preint_->dR_dbg_;
    auto dv_dbg = preint_->dV_dbg_;
    auto dp_dbg = preint_->dP_dbg_;
    auto dv_dba = preint_->dV_dba_;
    auto dp_dba = preint_->dP_dba_;

    // 估计值
    Vec3d vi = v1->estimate();
    Vec3d vj = v2->estimate();
    Vec3d pi = p1->estimate().translation();
    Vec3d pj = p2->estimate().translation();

    const SO3 dR = preint_->GetDeltaRotation(bg);
    const SO3 eR = SO3(dR).inverse() * R1T * R2;
    const Vec3d er = eR.log();
    const Mat3d invJr = SO3::jr_inv(eR);
    const Vec3d dv = preint_->GetDeltaVelocity(bg, ba);
    const Vec3d ev = r_vij(R1.matrix(), vj, vi, grav_, dt_, dv);
    const Vec3d dp = preint_->GetDeltaPosition(bg, ba);
    const Vec3d ep = r_pij(R1.matrix(), pj, pi, vi, dt_, grav_, dp);

    /// 雅可比矩阵
    /// 注意有3个index, 顶点的，自己误差的，顶点内部变量的
    /// 变量顺序：pose1(R1,p1), v1, bg1, ba1, pose2(R2,p2), v2
    /// 残差顺序：eR, ev, ep，残差顺序为行，变量顺序为列

    //       | R1 | p1 | v1 | bg1 | ba1 | R2 | p2 | v2 |
    //  vert | 0       | 1  | 2   | 3   | 4       | 5  |
    //  col  | 0    3  | 0  | 0   | 0   | 0    3  | 0  |
    //    row
    //  eR 0 |
    //  ev 3 |
    //  ep 6 |

    /************************************* 解析求导jacobian ****************************************/
    std::cout << "************** 解析求导 **************** \n" << std::endl;
    /// 残差对R1, 9x3
    _jacobianOplus[0].setZero();
    // dR/dR1, 4.42
    _jacobianOplus[0].block<3, 3>(0, 0) = -invJr * (R2.inverse() * R1).matrix();
    std::cout << "dR/dR1: \n" << _jacobianOplus[0].block<3, 3>(0, 0) << std::endl;
    // dv/dR1, 4.47
    _jacobianOplus[0].block<3, 3>(3, 0) = SO3::hat(R1T * (vj - vi - grav_ * dt_));
    std::cout << "dv/dR1: \n" << _jacobianOplus[0].block<3, 3>(3, 0) << std::endl;
    // dp/dR1, 4.48d
    _jacobianOplus[0].block<3, 3>(6, 0) = SO3::hat(R1T * (pj - pi - v1->estimate() * dt_ - 0.5 * grav_ * dt_ * dt_));
    std::cout << "dp/dR1: \n" << _jacobianOplus[0].block<3, 3>(6, 0) << std::endl;

    /// 残差对p1, 9x3
    // dp/dp1, 4.48a
    _jacobianOplus[0].block<3, 3>(6, 3) = -R1T.matrix();
    std::cout << "dp/dp1: \n" << _jacobianOplus[0].block<3, 3>(6, 3) << std::endl;

    /// 残差对v1, 9x3
    _jacobianOplus[1].setZero();
    // dv/dv1, 4.46a
    _jacobianOplus[1].block<3, 3>(3, 0) = -R1T.matrix();
    std::cout << "dv/dv1: \n" << _jacobianOplus[1].block<3, 3>(3, 0) << std::endl;
    // dp/dv1, 4.48c
    _jacobianOplus[1].block<3, 3>(6, 0) = -R1T.matrix() * dt_;
    std::cout << "dp/dv1: \n" << _jacobianOplus[1].block<3, 3>(6, 0) << std::endl;

    /// 残差对bg1
    _jacobianOplus[2].setZero();
    // dR/dbg1, 4.45
    _jacobianOplus[2].block<3, 3>(0, 0) = -invJr * eR.inverse().matrix() * SO3::jr((dR_dbg * dbg).eval()) * dR_dbg;
    std::cout << "dR/dbg1: \n" << _jacobianOplus[2].block<3, 3>(0, 0) << std::endl;
    // dv/dbg1
    _jacobianOplus[2].block<3, 3>(3, 0) = -dv_dbg;
    std::cout << "dv/dbg1: \n" << _jacobianOplus[2].block<3, 3>(3, 0) << std::endl;
    // dp/dbg1
    _jacobianOplus[2].block<3, 3>(6, 0) = -dp_dbg;
    std::cout << "dp/dbg1: \n" << _jacobianOplus[2].block<3, 3>(6, 0) << std::endl;

    /// 残差对ba1
    _jacobianOplus[3].setZero();
    // dv/dba1
    _jacobianOplus[3].block<3, 3>(3, 0) = -dv_dba;
    std::cout << "dv/dba1: \n" << _jacobianOplus[3].block<3, 3>(3, 0) << std::endl;
    // dp/dba1
    _jacobianOplus[3].block<3, 3>(6, 0) = -dp_dba;
    std::cout << "dp/dba1: \n" << _jacobianOplus[3].block<3, 3>(6, 0) << std::endl;

    /// 残差对pose2
    _jacobianOplus[4].setZero();
    // dR/dR2, 4.43
    _jacobianOplus[4].block<3, 3>(0, 0) = invJr;
    std::cout << "dR/dR2: \n" << _jacobianOplus[4].block<3, 3>(0, 0) << std::endl;
    // dp/dp2, 4.48b
    _jacobianOplus[4].block<3, 3>(6, 3) = R1T.matrix();
    std::cout << "dp/dp2: \n" << _jacobianOplus[4].block<3, 3>(6, 3) << std::endl;

    /// 残差对v2
    _jacobianOplus[5].setZero();
    // dv/dv2, 4,46b
    _jacobianOplus[5].block<3, 3>(3, 0) = R1T.matrix(); // OK
    std::cout << "dv/dv2: \n" << _jacobianOplus[5].block<3, 3>(3, 0) << std::endl;

    /************************************* 数值求导jacobian *****************************************/
    std::cout << "************** 数值求导 **************** \n" << std::endl;
    double eps = 1e-6;
    Vec3d inc0(eps, 0.0, 0.0);
    Vec3d inc1(0.0, eps, 0.0);
    Vec3d inc2(0.0, 0.0, eps);

    // dR/dR1, 4.42
    Mat3d dR_dR1;
    Vec3d Delta_er = r_Rij(dR.matrix(), R1.matrix() * ch4_Exp(inc0), R2.matrix());
    dR_dR1.block<3, 1>(0, 0) = (Delta_er - er) / eps;
    Delta_er = r_Rij(dR.matrix(), R1.matrix() * ch4_Exp(inc1), R2.matrix());
    dR_dR1.block<3, 1>(0, 1) = (Delta_er - er) / eps;
    Delta_er = r_Rij(dR.matrix(), R1.matrix() * ch4_Exp(inc2), R2.matrix());
    dR_dR1.block<3, 1>(0, 2) = (Delta_er - er) / eps;
    std::cout << "dR/dR1: \n" << dR_dR1 << std::endl;

    // dv/dR1, 4.47
    Mat3d dv_dR1;
    Vec3d Delta_ev = r_vij(R1.matrix() * ch4_Exp(inc0), vj, vi, grav_, dt_, dv);
    dv_dR1.block<3, 1>(0, 0) = (Delta_ev - ev) / eps;
    Delta_ev = r_vij(R1.matrix() * ch4_Exp(inc1), vj, vi, grav_, dt_, dv);
    dv_dR1.block<3, 1>(0, 1) = (Delta_ev - ev) / eps;
    Delta_ev = r_vij(R1.matrix() * ch4_Exp(inc2), vj, vi, grav_, dt_, dv);
    dv_dR1.block<3, 1>(0, 2) = (Delta_ev - ev) / eps;
    std::cout << "dv/dR1: \n" << dv_dR1 << std::endl;

    // dp/dR1, 4.48d
    Mat3d dp_dR1;
    Vec3d Delta_ep = r_pij(R1.matrix() * ch4_Exp(inc0), pj, pi, vi, dt_, grav_, dp);
    dp_dR1.block<3, 1>(0, 0) = (Delta_ep - ep) / eps;
    Delta_ep = r_pij(R1.matrix() * ch4_Exp(inc1), pj, pi, vi, dt_, grav_, dp);
    dp_dR1.block<3, 1>(0, 1) = (Delta_ep - ep) / eps;
    Delta_ep = r_pij(R1.matrix() * ch4_Exp(inc2), pj, pi, vi, dt_, grav_, dp);
    dp_dR1.block<3, 1>(0, 2) = (Delta_ep - ep) / eps;
    std::cout << "dp/dR1: \n" << dp_dR1 << std::endl;

    // dp/dp1, 4.48a
    Mat3d dp_dp1;
    Delta_ep = r_pij(R1.matrix(), pj, pi + inc0, vi, dt_, grav_, dp);
    dp_dp1.block<3, 1>(0, 0) = (Delta_ep - ep) / eps;
    Delta_ep = r_pij(R1.matrix(), pj, pi + inc1, vi, dt_, grav_, dp);
    dp_dp1.block<3, 1>(0, 1) = (Delta_ep - ep) / eps;
    Delta_ep = r_pij(R1.matrix(), pj, pi + inc2, vi, dt_, grav_, dp);
    dp_dp1.block<3, 1>(0, 2) = (Delta_ep - ep) / eps;
    std::cout << "dp/dp1: \n" << dp_dp1 << std::endl;

    // dv/dv1, 4.46a
    Mat3d dv_dv1;
    Delta_ev = r_vij(R1.matrix(), vj, vi + inc0, grav_, dt_, dv);
    dv_dv1.block<3, 1>(0, 0) = (Delta_ev - ev) / eps;
    Delta_ev = r_vij(R1.matrix(), vj, vi + inc1, grav_, dt_, dv);
    dv_dv1.block<3, 1>(0, 1) = (Delta_ev - ev) / eps;
    Delta_ev = r_vij(R1.matrix(), vj, vi + inc2, grav_, dt_, dv);
    dv_dv1.block<3, 1>(0, 2) = (Delta_ev - ev) / eps;
    std::cout << "dv/dv1: \n" << dv_dv1 << std::endl;

    // dp/dv1, 4.48c
    Mat3d dp_dv1;
    Delta_ep = r_pij(R1.matrix(), pj, pi, vi + inc0, dt_, grav_, dp);
    dp_dv1.block<3, 1>(0, 0) = (Delta_ep - ep) / eps;
    Delta_ep = r_pij(R1.matrix(), pj, pi, vi + inc1, dt_, grav_, dp);
    dp_dv1.block<3, 1>(0, 1) = (Delta_ep - ep) / eps;
    Delta_ep = r_pij(R1.matrix(), pj, pi, vi + inc2, dt_, grav_, dp);
    dp_dv1.block<3, 1>(0, 2) = (Delta_ep - ep) / eps;
    std::cout << "dp/dv1: \n" << dp_dv1 << std::endl;

    // dR/dbg1, 4.45
    Mat3d dR_dbg1;
    Delta_er = r_Rij(preint_->GetDeltaRotation(bg + inc0).matrix(), R1.matrix(), R2.matrix());
    dR_dbg1.block<3, 1>(0, 0) = (Delta_er - er) / eps;
    Delta_er = r_Rij(preint_->GetDeltaRotation(bg + inc1).matrix(), R1.matrix(), R2.matrix());
    dR_dbg1.block<3, 1>(0, 1) = (Delta_er - er) / eps;
    Delta_er = r_Rij(preint_->GetDeltaRotation(bg + inc2).matrix(), R1.matrix(), R2.matrix());
    dR_dbg1.block<3, 1>(0, 2) = (Delta_er - er) / eps;
    std::cout << "dR/dbg1: \n" << dR_dbg1 << std::endl;

    // dv/dbg1
    Mat3d dv_dbg1;
    Delta_ev = r_vij(R1.matrix(), vj, vi, grav_, dt_, preint_->GetDeltaVelocity(bg + inc0, ba));
    dv_dbg1.block<3, 1>(0, 0) = (Delta_ev - ev) / eps;
    Delta_ev = r_vij(R1.matrix(), vj, vi, grav_, dt_, preint_->GetDeltaVelocity(bg + inc1, ba));
    dv_dbg1.block<3, 1>(0, 1) = (Delta_ev - ev) / eps;
    Delta_ev = r_vij(R1.matrix(), vj, vi, grav_, dt_, preint_->GetDeltaVelocity(bg + inc2, ba));
    dv_dbg1.block<3, 1>(0, 2) = (Delta_ev - ev) / eps;
    std::cout << "dv/dbg1: \n" << dv_dbg1 << std::endl;

    // dp/dbg1
    Mat3d dp_dbg1;
    Delta_ep = r_pij(R1.matrix(), pj, pi, vi, dt_, grav_, preint_->GetDeltaPosition(bg + inc0, ba));
    dp_dbg1.block<3, 1>(0, 0) = (Delta_ep - ep) / eps;
    Delta_ep = r_pij(R1.matrix(), pj, pi, vi, dt_, grav_, preint_->GetDeltaPosition(bg + inc1, ba));
    dp_dbg1.block<3, 1>(0, 1) = (Delta_ep - ep) / eps;
    Delta_ep = r_pij(R1.matrix(), pj, pi, vi, dt_, grav_, preint_->GetDeltaPosition(bg + inc2, ba));
    dp_dbg1.block<3, 1>(0, 2) = (Delta_ep - ep) / eps;
    std::cout << "dp/dbg1: \n" << dp_dbg1 << std::endl;

    // dv/dba1
    Mat3d dv_dba1;
    Delta_ev = r_vij(R1.matrix(), vj, vi, grav_, dt_, preint_->GetDeltaVelocity(bg, ba + inc0));
    dv_dba1.block<3, 1>(0, 0) = (Delta_ev - ev) / eps;
    Delta_ev = r_vij(R1.matrix(), vj, vi, grav_, dt_, preint_->GetDeltaVelocity(bg, ba + inc1));
    dv_dba1.block<3, 1>(0, 1) = (Delta_ev - ev) / eps;
    Delta_ev = r_vij(R1.matrix(), vj, vi, grav_, dt_, preint_->GetDeltaVelocity(bg, ba + inc2));
    dv_dba1.block<3, 1>(0, 2) = (Delta_ev - ev) / eps;
    std::cout << "dv/dba1: \n" << dv_dba1 << std::endl;

    // dp/dba1
    Mat3d dp_dba1;
    Delta_ep = r_pij(R1.matrix(), pj, pi, vi, dt_, grav_, preint_->GetDeltaPosition(bg, ba + inc0));
    dp_dba1.block<3, 1>(0, 0) = (Delta_ep - ep) / eps;
    Delta_ep = r_pij(R1.matrix(), pj, pi, vi, dt_, grav_, preint_->GetDeltaPosition(bg, ba + inc1));
    dp_dba1.block<3, 1>(0, 1) = (Delta_ep - ep) / eps;
    Delta_ep = r_pij(R1.matrix(), pj, pi, vi, dt_, grav_, preint_->GetDeltaPosition(bg, ba + inc2));
    dp_dba1.block<3, 1>(0, 2) = (Delta_ep - ep) / eps;
    std::cout << "dp/dba1: \n" << dp_dba1 << std::endl;

    // dR/dR2, 4.43
    Mat3d dR_dR2;
    Delta_er = r_Rij(dR.matrix(), R1.matrix(), R2.matrix() * ch4_Exp(inc0));
    dR_dR2.block<3, 1>(0, 0) = (Delta_er - er) / eps;
    Delta_er = r_Rij(dR.matrix(), R1.matrix(), R2.matrix() * ch4_Exp(inc1));
    dR_dR2.block<3, 1>(0, 1) = (Delta_er - er) / eps;
    Delta_er = r_Rij(dR.matrix(), R1.matrix(), R2.matrix() * ch4_Exp(inc2));
    dR_dR2.block<3, 1>(0, 2) = (Delta_er - er) / eps;
    std::cout << "dR/dR2: \n" << dR_dR2 << std::endl;

    // dp/dp2, 4.48b
    Mat3d dp_dp2;
    Delta_ep = r_pij(R1.matrix(), pj + inc0, pi, vi, dt_, grav_, dp);
    dp_dp2.block<3, 1>(0, 0) = (Delta_ep - ep) / eps;
    Delta_ep = r_pij(R1.matrix(), pj + inc1, pi, vi, dt_, grav_, dp);
    dp_dp2.block<3, 1>(0, 1) = (Delta_ep - ep) / eps;
    Delta_ep = r_pij(R1.matrix(), pj + inc2, pi, vi, dt_, grav_, dp);
    dp_dp2.block<3, 1>(0, 2) = (Delta_ep - ep) / eps;
    std::cout << "dp/dp2: \n" << dp_dp2 << std::endl;

    // dv/dv2, 4,46b
    Mat3d dv_dv2;
    Delta_ev = r_vij(R1.matrix(), vj + inc0, vi, grav_, dt_, dv);
    dv_dv2.block<3, 1>(0, 0) = (Delta_ev - ev) / eps;
    Delta_ev = r_vij(R1.matrix(), vj + inc1, vi, grav_, dt_, dv);
    dv_dv2.block<3, 1>(0, 1) = (Delta_ev - ev) / eps;
    Delta_ev = r_vij(R1.matrix(), vj + inc2, vi, grav_, dt_, dv);
    dv_dv2.block<3, 1>(0, 2) = (Delta_ev - ev) / eps;
    std::cout << "dv/dv2: \n" << dv_dv2 << std::endl;
}

} // namespace sad