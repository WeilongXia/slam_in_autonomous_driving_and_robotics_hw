#include "ch5/nanoflann.hpp"
#include "common/point_cloud_utils.h"
#include "common/point_types.h"
#include "common/sys_utils.h"
#include <vector>

// 定义数据点结构
// struct Point
// {
//     double x, y, z;

//     Point(double x, double y, double z) : x(x), y(y), z(z)
//     {
//     }
// };

// 定义数据集适配器
class NanoPointCloud
{
  public:
    std::vector<sad::PointType> points;

    // 返回数据集中的点的数量
    inline size_t kdtree_get_point_count() const
    {
        return points.size();
    }

    // 返回给定索引处的点的指定维度的值
    inline double kdtree_get_pt(const size_t idx, const size_t dim) const
    {
        if (dim == 0)
            return points[idx].x;
        else if (dim == 1)
            return points[idx].y;
        else if (dim == 2)
            return points[idx].z;
        return 0.0;
    }

    // 估计数据集的边界框
    template <class BBOX> bool kdtree_get_bbox(BBOX &) const
    {
        return false;
    }
};

typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, NanoPointCloud>, NanoPointCloud, 3>
    NanoFlannKdTree;

void nanoflannSearch(NanoFlannKdTree &kdtree_nanoflann, const int &num_search,
                     std::vector<std::pair<size_t, size_t>> &matches, const NanoPointCloud &second)
{
    std::vector<unsigned int> result_idx(num_search);
    std::vector<double> result_dist(num_search);

    for (int i = 0; i < second.points.size(); i++)
    {
        double query_pt[3] = {second.points[i].x, second.points[i].y, second.points[i].z};
        kdtree_nanoflann.knnSearch(&query_pt[0], num_search, &result_idx[0], &result_dist[0]);
        for (int j = 0; j < result_idx.size(); ++j)
        {
            int m = result_idx[j];
            double d = result_dist[j];
            matches.push_back({m, i});
        }
    }
}