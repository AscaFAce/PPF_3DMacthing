#ifndef PPF_PUBLIC_H
#define PPF_PUBLIC_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <algorithm>
#include <map>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "MurmurHash3.h"


struct Pose3D
{
    Eigen::Affine3d pose;
    int votes;
};

struct Cluster
{
    std::vector<Pose3D> poses;
    int accu_votes;
};


static bool compare_x(Eigen::Vector3d& i, Eigen::Vector3d& j);
static bool compare_y(Eigen::Vector3d& i, Eigen::Vector3d& j);
static bool compare_z(Eigen::Vector3d& i, Eigen::Vector3d& j);

class ppf_public
{
public:
    ppf_public();

public:
    void samplePCPoisson(std::vector<Eigen::Vector3d> &pts, std::vector<Eigen::Vector3d> &nos,
                         double sample_step, std::vector<Eigen::Vector3d> &samplePC,
                         std::vector<Eigen::Vector3d> &samplePCNor);

    void computePPF(Eigen::Vector3d& pt1,  Eigen::Vector3d& nor1, Eigen::Vector3d& pt2,
                    Eigen::Vector3d& nor2, std::vector<double>& f);

    uint32_t hashPPF(std::vector<double> f, double angle_step, double distance_step);

    double computeAlpha(Eigen::Vector3d pt1, Eigen::Vector3d nor1,
                        Eigen::Vector3d pt2, Eigen::Vector3d nor2);

    double rangeMin(std::vector<Eigen::Vector3d> &pts, int dim);

    double rangeMax(std::vector<Eigen::Vector3d> &pts, int dim);

    bool checkDistance(std::vector<Eigen::Vector3d> &pts,std::vector<int>& neigh, int k, double rs);

    Eigen::Affine3d transformRT(Eigen::Vector3d pt,Eigen::Vector3d nor);

};

#endif // PPF_PUBLIC_H
