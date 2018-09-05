#ifndef PPF3D_PUBLIC_H
#define PPF3D_PUBLIC_H

#include "t_hash_int.hpp"
#include "stdafx.h"

class ppf3d_public
{
public:
    ppf3d_public();

    ppf3d_public(const double relativeSamplingStep, const double relativeDistanceStep=0.05, const double numAngles=30);

    //virtual ~ppf3d_public();

    void    setSearchParams(const double positionThreshold=-1, const double rotationThreshold=-1, const bool useWeightedClustering=false);

    void    trainM2E(std::vector<Eigen::Vector3f> M_pc, std::vector<Eigen::Vector3f> M_nr);

    void    matchM2E(std::vector<Eigen::Vector3f> S_pc, std::vector<Eigen::Vector3f> S_nr, std::vector<Pose_3DM2E>& results, const double relativeSceneSampleStep, const double relativeSceneDistance);

    void    loadPLYFile(const char* fileName, std::vector<Eigen::Vector3f> &PC, std::vector<Eigen::Vector3f> &NOR);

    void    savePLY(const std::string filename, std::vector<Eigen::Vector3f> &pts);




protected:

    double angle_stepM2E, angle_step_radiansM2E, distance_stepM2E;
    double sampling_step_relativeM2E, angle_step_relativeM2E, distance_step_relativeM2E;
    int num_ref_pointsM2E;

    hashtable_int* hash_tableM2E;
    THash* hash_nodesM2E;

    double position_threshold, rotation_threshold;
    bool use_weighted_avg, trained;

    int scene_sample_stepM2E;

    std::vector<Eigen::Vector3f> sample_pc, sample_nr;

    void    clearTrainingModels();

    void    samplePCM2E(std::vector<Eigen::Vector3f> pc, std::vector<Eigen::Vector3f> nr, std::vector<Eigen::Vector3f> &samoled_pc, std::vector<Eigen::Vector3f> &samoled_nr, float sampleStep, int weightByCenter);

    double  computeAlphaM2E(const Eigen::Vector3f& pt1, const Eigen::Vector3f& nor1, const Eigen::Vector3f& pt2);

    void    computeBBoxM2E(std::vector<Eigen::Vector3f>& pts_Model, Eigen::Vector2f& range_x, Eigen::Vector2f& range_y, Eigen::Vector2f& range_z);

    void    computePPFFeaturesM2E(const Eigen::Vector3f& pt1, const Eigen::Vector3f& nor1, const Eigen::Vector3f& pt2, const Eigen::Vector3f& nor2, Eigen::Vector4f& f);

    void    computetransformRTM2E(const Eigen::Vector3f& p1, const Eigen::Vector3f& n1,  Eigen::Matrix3f& R, Eigen::Vector3f& t);

    void    TNormalize3M2E(Eigen::Vector3f& v);

    void    aatoRM2E(const Eigen::Vector3f& axis, double angle, Eigen::Matrix3f& R);

    bool    matchPoseM2E(const Pose_3DM2E& sourcePose, const Pose_3DM2E& targetPose);

    void    clusterPosesM2E(std::vector<Pose_3DM2E>& poseList, int numPoses, std::vector<Pose_3DM2E> &finalPoses);



};

#endif // PPF3D_PUBLIC_H
