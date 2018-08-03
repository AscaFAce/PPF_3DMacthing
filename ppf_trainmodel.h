#ifndef PPF_TRAINMODEL_H
#define PPF_TRAINMODEL_H

#include "ppf_public.h"


class PPF_TrainModel
{
private:
    ppf_public m_ppf_public;

public:
    PPF_TrainModel();

    void trainModel(std::vector<Eigen::Vector3d>& pts, std::vector<Eigen::Vector3d>& nor,
                    std::vector<Eigen::Vector3d>& sampledPC, std::vector<Eigen::Vector3d>& sampledPCNor,
                    double& angle_step, double& distance_step,int& model_num,
                    std::map<int,std::vector<std::vector<double> > >& hash_table);


};

#endif // PPF_TRAINMODEL_H
