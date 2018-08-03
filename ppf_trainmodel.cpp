#include "ppf_trainmodel.h"



PPF_TrainModel::PPF_TrainModel()
{

}

void PPF_TrainModel::trainModel(std::vector<Eigen::Vector3d>& pts, std::vector<Eigen::Vector3d>& nor, std::vector<Eigen::Vector3d>& sampledPC, std::vector<Eigen::Vector3d>& sampledPCNor,double& angle_step, double& distance_step,int& model_num,std::map<int,std::vector<std::vector<double> > >& hash_table)
{
    //    double angle_step, distance_step;
    std::vector<double> range_x, range_y, range_z;
    // 求包围盒
    range_x.push_back(m_ppf_public.rangeMin(pts,1));
    range_x.push_back(m_ppf_public.rangeMax(pts,1));

    range_y.push_back(m_ppf_public.rangeMin(pts,2));
    range_y.push_back(m_ppf_public.rangeMax(pts,2));

    range_z.push_back(m_ppf_public.rangeMin(pts,3));
    range_z.push_back(m_ppf_public.rangeMax(pts,3));

    double dx = range_x.back() - range_x.front();
    double dy = range_y.back() - range_y.front();
    double dz = range_z.back() - range_z.front();

    double d = sqrt(dx*dx + dy*dy + dz*dz); //包围盒 对角线长度

    distance_step = d*0.04;
    angle_step    = (360/30)*M_PI/180;



    m_ppf_public.samplePCPoisson(pts,nor,0.04,sampledPC,sampledPCNor);


    std::cout << "sampled Point cloud size: " << sampledPC.size() << std::endl;
    model_num = sampledPC.size();
    //    std::map<int,std::vector<std::vector<double>>> hash_table;

    for (int i = 0; i < sampledPC.size(); ++i) {
        Eigen::Vector3d pt1 = sampledPC[i];
        Eigen::Vector3d nor1 = sampledPCNor[i];
        for (int j = 0; j < sampledPC.size(); ++j) {
            if(i!=j)
            {
                Eigen::Vector3d pt2 = sampledPC[j];
                Eigen::Vector3d nor2 = sampledPCNor[j];
                std::vector<double> f ;
                m_ppf_public.computePPF(pt1,nor1,pt2,nor2,f);
                uint32_t hash = m_ppf_public.hashPPF(f,angle_step,distance_step);
                double model_alpha = m_ppf_public.computeAlpha(pt1,nor1,pt2,nor2);
                std::vector<double> node;

                node.push_back(i);
                node.push_back(model_alpha);

                std::map<int,std::vector<std::vector<double> > >::iterator iter;
                iter = hash_table.find(hash);
                if (iter == hash_table.end())
                {
                    std::vector<std::vector<double> > nodes;
                    nodes.push_back(node);
                    hash_table[hash] = nodes;
                }else{
                    std::vector<std::vector<double> > nodes;
                    //                    std::cout << "got one " << std::endl;
                    nodes = hash_table[hash];
                    nodes.push_back(node);
                    hash_table[hash] = nodes;
                }
            }
        }
        if(i%10 == 0)
        {
            //std::cout << "trained: " << std::round(100*i/sampledPC.size()) << "%" << std::endl;
        }
    }
}

