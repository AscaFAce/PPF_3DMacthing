
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/features/normal_3d.h>
#include <pcl/registration/ppf_registration.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <vector>
#include <iostream>
#include <sys/time.h>
#include "ppf3d_public.h"

using namespace pcl;
using namespace std;


int main(int argc, char** argv)
{

    /******************   opencv_contrib 改版  ********************/


    struct timeval tstart, tend, mstart, mend, ticps, ticpe;
    double timer;


    string modelFileName = "/home/ds/Desktop/sample sets/test/mian_T-rex_high.ply";
    string sceneFileName = "/home/ds/Desktop/sample sets/test/rs1.ply";
    string SaveFileName  = "/home/ds/Desktop/sample sets/test/ppf_m2e_match.ply";


    std::vector<Eigen::Vector3f> M_pc, M_nr;
    std::vector<Eigen::Vector3f> S_pc, S_nr;

    ppf3d_public detector(0.05, 0.05);

    // Load a PLY File
    //Mat pcModel = detector.loadPLYSimple(modelFileName.c_str(), 1);
    detector.loadPLYFile(modelFileName.c_str(), M_pc, M_nr);

    // Read the scene
    //Mat pcScene = detector.loadPLYSimple(sceneFileName.c_str(), 1);
    detector.loadPLYFile(sceneFileName.c_str(), S_pc, S_nr);

    // Now train the model
    std::cout << "M2E Training..." << std::endl;

    gettimeofday(&tstart, NULL);
    detector.trainM2E(M_pc, M_nr);
    gettimeofday(&tend, NULL);
    timer= (tend.tv_sec-tstart.tv_sec)*1000 + (tend.tv_usec-tstart.tv_usec)/1000;
    std::cout << std::endl << "M2E Training complete in " << timer << " msec" << std::endl;

    // Match the model to the scene and get the pose
    std::cout << std::endl << "M2E Starting matching..." << std::endl;

    vector<Pose_3DM2E> resultsM2E;
    gettimeofday(&mstart, NULL);
    detector.matchM2E(S_pc, S_nr, resultsM2E, 1.0/10.0, 0.03);
    gettimeofday(&mend, NULL);
    timer= (mend.tv_sec-mstart.tv_sec)*1000 + (mend.tv_usec-mstart.tv_usec)/1000;
    std::cout << std::endl << "M2E Matching complete in " << timer << " msec" << std::endl;


    std::vector<Eigen::Vector3f> ptm_save;
    for (int i = 0; i < M_pc.size(); ++i)
    {
        Eigen::Vector4f v_tmp;
        v_tmp << M_pc[i] , 1;

        Eigen::Vector4f  pt_after_m = resultsM2E[0].pose*v_tmp;
        ptm_save.push_back(pt_after_m.head(3));
    }
    detector.savePLY(SaveFileName, ptm_save);

    /******************   界面显示  ********************/

    PointCloud<PointXYZ>::Ptr cloud_model(new PointCloud<PointXYZ>);
    io::loadPLYFile<PointXYZ>(modelFileName, *cloud_model);

    PointCloud<PointXYZ>::Ptr cloud_scene(new PointCloud<PointXYZ>);
    io::loadPLYFile<PointXYZ>(sceneFileName, *cloud_scene);

    PointCloud<PointXYZ>::Ptr cloud_match(new PointCloud<PointXYZ>);
    io::loadPLYFile<PointXYZ>(SaveFileName, *cloud_match);

    boost::shared_ptr<visualization::PCLVisualizer> viewer(new visualization::PCLVisualizer("3D viewer"));

    int v1(0), v2(0);

    viewer->createViewPort(0.0, 0.0, 0.5, 1.0, v1);
    viewer->setBackgroundColor(0, 0, 0, v1);
    viewer->addPointCloud<pcl::PointXYZ>(cloud_model,  "cloud model",v1);

    viewer->createViewPort(0.5, 0.0, 1.0, 1.0, v2);
    viewer->setBackgroundColor(0, 0, 0, v2);
    viewer->addPointCloud<pcl::PointXYZ>(cloud_scene, "cloud scene",v2);

    visualization::PointCloudColorHandlerCustom<PointXYZ> cloud2(cloud_match, 0, 255, 0);
    viewer->addPointCloud<pcl::PointXYZ>(cloud_match, cloud2, "cloud match", v2);

    while (!viewer->wasStopped())
    {
        viewer->spinOnce(100);
        boost::this_thread::sleep(boost::posix_time::microseconds(100000));
    }

    return 0;
}




