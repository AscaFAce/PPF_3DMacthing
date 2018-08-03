#include "opencv2/surface_matching.hpp"
#include <iostream>
#include "opencv2/surface_matching/ppf_helpers.hpp"
#include "opencv2/core/utility.hpp"
#include "opencv2/surface_matching/ppf_match_3d.hpp"
#include <pcl/features/ppf.h>
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
#include "utils.h"
#include <vector>



using namespace pcl;
using namespace std;
using namespace cv;
using namespace ppf_match_3d;



int main(int argc, char** argv)
{
//    vector<Eigen::Vector3d> pts, sampledPC, sampledPCNor;
//    vector<Eigen::Vector3d> nos;
//    vector<Eigen::Vector3d> pts_scene, nos_scene;

//    loadPLY("/home/ds/Desktop/ppf-master/data/mian_T-rex_high.ply", pts, nos);
//    loadPLY("/home/ds/Desktop/ppf-master/data/rs1.ply", pts_scene, nos_scene);


//    std::map< int,std::vector< std::vector<double> > > hash_table;
//    double angle_step = 0, distance_step = 0;
//    int model_num = 0;
//    trainModel(pts, nos, sampledPC, sampledPCNor, angle_step, distance_step, model_num, hash_table);

//    cout <<"hash table size: " <<  hash_table.size() << endl;

//    vector<Pose3D> poses_list;
//    match(pts_scene,nos_scene,sampledPC,sampledPCNor,hash_table,distance_step,angle_step,model_num,poses_list);
//    cout << poses_list.size() << endl;
//    sortPoses(poses_list);

//    vector<Cluster> clusters;
//    clusterPoses(poses_list,clusters,27.3,angle_step);
//    vector<Pose3D> ans_poses;
//    cout << "clusters: " << clusters.size() << endl;

//    averageClusters(clusters, ans_poses);
//    sortPoses(ans_poses);
//    cout << ans_poses[0].pose.rotation() << endl;
//    cout << ans_poses[0].pose.translation() << endl;
//    vector<Eigen::Vector3d> pts_save;

//    //    Eigen::MatrixXd TT(4,4);
//    //    TT << 0.986843,0.0105495,0.161341,-15.9453,-0.0349336,-0.960388,0.276469,185.037,0.157867,-0.278468,-0.947383,-1199.47,0,0,0,1;
//    //    std::cout << TT << std::endl;
//    for (int i = 0; i < pts.size(); ++i)
//    {
//        Eigen::Vector4d v_tmp;
//        v_tmp << pts[i] , 1;
//        Eigen::Vector4d pt_after = ans_poses[0].pose*v_tmp;
//        pts_save.push_back(pt_after.head(3));
//    }
    //savePLY("PPF_Matching.ply",pts_save);


//    PointCloud<PointXYZ>::Ptr cloud(new PointCloud<PointXYZ>);
//    io::loadPLYFile<PointXYZ>("/home/ds/Desktop/ppf-master/data/save.ply", *cloud);

//    visualization::CloudViewer viewer("Simple Cloud Viewer");//直接创造一个显示窗口
//    viewer.showCloud(cloud);//再这个窗口显示点云
//    while (!viewer.wasStopped())
//    {
//    }

//    return 0;


    /****************************************************/
    /**************　PLY文件法向量计算和添加　　**************/
    /****************************************************/

    /*
    string modelFileName = "/home/ds/Desktop/sample sets/bun_zipper.ply";
    string outputFileName = "/home/ds/Desktop/sample sets/bun_zipper_nor.ply";
    cv::Mat points, pointsAndNormals;

    cout << "Loading points\n";
    cv::ppf_match_3d::loadPLYSimple(modelFileName.c_str(), 1).copyTo(points);

    cout << "Computing normals\n";
    cv::Vec3d viewpoint(0, 0, 0);
    cv::ppf_match_3d::computeNormalsPC3d(points, pointsAndNormals, 2, true, viewpoint);

    std::cout << "Writing points\n";
    cv::ppf_match_3d::writePLY(pointsAndNormals, outputFileName.c_str());

    //the following function can also be used for debugging purposes
    //cv::ppf_match_3d::writePLYVisibleNormals(pointsAndNormals, outputFileName.c_str());


    PointCloud<PointXYZ>::Ptr cloud_model(new PointCloud<PointXYZ>);
    io::loadPLYFile<PointXYZ>("/home/ds/Desktop/sample sets/mian_T-rex_high.ply", *cloud_model);

    PointCloud<PointXYZ>::Ptr cloud_model_nor(new PointCloud<PointXYZ>);
    io::loadPLYFile<PointXYZ>("/home/ds/Desktop/sample sets/bun_zipper_nor.ply", *cloud_model_nor);


    boost::shared_ptr<visualization::PCLVisualizer> viewer(new visualization::PCLVisualizer("3D viewer"));
    int v1(0), v2(0), v3(0), v4(0);
    viewer->createViewPort(0.0, 0.0, 0.5, 1.0,v1);
    viewer->setBackgroundColor(0, 0, 0, v1);
    visualization::PointCloudColorHandlerCustom<PointXYZ> cloud0(cloud_model, 0, 0, 255);
    viewer->addPointCloud<pcl::PointXYZ>(cloud_model, cloud0, "cloud model",v1);


    viewer->createViewPort(0.5, 0.0, 1.0, 1.0,v2);
    viewer->setBackgroundColor(0, 0, 0, v2);
    visualization::PointCloudColorHandlerCustom<PointXYZ> cloud1(cloud_model_nor, 0, 255,0);
    viewer->addPointCloud<pcl::PointXYZ>(cloud_model_nor, cloud1, "cloud model_nor",v2);



    viewer->setPointCloudRenderingProperties(visualization::PCL_VISUALIZER_FONT_SIZE, 3, "cloud model");
    viewer->setPointCloudRenderingProperties(visualization::PCL_VISUALIZER_FONT_SIZE, 3, "cloud model_nor");
    viewer->addCoordinateSystem(1.0);

    //viewer->addPointCloudNormals<PointXYZ, Normal>(cloud_model, normals, 10, 0.05, "normals1", v1);
    while (!viewer->wasStopped())
    {
        viewer->spinOnce(100);
        boost::this_thread::sleep(boost::posix_time::microseconds(100000));
    }
    */


    string modelFileName = "/home/ds/Desktop/sample sets/chicken3.ply";
    string sceneFileName = "/home/ds/Desktop/sample sets/rs1_normals.ply";

    // Load a PLY File
    Mat pc = loadPLYSimple(modelFileName.c_str(), 1);


    // Now train the model
    cout << "Training..." << endl;
    int64 tick1 = cv::getTickCount();
    ppf_match_3d::PPF3DDetector detector(0.03, 0.05);
    detector.trainModel(pc);
    int64 tick2 = cv::getTickCount();
    cout << endl << "Training complete in " << (double)(tick2-tick1)/ cv::getTickFrequency() << " sec" << endl << "Loading model..." << endl;


    // Read the scene
    Mat pcTest = loadPLYSimple(sceneFileName.c_str(), 1);

    // Match the model to the scene and get the pose
    cout << endl << "Starting matching..." << endl;
    vector<Pose3DPtr> results;
    tick1 = cv::getTickCount();
    detector.match(pcTest, results, 1.0/10.0, 0.05);
    tick2 = cv::getTickCount();
    cout << endl << "PPF Elapsed Time " << (tick2-tick1)/cv::getTickFrequency() << " sec" << endl;


    // Get only first N results
    int N = 2;
    vector<Pose3DPtr>::const_iterator first = results.begin();
    vector<Pose3DPtr>::const_iterator last = results.begin() + N;
    vector<Pose3DPtr> resultsSub(first, last);

    // Create an instance of ICP
    ICP icp(200, 0.001f, 2.5f, 8);
    int64 t1 = cv::getTickCount();

    // Register for all selected poses
    cout << endl << "Performing ICP on " << N << " poses..." << endl;

    icp.registerModelToScene(pc, pcTest, resultsSub);
    int64 t2 = cv::getTickCount();
    cout << endl << "ICP Elapsed Time " << (t2-t1)/cv::getTickFrequency() << " sec" << endl;
    //cout << "Poses: " << endl;

    // debug first five poses
    for (size_t i=0; i<resultsSub.size(); i++)
    {
        Pose3DPtr result = resultsSub[i];
        cout << "Pose Result " << i << endl;
        result->printPose();
        if (i==1)
        {
            Mat pct = transformPCPose(pc, result->pose);
            writePLY(pct, "PPF_Matching.ply");
        }
    }


    PointCloud<PointXYZ>::Ptr cloud_model(new PointCloud<PointXYZ>);
    io::loadPLYFile<PointXYZ>(modelFileName, *cloud_model);

    PointCloud<PointXYZ>::Ptr cloud_scene(new PointCloud<PointXYZ>);
    io::loadPLYFile<PointXYZ>(sceneFileName, *cloud_scene);

    PointCloud<PointXYZ>::Ptr cloud_match(new PointCloud<PointXYZ>);
    io::loadPLYFile<PointXYZ>("PPF_Matching.ply", *cloud_match);


    boost::shared_ptr<visualization::PCLVisualizer> viewer(new visualization::PCLVisualizer("3D viewer"));
//    //visualization::CloudViewer viewer("Simple Cloud Viewer");//直接创造一个显示窗口
    int v1(0), v2(0), v3(0), v4(0);

    viewer->createViewPort(0.0, 0.0, 0.5, 0.5, v1);
    viewer->setBackgroundColor(0, 0, 0, v1);
    viewer->addPointCloud<pcl::PointXYZ>(cloud_model,  "cloud model",v1);


    viewer->createViewPort(0.5, 0.0, 1.0, 0.5, v2);
    viewer->setBackgroundColor(0, 0, 0, v2);
    viewer->addPointCloud<pcl::PointXYZ>(cloud_scene, "cloud scene",v2);

    //viewer->createViewPort(0.5, 0.5, 1.0, 1.0,v3);
    visualization::PointCloudColorHandlerCustom<PointXYZ> cloud2(cloud_match, 0, 255, 0);
    viewer->addPointCloud<pcl::PointXYZ>(cloud_match, cloud2, "cloud match", v2);
    //viewer->addCoordinateSystem(1.0);

    while (!viewer->wasStopped())
    {
        viewer->spinOnce(100);
        boost::this_thread::sleep(boost::posix_time::microseconds(100000));
    }

    return 0;

}
