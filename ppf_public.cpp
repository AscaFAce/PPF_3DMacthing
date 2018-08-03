#include "ppf_public.h"

ppf_public::ppf_public()
{

}

bool /*ppf_public::*/compare_x(Eigen::Vector3d& i, Eigen::Vector3d& j)
{
    return i.x() < j.x();
}

bool /*ppf_public::*/compare_y(Eigen::Vector3d& i, Eigen::Vector3d& j)
{
    return i.y() < j.y();
}

bool /*ppf_public::*/compare_z(Eigen::Vector3d& i, Eigen::Vector3d& j)
{
    return i.z() < j.z();
}


double ppf_public::rangeMax(std::vector<Eigen::Vector3d> &pts, int dim)
{

    if(dim == 1)
    {
        return (*std::max_element(pts.begin(),pts.end(),compare_x)).x();
    }else if(dim == 2)
    {
        return (*std::max_element(pts.begin(),pts.end(),compare_y)).y();
    }else if(dim == 3)
    {
        return (*std::max_element(pts.begin(),pts.end(),compare_z)).z();
    }else
    {
        std::cerr << "The dim of point cloud is error" << std::endl;
    }
}

double ppf_public::rangeMin(std::vector<Eigen::Vector3d> &pts, int dim)
{

    if(dim == 1)
    {
        return (*std::min_element(pts.begin(),pts.end(),compare_x)).x();
    }else if(dim == 2)
    {
        return (*std::min_element(pts.begin(),pts.end(),compare_y)).y();
    }else if(dim == 3)
    {
        return (*std::min_element(pts.begin(),pts.end(),compare_z)).z();
    }else
    {
        std::cerr << "The dim of point cloud is error" << std::endl;
    }
}


bool ppf_public::checkDistance(std::vector<Eigen::Vector3d> &pts,std::vector<int>& neigh, int k, double rs)
{
    for (int i = 0; i < neigh.size(); ++i) {
        if (std::pow((pts[neigh[i]].x() - pts[k].x()),2)+std::pow((pts[neigh[i]].y() - pts[k].y()),2)+std::pow((pts[neigh[i]].z() - pts[k].z()),2) < rs)
        {
            return false;
        }
    }
    return true;
}


void ppf_public::samplePCPoisson(std::vector<Eigen::Vector3d> &pts,std::vector<Eigen::Vector3d> &nos, double sample_step,
                     std::vector<Eigen::Vector3d> &samplePC,std::vector<Eigen::Vector3d> &samplePCNor)
{
    // poisson disc sampling
    std::vector<double> range_x, range_y, range_z;
    // 求包围盒
    range_x.push_back(rangeMin(pts,1));
    range_x.push_back(rangeMax(pts,1));

    range_y.push_back(rangeMin(pts,2));
    range_y.push_back(rangeMax(pts,2));

    range_z.push_back(rangeMin(pts,3));
    range_z.push_back(rangeMax(pts,3));

    double dx = range_x.back() - range_x.front();
    double dy = range_y.back() - range_y.front();
    double dz = range_z.back() - range_z.front();

    double d = sqrt(dx*dx + dy*dy + dz*dz); //包围盒 对角线长度
    double r = d*sample_step; //点点之间最小距离
    double rs = r*r;

    double boxsize = r/sqrt(3); //分割网格的边长

    int samples_in_dimx = static_cast<int>(floor(dx / boxsize));
    int samples_in_dimy = static_cast<int>(floor(dy / boxsize));
    int samples_in_dimz = static_cast<int>(floor(dz / boxsize));

    int map[samples_in_dimx][samples_in_dimy][samples_in_dimz]  = {0};

    for (int i = 0; i < pts.size(); ++i)
    {
        int x_cell = floor(samples_in_dimx*(pts[i].x()-range_x[0])/dx);
        int y_cell = floor(samples_in_dimy*(pts[i].y()-range_y[0])/dy);
        int z_cell = floor(samples_in_dimz*(pts[i].z()-range_z[0])/dz);

        std::vector<int> neigh;
        // 选取周围的5*5*5个box
        for (int j = 0; j < 5; ++j)
        {
            for (int k = 0; k < 5; ++k)
            {
                for (int l = 0; l < 5; ++l)
                {
                    int pos_x = std::max(x_cell-2+j,0);
                    pos_x = std::min(pos_x,samples_in_dimx-1);

                    int pos_y = std::max(y_cell-2+k,0);
                    pos_y = std::min(pos_y,samples_in_dimy-1);

                    int pos_z = std::max(z_cell-2+l,0);
                    pos_z = std::min(pos_z,samples_in_dimz-1);

                    if(map[pos_x][pos_y][pos_z]!=0)
                    {
                        neigh.push_back(map[pos_x][pos_y][pos_z]);
                    }
                }
            }
        }

        if (neigh.size() == 0)
        {
            map[x_cell][y_cell][z_cell] = i;
        }
        else if(checkDistance(pts,neigh,i,rs))
        {
            //查看距离
            map[x_cell][y_cell][z_cell] = i;
        }

    }

    std::vector<int> index;
    for (int m = 0; m < samples_in_dimx; ++m)
    {
        for (int i = 0; i < samples_in_dimy; ++i)
        {
            for (int j = 0; j < samples_in_dimz; ++j)
            {
                if (map[m][i][j]!=0)
                {
                    index.push_back(map[m][i][j]);
                }
            }
        }
    }
    for (int n = 0; n < index.size(); ++n)
    {
        samplePC.push_back(pts[index[n]]);
        samplePCNor.push_back(nos[index[n]].normalized());
    }
}


void ppf_public::computePPF(Eigen::Vector3d& pt1, Eigen::Vector3d& nor1, Eigen::Vector3d& pt2, Eigen::Vector3d& nor2,std::vector<double>& f)
{
    Eigen::Vector3d d = pt2 - pt1;
    double norm = d.norm();
//    d.normalize();

    f.push_back(atan2(nor1.cross(d).norm(), nor1.dot(d)));
    f.push_back(atan2(nor2.cross(d).norm(), nor2.dot(d)));
    f.push_back(atan2(nor1.cross(nor2).norm(), nor1.dot(nor2)));
    f.push_back(norm);
    /**
    double a[6] = {0}; double b[6] = {0};
    for (int i = 0; i < 3; ++i) {
        a[i] = pt1(i);
        b[i] = pt2(i);
        a[i+3] = nor1(i);
        b[i+3] = nor2(i);
    }

    double d[3]= {b[0]-a[0], b[1]-a[1], b[2]-a[2]};
    double norm=sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    d[0] /= norm;
    d[1] /= norm;
    d[2] /= norm;
    double cross1[3]= { a[4]*d[2]-a[5]*d[1]  ,  a[5]*d[0]-a[3]*d[2]  ,  a[3]*d[1]-a[4]*d[0] };
    double cross2[3]= { b[4]*d[2]-b[5]*d[1]  ,  b[5]*d[0]-b[3]*d[2]  ,  b[3]*d[1]-b[4]*d[0] };
    double cross3[3]= { a[4]*b[5]-a[5]*b[4]  ,  a[5]*b[3]-a[3]*b[5]  ,  a[3]*b[4]-a[4]*b[3] };
    f.push_back(atan2( sqrt(cross1[0]*cross1[0] + cross1[1]*cross1[1] + cross1[2]*cross1[2]) ,  d[0]*a[3]+d[1]*a[4]+d[2]*a[5] ));
    f.push_back(atan2( sqrt(cross2[0]*cross2[0] + cross2[1]*cross2[1] + cross2[2]*cross2[2]) ,  d[0]*b[3]+d[1]*b[4]+d[2]*b[5] ));
    f.push_back(atan2( sqrt(cross3[0]*cross3[0] + cross3[1]*cross3[1] + cross3[2]*cross3[2]) ,  a[3]*b[3]+a[4]*b[4]+a[5]*b[5] ));
    f.push_back(norm);
     **/

}


uint32_t ppf_public::hashPPF(std::vector<double> f, double angle_step, double distance_step)
{
    uint32_t *key = new uint32_t [4];             /* input scalar */
//    char key[4] = {0};
    key[0] = static_cast<uint32_t>(floor(f[0]/angle_step));
    key[1] = static_cast<uint32_t>(floor(f[1]/angle_step));
    key[2] = static_cast<uint32_t>(floor(f[2]/angle_step));
    key[3] = static_cast<uint32_t>(floor(f[3]/distance_step));


    uint32_t* outMatrix = new uint32_t [4];              /* output matrix */

    MurmurHash3_x86_32(key,16,42,outMatrix);

    return outMatrix[0];
//    return murmur3_32(key,16,42);
}

double ppf_public::computeAlpha(Eigen::Vector3d pt1,Eigen::Vector3d nor1,Eigen::Vector3d pt2,Eigen::Vector3d nor2)
{
    Eigen::Affine3d T = transformRT(pt1,nor1);
    Eigen::Vector3d mpt = T*pt2;
    double alpha = atan2(-mpt.z(),mpt.y());
    return alpha;
}

Eigen::Affine3d ppf_public::transformRT(Eigen::Vector3d pt,Eigen::Vector3d nor)
{
    /**
     * 把一个点法，归一到法向量在x轴上
     */

    double angle = acos(nor.x());
    Eigen::Vector3d axis(0,nor.z(),-nor.y());
    if(nor.y() == 0 && nor.z() == 0 )
    {
        axis.x() = 0;
        axis.y() = 1;
        axis.z() = 0;
    } else
    {
        axis.normalize();
    }
    Eigen::AngleAxisd angle_axis(angle,axis);
    Eigen::Matrix3d R = angle_axis.toRotationMatrix();
    Eigen::Vector3d t = -R*pt;
    Eigen::Affine3d T;
    T = R;
    T.translation() = t;
    return T;
}
