#include "ppf3d_public.h"
#include "hash_murmur.hpp"

#define EPS     1e-8

static const size_t PPF_LENGTH = 5;

static std::vector<std::string> split(const std::string &text, char sep)
{
    std::vector<std::string> tokens;
    std::size_t start = 0, end = 0;
    while ((end = text.find(sep, start)) != std::string::npos)
    {
        tokens.push_back(text.substr(start, end - start));
        start = end + 1;
    }
    tokens.push_back(text.substr(start));
    return tokens;
}

static bool pose_3DCompareM2E(const Pose_3DM2E& a, const Pose_3DM2E& b)
{
    return ( a.numVotes > b.numVotes);
}

static int Pose_3DClustersM2E(const ClusterM2E& a, const ClusterM2E& b)
{
  return ( a.accu_votes > b.accu_votes);
}

static KeyType hashPPFM2E(const Eigen::Vector4f& f, const double AngleStep, const double DistanceStep)
{
    KeyType *key = new KeyType[4];

    key[0] = (int)(f[0] / AngleStep);
    key[1] = (int)(f[1] / AngleStep);
    key[2] = (int)(f[2] / AngleStep);
    key[3] = (int)(f[3] / DistanceStep);

    KeyType hashKey = 0;

    murmurHash(key, 4*sizeof(int), 42, &hashKey);
    return hashKey;
}

static void rtToPoseM2E(const Eigen::Matrix3f& R, const Eigen::Vector3f& t, Eigen::Matrix4f& Pose)
{
    //Matx34d P;
    //hconcat(R, t, P);
    //vconcat(P, Matx14d(0, 0, 0, 1), Pose);
    Eigen::Matrix4f P;
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
        {
           P(i,j) = R(i,j);
        }
    }

    for(int i=0; i<3; i++)
    {
       P(i,3) = t(i,0);
    }

    P(3,0) = 0;
    P(3,1) = 0;
    P(3,2) = 0;
    P(3,3) = 1;

    Pose = P;
}

static bool compare_x(Eigen::Vector3f& i, Eigen::Vector3f& j)
{
    return i.x() < j.x();
}

static bool compare_y(Eigen::Vector3f& i, Eigen::Vector3f& j)
{
    return i.y() < j.y();
}

static bool compare_z(Eigen::Vector3f& i, Eigen::Vector3f& j)
{
    return i.z() < j.z();
}

static double rangeMax(std::vector<Eigen::Vector3f> &pts, int dim)
{

    if(dim == 1)
    {
        return (*max_element(pts.begin(),pts.end(),compare_x)).x();
    }
    else if(dim == 2)
    {
        return (*max_element(pts.begin(),pts.end(),compare_y)).y();
    }
    else if(dim == 3)
    {
        return (*max_element(pts.begin(),pts.end(),compare_z)).z();
    }
    else
    {
       std::cerr << "The dim of point cloud is error" << std::endl;
    }
}

static double rangeMin(std::vector<Eigen::Vector3f> &pts, int dim)
{

    if(dim == 1)
    {
        return (*min_element(pts.begin(),pts.end(),compare_x)).x();
    }
    else if(dim == 2)
    {
        return (*min_element(pts.begin(),pts.end(),compare_y)).y();
    }
    else if(dim == 3)
    {
        return (*min_element(pts.begin(),pts.end(),compare_z)).z();
    }
    else
    {
        std::cerr << "The dim of point cloud is error" << std::endl;
    }
}

static float medianF(float arr[], int n)
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ;
    high = n-1 ;
    median = (low + high) >>1;
    for (;;)
    {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1)
        {
            /* Two elements only */
            if (arr[low] > arr[high])
                std::swap(arr[low], arr[high]) ;
            return arr[median] ;
        }

        /* Find median of low, middle and high items; swap into position low */
        middle = (low + high) >>1;
        if (arr[middle] > arr[high])
            std::swap(arr[middle], arr[high]) ;
        if (arr[low] > arr[high])
            std::swap(arr[low], arr[high]) ;
        if (arr[middle] > arr[low])
            std::swap(arr[middle], arr[low]) ;

        /* Swap low item (now in position middle) into position (low+1) */
        std::swap(arr[middle], arr[low+1]) ;

        /* Nibble from each end towards middle, swapping items when stuck */
        ll = low + 1;
        hh = high;
        for (;;)
        {
            do
                ll++;
            while (arr[low] > arr[ll]) ;
            do
                hh--;
            while (arr[hh]  > arr[low]) ;

            if (hh < ll)
                break;

            std::swap(arr[ll], arr[hh]) ;
        }

        /* Swap middle item (in position low) back into correct position */
        std::swap(arr[low], arr[hh]) ;

        /* Re-set active partition */
        if (hh <= median)
            low = ll;
        if (hh >= median)
            high = hh - 1;
    }
}

static float getRejectionThreshold(float* r, int m, float outlierScale)
{
    float* t=(float*)calloc(m, sizeof(float));
    int i=0;
    float s=0, medR, threshold;

    memcpy(t, r, m*sizeof(float));
    medR = medianF(t, m);

    for (i=0; i<m; i++)
        t[i] = (float)fabs((double)r[i]-(double)medR);

    s = 1.48257968f * medianF(t, m);

    threshold = (outlierScale*s+medR);

    free(t);
    return threshold;
}

static void getUnitXRotationM2E(double angle, Eigen::Matrix3f& Rx)
{
  const double sx = sin(angle);
  const double cx = cos(angle);

  //Mat(Rx.eye()).copyTo(Rx);
  Rx(0, 0) = 1, Rx(1, 0) = 0, Rx(2, 0) = 0;
  Rx(0, 1) = 0, Rx(1, 1) = 1, Rx(2, 1) = 0;
  Rx(0, 2) = 0, Rx(1, 2) = 0, Rx(2, 2) = 1;

  Rx(1, 1) = cx;
  Rx(1, 2) = -sx;
  Rx(2, 1) = sx;
  Rx(2, 2) = cx;
}

static hashtable_int* getHashtable(int* data, size_t length, int numMaxElement)
{
    hashtable_int* hashtable = hashtableCreate(static_cast<size_t>(numMaxElement*2), 0);
    for (size_t i = 0; i < length; i++)
    {
        const KeyType key = (KeyType)data[i];
        hashtableInsertHashed(hashtable, key+1, reinterpret_cast<void*>(i+1));
    }

    return hashtable;
}



ppf3d_public::ppf3d_public()
{
    sampling_step_relativeM2E = 0.05;
    distance_step_relativeM2E = 0.05;
    scene_sample_stepM2E = (int)(1/0.04);
    angle_step_relativeM2E = 30;
    angle_step_radiansM2E = (360.0/angle_step_relativeM2E)*M_PI/180.0;
    angle_stepM2E = angle_step_radiansM2E;
    trained = false;

    hash_tableM2E = NULL;
    hash_nodesM2E = NULL;

    trained = false;
    setSearchParams();
}

ppf3d_public::ppf3d_public(const double relativeSamplingStep, const double relativeDistanceStep, const double numAngles)
{
    sampling_step_relativeM2E = relativeSamplingStep;
    distance_step_relativeM2E = relativeDistanceStep;
    angle_step_relativeM2E = numAngles;
    angle_step_radiansM2E = (360.0/angle_step_relativeM2E)*M_PI/180.0;
    angle_stepM2E = angle_step_radiansM2E;
    trained = false;

    hash_tableM2E = NULL;
    hash_nodesM2E = NULL;

    trained = false;

    setSearchParams();
}

//ppf3d_public::~ppf3d_public()
//{
//    clearTrainingModels();
//}

void    ppf3d_public::setSearchParams(const double positionThreshold, const double rotationThreshold, const bool useWeightedClustering)
{
    if (positionThreshold<0)
        position_threshold = sampling_step_relativeM2E;
    else
        position_threshold = positionThreshold;

    if (rotationThreshold<0)
        rotation_threshold = ((360/angle_stepM2E) / 180.0 * M_PI);
    else
        rotation_threshold = rotationThreshold;

    use_weighted_avg = useWeightedClustering; //false
}

void    ppf3d_public::TNormalize3M2E(Eigen::Vector3f& v)
{
    double norm = v.norm();
    if (norm > EPS)
    {
        v *= 1.0 / norm;
    }
}

void    ppf3d_public::aatoRM2E(const Eigen::Vector3f& axis, double angle, Eigen::Matrix3f& R)
{
    const double sinA = sin(angle);
    const double cosA = cos(angle);
    const double cos1A = (1 - cosA);
    uint i, j;

    R(0,0) = cosA;  R(0,1) = 0;     R(0,2) = 0;
    R(1,0) = 0;     R(1,1) = cosA;  R(1,2) = 0;
    R(2,0) = 0;     R(2,1) = 0;     R(2,2) = cosA;

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
        {
            if (i != j)
            {
                // Symmetry skew matrix
                R(i, j) += (((i + 1) % 3 == j) ? -1 : 1) * sinA * axis[3 - i - j];
            }
            R(i, j) += cos1A * axis[i] * axis[j];
        }
}

void    ppf3d_public::computetransformRTM2E(const Eigen::Vector3f& p1, const Eigen::Vector3f& n1,  Eigen::Matrix3f& R, Eigen::Vector3f& t)
{
    // dot product with x axis
    double angle = acos(n1[0]);

    // cross product with x axis
    Eigen::Vector3f axis(0, n1[2], -n1[1]);

    // we try to project on the ground plane but it's already parallel
    if (n1[1] == 0 && n1[2] == 0)
    {
        axis[1] = 1;
        axis[2] = 0;
    }
    else
    {
        TNormalize3M2E(axis);
    }

    aatoRM2E(axis, angle, R);
    t = -R * p1;
}

double  ppf3d_public::computeAlphaM2E(const Eigen::Vector3f& pt1, const Eigen::Vector3f& nor1, const Eigen::Vector3f& pt2)
{
    Eigen::Vector3f Tmg, mpt;
    Eigen::Matrix3f R;
    computetransformRTM2E(pt1,nor1, R, Tmg);
    mpt = Tmg + R*pt2;
    double alpha = atan2(-mpt[2],mpt[1]);

    if ( alpha != alpha)
    {
        return 0;
    }

    if (sin(alpha)*mpt[2]<0.0)
        alpha=-alpha;

    return (-alpha);
}

void    ppf3d_public::computeBBoxM2E(std::vector<Eigen::Vector3f>& pts_Model, Eigen::Vector2f& range_x, Eigen::Vector2f& range_y, Eigen::Vector2f& range_z)
{
    range_x[0] = rangeMin(pts_Model,1);
    range_x[1] = rangeMax(pts_Model,1);

    range_y[0] = rangeMin(pts_Model,2);
    range_y[1] = rangeMax(pts_Model,2);

    range_z[0] = rangeMin(pts_Model,3);
    range_z[1] = rangeMax(pts_Model,3);
}

void    ppf3d_public::samplePCM2E(std::vector<Eigen::Vector3f> pc, std::vector<Eigen::Vector3f> nr, std::vector<Eigen::Vector3f> &samoled_pc, std::vector<Eigen::Vector3f> &samoled_nr, float sampleStep, int weightByCenter)
{
    std::vector< std::vector<int> > map;

    int numSamplesDim = (int)(1.0/sampleStep);

    Eigen::Vector2f xrange, yrange, zrange;
    computeBBoxM2E(pc, xrange, yrange, zrange);

    float xr = xrange[1] - xrange[0];
    float yr = yrange[1] - yrange[0];
    float zr = zrange[1] - zrange[0];

    int numPoints = 0;

    map.resize((numSamplesDim+1)*(numSamplesDim+1)*(numSamplesDim+1));

    // OpenMP might seem like a good idea, but it didn't speed this up for me
    //#pragma omp parallel for
    for (int i=0; i<pc.size(); i++)
    {
        //const float* point = pc.ptr<float>(i);
        const Eigen::Vector3f pt = pc[i];
        //const Eigen::Vector3f no = nr[i];
        // quantize a point
        const int xCell =(int) ((float)numSamplesDim*(pt[0]-xrange[0])/xr);
        const int yCell =(int) ((float)numSamplesDim*(pt[1]-yrange[0])/yr);
        const int zCell =(int) ((float)numSamplesDim*(pt[2]-zrange[0])/zr);
        const int index = xCell*numSamplesDim*numSamplesDim+yCell*numSamplesDim+zCell;

        /*#pragma omp critical
        {*/
        map[index].push_back(i);
        //  }
    }

    for (unsigned int i=0; i<map.size(); i++)
    {
        numPoints += (map[i].size()>0);
    }

    //Mat pcSampled = Mat(numPoints, pc.cols, CV_32F);
    //int c = 0;

    for (unsigned int i=0; i<map.size(); i++)
    {
        double px=0, py=0, pz=0;
        double nx=0, ny=0, nz=0;

        std::vector<int> curCell = map[i];
        int cn = (int)curCell.size();
        if (cn>0)
        {
            if (weightByCenter)
            {
                int xCell, yCell, zCell;
                double xc, yc, zc;
                double weightSum = 0 ;
                zCell = i % numSamplesDim;
                yCell = ((i-zCell)/numSamplesDim) % numSamplesDim;
                xCell = ((i-zCell-yCell*numSamplesDim)/(numSamplesDim*numSamplesDim));

                xc = ((double)xCell+0.5) * (double)xr/numSamplesDim + (double)xrange[0];
                yc = ((double)yCell+0.5) * (double)yr/numSamplesDim + (double)yrange[0];
                zc = ((double)zCell+0.5) * (double)zr/numSamplesDim + (double)zrange[0];

                for (int j=0; j<cn; j++)
                {
                    const int ptInd = curCell[j];
                    //float* point = pc.ptr<float>(ptInd);
                    Eigen::Vector3f pt = pc[ptInd];
                    Eigen::Vector3f no = nr[ptInd];
                    const double dx = pt[0]-xc;
                    const double dy = pt[1]-yc;
                    const double dz = pt[2]-zc;
                    const double d = sqrt(dx*dx+dy*dy+dz*dz);
                    double w = 0;

                    if (d>EPS)
                    {
                        // it is possible to use different weighting schemes.
                        // inverse weigthing was just good for me
                        // exp( - (distance/h)**2 )
                        //const double w = exp(-d*d);
                        w = 1.0/d;
                    }

                    //float weights[3]={1,1,1};
                    px += w*(double)pt[0];
                    py += w*(double)pt[1];
                    pz += w*(double)pt[2];
                    nx += w*(double)no[0];
                    ny += w*(double)no[1];
                    nz += w*(double)no[2];

                    weightSum+=w;
                }
                px/=(double)weightSum;
                py/=(double)weightSum;
                pz/=(double)weightSum;
                nx/=(double)weightSum;
                ny/=(double)weightSum;
                nz/=(double)weightSum;
            }
            else
            {
                for (int j=0; j<cn; j++)
                {
                    const int ptInd = curCell[j];
                    //float* point = pc.ptr<float>(ptInd);
                    Eigen::Vector3f pt = pc[ptInd];
                    Eigen::Vector3f no = nr[ptInd];
                    px += (double)pt[0];
                    py += (double)pt[1];
                    pz += (double)pt[2];
                    nx += (double)no[0];
                    ny += (double)no[1];
                    nz += (double)no[2];
                }

                px/=(double)cn;
                py/=(double)cn;
                pz/=(double)cn;
                nx/=(double)cn;
                ny/=(double)cn;
                nz/=(double)cn;

            }

            //float *pcData = pcSampled.ptr<float>(c);
            Eigen::Vector3f samplept;
            Eigen::Vector3f sampleno;

            samplept[0]=(float)px;
            samplept[1]=(float)py;
            samplept[2]=(float)pz;

            // normalize the normals
            double norm = sqrt(nx*nx+ny*ny+nz*nz);

            if (norm>EPS)
            {
                sampleno[0]=(float)(nx/norm);
                sampleno[1]=(float)(ny/norm);
                sampleno[2]=(float)(nz/norm);
            }
            //#pragma omp atomic
            samoled_pc.push_back(samplept);
            samoled_nr.push_back(sampleno);
            //c++;

            curCell.clear();
        }
    }
    map.clear();
}

void    ppf3d_public::computePPFFeaturesM2E(const Eigen::Vector3f& pt1, const Eigen::Vector3f& nor1, const Eigen::Vector3f& pt2, const Eigen::Vector3f& nor2, Eigen::Vector4f& f)
{
    Eigen::Vector3f d = pt2 - pt1;
    double norm = d.norm();
    f[3] = norm;
    if (f[3] <= EPS)
        return;
    d *= 1.0 / f[3];

    f[0] = acos(nor1.dot(d));
    f[1] = acos(nor2.dot(d));
    f[2] = acos(nor1.dot(nor2));
}

void    ppf3d_public::trainM2E(std::vector<Eigen::Vector3f> M_pc, std::vector<Eigen::Vector3f> M_nr)
{
    // compute bbox
    Eigen::Vector2f xRange, yRange, zRange;
    computeBBoxM2E(M_pc, xRange, yRange, zRange);

    // compute sampling step from diameter of bbox
    float dx = xRange[1] - xRange[0];
    float dy = yRange[1] - yRange[0];
    float dz = zRange[1] - zRange[0];
    float diameter = sqrt ( dx * dx + dy * dy + dz * dz );

    float distanceStep = (float)(diameter * sampling_step_relativeM2E);

    //std::vector<Eigen::Vector3f> sample_pc, sample_nr;
    samplePCM2E(M_pc, M_nr, sample_pc, sample_nr, (float)sampling_step_relativeM2E,0);

    int size = sample_pc.size()*sample_pc.size();

    hashtable_int* hashTable = hashtableCreate(size, NULL);

    //int numPPF = sample_pc.size()*sample_pc.size();
    //ppfM2E = Mat(numPPF, PPF_LENGTH, CV_32FC1);

    // TODO: Maybe I could sample 1/5th of them here. Check the performance later.
    int numRefPoints = sample_pc.size();

    // pre-allocate the hash nodes
    hash_nodesM2E = (THash*)calloc(numRefPoints*numRefPoints, sizeof(THash));

    // TODO : This can easily be parallelized. But we have to lock hashtable_insert.
    // I realized that performance drops when this loop is parallelized (unordered
    // inserts into the hashtable
    // But it is still there to be investigated. For now, I leave this unparallelized
    // since this is just a training part.
    for (int i=0; i<numRefPoints; i++)
    {
        const Eigen::Vector3f p1 = sample_pc[i];
        const Eigen::Vector3f n1 = sample_nr[i];

        //printf("///////////////////// NEW REFERENCE ////////////////////////\n");
        for (int j=0; j<numRefPoints; j++)
        {
            // cannnot compute the ppf with myself
            if (i!=j)
            {
                const Eigen::Vector3f p2 = sample_pc[j];
                const Eigen::Vector3f n2 = sample_nr[j];

                Eigen::Vector4f f(0,0,0,0);
                computePPFFeaturesM2E(p1, n1, p2, n2, f);
                KeyType hashValue = hashPPFM2E(f, angle_step_radiansM2E, distanceStep);

                double alpha = computeAlphaM2E(p1, n1, p2);
                uint ppfInd = i*numRefPoints+j;

                THash* hashNode = &hash_nodesM2E[i*numRefPoints+j];
                hashNode->id = hashValue;
                hashNode->i = i;
                hashNode->ppfInd = ppfInd;
                hashNode->angle = (float)alpha;
                hashtableInsertHashed(hashTable, hashValue, (void*)hashNode);

                //cv::Vec4f matf(f(0,0), f(1,0), f(2,0), f(3,0));
                //Mat(matf).reshape(1, 1).convertTo(ppfM2E.row(ppfInd).colRange(0, 4), CV_32F);
                //ppfM2E.ptr<float>(ppfInd)[4] = (float)alpha;
            }
        }

        if(i%10 == 0)
        {
            std::cout << "trained: " << std::floor(100*i/numRefPoints) << "%" << std::endl;
        }
    }

    angle_stepM2E = angle_step_radiansM2E;
    distance_stepM2E = distanceStep;
    hash_tableM2E = hashTable;
    num_ref_pointsM2E = numRefPoints;
    trained = true;
}

void    ppf3d_public::matchM2E(std::vector<Eigen::Vector3f> S_pc, std::vector<Eigen::Vector3f> S_nr, std::vector<Pose_3DM2E>& results, const double relativeSceneSampleStep, const double relativeSceneDistance)
{
    scene_sample_stepM2E = (int)(1.0/relativeSceneSampleStep);

    int numAngles = (int) (floor (2 * M_PI / angle_stepM2E));
    float distanceStep = (float)distance_stepM2E;
    uint n = num_ref_pointsM2E;
    std::vector<Pose_3DM2E> poseList;
    int sceneSamplingStep = scene_sample_stepM2E;

    // compute bbox
    Eigen::Vector2f xRange, yRange, zRange;
    computeBBoxM2E(S_pc, xRange, yRange, zRange);

    // relativeSceneDistance = 0.04
    std::vector<Eigen::Vector3f> sampled_pc, sampled_nr;
    samplePCM2E(S_pc, S_nr, sampled_pc, sampled_nr, (float)relativeSceneDistance,0);

    // allocate the accumulator : Moved this to the inside of the loop
    /*#if !defined (_OPENMP)
     uint* accumulator = (uint*)calloc(numAngles*n, sizeof(uint));
    #endif*/

    poseList.reserve((sampled_pc.size()/sceneSamplingStep)+4);

#if defined _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < sampled_pc.size(); i += sceneSamplingStep)
    {
        uint refIndMax = 0, alphaIndMax = 0;
        uint maxVotes = 0;

        const Eigen::Vector3f p1  = sampled_pc[i];
        const Eigen::Vector3f n1  = sampled_nr[i];
        Eigen::Vector3f tsg(0,0,0);
        Eigen::Matrix3f Rsg = Eigen::Matrix3f::Zero(), RInv = Eigen::Matrix3f::Zero();

        uint* accumulator = (uint*)calloc(numAngles*n, sizeof(uint));
        computetransformRTM2E(p1, n1, Rsg, tsg);

        // Tolga Birdal's notice:
        // As a later update, we might want to look into a local neighborhood only
        // To do this, simply search the local neighborhood by radius look up
        // and collect the neighbors to compute the relative pose

        for (int j = 0; j < sampled_pc.size(); j ++)
        {
            if (i!=j)
            {
                const Eigen::Vector3f p2 = sampled_pc[j];
                const Eigen::Vector3f n2 = sampled_nr[j];

                Eigen::Vector3f p2t;
                double alpha_scene;

                Eigen::Vector4f f(0, 0, 0, 0);
                computePPFFeaturesM2E(p1, n1, p2, n2, f);
                KeyType hashValue = hashPPFM2E(f, angle_stepM2E, distanceStep);

                p2t = tsg + Rsg * Eigen::Vector3f(p2);

                alpha_scene=atan2(-p2t[2], p2t[1]);

                if ( alpha_scene != alpha_scene)
                {
                    continue;
                }

                if (sin(alpha_scene)*p2t[2]<0.0)
                    alpha_scene=-alpha_scene;

                alpha_scene=-alpha_scene;

                hashnode_i* node = hashtableGetBucketHashed(hash_tableM2E, (hashValue));

                while (node)
                {
                    THash* tData = (THash*) node->data;
                    int corrI = (int)tData->i;
                    int ppfInd = (int)tData->ppfInd;
                    //float* ppfCorrScene = ppfM2E.ptr<float>(ppfInd);
                    //double alpha_model = (float)ppfCorrScene[PPF_LENGTH-1];
                    double alpha_model = (float)tData->angle;
                    double alpha = alpha_model - alpha_scene;
                  /*
                  Tolga Birdal's note: Map alpha to the indices:
                  atan2 generates results in (-pi pi]
                  That's why alpha should be in range [-2pi 2pi]
                  So the quantization would be :
                  numAngles * (alpha+2pi)/(4pi)
                  */

                    //printf("%f\n", alpha);
                    int alpha_index = (int)(numAngles*(alpha + 2*M_PI) / (4*M_PI));

                    uint accIndex = corrI * numAngles + alpha_index;

                    accumulator[accIndex]++;
                    node = node->next;
                }
            }
        }

        // Maximize the accumulator
        for (uint k = 0; k < n; k++)
        {
            for (int j = 0; j < numAngles; j++)
            {
                const uint accInd = k*numAngles + j;
                const uint accVal = accumulator[ accInd ];
                if (accVal > maxVotes)
                {
                    maxVotes = accVal;
                    refIndMax = k;
                    alphaIndMax = j;
                }

#if !defined (_OPENMP)
                accumulator[accInd ] = 0;
#endif
            }
        }

        // invert Tsg : Luckily rotation is orthogonal:
        // Inverse = Transpose.
        // We are not required to invert

        Eigen::Vector3f tInv, tmg;
        Eigen::Matrix3f Rmg;
        RInv = Rsg.transpose();
        tInv = -RInv * tsg;

        Eigen::Matrix4f TsgInv;
        rtToPoseM2E(RInv, tInv, TsgInv);

        // TODO : Compute pose
        const Eigen::Vector3f pMax = sample_pc[refIndMax];
        const Eigen::Vector3f nMax = sample_nr[refIndMax];

        computetransformRTM2E(pMax, nMax, Rmg, tmg);

        Eigen::Matrix4f Tmg;
        rtToPoseM2E(Rmg, tmg, Tmg);

        // convert alpha_index to alpha
        int alpha_index = alphaIndMax;
        double alpha = (alpha_index*(4*M_PI))/numAngles-2*M_PI;

        // Equation 2:
        Eigen::Matrix4f Talpha;
        Eigen::Matrix3f R;
        Eigen::Vector3f t(0, 0, 0);
        getUnitXRotationM2E(alpha, R);
        rtToPoseM2E(R, t, Talpha);

        Eigen::Matrix4f rawPose = TsgInv * (Talpha * Tmg);

        Pose_3DM2E pose;
        pose.alpha = alpha;
        pose.modelIndex = refIndMax;
        pose.numVotes = maxVotes;
        pose.pose = rawPose;

#if defined (_OPENMP)
#pragma omp critical
#endif
        {
            poseList.push_back(pose);
        }

#if defined (_OPENMP)
        free(accumulator);
#endif
    }

    // TODO : Make the parameters relative if not arguments.
    //double MinMatchScore = 0.5;

    int numPosesAdded = sampled_pc.size()/sceneSamplingStep;

    std::cout << "numPoses = " << numPosesAdded << std::endl;

    clusterPosesM2E(poseList, numPosesAdded, results);
}

void    ppf3d_public::clusterPosesM2E(std::vector<Pose_3DM2E>& poseList, int numPoses, std::vector<Pose_3DM2E> &finalPoses)
{
    std::vector<ClusterM2E> poseClusters;

    finalPoses.clear();

    // sort the poses for stability
    std::sort(poseList.begin(), poseList.end(), pose_3DCompareM2E);

    for (int i=0; i<numPoses; i++)
    {
        Pose_3DM2E pose = poseList[i];
        bool assigned = false;

        // search all clusters
        for (size_t j=0; j<poseClusters.size() && !assigned; j++)
        {
            //const Pose_3D poseCenter = poseClusters[j]->poseList[0];
            const Pose_3DM2E poseCenter = poseClusters[j].poses[0];
            if (matchPoseM2E(pose, poseCenter))
            {
                //poseClusters[j]->addPose(pose);
                poseClusters[j].poses.push_back(pose);
                assigned = true;
                break;
            }
        }

        if (!assigned)
        {
            ClusterM2E poseCluster;
            poseCluster.poses.push_back(pose);
            poseClusters.push_back(poseCluster);
        }
    }

    // sort the clusters so that we could output multiple hypothesis
    std::sort(poseClusters.begin(), poseClusters.end(), Pose_3DClustersM2E);

    finalPoses.resize(poseClusters.size());

    // TODO: Use MinMatchScore

    if (use_weighted_avg)
    {
#if defined _OPENMP
#pragma omp parallel for
#endif
        // uses weighting by the number of votes
        for (int i=0; i<static_cast<int>(poseClusters.size()); i++)
        {
            // We could only average the quaternions. So I will make use of them here
            Eigen::Vector4f qAvg(0,0,0,0);
            Eigen::Vector3f tAvg(0,0,0);

            // Perform the final averaging
            ClusterM2E curCluster = poseClusters[i];
            std::vector<Pose_3DM2E> curPoses = curCluster.poses;
            int curSize = (int)curPoses.size();
            size_t numTotalVotes = 0;

            for (int j=0; j<curSize; j++)
                numTotalVotes += curPoses[j].numVotes;

            double wSum=0;

            for (int j=0; j<curSize; j++)
            {
                const double w = (double)curPoses[j].numVotes / (double)numTotalVotes;

                qAvg += w * curPoses[j].q;
                tAvg += w * curPoses[j].t;
                wSum += w;
            }

            tAvg *= 1.0 / wSum;
            qAvg *= 1.0 / wSum;


            //curPoses[0].updatePoseQuat(qAvg, tAvg);
            curPoses[0].numVotes= curCluster.accu_votes;

            finalPoses[i].pose =curPoses[0].pose;
            finalPoses[i].q =curPoses[0].q;
            finalPoses[i].t =curPoses[0].t;
            finalPoses[i].angle =curPoses[0].angle;
        }
    }
    else
    {
#if defined _OPENMP
#pragma omp parallel for
#endif
        for (int i=0; i<static_cast<int>(poseClusters.size()); i++)
        {
            // We could only average the quaternions. So I will make use of them here
            Eigen::Vector4f qAvg(0,0,0,0);
            Eigen::Vector3f tAvg(0,0,0);

            // Perform the final averaging
            ClusterM2E curCluster = poseClusters[i];
            std::vector<Pose_3DM2E> curPoses = curCluster.poses;
            const int curSize = (int)curPoses.size();

            for (int j=0; j<curSize; j++)
            {
                qAvg += curPoses[j].q;
                tAvg += curPoses[j].t;
            }

            tAvg *= 1.0 / curSize;
            qAvg *= 1.0 / curSize;

            //curPoses[0]->updatePoseQuat(qAvg, tAvg);
            curPoses[0].numVotes=curCluster.accu_votes;

            finalPoses[i].pose =curPoses[0].pose;
            finalPoses[i].q =curPoses[0].q;
            finalPoses[i].t =curPoses[0].t;
            finalPoses[i].angle =curPoses[0].angle;
        }
    }

    poseClusters.clear();
}

bool    ppf3d_public::matchPoseM2E(const Pose_3DM2E& sourcePose, const Pose_3DM2E& targetPose)
{
    // translational difference
    Eigen::Vector3f dv = targetPose.t - sourcePose.t;
    double dNorm = dv.norm();
    const double phi = fabs ( sourcePose.angle - targetPose.angle );

    return (phi<this->rotation_threshold && dNorm < this->position_threshold);
}

void    ppf3d_public::loadPLYFile(const char* fileName, std::vector<Eigen::Vector3f> &PC, std::vector<Eigen::Vector3f> &NOR)
{
    int numVertices = 0;
    int numCols = 3;
    int has_normals = 0;

    std::ifstream ifs(fileName);

    if (!ifs.is_open())
        std::cout << "Error opening input file: " << fileName << std::endl;

    std::string str;
    while (str.substr(0, 10) != "end_header")
    {
        std::vector<std::string> tokens = split(str,' ');
        if (tokens.size() == 3)
        {
            if (tokens[0] == "element" && tokens[1] == "vertex")
            {
                numVertices = atoi(tokens[2].c_str());
            }
            else if (tokens[0] == "property")
            {
                if (tokens[2] == "nx" || tokens[2] == "normal_x")
                {
                    has_normals = -1;
                    numCols += 3;
                }
                else if (tokens[2] == "r" || tokens[2] == "red")
                {
                    //has_color = true;
                    numCols += 3;
                }
                else if (tokens[2] == "a" || tokens[2] == "alpha")
                {
                    //has_alpha = true;
                    numCols += 1;
                }
            }
        }
        else if (tokens.size() > 1 && tokens[0] == "format" && tokens[1] != "ascii")
            std::cout << "Cannot read file, only ascii ply format is currently supported..." << std::endl;
        std::getline(ifs, str);
    }
    Eigen::Vector3f pt, nos;
    for (int i = 0; i < numVertices; i++)
    {
        int col = 0;
        for (; col < 6; ++col)
        {
            if(col < 3)
                ifs >> pt[col];
            else
                ifs >> nos[col-3];
        }
        for (; col < numCols; ++col)
        {
            float tmp;
            ifs >> tmp;
        }

        // normalize to unit norm
        double norm = sqrt(nos[0]*nos[0] + nos[1]*nos[1] + nos[2]*nos[2]);
        if (norm>0.00001)
        {
            nos[0]/=static_cast<float>(norm);
            nos[1]/=static_cast<float>(norm);
            nos[2]/=static_cast<float>(norm);
        }

        PC.push_back(pt);
        NOR.push_back(nos);
    }
}

void    ppf3d_public::savePLY(const std::string filename, std::vector<Eigen::Vector3f> &pts)
{
    std::ofstream file(filename.c_str(), std::ofstream::out);
    if(file.fail() == true)
    {
        std::cerr << filename << " cloud not be openned" << std::endl;
    }
    file << "ply" << std::endl;
    file << "format ascii 1.0" << std::endl;
    file << "comment VCGLIB generated" << std::endl;
    file << "comment VCGLIB generated" << std::endl;
    file << "element vertex " << pts.size() << std::endl;
    file << "property float x" << std::endl;
    file << "property float y" << std::endl;
    file << "property float z" << std::endl;
    file << "element face 0" << std::endl;
    file << "property list uchar int vertex_indices" << std::endl;
    file << "end_header" << std::endl;

    for (int i = 0; i < pts.size(); ++i)
    {
        file << std::setprecision(6)<< pts[i].x() << " " << pts[i].y() << " "<< pts[i].z() << /* " " << ptm[i].x() << " " << ptm[i].y() << " "<< ptm[i].z() << */std::endl;
    }

    file.close();

}

