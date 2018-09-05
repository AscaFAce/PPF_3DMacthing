#ifndef STDAFX_H
#define STDAFX_H

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <algorithm>
#include <map>
#include <stdint.h>


typedef uint KeyType;

typedef struct hashnode_i
{
  KeyType key;
  void *data;
  struct hashnode_i *next;
} hashnode_i ;


typedef struct HSHTBL_i
{
  size_t size;
  struct hashnode_i **nodes;
  size_t (*hashfunc)(uint);
} hashtable_int;

typedef struct THash
{
  int id;
  int i, ppfInd;
  float angle;
} THash;


typedef struct Pose_3DM2E
{
    double alpha;
    double residual;
    double angle;
    size_t modelIndex;
    size_t numVotes;
    Eigen::Matrix4f pose;
    Eigen::Vector3f t;
    Eigen::Vector4f q;

} Pose_3DM2E;



typedef struct ClusterM2E
{
    std::vector<Pose_3DM2E> poses;
    size_t accu_votes;
} ClusterM2E;



#endif // STDAFX_H
