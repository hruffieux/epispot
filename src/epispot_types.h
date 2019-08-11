#ifndef EPISPOT_TYPES_H_
#define EPISPOT_TYPES_H_

#include <RcppEigen.h>

// These typedefs have to be in a separate header file to be properly copied to RcppExports

typedef Eigen::ArrayXd Arr1D;
typedef Eigen::ArrayXXd Arr2D;
typedef Eigen::Map<Eigen::ArrayXd> MapArr1D;
typedef Eigen::Map<Eigen::ArrayXXd> MapArr2D;
typedef Eigen::Map<Eigen::VectorXd> MapVec;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;


#endif // EPISPOT_TYPES_H_
