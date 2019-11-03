
#ifndef QUATERNION_H
#define QUATERNION_H

#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include<fstream>
#include<cmath>
#include "math.h"
#include <cstdio>
// For exiting error 
#include <stdlib.h>

// For changing working directory 
#include <stdio.h>
#include <unistd.h>
#include <cstddef>

using std::vector;
using namespace std;
using namespace Eigen;
#include <string>
#include <sstream>
#include<map>


// Declare Quaternion class 
class quaternion {

public: 
	// Eigen vector to store the i,j,k components
	Vector3f imaginary; 
	// real component 
	float s; 

	// Matrix representation of the rotation 
	//Matrix4f rotationMatrix; 
	// Default constructor
	quaternion();
};

// Functions in Quaternion.cpp
Vector3f screenToNDC(int x_screen, int y_screen, int xres, int yres);
quaternion ComputeRotationQuaternion(int x_s, int y_s, int x_e, int y_e, int xres, int yres);
float computeRotationAngle(Vector3f &p_s, Vector3f &p_e);
quaternion quaternionProduct(quaternion A, quaternion B);
Matrix4d quaternionRotationMatrix(quaternion quat);
Vector3f computeRotationAxis(Vector3f &p_s, Vector3f &p_e);
quaternion IdentityQuaternion();
quaternion normalizeQuaternion(quaternion X);

#endif