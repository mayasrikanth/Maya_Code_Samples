#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include<fstream>
#include "transform.h"
#include<cmath>
#include "math.h"
// #include <cstdio>
#include<string>
#include <sstream>
using std::vector;
using namespace std;
using namespace Eigen; 

Matrix4f scalingMatrix(float x, float y, float z) {
	Matrix4f m; 
	m << x, 0, 0, 0, 
		 0, y, 0, 0, 
		 0, 0, z, 0,
		 0, 0, 0, 1; 
	return m; 
}
Matrix4f rotationMatrix(float x, float y, float z, float theta) {
	// ensure (x,y,z) is a unit vector
	// ASSUMING COUNTERCLOCKWISE ROTATION 
	Matrix4f m; 
	float x2_norm = pow(x,2.0);
	float y2_norm = pow(y, 2.0);
	float z2_norm = pow(z, 2.0);
	if(x2_norm + y2_norm  + z2_norm != 1) {
		float length = sqrt(x2_norm + y2_norm + z2_norm);
		x = x/length;
		y = y/length;
		z = z/length; 
	}
	float x2 = pow(x, 2.0);
	float y2 = pow(y, 2.0);
	float z2 = pow(z, 2.0);
	float test_norm = x2 + y2 + z2;
	// cout << "Normalized rotation vector: " << test_norm << "\n";
	float stheta = sin(theta);
	float ctheta = cos(theta);
	m << (x2 + (1 - x2)*ctheta), x*y*(1 - ctheta) - z*stheta, x*z*(1 -ctheta) + y*stheta, 0,
		 y*x*(1 - ctheta) + z*stheta, y2 + (1 - y2)*ctheta, y*z*(1-ctheta) - x*stheta, 0, 
		 z*x*(1 - ctheta) - y*stheta, z*y*(1 - ctheta) + x*stheta, z2 + (1 - z2)*ctheta, 0,
		 0, 0, 0, 1;
	return m; 
}
Matrix4f translationMatrix(float x, float y, float z) {
	Matrix4f m; 
	m << 1, 0, 0, x,
		 0, 1, 0, y,
		 0, 0, 1, z,
		 0, 0, 0, 1;
	// cout << "from translationMatrix: \n " << m << "\n"; 
	return m; 
}
/* Convert lines containing transformations to transformation matrix and
 * return result. */
Matrix4f getTransformationMatrix(std::vector<std::string> transformations) {
	// Vector to hold all transformation matrices
	std::vector<Matrix4f> matrices; 
	for(int i = 0; i < transformations.size(); i ++ ) {
		string c; // string to hold the transformation type 
		std::stringstream temp(transformations[i]); // line containing transformations 

		float x, y, z, theta; // doubles containing transformation info 
		temp >> c; 
		if(c == "t") {
			temp >> x >> y >> z; 
			// cout << "t " << x << " " << y << " " << z << "\n";
			matrices.push_back(translationMatrix(x, y, z));
		}
		else if (c == "r") {
			temp >> x >> y >> z >> theta; 
			// cout << "r " << x << " " << y << " " << z << " " << theta << "\n";
			matrices.push_back(rotationMatrix(x, y, z, theta));
		}
		else if (c == "s") {
			temp >> x >> y >> z; 
			// cout << "s " <<  x << " " << y << " " << z << "\n"; 
			matrices.push_back(scalingMatrix(x, y, z));
		}
	}

	Matrix4f res; 
	for(int i = 0; i < matrices.size(); i++) {
		if(i == 0) {
			res = matrices[i];
		}
		else {
			res = matrices[i] * res;
		}	 
	}
	return res; 	
}


