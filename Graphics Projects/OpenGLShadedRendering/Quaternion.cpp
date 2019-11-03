/* class that creates a quaterion  */
#include "Quaternion.h"
#include<cmath>
#include "math.h"
#include <Eigen/Dense>
using std::vector;
using namespace std;
using namespace Eigen;

Vector3f screenToNDC(int x_screen, int y_screen, int xres, int yres);
void ComputeRotationquaternion(int x_s, int y_s, int x_e, int y_e, int xres, int yres);
float computeRotationAngle(Vector3f &p_s, Vector3f &p_e);
quaternion quaternionProduct(quaternion R, quaternion L);
Matrix4d quaternionRotationMatrix(quaternion quat);
Vector3f computeRotationAxis(Vector3f &p_s, Vector3f &p_e);
quaternion Identityquaternion();
quaternion normalizeQuaternion(quaternion X);


// Default constructor : Default is identity quaternion 
quaternion::quaternion() {
	s = 1.0; 
	imaginary << 0.0, 0.0, 0.0;
}
// CALLED BY mouse_moved IF pressed is true in order to 
// "cache" all rotation data 
quaternion ComputeRotationQuaternion(int x_s, int y_s, int x_e, int y_e, int xres, int yres) {
	// Convert both points from screen -> NDC 
	Vector3f p_s = screenToNDC(x_s, y_s, xres, yres);
	Vector3f p_e = screenToNDC(x_e, y_e, xres, yres);

	// Retrieve unit vector for rotation 
	Vector3f unitRotation = computeRotationAxis(p_s, p_e);
	float theta = computeRotationAngle(p_s, p_e);
	// If y_s < y_e
	// theta *= -1; ss

	// Construct quaternion t
	quaternion rotationquaternion; 
	float real = cos((theta / 2.0)); 
	unitRotation *= (sin(theta / 2.0)); 

	rotationquaternion.s = real;
	rotationquaternion.imaginary = unitRotation; 
	// Returning quaternion to mouse-moved 
	return rotationquaternion; 
}

// Function that multiplies two quaternions and returns the product as 
// another quaternion (order: A x B)
quaternion quaternionProduct(quaternion A, quaternion B) {
	//cout << "A:  " << A.s << " \n" << A.imaginary << "\n";
	//cout << "B: " << B.s << " \n" << B.imaginary << "\n";
	float real = (A.s * B.s) - (A.imaginary.dot(B.imaginary));
	Vector3f imag = (A.s * B.imaginary) + (B.s *A.imaginary) 
	+ (A.imaginary.cross(B.imaginary));

	quaternion product; 
	product.s = real;
	product.imaginary = imag; 

	//cout << "Result of quaternionProduct: " << product.s << "  \n" 
	//<< product.imaginary << "\n"; 
	// Returning quaternion representing product of two 
	return product; 

}
quaternion normalizeQuaternion(quaternion X) {
	Vector3f imag = X.imaginary;
	float mag = sqrt(pow(X.s, 2.0) + pow(imag(0), 2.0) + pow(imag(1), 2.0) + pow(imag(2), 2.0));
	X.s = X.s / mag; 
	X.imaginary = X.imaginary * (1/mag);
	return X;

}
// Create and return a rotation matrix using the fields of the 
// rotationquaternion passed in. 
Matrix4d quaternionRotationMatrix(quaternion quat) {
	double q_s = quat.s; 
	Vector3f imag = quat.imaginary;
	// cout << "Real: " << q_s << "Imaginary: \n" << imag << "\n";
	double q_x = imag(0);
	double q_y = imag(1);
	double q_z = imag(2);

	double q_x2 = pow(q_x, 2.0);
	double q_y2 = pow(q_y, 2.0);
	double q_z2 = pow(q_z, 2.0);

	Matrix4d rotationMatrix; 
	rotationMatrix << (1 - 2*q_y2 - 2*q_z2), 2*(q_x*q_y - q_z*q_s), 2*(q_x*q_z + q_y*q_s), 0,
	2*(q_x*q_y + q_z*q_s), (1 - 2*q_x2 - 2*q_z2), 2*(q_y*q_z - q_x*q_s), 0,
	2*(q_x*q_z - q_y*q_s), 2*(q_y*q_z + q_x*q_s),  (1 - 2*q_x2 - 2*q_y2), 0, 
	0, 0, 0, 1;

	// cout << "RotationMatrix (from Quaternion: ) " << rotationMatrix << "\n";
	return rotationMatrix; 

}

// Returns identity quaternion, potentially very inefficient
quaternion IdentityQuaternion() {
	// returning identity quaternion 
	quaternion Identity;  // default constructor creates identity quaternion
	Identity.s = 1.0;
	Vector3f imagIdentity;
	imagIdentity << 0, 0, 0; 
	Identity.imaginary = imagIdentity; 
	return Identity; 
}

// Returns rotation angle (in radians)
float computeRotationAngle(Vector3f &p_s, Vector3f &p_e) {
	float numerator = p_s.dot(p_e);
	float mag_s = p_s.norm();
	float mag_e = p_e.norm();

	float arg = numerator * (1 / (mag_s * mag_e));

	// Returning arccos in RADIANS
	//float res = acos(arg);
	//float alternative = acos(1);
	float one = 1.0;
	float param = min(arg, one);
	float res = acos(param); 

	return res; 
}

// Return unit vector representing rotation 
Vector3f computeRotationAxis(Vector3f &p_s, Vector3f &p_e) {
	// p_start x p_end
	Vector3f rotation = p_s.cross(p_e);
	rotation.normalize();
	return rotation; 
}

// Returns an (x,y,z) vector representing coordinates mapped to 1x1x1 NDC sphere
Vector3f screenToNDC(int x_screen, int y_screen, int xres, int yres) {
	float x_ndc, y_ndc, z_ndc; 
	// y: [-1, 1] -> [800,0]
	y_ndc =  1 + ((-2 * (float)y_screen) / ((float)yres - 1)); 

	// x: [-1, 1] -> [0, 800]
	x_ndc = -1 + ((float)(x_screen)/(float)(xres - 1)) * 2;

	// Determine z-component 
	float x2 = pow(x_ndc, 2.0);
	float y2 = pow(y_ndc, 2.0);
	float sumxy = x2 + y2; 
	if(sumxy <= 1) {
		z_ndc = sqrt((1 - x2 - y2));
	}
	else if(sumxy > 1) {
		z_ndc = 0; 
	}
	Vector3f ndcPoint; 
	ndcPoint << x_ndc, y_ndc, z_ndc; 
	return ndcPoint; 
}

// Note, Draw_Scene will call glMultMatrixd, whose argument is a 
// pointer to a 4x4 column matrix. You will pass the final rotation
// matrix into this funciton. 



// Workflow: you have p'_s= (x'_s, y'_s) and 
// p'_e = (x'_e, y'_e). Both can be obtained in 
// the mouse_pressed function and both are screen coordinates. 


// (1) Map both p'_s and p'_e to NDC, ignoring z component. 
// Depending on whether the mapped points are on the NDC sphere,
// the z component is 0 or 1. 

// Represent these NDC-variant points as eigen Vector3f

// Having constructed these points, it is time to construct the quaternion. 

// (1) Find Rotation angle theta (using arc cos,dot product, magnitude)
		// REstrict to arccos(min(1,arg))
// (2) Find axis u to rotate about using cross product and NORMALIZE u. 


// Addition 
// Subtraction 
// Scalar multiplication 




// Big idea: you can multiply two quaternions together to get the most 
// updated rotation. 


