#include <iostream>
#include <string>
#include <sstream>
#include<vector> 
using std::vector; 
#include <Eigen/Dense>
using namespace std; 
using namespace Eigen; 
// Need to include makeScene.h in order to access 
// camera struct 




Matrix4f scalingMatrix(float x, float y, float z);

Matrix4f rotationMatrix(float x, float y, float z, float theta);

Matrix4f translationMatrix(float x, float y, float z);

Matrix4f getTransformationMatrix(std::vector<string> transformations); 

