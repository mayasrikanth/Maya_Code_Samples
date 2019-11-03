/* Header file for the openGL version of MakeScene.cpp */
#include<fstream>
#include <vector>
#include<iostream>

#ifndef MAKESCENEOPENGL_H
#define MAKESCENEOPENGL_H

// including eigen for header definitions 
#include <Eigen/Dense>
using namespace Eigen;
#include <iostream>

// for stringstream: 
#include<sstream>
// for chdir:
#include <unistd.h>
using std::vector;
using namespace std;

// Triple struct to store (x,y,z) for normals and vertices. 
struct Triple {
	float x;
	float y;
	float z;
};
	
// Holds point light attributes, store as float vectors. 
struct Light {
	// vector<float> position;
	float position[4];
	// vector<float> color; 
	float color[3];
	// Attenuation parameter
	float k; 
};

// Contains indices of vertices for face
struct Face {
	// Indices for vertices that comprise a face 
	int v1, v2, v3; 
	// Indices for surface normals that correspond to 
	// v1, v2, and v3
	int s1, s2, s3; 
}; 

// Contains camera information 
struct Camera {
	// (x,y,z) position of camera 
	vector<float> position; 
	// (x,y,z,theta) orientation of camera 
	vector<float> orientation;
	// Frustum parameters
	float n; 
	float f; 
	float l;
	float r;
	float t;
	float b; 
};

class Model
{
	/* Represents object given by .obj file */
public:
	// Material lighting properties
	float ambient_reflect[3];
	float diffuse_reflect[3];
	float specular_reflect[3];
	float shininess; 

	/* Model-related information. */ 
	// Vector of vertices 
	std::vector<Triple> vertices; 
	// Vector of normals 
	// Vector of faces where each face contains 3 vertices with corresponding 
	// surface normals
	std::vector<Face> faces;
	std::vector<Triple> normals; 

	/* Use the 3 vectors above to populate: */
	// Vector of all vertices in face-arrangement (3 points per face)
	std::vector<Triple> vertex_buffer;
	// Vector of all normals 
	std::vector<Triple> normal_buffer; 
	
	// Vector to hold transformation lines in reverse order 
	std::vector<string> modelTransformations; 

	/* File-related Information. */
	// Name of objfile to parse 
	std::string filename;  
	// Name of object 
	std::string objName;

	std::vector<bool> keepPoints;

	// Declaring non-default constructor 
	Model(std::string f, std::string name);
	// Function to populate the vertex_buffer and normal_buffer for this Model
	void ModelSetBuffers(); 
};

class Scene {
	/* Container for all models, camera, and 
	* transformations for a given scene. */
public:
	int xres;
	int yres;
	Camera cam; 
	// Contains all lights in the scene 
	vector<Light> sceneLights; 
	// Contains all models in the scene 
	vector<Model> sceneModels;
	// Mode: 1 = Phong Shading, 0 = Gourard Shading
	int shadingMode; 
};

// Declare all function headers for makeScene.cpp
/* Loads in the lights/ camera data from the file and stores in the lights vector/ camera struct
 * for the scene. */ 
void lightsCam(std::ifstream &currFile, Camera &cam, vector<Light> &lights);

/* Initializes model with vertex, normal, and face info. Also saves material 
 * lighting properties for the model: (ambient, diffuse, specular, shininess).  */
void initModel(vector<string> objectTransformations, std::string filename, int copyNum,
 std::string objName, vector<Model> &allObjects, float ambient[], float diffuse[], float specular[], 
 float shininess);

 /* Parse txt file containing camera, object filenames, and object transformations.
 * Create a scene object that contains the camera, and a vector of objects. 
 * Each fileInfo object should contain a list of vertices and faces. */ 
Scene setScene(std::ifstream &currFile, int xres, int yres, int shadingMode); 

#endif