#include "makeScene.h"
#include "transform.h"  // in order to apply transformations to file 

#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include<fstream>
#include<cmath>
#include "math.h"
#include <cstdio>
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

// ALL FUNCTION HEADERS FOR FILE: 
/* Function to parse file and populate lights and camera for the scene. */
void lightsCam(std::ifstream &currFile, Camera &cam, vector<Light> &lights); 
/* Function to parse entire text file and populate scene accordingly. */ 
Scene setScene(std::ifstream &currFile, int xres, int yres, int shadingMode); 
/* Function to create object and apply geometric transformation */
void initModel(Matrix4f transformMatrix, Matrix4f transformNormals, std::string filename, int copyNum,
 std::string objName, vector<Model>& allObjects, Vector3f ambient, Vector3f diffuse, Vector3f specular, 
 float shininess);
/* Returns matrix to transform points to camera space */
Matrix4f cameraSpaceMatrix(Camera &currCam); 
/* Returns matrix to transform points to perspective space. */ 
Matrix4f perspectiveMatrix(Camera &currCam);


/* Applies specified geometric transformation to all vertices in model, and adds model to list 
 * of scene models. Also transforms the surface normals of this object by the inverse transpose
 * of the geometric transformaion matrix. This function ALSO saves the matterial lighting 
 * properties for the model, namely (ambient, diffuse, specular, shininess) */
void initModel(Matrix4f transformMatrix, Matrix4f transformNormals, std::string filename, int copyNum,
 std::string objName, vector<Model> &allObjects, Vector3f ambient, Vector3f diffuse, Vector3f specular, 
 float shininess) {
	char buffer[50]; 
	std::string currObjName; 
	if(copyNum > 0) {
		int n = sprintf(buffer, "%s_copy%d", objName.c_str(), copyNum);
		std::string str(buffer);
		// Assigning char array buffer to the currObjName string 
		currObjName = str;
	}
	// Create objInfo object to represent the file 
	Model currObj = Model(filename, currObjName);
	// Save material lighting properties 
	currObj.ambient = ambient;
	currObj.diffuse = diffuse;
	currObj.specular = specular;
	currObj.shininess = shininess;

	// cout << "ambient: " << ambient << "\n";
	// Transform all vertex coordinates 
	// std::vector<Vertex> newVertices;
	//Vector3f vertdummyVector; 
	// Create a dummy vertex at origin and push it 
	//Vertex dummyV = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, vertdummyVector, 0, 0};

	for(int i = 0; i < currObj.vertices.size(); i ++) {
		Vertex &currVertex = currObj.vertices[i];
		// cout << "Before geometric: " << "\n" <<  currVertex.worldVec() << " and ";
		// in-place transformation (passing geometric transformation)
		currVertex.transformVertex(transformMatrix); 
		// For the sake of the PHONG MODEL, we also need to 
		// transform the x,y,z,w fields geometrically 
		currVertex.x = currVertex.currX;
		currVertex.y = currVertex.currY;
		currVertex.z = currVertex.currZ;
		currVertex.w = currVertex.currW; 
	}
	
	// Transform and normalize all surface normals using the transformNormals matrix 
	// (each model has a vector of surface normals)
	 
	for(int i = 0; i < currObj.normals.size(); i ++) {
		SurfaceNormal &currNormal = currObj.normals[i];
		// in-place transformation (passing rotation / scaling transformation)
		currNormal.transformNormal(transformNormals); 
	}


	// CHECKING IF THINGS WERE SAVED AFTER INIT MODEL
		for(int i = 1; i < currObj.vertices.size(); i ++) {
		Vertex &currVertex = currObj.vertices[i];
	}
	// Updating allObjects vector containing all objects in scene 
	allObjects.push_back(currObj);
}

/* Returns matrix to transform points to camera space */
Matrix4f cameraSpaceMatrix(Camera &currCam) {
	vector<float> pos = currCam.position;
	if(pos.size() != 3) {
		cout << "Error, incorrect dimension of camera position" << endl;
		exit(1);
	}
	// x, y, z
	Matrix4f translation = translationMatrix(currCam.position[0], \
		currCam.position[1], currCam.position[2]); 

	vector<float> rot = currCam.orientation;
	if(rot.size() != 4) {
		cout << "Error, incorrect dimension of camera rotation" << endl;
		exit(1); 
	}
	// x, y, z, theta
	Matrix4f rotation = rotationMatrix(currCam.orientation[0], \
		currCam.orientation[1], currCam.orientation[2], currCam.orientation[3]); 
	// Return 
	Matrix4f C = translation * rotation;
	Matrix4f cResult = C.inverse();
	// cout << "C inverse: " << "\n" << cResult << "\n";
	return C.inverse(); 

}
/* Loads in the lights/ camera data from the file and stores in the lights vector/ camera struct
 * for the scene. */ 
void LightsCam(std::ifstream &currFile, Camera &cam, vector<Light> &lights) {
	/* Stop parsing when you reach "objects:" */
	std::string line; 
	if(currFile) {
		while(getline(currFile, line)) {
			char delim = ' ';
			stringstream countstream(line);
			stringstream linestream(line);
			string a; 
			std:: string token;
			std::vector<string> numWords;
			// Determine number of words in the line 
			while(std::getline(countstream, token, delim)) {
				numWords.push_back(token);
			}
			/* CAMERA */ 
			// CASE 1: "camera:"
			if(numWords.size() == 1) {
				string word; 
				linestream >> word; 
				if(word.compare("camera:") == 0) {
					continue; 
				}
				else if (word.compare("objects:") == 0) {
					return;  // 
				}
			}
			// CASE 2: near, far, left, right, top, bottom
			else if(numWords.size() == 2) {
				string fieldName;
				linestream >> fieldName; 
				float val;
				linestream >> val;
				if(fieldName.compare("near") == 0) {
					cam.n = val;
					//cout<< "Near: " << cam.n << "\n";
				}
				else if(fieldName.compare("far") == 0) {
					cam.f = val;
					//cout<< "Far: " << cam.f << "\n";
				}
				else if(fieldName.compare("left") == 0) {
					cam.l = val; 
					//cout << "Left: " << cam.l << "\n";
				}
				else if(fieldName.compare("right") == 0) {
					cam.r = val;
					//cout << "Right: " << cam.r << "\n";
				}
				else if (fieldName.compare("top") == 0) {
					cam.t = val;
					//cout << "Top: " << cam.t << "\n";
				}
				else if (fieldName.compare("bottom") == 0) {
					cam.b = val; 
					//cout << "Bottom: " << cam.b << "\n"; 
				}

			}
			// CASE 3: position 
			else if(numWords.size() == 4) {
				string fieldName; 
				linestream >> fieldName;
				if(fieldName.compare("position") == 0) {
					float x, y, z;
					linestream >> x >> y >> z;
					std::vector<float> pos{x, y, z};
					cam.position = pos; 
				}
			}
			// CASE 4: ORIENTATION 
			else if(numWords.size() == 5) {
				string fieldName;
				linestream >> fieldName;
				if(fieldName.compare("orientation") == 0) {
					float x, y, z, theta; 
					linestream >> x >> y >> z >> theta;
					std::vector<float> orient{x, y, z, theta};
					cam.orientation = orient; 
				}
			}	
			/* LIGHTS*/
			else if(numWords.size() == 10) {
				string word; 
				linestream >> word;
				// Then parse line, store in light struct,
				// push onto vector
				if(word.compare("light") == 0) {
					float x, y, z, r, g, b, k;
					char ch; 
					linestream >> x >> y >> z >> ch >>
					r >> g >> b >> ch >> k;
					Light newLight; 
					// Assigning position
					Vector3f position;
					position << x, y, z;
					newLight.position = position;

					// Assigning color 
					Vector3f color; 
					color << r, g, b;
					newLight.color = color;

					// Assigning attenuation 
					newLight.k = k; 

				// Store in vector
					lights.push_back(newLight); 
				}
			}
		}
	}
}


/* Parse txt file containing camera, object filenames, and object transformations.
 * Create a scene object that contains the camera, and a vector of objects. 
 * Each fileInfo object should contain a list of vertices and faces. */ 
Scene setScene(std::ifstream &currFile, int xres, int yres, int shadingMode) {
	// Holds mapping objName: fileName
	std::map<std::string, std::string> fileMappings; 
	// Holds mapping objName : copy number 
	std::map<std::string, int> objCopies; 
	// Vector to store the current line of geometric transformations for the current object
	std::vector<string> currTrans; 
	// Vector to store the current line of rotation and scaling transformations for 
	// the surface normals of the current object 
	std::vector<string> currNormTrans; 
	// Holds the mapping objName : currentObjReference 
	std::vector<Model> allObjects; 
	// Holds name of object currently being processed 
	std::string currObjName = "null"; 
	// Struct to hold camera information 
	Camera currCam; 
	// Vector of lights to hold all light structs
	vector<Light> currLights; 
	// String to hold current line being processed 
	string line; 
	// Pass the file to loadCam to process camera information 
	LightsCam(currFile, currCam, currLights);
	// Variables for ambience, diffuse, specular, shininess
	// for the object currently being processed 
	Vector3f currAmbiance; 
	Vector3f currDiffuse; 
	Vector3f currSpecular; 
	float currShininess; 
	// Load remainder of file 
	if(currFile) {
		while(getline(currFile, line)) {
			char delim = ' ';
			stringstream countstream(line);
			stringstream linestream(line);
			string a; 
			std:: string token;
			std::vector<string> numWords;
			// Determine number of words in the line 
			while(std::getline(countstream, token, delim)) {
				numWords.push_back(token);
			}
			// CASE 1: "objects: ", or objName 
			if(numWords.size() == 1) {
				string objName; 
				linestream >> objName; 
				// the objects: should not be appearing 
			  	if(objName.compare("objects:") != 0) {
					// Update the objCopies value 
					if (objCopies.find(objName) == objCopies.end()) {
						// Not found
						cout << "Error with line-parsing method" << endl;
						exit(1); 
					}
					else {
					 // Save current object name 
					currObjName = objName; 
					// Update copy number 
					objCopies[objName] += 1; 
					}
			  	}
			}
			// CASE 2: obj_name file_name 
			if(numWords.size() == 2) {
				// First check for shininess 
				if(numWords[0] == "shininess") {
					// Parse line 
					string s; 
					float p;
					linestream >> s >> p; 
					currShininess = p; 
				}

				// Otherwise obj_name file_name
				else {
					// Otherwise do this 
				string objName; 
				string fileName; 
				linestream >> objName >> fileName; 
				// Add objName : fileName mapping 
				fileMappings[objName] = fileName; 
				// Initialize objCopies value
				objCopies[objName] = 0; 
				}
				
			}
			// CASE 3: transformation line, so append line to vector holding translations 
			// for object being processed. Either this OR the material light properties.
			// To verify which one, we have some boolean logic. 
			else if (numWords.size() >= 4) {

				if(numWords[0] == "ambient") {
					// parse line -> eigen vector
					string name; 
					float r, g, b;
					linestream >> name >> r >> g >> b;
					Vector3f amb;
					amb << r, g, b;
					currAmbiance = amb; 
				}

				else if(numWords[0] == "diffuse") {
					// parse line -> eigen vector 
					string name;
					float r, g, b; 
					linestream >> name >> r >> g >> b;
					Vector3f diff; 
					diff << r, g, b;
					currDiffuse = diff;
				}

				else if(numWords[0] == "specular") {
					// parse line -> eigen vectors 
					string name; 
					float r, g, b; 
					linestream >> name >> r >> g >> b;
					Vector3f spec; 
					spec << r, g, b; 
					currSpecular = spec; 
				}

				// Otherwise we have a translation vector. 
				else {
					currTrans.push_back(line); 
					// Append to currNormTrans if scaling or rotation 
					if(numWords[0] == "s" or numWords[0] == "r") {
						currNormTrans.push_back(line); 
					}
				}
			}
		
			// Matrix4f transformMatrix,std::string filename, int copyNum, std::string objName, std::vector<objInfo> &allObjects
			// CASE 4: empty line 
			else if (numWords.size() == 0) {
				// If currTrans.size() > 0, apply all transformations to currObj 
				if(currTrans.size() > 0) {
					// Retrieve transformation matrix 
					Matrix4f geometricTransformation = getTransformationMatrix(currTrans);
					Matrix4f geoTransNormals = getTransformationMatrix(currNormTrans); 
					// Retrieve the copy # we are on for this object 
					int copyNum = objCopies[currObjName];
					std::string filename = fileMappings[currObjName];
					// Apply all specified transformations to the vertices in this object
					initModel(geometricTransformation, geoTransNormals, filename, copyNum, currObjName, allObjects, 
						currAmbiance, currDiffuse, currSpecular, currShininess); 
					// Clear the transformation vectors, as this object is done processing 
					std::vector<string>().swap(currTrans);
					std::vector<string>().swap(currNormTrans);
					// reset currObjName to "null"
					currObjName = "null"; 
				}
			}
		}
	}

		// If currTrans.size() > 0, apply all transformations to currObj 
		if(currTrans.size() > 0) {
			// Retrieve transformation matrix 
			Matrix4f geometricTransformation = getTransformationMatrix(currTrans);
			Matrix4f geoTransNormals = getTransformationMatrix(currNormTrans); 
			// Retrieve the copy number we are on for this object 
			int copyNum = objCopies[currObjName];
			std::string filename = fileMappings[currObjName];
			// This function applies the geometric transformation to all vertices in the object,
			// as well as the appropriate transformation to all the surface normals in the object. 
			initModel(geometricTransformation, geoTransNormals, filename, copyNum, currObjName, allObjects, 
					currAmbiance, currDiffuse, currSpecular, currShininess); 
			// Clear the transformation vectors  
			std::vector<string>().swap(currTrans);
			std::vector<string>().swap(currNormTrans); 
			// reset currObjName to "null"
			currObjName = "null";
		}
		currFile.close(); 
	
	// NOW: create scene object, as you have all the information necessary 
	// to finish constructing the scene. 
	Scene currScene; 
	currScene.cam = currCam; 
	currScene.sceneLights = currLights;
	currScene.sceneModels = allObjects;
	// Assign shading mode 
	currScene.shadingMode = shadingMode;
	// Assign xres, yres
	currScene.xres = xres;
	currScene.yres = yres; 

	// COMPUTE AND STORE CAMERA TRANFORM MATRIX FOR SCENE 
	currScene.cameraTransMatrix = cameraSpaceMatrix(currCam); 

	// Right here, set the camera coordinates to cameraTransMatrix * origin point 
	Vector4f origin;
	origin << 0, 0, 0, 1.0;

	// Set this new point as the camera position (divide by homogenous component
	// and set the x,y,z as new camera position)
	// cout << "Old camera: " << currScene.cameraTransMatrix << "\n";  
	Matrix4f invCam = currScene.cameraTransMatrix;
	Matrix4f invCamera = invCam.inverse(); 
	Vector4f cameraCoords = invCamera * origin;
	float camX = cameraCoords(0);
	float camY = cameraCoords(1);
	float camZ = cameraCoords(2);
	float camW = cameraCoords(3);
	camX = camX / camW;
	camY = camY / camW;
	camZ = camZ / camW;

	// create vector
	vector<float> newCamCoords {camX, camY, camZ};
	// Set the field for camera 
	currScene.cam.position = newCamCoords;
	// COMPUTE AND STORE PERSPECTIVE TRANFORM MATRIX FOR SCENE 
	currScene.perspectiveTransMatrix = perspectiveMatrix(currCam);  
	// cout << "PerspectiveTransMatrix: \n" << currScene.perspectiveTransMatrix << "\n"; 
	// Return currScene to calling function (main) in shadedRenderer
	return currScene;  

}




/* Returns matrix to transform points to perspective space. */ 
Matrix4f perspectiveMatrix(Camera &currCam) {
	Matrix4f persMatrix; 
	float n = currCam.n;
	float r = currCam.r;
	float l = currCam.l;
	float f = currCam.f;
	float b = currCam.b; 
	float t = currCam.t; 

	persMatrix << (2*n)/(r - l), 0, (r + l)/(r - l), 0, 
				   0, (2*n)/(t - b), (t + b)/(t - b), 0,
				   0,  0, -(f + n)/(f - n), (-2 * f * n)/(f - n),
				   0,  0,  -1, 0; 
	return persMatrix; 
}







