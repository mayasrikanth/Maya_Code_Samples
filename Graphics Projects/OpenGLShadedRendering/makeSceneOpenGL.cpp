/* Create a Scene Object for the openGL shaded rendering
 * implementation. Pass this Scene Object to openGLRenderer's
 * main. This code includes Model-creation, file-parsing, line-parsing,
 * and scene creation.   */

#include "makeSceneOpenGL.h"
// No longer need to manually apply transforms  
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


/* Implementation of non-default constructor for Model class (located in makeSceneOpenGL.h) */
Model::Model(std::string f, std::string name){
		// FINDING THE CORRECT DIRECTORY
		// Save object name of form obj_copyx
		objName = name;
		std::string filename; 
		std::stringstream ss(f); 
		std::string tok; 
		vector<string> filepath; 
		// Were we given a file or a path to a file?
		while(std::getline(ss, tok, '/' )) {
			filepath.push_back(tok); 
		}
		// If we have a path, then we need to extract the directory and use chdir
		if(filepath.size() > 1) {
			std::size_t found = f.find_last_of("/");
			// Now change current working directory 
			std::string workingDir = f.substr(0, found); 
			const char *p = workingDir.c_str(); 
			chdir(p); 
			// Set the filename to everything from \-> end (excluding last slash)
			filename = f.substr(found + 1);
		}
		// We are given only a file
		else if(filepath.size() == 1) {
			filename = f;
		}
		// FILE PARSING 
		// Open object file for reading or writing 
		std::fstream infile (filename);
		vector<Triple> verts; 
		vector<Triple> norms;
		vector<Face> face;

		std::string line;

		// Create a dummy vertex and normal 
		Triple dummyV = {0.0, 0.0, 0.0};
		Triple dummyS = {0.0, 0.0, 0.0}; 

		verts.push_back(dummyV);
		norms.push_back(dummyS);

		char i;
		char ch; 
		float a, b, c; 
		// First parse first string, and from there decide how to 
		// parse the line depending on v, vn, or f
		if(infile) {
			while(getline(infile, line)) {
				char delim = ' ';
				stringstream countstream(line);
				stringstream linestream(line);
				string a; 
				std:: string token;
				std::vector<string> numWords;
				// Split line into strings to determine first string
				// in line 
				while(std::getline(countstream, token, delim)) {
					numWords.push_back(token);
				}
				if(numWords[0] == "v") {
					float x, y, z;
					string v;
					linestream >> v >> x >> y >> z;
					Triple newVert = {x, y, z};
					verts.push_back(newVert);
				}
				// SurfaceNormal
				else if(numWords[0] == "vn") {
					// Store for i, currI 
					float x, y, z; 
					string vn;
					linestream >> vn >> x >> y >> z;
					Triple newNorm = {x, y, z};
					norms.push_back(newNorm);
				}
				// Face with vertexIdx//normalIdx
				else if(numWords[0] == "f") {
					int v1, s1, v2, s2, v3, s3;
					char ch; 
					// First pair
					stringstream p1(numWords[1]);
					p1 >> v1 >> ch >> ch >> s1;
					// Second pair
					stringstream p2(numWords[2]);
					p2 >> v2 >> ch >> ch >> s2;
					// Third pair 
					stringstream p3(numWords[3]);
					p3 >> v3 >> ch >> ch >> s3;
					// Now create face 
					Face newFace = {0, 0, 0, 0, 0, 0};
					newFace.v1 = v1;
					newFace.v2 = v2;
					newFace.v3 = v3;
					newFace.s1 = s1;
					newFace.s2 = s2;
					newFace.s3 = s3;
					face.push_back(newFace);
				}
			}
		}
		// Saving vector fields 
		vertices = verts;
		faces = face; 
		normals = norms;
		filename = f; 	
		// Closing file after operations
		infile.close();
		// Populate buffers for this particular Model
		this->ModelSetBuffers();
	}
/* Implementation of function to populate the normal and vertex buffers. */
void Model::ModelSetBuffers() {
	// Ensure model vectors are populated, if not exit with error message 
	if (vertices.size() == 0 or faces.size() == 0 or normals.size() == 0) {
		cout << "The vertices, faces, or normals for this object are not populated!" << "\n";
		exit(EXIT_FAILURE);
	}
	vector<Triple> normbuff; 
	vector<Triple> vertbuff; 

	// Populate vertex buffer and normals buffer (v1 -> s1, v2 -> s2, v3 -> s3)
	for(int i = 0; i < faces.size(); i ++) {
		// Extract indices of normals and vertices 
		int v1 = faces[i].v1;
		int v2 = faces[i].v2;
		int v3 = faces[i].v3; 
		int s1 = faces[i].s1;
		int s2 = faces[i].s2;
		int s3 = faces[i].s3;

		// Extract vertices for the face
		Triple &vert1 = vertices[v1];
		Triple &vert2 = vertices[v2];
		Triple &vert3 = vertices[v3];

		vertbuff.push_back(vert1);
		vertbuff.push_back(vert2);
		vertbuff.push_back(vert3);

		// Extract corresponding normals for the face
		Triple &norm1 = normals[s1];
		Triple &norm2 = normals[s2];
		Triple &norm3 = normals[s3];

		normbuff.push_back(norm1);
		normbuff.push_back(norm2);
		normbuff.push_back(norm3);
	}
	// Now assign vertex_buffer and normal_buffer
	vertex_buffer = vertbuff;
	normal_buffer = normbuff; 
}


/* Initializes model with vertex, normal, and face info. Also saves material 
 * lighting properties for the model: (ambient, diffuse, specular, shininess).  */
void initModel(vector<string> objectTransformations, std::string filename, int copyNum,
 std::string objName, vector<Model> &allObjects, Vector3f &ambient, Vector3f &diffuse, Vector3f &specular, 
 float shininess) {
	char buffer[50]; 
	std::string currObjName; 
	if(copyNum > 0) {
		int n = sprintf(buffer, "%s_copy%d", objName.c_str(), copyNum);
		std::string str(buffer);
		// Assigning char array buffer to the currObjName string 
		currObjName = str;
	}
	// Create objInfo object to represent the file (should populate vertex_buffer
	// and normal_buffer for the Model)
	Model currObj = Model(filename, currObjName);
	// Reverse string vector containing transformations for this object
	reverse(objectTransformations.begin(), objectTransformations.end());
	// Assiging string vector of transformations 
	currObj.modelTransformations = objectTransformations;
	// Save material lighting properties (dereferencing pointers to get value)
	currObj.ambient_reflect[0] = ambient(0);
	currObj.ambient_reflect[1] = ambient(1);
	currObj.ambient_reflect[2] = ambient(2);

	currObj.diffuse_reflect[0] = diffuse(0);
	currObj.diffuse_reflect[1] = diffuse(1);
	currObj.diffuse_reflect[2] = diffuse(2);

	currObj.specular_reflect[0] = specular(0);
	currObj.specular_reflect[1] = specular(1);
	currObj.specular_reflect[2] = specular(2);

	currObj.shininess = shininess;

	// Updating allObjects vector containing all objects in scene 
	allObjects.push_back(currObj);
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
					float mag = sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0));
					x = x/mag;
					y = y/mag;
					z = z/mag;
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
					// Vector3f position;
					// position << x, y, z;
					newLight.position[0] = x;
					newLight.position[1] = y;
					newLight.position[2] = z;
					// float position[] = {x, y, z, 1.0}; 
					// newLight.position = position;

					// Assigning color 
					Vector3f color; 
					color << r, g, b;

					newLight.color[0] = r;
					newLight.color[1] = g;
					newLight.color[2] = b;

					//float color[] = {r, g, b};
					//newLight.color = *color;

					// Assigning attenuation 
					newLight.k = k; 

					// Store in vector
					lights.push_back(newLight); 
				}
			}
		}
	}
}

/* Set the scene. */
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
	Vector3f currAmbient; 
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
					//currAmbiance[0] = r;
					//currAmbiance[1] = g;
					//currAmbiance[2] = b;
					Vector3f amb;
					amb << r, g, b;
					currAmbient = amb; 
				}

				else if(numWords[0] == "diffuse") {
					// parse line -> eigen vector 
					string name;
					float r, g, b; 
					linestream >> name >> r >> g >> b;
					//currDiffuse[0] = r;
					//currDiffuse[1] = g;
					//currDiffuse[2] = b;
					Vector3f diff; 
					diff << r, g, b;
					currDiffuse = diff;
				}

				else if(numWords[0] == "specular") {
					// parse line -> eigen vectors 
					string name; 
					float r, g, b; 
					linestream >> name >> r >> g >> b;
					//currSpecular[0] = r;
					//currSpecular[1] = g;
					//currSpecular[2] = b; 
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
					// Retrieve the copy # we are on for this object 
					int copyNum = objCopies[currObjName];
					std::string filename = fileMappings[currObjName];
					// Initialize model and save tranformations 
					initModel(currTrans, filename, copyNum, currObjName, allObjects, 
						currAmbient, currDiffuse, currSpecular, currShininess); 
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
			// Retrieve the copy number we are on for this object 
			int copyNum = objCopies[currObjName];
			std::string filename = fileMappings[currObjName];
			// Initialize the model  
			initModel(currTrans, filename, copyNum, currObjName, allObjects, 
					currAmbient, currDiffuse, currSpecular, currShininess); 
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

	return currScene;  

}