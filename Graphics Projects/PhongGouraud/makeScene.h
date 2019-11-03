/* Holds all function declarations for makeScene.cpp */
#include<fstream>
#include <vector>
#include<iostream>

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

// Holds x,y,z components of surface normal 
struct SurfaceNormal {
	float x;
	float y;
	float z;
	// Function to put the most-recently transformed 
	// coordinates of surface normal into a eigen vector. 
	Vector3f Vectorize() {
		Vector3f vectorRep; 
		float xtemp = this->x;
		float ytemp = this->y;
		float ztemp = this->z;
		vectorRep << x, y, z; 
		return vectorRep; 
	}

	// Function to manually set the fields with an eigenvector
	void setNormal(Vector3f normal) {
		this->x = normal(0);
		this->y = normal(1);
		this->z = normal(2);
	}
	// Transform and normalize surface normal vector,
	// assuming w = 1. This function takes a transformation 
	// matrix that only accounts for rotation and scaling. 
	void transformNormal(Matrix4f &transformMatrix) {
		// Inverse normal 
		Matrix4f geoInv = transformMatrix.inverse(); 
		// Transpose
		Matrix4f invTrans = geoInv.transpose(); 
		// Now Transform 
		Vector4f currNorm; 
		currNorm << this->x , this->y, this->z, 1.0;  
		// Multiply with invTrans
		Vector4f transformedNormals = invTrans * currNorm;
		// Extract x,y,z componens 
		float normx = transformedNormals(0);
		float normy = transformedNormals(1);
		float normz = transformedNormals(2);
		// Technically the w component should be 1
		float normw = transformedNormals(3); 
		// Normalize by dividing by the magnitude (x,y,z) since the 
		// lighting and shading algorithms assume unit normals 
		float mag = sqrt(pow(normx, 2) + pow(normy, 2) + pow(normz, 2));
		normx = normx / mag;
		normy = normy / mag;
		normz = normz / mag; 
		// Re-assign the surface normal 
		this->x = normx; 
		this->y = normy;
		this->z = normz; 
	}


};
// Holds point light attributes, storing as eigenvectors to facilitate
// computation later on in lighting and shading algorithms  
struct Light {
	Vector3f position;
	Vector3f color; 
	// Attenuation parameter
	float k; 
};

struct Material {
	// r,g,b values for the material reflectance property
	float r;
	float g;
	float b;
};
// Contains vertex coordinates 
struct Vertex {
	// Original world-space coordinates transformed geometrically. 
	// No other transformations will be applied to these fields.  
	float x;
	float y;
	float z;
	float w;
	// Current transformation (this can be
	// world, geometric, camera, perspective, ndc)
	float currX;
	float currY;
	float currZ;
	float currW; 

	
	// Function to return world-coordinates in vector form
	// in-place transformation of a vertex
	Vector4f worldVec() {
		Vector4f vectorRep; 
		float xtemp = this->x;
		float ytemp = this->y;
		float ztemp = this->z;
		// Should be 1.0 if NDC
		float wtemp = this->w;
		vectorRep << xtemp, ytemp, ztemp, wtemp; 
		return vectorRep; 

	}
	void transformVertex(Matrix4f &transMatrix) {
		Vector4f currPoint;
		currPoint << this->currX, this->currY, this->currZ, 
			 this->currW; 
		Vector4f transPoint = transMatrix * currPoint;
		this->currX = transPoint(0);
		this->currY = transPoint(1);
		this->currZ = transPoint(2);
		this->currW = transPoint(3);	 
	}
	void divideByW(){
		float homog = this->currW;
		float xtemp = this->currX;
		float ytemp = this->currY;
		float ztemp = this->currZ;
		this->currX = (xtemp / homog);
		this->currY = (ytemp / homog);
		this->currZ = (ztemp / homog); 
		this->currW = 1.0; 
	}
	// Function to put the most-recently transformed 
	// coordinates into a eigen vector. 
	Vector4f Vectorize() {
		Vector4f vectorRep; 
		float xtemp = this->currX;
		float ytemp = this->currY;
		float ztemp = this->currZ;
		// Should be 1.0 if NDC
		float wtemp = this->currW;
		vectorRep << xtemp, ytemp, ztemp, wtemp; 
		return vectorRep; 
	}
	// Color value for world-space components 
	// Colors are of form [red, green, blue]
	Vector3f color; 
	// Finally, store screen coordinates 
	int xscreen;
	int yscreen;
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
	vector<float> position; 
	// Function to return eigen vector representation of camera
	// position 
	Vector3f vectorizePos() {
		Vector3f res; 
		float x = this->position[0];
		float y = this->position[1];
		float z = this->position[2]; 
		res << x, y, z;
		return res; 
	}
	vector<float> orientation;
	float n; 
	float f; 
	float l;
	float r;
	float t;
	float b; 
};
/* Contains color vector for pixel on the screen. */
struct Pixel {
	 int x;
	 int y; 
	vector<int> pixColor; 
};

class Model
{
	/* Represents object given by .obj file */
public:
	// Material lighting properties
	Vector3f ambient;
	Vector3f diffuse; 
	Vector3f specular; 
	float shininess; 

	// Vector of all original world-space positions of vertices as passed in,
	// before any transformation 
	std::vector<Vertex> vertices;
	// Vector of all normals 
	std::vector<SurfaceNormal> normals; 
	// Vector of faces where each face contains 3 vertices with corresponding 
	// surface normals
	std::vector<Face> faces;
	// Name of objfile to parse 
	std::string filename;  
	// Name of object 
	std::string objName;
	// NOW STORING CURRENT TRANSFORMATION OF VERTEX IN THE STRUCT ITSELF 


	// NOT IN USE AT THE MOMENT: 
	/* For each entry keepPoints[i], store a 1 to indicate that the point is 
	 * within range or a 0 to indicate a point is out of range of the 
	 * perspective cube. */
	std::vector<bool> keepPoints;

	Model(std::string f, std::string name){
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
		vector<Vertex> vert;
		vector<SurfaceNormal> surfNorm;
		vector<Face> face;
		std::string line;

		Vector3f vertdummyVector; 
		// Create a dummy vertex at origin and push it 
		Vertex dummyV = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, vertdummyVector, 0, 0};
		SurfaceNormal dummyS = {0.0, 0.0, 0.0}; 

		vert.push_back(dummyV);
		surfNorm.push_back(dummyS);

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
				// Use linestream from now on
				// Vertex
				// cout << line << "\n";
				if(numWords[0] == "v") {
					float x, y, z;
					string v;
					linestream >> v >> x >> y >> z;
					Vector3f dummyVector;
					// cout << "x: " << x << " y: " << y << " z: " << z << "\n";
					Vertex newVert = {x, y, z, 1.0, x, y, z, 1.0, dummyVector, 0, 0};
					vert.push_back(newVert);
				}
				// SurfaceNormal
				else if(numWords[0] == "vn") {
					// Store for i, currI 
					float x, y, z; 
					string vn;
					linestream >> vn >> x >> y >> z;
					SurfaceNormal newSurf = {x, y, z};
					surfNorm.push_back(newSurf);
				}
				// Face with vertexIdx//normalIdx
				else if(numWords[0] == "f") {
					// cout << line << " : ";
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
					//cout << v1 << "\\" << s1 << " " << v2 << "\\" << s2 << 
					//" " << v3 << "\\" << s3 << "\n";
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
		// Saving all fields 
		vertices = vert;
		faces = face; 
		// cout << "Size of faces from makeScene.h: " << faces.size() << "\n";
		normals = surfNorm; 
		filename = f; 	
		// Closing file after operations
		infile.close();
	}
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
	// Storing perspective matrix for scene 
	Matrix4f cameraTransMatrix; 
	// Storing camera transformation for scene 
	Matrix4f perspectiveTransMatrix; 
};

// Declare all function headers for makeScene.cpp
/* Loads in the lights/ camera data from the file and stores in the lights vector/ camera struct
 * for the scene. */ 
void lightsCam(std::ifstream &currFile, Camera &cam, vector<Light> &lights);

/* Returns matrix to transform points to camera space */
Matrix4f cameraSpaceMatrix(Camera &currCam);

/* Returns matrix to transform points to perspective space. */ 
Matrix4f perspectiveMatrix(Camera &currCam);

/* Receives reference to scene object and applies 
 * world-space -> camera space transformation to objects.
 * Updates the currPosiions value for each model. */
void cameraSpaceTransformation(Scene &sceneContents);

/* Applies camera space -> ndc space. Basically exact same as 
 * cameraSpaceTransformation. */
void perspectiveTransformation(Scene &sceneContents);

/* ndc -> screen. */
void screenTransformation(Scene &sceneContents, int xres, int yres);

/* Sequentially applies transforms to go from transformed world space -> perspecive space. */ 
void applyAllTransforms(Scene &sceneContents, int xres, int yres);

/* Applies specified transformation to all vertices in model, and adds model to list of scene models. */
void modelTransformation(Matrix4f transformMatrix, std::string filename, int copyNum, std::string objName, \
 vector<Model> &allObjects);

 /* Parse txt file containing camera, object filenames, and object transformations.
 * Create a scene object that contains the camera, and a vector of objects. 
 * Each fileInfo object should contain a list of vertices and faces. */ 
Scene setScene(std::ifstream &currFile, int xres, int yres, int shadingMode); 