#include "makeScene.h" // Need this to access all the structures and classes for makeScene.cpp
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
// For maximum float size (depth buffer)
#include <limits>

using std::vector;
using namespace std;
using namespace Eigen;
#include <string>
#include <sstream>
#include<map>

/* For testing purposes: */
void drawLine(Vertex &v0, Vertex &v1, int** grid, int xres, int yres);
/* A function to convert transform all the vertices in the scene to NDC
 * coordinates. */
void worldToNDC(Scene &currScene);

 /* Procedure to determine color field for each point. */
void Lighting(Scene &currScene, SurfaceNormal &currNorm, Vertex &currVertex, Vector3f &diffuse, Vector3f &ambient, 
	Vector3f &specular, float p, Vector3f &cameraCoords);
/* Preprocesses the vertices (i.e. sets their color value) before calling GouraudRasterize. */
void GouraudShading(Model &currModel, Face &currFace, Scene &currScene, float** depthbuffer, vector<vector<Pixel>> &screenColors);
 /* Performs rasterization, depth-buffering, and backface culling using Gouraud 
 * technique. */
void GouraudRasterize(Scene &currScene, Model &currModel, Face &currFace, float** depthbuffer, 
	vector<vector<Pixel>> &screenColors);
/* Return alpha, beta, gamma for the triangle. */
float computeAlpha(int x_a, int y_a, int x_b, int y_b, int x_c, int y_c, int x_p, int y_p);
float computeBeta(int x_a, int y_a, int x_b, int y_b, int x_c, int y_c, int x_p, int y_p);
float computeGamma(int x_a, int y_a, int x_b, int y_b, int x_c, int y_c, int x_p, int y_p);
// Returns true if the NDC coordinate is within the 2x2x2 perspective cube, false if 
// it isn't. 
bool checkNDCRange(Vector4f NDC); 
// in-place conversion of a single vertex screen coordinates to use in the Phong algorithm 
void NDCToScreenVertex(Vertex &currVertex, int xres, int yres);
/* Performs phong shading and rasterization. */
void PhongShading(Model &currModel, Face &currFace, Scene &currScene, float** depthbuffer, vector<vector<Pixel>> &screenColors);
/* Outputs screenColors grid to ppm file. */
void createPPM(vector<vector<Pixel>> &screenColors, int xres, int yres);

float computeAlpha(int x_a, int y_a, int x_b, int y_b, int x_c, int y_c, int x_p, int y_p) {
	// screen coordinates -> barycentric to determine alpha 
	float numerator = (y_b - y_c) * x_p + (x_c - x_b) * y_p + x_b * y_c - x_c * y_b;
	if(numerator != 0) {
		// cout << numerator << "\n";
	}
	float denominator = (y_b - y_c) * x_a + (x_c - x_b) * y_a + x_b * y_c - x_c * y_b;
	float ans = numerator / denominator; 
	return ans; 
}
float computeBeta(int x_a, int y_a, int x_b, int y_b, int x_c, int y_c, int x_p, int y_p) {
	float numerator = (y_a - y_c) * x_p + (x_c - x_a) * y_p + x_a * y_c - x_c * y_a;
	float denominator = (y_a - y_c) * x_b + (x_c - x_a) * y_b + x_a * y_c - x_c * y_a;
	float ans = numerator / denominator;
	return ans;
}
float computeGamma(int x_a, int y_a, int x_b, int y_b, int x_c, int y_c, int x_p, int y_p) {
	float numerator = (y_a - y_b) * x_p + (x_b - x_a) * y_p + x_a * y_b - x_b * y_a;
	float denominator = (y_a - y_b) * x_c + (x_b - x_a) * y_c + x_a * y_b - x_b * y_a;
	float ans = numerator / denominator; 
	return ans; 
}
// Returns true if the NDC coordinate is within the 2x2x2 perspective cube, false if 
// it isn't. 
bool checkNDCRange(Vector4f NDC) {
	float x = NDC(0);
	float y  = NDC(1);
	float z = NDC(2);
	if(x >= -1 and x <= 1 and y >= -1 and y <= 1 and z >= -1 and z <= 1) {
		return true;
	}
	else {
		return false; 
	}
}

/* Procedure to apply to one face at a time: Preprocesses the vertices (i.e. sets their color value) before 
 * calling GouraudRasterize. */
void GouraudShading(Model &currModel, Face &currFace, Scene &currScene, float** depthbuffer, vector<vector<Pixel>> &screenColors)
 {
 	// Obtain vectors 
 	Vector3f diffuse = currModel.diffuse;
 	Vector3f ambient = currModel.ambient; 
 	Vector3f specular = currModel.specular; 
 	float p = currModel.shininess; 
 	// Obtain the camera position in eigen vector form 
 	Vector3f cameraCoords = currScene.cam.vectorizePos(); 
	// Retrieve indices for vertices and normals 
	int v1 = currFace.v1;
	int v2 = currFace.v2;
	int v3 = currFace.v3;
	int s1 = currFace.s1;
	int s2 = currFace.s2;
	int s3 = currFace.s3;
	
	// Grab the vertices and normals from the currModel
	Vertex &vert_a = currModel.vertices[v1];
	Vertex &vert_b = currModel.vertices[v2];
	Vertex &vert_c = currModel.vertices[v3];
	// Extract surface normals
	SurfaceNormal &a_norm = currModel.normals[s1];
	SurfaceNormal &b_norm = currModel.normals[s2];
	SurfaceNormal &c_norm = currModel.normals[s3]; 
	// Call lighting model on each pair to perform in-place assignment of the color
	// fields for each vertex. 
	Lighting(currScene, a_norm, vert_a, diffuse, ambient, specular, p, cameraCoords);
	Lighting(currScene, b_norm, vert_b, diffuse, ambient, specular, p, cameraCoords);
	Lighting(currScene, c_norm, vert_c, diffuse, ambient, specular, p, cameraCoords);

	// void Lighting(Scene &currScene, SurfaceNormal &currNorm, Vertex &currVertex, Vector3f &diffuse, Vector3f &ambient, 
	// Vector3f &specular, float p, Vector3f &cameraCoords)
	GouraudRasterize(currScene, currModel, currFace, depthbuffer, screenColors); 

	// Now call Gouraud Rasterize on this face 
	// Call the lighting model 
	// void Lighting(Scene &currScene, SurfaceNormal &currNorm, Vertex &currVertex, Vector3f &diffuse, Vector3f &ambient, 
	// Vector3f &specular, float p, Vector3f &cameraCoords)
}

/* Rasterizes a single face. Backface culling and depth-buffering as well.*/
void GouraudRasterize(Scene &currScene, Model &currModel, Face &currFace, float** depthbuffer, 
	vector<vector<Pixel>> &screenColors) {
	// dimensions 
	int xres = currScene.xres;
	int yres = currScene.yres;
	// CURRENTLY IN NDC SPACE:
	// Retrieve indices for vertices and normals 
	int v1 = currFace.v1;
	int v2 = currFace.v2;
	int v3 = currFace.v3;
	if(v1 != v2 or v2 != v3 or v1 != v3) {
		// cout << "v1: " << v1 << " v2: " << v2 << " v3: " << v3 << "\n"; 
	}
	// Grab the vertices and normals from the currModel (the color fields 
	// are already assigned)
	Vertex &a = currModel.vertices[v1];
	Vertex &b = currModel.vertices[v2];
	Vertex &c = currModel.vertices[v3];
	// Remember that the wordl-> NDC transformation has been applied by the time 
	// GouraudRasterize is called 
	Vector4f NDC_a = a.Vectorize();
	Vector3f NDC_a3;
	NDC_a3 << NDC_a(0), NDC_a(1), NDC_a(2); 
	Vector4f NDC_b = b.Vectorize();
	Vector3f NDC_b3; 
	NDC_b3 << NDC_b(0), NDC_b(1), NDC_b(2);
	Vector4f NDC_c = c.Vectorize(); 
	Vector3f NDC_c3; 
	NDC_c3 << NDC_c(0), NDC_c(1), NDC_c(2);

	// Backface Culling: check cross product of NDC vectors (should be assigned already
	// in the vertex struct)
	// In-place NDC->Screen transformation (we stillhave access to geometric coords)

	NDCToScreenVertex(a, xres, yres);
	NDCToScreenVertex(b, xres, yres);
	NDCToScreenVertex(c, xres, yres); 


	// Check cross-product 
	Vector3f cb = NDC_c3 - NDC_b3; 
	Vector3f ab = NDC_a3 - NDC_b3; 
	Vector3f crossproduct = cb.cross(ab); 
	// check z-component < 0 
	if(crossproduct(2) < 0) {
		return;
	}
	int x_a, y_a, x_b, y_b, x_c, y_c; 
	x_a = a.xscreen;
	y_a = a.yscreen;
	x_b = b.xscreen;
	y_b = b.yscreen;
	x_c = c.xscreen;
	y_c = c.yscreen; 
	// Determine bounding box 
	int x_min = min({x_a, x_b, x_c});
	int x_max = max({x_a, x_b, x_c});
	int y_min = min({y_a, y_b, y_c});
	int y_max = max({y_a, y_b, y_c});

	// Obtain color vectors for each vertex 
	Vector3f color_a = a.color;
	Vector3f color_b = b.color;
	Vector3f color_c = c.color;
	// Now extract the individual R,G,B values for these
	float R_a, R_b, R_c, G_a, G_b, G_c, B_a, B_b, B_c;
	R_a = color_a[0];
	G_a = color_a[1];
	B_a = color_a[2];
	R_b = color_b[0];
	G_b = color_b[1];
	B_b = color_b[2];
	R_c = color_c[0];
	G_c = color_c[1];
	B_c = color_c[2]; 
	// Depth-buffering: Update buffer right before you fill the pixel grid 
	for(int x = x_min; x <= x_max; x++) {
		for(int y = y_min; y <= y_max; y++) {
			float alpha = computeAlpha(x_a, y_a, x_b, y_b, x_c, y_c, x, y);
			float beta = computeBeta(x_a, y_a, x_b, y_b, x_c, y_c, x, y);
			float gamma = computeGamma(x_a, y_a, x_b, y_b, x_c, y_c, x, y);
			Vector4f currNDC; 
			if(alpha >= 0 and alpha <= 1 and beta >= 0 and beta <= 1 and gamma >= 0 and gamma <= 1) {
				currNDC = (alpha * NDC_a) + (beta * NDC_b) + (gamma * NDC_c);
				// Remember: depthbuffer and pixelgrid are indexed by (y_value, x_value)
				if(checkNDCRange(currNDC) == true and !(currNDC(2) > depthbuffer[y][x])) {
					// Saving the z-component of the NDC vector in the depthbuffer 
					depthbuffer[y][x] = currNDC(2);
					// compute the color 
					float R = alpha * R_a + beta * R_b + gamma * R_c;
					float G = alpha * G_a + beta * G_b + gamma * G_c;
					float B = alpha * B_a + beta * B_b + gamma * B_c;
					// Save the color for this screen point, scaling up to max value 
					R = R * 255;
					G = G * 255;
					B = B * 255;
					if(R != 0 or G != 0 or B != 0) {
						// cout << R << " " << G << " " << B << " \n";
					}
					//Pixel &currPix = screenColors[y][x];
					Pixel currPix; 
					vector<int> currColor {int(R), int(G), int(B)};
					// vector<int> currColor {255, 0, 255};
					currPix.pixColor = currColor; 
					screenColors[y][x] = currPix;
				}
			}
		}
	}
}

/* Rasterizes a single face. Backface culling and depth-buffering as well. */
void PhongShading(Model &currModel, Face &currFace, Scene &currScene, float** depthbuffer, 
	vector<vector<Pixel>> &screenColors) {
	int xres = currScene.xres;
	int yres = currScene.yres;
	// Retrieve indices of vertices and normals 
	int v1, v2, v3, s1, s2, s3;
	v1 =  currFace.v1;
	v2 = currFace.v2;
	v3 = currFace.v3;
	s1 = currFace.s1;
	s2 = currFace.s2;
	s3 = currFace.s3;
	// Retrieve vertices 
	Vertex &a = currModel.vertices[v1];
	Vertex &b = currModel.vertices[v2];
	Vertex &c = currModel.vertices[v3]; 
	// Retrieve world coordinates (the geometric transformation applied 
	// in MakeScene.cpp's initModel
	Vector4f v_a = a.worldVec();
	Vector4f v_b = b.worldVec();
	Vector4f v_c = c.worldVec(); 
	// Retrieve normals 
	SurfaceNormal &norma = currModel.normals[s1];
	SurfaceNormal &normb = currModel.normals[s2];
	SurfaceNormal &normc = currModel.normals[s3];
	// Retrieve normal Vectors 
	Vector3f n_a = norma.Vectorize();
	Vector3f n_b = normb.Vectorize();
	Vector3f n_c = normc.Vectorize(); 
	// Recall that the x,y,z,w fields of the vertices still have the
	// world-space, geometrically transformed coordinates 
	NDCToScreenVertex(a, xres, yres);
	NDCToScreenVertex(b, xres, yres);
	NDCToScreenVertex(c, xres, yres);
	// Retrieve NDC vectors 
	Vector4f NDC_a = a.Vectorize();
	Vector3f NDC_a3; 
	NDC_a3 << NDC_a(0), NDC_a(1), NDC_a(2);
	Vector4f NDC_b = b.Vectorize();
	Vector3f NDC_b3; 
	NDC_b3 << NDC_b(0), NDC_b(1), NDC_b(2);
	Vector4f NDC_c = c.Vectorize(); 
	Vector3f NDC_c3; 
	NDC_c3 << NDC_c(0), NDC_c(1), NDC_c(2);
	// v_2 - v_1 
	// v_0 - v_1
	// Back-face culling
	Vector3f cb = NDC_c3 - NDC_b3; 
	Vector3f ab = NDC_a3 - NDC_b3; 
	Vector3f cross = cb.cross(ab); 
	// check z-component < 0 
	if(cross(2) < 0) {
		return;
	}

	// Set screen coordinates for the vertices 
	NDCToScreenVertex(a, xres, yres);
	NDCToScreenVertex(b, xres, yres);
	NDCToScreenVertex(c, xres, yres);
	// Obtain the pixel coordinates 
	int x_a, y_a, x_b, y_b, x_c, y_c;
	x_a = a.xscreen;
	y_a = a.yscreen;
	x_b = b.xscreen;
	y_b = b.yscreen;
	x_c = c.xscreen;
	y_c = c.yscreen;
	// Determine Bounding Box 
	int x_min = min({x_a, x_b, x_c});
	int x_max = max({x_a, x_b, x_c});
	int y_min = min({y_a, y_b, y_c});
	int y_max = max({y_a, y_b, y_c});
	// Obtain the scene and camera information 
	Vector3f &diffuse = currModel.diffuse;
	Vector3f &ambient = currModel.ambient;
	Vector3f &specular = currModel.specular;
	float p = currModel.shininess;
	Vector3f cameraCoords = currScene.cam.vectorizePos(); 

	// Rasterize 
	for(int x = x_min; x <= x_max; x++) {
		for(int y = y_min; y <= y_max; y++) {
			float alpha = computeAlpha(x_a, y_a, x_b, y_b, x_c, y_c, x, y);
			float beta = computeBeta(x_a, y_a, x_b, y_b, x_c, y_c, x, y);
			float gamma = computeGamma(x_a, y_a, x_b, y_b, x_c, y_c, x, y);
			if(alpha >= 0 and alpha <= 1 and beta >= 0 and beta <= 1 and gamma >= 0 and gamma <= 1) {
				Vector4f currNDC = (alpha * NDC_a) + (beta * NDC_b) + (gamma * NDC_c);
				// Remember: depthbuffer and pixelgrid are indexed by (y_value, x_value)
				if(checkNDCRange(currNDC) == true and (currNDC(2) <= depthbuffer[y][x])) {
					// Saving the z-component of the NDC vector in the depthbuffer 
					depthbuffer[y][x] = currNDC(2);
					// Compute interpolated surface normal 	
					Vector3f interpNorm = (alpha * n_a) + (beta * n_b) + (gamma * n_c);
					// Turn this sum into a SurfaceNormal Instance
					SurfaceNormal currNorm;
					currNorm.setNormal(interpNorm);
					// Compute interpolated world-coordinate for point 
					Vector4f interpVert = (alpha * v_a) + (beta * v_b) + (gamma * v_c);
					// Turn this sum into a Vertex Instance 
					Vertex currVert; 
					currVert.x = interpVert(0);
					currVert.y = interpVert(1);
					currVert.z = interpVert(2);
					currVert.w = interpVert(3); 
					// Call the lighting model given on this vertex-normal pair 
					Lighting(currScene, currNorm, currVert, diffuse, ambient, specular, p, cameraCoords); 
					// Lighting(currScene, a_norm, vert_a, diffuse, ambient, specular, p, cameraCoords);
					// Extract color field from this vertex 
					Vector3f pointColor = currVert.color; 
					// Scaling up by max value of 255
					pointColor *= 255;
					// Update the screenColors Pixel with this color 
					// Convert to regular vector 
					Pixel currPix;
					vector<int> currColor {int(pointColor(0)), int(pointColor(1)), int(pointColor(2))};
					currPix.pixColor = currColor; 
					screenColors[y][x] = currPix;
				}
			}
		}
	}
}

void worldToNDC(Scene &currScene) {	
	vector<Vertex> transVert; 
	Matrix4f &cameraMatrix = currScene.cameraTransMatrix;
	Matrix4f &perspectiveMatrix = currScene.perspectiveTransMatrix;
	for(int i = 0; i < currScene.sceneModels.size(); i ++) {
		Model &currModel = currScene.sceneModels[i];
		for(int j = 0; j < currModel.vertices.size(); j ++) {
			// in-place transformation world space -> camera space 
			Vertex &currVertex = currModel.vertices[j];
			currVertex.transformVertex(cameraMatrix);
			// in-place transformation camera space -> homogenous ndc
			currVertex.transformVertex(perspectiveMatrix);
			// in-place scaling down by homogenous component -> cartesian ndc  
			currVertex.divideByW();
			transVert.push_back(currVertex);

		}
		currModel.vertices = transVert; 
		vector<Vertex>().swap(transVert);
	}
}

// Assuming that the points are already in NDC space! Performs in-place 
// transformation from ndc to screen. This does not overwrite the 
// most recent coordinates, rather it saves the screen coordinates to 
// separate fields in the vertex. 
void NDCToScreen(Scene &currScene) {
	int xres = currScene.xres;
	int yres = currScene.yres;

	for(int i = 0; i < currScene.sceneModels.size(); i ++) {
		Model &currModel = currScene.sceneModels[i];
		for(int j = 0; j < currModel.vertices.size(); j ++) {
			Vertex &currVertex = currModel.vertices[j];
			float x_ndc = currVertex.currX;
			float y_ndc = currVertex.currY;
			float z_ndc = currVertex.currZ;
			// x_range goes from [-1, 1] -> [0, w - 1]
			int x_screen = ((xres - 1) / 2) * (x_ndc + 1);  // need to be integer
			// y_range goes from [-1, 1] -> [h - 1, 0]
			int y_screen = (-1 * (yres - 1) / 2) * (y_ndc - 1);  // needs to be integer 
			// Assign screen coordinates to the vertex 
			currVertex.xscreen = x_screen;
			currVertex.yscreen = y_screen;	
		}
	}

}

// in-place conversion of a single vertex screen coordinates to use in the Phong algorithm 
void NDCToScreenVertex(Vertex &currVertex, int xres, int yres) {
	float x_ndc = currVertex.currX;
	float y_ndc = currVertex.currY;
	float z_ndc = currVertex.currZ;
	// x_range goes from [-1, 1] -> [0, w - 1]
	int x_screen = ((xres - 1) / 2) * (x_ndc + 1);  // need to be integer
	// y_range goes from [-1, 1] -> [h - 1, 0]
	int y_screen = (-1 * (yres - 1) / 2) * (y_ndc - 1);  // needs to be integer 
	// Assign screen coordinates to the vertex 
	currVertex.xscreen = x_screen;
	currVertex.yscreen = y_screen;

}
/* Sets the color vector for each vertex in the model in worldspace. Operates on a single vertex
 * struct and a single corresponding normal struct. Note that NO OTHER transformation except 
 * worldspace geometric should be applied before calling this function. */
void Lighting(Scene &currScene, SurfaceNormal &currNorm, Vertex &currVertex, Vector3f &diffuse, Vector3f &ambient, 
	Vector3f &specular, float p, Vector3f &cameraCoords) {
	// Extract coordinates of surface unit normal  and place in eigenvector 
	Vector3f normPos = currNorm.Vectorize();
	// Extract coordinates of current vertex (should be geometrically transformed ONLY)
	Vector3f vertPos;
	vertPos << currVertex.x, currVertex.y, currVertex.z; 
	// Create vectors to hold diffuse and specular fields
	Vector3f diffuseSum; 
	diffuseSum << 0, 0, 0;
	Vector3f specularSum; 
	specularSum << 0, 0, 0;
	// compute e_direction 
	Vector3f e_direction = cameraCoords - vertPos;
	// in-place normalization 
	e_direction.normalize(); 
	// vector of 1's for component-wise minimum later on
	Vector3f ones;
	ones << 1, 1, 1;
	// For each light 
	for(int i = 0; i < currScene.sceneLights.size(); i++) {
		Light &currLight = currScene.sceneLights[i]; 
		// Storing copy of current light's position
		Vector3f l_p = currLight.position; 
		Vector3f l_c = currLight.color;
		Vector3f l_direction = l_p - vertPos;
		// in-place normalization
		l_direction.normalize();
		// Computing Attenuation 
		float k = currLight.k;
		// (position of light - vertex position)
		Vector3f diff = l_p - vertPos;
		// Obtains frobenous norm, square root of the sum of squares of the components
		float euclidDist = diff.norm();
		float attenuationParam = 1 / (1 + k * euclidDist);
		// Multiplying l_c by this scalar 
		l_c *= attenuationParam;
		// Finding scalar coefficient to multiply l_c (color vector) by 
		float normDir = normPos.dot(l_direction);
		float zero = 0;
		float diffuseCoeff = max(zero, normDir);
		// Compute diffuse for this light
		Vector3f l_diffuse = l_c * diffuseCoeff;
		diffuseSum += l_diffuse;
		// Finding scalar coefficient to multiply l_c by 
		Vector3f norm_el = e_direction + l_direction;
		// in-place normalize 
		norm_el.normalize();
		// compute dot product with normal vector 
		float dotProd = normPos.dot(norm_el);
		float zero1 = 0;
		float maxval = max(zero1, dotProd);
		float specCoeff = pow(maxval, p);

		Vector3f l_specular = l_c * specCoeff;
		specularSum += l_specular; 
	}
	// Now compute color for the point
	// Component-wise product of diffuseSum and diffuse property of material 
	Vector3f product1 = diffuseSum.array() * diffuse.array();
	product1 += ambient;
	// Component-wise product of specularSum and specular property of material
	Vector3f product2 = specularSum.array() * specular.array();
	Vector3f term2 = product1 + product2;
	// Component-wise minimum 
	Vector3f c = term2.cwiseMin(ones);
	if(c.isZero(0) == 0) { 
		// cout << c(0) << " " << c(1) << " " << c(2) << "\n";
	}
	// Assign to vertex color field 
	currVertex.color = c; 
}


// Call the appropriate shading algorithm to produce a shaded, rasterized image containing 
// all models in the scene.  
void allAlgos(Scene &currScene) {
	// Initialize the pixel grid and depth buffer, as these are 
	// common to both shading algorithms. 
	int xres = currScene.xres;
	int yres = currScene.yres; 
	// Depth buffer: holds min distance of each NDC point from the camera. 
	float** depthbuffer = new float*[yres];
	// Init vector of type pixel to hold colors for each screen coordinate 
	// contains vectors of type pixels to represent the rows. 
	vector<Pixel> v(xres);
	vector<vector<Pixel>> screenColors(yres, v); 
	for(int i = 0; i < yres; i ++) {
		depthbuffer[i] = new float[xres];
		vector<Pixel> row; 
	}
	// Init everything to 0 (i.e. the background color)
	for(int r = 0; r < yres; r ++) {
		for(int c = 0; c < xres; c++) {
			depthbuffer[r][c] = std::numeric_limits<float>::max();
			// Should be integers from [0,255]
			vector<int> initColor {0, 0, 0};
			Pixel currPix;
			currPix.pixColor = initColor; 
			currPix.x = 0;
			currPix.y = 0; 
			screenColors[r][c] = currPix; 
		}
	} 

	// For the Gourard algorithm, use the preset color for the 
	// 3 vertices, and interpolate using barycentric coordinates. 
	if(currScene.shadingMode == 0) {
		// Gouraurd shading
		// Transform all points from world-NDC space (note that we still have access
		// to world-space geometrically transformed coordinates within the vertex 
		// structs)
		worldToNDC(currScene); 
		// Call GouraudShading for each face 
		//cout << "After NDC, Assignment 2: " << "\n";
		for(int i = 0; i < currScene.sceneModels.size(); i ++) {
			Model &currModel = currScene.sceneModels[i];
			for(int j = 0; j < currModel.faces.size(); j++) {
				Face &currFace = currModel.faces[j];
				GouraudShading(currModel, currFace, currScene, depthbuffer, screenColors); 
			}
		}
		  createPPM(screenColors, xres, yres);
	}


	else if(currScene.shadingMode == 1) {
		// Transform all points to NDC (note that we still have access to world-space
		// geometrically-transformed coordinates within the vertex struct)
		worldToNDC(currScene);
		// Model &currModel, Face &currFace, Scene &currScene, float** depthbuffer, vector<Pixel> &screenColors
		for(int i = 0; i < currScene.sceneModels.size(); i++) {
			Model &currModel = currScene.sceneModels[i];
			for(int j = 0; j < currModel.faces.size(); j ++) {
				Face &currFace = currModel.faces[j]; 
				// PhongShading, function signature below: 
				// (Model &currModel, Face &currFace, Scene &currScene, float** depthbuffer, vector<Pixel> &screenColors)
				PhongShading(currModel, currFace, currScene, depthbuffer, screenColors);

			}
		}
		createPPM(screenColors, xres, yres);

		/* Testing the vertices positions. */
		// Grid is a pointer to an array of pointers 
		/*int** grid =  new int*[yres];
		for(int i = 0; i < yres; i ++) {
			grid[i] = new int[xres]; 
		}
		// Init everything to 0 (i.e. the background color)
		for(int r = 0; r < yres; r ++) {
			for(int c = 0; c < xres; c++) {
				grid[r][c] = 0;
			}
		}
		worldToNDC(currScene);
		NDCToScreen(currScene); 
		for(int i = 0; i < currScene.sceneModels.size(); i ++ ){
			Model &currModel = currScene.sceneModels[i];
			for(int j = 0; j < currModel.faces.size(); j ++) {
				Face &currFace = currModel.faces[j];
				int v1 = currFace.v1;
				int v2 = currFace.v2; 
				int v3 = currFace.v3; 

				Vertex &vert1 = currModel.vertices[v1];
				Vertex &vert2 = currModel.vertices[v2];
				Vertex &vert3 = currModel.vertices[v3];
				drawLine(vert1, vert2, grid, xres, yres);
				drawLine(vert2, vert3, grid, xres, yres);
				drawLine(vert1, vert3, grid, xres, yres);
			}
		}
		// Now outputting vertices. 
		cout << "P3" << "\n"; // ppm format header 	
		cout << xres << " " << yres << "\n"; 
		cout << 255 << "\n";
		// white: (255, 255, 255)
		// black: (0,0,0)
		for (int r = 0; r < yres; r++) {
			for (int c = 0; c < xres; c++) {
				if(grid[r][c] == 0){
					// then black
					cout << 0 << " " << 0 << " " << 0 << "\n";
				}
				else if(grid[r][c] == 1) {
					// then white 
					cout << 255 << " " << 255 << " " << 255 << "\n";
				}
			}
		}*/

		// now deleting array of pointers 
		/*for(int i = 0; i < yres; i ++) {
			delete[] grid[i];
		}
		// deleting pointer to array of pointers 
		delete[] grid;*/
	}
		
		// FINALLY: FREE DYNAMICALLY ALLOCATED MEMORY
		// now deleting array of pointers 
		for(int i = 0; i < yres; i ++) {
			delete[] depthbuffer[i];
		}
		// deleting pointer to array of pointers 
		delete[] depthbuffer;
}

// USED FOR Debugging Purposes 
void drawLine(Vertex &v0, Vertex &v1, int** grid, int xres, int yres) {
	int x0 = v0.xscreen;
	int y0 = v0.yscreen;
	int x1 = v1.xscreen;
	int y1 = v1.yscreen;
	
	// Attempting to plot line from (x0, y0) -> (x1, y1)

	// First check if both points are in range. If not,
	// do not draw the line. 
	if(x0 < 0 or x0 > (xres - 1) or x1 < 0 or x1 > (xres - 1)) {

		return;
	} 
	if(y0 < 0 or y0 > (yres - 1) or y1 < 0 or y1 > (xres - 1)) {
		return; 
	}
	// Check if horizontal Line
	if(y1 == y0) {
		// Loop from x1 to x0
		if(x0 > x1) {
			for(int i = x1; i <= x0; i ++) {
				grid[y0][i] = 1;
			}
		}
		// Loop from x0 to x1
		else if(x1 > x0) {
			for(int i = x0; i <= x1; i ++) {
				grid[y0][i] = 1;
			}
		}
		return; 
	}

	// Check if vertical Line
	if(x1 == x0)  {
		//cout <<"Vertical line between " << "(" << x0 << ", " << y0 << ") and";
		//cout << "(" << x1 << "," << y1 << ")\n" ;
		// Loop from y1 to y0
		if(y0 > y1) {
			for(int i = y1; i <= y0; i ++) {
				grid[i][x0] = 1;
			}
		}
		// Loop from y0 to y1
		else if (y1 > y0) {
			for(int i = y0; i <= y1; i ++) {
				grid[i][x0] = 1; 
			}
		}
		return; 
	}
	// Attempting to draw line from (x0, y0) to (x1, y1)
	int x_s, x_e, y_s, y_e;
	// Calculating slope 
	// double m = (y - y) / (x_e - x_s); 
	int err, dy, dx, y_bresen, x_bresen; 
	// Check whether to iterate through y or x
	dy = abs(y0 - y1);
	dx = abs(x0 - x1); 

	if(dy >= dx) {
		// INCREMENTING THROUGH Y, optional x movement
		// enforce y direction is low -> high 
		if(y1 < y0) {
			x_s = x1;
			y_s = y1;
			x_e = x0;
			y_e = y0;
		}
		else if(y0 < y1) {
			x_s = x0;
			y_s = y0;
			x_e = x1;
			y_e = y1;
		}
		// Now determine x direction (right or left)? 
		if(x_s > x_e) {
			// move left
			err = 0;
			x_bresen = x_s;
			dx = -1 * (x_e - x_s);
			dy = y_e - y_s;
			for(int y = y_s; y <= y_e; y++){
				grid[y][x_bresen] = 1;
				if((2 * (err + dx)) < dy) {
					err = err + dx;
				}
				else {
					err = err + dx - dy;
					x_bresen = x_bresen - 1;
				}
			}
		}
		else if (x_e > x_s) {
			// move right 
			err = 0;
			x_bresen = x_s;
			dx = x_e - x_s;
			dy = y_e - y_s;
			for(int y = y_s; y <= y_e; y++){
				grid[y][x_bresen] = 1;
				if((2 * (err + dx)) < dy) {
					err = err + dx;
				}
				else {
					err = err + dx - dy;
					x_bresen = x_bresen + 1;
				}
			}
		}
	}

	else if (dx > dy) {
		// INCREMENTING THROUGH X, optional y movement 
		// enforce x direction is left -> right 
		if(x0 < x1) {
			x_s = x0;
			x_e = x1;
			y_s = y0;
			y_e = y1;
		}
		else if (x1 < x0) {
			x_s = x1;
			y_s = y1;
			x_e = x0;
			y_e = y0;
		}
		// Now determine y direction (up or down?) 
		if(y_e > y_s) {
			// move up
			err = 0;
			y_bresen = y_s;
			dx = x_e - x_s;
			dy = y_e - y_s;
			for(int x = x_s; x <= x_e; x++) {
				grid[y_bresen][x] = 1;
				if((2 *(err + dy)) < dx) {
					err = err + dy;
				}
				else {
					err = err + dy - dx;
					y_bresen = y_bresen + 1;
				}
			}
		}
		else if (y_s > y_e) {
			// move down 
			err = 0;
			y_bresen = y_s;
			dx = x_e - x_s;
			dy = -1 * (y_e - y_s);
			for(int x = x_s; x <= x_e; x++) {
				grid[y_bresen][x] = 1;
				if((2 *(err + dy)) < dx) {
					err = err + dy;
				}
				else {
					err = err + dy - dx;
					y_bresen = y_bresen - 1;
				}
			}
		}
	}
}

void createPPM(vector<vector<Pixel>> &screenColors, int xres, int yres) {
	cout << "P3" << "\n"; // ppm format header 	
	cout << xres << " " << yres << "\n"; 
	cout << 255 << "\n"; 
	// white: (255, 255, 255)
	// black: (0,0,0)
	for (int r = 0; r < yres; r++) {
		for (int c = 0; c < xres; c++) {
			Pixel currPix = screenColors[r][c];
			// should be a length-3 vector 
			vector<int> color = currPix.pixColor;
			cout << color[0] << " " << color[1] << " " << color[2] << "\n";
		}
	}
}
// Establish correct working directory and call setScene in makeScene.cpp. 
// Input is of form ./shaded_renderer scene_file.txt xres yres mode
int main(int argc, char** argv) {
	// Grab resolution 
	int xres = atoi(argv[2]);
	int yres = atoi(argv[3]);
	int shadingMode = atoi(argv[4]); 
	std::string filename; 
	// Determining correct file path 
	std::string path = argv[1];
	std::stringstream ss(path); 
	std::string tok; 
	vector<string> filepath; 
	while(std::getline(ss, tok, '/' )) {
		filepath.push_back(tok); 
	}
	// Then we need to extract the directory and use chdir
	if(filepath.size() > 1) {
		std::size_t found = path.find_last_of("/");
		// Now change current working directory 
		std::string workingDir = path.substr(0, found); 
		const char *p = workingDir.c_str(); 
		chdir(p); 
		// Set the filename to everything from \-> end (excluding last slash)
		filename = path.substr(found + 1);
	}
	// Then we are given only a file
	else if(filepath.size() == 1) {
		filename = path; 
	}
	std::ifstream infile(filename); 
	// cout<< "Filename: " << filename << "\n";
	// cout<< "Filepath: " << path << "\n";
	// what is current working directory?
	char cwd[256];
	if(getcwd(cwd, sizeof(cwd)) == NULL) {
		cout << "error... with retrieving directory" << endl;
	}
	else {
		// cout<< "Current working directory is: " << string(cwd) << endl;
	}
	// Grab the fully populated scene object (is there a more efficientt
	// way to do this?)
	Scene currScene = setScene(infile, xres, yres, shadingMode);
	allAlgos(currScene); 
	return 0;
}



