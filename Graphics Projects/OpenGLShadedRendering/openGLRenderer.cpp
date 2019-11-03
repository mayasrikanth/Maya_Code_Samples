#include "makeSceneOpenGL.h" // Access all structures and classes for makeScene.cpp
#include "Quaternion.h" // Access all current quaternion classes 
// For accessing main functions, data structures, and variables 
// to allow OpenGL development. 
#include <GL/glew.h>
#include <GL/glut.h>

#include "math.h"
#define _USE_MATH_DEFINES // For accessing double pi, need to convert to float


#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include<fstream>
#include<cmath>

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
#include<map>   // For dictionary 

/* All function headers. */
void GLTransform(vector<string> transformations);
void init(void);
void reshape(int width, int height);
void display(void);
void init_lights();
void set_lights();
void draw_objects();
void mouse_pressed(int button, int state, int x, int y);
void mouse_moved(int x, int y);
void key_pressed(unsigned char key, int x, int y);
float deg2rad(float angle);
Matrix4d Get_Current_Rotation(); 
/* Declaring global scene object that will be populated by main.*/
Scene currScene; 

/* Global parameters for creating an interactive first-person 
 * camera view of the scene. */
int mouse_x, mouse_y;

bool drag = false; 
bool is_pressed = false;
bool wireframe_mode = false;

/* Declaring global variables to represent current and last (Quaternion) rotations. */
quaternion current_rotation; 
quaternion last_rotation; 
// Assuming nothing different needs to be done to account for 
// transformation of normals 

/* Specify OpenGL states */
void init(void){
	/* Initializing current_rotation and last_rotation to Identity. */
	current_rotation = IdentityQuaternion(); 
	last_rotation = IdentityQuaternion(); 
	/* Smooth shading. */
	glShadeModel(GL_SMOOTH);
	/* Back-face culling. */
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	/* Use depth-buffering. */
	glEnable(GL_DEPTH_TEST);
	/* Automatically normalize vectors before passing them into
	 * the normals array. */
	glEnable(GL_NORMALIZE);
	/* Enable vertex array and norrmal array functionality. */
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	/* Tell openGL that we will be modifying the Projection Matrix. */
	// Projection matrix is applied to points in camera space -> NDC
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); // setting Projection Matrix to identity
	/* Create perspective projection matrix using frustum parameters.
	 * Current main matrix (identity) is post-multiplied
	 * with the newly-created matrix to create the new projection 
	 * matrix. */
	glFrustum(currScene.cam.l, currScene.cam.r, 
			  currScene.cam.b, currScene.cam.t, 
			  currScene.cam.n, currScene.cam.f);
	/* Note that the ModelView Matrix is the matrix that OpenGL
	 * applies to untransformed points in world space. 
	 * Here, it will factor in the geometric 
	 * and camera space transformations.*/
	// From now-on, any commands called will modify the ModelView matrix
	glMatrixMode(GL_MODELVIEW);

	/* Initialize lights to represent our Point Light structures. */
	init_lights();
}

void reshape(int width, int height) {
	// Prevent width or height from ever becoming 0, 1x1 is smallest
	height = (height == 0) ? 1 : height;
	width = (width == 0) ? 1 : width; 
	// Specify how to convert from NDC -> screen given new dimensions, 
	// setting lower-left corner to be (0,0)
	glViewport(0, 0, width, height);

	// Tell OpenGL tha tour program window needs to be re-displayed
	// i.e. everything that was being displayed before resizing needs
	// to be re-rendered
	glutPostRedisplay();
}
Matrix4d Get_Current_Rotation() {
	// Compute current rotation 
	//cout << "Current rotation:  " << current_rotation.s << " \n" << current_rotation.imaginary 
	//<< "\n Last Rotation:  " << last_rotation.s << " \n" << last_rotation.imaginary << "\n"; s
	quaternion prod = quaternionProduct(current_rotation, last_rotation); 
	// cout << "Product: " << prod.s << " \n" << prod.imaginary << "\n";
	// Normalize this product 
	quaternion norm_prod = normalizeQuaternion(prod);
	// Retrieve rotation matrix
	Matrix4d currRot = quaternionRotationMatrix(norm_prod); 
	// cout << "Rotation Matrix: \n" << currRot << "\n";
	return currRot; 
}

/* Handle all processing of points in world and camera space for display. */
void display(void) {
	// Reset color buffer (black-out everything) and depth buffer 
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// glMatrixMode(GL_MODELVIEW);
	// Assuming we are currently in ModelView Matrix mode
	glLoadIdentity(); // Initializing ModelView to Identity matrix
	// Due to openGL's post-multiplication, we need to specify 
	// all transformations in reverse order. 
	// Specify inverse rotation of camera by its orientation angle and
	// about its orientation axis
	Camera &sceneCam = currScene.cam; // Grabbing scene camera
	vector<float> &camOrientation = sceneCam.orientation;
	// Arguments are (-theta, x, y, z) (normalized during parsing)
	float theta = camOrientation[3] * (180.0 / M_PI);

	glRotatef(-theta, camOrientation[0], 
			   camOrientation[1], camOrientation[2]);
	vector<float> &camPos = sceneCam.position;
	glTranslatef(-camPos[0], -camPos[1], -camPos[2]);

	// Current rotation 
	// quaternion norm_rotation = normalizeQuaternion(current_rotation); 
	// Get current rotation matrix 
	Matrix4d currRot = Get_Current_Rotation();
	//Matrix4d currRot = quaternionRotationMatrix(norm_rotation);
	// Do we need to transform a specific matrix ? (MatrixMode is ModelView)
	glMultMatrixd(currRot.data()); 


	// Set up all lights in specified position
	set_lights();

	// cout << "Calling draw_objects from display() " << "\n";
	// Once lights set, specify points and faces to draw.
	draw_objects();

	// Double-buffering for smoother user experience 
	glutSwapBuffers();
}

/* Enables lighting calculation during rendering proces. */
void init_lights() {
	// Enable lighting calculations during rendering process 
	// (Phong reflection model or lighting model)
	glEnable(GL_LIGHTING);
	int num_lights = currScene.sceneLights.size();
	vector<Light> &lights = currScene.sceneLights;
	for(int i = 0; i < num_lights; i ++) {
		// Associate each scene light with one of OpenGL's built in lights
		int light_id = GL_LIGHT0 + i;
		glEnable(light_id);
		// Set color of all components of light (assuming color 
		// component is an array)
		glLightfv(light_id, GL_AMBIENT, lights[i].color);
		glLightfv(light_id, GL_DIFFUSE, lights[i].color);
		glLightfv(light_id, GL_SPECULAR, lights[i].color);
		// Setting attenuation 
		glLightf(light_id, GL_QUADRATIC_ATTENUATION, 
			lights[i].k); 
	}
}

/* Position the lights (uses ModelView Matrix to position these lights
 * correctly in camera space. Assuming ModelView is up to date with 
 * correct geometric and camera transformations. */
void set_lights() {
	int num_lights = currScene.sceneLights.size();
	vector<Light> &lights = currScene.sceneLights;
	for(int i = 0; i < num_lights; i ++) {
		int light_id = GL_LIGHT0 + i;
		glLightfv(light_id, GL_POSITION, lights[i].position);
	}
}
/* Function to parse transformation line to apply to the ModelViewMatrix.
 * The vector passed into this function should contain transformations in 
 * the reverse order of their parsing order in SetScene's InitModel.  */
void GLTransform(vector<string> transformations) {
	// Create separate rotate, translate, transform procedure for calling 
	// opengl 
	for(int i = 0; i < transformations.size(); i ++) {
		string c; 
		stringstream temp(transformations[i]);
		float x,y,z,theta; 
		temp >> c; 
		// Translation
		if(c == "t") {
			temp >> x >> y >> z; 
			// cout << "Translate: " << x << " " << y << " " << z << "\n";
			glTranslatef(x, y, z);
		}
		// Rotation 
		else if(c == "r") {
			temp >> x >> y >> z >> theta; 
			theta = theta * (180.0 / M_PI);
			float mag = sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0));
			x = x/mag;
			y = y/mag;
			z = z/mag;
			// cout << "Rotate: " << x << " " << y << " " << z << " " << theta << "\n";
			glRotatef(theta, x, y, z);
		}
		// scaling 
		else if(c == "s") {
			temp >> x >> y >> z;
			// cout << "Scale: " << x << " " << y << " " << z << "\n";
			glScalef(x, y, z);
		}

	}
}

/* Tells OpenGL to render objects to the display screen. */
void draw_objects() {
	for(int i = 0; i < currScene.sceneModels.size(); i++) {
		// Push another copy of the current ModelViewMatrix onto stack
		glPushMatrix();
		// Access current model 
		Model &currModel = currScene.sceneModels[i];
		// Apply geometric transformations 
		GLTransform(currModel.modelTransformations);
		// Specify material properties to render 
		glMaterialfv(GL_FRONT, GL_AMBIENT, currModel.ambient_reflect);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, currModel.diffuse_reflect);
		glMaterialfv(GL_FRONT, GL_SPECULAR, currModel.specular_reflect);
		glMaterialf(GL_FRONT, GL_SHININESS, currModel.shininess);
		// Specify how to render geometry
		glVertexPointer(3, GL_FLOAT, 0, &currModel.vertex_buffer[0]);
		glNormalPointer(GL_FLOAT, 0, &currModel.normal_buffer[0]);
		// Render the object
		int buffer_size = currModel.vertex_buffer.size();
		if(!wireframe_mode) {
			// cout << "About to draw arrays" << "\n";
			glDrawArrays(GL_TRIANGLES, 0, buffer_size);
		}
		else {
			// Render lines instead of triangle surfaces 
			for(int j = 0; j < buffer_size; j += 3) {
				glDrawArrays(GL_LINE_LOOP, j, 3);
			}
		}
 		// Get back original Modelview Matrix befoer proceeding to 
 		// render the next model 
 		glPopMatrix();
 	}
}

/* Respond to mouse clicks and releases. */
void mouse_pressed(int button, int state, int x, int y) {
	// If the left-mouse button was clicked down, then...
    if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
    {
        /* Store the mouse position in our global variables.
         */
        mouse_x = x;
        mouse_y = y;
        
        /* Since the mouse is being pressed down, we set our 'is_pressed"
         * boolean indicator to true.
         */
        is_pressed = true;
    }
    /* If the left-mouse button was released up, then...
     */
    else if(button == GLUT_LEFT_BUTTON && state == GLUT_UP)
    {
        /* Mouse is no longer being pressed, so set our indicator to false.
         */
    	// Update current_rotation (quaternion)
    	quaternion product = quaternionProduct(current_rotation, last_rotation);
    	// Cache all previous rotations in last_rotation
    	last_rotation = product; 
    	// Reset current rotation to identity
    	current_rotation = IdentityQuaternion();
    	// glutPostRedisplay();
        is_pressed = false;
    }
}
/* Respond to click-drag. (i.e. mouse is being moved. */
void mouse_moved(int x, int y) {
	if(is_pressed) {
	// drag = true; 
    // Compute Quaternion Rotation matrix 
    // Grab current x and y dimensions 
    GLfloat V[4];
    glGetFloatv(GL_VIEWPORT, V);
    int xres = V[2];
    int yres = V[3];
    // Update global current_rotation quaternion
    current_rotation = ComputeRotationQuaternion(mouse_x, mouse_y, x, y, xres, yres);
    // Update mouse positions to current position
    // cout << "Current Rotation:  " << current_rotation.s << " \n" << current_rotation.imaginary << "\n";
    /* Tell OpenGL that it needs to re-render our scene with the new camera
     * angles.
     */
    glutPostRedisplay();

	}
	
}
// Degrees -> Radians
float deg2rad(float angle) {
	return angle * M_PI / 180.0;
}
/* Respond to key presses on keyboard. */
void key_pressed(unsigned char key, int x, int y)
{
    /* If 'q' is pressed, quit the program.
     */
    if(key == 'q')
    {
        exit(0);
    }
    /* If 't' is pressed, toggle our 'wireframe_mode' boolean to make OpenGL
     * render our cubes as surfaces of wireframes.
     */
    else if(key == 't')
    {
        wireframe_mode = !wireframe_mode;
        /* Tell OpenGL that it needs to re-render our scene with the cubes
         * now as wireframes (or surfaces if they were wireframes before).
         */
        glutPostRedisplay();
    }
    else
    {
        /* Use current change in the horizontal camera angle (ie. the
         * value of 'x_view_angle') to compute the correct changes in our x and
         * z coordinates in camera space as we move forward, backward, to the left,
         * or to the right.
         *
         * 'step_size' is an arbitrary value to determine how "big" our steps
         * are.
         *
         * We make the x and z coordinate changes to the camera position, since
         * moving forward, backward, etc is basically just shifting our view */
         

    }
}

/* Input in format: ./openGLRenderer filename.txt xres yres */
int main(int argc, char* argv[]) {

	int xres = atoi(argv[2]);
	int yres = atoi(argv[3]);
	/* Populate the scene. */
	int shadingMode = 0; // dummy shading mode, not used for openGL
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
	// Find current working directory
	char cwd[256];
	if(getcwd(cwd, sizeof(cwd)) == NULL) {
		cout << "error... with retrieving directory" << endl;
	}
	else {
		// cout<< "Current working directory is: " << string(cwd) << endl;
	}
	
	// Populate the global scene variable 
	currScene = setScene(infile, xres, yres, shadingMode);
	//init();

	 /* 'glutInit' intializes the GLUT (Graphics Library Utility Toolkit) library.
     *
     * 'glutInit' takes the 'main' function arguments as parameters. This is not
     * too important for us, but it is possible to give command line specifications
     * to 'glutInit' by putting them with the 'main' function arguments.
     */
    glutInit(&argc, argv);
    /* The following line of code tells OpenGL that we need a double buffer,
     * a RGB pixel buffer, and a depth buffer.
     */
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    /* The following line tells OpenGL to create a program window of size
     * 'xres' by 'yres'.
     */
    glutInitWindowSize(xres, yres);
    /* The following line tells OpenGL to set the program window in the top-left
     * corner of the computer screen (0, 0).
     */
    glutInitWindowPosition(0, 0);
    /* The following line tells OpenGL to name the program window 
     */
    glutCreateWindow("Test");
    
    /* Call our 'init' function...
     */
    init();
    /* Specify to OpenGL our display function.
     */
    glutDisplayFunc(display);
    /* Specify to OpenGL our reshape function.
     */
    glutReshapeFunc(reshape);
    /* Specify to OpenGL our function for handling mouse presses.
     */
    glutMouseFunc(mouse_pressed);
    /* Specify to OpenGL our function for handling mouse movement.
     */
    glutMotionFunc(mouse_moved);
    /* Specify to OpenGL our function for handling key presses.
     */
    glutKeyboardFunc(key_pressed);
    /* The following line tells OpenGL to start the "event processing loop". This
     * is an infinite loop where OpenGL will continuously use our display, reshape,
     * mouse, and keyboard functions to essentially run our program.
     */
    glutMainLoop();


}