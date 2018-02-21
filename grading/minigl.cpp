/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>

// #include<iostream>	// for testing 

using namespace std;

/**
 * Useful data types
 */
typedef mat<MGLfloat,4,4> mat4; //data structure storing a 4x4 matrix, see mat.h
typedef mat<MGLfloat,3,3> mat3; //data structure storing a 3x3 matrix, see mat.h
typedef vec<MGLfloat,4> vec4;   //data structure storing a 4 dimensional vector, see vec.h
typedef vec<MGLfloat,3> vec3;   //data structure storing a 3 dimensional vector, see vec.h
typedef vec<MGLfloat,2> vec2;   //data structure storing a 2 dimensional vector, see vec.h

struct Vertex{
	vec4 position;
	vec3 vertex_color;
};
struct Triangle{
	vector<Vertex> triangle_verticies;
};

vector<Triangle> triangles;

vec3 color;
vector<Vertex> verticies;
MGLpoly_mode shape;
MGLmatrix_mode matrix_mode;
int orthoFlag = 0;
mat4 projection;

/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}

/** Helper function to calculate area of a triangle **/
float Triangle_Area(Vertex a, Vertex b, Vertex c){
	// if(1){
	// float length_AB, length_AC, length_BC;
	
	// // I'm using Heron's formula to calculate the area. 
	// // Find side lengths:
	
	// float newA = b.position[0] - a.position[0];	// may cause seg fault, look out
	// float newB = b.position[1] - a.position[1];
	// float newC = b.position[2] - a.position[2];
	
	// float d = sqrt(pow(newA, 2) + pow(newB, 2) + pow(newC, 2));	//distance formula
	// length_AB = d;
	
	// newA = c.position[0] - a.position[0];
	// newB = c.position[1] - a.position[1];
	// newC = c.position[2] - a.position[2];
	
	// d = sqrt(pow(newA, 2) + pow(newB, 2) + pow(newC, 2));
	// length_AC = d;
	
	// newA = c.position[0] - b.position[0];
	// newB = c.position[1] - b.position[1];
	// newC = c.position[2] - b.position[2];
	
	// d = sqrt(pow(newA, 2) + pow(newB, 2) + pow(newC, 2));
	// length_BC = d;
	
	// float s = (length_AB + length_AC + length_BC) / 2;

	// }	
	// return sqrt(s * (s - length_AB) * (s - length_AC) * (s - length_BC));


	float ax = a.position[0];
	float ay = a.position[1];
	float bx = b.position[0];
	float by = b.position[1];
	float cx = c.position[0];
	float cy = c.position[1];
	return ax*(by - cy) + ay * (cx - bx) + (bx * cy - by * cx);
	
}

float Triangle_Area(Triangle tri){
	return Triangle_Area(tri.triangle_verticies.at(0), tri.triangle_verticies.at(1), tri.triangle_verticies.at(2));
}


/** Helper function that gets a triangle and rasterizes the triangle on screen by setting the colors in data*/
void Rasterize_Triangle(const Triangle& tri, int width, int height, MGLpixel* data)
{
	// Calculate pixel coordinates
	float Ax, Ay, Bx, By, Cx, Cy;
	Vertex v1, v2, v3;
	
	Ax = (tri.triangle_verticies.at(0).position[0] + 1) * 0.5 * width;
	Ay = (tri.triangle_verticies.at(0).position[1] + 1) * 0.5 * height;
	Bx = (tri.triangle_verticies.at(1).position[0] + 1) * 0.5 * width;
	By = (tri.triangle_verticies.at(1).position[1] + 1) * 0.5 * height;
	Cx = (tri.triangle_verticies.at(2).position[0] + 1) * 0.5 * width;
	Cy = (tri.triangle_verticies.at(2).position[1] + 1) * 0.5 * height;
	
	v1.position = vec4(Ax, Ay, 0, 1);
	v2.position = vec4(Bx, By, 0, 1);
	v3.position = vec4(Cx, Cy, 0, 1);
	
	// For each pixel on screen,
	for(int i = 0; i < width; i++){
		for(int j = 0; j < height; j++){
			// Calculate barycentric coordinates
			Triangle t1, t2, t3;
			
			Vertex currVertex;
			currVertex.position = vec4(i, j, 0, 1);
			
			t1.triangle_verticies.push_back(v1);
			t1.triangle_verticies.push_back(v2);
			t1.triangle_verticies.push_back(currVertex);
			
			t2.triangle_verticies.push_back(v3);
			t2.triangle_verticies.push_back(v1);
			t2.triangle_verticies.push_back(currVertex);
			
			t3.triangle_verticies.push_back(v2);
			t3.triangle_verticies.push_back(v3);
			t3.triangle_verticies.push_back(currVertex);
			
			float alpha, beta, gamma;
			
			alpha = Triangle_Area(t1)/Triangle_Area(tri);
			beta = Triangle_Area(t2)/Triangle_Area(tri);
			gamma = Triangle_Area(t3)/Triangle_Area(tri);
			
			// Decide if pixel is inside the triangle (0 <= a/b/g)
			if(alpha >= 0 && beta >= 0 && gamma >= 0){
				data[i+j*width] = Make_Pixel(255,255,255);
			}
		}
	}

	
}


/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{
	// Fill whole pixel data with black initially
	for(unsigned i = 0; i < width; i++){
		for(unsigned j = 0; j < height; j++){
			data[i+j*width] = Make_Pixel(0,0,0);
		}
	}
	// For each triangle in the triangle list call Rasterize_Triangle
	for(unsigned i = 0; i < triangles.size(); i++){
		Rasterize_Triangle(triangles.at(i), width, height, data);
	}
	
	// Clear triangle list
	while(!triangles.empty()){
		triangles.pop_back();
	}
}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
	shape = mode;
}


/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
	if(shape == MGL_TRIANGLES){
		for(unsigned i = 0; i < verticies.size(); ){
			Vertex v1 = verticies.at(i);
			i++;
			Vertex v2 = verticies.at(i);
			i++;
			Vertex v3 = verticies.at(i);
			i++;
			
			Triangle t;
			t.triangle_verticies.push_back(v1);
			t.triangle_verticies.push_back(v2);
			t.triangle_verticies.push_back(v3);
			
			triangles.push_back(t);
		}
	}
	else if(shape == MGL_QUADS){
		for(unsigned i = 0; i < verticies.size(); ){
			Vertex v1 = verticies.at(i);
			i++;
			Vertex v2 = verticies.at(i);
			i++;
			Vertex v3 = verticies.at(i);
			i++;
			Vertex v4 = verticies.at(i);
			i++;
			
			Triangle t1, t2;
			t1.triangle_verticies.push_back(v1);
			t1.triangle_verticies.push_back(v2);
			t1.triangle_verticies.push_back(v3);
			
			t2.triangle_verticies.push_back(v1);
			t2.triangle_verticies.push_back(v3);
			t2.triangle_verticies.push_back(v4);
			
			triangles.push_back(t1);
			triangles.push_back(t2);
		}
	}
	
	while(!verticies.empty()){
		verticies.pop_back();
	}
}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
	mglVertex3(x,y,0);
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
	Vertex temp;
	vec4 pos = vec4(x,y,z,1);
	temp.vertex_color = color;
	temp.position = pos;
	
	if(orthoFlag){
		temp.position = projection * temp.position;
	}
	
	verticies.push_back(temp);
	
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
	matrix_mode = mode;
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
	
}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
	// https://www.scratchapixel.com/lessons/3d-basic-rendering/perspective-and-orthographic-projection-matrix/orthographic-projection-matrix
	// https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/glOrtho.xml
	if(matrix_mode == MGL_PROJECTION){
		orthoFlag = 1;
		MGLfloat tx, ty, tz;
		tx = -(right + left)/(right - left);
		ty = -(top + bottom)/(top - bottom);
		tz = -(far + near)/(far - near);
		// projection = {{(2/(right - left)), 0, 0, tx, 0, (2/(top - bottom)), 0, ty, 0, 0, (-2/(far - near)), tz, 0, 0, 0, 1}};
		projection = {{(2/(right - left)), 0, 0, 0, 0, (2/(top - bottom)), 0, 0, 0, 0, (-2/(far - near)), 0, tx, ty, tz, 1}};
	}
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
	color = vec3(red,green,blue);
}

