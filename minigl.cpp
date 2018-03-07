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
// int orthoFlag = 0;
// int frustumFlag = 0;
// mat4 projection = {{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}};	// initialize to identity
mat4 projection = {{4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}};
mat4 modelview = {{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}};	// initialize to identity

// zbuffer: https://www.cs.utexas.edu/~fussell/courses/cs384g/lectures/Lecture9-Zbuffer_pipeline.pdf
// SLIDE 16
vector<vector<float> > zbuffer;
int ZMAX = 1;
int ZMIN = -1;
vector<mat4> projection_stack{projection};
vector<mat4> modelview_stack{modelview};

/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}

/** Helper function to calculate area of a triangle **/
float Triangle_Area(Vertex a, Vertex b, Vertex c){
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
	
	Triangle tri2;
	Vertex tri2_v1, tri2_v2, tri2_v3;
	tri2_v1.position = vec4(Ax, Ay, 0, 1);
	tri2_v2.position = vec4(Bx, By, 0, 1);
	tri2_v3.position = vec4(Cx, Cy, 0, 1);
	tri2.triangle_verticies.push_back(tri2_v1);
	tri2.triangle_verticies.push_back(tri2_v2);
	tri2.triangle_verticies.push_back(tri2_v3);
	
	v1.position = vec4(Ax, Ay, 0, 1);
	v2.position = vec4(Bx, By, 0, 1);
	v3.position = vec4(Cx, Cy, 0, 1);
	
	Triangle t1, t2, t3;
	t1.triangle_verticies.push_back(v1);
	t1.triangle_verticies.push_back(v2);
	t2.triangle_verticies.push_back(v3);
	t2.triangle_verticies.push_back(v1);
	t3.triangle_verticies.push_back(v2);
	t3.triangle_verticies.push_back(v3);
			
	Vertex currVertex;
	
	// For each pixel on screen;
	// Use bounding rectangle, more efficient
	// http://www.cs.ucr.edu/~shinar/courses/cs130/content/07-TriangleRasterization.pdf
	// slide 19
	
	int x_max, x_min, y_max, y_min;
	if(Ax >= Bx && Ax >= Cx){
		x_max = Ax;
	}
	else if(Bx >= Ax && Bx >= Cx){
		x_max = Bx;
	}
	else{
		x_max = Cx;
	}
	
	if(Ax <= Bx && Ax <= Cx){
		x_min = Ax;
	}
	else if(Bx <= Ax && Bx <= Cx){
		x_min = Bx;
	}
	else{
		x_min = Cx;
	}
	
	if(Ay >= By && Ay >= Cy){
		y_max = Ay;
	}
	else if(By >= Ay && By >= Cy){
		y_max = By;
	}
	else{
		y_max = Cy;
	}
	
	if(Ay <= By && Ay <= Cy){
		y_min = Ay;
	}
	else if(By <= Ay && By <= Cy){
		y_min = By;
	}
	else{
		y_min = Cy;
	}
	
	if(y_min < 0){
		y_min = 0;
	}
	if(x_min < 0){
		x_min = 0;
	}
	if(x_max >= width){
		x_max = width;
	}
	if(y_max >= height){
		y_max = height;
	}
	
	for(int i = x_min; i < x_max; i++){
		for(int j = y_min; j < y_max; j++){
			// Calculate barycentric coordinates
			
			currVertex.position = vec4(i, j, 0, 1);
			t1.triangle_verticies.push_back(currVertex);
			t2.triangle_verticies.push_back(currVertex);
			t3.triangle_verticies.push_back(currVertex);
			
			float alpha = Triangle_Area(t1)/Triangle_Area(tri2);
			float beta = Triangle_Area(t2)/Triangle_Area(tri2);
			float gamma = Triangle_Area(t3)/Triangle_Area(tri2);
			
			t1.triangle_verticies.pop_back();
			t2.triangle_verticies.pop_back();
			t3.triangle_verticies.pop_back();
			// Decide if pixel is inside the triangle (0 <= a/b/g)
			if(alpha >= 0 && beta >= 0 && gamma >= 0){

				float currentZ = gamma * tri.triangle_verticies.at(0).position[2] + 
						         beta *  tri.triangle_verticies.at(1).position[2] +
						         alpha * tri.triangle_verticies.at(2).position[2];
				currentZ = -currentZ;
				if(currentZ >= zbuffer.at(i).at(j) && currentZ <= ZMAX){
					int r,g,b;
					
					r = gamma * tri.triangle_verticies.at(0).vertex_color[0]
					  + beta * tri.triangle_verticies.at(1).vertex_color[0]
					  + alpha * tri.triangle_verticies.at(2).vertex_color[0];
					  
					g = gamma * tri.triangle_verticies.at(0).vertex_color[1]
					  + beta * tri.triangle_verticies.at(1).vertex_color[1]
					  + alpha * tri.triangle_verticies.at(2).vertex_color[1];
					  
					b = gamma * tri.triangle_verticies.at(0).vertex_color[2]
					  + beta * tri.triangle_verticies.at(1).vertex_color[2]
					  + alpha * tri.triangle_verticies.at(2).vertex_color[2];
						
					zbuffer.at(i).at(j) = currentZ;
					data[i+j*width] = Make_Pixel(r, g, b);
					// http://www.cs.ucr.edu/~shinar/courses/cs130/content/07-TriangleRasterization.pdf
					// SLIDE 17 (!!)
				}
			}
		}
		
	}

	
}

/* Helper function to return reference to the top of active stack (based on current matrix mode) */

mat4& top_of_active_matrix_stack(){
	if(matrix_mode == MGL_PROJECTION){
		return projection_stack.back();
	}
	else{
		return modelview_stack.back();
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
	// populate zbuffer with 1000 vectors initially
	zbuffer.resize(1000);
	
	// Fill whole pixel data with black initially
	for(unsigned i = 0; i < width; i++){
		zbuffer.at(i).resize(1000, ZMIN);	// initialize background to ZMIN far away from screen
		for(unsigned j = 0; j < height; j++){
			data[i+j*width] = Make_Pixel(0,0,0);
		}
	}
	// For each triangle in the triangle list call Rasterize_Triangle
	for(unsigned i = 0; i < triangles.size(); i++){
		Rasterize_Triangle(triangles.at(i), width, height, data);
	}
	
	// Clear triangle list
	triangles.clear();
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
	
	temp.position = projection_stack.back() * modelview_stack.back() * temp.position;
	
	temp.position /= temp.position[3];
	
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
	mat4 temp = top_of_active_matrix_stack();
	if(matrix_mode == MGL_PROJECTION){
		projection_stack.push_back(temp);
	}
	else{
		modelview_stack.push_back(temp);
	}
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
	if(matrix_mode == MGL_PROJECTION){
		projection_stack.pop_back();
	}
	else{
		modelview_stack.pop_back();
	}
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
	mat4 &temp = top_of_active_matrix_stack();
	temp = {{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}};
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
	mat4 &current_matrix = top_of_active_matrix_stack();
	mat4 temp_matrix;
	float a1 = matrix[0];
	float a2 = matrix[1];
	float a3 = matrix[2];
	float a4 = matrix[3];
	float a5 = matrix[4];
	float a6 = matrix[5];
	float a7 = matrix[6];
	float a8 = matrix[7];
	float a9 = matrix[8];
	float a10 = matrix[9];
	float a11 = matrix[10];
	float a12 = matrix[11];
	float a13 = matrix[12];
	float a14 = matrix[13];
	float a15 = matrix[14];
	float a16 = matrix[15];
	temp_matrix = mat4{{a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16}};
	current_matrix = temp_matrix;
}

/* Helper to multiply 2 4x4 matricies, return result */

mat4 mult_matricies(mat4 a, mat4 b){
	float a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15;
	
	a0 = (a.values[0] * b.values[0]) + (a.values[4] * b.values[1]) + (a.values[8] * b.values[2]) + (a.values[12] * b.values[3]);
	a1 = (a.values[1] * b.values[0]) + (a.values[5] * b.values[1]) + (a.values[9] * b.values[2]) + (a.values[13] * b.values[3]);
	a2 = (a.values[2] * b.values[0]) + (a.values[6] * b.values[1]) + (a.values[10] * b.values[2]) + (a.values[14] * b.values[3]);
	a3 = (a.values[3] * b.values[0]) + (a.values[7] * b.values[1]) + (a.values[11] * b.values[2]) + (a.values[15] * b.values[3]);
	
	a4 = (a.values[0] * b.values[4]) + (a.values[4] * b.values[5]) + (a.values[8] * b.values[6]) + (a.values[12] * b.values[7]);
	a5 = (a.values[1] * b.values[4]) + (a.values[5] * b.values[5]) + (a.values[9] * b.values[6]) + (a.values[13] * b.values[7]);
	a6 = (a.values[2] * b.values[4]) + (a.values[6] * b.values[5]) + (a.values[10] * b.values[6]) + (a.values[14] * b.values[7]);
	a7 = (a.values[3] * b.values[4]) + (a.values[7] * b.values[5]) + (a.values[11] * b.values[6]) + (a.values[15] * b.values[7]);
	
	a8 = (a.values[0] * b.values[8]) + (a.values[4] * b.values[9]) + (a.values[8] * b.values[10]) + (a.values[12] * b.values[11]);
	a9 = (a.values[1] * b.values[8]) + (a.values[5] * b.values[9]) + (a.values[9] * b.values[10]) + (a.values[13] * b.values[11]);
	a10 = (a.values[2] * b.values[8]) + (a.values[6] * b.values[9]) + (a.values[10] * b.values[10]) + (a.values[14] * b.values[11]);
	a11 = (a.values[3] * b.values[8]) + (a.values[7] * b.values[9]) + (a.values[11] * b.values[10]) + (a.values[15] * b.values[11]);
	
	a12 = (a.values[0] * b.values[12]) + (a.values[4] * b.values[13]) + (a.values[8] * b.values[14]) + (a.values[12] * b.values[15]);
	a13 = (a.values[1] * b.values[12]) + (a.values[5] * b.values[13]) + (a.values[9] * b.values[14]) + (a.values[13] * b.values[15]);
	a14 = (a.values[2] * b.values[12]) + (a.values[6] * b.values[13]) + (a.values[10] * b.values[14]) + (a.values[14] * b.values[15]);
	a15 = (a.values[3] * b.values[12]) + (a.values[7] * b.values[13]) + (a.values[11] * b.values[14]) + (a.values[15] * b.values[15]);
	
	mat4 ret = {{a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15}};
	return ret;
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
	mat4 &current_matrix = top_of_active_matrix_stack();
	mat4 temp_matrix;
	float a1 = matrix[0];
	float a2 = matrix[1];
	float a3 = matrix[2];
	float a4 = matrix[3];
	float a5 = matrix[4];
	float a6 = matrix[5];
	float a7 = matrix[6];
	float a8 = matrix[7];
	float a9 = matrix[8];
	float a10 = matrix[9];
	float a11 = matrix[10];
	float a12 = matrix[11];
	float a13 = matrix[12];
	float a14 = matrix[13];
	float a15 = matrix[14];
	float a16 = matrix[15];
	temp_matrix = mat4{{a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16}};
	current_matrix = mult_matricies(current_matrix, temp_matrix);
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
	mat4 &curr = top_of_active_matrix_stack();
	mat4 translateFactor = {{1,0,0,0,0,1,0,0,0,0,1,0,x,y,z,1}};
	
	curr = mult_matricies(curr, translateFactor);
	
}

/* Helper function: convert degrees (float) to radians (float) */

float convert_to_radians(float degrees){
	return degrees * 3.1415926535/180;
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
	// Normalize x,y,z
	
	float a = sqrt((x * x) + (y * y) + (z * z));
	
	x/= a;
	y/= a;
	z/= a;
	
	// (!!) function uses degrees as input, however cos and sin function use radians.
	
	mat4 &curr = top_of_active_matrix_stack();
	float c = cos(convert_to_radians(angle));
	float s = sin(convert_to_radians(angle));
	
	// https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/glRotate.xml
	
	mat4 rotateFactor;
	rotateFactor = {{x*x*(1 - c)+c, y*x*(1-c) + z*s, x*z*(1-c) - y*s, 0,
			x*y*(1-c)-z*s, y*y*(1-c)+c, y*z*(1-c)+x*s, 0,
			x*z*(1-c)+y*s, y*z*(1-c)-x*s, z*z*(1-c)+c, 0,
			0, 0, 0, 1}};
			
	curr = mult_matricies(curr, rotateFactor);		
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
	mat4 &curr = top_of_active_matrix_stack();
	mat4 scaleFactor = {{x, 0, 0, 0, 0, y, 0, 0, 0, 0, z, 0, 0, 0, 0, 1}};
	
	curr = mult_matricies(curr, scaleFactor);
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
	// if(matrix_mode == MGL_PROJECTION){
	// 	frustumFlag = 1;
	// 	float A,B,C,D;
	// 	A = (right + left)/(right - left);
	// 	B = (top + bottom)/(top - bottom);
	// 	C = (far + near)/(far - near);
	// 	D = (2 * far * near)/(far - near);
	// 	projection = {{(2 * near)/(right - left), 0, 0, 0, 0, (2 * near)/(top - bottom), 0, 0, A, B, C, -1, 0, 0, D, 0}};
	// }
	mat4 &temp = top_of_active_matrix_stack();
	float A,B,C,D;
	A = (right + left)/(right - left);
	B = (top + bottom)/(top - bottom);
	C = (far + near)/(far - near);
	D = (2 * far * near)/(far - near);
	temp = {{(2 * near)/(right - left), 0, 0, 0, 0, (2 * near)/(top - bottom), 0, 0, A, B, C, -1, 0, 0, D, 0}};
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
	// if(matrix_mode == MGL_PROJECTION){
	// 	orthoFlag = 1;
	// 	MGLfloat tx, ty, tz;
	// 	tx = -(right + left)/(right - left);
	// 	ty = -(top + bottom)/(top - bottom);
	// 	tz = -(far + near)/(far - near);
	// 	projection = {{(2/(right - left)), 0, 0, 0, 0, (2/(top - bottom)), 0, 0, 0, 0, (-2/(far - near)), 0, tx, ty, tz, 1}};
	// }
	ZMAX = far;
	ZMIN = near;
	mat4 &temp = top_of_active_matrix_stack();
	MGLfloat tx, ty, tz;
	tx = -(right + left)/(right - left);
	ty = -(top + bottom)/(top - bottom);
	tz = -(far + near)/(far - near);
	temp = {{(2/(right - left)), 0, 0, 0, 0, (2/(top - bottom)), 0, 0, 0, 0, (-2/(far - near)), 0, tx, ty, tz, 1}};
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
	color = vec3(red * 255, green * 255, blue * 255);
}

