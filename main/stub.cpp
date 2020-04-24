/******************************************************************************|
| CPSC 4050/6050 Computer Garphics Assignment 5, Daljit Singh Dhillon, 2020    |
| Reference:                                                                   |
|                                                                              |
| Some OpenGL setup code here including math_funcs, and gl_utils               |
| are from Angton Gerdelan and "Anton's OpenGL 4 Tutorials."                   |
| http://antongerdelan.net/opengl/                                             |
| Email: anton at antongerdelan dot net                                        |
| Copyright Dr Anton Gerdelan, Trinity College Dublin, Ireland.                |
|******************************************************************************/
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <assert.h>


#include <math.h>
#include <time.h>

#include "maths_funcs.h"   // Anton's maths functions.
#include "gl_utils.h"      // Anton's opengl functions and small utilities like logs
#include "stb_image.h"     // Sean Barrett's image loader with Anton's load_texture()

#define _USE_MATH_DEFINES
#define ONE_DEG_IN_RAD (2.0 * M_PI) / 360.0 // 0.017444444

#define STEPS_Y (128)
#define STEPS_THETA (64)
#define MESH_LENGTH ((STEPS_Y-1) * (STEPS_THETA) * 18)
#define UV_LENGTH ((STEPS_Y-1) * (STEPS_THETA) * 12)

#define SQ(x) ((x)*(x))

void drawCubicBezier(
	GLfloat *TRIANGLES, 
	GLfloat *NORMALS,
	GLfloat *UV,
	int size, int Y, int THETA);
void get_normal(
	float x, float y, float z,
	float ax,float ay,float az, 
	float bx, float by, float bz, 
	float *ox, float *oy, float *oz);
void normalize(float *x, float *y, float *z);


mat4 view_mat;
mat4 proj_mat;
mat4 model_mat;


int pointCount;


int DIFFUSE_SHADING_location;// = glGetUniformLocation (shader_programme, "DIFFUSE_SHADING");
int SPECULAR_LIGHTING_location;// = glGetUniformLocation (shader_programme, "SPECULAR_LIGHTING");
int DIFFUSE_MAP_location;// = glGetUniformLocation (shader_programme, "DIFFUSE_MAP");
int SPECULAR_EXPONENT_location;// = glGetUniformLocation (shader_programme, "SPECULAR_EXPONENT");
int LIGHT_POSITION_location;

int DIFFUSE = 0;
int DIFFUSEMAP = 0;
int SPEC = 0;
int SPECEX = 0;
float LIGHTPOS[3] = {0.0,0.0,5.0};

void loadSurfaceOfRevolution() 
{
/*------------------------------CREATE GEOMETRY-------------------------------*/
	// GLfloat vp[18];    // array of vertex points
	
	// //face 1, vertex 1
	// vp[0] = -1; //x
	// vp[1] = -1; //y
	// vp[2] = 0; //z
	// //face 1, vertex 2
	// vp[3] = 1; //x
	// vp[4] = -1; //y
	// vp[5] = 0; //z
	// //face 1, vertex 3
	// vp[6] = -1; //x
	// vp[7] =  1; //y
	// vp[8] =  0; //z
	
	// //face 2, vertex 1
	// vp[ 9] = -1; //x
	// vp[10] =  1; //y
	// vp[11] = 0; //z
	// //face 2, vertex 2
	// vp[12] =  1; //x
	// vp[13] = -1; //y
	// vp[14] = 0; //z
	// //face 2, vertex 3
	// vp[15] =  1; //x
	// vp[16] =  1; //y
	// vp[17] =  0; //z

	GLfloat vp[MESH_LENGTH];    // array of vertex points
	GLfloat normal[MESH_LENGTH];	// Array of surface normals
	GLfloat uv[UV_LENGTH]; // UV Texture mapping

	drawCubicBezier(&(vp[0]),&(normal[0]),&(uv[0]),MESH_LENGTH,STEPS_Y, STEPS_THETA);

	// VAO -- vertex attribute objects bundle the various things associated with vertices
	GLuint vao;
	glGenVertexArrays (1, &vao);   // generating and binding is common pattern in OpenGL
	glBindVertexArray (vao);       // basically setting up memory and associating it

	// VBO -- vertex buffer object to contain coordinates
	// MODIFY THE FOLLOWING BLOCK OF CODE APPRORIATELY FOR YOUR SURFACE OF REVOLUTION
	GLuint points_vbo;
	glGenBuffers(1, &points_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
	glBufferData(GL_ARRAY_BUFFER, MESH_LENGTH * sizeof (GLfloat), vp, GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray(0);

	// VBO -- normals -- needed for shading calcuations
	// ADD CODE TO POPULATE AND LOAD PER-VERTEX SURFACE NORMALS  
	// [HINT] Vertex normals are organized in same order as that for vertex coordinates
	GLuint vertex_vbo;
	glGenBuffers(1, &vertex_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vertex_vbo);
	glBufferData(GL_ARRAY_BUFFER, MESH_LENGTH * sizeof (GLfloat), normal, GL_STATIC_DRAW);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray(1);

    // VBO -- vt -- texture coordinates
	// ADD CODE TO POPULATE AND LOAD PER-VERTEX TEXTURE COORDINATES  
	// [HINT] texture coordinates are organized in same order as that for vertex coordinates
	// [HINT] there are two texture coordinates instead of three vertex coordinates for each vertex
	// GLuint texture_vbo;
	// glGenBuffers(1, &texture_vbo);
	// glBindBuffer(GL_ARRAY_BUFFER, texture_vbo);
	// glBufferData(GL_ARRAY_BUFFER, UV_LENGTH * sizeof (GLfloat), uv, GL_STATIC_DRAW);
	// glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, NULL);
	// glEnableVertexAttribArray(2);
}


	
void loadUniforms(GLuint shader_programme)
{	
/*---------------------------SET RENDERING DEFAULTS---------------------------*/

	// Choose vertex and fragment shaders to use as well as view and proj matrices.
	int model_mat_location = glGetUniformLocation (shader_programme, "model_mat");
	int view_mat_location  = glGetUniformLocation (shader_programme, "view_mat");
	int proj_mat_location  = glGetUniformLocation (shader_programme, "proj_mat");
	
	glUniformMatrix4fv (view_mat_location, 1, GL_FALSE, view_mat.m);
	glUniformMatrix4fv (proj_mat_location, 1, GL_FALSE, proj_mat.m);
	glUniformMatrix4fv (model_mat_location, 1, GL_FALSE, model_mat.m);
	
	// WRITE CODE TO LOAD OTHER UNIFORM VARIABLES LIKE FLAGS FOR ENABLING OR DISABLING CERTAIN FUNCTIONALITIES
	// uniform float DIFFUSE_SHADING, SPECULAR_LIGHTING, DIFFUSE_MAP;
	// uniform float SPECULAR_EXPONENT;
	// DIFFUSE_SHADING_location = glGetUniformLocation (shader_programme, "DIFFUSE_SHADING");
	// SPECULAR_LIGHTING_location = glGetUniformLocation (shader_programme, "SPECULAR_LIGHTING");
	// DIFFUSE_MAP_location = glGetUniformLocation (shader_programme, "DIFFUSE_MAP");
	// SPECULAR_EXPONENT_location = glGetUniformLocation (shader_programme, "SPECULAR_EXPONENT");
	// LIGHT_POSITION_location = glGetUniformLocation (shader_programme, "LIGHT_POSITION");
}

void drawSurfaceOfRevolution()
{
	// MODIFY THIS LINE OF CODE APPRORIATELY FOR YOUR SURFACE OF REVOLUTION
	glDrawArrays(GL_TRIANGLES, 0, MESH_LENGTH);
}
	
void keyboardFunction(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	// MODIFY THIS FUNCTION FOR KEYBOARD INTERACTIVITY
	//GLFW Reference Links:
	// Callback Example: https://www.glfw.org/docs/3.3/input_guide.html#input_key
	// List of Keys: https://www.glfw.org/docs/3.3/group__keys.html
	
 //    if (key == GLFW_KEY_A && action == GLFW_PRESS)
 //    {
	// 	printf("\nKey 'A' pressed.... \nDiffuse map %s\n", DIFFUSEMAP % 2 == 0 ? "OFF" : "ON");
	// 	DIFFUSEMAP += 1;
	// 	glUniform1i(DIFFUSE_MAP_location, DIFFUSEMAP);
	// }

	// if (key == GLFW_KEY_S && action == GLFW_PRESS)
	// {
	// 	printf("\nKey 'S' pressed.... \nSpecular %s\n", SPEC % 2 == 0 ? "OFF" : "ON");
 //        SPEC += 1;
 //        glUniform1i(SPECULAR_LIGHTING_location, SPEC);	
	// }
 
	// if (key == GLFW_KEY_D && action == GLFW_PRESS)
	// {
	// 	printf("\nKey 'D' pressed.... \nDiffuse %s\n", DIFFUSE % 2 == 0 ? "OFF" : "ON");
	// 	DIFFUSE += 1;
	// 	glUniform1i(DIFFUSE_SHADING_location, DIFFUSE);
	// }

	// if (key == GLFW_KEY_E && action == GLFW_PRESS)
	// {
	// 	printf("\nKey 'E' pressed.... \nIncreasing Specular Exponent: %d\n", SPECEX);
	// 	SPECEX += 1;
	// 	glUniform1i(SPECULAR_EXPONENT_location, SPECEX);
	// }

	// if (key == GLFW_KEY_RIGHT && action == GLFW_PRESS)
	// {
	// 	LIGHTPOS[0] -= 0.25;
	// 	printf("\nKey 'RIGHT' pressed.... \nMoving light positive x %f\n", LIGHTPOS[0]);
	// 	glUniform3fv(LIGHT_POSITION_location, 1, &(LIGHTPOS[0]));
	// }

	// if (key == GLFW_KEY_LEFT && action == GLFW_PRESS)
	// {
	// 	LIGHTPOS[0] += 0.25;
	// 	printf("\nKey 'LEFT' pressed.... \nMoving light negative x %f\n", LIGHTPOS[0]);
	// 	glUniform3fv(LIGHT_POSITION_location, 1, &(LIGHTPOS[0]));
	// }

	// if (key == GLFW_KEY_UP && action == GLFW_PRESS)
	// {
	// 	LIGHTPOS[1] -= 0.25;
	// 	printf("\nKey 'UP' pressed.... \nMoving light positive y %f\n", LIGHTPOS[1]);
	// 	glUniform3fv(LIGHT_POSITION_location, 1, &(LIGHTPOS[0]));
	// }

	// if (key == GLFW_KEY_DOWN && action == GLFW_PRESS)
	// {
	// 	LIGHTPOS[1] += 0.25;
	// 	printf("\nKey 'DOWN' pressed.... \nMoving light negative y %f\n", LIGHTPOS[1]);
	// 	glUniform3fv(LIGHT_POSITION_location, 1, &(LIGHTPOS[0]));
	// }

	// if (key == GLFW_KEY_Z && action == GLFW_PRESS)
	// {
	// 	LIGHTPOS[2] -= 0.25;
	// 	printf("\nKey 'UP' pressed.... \nMoving light positive Z %f\n", LIGHTPOS[2]);
	// 	glUniform3fv(LIGHT_POSITION_location, 1, &(LIGHTPOS[0]));
	// }

	// if (key == GLFW_KEY_X && action == GLFW_PRESS)
	// {
	// 	LIGHTPOS[2] += 0.25;
	// 	printf("\nKey 'DOWN' pressed.... \nMoving light negative Z %f\n", LIGHTPOS[2]);
	// 	glUniform3fv(LIGHT_POSITION_location, 1, &(LIGHTPOS[0]));
	// }

	// if (key == GLFW_KEY_R && action == GLFW_PRESS)
	// {
	// 	SPECEX -= 1;

	// 	if (SPECEX < 0)
	// 	{
	// 		SPECEX = 0;
	// 	}
	// 	printf("\nKey 'E' pressed.... \nDecreasing Specular Exponent: %d\n", SPECEX);

	// 	glUniform1i(SPECULAR_EXPONENT_location, SPECEX);
	// }

	if (GLFW_PRESS == glfwGetKey (g_window, GLFW_KEY_ESCAPE)) {
		// Close window when esacape is pressed
			glfwSetWindowShouldClose (g_window, 1);
	}

}



//=============================================================================
//=============================================================================
//=============================================================================
//=============================================================================

// Compute a cubic bezier using fixed control points
void drawCubicBezier(
	GLfloat *TRIANGLES, 
	GLfloat *NORMALS,
	GLfloat *UV,
	int size, int Y, int THETA)
{
	float t;
	float x,y,z;
	int incr = 0;

	// START From function createBezierSurface()
	int theta = 0;
	float r;
	// END From function createBezierSurface()

	// START Replace global variables from cg04

	// float BEZIER[STEPS_Y][3];
	float **BEZIER = (float**) calloc(STEPS_Y, sizeof(float*));
	for (int i = 0; i < STEPS_Y; i++)
	{
		BEZIER[i] = (float*) calloc(3, sizeof(float));
	}

	int X_DIM = 0, Y_DIM = 1, Z_DIM = 2;

	float ptX[4] = { 0.0, 1.0, 1.0, 0.0 };
	float ptY[4] = { 1.0, 0.5,-0.5,-1.0 };

	const float PI = 3.1415926;
	float THETA_INCR = 2.0 * PI / (float) STEPS_THETA;

	// float SURFACE[STEPS_Y][STEPS_THETA][3];
	float ***SURFACE = (float***) calloc(STEPS_Y,sizeof(float**));
	for (int i = 0; i < STEPS_Y; i++)
	{
		SURFACE[i] = (float**) calloc(STEPS_THETA, sizeof(float**));

		for (int j = 0; j < STEPS_THETA; j++)
		{
			SURFACE[i][j] = (float*) calloc(3, sizeof(float));
		}
	}

	// END Replace global variables from cg04

	// Each vertex in rows 0 - STEPS_Y-1 creates:
	// 	 2 triangles
	// 	 	1 triangle = 9 XYZ values
	// 	 
	// int triangle_vertex_size = STEPS_THETA * (STEPS_Y - 1) * 2 * 9;
	// float **TRIANGLES = (float*) calloc(triangle_vertex_size,sizeof(float));
	int triangle_index = 0, normal_index = 0, uv_index = 0;
	float u_val = 0.0; // X direction, right to 1.0, upper right at (1.0,1.0) = (u,v)
	float v_val = 1.0; // Y direction, down to 0.0, lower left at (0.0,0.0) = (u,v)
	float u_step = 1.0 / (float) STEPS_THETA;
	float v_step = 1.0 / (float) STEPS_Y;

	// Draw STEPS_Y points
	for (t = 0.0; t < 1.0; t += 1.0 / (float) STEPS_Y)
	{
		// Cubic bezier equaiton for both coordinates
		x = pow(1.0-t,3.0) * ptX[0] + 
			3.0 * t * pow(1.0-t,2)*ptX[1] + 
			3.0 * pow(t,2) * (1.0-t) * ptX[2] + 
			pow(t,3.0) * ptX[3]; 

		y = pow(1.0-t,3.0) * ptY[0] + 
			3.0 * t * pow(1.0-t,2)*ptY[1] + 
			3.0 * pow(t,2) * (1.0-t) * ptY[2] + 
			pow(t,3.0) * ptY[3]; 	

		BEZIER[incr][X_DIM] = x;
		BEZIER[incr][Y_DIM] = y;
		BEZIER[incr][Z_DIM] = 0.0;
		incr++;
	}


	// Create surface of revolution
	for (int y = 0; y < STEPS_Y; y++)
	{
		theta = 0;
		r = BEZIER[y][X_DIM] * BEZIER[y][X_DIM];

		for (float t = THETA_INCR; t < 2.0 * PI; t += THETA_INCR)
		{
			x = r * cos(t);
			z = r * sin(t);

			SURFACE[y][theta][X_DIM] = x;
			SURFACE[y][theta][Y_DIM] = BEZIER[y][Y_DIM];
			SURFACE[y][theta][Z_DIM] = z;

			theta++;
		}
	}

	float nx,ny,nz;

	for (int row = 0; row < STEPS_Y-1; row++)
	{
		for (int theta_incr = 0; theta_incr < STEPS_THETA; theta_incr++)
		{
			// Bottom Triangle
			// 3
			// 1	2
			TRIANGLES[triangle_index] = SURFACE[row+1][theta_incr][X_DIM]; triangle_index++;
			TRIANGLES[triangle_index] = SURFACE[row+1][theta_incr][Y_DIM]; triangle_index++;
			TRIANGLES[triangle_index] = SURFACE[row+1][theta_incr][Z_DIM]; triangle_index++;

			TRIANGLES[triangle_index] = SURFACE[row+1][(theta_incr+1) % STEPS_THETA][X_DIM]; triangle_index++;
			TRIANGLES[triangle_index] = SURFACE[row+1][(theta_incr+1) % STEPS_THETA][Y_DIM]; triangle_index++;
			TRIANGLES[triangle_index] = SURFACE[row+1][(theta_incr+1) % STEPS_THETA][Z_DIM]; triangle_index++;	
			
			TRIANGLES[triangle_index] = SURFACE[row][theta_incr][X_DIM]; triangle_index++;
			TRIANGLES[triangle_index] = SURFACE[row][theta_incr][Y_DIM]; triangle_index++;
			TRIANGLES[triangle_index] = SURFACE[row][theta_incr][Z_DIM]; triangle_index++;		

			// NORMALS
			get_normal(
				SURFACE[row+1][theta_incr][X_DIM], SURFACE[row+1][theta_incr][Y_DIM], SURFACE[row+1][theta_incr][Z_DIM],
				SURFACE[row+1][(theta_incr+1) % STEPS_THETA][X_DIM], SURFACE[row+1][(theta_incr+1) % STEPS_THETA][Y_DIM], SURFACE[row+1][(theta_incr+1) % STEPS_THETA][Z_DIM],
				SURFACE[row][theta_incr][X_DIM], SURFACE[row][theta_incr][Y_DIM], SURFACE[row][theta_incr][Z_DIM],
				&nx,&ny,&nz);

			NORMALS[normal_index] = nx; normal_index++;
			NORMALS[normal_index] = ny; normal_index++;
			NORMALS[normal_index] = nz; normal_index++;

			get_normal(
				SURFACE[row+1][(theta_incr+1) % STEPS_THETA][X_DIM], SURFACE[row+1][(theta_incr+1) % STEPS_THETA][Y_DIM], SURFACE[row+1][(theta_incr+1) % STEPS_THETA][Z_DIM],
				SURFACE[row][theta_incr][X_DIM], SURFACE[row][theta_incr][Y_DIM], SURFACE[row][theta_incr][Z_DIM],
				SURFACE[row+1][theta_incr][X_DIM], SURFACE[row+1][theta_incr][Y_DIM], SURFACE[row+1][theta_incr][Z_DIM],
				&nx,&ny,&nz);

			NORMALS[normal_index] = nx; normal_index++;
			NORMALS[normal_index] = ny; normal_index++;
			NORMALS[normal_index] = nz; normal_index++;

			get_normal(
				SURFACE[row][theta_incr][X_DIM], SURFACE[row][theta_incr][Y_DIM], SURFACE[row][theta_incr][Z_DIM],
				SURFACE[row+1][theta_incr][X_DIM], SURFACE[row+1][theta_incr][Y_DIM], SURFACE[row+1][theta_incr][Z_DIM],
				SURFACE[row+1][(theta_incr+1) % STEPS_THETA][X_DIM], SURFACE[row+1][(theta_incr+1) % STEPS_THETA][Y_DIM], SURFACE[row+1][(theta_incr+1) % STEPS_THETA][Z_DIM],
				&nx,&ny,&nz);

			NORMALS[normal_index] = nx; normal_index++;
			NORMALS[normal_index] = ny; normal_index++;
			NORMALS[normal_index] = nz; normal_index++;

			// UV Coords
			UV[uv_index] = u_val; uv_index++;
			UV[uv_index] = v_val - v_step; uv_index++; // row + 1

			UV[uv_index] = u_val + u_step; uv_index++; // row & col + 1
			UV[uv_index] = v_val - v_step; uv_index++;

			UV[uv_index] = u_val; uv_index++;
			UV[uv_index] = v_val; uv_index++;




			// Top triangle
			// 1	3
			// 		2
			TRIANGLES[triangle_index] = SURFACE[row][theta_incr][X_DIM]; triangle_index++;
			TRIANGLES[triangle_index] = SURFACE[row][theta_incr][Y_DIM]; triangle_index++;
			TRIANGLES[triangle_index] = SURFACE[row][theta_incr][Z_DIM]; triangle_index++;

			TRIANGLES[triangle_index] = SURFACE[row+1][(theta_incr+1) % STEPS_THETA][X_DIM]; triangle_index++;
			TRIANGLES[triangle_index] = SURFACE[row+1][(theta_incr+1) % STEPS_THETA][Y_DIM]; triangle_index++;
			TRIANGLES[triangle_index] = SURFACE[row+1][(theta_incr+1) % STEPS_THETA][Z_DIM]; triangle_index++;	
			
			TRIANGLES[triangle_index] = SURFACE[row][(theta_incr+1) % STEPS_THETA][X_DIM]; triangle_index++;
			TRIANGLES[triangle_index] = SURFACE[row][(theta_incr+1) % STEPS_THETA][Y_DIM]; triangle_index++;
			TRIANGLES[triangle_index] = SURFACE[row][(theta_incr+1) % STEPS_THETA][Z_DIM]; triangle_index++;

			get_normal(
				SURFACE[row][theta_incr][X_DIM], SURFACE[row][theta_incr][Y_DIM], SURFACE[row][theta_incr][Z_DIM],
				SURFACE[row+1][(theta_incr+1) % STEPS_THETA][X_DIM], SURFACE[row+1][(theta_incr+1) % STEPS_THETA][Y_DIM], SURFACE[row+1][(theta_incr+1) % STEPS_THETA][Z_DIM],
				SURFACE[row][(theta_incr+1) % STEPS_THETA][X_DIM], SURFACE[row][(theta_incr+1) % STEPS_THETA][Y_DIM], SURFACE[row][(theta_incr+1) % STEPS_THETA][Z_DIM],
				&nx,&ny,&nz);

			NORMALS[normal_index] = nx; normal_index++;
			NORMALS[normal_index] = ny; normal_index++;
			NORMALS[normal_index] = nz; normal_index++;

			get_normal(
				SURFACE[row+1][(theta_incr+1) % STEPS_THETA][X_DIM], SURFACE[row+1][(theta_incr+1) % STEPS_THETA][Y_DIM], SURFACE[row+1][(theta_incr+1) % STEPS_THETA][Z_DIM],
				SURFACE[row][(theta_incr+1) % STEPS_THETA][X_DIM], SURFACE[row][(theta_incr+1) % STEPS_THETA][Y_DIM], SURFACE[row][(theta_incr+1) % STEPS_THETA][Z_DIM],
				SURFACE[row][theta_incr][X_DIM], SURFACE[row][theta_incr][Y_DIM], SURFACE[row][theta_incr][Z_DIM],
				&nx,&ny,&nz);

			NORMALS[normal_index] = nx; normal_index++;
			NORMALS[normal_index] = ny; normal_index++;
			NORMALS[normal_index] = nz; normal_index++;

			get_normal(
				SURFACE[row][(theta_incr+1) % STEPS_THETA][X_DIM], SURFACE[row][(theta_incr+1) % STEPS_THETA][Y_DIM], SURFACE[row][(theta_incr+1) % STEPS_THETA][Z_DIM],
				SURFACE[row][theta_incr][X_DIM], SURFACE[row][theta_incr][Y_DIM], SURFACE[row][theta_incr][Z_DIM],
				SURFACE[row+1][(theta_incr+1) % STEPS_THETA][X_DIM], SURFACE[row+1][(theta_incr+1) % STEPS_THETA][Y_DIM], SURFACE[row+1][(theta_incr+1) % STEPS_THETA][Z_DIM],
				&nx,&ny,&nz);

			NORMALS[normal_index] = nx; normal_index++;
			NORMALS[normal_index] = ny; normal_index++;
			NORMALS[normal_index] = nz; normal_index++;

			// UV Coords
			UV[uv_index] = u_val; uv_index++;
			UV[uv_index] = v_val; uv_index++;

			UV[uv_index] = u_val + u_step; uv_index++; // row & col + 1
			UV[uv_index] = v_val - v_step; uv_index++;

			UV[uv_index] = u_val + u_step; uv_index++; // col + 1
			UV[uv_index] = v_val; uv_index++;

			u_val += u_step; // Increase base u value
		}

		v_val -= v_step;
		u_val = 0.0;
	}

	for (int i = 0; i < STEPS_Y; i++)
	{
		for (int j = 0; j < STEPS_THETA; j++)
		{
			free(SURFACE[i][j]);
		}

		free(SURFACE[i]);
		free(BEZIER[i]);
	}
	free(SURFACE);
	free(BEZIER);
}


void get_normal(
	float x, float y, float z,
	float ax,float ay,float az, 
	float bx, float by, float bz, 
	float *ox, float *oy, float *oz)
{
	float tax = ax - x;
	float tay = ay - y;
	float taz = az - z;
	float tbx = bx - x;
	float tby = by - y;
	float tbz = bz - z;

	*ox = tay * tbz - taz * tby;
	*oy = taz * tbx - tax * tbz;
	*oz = tax * tby - tay * tbx;

	normalize(ox,oy,oz);
}

void normalize(float *x, float *y, float *z)
{
	float tx = *x, ty = *y, tz = *z;

	float norm = sqrt(SQ(tx) + SQ(ty) + SQ(tz));

	// avoid NaN
	if (norm == 0.0)
	{
		norm = 0.00000001;
	}

	*x = tx / norm;
	*y = ty / norm;
	*z = tz / norm;
}








