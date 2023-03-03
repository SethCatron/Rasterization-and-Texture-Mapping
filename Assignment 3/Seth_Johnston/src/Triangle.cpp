#include "Triangle.h"
#include <GL/glew.h>
#include <glm/gtc/type_ptr.hpp>
#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

// A function clamping the input values to the lower and higher bounds
#define CLAMP(in, low, high) ((in) < (low) ? (low) : ((in) > (high) ? (high) : in))

Triangle::Triangle()
{
	v[0] = glm::vec3(0.0f, 0.0f, 0.0f);
	v[1] = glm::vec3(0.0f, 0.0f, 0.0f);
	v[2] = glm::vec3(0.0f, 0.0f, 0.0f);

	c[0] = glm::vec3(0.0f, 0.0f, 0.0f);
	c[1] = glm::vec3(0.0f, 0.0f, 0.0f);
	c[2] = glm::vec3(0.0f, 0.0f, 0.0f);

	t[0] = glm::vec2(0.0f, 0.0f);
	t[1] = glm::vec2(0.0f, 0.0f);
	t[2] = glm::vec2(0.0f, 0.0f);
}

Triangle::Triangle(glm::vec3& v0, glm::vec3& v1, glm::vec3& v2)
{
	v[0] = v0;
	v[1] = v1;
	v[2] = v2;

	c[0] = glm::vec3(1.0f, 1.0f, 1.0f);
	c[1] = glm::vec3(1.0f, 1.0f, 1.0f);
	c[2] = glm::vec3(1.0f, 1.0f, 1.0f);

	t[0] = glm::vec2(0.0f, 0.0f);
	t[1] = glm::vec2(0.0f, 0.0f);
	t[2] = glm::vec2(0.0f, 0.0f);

};

Triangle::Triangle(glm::vec3& v0, glm::vec3& v1, glm::vec3& v2, glm::vec2& t0, glm::vec2& t1, glm::vec2& t2)
{
	v[0] = v0;
	v[1] = v1;
	v[2] = v2;

	t[0] = t0;
	t[1] = t1;
	t[2] = t2;

	c[0] = glm::vec3(1.0f, 1.0f, 1.0f);
	c[1] = glm::vec3(1.0f, 1.0f, 1.0f);
	c[2] = glm::vec3(1.0f, 1.0f, 1.0f);
};

void Triangle::setColorMode0() {
	float redcolor = ((float)rand()/RAND_MAX);
	float greencolor = ((float)rand()/RAND_MAX);
	float bluecolor = ((float)rand()/RAND_MAX);
	this->c[0] = glm::vec3(redcolor, greencolor, bluecolor);
	this->c[1] = glm::vec3(redcolor, greencolor, bluecolor);
	this->c[2] = glm::vec3(redcolor, greencolor, bluecolor);
}

void Triangle::setColorMode1() {
	// random color to each vertex
	float redcolor1 = ((float)rand()/RAND_MAX);
	float greencolor1 = ((float)rand()/RAND_MAX);
	float bluecolor1 = ((float)rand()/RAND_MAX);
	float redcolor2 = ((float)rand()/RAND_MAX);
	float greencolor2 = ((float)rand()/RAND_MAX);
	float bluecolor2 = ((float)rand()/RAND_MAX);
	float redcolor3 = ((float)rand()/RAND_MAX);
	float greencolor3 = ((float)rand()/RAND_MAX);
	float bluecolor3 = ((float)rand()/RAND_MAX);
	this->c[0] = glm::vec3(redcolor1, greencolor1, bluecolor1);
	this->c[1] = glm::vec3(redcolor2, greencolor2, bluecolor2);
	this->c[2] = glm::vec3(redcolor3, greencolor3, bluecolor3);
}

float globalmaxz = 0.0;
float globalminz = FLT_MAX;
std::vector<float> minvec = {};
std::vector<float> maxvec = {};
void Triangle::setColorMode2() {

	
	float localminz = FLT_MAX;
	float localmaxz = 0.0;

	float zval0 = this->v[0].z; 
	float zval1 = this->v[1].z;
	float zval2 = this->v[2].z;

	std::vector<float> localzvals = {zval0, zval1, zval2};

	for (int i = 0; i < localzvals.size(); i++) {
		if (localzvals[i] < localminz) {
			localminz = localzvals[i];
		} 
		if (localzvals[i] > localmaxz) {
			localmaxz = localzvals[i];
		}
	}

	maxvec.push_back(localmaxz);
	minvec.push_back(localminz);

	for (int i = 0; i < minvec.size(); i++) {
		if (minvec[i] < globalminz) {
			globalminz = minvec[i];
		} 
	}
	
	for (int i = 0; i < maxvec.size(); i++) {
		if (maxvec[i] > globalmaxz) {
			globalmaxz = maxvec[i];
		} 
	}

	// find max - min = global difference
	// find local max - min = local z
	// 1 - (localz / global z) = red float value
	float globaldiff = globalmaxz - globalminz;
	float localdiff = localmaxz - localminz;

	float calc = 1 - (localdiff / globaldiff);

	this->c[0] = glm::vec3(calc, 0, 0);
	this->c[1] = glm::vec3(calc, 0, 0);
	this->c[2] = glm::vec3(calc, 0, 0);
}

// 1D Linear Interpolation
glm::vec3 Triangle::lerpcalc(glm::vec3 v0, glm::vec3 v1, float s) {
	return (v0 + (s * (v1 - v0)));
}

glm::vec3 Triangle::bilinearcalc(glm::vec2 textureCoordinate, std::vector<float*> texture, int width, int d = 0) {
	
	glm::vec2 topRight = glm::vec2(floor(textureCoordinate[0] + 0.5), floor(textureCoordinate[1] + 0.5));
	glm::vec3 u11;
	glm::vec3 u01;
	glm::vec3 u10;
	glm::vec3 u00;

    u11[0] = texture[d][(int)(topRight[1] * width * 3 + topRight[0]* 3 + 0)];
	u11[1] = texture[d][(int)(topRight[1] * width * 3 + topRight[0]* 3 + 1)];
	u11[2] = texture[d][(int)(topRight[1] * width * 3 + topRight[0]* 3 + 2)];

	u01[0] = texture[d][(int)((topRight[1] - 1) * width * 3 + topRight[0]* 3 + 0)];
	u01[1] = texture[d][(int)((topRight[1] - 1) * width * 3 + topRight[0]* 3 + 1)];
	u01[2] = texture[d][(int)((topRight[1] - 1) * width * 3 + topRight[0]* 3 + 2)];

	u10[0] = texture[d][(int)((topRight[1]) * width * 3 + (topRight[0] - 1) * 3 + 0)];
	u10[1] = texture[d][(int)((topRight[1]) * width * 3 + (topRight[0] - 1) * 3 + 1)];
	u10[2] = texture[d][(int)((topRight[1]) * width * 3 + (topRight[0] - 1) * 3 + 2)];

	u00[0] = texture[d][(int)((topRight[1] - 1) * width * 3 + (topRight[0] - 1) * 3 + 0)];
	u00[1] = texture[d][(int)((topRight[1] - 1) * width * 3 + (topRight[0] - 1) * 3 + 1)];
	u00[2] = texture[d][(int)((topRight[1] - 1) * width * 3 + (topRight[0] - 1) * 3 + 2)];

    // get distances to bottom left from p
    float t = textureCoordinate[1] - (topRight[1] - 1) - 0.5;
    float s = textureCoordinate[0] - (topRight[0] - 1) - 0.5;

    // helper lerps
    glm::vec3 u0 = lerpcalc(u00, u10, s);
    glm::vec3 u1 = lerpcalc(u01, u11, s);

	// final vertical lerp calculation
    return lerpcalc(u0, u1, t);
}

float Triangle::lcalc(float dudx, float dudy, float dvdx, float dvdy) {
	float x = sqrt(pow(dudx, 2) + pow(dvdx, 2));
	float y = sqrt(pow(dudy, 2) + pow(dvdy, 2));
	float max = x;
	if (y > x) {
		max = y;
	}
	return max;
}


// Rendering the triangle using OpenGL
void Triangle::RenderOpenGL(glm::mat4 &modelViewMatrix, glm::mat4 &projectionMatrix, bool isTextured)
{

	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixf(glm::value_ptr(modelViewMatrix));

	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(glm::value_ptr(projectionMatrix));
	
	// For textured object
	if (isTextured)
	{
		glEnable(GL_TEXTURE_2D);

		// Avoid modulating the texture by vertex color
		glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

		glBegin(GL_TRIANGLES);

			glTexCoord2f(t[0].x, t[0].y);
			glVertex3f(v[0].x, v[0].y, v[0].z);

			glTexCoord2f(t[1].x, t[1].y);
			glVertex3f(v[1].x, v[1].y, v[1].z);

			glTexCoord2f(t[2].x, t[2].y);
			glVertex3f(v[2].x, v[2].y, v[2].z);

		glEnd();

		glDisable(GL_TEXTURE_2D);


	}
	// For object with only vertex color
	else
	{
		glBegin(GL_TRIANGLES);

			glColor3f(c[0].x, c[0].y, c[0].z);
			glVertex3f(v[0].x, v[0].y, v[0].z);

			glColor3f(c[1].x, c[1].y, c[1].z);
			glVertex3f(v[1].x, v[1].y, v[1].z);

			glColor3f(c[2].x, c[2].y, c[2].z);
			glVertex3f(v[2].x, v[2].y, v[2].z);
		
		glEnd();
	}

}

// Render the triangle on CPU
// Inside this function you have access to the 3 vertices of the triangle for which this function is called.
// For the first step, you need to pass the modelview and projection matrices from the Display code into RenderCPU code and then perform the transformations

void Triangle::RenderCPU(glm::mat4 &modelViewMatrix, glm::mat4 &projectionMatrix, float (&pixelColor)[1024][1024][3], float (&depth)[1024][1024], bool activeTexture, int width, int height, std::vector<float*> texture, int textureMode)
{
	// mat3 into homogenous mat4 vectors
	glm::vec4 v0(this->v[0], 1.0);
	glm::vec4 v1(this->v[1], 1.0);
	glm::vec4 v2(this->v[2], 1.0);
	
	// Must create the viewport matrix as identity based on screen resolution
	glm::mat4 viewport(1.0f);
	viewport[0][0] = 512.0;
	viewport[1][1] = 512.0;
	viewport[3][0] = 512.0;
	viewport[3][1] = 512.0;

	// Transform the Triangle to the screen  -- transformations done in homogeneous coordinate (4x4 matrices) --> PVM * each vertex = ndc  --> P=ortho, V=lookat, M=perspective
		
		// first need to apply model view projection transformation to bring them to normalize device coordinate (NDC)
		v0 = modelViewMatrix * v0;
		v1 = modelViewMatrix * v1;
		v2 = modelViewMatrix * v2;

		// used for perspective interpolation for the texture P2
		float zv0 = v0[2];
		float zv1 = v1[2];
		float zv2 = v2[2];

		// multiply each vertex by projection matrix and then normalize by dividing by the w value of the projection matrix which corresponds to projectionMatrix[3]
		v0 = projectionMatrix * v0; // normalize this output
		v0 = v0 / v0[3];
		v1 = projectionMatrix * v1; // normalize this output
		v1 = v1 / v1[3];
		v2 = projectionMatrix * v2; // normalize this output
		v2 = v2 / v2[3];

		// Finally you apply viewport transformation to go from NDC to screen space.
		// bound checking to multiply by viewport 
		v0 = viewport * v0;
		v1 = viewport * v1;
		v2 = viewport * v2;

	// Rasterize the Triangle
		
		// loop over the pixels on the screen
		int minboundx = INT32_MAX;
		int maxboundx = 0.0;
		int minboundy = INT32_MAX;
		int maxboundy = 0.0;

		std::vector<float> xvals = {v0[0], v1[0], v2[0]};
		std::vector<float> yvals = {v0[1], v1[1], v2[1]};

		// xmin and xmax finder
		for (int i = 0; i < xvals.size(); i++) {
			if (truncf(xvals[i]) < minboundx) {
				minboundx = truncf(xvals[i]);
			}
			if (ceil(xvals[i]) > maxboundx) {
				maxboundx = ceil(xvals[i]);
			}
			if (truncf(yvals[i]) < minboundy) {
				minboundy = truncf(yvals[i]);
			}
			if (ceil(yvals[i]) > maxboundy) {
				maxboundy = ceil(yvals[i]);
			}
		}

		glm::vec3 c0(this->c[0]);
		glm::vec3 c1(this->c[1]);
		glm::vec3 c2(this->c[2]);

		glm::vec2 t0(this->t[0]);
		glm::vec2 t1(this->t[1]);
		glm::vec2 t2(this->t[2]);
		

		float alpha = 0.0;
		float beta = 0.0;
		float gamma = 0.0;
		
		#define CLAMP(in, low, high) ((in) < (low) ? (low) : ((in) > (high) ? (high) : in))
		maxboundx = CLAMP(maxboundx, 0, 1024 - 1);
		maxboundy = CLAMP(maxboundy, 0, 1024 - 1);
		minboundx = CLAMP(minboundx, 0, 1024 - 1);
		minboundy = CLAMP(minboundy, 0, 1024 - 1);


		for (int x = minboundx; x < maxboundx; x++) {
			for (int y = minboundy; y < maxboundy; y++) {
				// x
				float xa = v0[0];
				float xb = v1[0];
				float xc = v2[0];
				float xp = (float)x + 0.5;
				// y
				float ya = v0[1];
				float yb = v1[1]; 
				float yc = v2[1];  
				float yp = (float)y + 0.5;
				// alpha beta gamma calculations
				alpha = (((-1)*(xp-xb)*(yc-yb)) + ((yp-yb)*(xc-xb))) / (((-1)*(xa-xb)*(yc-yb)) + ((ya-yb)*(xc-xb)));
				// cout << "Alpha " << alpha << endl;
				beta = (((-1)*(xp-xc)*(ya-yc)) + ((yp-yc)*(xa-xc))) / (((-1)*(xb-xc)*(ya-yc)) + ((yb-yc)*(xa-xc)));
				gamma = 1.0 - alpha - beta;

				
				if (alpha <= 1.0 && alpha >= 0 && beta <= 1.0 && beta >= 0 && gamma <= 1.0 && gamma >= 0) {
					float pixelDepth = v0[2]*alpha + v1[2]*beta + v2[2]*gamma;
					if (pixelDepth < depth[y][x]) {
						depth[y][x] = pixelDepth;
						if (activeTexture == true) {

////////////////////////////// NEAREST NEIGHBOR ////////////////////////////////////
							if (textureMode == 0) {

////////////////////////// COMPUTE INTERPOLATED TEXTURE COORDINATE //////////////

								// calculate z inverses
								float zv0inv = 1.0/zv0;
								float zv1inv = 1.0/zv1;
								float zv2inv = 1.0/zv2;

								// divide t values by z to prepare for calculations
								glm::vec2 q0 = t0 / zv0;
								glm::vec2 q1 = t1 / zv1;
								glm::vec2 q2 = t2 / zv2;
								
								// interpolate z-inverses to get a float value
								float zint = zv0inv*alpha + zv1inv*beta + zv2inv*gamma;
								zint = 1/zint;

								// interpolate Qsca values 
								glm::vec2 qint = q0*alpha + q1*beta + q2*gamma;

								// create final cacluation vec2 which divides Qsca vec2 by zint interpolated float
								glm::vec2 final = qint*zint;

////////////////////////////////// WRAPAROUND CALCULATIONS ////////////////////////////
								for (int i = 0; i <= 1; i++) {
									// wrap around texture value > 1
									if (final[i] >= 1) {
										float val = final[i] - (int)final[i];
										final[i] = val;
									}
									// wrap around texture value < 0
									else if (final[i] < 0) {
										float val = final[i] - (int)final[i];
										final[i] = 1.0 + val;
									} 
								}

////////////////////////////////// SCALING CALCULATIONS ////////////////////////////////////
								final[0] = floor(final[0] * width);
								final[1] = floor(final[1] * height);
								// looking up the corresponding color in the texture image and assigning the interpolated texture colors to the pixel colors
								pixelColor[y][x][0] = texture[0][((int)final[0] * 3) + ((int)final[1] * width * 3) + 0];  // red
								pixelColor[y][x][1] = texture[0][((int)final[0] * 3) + ((int)final[1] * width * 3) + 1];  // green
								pixelColor[y][x][2] = texture[0][((int)final[0] * 3) + ((int)final[1] * width * 3) + 2];  // blue
						
////////////////////////////// BILINEAR INTERPOLATION //////////////////////////////
							} else if (textureMode == 1) {
								////////////////////////// COMPUTE INTERPOLATED TEXTURE COORDINATE //////////////

								// calculate z inverses
								float zv0inv = 1.0/zv0;
								float zv1inv = 1.0/zv1;
								float zv2inv = 1.0/zv2;

								// divide t values by z to prepare for calculations
								glm::vec2 q0 = t0 / zv0;
								glm::vec2 q1 = t1 / zv1;
								glm::vec2 q2 = t2 / zv2;
								
								// interpolate z-inverses to get a float value
								float zint = zv0inv*alpha + zv1inv*beta + zv2inv*gamma;
								zint = 1/zint;

								// interpolate Qsca values 
								glm::vec2 qint = q0*alpha + q1*beta + q2*gamma;

								// create final cacluation vec2 which divides Qsca vec2 by zint interpolated float
								glm::vec2 final = qint*zint;

////////////////////////////////// WRAPAROUND CALCULATIONS ////////////////////////////
								for (int i = 0; i <= 1; i++) {
									// wrap around texture value > 1
									if (final[i] >= 1) {
										float val = final[i] - (int)final[i];
										final[i] = val;
									}
									// wrap around texture value < 0
									else if (final[i] < 0) {
										float val = final[i] - (int)final[i];
										final[i] = 1.0 + val;
									} 
								}

////////////////////////////////// SCALING CALCULATIONS ////////////////////////////////////
								final[0] = floor(final[0] * width);
								final[1] = floor(final[1] * height);
								// when calculating which 4 pixels to look for you add .5 before flooring it (looking at the center of the pixels)
								int d = 0;
								glm::vec3 finalColor = bilinearcalc(final, texture, width, d);
								pixelColor[y][x][0] = finalColor[0];
								pixelColor[y][x][1] = finalColor[1];
								pixelColor[y][x][2] = finalColor[2];
								// mix all 4 colors no using the coordinates, but i am mixing the colors!!! nneed a way to sample the color of those extures, given those tx and ty how do i get the colors of the texel? bottom left, top left, bottom right, top right texels
								
////////////////////////////////// MIPMAPPING //////////////////////////////////////////
							} else if (textureMode == 2) {
								
////////////////////////////////// COMPUTE INTERPOLATED TEXTURE COORDINATE //////////////

								// calculate z inverses
								float zv0inv = 1.0/zv0;
								float zv1inv = 1.0/zv1;
								float zv2inv = 1.0/zv2;

								// divide t values by z to prepare for calculations
								glm::vec2 q0 = t0 / zv0;
								glm::vec2 q1 = t1 / zv1;
								glm::vec2 q2 = t2 / zv2;
								
								// interpolate z-inverses to get a float value
								float zint = zv0inv*alpha + zv1inv*beta + zv2inv*gamma;
								zint = 1/zint;

								// interpolate Qsca values 
								glm::vec2 qint = q0*alpha + q1*beta + q2*gamma;

								// create final cacluation vec2 which divides Qsca vec2 by zint interpolated float
								glm::vec2 final = qint*zint;

////////////////////////////////// SCALING CALCULATIONS ////////////////////////////////////
								final[0] = floor(final[0] * width);
								final[1] = floor(final[1] * height);

								// calculate L value with 2D vector calculations -- compute after scaling
								glm::vec2 addx = final;
								addx[0] = final[0]+1;
								glm::vec2 newx = addx - final;
								float dudx = newx[0];
								float dudy = newx[1];

								glm::vec2 addy = final;
								addy[1] = final[1]+1;
								glm::vec2 newy = addy - final;
								float dvdx = newy[0];
								float dvdy = newy[1];
								float newL = lcalc(dudx, dudy, dvdx, dvdy); // how many pixels are in the square

								// FIND MIPMAP LEVEL D
								float d = log2(newL);

								// ROUND D VALUE with ceil and floor to get two D values
								int upperd = (int)ceil(d);
								int lowerd = (int)floor(d);

								// do bilinear for each D using scaled final values resulting in C1 and C2 floats for each bilinear operation
								glm::vec3 C1 = bilinearcalc(final, texture, width, lowerd);
								glm::vec3 C2 = bilinearcalc(final, texture, width, upperd);

								// linear interpolation on each C value to get a final C which is used in the r g b value assignments for the pixelColor buffer
								float C1int = alpha*C1[0] + beta*C1[1] + gamma*C1[2];
								float C2int = alpha*C2[0] + beta*C2[1] + gamma*C2[2];
								float Cfinal = (int)(C1int/C2int);

								// pixelColor assignments
								pixelColor[y][x][0] = texture[Cfinal][((int)final[0] * 3) + ((int)final[1] * width * 3) + 0];
								pixelColor[y][x][1] = texture[Cfinal][((int)final[0] * 3) + ((int)final[1] * width * 3) + 1];
								pixelColor[y][x][2] = texture[Cfinal][((int)final[0] * 3) + ((int)final[1] * width * 3) + 2];

////////////////////////////////// WRAPAROUND CALCULATIONS ////////////////////////////
								for (int i = 0; i <= 1; i++) {
									// wrap around texture value > 1
									if (final[i] >= 1) {
										float val = final[i] - (int)final[i];
										final[i] = val;
									}
									// wrap around texture value < 0
									else if (final[i] < 0) {
										float val = final[i] - (int)final[i];
										final[i] = 1.0 + val;
									} 
								}


							}

///////////////////////// // Interpolate the color of the Pixel using Barycentric Coordinates ////////////////////
						} else { 
							pixelColor[y][x][0] = (alpha*c0[0]) + (beta*c1[0]) + (gamma*c2[0]);
							// depth[y][x] = v0[2];
						
							pixelColor[y][x][1] = (alpha*c0[1]) + (beta*c1[1]) + (gamma*c2[1]);
							// depth[y][x] = v1[2];
						
							pixelColor[y][x][2] = (alpha*c0[2]) + (beta*c1[2]) + (gamma*c2[2]);
							// depth[y][x] = v2[2];
						}
					}
				}
			}
		}
}