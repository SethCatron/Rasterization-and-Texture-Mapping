#pragma once

#include <stdlib.h>
#include <math.h>
#include <vector>

#include <glm/glm.hpp>


class Triangle {
	private:
		glm::vec3 v[3];		// Triangle vertices
		glm::vec3 c[3];		// Vertex color
		glm::vec2 t[3];		// Texture coordinates

	public:

		// Default constructor
		Triangle();

		// Constructor without texture coordinates
		Triangle(glm::vec3& v0, glm::vec3& v1, glm::vec3& v2);

		// Constructor with texture coordinates
		Triangle(glm::vec3& v0, glm::vec3& v1, glm::vec3& v2, glm::vec2& t0, glm::vec2& t1, glm::vec2& t2);

		void setColorMode0();
		void setColorMode1();
		void setColorMode2();

		glm::vec3 bilinearcalc(glm::vec2 textureCoordinate, std::vector<float*> texture, int width, int d);
		glm::vec3 lerpcalc(glm::vec3 v0, glm::vec3 v1, float s);
		float lcalc(float dudx, float dudy, float dvdx, float dvdy);
		float dcalc();
		// Rendering the triangle using OpenGL
		void RenderOpenGL(glm::mat4 &modelViewMatrix, glm::mat4 &projectionMatrix, bool textureMode);

		// Rendering the triangle using CPU
		void RenderCPU(glm::mat4 &modelViewMatrix, glm::mat4 &projectionMatrix,  float (&pixelColor)[1024][1024][3], float (&depth)[1024][1024], bool activeTexture, int width, int height, std::vector<float*> texture, int textureMode);
};
