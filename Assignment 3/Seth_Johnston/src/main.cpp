#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "stb_image_resize.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
#include "Triangle.h"

using namespace std;

#define WINDOW_WIDTH 1024
#define WINDOW_HEIGHT 1024

GLFWwindow *window;
bool lButtonPressed;
bool rButtonPressed;



float color[WINDOW_HEIGHT][WINDOW_WIDTH][3];  // draw the current triangle onto this
float depth[WINDOW_HEIGHT][WINDOW_WIDTH];  // use to implemeny z-buffer


std::vector<Triangle> triangleVector;
std::vector<float*> texture; // stores the texture image at multiple scales

bool isOpenGL = true;
bool isTextured = false;
float eyeDistance = 5.0f;
int textureMode = 0;
int colorMode = 0;
float angle = 0;

std::string mainName = "Assignment3 - Seth Johnston";

int texWidth, texHeight;

GLuint texID;

void ClearFrameBuffer()
{
	memset(&color[0][0][0], 0.0f, sizeof(float) * WINDOW_WIDTH * WINDOW_HEIGHT * 3);
}

void Display()
{	
	glm::mat4 projectionMatrix = glm::perspective(glm::radians(60.0f), float(WINDOW_WIDTH) / float(WINDOW_HEIGHT), 0.1f, 100.0f);
	glm::mat4 modelViewMatrix = glm::lookAt(eyeDistance * glm::vec3(0.0f, 0.0f, 1.0f), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f)) * glm::rotate(glm::mat4(1.0f), angle, glm::vec3(0.0f, 1.0f, 0.0f));

	if (isOpenGL)
	{
		if (isTextured)
			glBindTexture(GL_TEXTURE_2D, texID);
		
		for (int i = 0; i < triangleVector.size(); i++)
			triangleVector[i].RenderOpenGL(modelViewMatrix, projectionMatrix, isTextured);  // drawing using OpenGL
		
		if (isTextured)
			glBindTexture(GL_TEXTURE_2D, 0);
	}
	else
	{
		for (int x = 0; x < WINDOW_WIDTH; x++) {
            for (int y = 0; y < WINDOW_HEIGHT; y++) {
                depth[y][x] = std::numeric_limits<float>::max();
            }
        }
		for (int i = 0; i < triangleVector.size(); i++) 
			triangleVector[i].RenderCPU(modelViewMatrix, projectionMatrix, color, depth, isTextured, texWidth, texHeight, texture, textureMode);

		glDrawPixels(WINDOW_WIDTH, WINDOW_HEIGHT, GL_RGB, GL_FLOAT, &color[0][0][0]);
		ClearFrameBuffer();
	}

	glFlush();
}

// Keyboard character callback function
// float minz = FLT_MAX;
// float maxz = 0.0;

void CharacterCallback(GLFWwindow* lWindow, unsigned int key)
{
	switch (key) 
	{
	case '0':  // all vertices in each triangle to same color to rand() floats color
		colorMode = 0;
		// loop thorugh all of the triangles in the TraingleVector and assign all vertices per triangle to the same random float value
		for (int i = 0; i < triangleVector.size(); i++) {
			triangleVector[i].setColorMode0();
		}
		break;
	case '1':
		colorMode = 1;
		for (int i = 0; i < triangleVector.size(); i++) {
			triangleVector[i].setColorMode1();
		}
		break;
	case '2':
		colorMode = 2;
		for (int i = 0; i < triangleVector.size(); i++) {
			triangleVector[i].setColorMode2();
		}
		break;
	case 'w':  // move closer to object
		eyeDistance *= (1 - 0.05);
		break;
	case 's':  // move further from object
		eyeDistance *= (1 + 0.05);
		break;
	case 'a':
		angle -= 0.03;
		break;
	case 'd':
		angle += 0.03;
		break;
	case ' ':
		isOpenGL = !isOpenGL;
		break;
	case 't':
	{
		if (!texture.empty())
			isTextured = !isTextured;
		break;
	}
	// The one you will be using for this part is the original image which is stored in the first element of the vector, i.e., texture[0]
	case 'n':
		glBindTexture(GL_TEXTURE_2D, texID);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glBindTexture(GL_TEXTURE_2D, 0);
		textureMode = 0;
		break;
	case 'l':
		glBindTexture(GL_TEXTURE_2D, texID);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glBindTexture(GL_TEXTURE_2D, 0);
		textureMode = 1;
		break;
	case 'm':
		glBindTexture(GL_TEXTURE_2D, texID);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glBindTexture(GL_TEXTURE_2D, 0);
		textureMode = 2;
		break;
	case 'q':
		glfwSetWindowShouldClose(window, GLFW_TRUE);
		break;
	default:
		break;
	}
}

// Create a vector of triangles. Considers the texture coordinates if they are available.
void CreateTriangleVector(std::vector<glm::vec3> &vertices, std::vector<glm::vec2>& texCoords)
{
	for (int i = 0; i < vertices.size() / 3; i++)
	{
		Triangle myTriangle;

		if (texCoords.empty())
			myTriangle = Triangle(vertices[i * 3 + 0], vertices[i * 3 + 1], vertices[i * 3 + 2]);
		else
			myTriangle = Triangle(vertices[i * 3 + 0], vertices[i * 3 + 1], vertices[i * 3 + 2], 
								texCoords[i * 3 + 0], texCoords[i * 3 + 1], texCoords[i * 3 + 2]);

		triangleVector.push_back(myTriangle);
	}
}

// Load the geometry and texture coordinates if available
void LoadModel(char* name, std::vector<glm::vec3> &vertices, std::vector<glm::vec2>& texCoords)
{
	// Taken from Shinjiro Sueda with slight modification
	std::string meshName(name);
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string errStr;
	bool rc = tinyobj::LoadObj(&attrib, &shapes, &materials, &errStr, meshName.c_str());
	if (!rc) {
		std::cerr << errStr << std::endl;
	}
	else {
		// Some OBJ files have different indices for vertex positions, normals,
		// and texture coordinates. For example, a cube corner vertex may have
		// three different normals. Here, we are going to duplicate all such
		// vertices.
		// Loop over shapes
		for (size_t s = 0; s < shapes.size(); s++) {
			// Loop over faces (polygons)
			size_t index_offset = 0;
			for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
				size_t fv = shapes[s].mesh.num_face_vertices[f];
				// Loop over vertices in the face.
				for (size_t v = 0; v < fv; v++) {
					// access to vertex
					tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
					vertices.push_back(glm::vec3(attrib.vertices[3 * idx.vertex_index + 0],
												 attrib.vertices[3 * idx.vertex_index + 1],
												 attrib.vertices[3 * idx.vertex_index + 2]));
					if (!attrib.texcoords.empty()) {
						texCoords.push_back(glm::vec2(attrib.texcoords[2 * idx.texcoord_index + 0],
							attrib.texcoords[2 * idx.texcoord_index + 1]));
					}
				}
				index_offset += fv;
			}
		}
	}
}

// Load texture and create downsampled versions of it for mipmapping
void LoadTexture(char* name)
{
	std::string texName(name);
	int c;
	stbi_set_flip_vertically_on_load(true);
	stbi_hdr_to_ldr_gamma(1.0f);
	float* image = stbi_loadf(texName.c_str(), &texWidth, &texHeight, &c, 0);
	
	if (!image)
		std::cerr << texName << " not found" << std::endl;
	else if (c != 3)
		std::cerr << texName << " must have 3 channels (RGB)" << std::endl;
	else if ((texWidth % 2) != 0 || (texHeight % 2) != 0)
		std::cerr << " must be a power of 2" << std::endl;
	else
		texture.push_back(image);

	int length = std::min(texWidth, texHeight);
	int numLevels = log2(length);
	
	float** downImages = new float* [numLevels];
	for (int i = 0; i < numLevels; i++)
		downImages[i] = new float[texWidth * texHeight * c];

	for (int i = 0; i < numLevels; i++)
	{
		int curWidth = texWidth / pow(2, i + 1);
		int curHeight = texHeight / pow(2, i + 1);
		float* temp = new float[curWidth * curHeight * c];
		stbir_resize_float(image, texWidth, texHeight, 0, temp, curWidth, curHeight, 0, c);
		stbir_resize_float(temp, texWidth / pow(2, i + 1), texHeight / pow(2, i + 1), 0, downImages[i], texWidth, texHeight, 0, c);
		texture.push_back(downImages[i]);
		stbi_image_free(temp);
	}


	if (!texture.empty())
	{
		glGenTextures(1, &texID);
		glBindTexture(GL_TEXTURE_2D, texID);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texWidth, texHeight, 0, GL_RGB, GL_FLOAT, texture[0]);
		glGenerateMipmap(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 0);
	}
		
}

std::string WindowTitle(std::string mainName)
{
	std::string hardwareName;
	if (isOpenGL)
		hardwareName = " - GPU";
	else
		hardwareName = " - CPU";

	std::string textureMethod;
	if (textureMode == 0)
		textureMethod = " - Nearest";
	else if (textureMode == 1)
		textureMethod = " - Bilinear";
	else if (textureMode == 2)
		textureMethod = " - Mipmap";

	std::string colorMethod;
	if (textureMode == 0)
		colorMethod = " - Mode 0";
	else if (textureMode == 1)
		colorMethod = " - Mode 1";
	else if (textureMode == 2)
		colorMethod = " - Mode 2";

	if (isTextured)
		return (mainName + hardwareName + std::string(" - Textured") + textureMethod);
	else
		return (mainName + hardwareName + std::string(" - Colored") + colorMethod);
}

void Init()
{
	glfwInit();
	glfwWindowHint(GLFW_COCOA_RETINA_FRAMEBUFFER, GL_FALSE);
	window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, WindowTitle(mainName).c_str(), NULL, NULL);
	glfwMakeContextCurrent(window);
	glfwSetCharCallback(window, CharacterCallback);
	glewExperimental = GL_TRUE;
	glewInit();
	glViewport(0, 0, WINDOW_WIDTH, WINDOW_HEIGHT);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glEnable(GL_DEPTH_TEST);

	ClearFrameBuffer();

	std::vector<glm::vec3> vertices;
	std::vector<glm::vec2> texCoords;

	char* userInput = new char[256];
	char* userInput2 = new char[256];
	cout << "Please Enter Object File Path for the Model: ";
	cin >> userInput;
	cout << "Please Enter Object File Path for the Texture: ";
	cin >> userInput2;
	// cout << "Please Enter Object File Path for the Texture: ";
	// cin2 >> userInput2;
	LoadModel("../resources/sphere.obj", vertices, texCoords);  // loads the model; reads the vertices of triangles from the desired object file, and writes them into the vertices vector
	
	if (!texCoords.empty())
	{
		LoadTexture("../resources/earth.jpg");  // loads the texture  -- earth
		if (texture.empty())
			isTextured = false;
	}
	else
		isTextured = false;
		
	CreateTriangleVector(vertices, texCoords); // creates an instance of the triangle class for each group of three vertices in the vertices vector, and pushes into triangleVector vector
	
}



int main()
{	
	// "../resources/bunny.obj"  -- loadmodel
	// "../resources/duck.obj"
	// "../resources/sphere.obj"
	//  "../resources/earth.jpg" -- loadtexture

	Init();
	while ( glfwWindowShouldClose(window) == 0) 
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		Display();
		glfwSwapBuffers(window);
		glfwPollEvents();
		glfwSetWindowTitle(window, WindowTitle(mainName).c_str());
	}

	glfwTerminate();
	return 0;
}