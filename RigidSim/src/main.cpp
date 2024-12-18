// ----------------------------------------------------------------------------
// main.cpp
//
//  Created on: 13 Dec 2020
//      Author: Kiwon Um
//        Mail: kiwon.um@telecom-paris.fr
//
// Description: Lab - Rigid Body Simulator (DO NOT distribute!)
//
// Copyright 2020-2024 Kiwon Um
//
// The copyright to the computer program(s) herein is the property of Kiwon Um,
// Telecom Paris, France. The program(s) may be used and/or copied only with
// the written permission of Kiwon Um or in accordance with the terms and
// conditions stipulated in the agreement/contract under which the program(s)
// have been supplied.
// ----------------------------------------------------------------------------

#define _USE_MATH_DEFINES

#include <glad/gl.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <glm/ext/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ios>
#include <vector>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>
#include <exception>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// #include "Error.h" // OpenGL 4.3 or later
#include "ShaderProgram.h"
#include "Camera.h"
#include "Mesh.h"

#include "RigidSolver.hpp"

// window parameters
GLFWwindow *g_window = nullptr;
int g_windowWidth = 1024;
int g_windowHeight = 768;

// pointer to the current camera model
std::shared_ptr<Camera> g_cam;

// camera control variables
float g_meshScale = 1.0; // to update based on the mesh size, so that navigation runs at scale
bool g_rotatingP = false;
bool g_panningP = false;
bool g_zoomingP = false;
double g_baseX = 0.0, g_baseY = 0.0;
glm::vec3 g_baseTrans(0.0);
glm::vec3 g_baseRot(0.0);

// timer
float g_appTimer = 0.0;
float g_appTimerLastClockTime;
bool g_appTimerStoppedP = true;

// textures
unsigned int g_availableTextureSlot = 0;

GLuint g_albedoTex;
unsigned int g_albedoTexOnGPU;

GLuint g_normalTex;
unsigned int g_normalTexOnGPU;

GLuint loadTextureFromFileToGPU(const std::string &filename)
{
  int width, height, numComponents;
  // Loading the image in CPU memory using stbd_image
  unsigned char *data = stbi_load(
    filename.c_str(),
    &width,
    &height,
    &numComponents, // 1 for a 8 bit greyscale image, 3 for 24bits RGB image, 4 for 32bits RGBA image
    0);

  // Create a texture in GPU memory
  GLuint texID;
  glGenTextures(1, &texID);
  glBindTexture(GL_TEXTURE_2D, texID);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  // Uploading the image data to GPU memory
  glTexImage2D(
    GL_TEXTURE_2D,
    0,
    (numComponents == 1 ? GL_RED : numComponents == 3 ? GL_RGB : GL_RGBA), // For greyscale images, we store them in the RED channel
    width,
    height,
    0,
    (numComponents == 1 ? GL_RED : numComponents == 3 ? GL_RGB : GL_RGBA), // For greyscale images, we store them in the RED channel
    GL_UNSIGNED_BYTE,
    data);

  // Generating mipmaps for filtered texture fetch
  glGenerateMipmap(GL_TEXTURE_2D);

  // Freeing the now useless CPU memory
  stbi_image_free(data);
  glBindTexture(GL_TEXTURE_2D, 0); // unbind the texture
  return texID;
}

struct Light {
  glm::vec3 position;
  glm::vec3 color;
  float intensity;
};


struct Scene {
  Light light;

  RigidSolver solver = RigidSolver(nullptr, Vec3f(0, -0.98, 0));
  std::shared_ptr<BodyAttributes> rigidAtt = nullptr;

  // meshes
  std::shared_ptr<Mesh> rigid = nullptr;
  std::shared_ptr<Mesh> plane = nullptr;

  // transformation matrices
  glm::mat4 rigidMat = glm::mat4(1.0);
  glm::mat4 planeMat = glm::mat4(1.0);
  glm::mat4 floorMat = glm::mat4(1.0);

  glm::vec3 scene_center = glm::vec3(0);
  float scene_radius = 1.f;

  // shaders to render the meshes
  std::shared_ptr<ShaderProgram> mainShader;

  // useful for debug
  bool saveScreenShot = false;
  int savedCnt = 0;

  void resetSim()
  {
    *rigidAtt = Box(.1f, .1f, .1f);
    solver.init(rigidAtt.get());
    rigidMat = glm::mat4(1.0);
  }

  void render()
  {
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glViewport(0, 0, g_windowWidth, g_windowHeight);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // clear the color and z buffers.

    //glDisable(GL_CULL_FACE);    // or
    glCullFace(GL_BACK);

    mainShader->use();

    // camera
    mainShader->set("camPos", g_cam->getPosition());
    mainShader->set("viewMat", g_cam->computeViewMatrix());
    mainShader->set("projMat", g_cam->computeProjectionMatrix());

    // light
    mainShader->set(std::string("lightSrc.position"), light.position);
    mainShader->set(std::string("lightSrc.color"), light.color);
    mainShader->set(std::string("lightSrc.intensity"), light.intensity);

    // back-wall
    mainShader->set("material.albedo", glm::vec3(0.29, 0.51, 0.82)); // default value if the texture was not loaded
    mainShader->set("material.albedoTexLoaded", 0);
    mainShader->set("material.normalTex", (int)g_normalTexOnGPU);
    mainShader->set("material.normalTexLoaded", 1);
    mainShader->set("modelMat", planeMat);
    mainShader->set("normMat", glm::mat3(glm::inverseTranspose(planeMat)));
    plane->render();

    // floor
    mainShader->set("material.albedo", glm::vec3(0.8, 0.8, 0.9));
    mainShader->set("material.albedoTexLoaded", 0);
    mainShader->set("material.normalTexLoaded", 0);
    mainShader->set("modelMat", floorMat);
    mainShader->set("normMat", glm::mat3(glm::inverseTranspose(floorMat)));
    plane->render();

    // rigid
    mainShader->set("material.albedo", glm::vec3(1, 0.71, 0.29));
    mainShader->set("material.albedoTex", (int)g_albedoTexOnGPU);
    mainShader->set("material.albedoTexLoaded", 1);
    mainShader->set("material.normalTexLoaded", 0);
    mainShader->set("modelMat", rigidMat);
    mainShader->set("normMat", glm::mat3(glm::inverseTranspose(rigidMat)));
    rigid->render();

    mainShader->stop();

    if(saveScreenShot) {
      std::stringstream fpath;
      fpath << "s" << std::setw(4) << std::setfill('0') << savedCnt++ << ".tga";

      std::cout << "Saving file " << fpath.str() << " ... " << std::flush;
      const short int w = g_windowWidth;
      const short int h = g_windowHeight;
      std::vector<int> buf(w*h*3, 0);
      glReadPixels(0, 0, w, h, GL_BGR, GL_UNSIGNED_BYTE, &(buf[0]));

      FILE *out = fopen(fpath.str().c_str(), "wb");
      short TGAhead[] = {0, 2, 0, 0, 0, 0, w, h, 24};
      fwrite(&TGAhead, sizeof(TGAhead), 1, out);
      fwrite(&(buf[0]), 3*w*h, 1, out);
      fclose(out);
      saveScreenShot = false;

      std::cout << "Done" << std::endl;
    }
  }
};

Scene g_scene;

void printHelp()
{
  std::cout <<
    "> Help:" << std::endl <<
    "    Mouse commands:" << std::endl <<
    "    * Left button: rotate camera" << std::endl <<
    "    * Middle button: zoom" << std::endl <<
    "    * Right button: pan camera" << std::endl <<
    "    Keyboard commands:" << std::endl <<
    "    * H: print this help" << std::endl <<
    "    * P: toggle simulation" << std::endl <<
    "    * R: reset simulation" << std::endl <<
    "    * S: save a screenshot" << std::endl <<
    "    * W: wireframe rendering" << std::endl <<
    "    * F: surface rendering" << std::endl <<
    "    * ESC: quit the program" << std::endl;
}

// Executed each time the window is resized. Adjust the aspect ratio and the rendering viewport to the current window.
void windowSizeCallback(GLFWwindow *window, int width, int height)
{
  g_windowWidth = width;
  g_windowHeight = height;
  g_cam->setAspectRatio(static_cast<float>(width)/static_cast<float>(height));
  glViewport(0, 0, (GLint)width, (GLint)height); // Dimension of the rendering region in the window
}

// Executed each time a key is entered.
void keyCallback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
  if(action == GLFW_PRESS && key == GLFW_KEY_H) {
    printHelp();
  } else if(action == GLFW_PRESS && key == GLFW_KEY_R) {
    g_scene.resetSim();
  } else if(action == GLFW_PRESS && key == GLFW_KEY_S) {
    g_scene.saveScreenShot = true;
  } else if(action == GLFW_PRESS && key == GLFW_KEY_P) {
    g_appTimerStoppedP = !g_appTimerStoppedP;
    if(!g_appTimerStoppedP)
      g_appTimerLastClockTime = static_cast<float>(glfwGetTime());
  } else if(action == GLFW_PRESS && key == GLFW_KEY_W) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  } else if(action == GLFW_PRESS && key == GLFW_KEY_F) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  } else if(action == GLFW_PRESS && key == GLFW_KEY_ESCAPE) {
    glfwSetWindowShouldClose(window, true); // Closes the application if the escape key is pressed
  }
}

// Called each time the mouse cursor moves
void cursorPosCallback(GLFWwindow *window, double xpos, double ypos)
{
  int width, height;
  glfwGetWindowSize(window, &width, &height);
  const float normalizer = static_cast<float>((width + height)/2);
  const float dx = static_cast<float>((g_baseX - xpos) / normalizer);
  const float dy = static_cast<float>((ypos - g_baseY) / normalizer);
  if(g_rotatingP) {
    const glm::vec3 dRot(-dy*M_PI, dx*M_PI, 0.0);
    g_cam->setRotation(g_baseRot + dRot);
  } else if(g_panningP) {
    g_cam->setPosition(g_baseTrans + g_meshScale*glm::vec3(dx, dy, 0.0));
  } else if(g_zoomingP) {
    g_cam->setPosition(g_baseTrans + g_meshScale*glm::vec3(0.0, 0.0, dy));
  }
}

// Called each time a mouse button is pressed
void mouseButtonCallback(GLFWwindow *window, int button, int action, int mods)
{
  if(button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
    if(!g_rotatingP) {
      g_rotatingP = true;
      glfwGetCursorPos(window, &g_baseX, &g_baseY);
      g_baseRot = g_cam->getRotation();
    }
  } else if(button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
    g_rotatingP = false;
  } else if(button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
    if(!g_panningP) {
      g_panningP = true;
      glfwGetCursorPos(window, &g_baseX, &g_baseY);
      g_baseTrans = g_cam->getPosition();
    }
  } else if(button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_RELEASE) {
    g_panningP = false;
  } else if(button == GLFW_MOUSE_BUTTON_MIDDLE && action == GLFW_PRESS) {
    if(!g_zoomingP) {
      g_zoomingP = true;
      glfwGetCursorPos(window, &g_baseX, &g_baseY);
      g_baseTrans = g_cam->getPosition();
    }
  } else if(button == GLFW_MOUSE_BUTTON_MIDDLE && action == GLFW_RELEASE) {
    g_zoomingP = false;
  }
}

void initGLFW()
{
  // Initialize GLFW, the library responsible for window management
  if(!glfwInit()) {
    std::cerr << "ERROR: Failed to init GLFW" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Before creating the window, set some option flags
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);

  // Create the window
  g_window = glfwCreateWindow(g_windowWidth, g_windowHeight, "A Simple Rigid Body Simulator", nullptr, nullptr);
  if(!g_window) {
    std::cerr << "ERROR: Failed to open window" << std::endl;
    glfwTerminate();
    std::exit(EXIT_FAILURE);
  }

  // Load the OpenGL context in the GLFW window using GLAD OpenGL wrangler
  glfwMakeContextCurrent(g_window);

  // not mandatory for all, but MacOS X
  glfwGetFramebufferSize(g_window, &g_windowWidth, &g_windowHeight);

  // Connect the callbacks for interactive control
  glfwSetWindowSizeCallback(g_window, windowSizeCallback);
  glfwSetKeyCallback(g_window, keyCallback);
  glfwSetCursorPosCallback(g_window, cursorPosCallback);
  glfwSetMouseButtonCallback(g_window, mouseButtonCallback);
}

void clear();
void exitOnCriticalError(const std::string &message)
{
  std::cerr << "> [Critical error]" << message << std::endl;
  std::cerr << "> [Clearing resources]" << std::endl;
  clear();
  std::cerr << "> [Exit]" << std::endl;
  std::exit(EXIT_FAILURE);
}

void initOpenGL()
{
  // Load extensions for modern OpenGL
  if(!gladLoadGL(glfwGetProcAddress))
    exitOnCriticalError("[Failed to initialize OpenGL context]");

  // supported from OpenGL 4.3
  // glEnable(GL_DEBUG_OUTPUT);    // Modern error callback functionality
  // glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS); // For recovering the line where the error occurs, set a debugger breakpoint in DebugMessageCallback
  // glDebugMessageCallback(debugMessageCallback, 0); // Specifies the function to call when an error message is generated.

  glCullFace(GL_BACK); // Specifies the faces to cull (here the ones pointing away from the camera)
  glEnable(GL_CULL_FACE); // Enables face culling (based on the orientation defined by the CW/CCW enumeration).
  glDepthFunc(GL_LESS);   // Specify the depth test for the z-buffer
  glEnable(GL_DEPTH_TEST);      // Enable the z-buffer test in the rasterization
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // specify the background color, used any time the framebuffer is cleared

  // Loads and compile the programmable shader pipeline
  try {
    g_scene.mainShader = ShaderProgram::genBasicShaderProgram("src/vertexShader.glsl", "src/fragmentShader.glsl");
    g_scene.mainShader->stop();
  } catch(std::exception &e) {
    exitOnCriticalError(std::string("[Error loading shader program]") + e.what());
  }
}

void initScene()
{
  // Init camera
  int width, height;
  glfwGetWindowSize(g_window, &width, &height);
  g_cam = std::make_shared<Camera>();
  g_cam->setAspectRatio(static_cast<float>(width)/static_cast<float>(height));

  // Load meshes in the scene
  {
    g_scene.rigid = std::make_shared<Mesh>();
    g_scene.rigid->addBox(.1f, .1f, .1f);
    g_scene.rigid->init();

    // for the solver
    g_scene.rigidAtt = std::make_shared<Box>(.1f, .1f, .1f);
    g_scene.solver.init(g_scene.rigidAtt.get());

    g_scene.plane = std::make_shared<Mesh>();
    g_scene.plane->addPlane();
    g_scene.plane->init();
    g_scene.planeMat = glm::translate(glm::mat4(1.0), glm::vec3(0, 0, -1.0));
    g_scene.floorMat = glm::translate(glm::mat4(1.0), glm::vec3(0, -1.0, 0))*
      glm::rotate(glm::mat4(1.0), (float)(-0.5f*M_PI), glm::vec3(1.0, 0.0, 0.0));
  }

  // Load textures
  {
    g_albedoTex = loadTextureFromFileToGPU("data/dice.png");
    g_normalTex = loadTextureFromFileToGPU("data/normal.png");
  }

  // Setup textures on the GPU
  {
    g_albedoTexOnGPU = g_availableTextureSlot;
    ++g_availableTextureSlot;
    glActiveTexture(GL_TEXTURE0 + g_albedoTexOnGPU);
    glBindTexture(GL_TEXTURE_2D, g_albedoTex);
  }
  {
    g_normalTexOnGPU = g_availableTextureSlot;
    ++g_availableTextureSlot;
    glActiveTexture(GL_TEXTURE0 + g_normalTexOnGPU);
    glBindTexture(GL_TEXTURE_2D, g_normalTex);
  }

  // Setup light
  g_scene.light.position = glm::vec3(0.0, 1.0, 1.0);
  g_scene.light.color = glm::vec3(1.0, 1.0, 1.0);
  g_scene.light.intensity = 1.0f;

  // Adjust the camera to the mesh
  g_meshScale = g_scene.scene_radius;
  g_cam->setPosition(g_scene.scene_center + glm::vec3(0.0, 0.0, 3.0*g_meshScale));
  g_cam->setNear(g_meshScale/100.f);
  g_cam->setFar(6.0*g_meshScale);
}

void init()
{
  initGLFW();                   // Windowing system
  initOpenGL();                 // OpenGL Context and shader pipeline
  initScene();                  // Actual scene to render
}

void clear()
{
  g_cam.reset();
  g_scene.rigid.reset();
  g_scene.plane.reset();
  g_scene.mainShader.reset();
  glfwDestroyWindow(g_window);
  glfwTerminate();
}

// The main rendering call
void render()
{
  g_scene.render();
}

// Update any accessible variable based on the current time
void update(const float currentTime)
{
  if(!g_appTimerStoppedP) {
    // Animate any entity of the program here
    const float dt = currentTime - g_appTimerLastClockTime;
    g_appTimerLastClockTime = currentTime;
    g_appTimer += dt;
    // <---- Update here what needs to be animated over time ---->

    g_scene.solver.step(std::min(dt, 0.017f)); // solve for the next step; avoid any chances of too large time step
    g_scene.rigidMat = g_scene.rigidAtt->worldMat(); // update position/orientation for rendering
  }
}

int main(int argc, char **argv)
{
  init();
  while(!glfwWindowShouldClose(g_window)) {
    update(static_cast<float>(glfwGetTime()));
    render();
    glfwSwapBuffers(g_window);
    glfwPollEvents();
  }
  clear();
  std::cout << " > Quit" << std::endl;
  return EXIT_SUCCESS;
}
