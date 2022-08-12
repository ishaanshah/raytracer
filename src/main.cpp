// #define DEBUG_SINGLE
#include <GLFW/glfw3.h>
#include <glad/glad.h>
#include <iostream>
#include <vector>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include "../include/shader.hpp"

void processInput(GLFWwindow *window);
void render(Shader &shader, float vertices[], unsigned int VAO,
            unsigned int vertLen, unsigned int samples);
unsigned int createFramebuffer(unsigned int *texture);
void saveImage(char *filepath, GLFWwindow *w);

// settings
unsigned int SCR_SIZE = 0;
unsigned int SAMPLES = 0;
char *PATH = NULL;

int main(int argc, char *argv[]) {
  // Validate command line arguments
  if (argc != 4) {
    std::cout << "USAGE::raytracer <NUM_SAMPLES> <IMAGE_SIZE> <IMAGE_PATH>"
              << std::endl;
    return 1;
  }
#ifdef DEBUG_SINGLE
  SAMPLES = 1e9;
#else
  SAMPLES = atoi(argv[1]);
#endif
  SCR_SIZE = atoi(argv[2]);
  PATH = argv[3];

  // glfw: initialize and configure
  // ------------------------------
  glfwInit();
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  srand(time(NULL));
#ifdef __APPLE__
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
  // glfw window creation
  // --------------------
  GLFWwindow *window =
      glfwCreateWindow(SCR_SIZE, SCR_SIZE, "Ray Tracer", NULL, NULL);
  if (window == NULL) {
    std::cout << "Failed to create GLFW window" << std::endl;
    glfwTerminate();
    return -1;
  }
  glfwMakeContextCurrent(window);

  // glad: load all OpenGL function pointers
  // ---------------------------------------
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
  //
  // Build and compile our shader zprogram
  // Note: Paths to shader files should be relative to location of executable
  Shader shader("../shaders/vert.glsl", "../shaders/frag.glsl");
  Shader copyShader("../shaders/vert.glsl", "../shaders/copy.glsl");
  Shader dispShader("../shaders/vert.glsl", "../shaders/disp.glsl");

  float vertices[] = {
      -1, -1, -1, +1, +1, +1, -1, -1, +1, +1, +1, -1,
  };

  unsigned int VBO, VAO;
  glGenVertexArrays(1, &VAO);
  glGenBuffers(1, &VBO);

  glBindVertexArray(VAO);

  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

  // Position attribute
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void *)0);
  glEnableVertexAttribArray(0);

  // Create two needed framebuffers
  unsigned int fbTexture1;
  unsigned int fb1 = createFramebuffer(&fbTexture1);
  unsigned int fbTexture2;
  unsigned int fb2 = createFramebuffer(&fbTexture2);

  unsigned int samples = 0;

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, fbTexture1);
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D, fbTexture2);

  // Store start time
  double t0 = glfwGetTime();
  while (!glfwWindowShouldClose(window) && samples <= SAMPLES) {
    // input
    // -----
    processInput(window);

    // Render pass on fb1
    glBindFramebuffer(GL_FRAMEBUFFER, fb1);
    render(shader, vertices, VAO, sizeof(vertices), samples);

    // Copy to fb2
    glBindFramebuffer(GL_FRAMEBUFFER, fb2);
    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    unsigned int ID = copyShader.ID;
    copyShader.use();
    copyShader.setInt("fb", 0);
    glUniform2f(glGetUniformLocation(ID, "resolution"), SCR_SIZE, SCR_SIZE);
    glBindVertexArray(VAO);
    glDrawArrays(GL_TRIANGLES, 0, sizeof(vertices) / 3);
    samples += 1;

    // Render to screen
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    dispShader.use();
    ID = dispShader.ID;

    // Render from fb2 texture
    dispShader.setFloat("exposure", 5);
    dispShader.setInt("screenTexture", 0);

    glUniform2f(glGetUniformLocation(ID, "resolution"), SCR_SIZE, SCR_SIZE);
    glBindVertexArray(VAO);
    glDrawArrays(GL_TRIANGLES, 0, sizeof(vertices) / 3);

    // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved
    // etc.)
    // -------------------------------------------------------------------------------
    glfwSwapBuffers(window);
    glfwPollEvents();

    std::cout << "Progress: " << samples << "/" << SAMPLES << " samples" << '\r';
  }

  // Write to file
  std::cout << "INFO::Output image written to " << PATH << std::endl;
  saveImage(PATH, window);

  std::cout << "INFO::Time taken: " << glfwGetTime() - t0 << "s"
            << std::endl;

  // Deallocate all resources once they've outlived their purpose
  glDeleteVertexArrays(1, &VAO);
  glDeleteBuffers(1, &VBO);
  glDeleteFramebuffers(1, &fb1);
  glDeleteFramebuffers(1, &fb2);
  glDeleteTextures(1, &fbTexture1);
  glDeleteTextures(1, &fbTexture2);

  // glfw: terminate, clearing all previously allocated GLFW resources.
  // ------------------------------------------------------------------
  glfwTerminate();
  return 0;
}

// process all input: query GLFW whether relevant keys are pressed/released this
// frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow *window) {
  if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
    glfwSetWindowShouldClose(window, true);
  }
}

void render(Shader &shader, float vertices[], unsigned int VAO,
            unsigned int vertLen, unsigned int samples) {
  // render
  // ------
  glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);

  // render the shader
  shader.use();

  // Set uniforms
  unsigned int ID = shader.ID;

  shader.setInt("prevFrame", 1);
#ifndef DEBUG_SINGLE
  shader.setFloat("acc", float(samples) / float(samples + 1));
#else
  shader.setFloat("acc", 0);
#endif

  shader.setInt("samples", samples);
  shader.setInt("numBounces", 8);
  shader.setFloat("time", glfwGetTime());
  shader.setInt("seedInit", rand());

  // Light properties
  shader.setFloat("intensity", 1.0);

  // Camera properties
  shader.setFloat("focalDistance", 2);
  glUniform2f(glGetUniformLocation(ID, "resolution"), SCR_SIZE, SCR_SIZE);

  // Checkerboard
  shader.setFloat("checkerboard", 2);

  glBindVertexArray(VAO);
  glDrawArrays(GL_TRIANGLES, 0, vertLen / 3);
}

unsigned int createFramebuffer(unsigned int *texture) {
  // Create a framebuffer to write output to
  unsigned int fb;
  glGenFramebuffers(1, &fb);
  glBindFramebuffer(GL_FRAMEBUFFER, fb);

  // Create a texture to write to
  glGenTextures(1, texture);
  glBindTexture(GL_TEXTURE_2D, *texture);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, SCR_SIZE, SCR_SIZE, 0, GL_RGBA,
               GL_UNSIGNED_BYTE, NULL);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glBindTexture(GL_TEXTURE_2D, 0);

  // Attach texture to framebuffer
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
                         *texture, 0);

  // Check if framebuffer is ready to be written to
  if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
    std::cout << "ERROR::FRAMEBUFFER:: Framebuffer is not complete!"
              << std::endl;
  }

  return fb;
}

void saveImage(char *filepath, GLFWwindow *w) {
  int width, height;
  glfwGetFramebufferSize(w, &width, &height);
  GLsizei nrChannels = 3;
  GLsizei stride = nrChannels * width;
  stride += (stride % 4) ? (4 - stride % 4) : 0;
  GLsizei bufferSize = stride * height;
  std::vector<char> buffer(bufferSize);
  glPixelStorei(GL_PACK_ALIGNMENT, 4);
  glReadBuffer(GL_FRONT);
  glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, buffer.data());
  stbi_flip_vertically_on_write(true);
  stbi_write_png(filepath, width, height, nrChannels, buffer.data(), stride);
}
