#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>

#include "../include/shader.hpp"

void processInput(GLFWwindow *window);
void render(
  Shader &shader,
  float vertices[],
  unsigned int VAO,
  unsigned int vertLen,
  unsigned int texture,
  unsigned int samples
);
unsigned int createFramebuffer(unsigned int *texture);

// settings
const unsigned int SCR_SIZE = 720;

int main()
{
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
    GLFWwindow* window = glfwCreateWindow(SCR_SIZE, SCR_SIZE, "Ray Tracer", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }    
    //
    // Build and compile our shader zprogram
    // Note: Paths to shader files should be relative to location of executable
    Shader shader("../shaders/vert.glsl", "../shaders/frag.glsl");
    Shader dispShader("../shaders/vert.glsl", "../shaders/disp.glsl");

    float vertices[] = {
      -1, -1,
      -1, +1,
      +1, +1,
      -1, -1,
      +1, +1,
      +1, -1,
    };

    unsigned int VBO, VAO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    // Position attribute
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // Create two needed framebuffers
    unsigned int fbTexture1;
    unsigned int fb1 = createFramebuffer(&fbTexture1);
    unsigned int fbTexture2;
    unsigned int fb2 = createFramebuffer(&fbTexture2);

    unsigned int samples = 0;

    // Render first pass on fb1
    glBindFramebuffer(GL_FRAMEBUFFER, fb1);
    render(shader, vertices, VAO, sizeof(vertices), fbTexture1, samples);
    samples += 1;

    while (!glfwWindowShouldClose(window)) {
      // input
      // -----
      processInput(window);

      // Render to fb2
      glBindFramebuffer(GL_FRAMEBUFFER, fb2);
      render(shader, vertices, VAO, sizeof(vertices), fbTexture1, samples);
      samples += 1;

      // Render to fb1
      glBindFramebuffer(GL_FRAMEBUFFER, fb1);
      render(shader, vertices, VAO, sizeof(vertices), fbTexture2, samples);
      samples += 1;

      // Render to screen
      glBindFramebuffer(GL_FRAMEBUFFER, 0);
      glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
      glClear(GL_COLOR_BUFFER_BIT);

      dispShader.use();
      unsigned int ID = dispShader.ID;

      // Render from fb2 texture
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, fbTexture2);
      dispShader.setFloat("exposure", 2.2);
      dispShader.setInt("screenTexture", 0);
      glUniform2f(glGetUniformLocation(ID, "resolution"), SCR_SIZE, SCR_SIZE);
      glBindVertexArray(VAO);
      glDrawArrays(GL_TRIANGLES, 0, sizeof(vertices) / 3);

      // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
      // -------------------------------------------------------------------------------
      glfwSwapBuffers(window);
      glfwPollEvents();
    }

    // Deallocate all resources once they've outlived their purpose
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteFramebuffers(1, &fb1);
    glDeleteFramebuffers(1, &fb2);
    glDeleteTextures(1, &fbTexture1);
    glDeleteTextures(1, &fbTexture2);

    std::cout << "Total samples: " << samples;

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    return 0;
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow *window)
{
  if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
    glfwSetWindowShouldClose(window, true);
  }
}

void render(
  Shader &shader,
  float vertices[],
  unsigned int VAO,
  unsigned int vertLen,
  unsigned int texture,
  unsigned int samples
) {
  // render
  // ------
  glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);

  // render the shader
  shader.use();

  // Set uniforms
  unsigned int ID = shader.ID;

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, texture);
  shader.setInt("prevFrame", 0);
  shader.setFloat("acc", float(samples) / float(samples+1));

  shader.setInt("samples", samples);
  shader.setInt("numBounces", 8);
  shader.setFloat("time", glfwGetTime());
  shader.setInt("seedInit", rand());

  // Light properties
  shader.setFloat("intensity", 1.0);

  // Camera properties
  shader.setFloat("focalDistance", 2);
  glUniform2f(glGetUniformLocation(ID, "resolution"), SCR_SIZE, SCR_SIZE);

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
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, SCR_SIZE, SCR_SIZE, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

  // Attach texture to framebuffer
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, fb, 0);

  // Check if framebuffer is ready to be written to
  if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
    std::cout << "ERROR::FRAMEBUFFER:: Framebuffer is not complete!" << std::endl;
  }

  return fb;
}
