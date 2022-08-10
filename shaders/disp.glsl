#version 330 core

#define DEBUG_COLOR
out vec4 fragColor;
uniform float exposure;

uniform sampler2D screenTexture;
uniform vec2 resolution;

void main() {
  vec3 op = texture(screenTexture, gl_FragCoord.xy / resolution.xy).xyz;
#ifdef DEBUG_COLOR
  fragColor = vec4(op, 1);
#else
  fragColor = vec4(pow(exposure * op, vec3(1.0 / 2.2)), 1);
#endif
}
