// For debugging progressive path tracer
#version 330 core

out vec4 fragColor;
  
uniform sampler2D fb1;
uniform sampler2D fb2;

uniform vec2 resolution;

void main() { 
  vec3 op1 = texture(fb1, gl_FragCoord.xy / resolution.xy).xyz;
  vec3 op2 = texture(fb2, gl_FragCoord.xy / resolution.xy).xyz;
  fragColor = op2-op1;
}
