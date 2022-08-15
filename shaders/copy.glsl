// Copy to intermediate buffer
#version 330 core

out vec4 fragColor;

uniform sampler2D fb;

uniform vec2 resolution;

void main() { fragColor = texture(fb, gl_FragCoord.xy / resolution.xy); }
