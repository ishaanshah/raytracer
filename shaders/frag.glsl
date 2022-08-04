#version 330

out vec4 fragColor;

uniform int samples;

// Light properties
uniform float intensity;

// Camera properties
uniform vec2  resolution;
uniform float focalDistance;
uniform float exposure;

// Constants
const float pi = 3.14159265;
const float inf = 99999999.0;
const float eps = 1e-6;

// Primitives
struct Ray
{
  vec3 origin;
  vec3 dir;
};

struct Sphere {
  vec3 center;
  vec3 radius;
};

struct Cuboid {
  vec3 size;
  vec3 center;
};

struct Rectangle {
  vec3 position; // Center of rectangle
  vec3 normal;
  float height;
  float width;
};

struct Triangle {
  vec3 vert1;
  vec3 vert2;
  vec3 vert3;
  vec3 normal;
};

struct Material {
  int type; // 0 for reflective, 1 for refractive, 2 for emmisive
  float roughness;
  vec3 color;
};

// Wals
vec4[5] walls = vec4[5](
  vec4(0, 1, 0, -5),
  vec4(0, -1, 0, -5),
  vec4(1, 0, 0, -5),
  vec4(-1, 0, 0, -5),
  vec4(0, 0, 1, -5)
);
Material[5] wallsMat = Material[5](
  Material(0, 1, vec3(0.5, 0.0, 0.5)),
  Material(0, 1, vec3(0.0, 0.5, 0.5)),
  Material(0, 1, vec3(0.0, 1.0, 0.0)),
  Material(0, 1, vec3(1.0, 0.0, 0.0)),
  Material(0, 1, vec3(0.5, 0.5, 0.5))
);

// Light
vec4 light = vec4(0, -1, 0, -5),

// Intersection
bool intersectPlane(Ray ray, vec4 plane, out float t)
{
  t = -dot(plane, vec4(ray.origin, 1.0)) / dot(plane.xyz, ray.dir);
  return t > 0.0;
}

// Camera is placed at (0, 0, 0) looking at (0, 0, -1)
Ray genRay()
{
  Ray ray;

  vec2 xy = 2.0*gl_FragCoord.xy/resolution - vec2(1.0);

  // TODO: evaluate the use of this
  ray.dir = normalize(vec3(xy, -focalDistance));
  // float ft = focalDistance/ray.dir.z;
  float ft = -1;
  vec3 pFocus = ray.dir*ft;

  // TODO: Antialiasing by subsampling pixel
  ray.origin = vec3(0);
  ray.dir = normalize(pFocus - ray.origin);

  return ray;
}

bool intersect(Ray ray, out float t, out Material mat) {
  bool intersect = false;
  t = inf; 
  float tmp = 0;

  for (int i = 0; i < 5; i += 1) {
    if (intersectPlane(ray, walls[i], tmp)) {
      intersect = true;
      if (tmp < t) {
        t = tmp;
        mat = wallsMat[i];
      }
    }
  }

  return intersect;
}

vec3 rotation_y(vec3 v, float a)
{
    vec3 r;
    r.x =  v.x*cos(a) + v.z*sin(a);
    r.y =  v.y;
    r.z = -v.x*sin(a) + v.z*cos(a);
    return r;
}

vec3 rayTrace(Ray ray) {
  Material mat;
  float t;
  if (intersect(ray, t, mat)) {
    return mat.color;
  }
  return vec3(0);
}

void main()
{
  Ray ray = genRay();

  fragColor = vec4(rayTrace(ray), 1);
}
