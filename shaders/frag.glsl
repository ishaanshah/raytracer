#version 330

#define DEBUG_NORMALS

out vec4 fragColor;

// Render settings
uniform int samples;
uniform int numBounces;
uniform float time;

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

int seed = 0;
int flatIdx = 0;

// Primitives
struct Ray
{
  vec3 origin;
  vec3 dir;
};

struct Sphere {
  vec3 center;
  float radius;
};

struct Cuboid {
  vec3 size;
  vec3 center;
};

struct Rectangle {
  vec4 plane;
  vec3 center;
  vec3 dirX;
  vec3 dirY;
  float lenX;
  float lenY;
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

// Walls
vec4[5] walls = vec4[5](
  vec4(0, 1, 0, 5),
  vec4(0, -1, 0, 5),
  vec4(1, 0, 0, 5),
  vec4(-1, 0, 0, 5),
  vec4(0, 0, 1, 5)
);
Material[5] wallsMat = Material[5](
  Material(0, 1, vec3(0.1, 0.1, 0.1)),
  Material(0, 1, vec3(0.0, 0.0, 0.5)),
  Material(0, 1, vec3(0.0, 1.0, 0.0)),
  Material(0, 1, vec3(1.0, 0.0, 0.0)),
  Material(0, 1, vec3(0.5, 0.5, 0.5))
);

// Sphere
Sphere sphere = Sphere(vec3(1.5, -3.5, 3), 1.5);
Material sphereMat = Material(1, 0, vec3(0.3, 0.1, 0.4));

// Cuboid
Cuboid cuboid = Cuboid(vec3(2), vec3(-1.5, -4, 1));
Material cuboidMat = Material(0, 0.5, vec3(0.5, 0.5, 0.0));

// Light
Rectangle light = Rectangle(
  vec4(0, -1, 0, 5-eps),
  vec3(0, 5-eps, 0),
  vec3(1, 0, 0),
  vec3(0, 0, 1),
  5, 5
);
Material lightMat = Material(2, 0, vec3(1));

// Camera is placed at (0, 0, 0) looking at (0, 0, -1)
Ray genRay() {
  Ray ray;

  vec2 xy = 2.0*gl_FragCoord.xy/resolution - vec2(1.0);

  // TODO: Antialiasing by subsampling pixel
  ray.origin = vec3(0, 0, 15);
  ray.dir = normalize(vec3(xy, 15-focalDistance) - ray.origin);

  return ray;
}

// Random number generator
void encryptTea(inout uvec2 arg) {
	uvec4 key = uvec4(0xa341316c, 0xc8013ea4, 0xad90777d, 0x7e95761e);
	uint v0 = arg[0], v1 = arg[1];
	uint sum = 0u;
	uint delta = 0x9e3779b9u;

	for(int i = 0; i < 32; i++) {
		sum += delta;
		v0 += ((v1 << 4) + key[0]) ^ (v1 + sum) ^ ((v1 >> 5) + key[1]);
		v1 += ((v0 << 4) + key[2]) ^ (v0 + sum) ^ ((v0 >> 5) + key[3]);
	}
	arg[0] = v0;
	arg[1] = v1;
}

vec2 getRandom() {
  uvec2 arg = uvec2(flatIdx, seed++);
  encryptTea(arg);
  return fract(vec2(arg) / vec2(0xffffffffu));
}

// Intersection
bool intersectPlane(Ray ray, vec4 plane, out float t, out vec3 normal)
{
  t = -dot(plane, vec4(ray.origin, 1.0)) / dot(plane.xyz, ray.dir);
  normal = plane.xyz;
  return t > 0.0;
}

bool intersectRectangle(Ray ray, Rectangle rect, out float t, out vec3 normal) {
  // Check if on correct side of plane
  if (!intersectPlane(ray, rect.plane, t, normal)) {
    return false;
  }
  // Check if lies inside rectangle
  vec3 newPos = ray.origin + ray.dir*t;
  newPos = newPos - rect.center;
  float projX = dot(newPos, rect.dirX);
  float projY = dot(newPos, rect.dirY);
  if (abs(projX) < rect.lenX / 2 && abs(projY) < rect.lenY / 2) {
    return true;
  }
  return false;
}

bool intersectCuboid(Ray ray, Cuboid cuboid, out float t, out vec3 normal) {
  // Check cuboid 
  return true;
}

bool intersectSphere(Ray ray, Sphere sphere, out float t, out vec3 normal) {
  vec3 oc = ray.origin - sphere.center;
  float b = dot(oc, ray.dir);
  float c = dot(oc, oc) - sphere.radius*sphere.radius;
  float disc = b*b - c;
  if(disc < 0.0) {
    return false;
  }
  disc = sqrt(disc);
  t = min(-b-disc, -b+disc);
  normal = normalize(ray.origin + ray.dir*t - sphere.center);
  return true;
}

// TODO: BVH
bool intersect(Ray ray, bool shadow, out float t, out Material mat, out vec3 normal) {
  bool intersect = false;
  t = inf; 
  float tmpT = 0;
  vec3 tmpNorm;

  if (intersectSphere(ray, sphere, tmpT, tmpNorm)) {
    intersect = true;
    if (tmpT < t) {
      t = tmpT;
      mat = sphereMat;
      normal = tmpNorm;
    }
  }

  // Light doesn't cast shadow
  if (!shadow) {
    if (intersectRectangle(ray, light, tmpT, tmpNorm)) {
      intersect = true;
      if (tmpT < t) {
        t = tmpT;
        mat = lightMat;
        normal = tmpNorm;
      }
    }
  }

  // Walls don't cast shadows
  if (!shadow) {
    for (int i = 0; i < 5; i += 1) {
      if (intersectPlane(ray, walls[i], tmpT, tmpNorm)) {
        intersect = true;
        if (tmpT < t) {
          t = tmpT;
          mat = wallsMat[i];
          normal = tmpNorm;
        }
      }
    }
  }

  return intersect;
}

// Reflective GGX
float D(vec3 N, vec3 H, float roughness) {
  float a2 = roughness * roughness;
  float NdotH = dot(N, H);
  float denom = ((NdotH*NdotH)*(a2*a2-1) + 1);
  return a2 / (pi * denom * denom);
}

float G(vec3 N, vec3 V, vec3 L, float roughness) {
  float NdotV = max(dot(N, V), 0);
  float NdotL = max(dot(N, L), 0);
  float k = roughness * roughness / 2;
  float G_V = NdotV / (NdotV*(1.0 - k) + k);
  float G_L = NdotL / (NdotL*(1.0 - k) + k);

  return G_V * G_L;
}

vec3 F(vec3 look, vec3 H, vec3 F0) {
  float LdotH = dot(look, H);
  return F0 + (1.0 - F0) * pow(1.0 - LdotH, 5.0);
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
  vec3 normal;
  vec3 Le = vec3(0);
  float tp = 1;

  #ifdef DEBUG_LOCATION
  if (intersect(ray, false, t, mat, normal)) {
    return mat.color;
  }
  #endif

  #ifdef DEBUG_NORMALS
  if (intersect(ray, false, t, mat, normal)) {
    return normalize(normal + vec3(1));
  }
  #endif

  for (int i = 0; i < numBounces; i++) {
    if (!intersect(ray, false, t, mat, normal)) {
      break;
    }
  }

  return tp*Le;
}

void main() {
  flatIdx = int(gl_FragCoord.x * resolution.x + gl_FragCoord.y);
  Ray ray = genRay();

  fragColor = vec4(rayTrace(ray), 1);
}
