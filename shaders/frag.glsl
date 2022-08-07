#version 330

// #define DEBUG_LOCATION
// #define DEBUG_NORMALS
#define USE_NEE
#define USE_VNDF

out vec4 fragColor;

// Previous frame
uniform sampler2D prevFrame;
uniform float acc;

// Render settings
uniform int numBounces;
uniform float time;
uniform int seedInit;

// Light properties
uniform float intensity;

// Camera properties
uniform vec2 resolution;
uniform float focalDistance;

// Constants
const float pi = 3.14159265;
const float inf = 99999999.0;
const float eps = 1e-6;

int seed = 0;
int flatIdx = 0;

// Primitives
struct Ray {
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
  // normal = cross(dirx, diry)
  // d = dot(normal, center)
  vec4 plane;  // vec4(normal, d)
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
  int type;  // 0 for reflective, 1 for refractive, 2 for emmisive
  float roughness;
  vec3 color;
};

// Walls
vec4[5] walls = vec4[5](vec4(0, 1, 0, 5), vec4(0, -1, 0, 5), vec4(1, 0, 0, 5),
                        vec4(-1, 0, 0, 5), vec4(0, 0, 1, 5));
Material[5] wallsMat = Material[5](Material(0, 0.2, vec3(0.5, 0.5, 0.5)),
                                   Material(0, 0.2, vec3(0.5, 0.5, 0.5)),
                                   Material(0, 0.2, vec3(0.0, 1.0, 0.0)),
                                   Material(0, 0.2, vec3(1.0, 0.0, 0.0)),
                                   Material(0, 0.2, vec3(0.5, 0.5, 0.5)));

// Sphere
Sphere sphere = Sphere(vec3(1.5, -3.5, 3), 1.5);
Material sphereMat = Material(1, 1.0, vec3(1.0, 1.0, 1.0));

// Cuboid
Cuboid cuboid = Cuboid(vec3(2), vec3(-1.5, -4, 1));
Material cuboidMat = Material(0, 0.5, vec3(0.5, 1.0, 0.0));

// Light
Rectangle lightRect = Rectangle(vec4(0, -1, 0, 5 - eps), vec3(0, 5 - eps, 0),
                                vec3(1, 0, 0), vec3(0, 0, 1),
                                // 100, 100
                                5, 5);
Material lightMat = Material(2, 0, vec3(3));

// Random number generator
void encryptTea(inout uvec2 arg) {
  uvec4 key = uvec4(0xa341316c, 0xc8013ea4, 0xad90777d, 0x7e95761e);
  uint v0 = arg[0], v1 = arg[1];
  uint sum = 0u;
  uint delta = 0x9e3779b9u;

  for (int i = 0; i < 32; i++) {
    sum += delta;
    v0 += ((v1 << 4) + key[0]) ^ (v1 + sum) ^ ((v1 >> 5) + key[1]);
    v1 += ((v0 << 4) + key[2]) ^ (v0 + sum) ^ ((v0 >> 5) + key[3]);
  }
  arg[0] = v0;
  arg[1] = v1;
}

float clampDot(vec3 a, vec3 b, bool zero) {
  return max(dot(a, b), zero ? 0 : eps);
}

vec2 getRandom() {
  uvec2 arg = uvec2(flatIdx, seed++);
  encryptTea(arg);
  return fract(vec2(arg) / vec2(0xffffffffu));
}

mat3 constructONBFrisvad(vec3 normal) {
  mat3 ret;
  ret[1] = normal;
  if (normal.z < -0.999805696) {
    ret[0] = vec3(0.0, -1.0, 0.0);
    ret[2] = vec3(-1.0, 0.0, 0.0);
  } else {
    float a = 1.0 / (1.0 + normal.z);
    float b = -normal.x * normal.y * a;
    ret[0] = vec3(1.0 - normal.x * normal.x * a, b, -normal.x);
    ret[2] = vec3(b, 1.0 - normal.y * normal.y * a, -normal.y);
  }
  return ret;
}

// Camera is placed at (0, 0, 0) looking at (0, 0, -1)
Ray genRay(vec2 rng) {
  Ray ray;

  vec2 xy = 2.0 * gl_FragCoord.xy / resolution - vec2(1.0);

  ray.origin = vec3(0, 0, 15);
  ray.dir = normalize(
      vec3(xy + rng.x * dFdx(xy) + rng.y * dFdy(xy), 15 - focalDistance) -
      ray.origin);

  return ray;
}

// Intersection
bool intersectPlane(Ray ray, vec4 plane, out float t, out vec3 normal) {
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
  vec3 newPos = ray.origin + ray.dir * t;
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
  float c = dot(oc, oc) - sphere.radius * sphere.radius;
  float disc = b * b - c;
  if (disc < 0.0) {
    return false;
  }
  disc = sqrt(disc);
  t = min(-b - disc, -b + disc);
  normal = normalize(ray.origin + ray.dir * t - sphere.center);
  return true;
}

// TODO: BVH
bool intersect(Ray ray, bool shadow, out float t, out Material mat,
               out vec3 normal) {
  bool intersect = false;
  t = inf;
  float tmpT;
  vec3 tmpNorm;

  // if (intersectSphere(ray, sphere, tmpT, tmpNorm)) {
  //   intersect = true;
  //   if (tmpT < t) {
  //     t = tmpT;
  //     mat = sphereMat;
  //     normal = tmpNorm;
  //   }
  // }

  if (intersectRectangle(ray, lightRect, tmpT, tmpNorm)) {
    intersect = true;
    if (tmpT < t) {
      t = tmpT;
      mat = lightMat;
      normal = tmpNorm;
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
  float NdotH = clampDot(N, H, true);
  float denom = ((NdotH * NdotH) * (a2 * a2 - 1) + 1);
  return a2 / (pi * denom * denom);
}

float G(vec3 N, vec3 V, vec3 L, float roughness) {
  float NdotV = clampDot(N, V, false);
  float NdotL = clampDot(N, L, false);
  float k = roughness * roughness / 2;
  float G_V = NdotV / (NdotV * (1.0 - k) + k);
  float G_L = NdotL / (NdotL * (1.0 - k) + k);

  return G_V * G_L;
}

vec3 F(vec3 L, vec3 H, vec3 F0) {
  float LdotH = clampDot(L, H, true);
  return F0 + (1.0 - F0) * pow(1.0 - LdotH, 5.0);
}

#ifdef USE_VNDF
// Adapted from https://jcgt.org/published/0007/04/01/
vec3 sampleVndf(vec3 V, float roughness, vec2 rng, mat3 onb) {
  V = vec3(dot(V, onb[0]), dot(V, onb[2]), dot(V, onb[1]));

  // Transforming the view direction to the hemisphere configuration
  V = normalize(vec3(roughness * V.x, roughness * V.y, V.z));

  // Orthonormal basis (with special case if cross product is zero)
  float lensq = V.x * V.x + V.y * V.y;
  vec3 T1 =
      lensq > 0. ? vec3(-V.y, V.x, 0) * inversesqrt(lensq) : vec3(1, 0, 0);
  vec3 T2 = cross(V, T1);

  // Parameterization of the projected area
  float r = sqrt(rng.x);
  float phi = 2.0 * pi * rng.y;
  float t1 = r * cos(phi);
  float t2 = r * sin(phi);
  float s = 0.5 * (1.0 + V.z);
  t2 = (1.0 - s) * sqrt(1.0 - t1 * t1) + s * t2;

  // Reprojection onto hemisphere
  vec3 Nh = t1 * T1 + t2 * T2 + sqrt(max(0.0, 1.0 - t1 * t1 - t2 * t2)) * V;
  // Transforming the normal back to the ellipsoid configuration
  vec3 Ne = normalize(vec3(roughness * Nh.x, max(0.0, Nh.z), roughness * Nh.y));
  return normalize(onb * Ne);
}

float sampleVndfPdf(vec3 V, vec3 H, float D) {
  float VdotH = clampDot(V, H, false);
  return VdotH > 0 ? D / (4 * VdotH) : 0;
}

#else

vec3 sampleNdf(float roughness, vec2 rng, mat3 onb) {
  // GGX NDF sampling
  float a2 = roughness * roughness;
  float cosThetaH = sqrt(max(0.0, (1.0 - rng.x) / ((a2 - 1.0) * rng.x + 1.0)));
  float sinThetaH = sqrt(max(0.0, 1.0 - cosThetaH * cosThetaH));
  float phiH = rng.y * pi * 2.0;

  // Get our GGX NDF sample (i.e., the half vector)
  vec3 H = vec3(sinThetaH * cos(phiH), cosThetaH, sinThetaH * sin(phiH));
  return normalize(onb * H);
}

float sampleNdfPdf(vec3 L, vec3 H, vec3 N, float D) {
  return D * clampDot(N, H, true) / (4 * clampDot(H, L, false));
}
#endif

// Light sampling
vec3 sampleLight(Rectangle light, vec2 rng) {
  float x = (rng.x - 0.5) * light.lenX;
  float y = (rng.y - 0.5) * light.lenY;
  vec3 pos = (light.dirX * x + light.dirY * y) + light.center;
  return pos;
}

float sampleLightPdf(Rectangle light) { return 1 / (light.lenX * light.lenY); }

// Convert from area measure to angle measure
float pdfA2W(float pdf, float dist2, float cos_theta) {
  float abs_cos_theta = abs(cos_theta);
  if (abs_cos_theta < 1e-8) {
    return 0.0;
  }

  return pdf * dist2 / abs_cos_theta;
}

vec3 rotationY(vec3 v, float a) {
  vec3 r;
  r.x = v.x * cos(a) + v.z * sin(a);
  r.y = v.y;
  r.z = -v.x * sin(a) + v.z * cos(a);
  return r;
}

vec3 rayTrace(Ray ray) {
  Material mat;
  float t;
  vec3 normal;
  vec3 contrib = vec3(0);
  vec3 tp = vec3(1);
  vec3 pos = ray.origin;

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

  bool hit = intersect(ray, false, t, mat, normal);
  if (!hit) {
    return vec3(0);
  }

  if (mat.type == 2) {
    return mat.color;
  }

  for (int i = 0; i < numBounces; i++) {
    mat3 onb = constructONBFrisvad(normal);
    vec3 V = -ray.dir;
// Next event estimation
#ifdef USE_NEE
    {
      float lightPdf = sampleLightPdf(lightRect);
      vec3 pos = sampleLight(lightRect, getRandom());
      vec3 newPos = ray.origin + t * ray.dir;
      vec3 L = normalize(pos - newPos);
      float dist = distance(pos, newPos);
      dist = dist * dist;

      float tTmp;
      Material tMat;
      vec3 tNormal;
      bool hit =
          intersect(Ray(ray.origin + eps * L, L), true, tTmp, tMat, tNormal);
      if (hit && tMat.type == 2) {
        vec3 H = normalize(V + L);

        float d = D(normal, H, mat.roughness);
        float g = G(normal, V, L, mat.roughness);
        vec3 f = F(L, H, mat.color);
        vec3 brdf = d * g * f;
        brdf = brdf /
               (4 * clampDot(normal, H, false) * clampDot(normal, L, false));

#ifdef USE_VNDF
        float brdfPdf = sampleVndfPdf(V, H, d);
#else
        float brdfPdf = sampleNdfPdf(V, H, normal, d);
#endif

        float lightPdfW = pdfA2W(lightPdf, dist, dot(-L, tNormal));

        float misW = lightPdfW / (lightPdfW + brdfPdf);

        contrib += misW * tMat.color * tp * brdf / lightPdfW;
      }
    }
#endif
    {
      Material nextMat;
      vec3 nextNormal;

#ifdef USE_VNDF
      vec3 H = sampleVndf(V, mat.roughness, getRandom(), onb);
#else
      vec3 H = sampleNdf(mat.roughness, getRandom(), onb);
#endif

      vec3 L = normalize(reflect(ray.dir, H));

      // New ray
      ray = Ray(ray.origin + t * ray.dir + L * eps, L);

      // Didn't hit anything
      if (!intersect(ray, false, t, nextMat, nextNormal)) {
        break;
      }

      float d = D(normal, H, mat.roughness);
      float g = G(normal, V, L, mat.roughness);
      vec3 f = F(L, H, mat.color);
      vec3 brdf = d * g * f;
      brdf =
          brdf / (4 * clampDot(normal, H, false) * clampDot(normal, L, false));

#ifdef USE_VNDF
      float brdfPdf = sampleVndfPdf(V, H, d);
#else
      float brdfPdf = sampleNdfPdf(L, H, normal, d);
#endif

      if (nextMat.type == 2) {
#ifdef USE_NEE
        float lightPdf = sampleLightPdf(lightRect);
        sampleLight(lightRect, getRandom());
        float lightPdfW = pdfA2W(lightPdf, dot(L, L), dot(-L, nextNormal));
        float misW = brdfPdf / (lightPdfW + brdfPdf);
#else
        float misW = 1;
#endif

        contrib += misW * nextMat.color * tp * brdf / brdfPdf;
        break;
      }

      tp *= clampDot(normal, L, true) * brdf / brdfPdf;
      mat = nextMat;
      normal = nextNormal;
    }
  }

  return contrib;
}

void main() {
  seed = seedInit;
  flatIdx = int(dot(gl_FragCoord.xy, vec2(1, 4096)));

  // Generate sample
  Ray ray = genRay(getRandom());
  vec3 op = rayTrace(ray);

  // Get previous data
  vec3 prev = texture(prevFrame, gl_FragCoord.xy / resolution.xy).rgb;

  // Mix it to produce new frame
  fragColor = vec4(mix(op, prev, acc), 1.0);
}
