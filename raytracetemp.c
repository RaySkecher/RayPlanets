#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Path tracing settings
#define MAX_BOUNCES 3
#define FOV 70.0
#define LIGHT_INTENSITY 2.0
#define BRIGHTNESS_SHIFT 4
#define NUM_SAMPLES 16

// Animation settings
#define ANIMATION_TIME 0.0f  // Current time (can be modified for different frames)
#define PLANET1_ORBIT_RADIUS 2.5f
#define PLANET1_ORBIT_SPEED 1.0f
#define PLANET2_ORBIT_RADIUS 1.6f  
#define PLANET2_ORBIT_SPEED 0.6f

// Scene dimensions
#define NUM_SPHERES 3
#define NUM_RINGS 5

#define HEIGHT 256
#define WIDTH 256

// Fixed-point math settings — 16-bit total (4 integer + 12 fractional)
typedef int16_t fp_t;
#define FRAC_BITS 12
#define ONE (1 << FRAC_BITS)
// Convert a floating-point value to 4.12 fixed-point (rounded)
#define F(x) ((fp_t)((x) * ONE + ((x) >= 0 ? 0.5f : -0.5f)))
#define I(x) ((x) >> FRAC_BITS)

// Handy constants
#define FP_EPS ((fp_t)1)             // ≈ 0.00024 in real units
#define FP_INF 0x7FFFFFFF            // Large “infinite” distance sentinel

//Textures
#define COLOR_LAYERS 3


// Simple vector struct
typedef struct {
    fp_t x, y, z;
} Vec3;

typedef struct {
    uint8_t r, g, b;
} Color;

// Ray
typedef struct {
    Vec3 orig, dir;
} Ray;

// Material
typedef struct {
    Vec3 color;
    int is_light;
} Material;

// Sphere
typedef struct {
    Vec3 center;
    fp_t radius;
    Material material;
} Sphere;

// Structure to hold intersection results
typedef struct {
    int32_t t;            // keep 32-bit for extra head-room during comparisons
    int hit_index;
    int hit_type;        // 0 = sphere, 1 = ring
    int hit;      // 1 if an object was hit, 0 otherwise
} Intersection;

//Ring
typedef struct {
    Vec3 center;
    Vec3 normal;
    fp_t inner_radius;
    fp_t outer_radius;
    //fp_t ellipse_ratio;
    Material material;
} Ring;

typedef struct {
    float x, y, z;
} Vec3f;

#define UNIT_VECTOR_LUT_SIZE 128
// 64 pre-computed unit vectors stored at 8.8 precision; they will be left-shifted
// at runtime to 4.12 to keep the table small and readable.
static const Vec3 g_unit_vector_lut[UNIT_VECTOR_LUT_SIZE] = {
    {.x = 148,  .y = 148,  .z = 148 }, {.x = -148, .y = -148, .z = -148},
    {.x = 210,  .y = 0,    .z = 148 }, {.x = -210, .y = 0,    .z = -148},
    {.x = 0,    .y = 210,  .z = 148 }, {.x = 0,    .y = -210, .z = -148},
    {.x = 148,  .y = 210,  .z = 0   }, {.x = -148, .y = -210, .z = 0   },
    {.x = 256,  .y = 0,    .z = 0   }, {.x = -256, .y = 0,    .z = 0   },
    {.x = 0,    .y = 256,  .z = 0   }, {.x = 0,    .y = -256, .z = 0   },
    {.x = 0,    .y = 0,    .z = 256 }, {.x = 0,    .y = 0,    .z = -256},
    {.x = 182,  .y = 182,  .z = 0   }, {.x = -182, .y = -182, .z = 0   },
    {.x = 182,  .y = 0,    .z = 182 }, {.x = -182, .y = 0,    .z = -182},
    {.x = 0,    .y = 182,  .z = 182 }, {.x = 0,    .y = -182, .z = -182},
    {.x = 102,  .y = 102,  .z = 226 }, {.x = -102, .y = -102, .z = -226},
    {.x = 102,  .y = 226,  .z = 102 }, {.x = -102, .y = -226, .z = -102},
    {.x = 226,  .y = 102,  .z = 102 }, {.x = -226, .y = -102, .z = -102},
    {.x = 128,  .y = 128,  .z = 182 }, {.x = -128, .y = -128, .z = -182},
    {.x = 128,  .y = 182,  .z = 128 }, {.x = -128, .y = -182, .z = -128},
    {.x = 182,  .y = 128,  .z = 128 }, {.x = -182, .y = -128, .z = -128},
    {.x =  48,  .y = 252,  .z =  39 }, {.x =  -48, .y = -252, .z =  -39},
    {.x = 252,  .y =  39,  .z =  48 }, {.x = -252, .y =  -39, .z =  -48},
    {.x =  39,  .y =  48,  .z = 252 }, {.x =  -39, .y =  -48, .z = -252},
    {.x = 156,  .y = 195,  .z =   5 }, {.x = -156, .y = -195, .z =   -5},
    {.x = 195,  .y =   5,  .z = 156 }, {.x = -195, .y =   -5, .z = -156},
    {.x =   5,  .y = 156,  .z = 195 }, {.x =   -5, .y = -156, .z = -195},
    {.x =  70,  .y =  70,  .z = 238 }, {.x =  -70, .y =  -70, .z = -238},
    {.x =  70,  .y = 238,  .z =  70 }, {.x =  -70, .y = -238, .z =  -70},
    {.x = 238,  .y =  70,  .z =  70 }, {.x = -238, .y =  -70, .z =  -70},
    {.x = 217,  .y =  86,  .z =  98 }, {.x = -217, .y =  -86, .z =  -98},
    {.x =  86,  .y =  98,  .z = 217 }, {.x =  -86, .y =  -98, .z = -217},
    {.x =  98,  .y = 217,  .z =  86 }, {.x =  -98, .y = -217, .z =  -86},
    {.x = 239,  .y =  48,  .z =  78 }, {.x = -239, .y =  -48, .z =  -78},
    {.x =  48,  .y =  78,  .z = 239 }, {.x =  -48, .y =  -78, .z = -239},
    {.x =  78,  .y = 239,  .z =  48 }, {.x =  -78, .y = -239, .z =  -48},
    {.x = 180, .y =  34, .z = 193}, {.x = -180, .y =  34, .z = -193},
    {.x =  34, .y = 193, .z = 180}, {.x =  -34, .y = -193, .z = -180},
    {.x = 193, .y = 180, .z =  34}, {.x = -193, .y = -180, .z =  -34},
    {.x =  45, .y = 210, .z = 160}, {.x =  -45, .y = -210, .z = -160},
    {.x = 210, .y = 160, .z =  45}, {.x = -210, .y = -160, .z =  -45},
    {.x = 160, .y =  45, .z = 210}, {.x = -160, .y =  -45, .z = -210},
    {.x = 125, .y = 230, .z =  70}, {.x = -125, .y = -230, .z =  -70},
    {.x = 230, .y =  70, .z = 125}, {.x = -230, .y =  -70, .z = -125},
    {.x =  70, .y = 125, .z = 230}, {.x =  -70, .y = -125, .z = -230},
    {.x = 190, .y = 200, .z =  20}, {.x = -190, .y = -200, .z =  -20},
    {.x = 200, .y =  20, .z = 190}, {.x = -200, .y =  -20, .z = -190},
    {.x =  20, .y = 190, .z = 200}, {.x =  -20, .y = -190, .z = -200},
    {.x =  95, .y = 215, .z = 140}, {.x =  -95, .y = -215, .z = -140},
    {.x = 215, .y = 140, .z =  95}, {.x = -215, .y = -140, .z =  -95},
    {.x = 140, .y =  95, .z = 215}, {.x = -140, .y =  -95, .z = -215},
    {.x =  60, .y = 180, .z = 230}, {.x =  -60, .y = -180, .z = -230},
    {.x = 180, .y = 230, .z =  60}, {.x = -180, .y = -230, .z =  -60},
    {.x = 230, .y =  60, .z = 180}, {.x = -230, .y =  -60, .z = -180},
    {.x = 130, .y = 200, .z = 150}, {.x = -130, .y = -200, .z = -150},
    {.x = 200, .y = 150, .z = 130}, {.x = -200, .y = -150, .z = -130},
    {.x = 150, .y = 130, .z = 200}, {.x = -150, .y = -130, .z = -200},
    {.x =  85, .y = 145, .z = 240}, {.x =  -85, .y = -145, .z = -240},
    {.x = 145, .y = 240, .z =  85}, {.x = -145, .y = -240, .z =  -85},
    {.x = 240, .y =  85, .z = 145}, {.x = -240, .y =  -85, .z = -145},
    {.x = 170, .y = 100, .z = 220}, {.x = -170, .y = -100, .z = -220},
    {.x = 100, .y = 220, .z = 170}, {.x = -100, .y = -220, .z = -170},
    {.x = 220, .y = 170, .z = 100}, {.x = -220, .y = -170, .z = -100},
    {.x =  66, .y = 210, .z = 195}, {.x =  -66, .y = -210, .z = -195},
    {.x = 210, .y = 195, .z =  66}, {.x = -210, .y = -195, .z =  -66},
    {.x = 195, .y =  66, .z = 210}, {.x = -195, .y =  -66, .z = -210},
    {.x = 110, .y = 180, .z = 210}, {.x = -110, .y = -180, .z = -210}
};

static Sphere g_spheres[NUM_SPHERES] = {
    {.center = {.x = F(0.7), .y = F(0.75), .z = F(0.1)}, .radius = F(0.2), .material = {.color = {.x = F(0.8), .y = F(0.6), .z = F(0.3)}, .is_light = 0}}, // Planet 1 - orange/brown
    {.center = {.x = F(0.0), .y = F(0.75), .z = F(-0.4)}, .radius = F(0.4), .material = {.color = {.x = F(LIGHT_INTENSITY), .y = F(LIGHT_INTENSITY), .z = F(LIGHT_INTENSITY)}, .is_light = 1}},  // Sun - white light source
    {.center = {.x = F(-0.8), .y = F(0.75), .z = F(-0.7)}, .radius = F(0.15), .material = {.color = {.x = F(0.4), .y = F(0.6), .z = F(0.9)}, .is_light = 0}}, // Planet 2 - blue (smaller)
   
};

static const Vec3 gradient_ray_perlin[12] = {
    {1, 1, 0},{-1, 1, 0},{1 ,-1, 0},{-1, -1, 0},
    {1, 0, 1},{-1, 0, 1},{1 ,0, -1},{-1, 0, -1},
    {0, 1, 1},{0, -1, 1},{0, 1, -1},{0, -1, -1},
    {1, 1, 0},{-1, 1, 0},{1 ,-1, 0},{-1, -1, 0},
};



// Fixed-point multiplication
// 8.8 fixed-point multiply (fp_t × fp_t → 32-bit)

static int32_t mul(int32_t a, int32_t b) {
    return (int32_t)(((int64_t)a * (int64_t)b) >> FRAC_BITS);
}

// Fixed-point division
int32_t div_fp(int32_t a, int32_t b) {
    if (b == 0) return 0;
    //int32_t X0 = 
    return (int32_t)(((int64_t)a << FRAC_BITS) / b);
}


// Fast inverse square root for fixed-point numbers
static int32_t inv_sqrt_fp(int32_t x) {
    if (x <= 0) return 0;
    float x_f = (float)x / (float)ONE;
    union { float f; uint32_t i; } u;
    u.f = x_f;
    u.i = 0x5f3759df - (u.i >> 1);
    u.f = u.f * (1.5f - 0.5f * x_f * u.f * u.f);
    return (int32_t)(u.f * (float)ONE);
}

// Fixed-point square root
int32_t sqrt_fp(int32_t n) {
    if (n <= 0) return 0;
    return mul(n, inv_sqrt_fp(n));
}

/* Fixed‑point → byte, with brightness boost */
static inline int fp_to_u8(fp_t v)
{
    /* scale by 2^BRIGHTNESS_SHIFT, then normalise to 0…255 */
    int32_t disp = v << BRIGHTNESS_SHIFT;
    int val = ((int64_t)disp * 255) >> FRAC_BITS;
    return val > 255 ? 255 : (val < 0 ? 0 : val);
}

// Random number generation
static uint32_t rand_state = 12345;
static uint32_t rand_u32() { // maybe make it generate 3 randon numbers each clock cycle?
    // xorshift
    rand_state ^= rand_state << 13;
    rand_state ^= rand_state >> 17;
    rand_state ^= rand_state << 5;
    return rand_state;
}

int32_t rand_fp() {
    return (int32_t)((uint64_t)rand_u32() * ONE >> 32);
}


Vec3 random_unit_vector() {
    uint32_t r_val = rand_u32();
    int lut_idx = r_val & 0x7F;
    Vec3 base = g_unit_vector_lut[lut_idx];
    // Convert from 8.8 → 4.12 by left-shifting 4 bits.
    return (Vec3){ (fp_t)(base.x << (FRAC_BITS - 8)),
                   (fp_t)(base.y << (FRAC_BITS - 8)),
                   (fp_t)(base.z << (FRAC_BITS - 8)) };
}

// Simple hash function for Vec3 to create deterministic stars
static uint32_t hash_vec3(Vec3 v) {
    uint32_t h = (uint32_t)v.x;
    h = (h * 73856093) ^ ((uint32_t)v.y * 19349663);
    h = (h * 83492791) ^ ((uint32_t)v.z * 101390427);
    return h;
}

// Function to get a deterministic pseudo-random fixed-point number from a seed [0, ONE)
static fp_t get_deterministic_rand_fp(uint32_t seed) {
    uint32_t s = seed;
    s ^= s << 13;
    s ^= s >> 17;
    s ^= s << 5;
    // Take the most significant FRAC_BITS for the fractional part.
    // This effectively scales the 32-bit random number to the [0, ONE) fixed-point range.
    return (fp_t)(s >> (32 - FRAC_BITS));
}

// Forward declarations
int is_on_light(Vec3 p);
int32_t intersect_sphere(Ray r, Sphere s);
Intersection intersect_scene(Ray r);
int32_t intersect_ring(Ray r, Ring ring);


// Vector operations
Vec3 vec_add(Vec3 a, Vec3 b) { return (Vec3){a.x + b.x, a.y + b.y, a.z + b.z}; }
Vec3 vec_sub(Vec3 a, Vec3 b) { return (Vec3){a.x - b.x, a.y - b.y, a.z - b.z}; }
int32_t vec_dot(Vec3 a, Vec3 b) { return mul(a.x, b.x) + mul(a.y, b.y) + mul(a.z, b.z); }
Vec3 vec_mul(Vec3 a, Vec3 b) { return (Vec3){ (fp_t)mul(a.x, b.x), (fp_t)mul(a.y, b.y), (fp_t)mul(a.z, b.z)}; }
Vec3 vec_scale(Vec3 v, int32_t s) {
    fp_t s_fp = (fp_t)s; // clamp / cast the scalar into fp_t range
    return (Vec3){ (fp_t)mul(v.x, s_fp), (fp_t)mul(v.y, s_fp), (fp_t)mul(v.z, s_fp)};
}
int32_t vec_len_sq(Vec3 v) { return vec_dot(v, v); }
Vec3 vec_norm(Vec3 v) {
    int32_t len_sq = vec_len_sq(v);
    int32_t inv_len = inv_sqrt_fp(len_sq);
    return vec_scale(v, inv_len);
}

Vec3f vec3_to_float(Vec3 v) {
    return (Vec3f){
        (float)v.x / ONE,
        (float)v.y / ONE,
        (float)v.z / ONE
    };
}

// Rotation matrix operations for orbital motion
Vec3 rotate_y(Vec3 v, float angle) {
    float cos_a = cosf(angle);
    float sin_a = sinf(angle);
    
    return (Vec3){
        F(cos_a * I(v.x) + sin_a * I(v.z)),
        v.y,
        F(-sin_a * I(v.x) + cos_a * I(v.z))
    };
}

// Calculate orbital position for a planet
Vec3 get_orbital_position(Vec3 sun_center, float orbit_radius, float orbit_speed, float time) {
    float angle = orbit_speed * time;
    Vec3 offset = {F(orbit_radius), F(0), F(0)};
    Vec3 rotated_offset = rotate_y(offset, angle);
    return vec_add(sun_center, rotated_offset);
}

// Update planet positions based on current time
void update_planet_positions(float time) {
    Vec3 sun_pos = g_spheres[1].center; // Sun is at index 1
    
    // Update Planet 1 (index 0)
    g_spheres[0].center = get_orbital_position(sun_pos, PLANET1_ORBIT_RADIUS, PLANET1_ORBIT_SPEED, time);
    
    // Update Planet 2 (index 2) 
    g_spheres[2].center = get_orbital_position(sun_pos, PLANET2_ORBIT_RADIUS, PLANET2_ORBIT_SPEED, time);
}

Ring g_rings[NUM_RINGS];

void update_ring_positions() {
    for (int i = 0; i < 2; ++i) g_rings[i].center = g_spheres[2].center;
    for (int i = 2 ; i < NUM_RINGS ; ++i ) g_rings[i].center = g_spheres[0].center;
}

void setup_rings() {
    Vec3 jupiter_center = g_spheres[2].center;  // Position of Jupiter
    Vec3 orange_center = g_spheres[0].center;
    // Manually customized ring sizes and colors
    fp_t inner_radii[NUM_RINGS] = {F(0.25), F(0.35), F(0.35)};  // > planet radius (0.15)
    fp_t outer_radii[NUM_RINGS] = {F(0.3), F(0.4), F(0.5)};
    //fp_t ellipse_ratios[NUM_RINGS] = {F(1.0), F(0.95)};
    Vec3 colors[NUM_RINGS] = {
        {F(0.4), F(0.4), F(0.5)},
        {F(0.2), F(0.2), F(0.3)},
        {F(0.5), F(0.3), F(0.1)},
    };

    //blue planet
    for (int i = 0; i < 2; ++i) {
         g_rings[i].center = jupiter_center;

        // All rings lie in roughly the same flat orbital plane
        g_rings[i].normal = vec_norm((Vec3){F(0.1), F(1), F(0.1)});

        g_rings[i].inner_radius = inner_radii[i];
        g_rings[i].outer_radius = outer_radii[i];
        //g_rings[i].ellipse_ratio = ellipse_ratios[i];

        g_rings[i].material = (Material){
            .color = colors[i],
            .is_light = 0
        };
    }

    //orange planet
    for (int i = 2; i < 4; ++i) {
         g_rings[i].center = orange_center;

        // All rings lie in roughly the same flat orbital plane
        g_rings[i].normal = vec_norm((Vec3){F(-0.5), F(1), F(0)});

        g_rings[i].inner_radius = inner_radii[i];
        g_rings[i].outer_radius = outer_radii[i];

        g_rings[i].material = (Material){
            .color = colors[i],
            .is_light = 0
        };
    }
}

// Returns an Intersection result.
Intersection intersect_scene(Ray r) {
    Intersection result = {.t = FP_INF, .hit = 0, .hit_index = -1, .hit_type = -1};
    int32_t closest_t = FP_INF;

    // Check spheres
    for (int i = 0; i < NUM_SPHERES; ++i) {
        int32_t t = intersect_sphere(r, g_spheres[i]);
        if (t < closest_t) {
            closest_t = t;
            result.t = t;
            result.hit_index = i;
            result.hit_type = 0;  // sphere
            result.hit = 1;
        }
    }

    // Check rings
    for (int i = 0; i < NUM_RINGS; ++i) {
        int32_t t = intersect_ring(r, g_rings[i]);
        if (t < closest_t) {
            closest_t = t;
            result.t = t;
            result.hit_index = i;
            result.hit_type = 1;  // ring
            result.hit = 1;
        }
    }

    return result;
}
// Ray-sphere intersection
int32_t intersect_sphere(Ray r, Sphere s) {
    //#pragma HLS ALLOCATION function instances=mul limit=2
    Vec3 oc = vec_sub(r.orig, s.center);
    int32_t a = vec_dot(r.dir, r.dir);
    int32_t b = 2 * vec_dot(oc, r.dir);
    int32_t c = vec_dot(oc, oc) - mul(s.radius, s.radius);
    int32_t discriminant = mul(b, b) - 4 * mul(a, c);
    if (discriminant < 0) return FP_INF;
    
    int32_t sqrt_d = sqrt_fp(discriminant);
    // Pre-compute 1/(2a) once.
    int32_t inv_2a = div_fp(ONE, 2 * a);

    int32_t t  = mul(-b - sqrt_d, inv_2a);
    if (t > FP_EPS) return t;

    int32_t t2 = mul(-b + sqrt_d, inv_2a);
    if (t2 > FP_EPS) return t2;
    
    return FP_INF;
}

int32_t intersect_ring(Ray r, Ring ring) {
    // Ring plane intersection (dot(n, P - C) = 0)
    int32_t denom = vec_dot(ring.normal, r.dir);
    if (abs(denom) < F(0.0001)) return FP_INF;  // Parallel

    int32_t t = div_fp(vec_dot(ring.normal, vec_sub(ring.center, r.orig)), denom);
    if (t < 0) return FP_INF;

    Vec3 p = vec_add(r.orig, vec_scale(r.dir, t));
    Vec3 rel = vec_sub(p, ring.center);
    int32_t dist2 = vec_dot(rel, rel);  // use circular distance

    if (dist2 < mul(ring.inner_radius, ring.inner_radius) ||
        dist2 > mul(ring.outer_radius, ring.outer_radius))
        return FP_INF;

    return t;
}

Vec3 randomGradient(int32_t x, int32_t y , int32_t z){
    // Simple hash function that produces values 0-255
    int32_t hash_value = (int32_t)((x * 16197654 + y * 31337764 + z * 6971234) & 0xFF
    );
    
    //Index to 16;
    int idx = hash_value & 15;
    
    // Return the corresponding gradient vector
    return gradient_ray_perlin[idx];
}

float Fade(float t) {
	return ((6*t - 15)*t + 10)*t*t*t; //smoothing function
}

float Lerp(float t, float a1, float a2) {
	return a1 + t*(a2-a1);
}

float noise(Vec3 hit_point){
    //Point 0,0,0

    Vec3f hit_pointf = vec3_to_float(hit_point);

    int32_t x0 = (int32_t)(floor(hit_pointf.x));
    int32_t y0 = (int32_t)floor(hit_pointf.y);
    int32_t z0 = (int32_t)floor(hit_pointf.z);

    //Point 1,1,0 
    int32_t x1 = (x0 + 1);
    int32_t y1 = (y0 + 1);
    int32_t z1 = (z0 + 1);

    //Interpolation weights
    float xf = hit_pointf.x - x0; 
    float yf = hit_pointf.y - y0; 
    float zf = hit_pointf.z - z0; 

    //Smoothing function
    float u = Fade(xf);
    float v = Fade(yf);
    float w = Fade(zf);

    //Get vector for each of the eight other corners
    Vec3f c000 = vec3_to_float(randomGradient(x0, y0, z0));
    Vec3f c100 = vec3_to_float(randomGradient(x1, y0, z0));
    Vec3f c010 = vec3_to_float(randomGradient(x0, y1, z0));
    Vec3f c110 = vec3_to_float(randomGradient(x1, y1, z0));
    Vec3f c001 = vec3_to_float(randomGradient(x0, y0, z1));
    Vec3f c101 = vec3_to_float(randomGradient(x1, y0, z1));
    Vec3f c011 = vec3_to_float(randomGradient(x0, y1, z1));
    Vec3f c111 = vec3_to_float(randomGradient(x1, y1, z1));

    // Compute dot products 
    float d000 = c000.x * xf + c000.y * yf + c000.z * zf;
    float d100 = c100.x * (xf-1.0f) + c100.y * yf + c100.z * zf;
    float d010 = c010.x * xf + c010.y * (yf-1.0f) + c010.z * zf;
    float d110 = c110.x * (xf-1.0f) + c110.y * (yf-1.0f) + c110.z * zf;
    float d001 = c001.x * xf + c001.y * yf + c001.z * (zf-1.0f);
    float d101 = c101.x * (xf-1.0f) + c101.y * yf + c101.z * (zf-1.0f);
    float d011 = c011.x * xf + c011.y * (yf-1.0f) + c011.z * (zf-1.0f);
    float d111 = c111.x * (xf-1.0f) + c111.y * (yf-1.0f) + c111.z * (zf-1.0f);

    //Interpolate along x axis 
    float x00 = Lerp(u, d000, d100);
    float x10 = Lerp(u, d010, d110);
    float x01 = Lerp(u, d001, d101);
    float x11 = Lerp(u, d011, d111);

    // Interpolate along y axis (2 lerps)
    float yt0 = Lerp(v, x00, x10);
    float yt1 = Lerp(v, x01, x11);

    // Final interpolation along z axis
    return Lerp(w, yt0, yt1); //Return float [-1,1]
}


float layerdPerlin_noise(Vec3 positiion, int layers){
    float frequency = 1.0f; //How much detail
    float amplitude = 1.0f; //How sharp are the details
    float total = 0.0f;
    float max_val = 0.0f;
    if (layers != 0){
        for (int i = 0; i < layers; i++){
            Vec3 scaled_position = {positiion.x * frequency, positiion.y * frequency, positiion.z * frequency};
            float noise_level = noise(scaled_position) * amplitude; //convets to float
            //Adding the values 
            total = total + noise_level;
            max_val = max_val + amplitude;
            //Scale down the sharpness of each layer
            amplitude = amplitude * 0.5f;
            frequency = frequency * 2.0f;
        }
    }
    else{
        return (noise(positiion)+1.0f)*0.5f;
    }
    return ((total/max_val)+1.0f)*0.5f;
}



Vec3 convert_color(Vec3 colorRGB){
    return (Vec3){F(colorRGB.x/255.0),F(colorRGB.y/255.0),F(colorRGB.z/255.0)};
}

Vec3 getPlanetColor(Vec3 hit_point, Vec3 hit_normal, int hit_object_index){
    int32_t lattitude = abs(hit_normal.y);
    int32_t terrain_noise_fp = F(layerdPerlin_noise(hit_point, COLOR_LAYERS));
    Vec3 color[COLOR_LAYERS];

    for (int i = 0; i < COLOR_LAYERS; i++){
        color[i] = g_spheres[hit_object_index].material.color;
    }

    //Terrian_boundries 
    int32_t ice_end = F(0.8) + mul(terrain_noise_fp, F(0.1));
    int32_t forest_start = F(0.8) + mul(terrain_noise_fp, F(0.1));
    int32_t forest_end = F(0.4) + mul(terrain_noise_fp, F(0.1));
    int32_t desert_start = F(0.4) + mul(terrain_noise_fp, F(0.1));

    
    //Colors
    switch (hit_object_index)
    {
    case 0: //outer planet
        color[0] = convert_color((Vec3){153, 41, 234}); //rgb(153,41,234);
        color[1] = convert_color((Vec3){204, 102, 218}); //rgb(204, 102, 218);
        color[2] = convert_color((Vec3){250, 235, 146}); //rgb(250, 235, 146);
        break;
    case 2: //inner planet
        color[0]= convert_color((Vec3){216,174,109}); //rgb(216,174,109);
        color[1] = convert_color((Vec3){255,225,171}); //rgb(255,225,171);
        color[2] = convert_color((Vec3){219,181,124}); //rgb(219,181,124);
    default:
        break;
    }

    //Caculate colors
    //F(value = percent of surface)
    if (lattitude >= ice_end){
        return color[0];
    }
    else if (lattitude >= forest_end && lattitude < forest_start){
        return color[1];
    }
    else if (lattitude < desert_start){
        return color[2];
    }
    else return (Vec3){255,255,255};
}


Color trace_path(int16_t x, int16_t y) {
    //#pragma HLS ALLOCATION function instances=mul limit=32
    //#pragma HLS ALLOCATION function instances=div_fp limit=8
    //#pragma HLS ALLOCATION function instances=intersect_plane limit=1
    //#pragma HLS ALLOCATION function instances=intersect_sphere limit=1
    #pragma HLS bind_storage variable=g_unit_vector_lut type=rom_1p

    Ray cam = {{F(0), F(1.2), F(3)}, {F(0), F(0), F(-1)}};
    float fov_rad = FOV * M_PI / 180.0;
    float fov_scale = tan(fov_rad / 2.0);
    //float aspect_ratio = (float)WIDTH / HEIGHT;

    float sx_ndc = (2.0f * (x) / WIDTH) - 1.0f;
    float sy_ndc = 1.0f - (2.0f * (y) / HEIGHT);

    Ray r;

    int32_t acc_r = 0, acc_g = 0, acc_b = 0;

    for (int sample = 0; sample < NUM_SAMPLES; sample++) {
        r = cam;
        r.dir.x = F(sx_ndc * fov_scale);
        r.dir.y = F(sy_ndc * fov_scale);
        r.dir = vec_norm(r.dir);
        Vec3 path_color = {F(0), F(0), F(0)};
        Vec3 path_attenuation = {ONE, ONE, ONE};

        for (int b = 0; b < MAX_BOUNCES; ++b) {
            #pragma HLS pipeline off
            Intersection inter = intersect_scene(r);

            if (!inter.hit) {
                // deep_purple_blue: now approximates ~RGB(40, 20, 100)
                // light_purple_blue: now approximates ~RGB(100, 80, 200)
                //Both values adjusted for BRIGHTNESS_SHIFT
                Vec3 deep_purple_blue = {F(0.0098), F(0.0049), F(0.0245)}; 
                Vec3 light_purple_blue = {F(0.0245), F(0.0196), F(0.049)};

                // Get interpolation factor based on ray's Y-direction (vertical position)
                // Map r.dir.y from [-ONE, ONE] to [0, ONE]
                // r.dir.y is in fixed-point, range approx [-1.0, 1.0]
                // To map to [0.0, 1.0], we do (r.dir.y + 1.0) / 2.0
                fp_t interp_factor_y = div_fp(r.dir.y + ONE, F(2.0));

                // Clamp interpolation factor to ensure it's within [0, ONE]
                if (interp_factor_y < 0) interp_factor_y = 0;
                if (interp_factor_y > ONE) interp_factor_y = ONE;

                Vec3 base_sky_color;
                base_sky_color.x = deep_purple_blue.x + mul(light_purple_blue.x - deep_purple_blue.x, interp_factor_y);
                base_sky_color.y = deep_purple_blue.y + mul(light_purple_blue.y - deep_purple_blue.y, interp_factor_y);
                base_sky_color.z = deep_purple_blue.z + mul(light_purple_blue.z - deep_purple_blue.z, interp_factor_y);

                Vec3 sky_color = base_sky_color; // Start with the interpolated color

                // Add stars based on ray direction (deterministic)
                uint32_t star_seed = hash_vec3(r.dir);
                fp_t star_prob = get_deterministic_rand_fp(star_seed);

                // Adjust this threshold for star density (smaller value = fewer stars)
                fp_t star_threshold = F(0.0003f); // Fixed-point threshold

                if (star_prob < star_threshold) {
                    fp_t intensity_factor = div_fp(star_prob, star_threshold);
                    fp_t one_minus_factor = ONE - intensity_factor;
                    fp_t star_intensity_fp = F(0.5) + mul(one_minus_factor, F(0.5));

                    // Add star intensity to the interpolated sky color
                    sky_color.x = sky_color.x + star_intensity_fp;
                    sky_color.y = sky_color.y + star_intensity_fp;
                    sky_color.z = sky_color.z + star_intensity_fp;
                }

                path_color = vec_add(path_color, vec_mul(path_attenuation, sky_color));
                path_attenuation = (Vec3){F(0), F(0), F(0)}; // Stop path tracing
                break; // Path ends in sky
            }

            int32_t t = inter.t;
            int hit_object_index = inter.hit_index;

            Vec3 hit_point = vec_add(r.orig, vec_scale(r.dir, t));
            Vec3 hit_normal;
            Material mat;

           if (inter.hit_type == 0) { // Sphere
                mat = g_spheres[hit_object_index].material;
                hit_normal = vec_norm(vec_sub(hit_point, g_spheres[hit_object_index].center));
                if (!mat.is_light) {
                    mat.color = getPlanetColor(hit_point,hit_normal, hit_object_index);
               } 
            } else if (inter.hit_type == 1) { // Ring
                mat = g_rings[hit_object_index].material;
                hit_normal = g_rings[hit_object_index].normal;
                if (vec_dot(hit_normal, r.dir) > 0) {
                    hit_normal.x = -hit_normal.x;
                    hit_normal.y = -hit_normal.y;
                    hit_normal.z = -hit_normal.z;
                }
            }

            Material surface_mat = mat;
            if (mat.is_light) { // If we hit the ceiling plane
                path_color = vec_add(path_color, vec_mul(path_attenuation, mat.color));
                path_attenuation = (Vec3){F(0), F(0), F(0)};
                break;
            }

            // Check if it is in a shadow
            // Shadow ray to the light source (sphere 2)
            Vec3 light_center = g_spheres[1].center; // Assuming sphere 2 is the light
            Vec3 light_point = vec_add(light_center, vec_scale(random_unit_vector(), g_spheres[1].radius));
            Vec3 light_vec = vec_sub(light_point, hit_point);
            int32_t dist_sq = vec_len_sq(light_vec);
            Vec3 light_dir = vec_norm(light_vec);

            
            Ray shadow_ray = {vec_add(hit_point, vec_scale(hit_normal, F(0.01))), light_dir};
            int occluded = 0;
            for (size_t i = 0; i < NUM_SPHERES; ++i) {
                if (i == hit_object_index || i == 1) continue; // Skip the hit object and light source
                int32_t shadow_t = intersect_sphere(shadow_ray, g_spheres[i]);
                if (shadow_t < FP_INF && mul(shadow_t, shadow_t) < dist_sq) {
                    occluded = 1;
                    break;
                }
            }

            for (size_t i = 0; i < NUM_RINGS; ++i) {
                if (i == hit_object_index) continue; // Skip the hit ring
                int32_t shadow_t = intersect_ring(shadow_ray, g_rings[i]);
                if (shadow_t < FP_INF && mul(shadow_t, shadow_t) < dist_sq) {
                    occluded = 1;
                    break;
                }
            }

            if (!occluded) {
                // if it is NOT in a shadow, calculate the direct light contribution
                int32_t cos_theta = vec_dot(hit_normal, light_dir);
                Vec3 light_normal = vec_norm(vec_sub(light_point, light_center));;
                int32_t cos_alpha = vec_dot(light_normal, light_dir);

                if (cos_theta > 0 && cos_alpha > 0) { // if the light and surface are facing each other
                    Material light_mat = g_spheres[1].material;
                    int32_t light_area = mul(g_spheres[1].radius, g_spheres[1].radius); // Approximate light area
                    int32_t geom_term_num = mul(cos_theta, cos_alpha);
                    int32_t geom_term = div_fp(geom_term_num, dist_sq);

                    Vec3 direct_light = vec_mul(path_attenuation, mat.color);
                    direct_light = vec_mul(direct_light, light_mat.color);
                    direct_light = vec_scale(direct_light, geom_term);
                    direct_light = vec_scale(direct_light, light_area);
                    // Divide by PI for diffuse BRDF
                    direct_light = vec_scale(direct_light, F(0.8)); 
                    path_color = vec_add(path_color, direct_light);
                }
            }

            // Attenuate path for next bounce (indirect light)
            path_attenuation = vec_mul(path_attenuation, surface_mat.color);
            
            // New random direction for bounced ray
            Vec3 random_dir = random_unit_vector();
            Vec3 bounce_dir = vec_add(hit_normal, random_dir);

            r.orig = vec_add(hit_point, vec_scale(hit_normal, F(0.01)));
            r.dir = vec_norm(bounce_dir);
        }
        acc_r += path_color.x;
        acc_g += path_color.y;
        acc_b += path_color.z;
    }

    return (Color){fp_to_u8((fp_t)(acc_r / NUM_SAMPLES)), fp_to_u8((fp_t)(acc_g / NUM_SAMPLES)), fp_to_u8((fp_t)(acc_b / NUM_SAMPLES))};
}

// New top-level function for Vitis HLS
void render_frame(float animation_time, Color output_buffer[HEIGHT][WIDTH]) {
    // HLS INTERFACE pragmas define how arguments map to hardware ports
    #pragma HLS INTERFACE s_axilite port=animation_time bundle=control
    #pragma HLS INTERFACE m_axi port=output_buffer offset=slave bundle=gmem
    #pragma HLS INTERFACE s_axilite port=return bundle=control

    // Initialize ring properties once per frame (or once if animation_time is fixed)
    setup_rings(); 

    // Update planet and ring positions for current animation time
    update_planet_positions(animation_time);
    update_ring_positions(); 
    
    // Render each pixel
    #pragma HLS DATAFLOW
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            #pragma HLS pipeline II=1
            output_buffer[y][x] = trace_path(x, y);
        }
    }
}

//COMMENT OUT main IF YOU'RE RUNNING THIS ON VITIS!

int main(int argc, char *argv[]) {
    // Parse animation time from command line (default to 0.0)
    float animation_time = ANIMATION_TIME;
    if (argc > 1) {
        animation_time = atof(argv[1]);
    }

    // Update planet positions for current animation time
    update_planet_positions(animation_time);
    setup_rings();
    update_ring_positions();

    // Create PPM image header
    printf("P3\n");
    printf("%d %d\n", WIDTH, HEIGHT);
    printf("255\n"); 
    
    // Render each pixel
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            Color pixel = trace_path(x, y);
            printf("%d %d %d ", pixel.r, pixel.g, pixel.b);
        }
        printf("\n");
        
        // Progress indicator to stderr
        if (y % 32 == 0) {
            fprintf(stderr, "Progress: %.1f%%\n", (float)y / HEIGHT * 100.0);
        }
    }
    
    fprintf(stderr, "Rendering complete!\n");
    return 0;
} 
