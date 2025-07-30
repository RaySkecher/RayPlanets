//testbench for raytracetemp
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

//#include "trace_path.h"

typedef struct {
    uint8_t r, g, b;
} Color;

extern Color trace_path(int x, int y);

/* Copy widths from the kernel */
#define WIDTH   256
#define HEIGHT  256
#define FRAC_BITS 12
#define ONE (1 << FRAC_BITS)
#define ANIMATION_TIME 0.0f

//function declare

extern Color trace_path(int x, int y);
extern void update_planet_positions(float t);
extern void setup_rings(void);
extern void update_ring_positions(void);


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
