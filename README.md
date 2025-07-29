# RayPlanets Raytracer

A fixed-point path tracing raytracer written in C that renders a scene with spheres using Monte Carlo sampling.

## Features

- **Fixed-point arithmetic**: Uses 4.12 fixed-point math for precision and performance
- **Path tracing**: Implements Monte Carlo path tracing for realistic lighting
- **Multiple sampling**: 16 samples per pixel for anti-aliasing and noise reduction
- **Direct lighting**: Calculates direct illumination from light sources
- **Indirect lighting**: Simulates light bounces for global illumination
- **Shadow calculation**: Implements shadow rays for realistic shadows

## Scene Description

The current scene contains:

- 3 spheres positioned in 3D space
- One light-emitting sphere (white) acting as the light source
- Two yellow/orange spheres as objects in the scene
- Camera positioned at (0, 0.8, 2) looking down the negative Z axis

## Quick Start

### One-command build and run (recommended)

```bash
make run
```

### Step-by-step commands

```bash
# Test compilation
make test

# Build the raytracer
make

# Render the image
make render

# View the rendered image
make view

# Convert to PNG (requires ImageMagick)
make png

# Clean build artifacts
make clean
```

### Manual compilation

```bash
# Compile
gcc -std=c99 -Wall -Wextra -O3 -march=native -ffast-math -o raytracer raytracer.c -lm

# Render
./raytracer > output.ppm

# View (macOS)
open output.ppm
```

## Output

The raytracer outputs a PPM (Portable Pixmap) image file. PPM is a simple image format that can be viewed by most image viewers and easily converted to other formats:

```bash
# Convert to PNG (requires ImageMagick)
convert output.ppm output.png

# Convert to JPEG (requires ImageMagick)
convert output.ppm output.jpg
```

## Configuration

Key parameters can be modified in `raytracer.c`:

- `HEIGHT`, `WIDTH`: Image resolution (currently 256x256)
- `NUM_SAMPLES`: Samples per pixel for anti-aliasing (currently 16)
- `MAX_BOUNCES`: Maximum light bounces (currently 3)
- `FOV`: Field of view in degrees (currently 60°)
- `LIGHT_INTENSITY`: Brightness of light sources

## Technical Details

### Fixed-Point Math

- Uses 16-bit fixed-point numbers with 4 integer and 12 fractional bits
- Provides good precision while maintaining performance
- Custom multiply, divide, and square root functions

### Rendering Algorithm

1. Cast rays from camera through each pixel
2. Find intersection with scene geometry (spheres)
3. Calculate direct lighting from light sources
4. Trace indirect lighting through random bounces
5. Accumulate color contributions
6. Average multiple samples per pixel

### Performance

- Optimized for compilation with `-O3 -march=native -ffast-math`
- Uses lookup tables for random unit vectors
- Rendering time: ~1-5 minutes depending on your system

## Viewing the Output

The rendered image will be saved as `output.ppm`. You can view it with:

- **macOS**: Preview (opens automatically with `open output.ppm`)
- **Windows**: Most image viewers support PPM
- **Linux**: GIMP, ImageMagick, or most image viewers
- **Web browsers**: Many browsers can display PPM files directly

## Dependencies

- GCC or compatible C compiler
- Standard C math library (`-lm`)
- No external graphics libraries required

## Future Goals

Originally planned for FPGA implementation with hopes to incorporate:

- Rotation matrices
- Textures
- Interactive camera
- Moons orbitting planets
- Stars/skybox

## Troubleshooting

**Long rendering time**: This is normal for path tracing. Reduce `NUM_SAMPLES` or image resolution for faster rendering.

**Black image**: Check that light sources are properly configured and visible from the camera.

**Compilation errors**: Ensure you have a C99-compatible compiler and the math library linked (`-lm`).

## License

Open source - feel free to modify and experiment with the code!
