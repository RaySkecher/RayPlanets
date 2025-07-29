# RayPlanets Raytracer

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



**Compilation errors**: Ensure you have a C99-compatible compiler and the math library linked (`-lm`).

## License

Open source - feel free to modify and experiment with the code!
