# RayPlanets Raytracer

### One-command build and run (recommended)

```bash
make
```

### Step-by-step commands

```bash
# Clean build artifacts
make clean
```

### Manual compilation

```bash
# Compile
gcc -std=c99 -Wall -Wextra -O3 -march=native -ffast-math -o raytracer raytracer.c -lm

# Render
./raytracer <time> > output.ppm

# View (macOS)
open output.ppm
```




