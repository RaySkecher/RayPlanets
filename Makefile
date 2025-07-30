CC = gcc
CFLAGS = -std=c99 -Wall -Wextra -O3 -march=native -ffast-math
TARGET = build/raytracer1
SOURCE = raytracetemp.c

TOTAL_TIME = $(shell echo "scale=8; 3 * 3.141592653589793" | bc -l)
TIME_STEP = $(shell echo "scale=8; $(TOTAL_TIME) / 12" | bc -l)
GIF_DELAY = 20

.PHONY: all clean gif

all: $(TARGET)
	@mkdir -p build
	@for i in $$(seq 0 18); do \
		frame_num=$$(printf "%02d" $$(($$i+1))); \
		echo "rendering..$$frame_num"; \
		time=$$(echo "scale=8; $$i * $(TIME_STEP)" | bc -l); \
		./$(TARGET) $$time > build/frame$$frame_num.ppm; \
	done

$(TARGET): $(SOURCE)
	@mkdir -p build
	$(CC) $(CFLAGS) -o $(TARGET) $(SOURCE) -lm

gif: all
	@convert -delay $(GIF_DELAY) build/frame*.ppm -loop 0 build/animation.gif
	@if command -v open >/dev/null 2>&1; then \
		open build/animation.gif; \
	elif command -v xdg-open >/dev/null 2>&1; then \
		xdg-open build/animation.gif; \
	elif command -v start >/dev/null 2>&1; then \
		start build/animation.gif; \
	else \
		echo "GIF created at build/animation.gif"; \
	fi

clean:
	rm -rf build/
