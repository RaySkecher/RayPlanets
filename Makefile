CC = gcc
CFLAGS = -std=c99 -Wall -Wextra -O3 -march=native -ffast-math
TARGET = raytracer
SOURCE = raytracer.c

.PHONY: all clean gif

all: $(TARGET)
	./$(TARGET) 0.0 > frame1.ppm
	./$(TARGET) 2.0 > frame2.ppm
	./$(TARGET) 4.0 > frame3.ppm

$(TARGET): $(SOURCE)
	$(CC) $(CFLAGS) -o $(TARGET) $(SOURCE) -lm

gif: all
	convert -delay 50 frame1.ppm frame2.ppm frame3.ppm -loop 0 animation.gif

clean:
	rm -f $(TARGET) frame1.ppm frame2.ppm frame3.ppm animation.gif
