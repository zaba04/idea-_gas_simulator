# Makefile dla macOS i Linux
SOURCES = animateAr.cpp
CC = g++

# Wykrywanie systemu
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
    # macOS
    CFLAGS = -pedantic -Wall -O2 -Wno-deprecated-declarations
    CLINK = -framework GLUT -framework OpenGL -framework Cocoa
else
    # Linux
    CFLAGS = -pedantic -Wall -O2
    CLINK = -lglut -lGL -lGLU -lGLEW -lm
endif

all: $(SOURCES)
	@echo "🔨 Compiling Enhanced Argon Visualizer for $(UNAME_S)..."
	$(CC) $(CFLAGS) $(SOURCES) $(CLINK) -o animateAr
	@echo "✅ Build successful!"
	@echo "📦 Run with: ./animateAr params xyz.dat out.dat"
	
clean:
	@echo "🧹 Cleaning..."
	rm -f *.o animateAr
	@echo "✨ Done!"

run: all
	@echo "🚀 Running with default files..."
	./animateAr params xyz.dat out.dat

.PHONY: all clean run