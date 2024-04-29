# Makefile

# Define the default target. When you run "make" without specifying a target, it will run the first target in the file.
all: build_ext

# Define the build_ext target
build_ext:
	python setup.py build_ext -i

# Add a clean target for removing the build artifacts
clean:
	rm -rf build
	rm -rf *.so
	rm -rf *.pyd
	rm -rf *.c
	rm -rf *.cpp

.PHONY: all build_ext clean

