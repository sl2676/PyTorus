# Define the default target
all: build_ext

# Define the build_ext target
build_ext:
	CFLAGS=`python -m pybind11 --includes` python setup.py build_ext -i

# Add a clean target for removing the build artifacts
clean:
	rm -rf build
	rm -rf *.so
	rm -rf *.pyd
	rm -rf *.c
	rm -rf *.cpp

.PHONY: all build_ext clean

