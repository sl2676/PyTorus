from setuptools import setup, Extension
import pybind11

cpp_args = ['-std=c++17', '-stdlib=libc++', '-Wno-register', '-mmacosx-version-min=10.14', '-Wall']

ext_modules = [
    Extension(
        "PytTorus._torus",
        ["wrap.cpp", "./Torus/src/utils/WDMath.cc", "./Torus/src/utils/Compress.cc", "Torus/src/utils/Err.cc", "Torus/src/utils/CHB.cc"],
        include_dirs=[
            pybind11.get_include(),
            "./Torus/src/utils",  # Ensure correct include path
            "./Torus/include",  # Add if headers are in a separate include folder
        ],
        language="c++",
        extra_compile_args=cpp_args,
    ),
]

setup(
    name="PytTorus",
    version="0.0.2",
    author="Sean Ly",
    author_email="seanly1101@gmail.com",
    description="Torus library for orbital mechanics",
    ext_modules=ext_modules,
    packages=["PytTorus"],
    package_dir={"PytTorus": "PytTorus"},
    python_requires=">=3.7",
    install_requires=["pybind11"],
    package_data={"": ["Torus/src/utils/*.h", "Torus/src/utils/*.templates", "Torus/src/utils/*.cpp"]},  # Ensure headers are included in the package
    include_package_data=True,
)

