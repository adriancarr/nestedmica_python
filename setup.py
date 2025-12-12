from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension(
        "nestedmica.model.dp_likelihood",
        ["nestedmica/model/dp_likelihood.pyx"],
        include_dirs=[np.get_include()],
        extra_compile_args=["-O3", "-ffast-math"],
    ),
    Extension(
        "nestedmica.model.cython_model",
        ["nestedmica/model/cython_model.pyx"],
        include_dirs=[np.get_include()],
        extra_compile_args=["-O3", "-ffast-math"],
    ),
]

setup(
    name="nestedmica",
    version="1.2.0",
    description="High-performance Nested Sampling Motif Discovery",
    ext_modules=cythonize(
        extensions,
        compiler_directives={
            "language_level": "3",
            "boundscheck": False,
            "wraparound": False,
            "cdivision": True,
        }
    ),
)
