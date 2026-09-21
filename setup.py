import os

from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

compile_args = ["-O3", "-fopenmp", "-std=c++17"]
if os.environ.get("CHAAI_NATIVE", "1") != "0":
    compile_args.append("-march=native")

ext_modules = [
    Pybind11Extension(
        f"chaai._core.{name}",
        [f"src/core/{name}.cpp"],
        include_dirs=["src"],
        extra_compile_args=compile_args,
        extra_link_args=["-fopenmp"],
    )
    for name in ("genecall", "dist", "search")
]

setup(ext_modules=ext_modules, cmdclass={"build_ext": build_ext})
