from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

bjet_core_ext_module = Extension(
    "bjet_core._bj_core",
    sources=[
        "bjet_core/bj_core.i",
        "bjet_core/bj_core.cpp",
        "bjet_core/processes_supp_core.cpp",
    ],
    swig_opts=["-c++"],
)
setup(
    cmdclass={"build_ext": build_ext},
    ext_modules=[bjet_core_ext_module],
    # py_modules=["bj_core"],
    packages=["_bj_core", "blazar_mcmc"],
    # package_dir={"bjet_core": "bjet_core", "blazar_mcmc": "blazar_mcmc"},
    include_package_data=True,
)