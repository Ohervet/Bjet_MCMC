from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

bjet_core_ext_module = Extension(
    "bjet_core._bj_core",
    sources=[
        "bjet_core/bj_core.i",
        "bjet_core/bj_core.cpp",
        "bjet_core/processes_supp_core.cpp",
    ],
    swig_opts=[
        "-c++",
        "-py3",
        "-outdir",
        "bjet_core",
    ],
)

setup(
    cmdclass={"build_ext": build_ext},
    ext_modules=[bjet_core_ext_module],
    packages=[
        "bjet_core",
        "bjet_mcmc",
    ],  # these are the names of the packages to import in a python session.
    include_package_data=True,
)
