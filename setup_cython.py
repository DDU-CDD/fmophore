from setuptools import setup, Extension
from Cython.Build import cythonize
import os

# List of Python modules to compile
modules = [
    "FMOPhore/FMOPhore_run.py",
    "FMOPhore/modules/FMOPhore_merge.py",
    "FMOPhore/modules/FMOPhore_download_pdb.py",
    "FMOPhore/modules/FMOPhore_prep.py",
    "FMOPhore/modules/FMOPhore_chain_correction.py",
    "FMOPhore/modules/FMOPhore_align.py",
    "FMOPhore/modules/FMOPhore_sep_LIGs.py",
    "FMOPhore/modules/FMOPhore_cutoff.py",
    "FMOPhore/modules/FMOPhore_split.py",
    "FMOPhore/modules/FMOPhore_ph4.py",
    "FMOPhore/modules/FMOPhore_Fragmention_com.py",
    "FMOPhore/modules/FMOPhore_complex_analysis.py",
    "FMOPhore/modules/FMOPhore_analysis.py",
    "FMOPhore/modules/FMOPhore_library_analysis.py",
    "FMOPhore/QM_run/run_DFTB.py",
    "FMOPhore/QM_run/run_MP2.py",
]

# Create a list of Extension objects
extensions = [
    Extension(
        module.replace("/", ".").replace(".py", ""),
        [module],
    )
    for module in modules
]

setup(
    name="FMOPhore",
    ext_modules=cythonize(extensions, compiler_directives={"language_level": "3"}),
)