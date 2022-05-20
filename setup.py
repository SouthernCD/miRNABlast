# coding utf8
import setuptools
from mirnablast.versions import get_versions
import platform
import sys

# system info
if platform.system() != "Linux":
    sys.stdout.write("Error: currently miRNABlast is not supported for " + platform.system() + "! ")
    exit()

setuptools.setup(
    name="miRNABlast",
    version=get_versions(),
    author="Yuxing Xu",
    author_email="xuyuxing@mail.kib.ac.cn",
    license='MIT',
    description="Use mature query miRNA sequences to BLAST the genome, fold the flanking sequences of the subjects to find possible homologous miRNA locus",
    url="https://github.com/SouthernCD/miRNABlast",

    packages=setuptools.find_packages(),

    install_requires=[
        "ToolBiox>=0.0.1",
    ],

    python_requires='>=3.5',

    entry_points={'console_scripts': ['miRNABlast = mirnablast.cli:main']},

    classifiers=[
            "Development Status :: 5 - Production/Stable",
            "Environment :: Console",
            "Intended Audience :: Science/Research",
            "Natural Language :: English",
            "Operating System :: Unix",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: Implementation :: PyPy",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]

)
