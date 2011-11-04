import distribute_setup
distribute_setup.use_setuptools()
from setuptools import setup, find_packages

import tapir

if __name__ == '__main__':
    setup(
        name="tapir",
        version = tapir.__version__,
        description="phylogenetic informativeness calculator for marker evaluation (tapir)",
        author="Brant Faircloth, Jonathan Chang, Mike Alfaro",
        author_email="brant.faircloth+tapir@gmail.com",
        url="http://github.com/faircloth-lab/tapir/",
        license="http://www.opensource.org/licenses/BSD-3-Clause",
        classifiers=[
            "Development Status :: 4 - Beta",
            "Environment :: Console",
            "Operating System :: OS Independent",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: BSD License",
            "Programming Language :: Python",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
             ],
        requires=["dendropy","numpy(>=1.3)","scipy(>=0.9.0)"],
        long_description=open("README.rst").read(),
        scripts=["bin/tapir_compute.py",
                "bin/tapir_plot.py",
                "bin/tapir_matplot.py",
                "bin/tapir_compare.py"
            ],
        packages=[
            'tapir',
            'tapir.tests',
            ],
        package_data = {
            '':['*.txt'],
            'tapir':['tests/test-data/*', 'data/*']
            },
        include_package_data = True,
        test_suite = 'tapir'
    )

