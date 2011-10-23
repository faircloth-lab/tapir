import distribute_setup
distribute_setup.use_setuptools()
from setuptools import setup, find_packages


if __name__ == '__main__':
    setup(
        name="picme",
        version="1.0",
        description="phylogenetic informativeness calculator for marker evaluation",
        author="Brant Faircloth, Jonathan Chang, Mike Alfaro",
        author_email="brant.faircloth+picme@gmail.com",
        url="http://github.com/faircloth-lab/picme/",
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
        requires=["dendropy","matplotlib","numpy","scipy","rpy2"],
        long_description=open("README.rst").read(),
        packages=[],
        scripts=["bin/picme_compute.py",
                "bin/picme_plot.py",
                "bin/picme_matplot.py"
                ],
        )

