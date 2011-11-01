import distribute_setup
distribute_setup.use_setuptools()
from setuptools import setup, find_packages

import picme

if __name__ == '__main__':
    setup(
        name="picme",
        version = picme.__version__,
        description="phylogenetic informativeness calculator for marker evaluation (PICME)",
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
        requires=["dendropy","numpy(>=1.3)","scipy(>=0.9.0)"],
        long_description=open("README.rst").read(),
        scripts=["bin/picme_compute.py",
                "bin/picme_plot.py",
                "bin/picme_matplot.py"
            ],
        packages=[
            'picme',
            'picme.tests',
            ],
        package_data = {
            '':['*.txt'],
            'picme':['tests/test-data/*', 'data/*']
            },
        include_package_data = True,
        test_suite = 'picme'
    )

