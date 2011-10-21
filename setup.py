import distribute_setup
distribute_setup.use_setuptools()
from setuptools import setup, find_packages


if __name__ == '__main__':
    setup(
        name='pd-ht',
        version="1.0",
        description="Parallel computation of phylogenetic informativeness",
        author="Brant Faircloth, Jonathan Chang, Mike Alfaro",
        author_email="brant.faircloth+pdht@gmail.com",
        url="http://github.com/faircloth-lab/pd-ht/",
        license="http://www.opensource.org/licenses/BSD-3-Clause",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Operating System :: OS Independent',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
             ],
        requires=['dendropy','matplotlib','numpy','scipy'],
        long_description=open('README.rst').read(),
        packages=[],
        scripts=['bin/estimate_p_i.py',
                'bin/create_plots.py',
                ],
        )

