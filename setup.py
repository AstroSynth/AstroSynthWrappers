from setuptools import setup
from os import path

HERE = path.abspath(path.dirname(__file__))

setup(name='AstroSynthWrappers',
      version='0.1.4',
      description='Data Wrappers meant to mimic the astroSynth data interface',
      url='https://github.com/tboudreaux/AstroSynthWrappers.git',
      author='Thomas Boudreaux',
      author_email='thomas@boudreauxmail.com',
      license='MIT',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3'
      ],
      install_requires=[
          'numpy>=1.14.1',
          'pandas>=0.20.3',
          'tqdm>=4.19.4',
          'scipy>=1.0.0',
          'astropy>=2.0.2',
          'pymongo>=3.5.1'
      ],
      packages=['AstroSynthWrappers', 'astroSynthWrappers.PTF'],
      zip_safe=False)