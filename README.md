# AstroSynthWrappers
Eventually a set of data wrappers for different data sets that are meant to emulate the interface of astroSynth, so that all code written to work with astroSynth can also work with these data sets. Currently this project has only focued on building a wrapper around Palomar Transient Factory Data.

# NOTE
This repository, while publically avalible is not nessisarily intended for use by anyone other than me, that is all to say that this code will likely not work for whatever you may try and do. If if does great, all the better, but I wrote it with some pretty narrow goals and unusual assumtptions. While these were valid in my case they are most likely not valid in yours. Anyways thats all there, have a great day!

## Prerequisits
1) MongoDB installed
2) UNIX system
3) MongoDB should have PTF data stored in a collection

## Installation
If you really want to install this there are two ways. If you are new to python (which seams unlikely if you are trying to install this wierd module that you probably dug up after hours of searching for one very particular thing) then I recommend pip
```bash
pip install astroSynthWrappers
```
I will try to keep the pip listing up to date. However, If you want to be sure you have the most up to data version clone this repository and then run
```bash
python setup.py install
```
