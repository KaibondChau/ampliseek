import ampliseek
from setuptools import setup

setup(
    name='AmpliSeek',
    version=ampliseek.__version__,
    author=" Kevin Chau",
    author_email=" kevin.chau@ndm.ox.ac.uk",
    description=' AmpliSeq amplicon mapper tool',
    url=' https://github.com/KaibondChau/ampliseek',
    license='LICENSE',
    python_requires='>=3.6',
    packages=['ampliseek'],
    install_requires=['biopython'],
    entry_points={'console_scripts':['ampliseek=ampliseek.ampliseek:main']},
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English'])