from setuptools import setup, find_packages

# I do not have a readme yet
#with open("README.md", "r") as fh:
#    long_description = fh.read()

setup(
    name='flux',
    version='-8.9.9',
    author='Lukas Finkbeiner, C. D. Nunhokee, Aaron Parsons',
    author_email='lfinkbeiner@berkeley.edu',
    url='https://github.com/3276908917/HERA',
    packages=find_packages(),
    package_data={'': ['gleam_with_alpha.txt', 'ant_dict.pk', 'J.npy']},
    include_package_data=True,
    install_requires=[], #I definitely need to come back and fix this
)
