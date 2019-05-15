import setuptools

with open('molgif/_version.py', 'r') as fid:
    exec(fid.read())

setuptools.setup(name='molgif',
                 version=__version__,
                 author='Michael Cowan',
                 url='https://www.github.com/michael-cowan/molgif',
                 packages=['molgif'],
                 python_requires='>=2.7',
                 install_requires=['matplotlib',
                                   'ase'])
