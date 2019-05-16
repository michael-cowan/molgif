import setuptools

with open('molgif/_version.py', 'r') as fid:
    exec(fid.read())

with open('README.md', 'r') as readme:
    # ignore gifs
    description = ''.join([i for i in readme.readlines()
                           if not i.startswith('![')])

setuptools.setup(name='molgif',
                 version=__version__,
                 author='Michael Cowan',
                 url='https://www.github.com/michael-cowan/molgif',
                 description="creates smooth gifs of rotating molecules",
                 long_description=description,
                 long_description_content_type='text/markdown',
                 packages=['molgif'],
                 python_requires='>=2.7',
                 install_requires=['matplotlib',
                                   'ase'])
