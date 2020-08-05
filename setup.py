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
                 packages=setuptools.find_packages(),
                 entry_points={
                    'console_scripts': ['molgif=molgif.command_line:cli'],
                 },
                 python_requires='>=3.5',
                 install_requires=['matplotlib',
                                   'numpy>=1.17.2',
                                   'pillow',
                                   'ase>=3.17.0',
                                   'click>=7'])
