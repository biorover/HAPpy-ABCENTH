from setuptools import setup, find_packages
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='HAPpy-ABCENTH',
    version='1.0.1',
    packages=[
        'HAPpy',
        'ABCENTH'
    ],
    install_requires=[
        'ete3',
        'numpy',
        'intervaltree'
    ],
    entry_points = {
        'console_scripts': [
            'thammerin = HAPpy.thammerin:main',
            'ABCENTH = ABCENTH.__main__:main',
            'HAPpy = HAPpy.__main__:main'
        ]
    },
    include_package_data = True,
    long_description=long_description,
    long_description_content_type='text/markdown'
)
