from setuptools import setup, find_packages

setup(
    name='HAPpy-ABCENTH',
    version='1.0',
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
    include_package_data = True
)
