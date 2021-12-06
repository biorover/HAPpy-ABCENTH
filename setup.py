from setuptools import setup

setup(
    name='HAPpy-ABCENTH',
    version='0.2.1',
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
)
