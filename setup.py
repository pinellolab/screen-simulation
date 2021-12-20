from setuptools import setup

setup(
    name='screen-simulation',
    version='0.0.1',    
    description='Simulate bulk CRISPR screen with sorting schemes and reproter/target edits',
    url='https://github.com/pinellolab/screen-simulation.git',
    author=['Jayoung Ryu', 'Basheer Becerra', 'Lucas Ferreira'],
    author_email=['jayoung_ryu@g.harvard.edu', 'basheer_becerra@dfci.harvard.edu', 'lferreiradasilva@mgh.harvard.edu'],
    license='MIT',
    packages=['screen_simulation'],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research', 
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 3',
    ],
)
