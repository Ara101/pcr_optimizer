from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read() 

INSTALL_REQUIRES = []

def doSetup(install_requires):
    setup(
        name='PCR_protocal_optimizer',
        version='0.1',
        author=' Lily Torp and K. Lionel Tukei ',
        author_email='ltorp3@uw.edu',
        url='https://github.com/Ara101/PCR_Optimization_class.git',
        description='A funtion PCR protocal optimizations',
        long_description=long_description,
        long_description_content_type='text/markdown',
        packages=['PCR_Optimization'],
        package_dir={'PCR_Optimization':
            'PCR_Optimization'},
        install_requires=install_requires,
        include_package_data=True,
    )

if __name__ == '__main__':
  doSetup(INSTALL_REQUIRES)
