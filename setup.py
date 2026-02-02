from setuptools import setup, find_packages

setup(
    name='metalotl',
    use_scm_version=True,
    setup_requires=['setuptools-scm'],
    install_requires=[
        'shiny',
        'shinywidgets',
        'plotly',
        'glasbey',
        'scanpy',
        'pandas',
        'setuptools-scm'
    ],
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    python_requires='>=3.12',
)
