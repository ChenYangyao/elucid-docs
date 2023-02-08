Documentation for ELUCID Project
================================

![ELUCID reconstructed density field](docs/imgs/CS-fields-and-SDSS-galaxies-small.jpg)

This repository hosts the documentation written for the ELUCID project (ses Wang et al. 2016, [arXiv:1608.01763](https://arxiv.org/abs/1608.01763) for a review).
The project description, data access description, data specification, and 
other related topics will be presented. 

The documentation is still under construction. If interested in either the 
data or the usage of it, please contact Yangyao Chen [yangyaochen.astro@foxmail.com](mailto:yangyaochen.astro@foxmail.com)) for details.


# Join Us

## How to Modify the documentation

Install a Sphinx environment. We prefer to use conda:

```bash
$ conda create --name sphinx
$ conda activate sphinx
$ conda install sphinx
$ pip install sphinx-rtd-theme
$ pip install sphinx-book-theme
```
This creates and enters a new conda environment and installs the dependencies 
``sphinx``, ``sphinx-rtd-theme`` and ``sphinx-book-theme``.

Clone/fork this repository, checkout to the ``dev`` branch, and edit the 
documentation files, e.g., ``index.rst`` or ``docs/data/index.rst``.

Compile the documentation by ``make html``. You can see the generated html files 
under ``_build/html`` directory.

Add and commit your changes (do not include ``_build`` and other private data), 
and push into the dev branch (included contributors), 
or create a pull request (external contributors).

The administrator will periodically merge the changes into the main branch.