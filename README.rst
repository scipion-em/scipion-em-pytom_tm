=============================
Scipion plugin for PYTOM_TM
=============================

.. image:: https://img.shields.io/pypi/v/scipion-em-pytom_tm.svg
        :target: https://pypi.python.org/pypi/scipion-em-pytom_tm
        :alt: PyPI release

.. image:: https://img.shields.io/pypi/l/scipion-em-pytom_tm.svg
        :target: https://pypi.python.org/pypi/scipion-em-pytom_tm
        :alt: License

.. image:: https://img.shields.io/pypi/pyversions/scipion-em-pytom_tm.svg
        :target: https://pypi.python.org/pypi/scipion-em-pytom_tm
        :alt: Supported Python versions

.. image:: https://img.shields.io/pypi/dm/scipion-em-pytom_tm
        :target: https://pypi.python.org/pypi/scipion-em-pytom_tm
        :alt: Downloads

This plugin provides a wrapper around the program `pytom-match-pick <https://sbc-utrecht.github.io/pytom-match-pick/>`_
(GPU-accelerated template matching for cryo-electron tomography, part of the PyTom suite) to use it within
`Scipion <https://scipion-em.github.io/docs/release-3.0.0/index.html>`_ framework.

Installation
------------

You will need to use `3.0 <https://scipion-em.github.io/docs/release-3.0.0/docs/scipion-modes/how-to-install.html>`_
version of Scipion to run these protocols. To install the plugin, you have two options:


a) Stable version:

.. code-block::

    scipion3 installp -p scipion-em-pytom_tm

b) Developer's version

    * download the repository from github:

    .. code-block::

        git clone -b devel https://github.com/scipion-em/scipion-em-pytom_tm.git

    * install:

    .. code-block::

        scipion3 installp -p /path/to/scipion-em-pytom_tm --devel

To check the installation, simply run the following Scipion test for the plugin:

    .. code-block::

        scipion3 tests pytom_tm.tests.tests_pytom_tm.TestPytomTM

To check the installation, simply run one of the tests. A complete list of tests can be displayed by executing

    .. code-block::

        scipion3 tests --grep pytom_tm --show

Protocols
-----------

* **Template matching** : Generates score and angular maps from the introduced tomograms using a reference template and mask, with GPU-accelerated exhaustive angular search.
* **Extract candidates** : Extracts particle coordinates (and orientations) from the score maps produced by template matching with pytom-match-pick.

Latest plugin versions
----------------------

If you want to check the latest version and release history go to `CHANGES <https://github.com/scipion-em/scipion-em-pytom_tm/blob/master/CHANGES.txt>`_

References
----------

1. Chaillet, M. L., Roet, S., Veltkamp, R. C., & Förster, F. (2025). pytom-match-pick: A tophat-transform constraint for automated classification in template matching. Journal of Structural Biology: X, 11, 100125. https://doi.org/10.1016/j.yjsbx.2025.100125
2. Chaillet, M. L., van der Schot, G., Gubins, I., Roet, S., Veltkamp, R. C., & Förster, F. (2023). Extensive angular sampling enables the sensitive localization of macromolecules in electron tomograms. International Journal of Molecular Sciences, 24(17), 13375. https://doi.org/10.3390/ijms241713375