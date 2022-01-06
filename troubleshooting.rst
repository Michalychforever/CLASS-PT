Troubleshooting
================

Listed below are common problems found when installing CLASS-PT and their resolutions.

Segmentation Faults when running ``classy``
-------------------------------------------

We have found that the OpenBLAS library conflicts with the library Intel MKL which is used in ``numpy`` version 1.16 and higher on some machines. This incompatibility makes ``classy`` crash with “segmentation fault” even though the code can be executed with a ``.ini`` file without any errors.  If this is the case on the user’s computer, an easy fix is to use a ``numpy`` version below 1.16.  We plan to resolve this issue in future releases.

Segmentation Faults when running MontePython with CLASS-PT
------------------------------------------------------------

When using some newer versions of Python, segmentation faults can occur when interfacing CLASS-PT with MontePython, even though there are no issues when ``classy`` is run on its own. This is related to the above problem, a sourced by an internal conflict of ``classy`` and the latest ``scipy`` modules. This can be straightforwardly fixed by ensuring that ``classy`` is imported before ``scipy`` in MontePython. Practically, the line ``from classy import Class`` should be inserted at the top of the ``montepython/Montepython.py`` script.
