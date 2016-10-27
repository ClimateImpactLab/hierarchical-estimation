
Working with SphinxScript documentation
---------------------------------------

`SphinxScript <https://github.com/ClimateImpactLab/sphinxscript>`_ uses 
`sphinx <http://www.sphinx-doc.org/>`_ to automatically generate documentation 
sites from `ReStructured Text (rst) <http://docutils.sourceforge.net/docs/user/rst/quickstart.html>`_ 
formatted headers in the scripts found in this directory. See the 
`SphinxScript Documentation <hhttp://sphinxscript.readthedocs.io/en/latest/>`_ 
for a list of available language parsers and tips on usage.


Documentation Best Practices at the Climate Impact Lab
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All files should have the following elements:

* A one-liner description giving someone a quick sense of the file's purpose
* Author
* Contact info (email)
* Last modified date
* Team
* Team Lead
* Input variables
* Output variables
* Source information
* A long, specific description, giving usage notes, caveats/warnings, source info, etc.

An example header (using stata syntax) is provided here:

.. literalinclude:: examples/header.rst
    :linenos:


Examples
~~~~~~~~~

.. toctree::
    :maxdepth: 1

    examples/stata
    examples/r
    examples/matlab
    examples/python
    examples/julia