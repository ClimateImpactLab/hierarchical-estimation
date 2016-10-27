
R Files
-------

R files can have comment headers defined using unassigned strings:

.. code-block:: r

    "
    This is a valid r comment

    :author: My Name
    :contact: my_email@gmail.com
    "

or using a block of start-of-line comments: 

.. code-block:: r

    # This is also a valid comment
    # :author: My Name
    # :contact: my_email@gmail.com

Example R file
~~~~~~~~~~~~~~

A template for a good R script can be found in :download:`docs/examples/r.R <r.R>`, reproduced here:


.. literalinclude:: r.R
    :linenos:
    :language: r