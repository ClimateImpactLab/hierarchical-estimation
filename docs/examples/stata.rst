
Stata Files
-----------

Stata files can have comment headers defined using block comments:

.. code-block:: stata

    /*
    This is a valid stata comment

    :author: My Name
    :contact: my_email@gmail.com
    */

or using a block of start-of-line comments: 

.. code-block:: stata

    * This is also a valid comment
    * :author: My Name
    * :contact: my_email@gmail.com


Example Stata File
~~~~~~~~~~~~~~~~~~

A template for a good stata script can be found in :download:`docs/examples/stata.do <stata.do>`, reproduced here:


.. literalinclude:: stata.do
    :linenos:
    :language: stata