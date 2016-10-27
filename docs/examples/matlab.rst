
Matlab Files
------------

Matlab files can have comment headers defined using block comments:

.. code-block:: matlab

    %{
    This is a valid matlab comment

    :author: My Name
    :contact: my_email@gmail.com
    %}

or using a block of start-of-line comments: 

.. code-block:: matlab

    % This is also a valid comment
    % :author: My Name
    % :contact: my_email@gmail.com


Example Matlab File
~~~~~~~~~~~~~~~~~~~

A template for a good matlab script can be found in :download:`docs/examples/matlab.m <matlab.m>`, reproduced here:


.. literalinclude:: matlab.m
    :linenos:
    :language: matlab