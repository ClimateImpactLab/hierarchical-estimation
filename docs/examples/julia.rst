
Julia Files
------------

Julia files can have comment headers defined using block quotes:

.. code-block:: julia

    #=
    This is a valid julia comment block

    :author: My Name
    :contact: my_email@gmail.com
    =#

or using a block of start-of-line comments: 

.. code-block:: julia

    # This is also a valid comment
    # :author: My Name
    # :contact: my_email@gmail.com

Example Julia File
~~~~~~~~~~~~~~~~~~~

A template for a good julia script can be found in :download:`docs/examples/julia.py <julia.jl>`, reproduced here:

.. literalinclude:: julia.jl
    :linenos:
    :language: julia
