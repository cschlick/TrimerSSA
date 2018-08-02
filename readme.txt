There is no way to run this yet without editing options in the code. I either edit and run main.py directly from in IDE, or use the ipython notebook.

The only uncommon dependencies are transformations.py and my PDBModule. Both of these could be replaced at some point.



There is sometimes a bug then ends with:

    assert (len(adding_edge.trimers) == 1)
AssertionError

It happens randomly and I haven't investigated.