Because the Assimulo dependency is difficult to build on third-party content providers (travis, RTD, etc), I've opted to manually control the examples in this directory. Here are the steps to keep the examples up to date:

1. Run each notebook independently, so that the results (figures, etc) it includes are pre-rendered and self-contained
2. Execute the supplied Makefile to convert ipynb -> rst and generate hte rendered figures
3. Commit the newly generated *_files/ folders and all new resources

In the future, this can be simplified by packaging the entire parcel model install. However, that will require either (a) re-opening the ability to use third-party ODE solvers from SciPy, or (b) packaging Assimulo and Sundials into conda for easy access on a third-party server.