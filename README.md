# README for ANN VLE Solver

# About VLE Solvers
Python based VLE solvers (VTFlash, UVFlash) which use the Peng Robinson equation of state have been added. Species Data can be added to the SpeciesData.py file.

# About ANN C++ Coupling
The files required to couple a trained ANN with a given C++ file are provided in the ANN directory. To run the given test case in C++_testing use the following commands - 
cd ANN/C++_testing
make
./runANN

# Citation
If using the solvers, please cite the following two journal articles: [ANN aided VLE](https://doi.org/10.1063/5.0219323) and  [ISAT aided VLE](https://doi.org/10.1016/j.jcp.2024.112752)

For BibTex users:
```bibtex
@article{Srinivasan_ANN_VLE,
    author = {Srinivasan, Navneeth and Yang, Suo},
    title = "{Artificial neural network aided vapor–liquid equilibrium model for multi-component high-pressure transcritical flows with phase change}",
    journal = {Physics of Fluids},
    volume = {36},
    pages = {083328},
    year = {2024},
    doi = {https://doi.org/10.1063/5.0219323}
}

@article{Zhang_ISAT_VLE,
    author = {Hongyuan Zhang and Navneeth Srinivasan and Suo Yang},
    title = {In situ adaptive tabulation of vapor-liquid equilibrium solutions for multi-component high-pressure transcritical flows with phase change},
    journal = {Journal of Computational Physics},
    volume = {500},
    pages = {112752},
    year = {2024},
    doi = {https://doi.org/10.1016/j.jcp.2024.112752}
}
