# Genetic Structure of Baculovirus Populations Revealed by Amplicon-Based SNP Analysis

This repository accompanies the following publication:

-   **Oehlmann, C., Fan, J., Rihlmann, M., Lemme, H., Müller, J., Ruoff, B., Wennmann, J.T., Jehle, J.A. (submitted to Virus Evolution).**

### Background

The underlying mathematical model is briefly described in the PDF `report_LdMNPV.pdf` in the `doc/` folder. You can also click [here](https://github.com/wennj/ldmnpv-population-structure/blob/main/doc/report_LdMNPV.pdf) to open it directly.

### How to run the code

1.  Set up Gurobi: Go to the Gurobi website (www.gurobi.com), download and install Gurobi. Create an account and get yourself a license (free for academics). Validate your license to your machine with the grbgetkey command.

    There are various tutorials online on how to set up Gurobi for different operating systems.

2.  Navigate to `./src/` and run `python3 main.py`

3.  You have to select the dataset to use in the first lines of `./src/main.py`

The dataset minimal_example provides an easy to understand example which can also be found in the documentation. The dataset original contains the original data. By execution the program, you should be able to reproduce the files in the subdirectory `results`.
