To run the code do the following:

Set up Gurobi:
- Go to the Gurobi website (www.gurobi.com), download and install Gurobi
- Create an account and get yourself a license (free for academics)
- Validate your license to your machine with the grbgetkey command

There are various tutorials online on how to set up Gurobi for different operating systems.

Navigate to ./src/ and run "python3 main.py"

You have to select the dataset to use in the first lines of ./src/main.py

The underlying mathematical model is briefly described in the pdf file doc/report_LdMNPV.pdf

The dataset minimal_example provides an easy to understand example which can also be found in the documentation. The dataset original contains the original data. By execution the program, you should be able to reproduce the files in the subdirectory "results".

