import logging
import genotypesanalysis as gt


logging.getLogger().setLevel(logging.INFO)



####################################### Parameters ################################################
dataset = "original"
#dataset = "minimal_example"

isolatespath = "../datasets/"+dataset+"/input.csv"
basispath = "../datasets/"+dataset+"/basis.txt"

writeanalysispath = "../datasets/"+dataset+"/analysis.csv"
header_analysis = ["isolate ","A1 ", "A2 ", "A3 ", "A4 ", "B1 ", "B2 ", "B3 ", "C1 ", "maximal error ", "large error entris (position: entry)"]

# for dataset minimal example
header_analysis = ["isolate ","B1 ", "B2 ", "maximal error ", "large error entris (position: entry)"]
writeinterpolationpath = "../datasets/"+dataset+"/interpolation.csv"
writeerrormatrixpath = "../datasets/"+dataset+"/errormatrix.csv"

# All entries of error matrices larger than epsilon are listed explicitly in main main analysis file
epsilon = 0.1
###################################################################################################


####################################### Read Input ################################################
logging.info("Begin reading input data")

numb_pos, numb_isolates, isolates, isolate_dict = gt.readIsolates(isolatespath)
Basis_min, numb_basis_min = gt.readBasisInput(basispath)

logging.info("Reading input completed")
###################################################################################################


############################ Analyse distribution for given basistypes ############################
logging.info("Starting analysis distribution for given genotypes")

composition, error_matrix, max_error = gt.solveForGivenBasis(
    isolate_dict, 
    numb_basis_min, 
    numb_pos, 
    Basis_min
    )

logging.info("analysis distribution for given genotypes finished")
##################################################################################################


####################################### Write output files #######################################
logging.info("Begin writing output data")

gt.writeOutput(
    composition, 
    error_matrix, 
    max_error, 
    writeanalysispath, 
    epsilon,
    header_analysis
    )

gt.writeErrorMatrix(
    writeerrormatrixpath, 
    error_matrix
    )


gt.writeInterpolation(
    composition, 
    Basis_min, 
    numb_pos, 
    numb_basis_min, 
    writeinterpolationpath
    )

logging.info("Finished writing output data")
###################################################################################################
