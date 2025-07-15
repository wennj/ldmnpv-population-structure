import numpy as np
import gurobipy as gp
import logging
import csv

from typing import Tuple, Dict, List

from gurobipy import GRB


def convertArrayToString(input_array: np.ndarray) -> str:
    """Converts an array of size nx4 with 0-1-values and 
    exactly one 1 per row into the corresponding n-bit string
    of A,C,G,T.

    :param input_array: the array representation of a basis genotype
    :type input_array: np.ndarray
    :return: the string representation of the basis genotype
    :rtype: str
    """
    s = ''
    for i in range(len(input_array)):
        letter = np.argmax(input_array[i,:])
        if letter == 0:
            s += 'A'
        if letter == 1:
            s += 'C'
        if letter == 2:
            s += 'G'
        if letter == 3:
            s += 'T'
    return s


def convertStringToArray(string: str) -> np.ndarray:
    """Converts an n-bit string of A,C,G,T into the corresponding
    array representation with 0-1-values and exactly one 1 per row.

    :param string: the string representation of the basis genotype
    :type string: str
    :return: the array representation of the basis genotype
    :rtype: np.ndarray
    """
    array_out = np.zeros((len(string),4))
    j = 0
    for c in string:
        if c=='A':
            array_out[j,0] = 1
        if c=='C':
            array_out[j,1] = 1
        if c=='G':
            array_out[j,2] = 1
        if c=='T':
            array_out[j,3] = 1
        j += 1
    return array_out


def convertSetToArray(theset):
    numb_basis = len(theset)
    numb_pos = len(next(iter(theset)))
    Basis = np.zeros((numb_basis,numb_pos,4))
    idx_b = 0
    for s in theset:
        Basis[idx_b,:,:] = convertStringToArray(s)
        idx_b += 1
    return Basis


def readIsolates(path: str) -> Tuple[int, int, np.ndarray, Dict[str, np.ndarray]]:
    """Method to read a set of isolates from a given file.

    :param path: the file to read
    :type path: str
    :return: number n of SNP-positions per isolate, number m of isolates, 
        isolates as mxnx4 array, isolates as dictionary indexed by isolate name.
    :rtype: Tuple[int, int, np.ndarray, Dict[str, np.ndarray]]
    """
    lines = []
    with open(path, "r", encoding='utf-8-sig') as f:
        for line in f.readlines():
            lines.append(line)

    isolate_dict = {}
    for line in lines:
        splitted = line.split(";")
        splitted = [s.rstrip("\n") for s in splitted]
        splitted[0] = splitted[0].lstrip("\ufeff")
        try:
            isolate_dict[splitted[0]] += splitted[1:]
        except KeyError:
            isolate_dict[splitted[0]] = splitted[1:]
    numb_pos = int(len(isolate_dict[lines[0].split(";")[0].lstrip("\ufeff")])/4)
    for key in isolate_dict:
        isolate_dict[key] = np.reshape(np.array(isolate_dict[key],dtype="int64"), (numb_pos,4))
        isolate_dict[key] = np.divide(isolate_dict[key], np.reshape(np.sum(isolate_dict[key],1), (numb_pos,1)))
    numb_isolates = len(isolate_dict)
    isolates = np.zeros((numb_isolates,numb_pos,4))
    i = 0
    for key in isolate_dict:
        isolates[i,:,:] = isolate_dict[key]
        i += 1
    return numb_pos, numb_isolates, isolates, isolate_dict


def readBasisInput(path: str) -> Tuple[np.ndarray, int]:
    """Method to read the given basis genotypes.

    :param path: file to read the basis from
    :type path: str
    :return: The read basis as #basistypes x #SNP-positions x 4 array and number of basis types.
    :rtype: Tuple[np.ndarray, int]
    """
    if path == "":
        return 0, 0
    else: 
        lines = []
        with open(path, "r") as f:
            for line in f.readlines():
                lines.append(line)
        data = {}
        for line in lines:
            splitted = line.split(':')
            data[int(splitted[0])] = splitted[1].rstrip("\n")

        numb_pos = len(data[1])                             # number of SNP-positions
        numb_basis = len(data)                              # number of basistypes

        Basis = np.zeros((numb_basis,numb_pos,4))

        i = 0
        for key in data:
            string = data[key]
            j = 0
            for c in string:
                if c=='A':
                    Basis[i,j,0] = 1
                if c=='C':
                    Basis[i,j,1] = 1
                if c=='G':
                    Basis[i,j,2] = 1
                if c=='T':
                    Basis[i,j,3] = 1
                j += 1
            i += 1
        return Basis, numb_basis


def writeBasis(path: str, basis_set: set[str]):
    """Method to write a set of basis types to a file.

    :param path: the file to write the basis to
    :type path: str
    :param basis_set: set of basis types given as set of A,C,G,T-strings.
    :type basis_set: set[str]
    """
    with open(path, "w", newline='') as f:
        basis_idx = 1
        for basis_string in basis_set:
            f.write(f"{basis_idx}:"+basis_string+"\n")
            basis_idx += 1


def solveForGivenBasis(isolates_dict: Dict[str, np.ndarray], numb_basis: int, numb_pos: int, Basis: np.ndarray) -> Tuple[Dict[str, List[float]], Dict[str, np.ndarray], Dict[str, float]]:
    """Solve the LP formulation for all isolates and the given basis.

    :param isolates_dict: dictionary of isolates as numpy arrays indexed by isolate name
    :type isolates_dict: Dict[str, np.ndarray]
    :param numb_basis: number ob basis genotypes
    :type numb_basis: int
    :param numb_pos: number of SNP-positions
    :type numb_pos: int
    :param Basis: representation of the given basis genotypes as numb_basis x numb_pos x 4 array
    :type Basis: np.ndarray
    :param writeanalysepath: file to write the analysis to
    :return: Dictionary of the computed composition of each isolate, dictionary of error matrices and dictionary of largest error matrix entry
    :rtype: Tuple[Dict[str, List[float]], Dict[str, np.ndarray], Dict[str, float]]
    """
    composition = {}
    error_matrix = {}
    z_out = {}
    max_err = {}
    x_sum = {}

    for isolate in isolates_dict:

        logging.info(f"Calculation isolate {isolate}")

        sample_data = isolates_dict[isolate]

        mod = gp.Model("Virus_Variants")
        mod.setParam('LogToConsole',0)
        z = mod.addMVar(shape = (numb_pos, 4), vtype = GRB.CONTINUOUS, name = "z")
        f = mod.addMVar(lb = -float('inf'), shape = (numb_pos, 4), vtype = GRB.CONTINUOUS, name = "f")
        x = mod.addMVar(lb = 0.0, ub = 1.0, shape = numb_basis, vtype = GRB.CONTINUOUS, name = "x")
        mod.setObjective(z.sum(), GRB.MINIMIZE)
        for i in range(4):
            for j in range(numb_pos):
                mod.addConstr(-f[j,i] - z[j,i] <= 0, name = "linearization+"+str(i)+str(j))
                mod.addConstr(f[j,i] - z[j,i] <= 0, name = "linearization-"+str(i)+str(j))
                mod.addConstr( Basis[:,j,i] @ x + f[j,i] == sample_data[j,i] , name="Constr"+str(j)+str(i))
        mod.optimize()

        composition[isolate] = [round(a,4) for a in x.X]
        error_matrix[isolate] = np.array(f.X)
        z_out[isolate] = np.array(z.X)
        max_err[isolate] = np.max(z_out[isolate])
        x_sum[isolate] = sum(composition[isolate])

        logging.info(f"{composition[isolate]}")

        logging.info(f"sum of x: {sum(composition[isolate]):.3f}")
        logging.info(f"max error: {np.max(z_out[isolate]):.4f}")

    return composition, error_matrix, max_err


def writeOutput(composition: Dict[str, List[float]], errormatrix: Dict[str, np.ndarray], max_err: Dict[int, float], path: str, epsilon: float, header: List[str]):
    """Method to write the calculated compositions of the isolates

    :param composition: dictionary of computed compositions indexed by isolate name
    :type composition: Dict[str, List[float]]
    :param errormatrix: the computed error matrix
    :type errormatrix: Dict[str, np.ndarray]
    :param max_err: dictionary of largest entry of error matrix for each isolate
    :type max_err: Dict[int, float]
    :param path: path to write to
    :type path: str
    :param epsilon: epsilon, to consider error entries as large
    :type epsilon: float
    :param header: header for the file
    :type header: List[str]
    """
    with open(path, "w", newline='') as f:
        writer = csv.writer(f, delimiter=";")
        writer.writerow(header)
        for isolate in composition:
            large_Error_String = ""
            for i in range(np.shape(errormatrix[isolate])[0]):
                if max([abs(errormatrix[isolate][i,j]) for j in range(4)]) > epsilon:
                    large_Error_String += f"{i+1}: {round(max([abs(errormatrix[isolate][i,j]) for j in range(4)]),4)}, "
            line = [isolate] + composition[isolate] + [round(max_err[isolate], 4), large_Error_String]
            writer.writerow(line)


def writeInterpolation(composition: Dict[str, List[float]], Basis: np.ndarray, numb_pos: int, numb_basis: int, path: str):
    """Method to write an interpolation of the samples with the computed composition

    :param composition: the computed composition
    :type composition: Dict[str, List[float]]
    :param Basis: the given basis genotypes
    :type Basis: np.ndarray
    :param numb_pos: number of SNP-positions
    :type numb_pos: int
    :param numb_basis: number of basis genotypes
    :type numb_basis: int
    :param path: path to write the interpolation to
    :type path: str
    """
    with open(path, "w", newline='') as f:
        writer = csv.writer(f, delimiter = ";")
        writer.writerow(["isolate", "A rel", "C rel", "G rel", "T rel"])
        for isolate in composition:
            scale = sum(composition[isolate])
            interpolation = np.zeros((numb_pos,4))
            for i in range(numb_basis):
                interpolation += composition[isolate][i]*Basis[i,:,:]
            interpolation = interpolation/scale
            for j in range(numb_pos):
                writer.writerow([isolate,interpolation[j,0],interpolation[j,1],interpolation[j,2],interpolation[j,3]])


def writeErrorMatrix(path: str, errordict: Dict[str, np.ndarray]):
    """Method to write the error matrix to a file.

    :param path: the path to write the error matrix
    :type path: str
    :param errordict: error matrix given as dictionnary of n x 4 arrays indexed by isolate name.
    :type errordict: Dict[str, np.ndarray]
    """
    with open(path, "w", newline='') as f:
        writer = csv.writer(f, delimiter=";")
        writer.writerow(["isloate", "A", "C", "G", "T"])
        for isolate in errordict:
            errormatrix = errordict[isolate]
            for i in range(np.shape(errormatrix)[0]):
                writer.writerow([isolate, round(errormatrix[i,0],4), round(errormatrix[i,1],4), round(errormatrix[i,2],4), round(errormatrix[i,3],4)])
    