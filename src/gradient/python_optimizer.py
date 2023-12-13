#!user/bin/python
from mpi4py import MPI
import numpy as np
from scipy.optimize import minimize
import sys

# emlation of enums
def enum(**keys):
    return type('ENUM',(), keys)


def generate_communicators():
    global globalComm, externalServerRank, pythonComm

    # COMM_WORLD contains the total procs, called by python and external program(s)
    globalComm = MPI.COMM_WORLD
    assert globalComm.size>1, "At least two procs are required"
    assert globalComm.rank!=0, "Python process is not master"

    # Split external comm from globalComm
    externalComm = globalComm.Split(color=MPI.UNDEFINED, key=0)
    assert externalComm == MPI.COMM_NULL, "Python process is not participateed with the external program"

    # Split Python comm from globalComm
    externalServerRank = 0 # Rank of the external proc in the pythonComm
    pythonCommRank = 1 # Rank of the python proc in the pythonComm
    pythonComm = globalComm.Split(color=0, key=pythonCommRank)
    assert pythonComm.size == 2, "Exactly 2 procs in pythonComm is required!"
    assert pythonComm.rank == pythonCommRank, "Python process must NOT be the master process!"    

    return


def get_options(outputFile):
    global nOptions, options, optMethod, tolerance, minF, nDim, xBnd, x0, nCons, consType

    # Optimization options, assigned the same numbers with the external server
    iOption = enum(LIBRARY=0, OPT_METHOD=1, MAX_ITR=2, TOLERANCE=3, IS_MIN=4, DISPLAY=5)

    # True/False Dictionary
    numbers = list(range(0,2))
    false_true = ['False', 'True']
    num2logic = dict(zip(numbers,false_true))

    # Optimizer method dictionary, assigned the same numbers with the external server 
    numbers = list(range(0,6))
    optMethodName = ['CG', 'NEWTON-CG', 'TNC', 'BFGS', 'L-BFGS-B', 'SLSQP']
    optMethodDic = dict(zip(numbers,optMethodName))

    status = MPI.Status()

    nOptions = np.array([0], dtype='i') # number of options that is requested by external program
    # Get nOptions from external server
    pythonComm.Recv(nOptions, source=externalServerRank, tag=MPI.ANY_TAG, status=status)
    # Check if a stop tag was sent by external server
    if status.tag == TAG.STOP:
        raise RuntimeError("Unexpected stop was sent by external server!")
    assert status.tag == TAG.ARGSIZE, "nOptions must have been sent by tag ARGSIZE, unexpected tag= "+str(status.tag)+" was received!"

    options = np.empty(nOptions[0]) # options that is requested by external server
    # Get options from external server
    pythonComm.Recv(options, source=externalServerRank, tag=MPI.ANY_TAG, status=status)
    if status.tag == TAG.STOP:
        raise RuntimeError("Unexpected stop was sent by external server!")
    assert status.tag == TAG.ARG, "Options must have been sent by tag ARG, unexpected tag= "+str(status.tag)+" was received!"
    # Send DONE tag to external server
    pythonComm.Send(nOptions, dest=externalServerRank, tag=TAG.DONE)

    # Translate optimization Method
    optMethod=optMethodDic.get(int(options[iOption.OPT_METHOD]))
    assert optMethod != None, "Unknown optimzer method code: "+str(int(options[iOption.OPT_METHOD]))

    # Translate min/max
    minF=num2logic.get(int(options[iOption.IS_MIN]))
    assert minF != None, "Unknown min/max option code: "+str(int(options[iOption.IS_MIN]))

    # Translate display option
    #display=num2logic.get(int(options[iOption.DISPLAY]))
    #assert display != None, "Unknown display option code: "+str(int(options[iOption.OPT_METHOD]))
    display = int(options[iOption.DISPLAY])
    assert display >= 0 and display <= 1, "Unknown display option number: "+str(int(options[iOption.OPT_METHOD]))

    # Save tolerance and maximum number of iterations
    tolerance = options[iOption.TOLERANCE]
    maxItr    = int(options[iOption.MAX_ITR])

    nDim = np.array([0], dtype='i') # number of parameter dimensions
    # Get nDim from external server
    pythonComm.Recv(nDim, source=externalServerRank, tag=MPI.ANY_TAG, status=status)
    if status.tag == TAG.STOP:
        raise RuntimeError("Unexpected stop was sent by external server!")
    assert status.tag == TAG.ARGSIZE, "nDim must have been sent by tag ARGSIZE, unexpected tag= "+str(status.tag)+" was received!"

    xVector = np.empty(nDim[0]*3)
    # Get xVector from external server
    pythonComm.Recv(xVector, source=externalServerRank, tag=MPI.ANY_TAG, status=status)
    if status.tag == TAG.STOP:
        raise RuntimeError("Unexpected stop was sent by external server!")
    assert status.tag == TAG.ARG, "xVector must have been sent by tag ARG, unexpected tag= "+str(status.tag)+" was received!"

    # Send DONE tag to external server
    pythonComm.Send(nDim, dest=externalServerRank, tag=TAG.DONE)

    # initial value of parameters
    x0   = np.empty(nDim[0])
    x0   = xVector[0:nDim[0]]
    # Bound of parameters
    xBnd = [(None, None)] * nDim[0]
    for i in range(0,nDim[0]):
        xBnd[i] = (xVector[nDim[0] + i], xVector[2*nDim[0]+ i])

    nCons = np.array([0], dtype='i') # number of constraints
    # Get nCons from external server
    pythonComm.Recv(nCons, source=externalServerRank, tag=MPI.ANY_TAG, status=status)
    if status.tag == TAG.STOP:
        raise RuntimeError("Unexpected stop was sent by external server!")
    assert status.tag == TAG.ARGSIZE, "nCons must have been sent by tag ARGSIZE, unexpected tag= "+str(status.tag)+" was received!"

    consTypeBool = np.empty(nCons[0], dtype='i') # Constraint Type in boolian
    # Get consType from external server
    pythonComm.Recv(consTypeBool, source=externalServerRank, tag=MPI.ANY_TAG, status=status)
    if status.tag == TAG.STOP:
        raise RuntimeError("Unexpected stop was sent by external server!")
    assert status.tag == TAG.ARG, "consType must have been sent by tag ARG, unexpected tag= "+str(status.tag)+" was received!"

    # Send DONE tag to external server
    pythonComm.Send(nCons, dest=externalServerRank, tag=TAG.DONE)

    # Equality Dictionary
    numbers = list(range(0,2))
    equalityName = ['ineq', 'eq']
    equalityDict = dict(zip(numbers,equalityName))

    consType = np.empty(nCons[0], dtype=type(equalityName)) # Constraint Type
    for i in range(0, nCons[0]):
        consType[i] = equalityDict.get(consTypeBool[i])
        assert consType[i] != None, "Unknown constraint type number: "+str(consTypeBool[i])
        assert consType[i] == 'eq', "Not yet implemented for inequality constraints!"


    outputFile.write("\nOptimization method: "+optMethod+"\n")
    outputFile.write("\nMax # iterations: "+str(maxItr)+"\n")
    outputFile.write("\nOptimization tolerance: "+str(tolerance)+"\n")
    outputFile.write("\nMinimize: "+minF+"\n")
    outputFile.write("\nDisplay: "+str(display)+"\n")
    outputFile.write("\nnDim: "+str(nDim[0])+"  xBnd: "+str(xBnd)+"\n")
    outputFile.write("\nnDim: "+str(nDim[0])+"  x0: "+str(x0)+"\n")
    outputFile.write("\nnCons: "+str(nCons[0])+"  consType: "+str(consType)+"\n")
    outputFile.flush()

    return {'maxiter':maxItr, 'disp': display}


# Write some information/errors in an outhput file
outputFile = open('python.out', 'w')

def Main():
    global TAG, Library, OptMethod, iOption

    # Communication tags, assigned the same numbers with the external server
    TAG = enum(DONE=0, STOP=1, ARGSIZE=2, ARG=3, CALL_OBJ=4, RESULT_OBJ=5, CALL_OBJ_GRAD=6,
                   RESULT_OBJ_GRAD=7, CALL_CONS=8, RESULT_CONS=9, CALL_CONS_GRAD=10,
                   RESULT_CONS_GRAD=11, CALLBACK = 12)

    # Python Library, assigned the same numbers with the external server
    Library = enum(SCIPY = 0)

    generate_communicators()
    outputFile.write("\nrank/size of globalComm :"+str(globalComm.rank)+"/"+str(globalComm.size)+"\n")
    outputFile.write("\nrank/size of pythonComm :"+str(pythonComm.rank)+"/"+str(pythonComm.size)+"\n")
    outputFile.flush()
    
    options = get_options(outputFile)

    # Calculate objective function by calling the external server
    def call_fObj(x):
        status = MPI.Status()
        pythonComm.Send(x, dest=externalServerRank, tag=TAG.CALL_OBJ)
        f = np.empty(1)
        pythonComm.Recv(f, source=externalServerRank, tag=MPI.ANY_TAG, status=status)
        if status.tag == TAG.STOP:
            raise RuntimeError("Unexpected stop was sent by external server!")
        assert status.tag == TAG.RESULT_OBJ, "call_fObj: Expected tag = RESULT_OBJ, got tag = "+str(status.tag)+"!"

        return f[0]


    # Calculate gradient of objective function by calling the external server
    def call_fObj_grad(x):
        status = MPI.Status()
        pythonComm.Send(x, dest=externalServerRank, tag=TAG.CALL_OBJ_GRAD)
        gradF = np.empty(nDim)
        pythonComm.Recv(gradF, source=externalServerRank, tag=MPI.ANY_TAG, status=status)
        if status.tag == TAG.STOP:
            raise RuntimeError("Unexpected stop was sent by external server!")
        assert status.tag == TAG.RESULT_OBJ_GRAD, "call_fObj_grad: Expected tag = RESULT_OBJ_GRAD, got tag = "+str(status.tag)+"!"

        return gradF

    # Calculate constraint by calling the external server
    def call_fCons(x):
        status = MPI.Status()
        # Send x and receive constraint function
        pythonComm.Send(x, dest=externalServerRank, tag=TAG.CALL_CONS)
        f = np.empty(nCons[0])
        pythonComm.Recv(f, source=externalServerRank, tag=MPI.ANY_TAG, status=status)
        if status.tag == TAG.STOP:
            raise RuntimeError("Unexpected stop was sent by external server!")
        assert status.tag == TAG.RESULT_CONS, "call_fCons: Expected tag = RESULT_CONS, got tag = "+str(status.tag)+"!"
        return f

    # Calculate constraint gradient by calling the external server
    def call_fCons_grad(x):
        status = MPI.Status()
        pythonComm.Send(x, dest=externalServerRank, tag=TAG.CALL_CONS_GRAD)
        gradFVec = np.empty([nCons[0]*nDim[0]])
        pythonComm.Recv(gradFVec, source=externalServerRank, tag=MPI.ANY_TAG, status=status)
        if status.tag == TAG.STOP:
            raise RuntimeError("Unexpected stop was sent by external server!")
        assert status.tag == TAG.RESULT_CONS_GRAD, "call_back: Expected tag = RESULT_CONS_GRAD, got tag = "+str(status.tag)+"!"

        # Substitude gradFVec into gradF
        gradF = np.empty([nCons[0],nDim[0]])
        for i in range(0, nCons[0]):
            for j in range(0, nDim[0]):
                index = i*nDim[0] + j
                gradF[i,j] = gradFVec[index]

        return gradF

    # Storing the major iteration information
    def callback(x):
        status = MPI.Status()
        pythonComm.Send(x, dest=externalServerRank, tag=TAG.CALLBACK)
        f = np.empty(1)
        pythonComm.Recv(f, source=externalServerRank, tag=MPI.ANY_TAG, status=status)
        if status.tag == TAG.STOP:
            raise RuntimeError("Unexpected stop was sent by external server!")
        assert status.tag == TAG.RESULT_OBJ, "call_fObj: Expected tag = RESULT_OBJ, got tag = "+str(status.tag)+"!"

    while True:
        status = MPI.Status()
        temp = np.array([0],dtype='i')
        pythonComm.Recv(temp, source=externalServerRank, tag=MPI.ANY_TAG, status=status)
        if status.tag == TAG.STOP:
            break
        assert status.tag == TAG.ARGSIZE, "Expected tag ARGSIZE, unexpected tag= "+str(status.tag)+" was received!"
        assert temp == nDim[0] + 1, "Expected a value="+str(nDim[0] + 1)+" from external server, got="+str(temp)

        # calculate the min value
        if (nCons[0] == 0):
            min_result = minimize(fun=call_fObj, x0=x0, method=optMethod, jac=call_fObj_grad,
                                    bounds = xBnd, tol = tolerance, callback = callback, options = options)
        else:
            min_result = minimize(fun=call_fObj, x0=x0, method=optMethod, jac=call_fObj_grad,
                                  bounds = xBnd, constraints={'type': 'eq', 'fun': call_fCons, 'jac': call_fCons_grad},
                                  tol = tolerance, callback = callback, options = options)


        # Send final result to the external server
        pythonComm.Send(min_result.x, dest=externalServerRank, tag=TAG.DONE)
        outputFile.write("\nResult: "+str(min_result.x)+"\n")
        outputFile.flush()

    return
    

try:
    Main()
except RuntimeError as exc:
    message = "\nPYTHON: Runtime error: "+str(exc)+"\nAborting MPI\n\n"
    sys.stderr.write(message)
    outputFile.write(message)
    outputFile.flush()
    MPI.COMM_WORLD.Abort(1)    
except AssertionError as exc:
    message = "\nPYTHON: Assertion error: "+str(exc)+"\nAborting MPI\n\n"
    sys.stderr.write(message)
    outputFile.write(message)
    outputFile.flush()
    MPI.COMM_WORLD.Abort(1)  
#except Exception as exc:
#    message = "\nPYTHON: Unexpected error: "+str(exc)+"\nAborting MPI\n\n"
#    sys.stderr.write(message)
#    outputFile.write(message)
#    outputFile.flush()
#    MPI.COMM_WORLD.Abort(1)    

MPI.Finalize()
outputFile.close()
sys.exit(0)
