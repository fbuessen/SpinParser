#set up pythonpath and import modules
import sys
pythonpath = sys.argv[1]
sys.path.append(pythonpath)

import numpy as np
import spinparser.obs as o

#read filename to use for testing
file = sys.argv[2]

#run tests
##test getLatticeBasis
reference = [[0.0,0.0,0.0],[1.0,0.0,0.0]]
data = o.getLatticeBasis(file)
np.isclose(data, reference).all() or sys.exit("Test getLatticeBasis failed.")
print("Test getLatticeBasis passed.")

##test getLatticePrimitives
reference = [[1.5,0.8660254,0.0],[1.5,-0.8660254,0.0],[0.0,0.0,1.0]]
data = o.getLatticePrimitives(file)
np.isclose(data, reference).all() or sys.exit("Test getLatticePrimitives failed.")
print("Test getLatticePrimitives passed.")

##test getLatticeSites
reference = [[0.0, 0.0, 0.0], [-0.5, 0.8660253882408142, 0.0], [0.0, 1.7320507764816284, 0.0], [-0.5, 2.598076105117798, 0.0], [-2.0, 0.0, 0.0], [1.5, 0.8660253882408142, 0.0], [1.5, -0.8660253882408142, 0.0], [-1.5, -0.8660253882408142, 0.0], [0.0, -1.7320507764816284, 0.0], [-1.5, 0.8660253882408142, 0.0], [-0.5, -0.8660253882408142, 0.0], [1.0, 0.0, 0.0], [2.5, 0.8660253882408142, 0.0], [1.0, 1.7320507764816284, 0.0], [2.5, -0.8660253882408142, 0.0], [1.0, -1.7320507764816284, 0.0], [-2.0, -1.7320507764816284, 0.0], [-0.5, -2.598076105117798, 0.0], [-2.0, 1.7320507764816284, 0.0]]
data = o.getLatticeSites(file, [0.0,0.0,0.0], verbose=False)
np.isclose(data, reference).all() or sys.exit("Test getLatticeSites failed.")

reference = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.5, 0.8660253882408142, 0.0], [1.5, -0.8660253882408142, 0.0], [-0.5, -0.8660253882408142, 0.0], [-0.5, 0.8660253882408142, 0.0], [2.5, 0.8660253882408142, 0.0], [1.0, 1.7320507764816284, 0.0], [2.5, -0.8660253882408142, 0.0], [1.0, -1.7320507764816284, 0.0], [-1.5, -0.8660253882408142, 0.0], [0.0, -1.7320507764816284, 0.0], [-1.5, 0.8660253882408142, 0.0], [0.0, 1.7320507764816284, 0.0], [3.0, 1.7320507764816284, 0.0], [3.0, 0.0, 0.0], [1.5, 2.598076105117798, 0.0], [3.0, -1.7320507764816284, 0.0], [1.5, -2.598076105117798, 0.0]]
data = o.getLatticeSites(file, [1.0,0.0,0.0], verbose=False)
np.isclose(data, reference).all() or sys.exit("Test getLatticeSites failed.")
print("Test getLatticeSites passed.")

##test getCorrelation
reference = [0.0, -9.50785e-06, -4.39011e-05, -0.000116698, -0.000250599, -0.000482914, -0.000874231, -0.0015225, -0.00257785, -0.00429219, -0.00706185, -0.0115148, -0.0156193, -0.018661, -0.0225649, -0.0277352, -0.0347753, -0.0445845, -0.0585239, -0.0786823, -0.107199, -0.148879, -0.203808, -0.235261, -0.273202, -0.318908, -0.37388, -0.440613, -0.522831, -0.624065, -0.731935, -0.828863, -0.911171, -0.974795]
data = o.getCorrelation(file, cutoff="all", site=[1.0, 0.0, 0.0], reference=[0.0, 0.0, 0.0], component='all', verbose=False).flatten()
np.isclose(data, reference, atol=1e-16).all() or sys.exit("Test getCorrelation failed.")

reference = [0.0, -3.16928e-06, -1.46337e-05, -3.88993e-05, -8.35329e-05, -0.000160971, -0.00029141, -0.000507499, -0.000859284, -0.00143073, -0.00235395, -0.00383827, -0.00520645, -0.00622034, -0.00752163, -0.00924508, -0.0115918, -0.0148615, -0.019508, -0.0262274, -0.035733, -0.0496263, -0.0679361, -0.0784205, -0.0910672, -0.106303, -0.124627, -0.146871, -0.174277, -0.208022, -0.243978, -0.276288, -0.303724, -0.324932]
data = o.getCorrelation(file, cutoff="all", site=[1.0, 0.0, 0.0], reference=[0.0, 0.0, 0.0], component='XX', verbose=False).flatten()
np.isclose(data, reference, atol=1e-16).all() or sys.exit("Test getCorrelation failed.")

reference = [0.686557, -0.324932, 0.0877276, -0.0207016, -0.0392286, 0.0877276, 0.0877276, 0.0877276, 0.0877276, 0.0877276, -0.324932, -0.324932, -0.0207016, -0.0392286, -0.0207016, -0.0392286, -0.0207016, -0.0207016, -0.0207016]
data = o.getCorrelation(file, cutoff=0.309031, site="all", reference=[0.0, 0.0, 0.0], component='ZZ', verbose=False).flatten()
np.isclose(data, reference, atol=1e-16).all() or sys.exit("Test getCorrelation failed.")

reference = [-0.276288, 0.0681969, -0.00924508, 0.000463789]
data = o.getCorrelation(file, cutoff=[0.381520, 2.058910], site=[[0.0, 0.0, 0.0],[2.5, 0.866025, 0.0]], reference=[1.0, 0.0, 0.0], component='ZZ', verbose=False).flatten()
np.isclose(data, reference, atol=1e-16).all() or sys.exit("Test getCorrelation failed.")
print("Test getCorrelation passed.")

##test getStructureFactor
reference = [0.0266215, 4.56356]
data = o.getStructureFactor(file, [[0.0,0.0,0.0],[4.188790,0.0,0.0]], cutoff=0.381520, verbose=False).flatten()
np.isclose(data, reference, atol=1e-16).all() or sys.exit("Test getStructureFactor failed.")
print("Test getStructureFactor passed.")

#success
sys.exit(0)