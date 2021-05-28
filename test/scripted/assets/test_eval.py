import sys
import h5py
import numpy as np

len(sys.argv) < 2 and sys.exit("Usage: test_eval.py FILE|OBJECT file1 [object1] file2 [object2]")
if sys.argv[1] == "FILE":
    len(sys.argv) == 4 or sys.exit("Usage: test_eval.py FILE file1 file2")
elif sys.argv[1] == "OBJECT":
    len(sys.argv) == 6 or sys.exit("Usage: test_eval.py OBJECT file1 object1 file2 object2")
else:
    sys.exit("Usage: test_eval.py FILE|OBJECT file1 [object1] file2 [object2]")

#define recursive comparison of HDF5 elements
def compare(obj1, obj2, tolerance = 1e-5):
    if isinstance(obj1, h5py.Group):
        #compare attributes
        len(obj1.attrs) == len(obj2.attrs) or sys.exit("Number of attributes differs in object %s" % obj1.name)
        if len(obj1.attrs) > 0:
            for k in obj1.attrs:
                eps = np.max(np.abs(obj1.attrs[k]-obj2.attrs[k]))
                eps < tolerance or sys.exit("Deviation found in attribute %s/%s (deviation: %.14f)" % (obj1.name,k,eps))
        #compare members
        len(obj1) == len(obj2) or sys.exit("Number of datasets differs in object %s" % obj1.name)
        for k in obj1.keys():
            compare(obj1[k], obj2[k])
    elif isinstance(obj1, h5py.Dataset):
        #compare dataset
        eps = np.max(np.abs(obj1[:]-obj2[:]))
        eps < tolerance or sys.exit("Deviation found in dataset %s (deviation: %.14f)" % (obj1.name,eps))

#run comparisons
if sys.argv[1] == "FILE":
    with h5py.File(sys.argv[2],"r") as f1, h5py.File(sys.argv[3],"r") as f2:
        compare(f1, f2)
elif sys.argv[1] == "OBJECT":
    with h5py.File(sys.argv[2],"r") as f1, h5py.File(sys.argv[4],"r") as f2:
        d1 = f1[sys.argv[3]]
        d2 = f2[sys.argv[5]]
        compare(d1,d2)

#success
sys.exit(0)