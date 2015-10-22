
import sys
import subprocess

import pandas
from pandas import DataFrame, Series

nCPU = 2
BLOCKSIZE = 96
NB = 256
MB = 256

def performance_test(name, sizes, mB, nB, nVP, nWRK):
    cmd = ["./mmperf",
           "-T", "%s" % name,
           "-L", "%s" % ",".join(map(lambda x: "%d"%x, sizes)),
           "-W", "%d" % nWRK,
           "-NB", "%d" % nB,
           "-MB", "%d" % mB,
           "-H", "%d" % nVP,
           "-g"] 

    perfdata = subprocess.check_output(cmd).strip()
    return eval(perfdata)
    
name = "Gemm"
filename = "Gemm.dat"
if len(sys.argv[1:]) > 0:
    name = sys.argv[1]
    filename = name + ".dat"
if len(sys.argv[1:]) > 1:
    filename = sys.argv[2]

if len(sys.argv[1:]) > 2:
    (NB, MB, VP, nW) = map(lambda x: int(x), sys.argv[3].split(':'))
else:
    NB = 256
    MB = 256
    VP = 96
    nW = 2

if len(sys.argv[1:]) > 3:
    sizes = map(lambda x: int(x), sys.argv[4:])
else:
    sizes = [10, 20, 40, 60, 80, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300]

print "running with params mB=%d, nB=%d, nVP=%d nWRK=%d sizes: %s" %(MB, NB, VP, nW, sizes)
pdata = performance_test(name, sizes, NB, MB, VP, nW)
series = Series(pdata)
print "GFlops/s: %s" % name
print series

if filename[0] != "-":
    pandas.save(series, filename)
