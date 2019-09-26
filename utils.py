import numpy
import scipy.special
import pandas as pd
import subprocess


def memory_expensive_random_choose(l, k, n=5000):
    """
    Randomly sample n k-tuples from l without returning 
    """
    l = list(l)
    n = min(int(scipy.special.comb(len(l), k, repetition=False)), n)
    cvs = set()
    while len(cvs) < n:
        elm = numpy.random.choice(l,k,replace=False)
        elm = tuple(elm)
        elm = frozenset(elm)
        if elm in cvs:
            continue
        cvs.add(elm)
        yield elm


def to_dict_no_nans(df):
    d = dict((k, v.dropna().to_dict()) for k, v in pd.compat.iteritems(df))
    return d


def subdict(bigdict, subset_keys):
    intersection_keys = bigdict.keys() & set(subset_keys)
    return {k: bigdict[k] for k in intersection_keys}


def run_command(cmd, get_output=False):
    "run command"

    if get_output == False:

        process = subprocess.Popen(cmd)

        process.wait()

    else:

        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        # wait for the process to terminate
        stdout, stderr = process.communicate()

        stdout = stdout.decode('UTF-8').strip()
        stderr = stderr.decode('UTF-8').strip()

        return stdout, stderr, process.returncode
