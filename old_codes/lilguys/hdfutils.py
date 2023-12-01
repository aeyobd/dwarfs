import h5py
import numpy as np
import re
from os import path

def make_default_header(N, mass):
    Npart = 2

    header = {}
    header["NumPart_ThisFile"] = (0, N)
    header["NumPart_Total"] = (0, N)
    header["MassTable"] = (0., mass)
    header["Time"] = 0.
    header["Redshift"] = 0.
    header["BoxSize"] = 350.
    header["NumFilesPerSnapshot"] = 1
    return header


def get_h5_header(h5_f):
    h5_header = h5_f["Header"].attrs
    header = {}
    for key, val in h5_header.items():
        if is_bytes(val):
            val = val.decode()
        header[key] = val

    return header


def get_h5_vector(h5_f, key):
    return np.array(h5_f["PartType1/" + key][()])

def set_h5_vector(h5_f, key, val):
    if key in h5_f["PartType1"].keys():
        del h5_f["PartType1/" + key]
    dtype, shape = get_dtype_shape(val)
    h5_f.create_dataset("PartType1/" + key, data=val, dtype=dtype)


def set_h5_header(h5_f, header):
    for key, val in header.items():
        set_h5_header_attr(h5_f, key, val)


def set_h5_header_attr(h5_f, key, val):
    attrs = h5_f["Header"].attrs

    if key in attrs.keys():
        attrs.modify(key, val)
    else:
        dtype, shape = get_dtype_shape(val)
        if is_bytes(val):
            val = val.decode()
        attrs.create(key, val, shape, dtype=dtype)



def get_dtype_shape(val):
    if is_array(val):
        return get_arr_dtype_shape(val)
    else:
        shape = 1
        dtype = get_dtype(val)
        return dtype, shape


def get_arr_dtype_shape(arr):
    if isinstance(arr, np.ndarray):
        shape = arr.shape
        if len(shape) == 1:
            shape = shape[0]
        el = arr.flatten()[0]
        dtype = get_dtype(el)
        return dtype, shape

    elif isinstance(arr, (tuple, list)):
        shape = len(arr)
        el = arr[0]
        dtype = get_dtype(el)
        return dtype, shape
    elif isinstance(arr, (str, np.str_)):
        str_type = h5py.string_dtype(encoding='utf-8')
        return str_type, None
    elif is_bytes(arr):
        str_type = h5py.string_dtype(encoding='ascii')
        return str_type, None
    else:
        raise NotImplementedError("not implemented", type(arr), arr)


def get_dtype(val):
    if isinstance(val, (float, np.floating)):
        return "double"
    if isinstance(val, (int, np.integer)):
        return "int"
    else:
        raise NotImplementedError(type(val))



def is_array(val):
    return hasattr(val, "__len__")

def is_bytes(val):
    return isinstance(val, (bytes, np.bytes_))





def get_epsilon(directory):
    filename = path.join(directory + "/parameters-usedvalues")
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("SofteningComovingClass0"):
                match = re.findall(r"\d*\.?\d+", line)[-1]
                return float(match)
