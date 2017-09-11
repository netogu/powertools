import numpy as np


def find_decade(res_val):
    if res_val < 0:
        return 'error'

    db = np.log(res_val) / np.log(10)

    if db < 1:
        dec = 1
    elif db < 2:
        dec = 10
    elif db < 3:
        dec = 100
    elif db < 4:
        dec = 1000
    elif db < 5:
        dec = 10000
    elif db < 6:
        dec = 100000
    elif db < 7:
        dec = 1000000

    return dec, db


def find_closest(val, array):
    array_min = np.abs(array - val)
    min_index = array_min.argmin()
    return min_index


def e96_std_value(res_val):
    ''' Returns closes 1% comercial value and lower and upper closest alternatives'''

    decade, db = find_decade(res_val)
    index = np.arange(0, 97)
    vals = np.round(10 ** (index / 96.0), 2) * decade
    min_index = find_closest(res_val, vals)
    comm_val = vals[min_index]

    print("\n--------- . . . . --------------- . . . . -------------")
    print("\nClosest 1% Resistance for {:.1f} is {:.1f} with {:.1f}% of error".format(res_val, comm_val, np.abs(
        (comm_val - res_val) / res_val * 100)))
    print("\n--------- . . . . --------------- . . . . -------------")

    return comm_val, vals[min_index - 1], vals[min_index + 1]
