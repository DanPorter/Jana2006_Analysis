"""
General Functions

Taken from Dans_Diffraction/functions_general.py

Version 0.1 12/8/19

Dan Porter
August 2019
"""

import numpy as np


# from Dans_Diffraction/functions_general.py
def stfm(val, err):
    """
    Create standard form string from value and uncertainty"
     str = stfm(val,err)
     Examples:
          '35.25 (1)' = stfm(35.25,0.01)
          '110 (5)' = stfm(110.25,5)
          '0.0015300 (5)' = stfm(0.00153,0.0000005)
          '1.56(2)E+6' = stfm(1.5632e6,1.53e4)

    Notes:
     - Errors less than 0.01% of values will be given as 0
     - The maximum length of string is 13 characters
     - Errors greater then 10x the value will cause the value to be rounded to zero
    """

    # Determine the number of significant figures from the error
    if err == 0. or val / float(err) >= 1E5:
        # Zero error - give value to 4 sig. fig.
        out = '{:1.5G}'.format(val)
        if 'E' in out:
            out = '{}(0)E{}'.format(*out.split('E'))
        else:
            out = out + ' (0)'
        return out
    elif np.log10(np.abs(err)) > 0.:
        # Error > 0
        sigfig = np.ceil(np.log10(np.abs(err))) - 1
        dec = 0.
    elif np.isnan(err):
        # nan error
        return '{} (-)'.format(val)
    else:
        # error < 0
        sigfig = np.floor(np.log10(np.abs(err)) + 0.025)
        dec = -sigfig

    # Round value and error to the number of significant figures
    rval = round(val / (10. ** sigfig)) * (10. ** sigfig)
    rerr = round(err / (10. ** sigfig)) * (10. ** sigfig)
    # size of value and error
    pw = np.floor(np.log10(np.abs(rval)))
    pwr = np.floor(np.log10(np.abs(rerr)))

    max_pw = max(pw, pwr)
    ln = max_pw - sigfig  # power difference

    if np.log10(np.abs(err)) < 0:
        rerr = err / (10. ** sigfig)

    # Small numbers - exponential notation
    if max_pw < -3.:
        rval = rval / (10. ** max_pw)
        fmt = '{' + '0:1.{:1.0f}f'.format(ln) + '}({1:1.0f})E{2:1.0f}'
        return fmt.format(rval, rerr, max_pw)

    # Large numbers - exponential notation
    if max_pw >= 4.:
        rval = rval / (10. ** max_pw)
        rerr = rerr / (10. ** sigfig)
        fmt = '{' + '0:1.{:1.0f}f'.format(ln) + '}({1:1.0f})E+{2:1.0f}'
        return fmt.format(rval, rerr, max_pw)

    fmt = '{' + '0:0.{:1.0f}f'.format(dec + 0) + '} ({1:1.0f})'
    return fmt.format(rval, rerr)


# From Dans_Diffraction/functions_crystallography.py
def gen_sym_pos(sym_ops, x, y, z):
    """
    Generate positions from symmetry operations
    Usage:
      uvw = gen_sym_pos(sym_ops,x,y,z)
      sym_ops = [n*'x,y,z'] array of string symmetry operations
      x,y,z = fractional coordinates of atomic posiiton to be modified by symmetry
      uvw = [[nx3]] array of symmetry defined factional coordinates [u,v,w]

    E.G.
      uvw = gen_sym_pos(['x,y,z','y,-x,z+1/2'],0.1,0.2,0.3)
      uvw >> [[0.1,0.2,0.3] , [0.2,-0.1,0.8]]
    """
    uvw = np.zeros([len(sym_ops), 3])
    for n in range(len(sym_ops)):
        sym = sym_ops[n]
        sym = sym.lower()
        # Evaluate string symmetry operation in terms of x,y,z
        sym = sym.replace('/', './')
        sym = sym.strip('\"\'')
        out = eval(sym)
        uvw[n] = np.array(out[0:3]) + 0.0  # add zero to remove -0.0 values
    return uvw