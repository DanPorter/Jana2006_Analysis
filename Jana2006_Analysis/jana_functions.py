"""
Python functions to load and display refinement data from Jana2006

Taken from script: Python/Jana2006_functions.py

Version 0.1 12/8/19

Dan Porter
August 2019
"""

import sys, os, re, time
import numpy as np

from . import general_functions as dgp

cf = os.path.dirname(__file__)
ini_dir = os.path.expanduser('~')

pi = np.pi  # mmmm tasty Pi
e = 1.6021733E-19  # C  electron charge
h = 6.62606868E-34  # Js  Plank consant
c = 299792458  # m/s   Speed of light
u0 = 4 * pi * 1e-7  # H m-1 Magnetic permeability of free space
me = 9.109e-31  # kg Electron rest mass
Na = 6.022e23  # Avagadro's No
A = 1e-10  # m Angstrom
Cu = 8.048  # Cu-Ka energy, keV
Mo = 17.4808  # Mo-Ka energy, keV

"--------------------------------------------------------------------"
"----------------------------File Readers----------------------------"
"--------------------------------------------------------------------"


def readm40(filename=None, debug=False):
    """
    Read Jana .m40 file
    Return atomic positions, ADPs and errors
    return crys, error
    crys, error are dicts with keys:
    'name'
    'atom position'
    'atom type'
    'atom label'
    'occupancy'
    'uaniso'
    """

    # Get file name
    (dirName, filetitle) = os.path.split(filename)
    (fname, Ext) = os.path.splitext(filetitle)

    # Open file
    with open(filename) as file:
        lines = file.readlines()
    if debug: print('%s'%filename)

    # Find first line
    for n in range(len(lines)):
        if len(re.findall('\d+', lines[n])) == 4:
            fst = n
            break
    for n in range(fst, len(lines)):
        if len(re.findall('[a-zA-Z_]+', lines[n])) > 0:
            fst_atom = n
            break

    natoms = int(lines[fst].split()[0])
    scale = float(lines[fst + 1].split()[0])
    extinction = float(lines[fst + 4].split()[0])
    if debug: print('Atoms: %2d %s %s'%(natoms, scale, extinction))

    # Atoms
    atom_name = ['' for x in range(natoms)]
    atom_type = np.zeros(natoms)
    atom_Utype = np.zeros(natoms)
    atom_pos = np.zeros([natoms, 3])
    atom_occ = np.zeros(natoms)
    atom_U = np.zeros([natoms, 6])
    extraline = 0
    for N in range(natoms):
        atomline = fst_atom + 2 * N + extraline
        ln1 = lines[atomline].split()  # Ru1       3  2     0.166667 0.000000 0.000000 0.000000
        ln2 = lines[atomline + 1]  # 0.000970 0.000970 0.003098 0.000485 0.000000 0.000000      0000101000
        vals1 = re.findall('.\d+\.\d+', lines[atomline])
        if debug:
            print('\n%s\n%s' % (lines[atomline].strip(), lines[atomline + 1].strip()))

        atom_name[N] = ln1[0]
        atom_type[N] = int(ln1[1])
        atom_Utype[N] = int(ln1[2])
        atom_pos[N, :] = [float(x) for x in vals1[1:]]
        atom_occ[N] = float(vals1[0])
        atom_U[N, :] = [float(ln2[x:x + 9]) for x in range(0, 6 * 9, 9)]

        if debug: print(atom_name[N], atom_type[N], atom_Utype[N], atom_pos[N, :], atom_occ[N], atom_U[N, :])

        # atoms with anharmonic Uaniso values have 2 extra lines
        if atom_Utype[N] > 2:
            extraline += 2

    "---------------s.u. block---------------"
    sublock = fst_atom + 2 * (natoms) + extraline + 1

    Dscale = float(lines[sublock].split()[0])
    Dextinction = float(lines[sublock + 3].split()[0])
    if debug: print('\ns.u.block: %s, %s' % (Dscale, Dextinction))

    # Atom Errors
    Datom_pos = np.zeros([natoms, 3])
    Datom_occ = np.zeros(natoms)
    Datom_U = np.zeros([natoms, 6])
    extraline = 0
    for N in range(natoms):
        atomline = sublock + (fst_atom - fst - 1) + 2 * N + extraline
        ln1 = lines[atomline].split()  # Ru1       3  2     0.166667 0.000000 0.000000 0.000000
        ln2 = lines[atomline + 1]  # 0.000970 0.000970 0.003098 0.000485 0.000000 0.000000      0000101000

        if debug:
            print('\n%s\n%s' % (lines[atomline].strip(), lines[atomline + 1].strip()))

        Datom_pos[N, :] = [float(x) for x in ln1[-3:]]
        Datom_occ[N] = float(ln1[1])
        Datom_U[N, :] = [float(ln2[x:x + 9]) for x in range(0, 6 * 9, 9)]

        if debug: print(atom_name[N], Datom_pos[N, :], Datom_occ[N], Datom_U[N, :])

        if atom_Utype[N] > 2:
            extraline += 2

    crys = {}
    crys['name'] = fname
    crys['atom position'] = atom_pos
    crys['atom type'] = atom_type
    crys['atom label'] = atom_name
    crys['occupancy'] = atom_occ
    crys['uaniso'] = atom_U

    error = {}
    error['name'] = fname
    error['atom position'] = Datom_pos
    error['atom type'] = atom_type
    error['atom label'] = atom_name
    error['occupancy'] = Datom_occ
    error['uaniso'] = Datom_U
    return crys, error


def readref(filename=None):
    """
    Read Jana .ref file
    Return R factors and GoF
    Return: GOFobs, GOFall, Robs, Rall, Rwobs, Rwall
    """

    # Open file
    with open(filename) as file:
        txt = file.read()

    GOFobs = float(re.findall('\d+.\d+', re.findall('GOF.obs\D+\d+.\d+', txt)[0])[0])
    GOFall = float(re.findall('\d+.\d+', re.findall('GOF.all\D+\d+.\d+', txt)[0])[0])
    Robs = float(re.findall('\d+.\d+', re.findall('R.obs\D+\d+.\d+', txt)[0])[0])
    Rall = float(re.findall('\d+.\d+', re.findall('R.all\D+\d+.\d+', txt)[0])[0])
    Rwobs = float(re.findall('\d+.\d+', re.findall('wR.obs\D+\d+.\d+', txt)[0])[0])
    Rwall = float(re.findall('\d+.\d+', re.findall('wR.all\D+\d+.\d+', txt)[0])[0])

    return GOFobs, GOFall, Robs, Rall, Rwobs, Rwall


def readm50(filename=None):
    """
    Read Jana .m50 file
    Return symmetry operations, space group and lattice parameters
    return lattpar, SG, SGN, sym_cen
    """

    # Get file name
    (dirName, filetitle) = os.path.split(filename)
    (fname, Ext) = os.path.splitext(filetitle)

    # Open file
    with open(filename) as file:
        txt = file.read()

    lat_txt = re.findall('cell [-0-9.]+ [-0-9.]+ [-0-9.]+ [-0-9.]+ [-0-9.]+ [-0-9.]+', txt)
    lat_typ = re.findall('lattice [A-Z]', txt)[0].split()[1]
    grp_txt = re.findall('spgroup \S+ \d+', txt)[0].split()
    sym_txt = re.findall('[\s-]+[xyz][\s-]+[xyz][\s-]+[xyz]', txt)
    sym_txt = [x.strip() for x in sym_txt]
    sym_txt = [x.replace(' ', ',') for x in sym_txt]

    # Cell centering
    sym_cen = []
    cen_fmt = ['{},{},{}']  # P
    if lat_typ == 'A':
        cen_fmt += ['{},{}+1/2,{}+1/2']  # A
    elif lat_typ == 'B':
        cen_fmt += ['{}+1/2,{},{}+1/2']  # B
    elif lat_typ == 'C':
        cen_fmt += ['{}+1/2,{}+1/2,{}']  # C
    elif lat_typ == 'I':
        cen_fmt += ['{}+1/2,{}+1/2,{}+1/2']  # I
    elif lat_typ == 'F':
        cen_fmt += ['{}+1/2,{}+1/2,{}']  # F
        cen_fmt += ['{},{}+1/2,{}+1/2']  # F
        cen_fmt += ['{}+1/2,{},{}+1/2']  # F
    elif lat_typ == 'R':
        cen_fmt += ['{}+2/3,{}+1/3,{}+1/3']  # R
        cen_fmt += ['{}+1/3,{}+2/3,{}+2/3']  # R

    for sym in sym_txt:
        x, y, z = sym.split(',')
        for fmt in cen_fmt:
            sym_cen += [fmt.format(x, y, z)]

    lattpar = [float(x) for x in lat_txt[0].replace('cell', '').split()]
    SG = grp_txt[1]
    SGN = int(grp_txt[2])
    return lattpar, SG, SGN, sym_cen


"--------------------------------------------------------------------"
"------------------------------Refinements---------------------------"
"--------------------------------------------------------------------"


def refine(filename=None, notes=''):
    """
    Display last Jana2006 refinement as a well formated text, and save this to a file
    """

    # Get file name
    (dirName, filetitle) = os.path.split(filename)
    (fname, Ext) = os.path.splitext(filetitle)

    # Read m40
    crys, err = readm40(os.path.join(dirName, fname + '.m40'))
    # Read ref file
    GOFobs, GOFall, Robs, Rall, Rwobs, Rwall = readref(os.path.join(dirName, fname + '.ref'))
    # Refinement time
    refDate = time.ctime(os.path.getmtime(filename))
    # Read m50
    LP, space_group, SGn, sym = readm50(os.path.join(dirName, fname + '.m50'))

    # Fractional Occupancy
    natom = len(crys['atom label'])
    nsymm = float(len(sym))
    occ_frac = np.ones(natom)
    for N in range(natom):
        all_pos = dgp.gen_sym_pos(sym, *crys['atom position'][N])
        occ_frac[N] = nsymm / len(all_pos)

    # Create Note file
    fout = 'Refinement ' + time.strftime('%Y %m%b %d') + '.txt'
    out = open(os.path.join(dirName, fout), 'a')

    out.write('--------------' + refDate + '--------------\n')
    out.write(filename[:-3] + '\n')
    out.write(notes + '\n')
    out.write('GoF = {:5.2f} R = {:5.2f} Rw = {:5.2f}\n'.format(GOFall, Rall, Rwall))
    print('--------------' + refDate + '--------------')
    print(notes)
    print('GoF = {:5.2f} R = {:5.2f} Rw = {:5.2f}'.format(GOFall, Rall, Rwall))

    # Display each atom
    for N in range(natom):
        x = dgp.stfm(crys['atom position'][N, 0], err['atom position'][N, 0])
        y = dgp.stfm(crys['atom position'][N, 1], err['atom position'][N, 1])
        z = dgp.stfm(crys['atom position'][N, 2], err['atom position'][N, 2])
        o = dgp.stfm(crys['occupancy'][N] * occ_frac[N], err['occupancy'][N] * occ_frac[N])
        u11 = dgp.stfm(crys['uaniso'][N, 0], err['uaniso'][N, 0])
        u22 = dgp.stfm(crys['uaniso'][N, 1], err['uaniso'][N, 1])
        u33 = dgp.stfm(crys['uaniso'][N, 2], err['uaniso'][N, 2])
        u12 = dgp.stfm(crys['uaniso'][N, 3], err['uaniso'][N, 3])
        u13 = dgp.stfm(crys['uaniso'][N, 4], err['uaniso'][N, 4])
        u23 = dgp.stfm(crys['uaniso'][N, 5], err['uaniso'][N, 5])
        fmt = '{:8s} x:{:12s} y:{:12s} z:{:12s} occ:{:12s}\n'
        fmt += 'U11:{:12s} U22:{:12s} U33:{:12s} U12:{:12s} U13:{:12s} U23:{:12s}'
        print(fmt.format(crys['atom label'][N], x, y, z, o, u11, u22, u33, u12, u13, u23))
        out.write(fmt.format(crys['atom label'][N], x, y, z, o, u11, u22, u33, u12, u13, u23))

    if np.any(crys['uaniso'] < -3 * err['uaniso']):
        print('***{} Negative ADPs***'.format(np.sum(crys['uaniso'] < -3 * err['uaniso'])))
        out.write('***{} Negative ADPs***\n'.format(np.sum(crys['uaniso'] < -3 * err['uaniso'])))
    if np.any(crys['uaniso'] > 0.1):
        print('***{} Large ADPs***'.format(np.sum(crys['uaniso'] > 0.1)))
        out.write('***{} Large ADPs***\n'.format(np.sum(crys['uaniso'] > 0.1)))
    if np.any(crys['occupancy'] < -err['occupancy']):
        print('***{} Negative Occupancies***'.format(np.sum(crys['occupancy'] < -3 * err['occupancy'])))
        out.write('***{} Negative Occupancies***\n'.format(np.sum(crys['occupancy'] < -3 * err['occupancy'])))
    if np.any(crys['occupancy'] * np.array(occ_frac) > 1):
        print('***{} Occupancies > 1***'.format(np.sum(crys['occupancy'] * np.array(occ_frac) > 1)))
        out.write('***{} Occupancies > 1***\n'.format(np.sum(crys['occupancy'] * np.array(occ_frac) > 1)))
    out.write('\n')
    out.close()


def refinementtable(filename=None, occ_frac=None):
    """
    Generate a latex table of the refinement results
    """

    # Get file name
    (dirName, filetitle) = os.path.split(filename)
    (fname, Ext) = os.path.splitext(filetitle)

    # Read m40
    crys, err = readm40(os.path.join(dirName, fname + '.m40'))
    # Read ref file
    GOFobs, GOFall, Robs, Rall, Rwobs, Rwall = readref(os.path.join(dirName, fname + '.ref'))
    # Refinement time
    refDate = time.ctime(os.path.getmtime(filename))
    # Read m50
    LP, space_group, SGn, sym = readm50(os.path.join(dirName, fname + '.m50'))

    # Fractional Occupancy
    natom = len(crys['atom label'])
    nsymm = float(len(sym))
    occ_frac = np.ones(natom)
    for N in range(natom):
        all_pos = dgp.gen_sym_pos(sym, *crys['atom position'][N])
        occ_frac[N] = nsymm / len(all_pos)

    print('\\begin{table}[htp]')
    print('    \\centering')
    print('       \\begin{tabular}{c|c|ccccc}')
    print('             & Site & x & y & z & Occ. & U$_{iso}$ \\\\ \hline')
    for N in range(natom):
        if err['atom position'][N, 0] > 0:
            x = dgp.stfm(crys['atom position'][N, 0], err['atom position'][N, 0])
        else:
            x = str(crys['atom position'][N, 0])
        if err['atom position'][N, 1] > 0:
            y = dgp.stfm(crys['atom position'][N, 1], err['atom position'][N, 1])
        else:
            y = str(crys['atom position'][N, 1])
        if err['atom position'][N, 2] > 0:
            z = dgp.stfm(crys['atom position'][N, 2], err['atom position'][N, 2])
        else:
            z = str(crys['atom position'][N, 2])
        if err['occupancy'][N] > 0:
            o = dgp.stfm(crys['occupancy'][N] * occ_frac[N], err['occupancy'][N] * occ_frac[N])
        else:
            o = str(crys['occupancy'][N] * occ_frac[N])

        # This only works for cubic crystals!!!
        if crys['uaniso'][N, 1] == 0:
            Uiso = crys['uaniso'][N, 0]
            dU = err['uaniso'][N, 0]
        else:
            Uiso = np.mean(crys['uaniso'][N, :3])
            dU = np.sqrt(np.sum(err['uaniso'][N, :3] ** 2)) / 3
        Uiso = dgp.stfm(Uiso, dU)
        nsite = int(nsymm / occ_frac[N])
        print('        {} & ${}a$ & {} & {} & {} & {} & {} \\\\'.format(crys['atom label'][N], nsite, x, y, z, o, Uiso))
    print('        \\end{tabular}')
    print('    \\caption{{Refinement of sample with R$_w$ = {}\%.}}'.format(Rwall))
    print('    \\label{tab:}')
    print('\\end{table}')

