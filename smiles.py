"""
  It runs prereduction depending on the given reduction
   parameters and options. 

  It utilizes two external libraries: OpenBabel and UManSysProp

"""

import os
import sys
import math

from utils import isfloat
from parameters import path_to_umansysprop

# OpenBabel (to be installed)
from openbabel import pybel, OBConversion, OBMol, OBMolAtomIter, OBAtomBondIter

# UManSysProp (to be downloaded)
if not os.path.exists(path_to_umansysprop):
    raise ValueError(
        "The path to UManSysProp is not found. Please update the path in the parameter file."
    )
sys.path.append(path_to_umansysprop)
from umansysprop import (
    groups,
    boiling_points,
    vapour_pressures,
    activity_coefficient_models,
)

mthVP = ["nannoolal", "myrdal_and_yalkowsky"]
mthBP = ["nannoolal", "stein_and_brown", "joback_and_reid"]

# Basic SMART structures
keyels = {
    "C": "[#6]",  # carbon
    "H": "[H]",  # hydrogen
    "O": "[O]",  # oxygen 6
    "N": "[N]",  # nitrogen 7
    #'S':'[S]'
}
key_pybel = {}
for key, val in keyels.items():
    key_pybel[key] = pybel.Smarts(val)


def get_pybelmol(mol):
    """get a molecule to pybel format"""
    if isinstance(mol, pybel.Molecule):
        return mol
    # if input a SMILES string
    elif isinstance(mol, str):
        try:
            mol = pybel.readstring("smi", mol)
            return mol
        except IOError:
            raise ValueError("get_pybelmol: not a smi string: ", mol)
    else:
        print(mol)
        raise ValueError(
            "get_pybelmol: not recognize data type, not mol and not string.", mol
        )


def get_functional_groups(smiles):
    """return a dictionary of the key functional groups read by pybel
    the key functional group can be modified in Parameter.py"""

    pymol = get_pybelmol(smiles)

    functionalGroups = {}
    # for key,val in smarts_pybel.items():
    for key, val in key_pybel.items():
        tmp = val.findall(pymol)
        if tmp != []:
            functionalGroups[key] = len(tmp) * 1.0
        else:
            functionalGroups[key] = 0.0
    return functionalGroups


def get_organic_ratios_and_du(smiles):
    """Get the number of key chemical elements in the given molecule,
    return OM/OC mass ratio, H/C, O/C, N/C atomic ratios and the degree of unsaturation
    """
    # For a compound with formula CaHbNcOdXe where X is F, Cl, Br or I, the degree of unsaturation is given by:
    # degree of unsaturation = 1/2 (2 + 2C + N - H - X)
    pymol = get_pybelmol(smiles)
    tmp = {}
    for key, val in key_pybel.items():
        tmp[key] = float(len(val.findall(pymol)))
    return {
        "OM/OC": round(pymol.molwt / (tmp["C"] * 12.0), 3),
        "H/C": round(tmp["H"] * 1.0 / tmp["C"], 3),
        "O/C": round(tmp["O"] * 1.0 / tmp["C"], 3),
        "N/C": round(tmp["N"] * 1.0 / tmp["C"], 3),
    }, int(0.5 * (2 + 2 * tmp["C"] - tmp["H"] + tmp["N"]))


# Not finished !!!
def get_toSOAPStructure(mol):
    """get SSH-aerosol strucutre of input compounds.
    input: one SMILES strucutre -> pymol
    output: a dict of SSH-aerosol functional groups {index: number,...}"""
    m = groups.to_SSH(mol)
    # transfer to SOAP format
    mkeys = list(sorted(m.keys()))
    return [mkeys, [float(m[i]) for i in mkeys]]


def get_saturation_vapor_pressure(mol, Type="evap", temperature=298):
    """return (Psat,T) Psat at the certain temperature T in unit: atm, K"""

    pymol = get_pybelmol(mol)

    # Boiling points [(K)] bp
    # Vapour pressures [log10 (atm) at a specific temperature] vp

    if Type == "evap":
        vp = vapour_pressures.evaporation(pymol, temperature)
    elif Type == "evap2":
        vp = vapour_pressures.evaporation2(pymol, temperature)
    elif Type == "simpol":
        vp = vapour_pressures.simpol(pymol, temperature)
    elif "VP" in Type and "BP" in Type:
        # example VP0BP0
        bp = eval("boiling_points." + mthBP[int(Type[5])] + "(pymol)")
        vp = eval(
            "vapour_pressures." + mthVP[int(Type[2])] + "(pymol, temperature, bp)"
        )
    else:
        vp = []
        # obtain vp
        for i in ["evaporation", "evaporation2", "simpol"]:
            # print mol, i, eval('vapour_pressures.'+i+'(pymol, temperature)')
            vp.append(eval("vapour_pressures." + i + "(pymol, temperature)"))

        for i in mthBP:
            bp = eval("boiling_points." + i + "(pymol)")
            for j in mthVP:
                # print(i,j,eval('vapour_pressures.'+j+'(pymol, temperature, bp)'))
                vp.append(eval("vapour_pressures." + j + "(pymol, temperature, bp)"))

        # compute
        if Type == "ave":
            num = 0
            for i in vp:
                num += 10**i
            return num / (len(vp) * 1.0)

        elif Type == "ave_log10":
            return 10 ** (sum(vp) / (len(vp) * 1.0))

        elif Type == "ave_log10_3":
            vp.sort()
            return 10 ** (sum(vp[0:3]) / 3.0)

        elif Type == "ave_log10_3_2":
            vp.sort()
            return 10 ** (sum(vp[0:3]) / 3.0) / 2.0

        elif Type == "min" or Type == "max":
            return 10 ** eval(Type + "(vp)")

        else:
            sys.exit("MP: not recognize saturation vapor pressure type")

    return 10**vp


def get_enthalpy_vaporization(mol, Type="evap", temp1=298.0, temp2=308.0):
    """return Hvap in unit J.mol-1"""
    if temp1 != temp2:
        psat1 = get_saturation_vapor_pressure(mol, temperature=temp1, Type=Type)
        psat2 = get_saturation_vapor_pressure(mol, temperature=temp2, Type=Type)

        #  Clausius-Clapeyron equation
        # print(psat1,temp1,psat2,temp2,1/temp1-1/temp2)
        return abs(math.log(psat2 / psat1) * 8.31446 / (1.0 / temp1 - 1.0 / temp2))

    else:
        sys.exit("Hvap: temp1 = temp2, " + str(temp1))


# for gamma
water = get_pybelmol("O")


def get_activity_coefficient_inf(compound, temperature=298.0):
    """This routine compute the activity coefficients at infinite dilution GAMMAinf
    given by UNIFAC for each species and compute the Henry's law constant from
    the saturation vapour pressure
    H = 1000*760/(18.0*GAMMAinf*Psat)"""

    item = get_pybelmol(compound)
    organic = {item: 1, water: 1e50}  # 1e10

    dic, m = activity_coefficient_models.calculate_activities_org(organic, temperature)

    dic = {}
    for c in m:
        dic[c.compound] = c.activity_coefficient

    # check water
    if dic[water] != 1.0:
        sys.exit(
            "water is not infinite: {:f}, gamma_inf: {:f}".format(dic[water], dic[item])
        )
    else:
        return dic[item]


def get_activity_coefficient_org(compounds, concs, temperature=298.0):
    """input: species name and pymol in a list"""
    organics = {}
    # add conc.
    for i in range(len(compounds)):
        organics[compounds[i].pymol] = concs[i]
    # compute
    dic, m = activity_coefficient_models.calculate_activities_org(organics, temperature)
    # out info
    out = []
    for i in compounds:
        tag = 0
        for c in m:
            if i.pymol == c.compound:
                out.append(c.activity_coefficient)
                # print(i.name,c.activity_coefficient)
                tag = 1
                break
        if not tag:
            sys.exit("activity_coefficient_org: not find gamma for species " + i.name)
    return out


def get_partitioning_coefficient(compund, temperature=298.0, Mavg=200):
    Kp = 8.314 * temperature / (compund.gamma * compund.psat_atm * Mavg)
    # print 'Kp',Kp,'psat_atm',compund.psat_atm
    return Kp
    # return 8.314*temperature/(self.gamma*self.psat_atm*Mavg)


def get_dHvap_simpol(SOAPStructure):
    simpol1 = {
        "bo": -9.0677e02,
        "C": -2.3229e02,
        "C=C": 5.9336e01,
        "OH": -8.7901e02,
        "aldehyde": -5.2267e02,
        "ketone": 1.9917e01,
        "COOH": -1.1963e03,
        "nitrate": -7.8253e02,
        "peroxide": 4.4567e02,
        "hydroperoxide": -7.9762e02,
        "aromatic ring": -1.3635e02,
        "ether": -2.2814e02,
        "phenol": -4.2981e02,
        "nitrophenol": 2.8685e02,
    }
    simpol0 = {
        "C": [0, 1, 2, 3],
        "C=C": [16, 17, 18, 19, 20],
        "OH": [26],
        "aldehyde": [31],
        "ketone": [29, 30],
        "COOH": [37],
        "nitrate": [39, 40, 41],
        "peroxide": [45, 46, 47, 48, 49, 50, 51, 52, 53],
        "hydroperoxide": [42, 43, 44],
        "aromatic ring": [21, 22],
        "ether": [34, 35, 36],
        "phenol": [28],
        "nitrophenol": [38],
    }
    ksim = list(simpol0.keys())
    dH = 0.0
    for i in SOAPStructure:
        tag = 0
        for k in ksim:
            if i in simpol0[k]:
                dH += SOAPStructure[i] * simpol1[k]
                tag = 1
                break
        if tag == 0:
            print("!!!!not find strucutre: ", i)
    if dH != 0.0:
        return -1 * (dH + simpol1["bo"]) * (2.303 * 8.314) / 1000
    else:
        return 0.0


def string_to_smiles(string, to_canonical=False):
    """Convert input SMILES to its canonical SMILES format"""

    # initialize Open Babel objects
    obConversion = OBConversion()
    obConversion.SetInAndOutFormats("smi", "can")

    mol = OBMol()
    is_valid = obConversion.ReadString(mol, string)

    # check if input SMILES is valid
    if not is_valid:
        print(f"Error: Invalid SMILES string: {string}")
        return False

    if to_canonical:  # convert to canonical SMILES
        smiles = obConversion.WriteString(mol).strip()
    else:
        smiles = string

    return smiles


def gecko_to_smiles(gecko, tag_canonical=True):
    """Convert input GECKO-A format to its canonical (optional) SMILES structure"""

    # convert functional groups (this might need to be expanded for all functional groups)
    functional_groups = {
        # C
        "Cd": "C",
        "c": "c",
        # H
        "H4": "",
        "H3": "",
        "H2": "",
        "H": "",
        # C=O
        "CHO": "C(=O)",
        "CO": "C(=O)",  #'CO' may overlap with COO!!!
        "C1O": "C1(=O)",
        "C2O": "C2(=O)",
        # NO2
        # "NO2": "[N+](=O)[O-]",
        "NO2": "N(=O)=O",
        # O
        "OH": "O",
        # Crieege radicals
        "EOO.": "OO.",
        "ZOO.": "OO.",
        # handwritten
        "#": "",  # need to test !!!
        "#mm": "",
        "-": "",
    }

    gecko_smiles = {}  # specific smiles

    # check if gecko in certain format
    if gecko in gecko_smiles.keys():
        return gecko_smiles[gecko]

    # rank gecko groups from the longest to shortest
    gecko_group = sorted(functional_groups.keys(), key=len, reverse=True)

    # initial check
    if "#" in gecko:
        print("Found handwritten compound: ", gecko)

    # replace gecko groups step by step
    for i in gecko_group:
        if i in gecko:
            gecko = gecko.replace(i, functional_groups[i])

    # check string and convert to canonical SMILES if need
    return string_to_smiles(gecko, tag_canonical)


def gecko_check_group(fgrp):
    """read the functional group code in GECKO and output as a dict"""
    chemical_groups = {
        "A": "acid",
        "D": "aldehyde",
        "E": "ether",
        "EE": "peroxide",
        "G": "peroxy-acid",
        "H": "hydroperoxy",
        "K": "ketone",
        "N": "nitrate",
        "O": "hydroxy",
        "P": "PAN",
        "R": "aromatic ring",
        "T": "aliphatic ring",
        "U": "unsaturated (double bond)",
        "V": "nitrite",
        "X": "ketene",
        "1": "alkoxy",
        "2": "peroxy",
        "3": "acyl RO2",
        "4": "Criegge",
    }

    return {}  # !!! not finished


def mechgen_to_smiles(mechgen, tag_canonical=True):
    """Convert input mechgen format to its canonical (optional) SMILES structure"""

    # convert functional groups (this might need to be expanded for all functional groups)
    functional_groups = {
        # C
        "aC": "c",  # aromatic
        "pC": "C",  # allylic
        # H
        "H4": "",
        "H3": "",
        "H2": "",
        "H": "",
        # Single bonds
        "-": "",
        "-^": "/",  # Cis and trans isomerizatio
        "-v": "\\",
        # C=O
        "CHO": "C(=O)",
        "CO": "C(=O)",  #'CO' may overlap with COO!!!
        "C1O": "C1(=O)",
        "C2O": "C2(=O)",
        # NO2
        "NO2": "N(=O)=O",
        # O
        "OH": "O",
        # radicals
        "C[.]": "[C]",
        "CH[.]": "[CH]",
        "CO[.]": "C[O]",
        "[O.]": "[O]",
        "[OO.]": "O[O]",
        "[OO]": "O[O]",
        # Syn/Anti isomerization not used
    }

    mechgen_smiles = {}  # specific smiles

    # check if contain non-smiles structure
    for i in ["syn", "anti"]:
        if i in mechgen:
            raise NameError("Contain non-smiles structure: ", mechgen)

    # check if mechgen in certain format
    if mechgen in mechgen_smiles.keys():
        return mechgen_smiles[mechgen]

    # rank mechgen groups from the longest to shortest
    mechgen_group = sorted(functional_groups.keys(), key=len, reverse=True)

    # replace mechgen groups step by step
    for i in mechgen_group:
        if i in mechgen:
            mechgen = mechgen.replace(i, functional_groups[i])

    # check if monocyclic "*" followed by a number # Must be after C !!!
    stars = [
        i
        for i in range(len(mechgen) - 1)
        if (mechgen[i] == "*" and not isfloat(mechgen[i + 1]))
    ]

    if stars != []:
        mech1 = ""  # new string
        # add number 1 to the closed 1
        for i, s in enumerate(stars):
            # get subline [s0:s1]
            if i == 0:
                subline = mechgen[0:s]
            elif i == len(stars) - 1:
                subline = mechgen[s : len(stars)]
            else:
                subline = mechgen[stars[i - 1] : s]

            # found C index in subline
            ind = -1  # initial index
            for j in range(len(subline) - 1, -1, -1):  # reverse
                if subline[j] == "C":
                    ind = j + 1  # found index for adding number
                    break
            # check if find C
            if ind == -1:
                raise ValueError("Not found C in the subline: ", subline, i, s, mechgen)
            # add number 1 at location ind
            mech1 += f"{subline[0:ind]}1{subline[ind:]}"

        # update string
        mechgen = mech1

    # remove * if contains any
    mechgen = mechgen.replace("*", "")

    # check string and convert to canonical SMILES if need
    return string_to_smiles(mechgen, tag_canonical)


def smiles_to_gecko(smiles):
    """Not finished!!!!"""

    # check input SMILES
    is_valid = string_to_smiles(smiles, False)  # no convert

    if not is_valid:
        print(f"In put invalid SMILES string: {smiles}")
        return False

    # start conversion
    gecko = ""  # init
    ind = 0

    # check atom one by one
    for atom in OBMolAtomIter(mol):
        iout = atom.GetType()
        atomic = atom.GetAtomicNum()
        bonds = OBAtomBondIter(atom)

        # check No.hydrogen
        n = atom.GetImplicitHCount()  # No. hydrogen
        if n > 0:
            iout = +f"H{n}"

        # iterate over its bonds/neighbors
        for bond in bonds:
            # Get the neighboring atom
            nbr_atom = bond.GetNbrAtom(atom)

            # C=C
            if bond.GetBondOrder() == 2:
                atom_type = atom_type.replace("C", "Cd")

        if iout == "":
            iout = atom.GetType()

        gecko += iout
        ind += 1

    return gecko


if __name__ == "__main__":

    if 0:
        soap = {6: 2.0, 29: 1.0, 30: 1.0, 31: 1.0, 43: 3.0}
        print(get_dHvap_simpol(soap))
    if 0:  # test Psat calculation
        mol = "C(=O)C(OO)C(OO)C(OO)C(=O)C(=O)C"
        for mode in [
            "VP0BP0",
            "VP0BP1",
            "VP0BP2",
            "VP1BP0",
            "VP1BP1",
            "VP1BP2",
            "simpol",
            "evap",
            "evap2",
        ]:
            psat = get_saturation_vapor_pressure(mol, Type=mode)
            print(f"Mol: {mol}, Type: {mode}, Psat_atm at 298K: {psat:6.3E}")

    if 0:
        mol = [
            "O=CCCC(=C)C(=O)CC(C)(C)O[N+](=O)[O-]",  #
            "O=CCC(C)(C)C(O)CCC(=O)C",
            "OCCC(=C)C(=O)CC(C)(C)C(O)CC=O",
            "CC(=O)CCC1C(CC1(C)C)C(=O)OO[N+](=O)[O-]",  # C1011PAN
            "OOc1c(C)cccc1O",
        ]
        for i in mol:
            print(i, toSOAPStructure(get_pybelmol(i)))

    if 0:  # test GECKO to SMILES
        mol = ["CO(OONO2)CH2CH(OH)CH2CH2CH2(OO.)", "CO(OONO2)CH2CH(OH)CH2CH2CHO"]
        for i in mol:
            print(f"From GECKO format: {i}\n To SMILES format: {gecko_to_smiles(i)}\n")

    if 0:  # test MechGen to SMILES
        mol = [
            "CH2=C*1-CH*2-CH2-CH2-CH(CH2*2)-C*1(CH3)-CH3",
            "CH3-C(CH3)(CO-CH2-OH)-C(ONO2)(CH2-CH2-O-OH)-CH2-CO-O-OH",
        ]
        for i in mol:
            print(
                f"From MechGen format: {i}\n To SMILES format: {mechgen_to_smiles(i)}\n"
            )

    if 0:  # test SMILES to GECKO
        mol = "CCCO"
        print(f"{smiles_to_gecko(mol)}")
