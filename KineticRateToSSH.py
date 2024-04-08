"""
  KineticMCMtoSSH.py converts the knietic rate contants from MCM format to a format can be executed by SSH-aerosol. 
  This is a part of the GENOA algorithm.
"""

from Convert import isfloat

# GECKO
gecko_to_ssh_keywords = {
 'HV':      'PHOT ',
 'TBODY':   'TB M ',
 'PERO':    'RO2 ', # MEPERO? 01-09
 'EXTRA':   'EXTRA ',
 'OXYGEN':  'TB O2 ',
 'FALLOFF': 'FALLOFF ',
 'ISOM':    'EXTRA 200 '
 }

# MCM
mcm_simple_rate_coefs = {'KRO2NO':1, 'KRO2HO2':2, 'KAPHO2':3, 'KAPNO':4,
                         'KRO2NO3':5, 'KNO3AL':6, 'KDEC':7, 'KROPRIM':8,
                         'KROSEC':9, 'KCH3O2':10, 'K298CH3O2':11, 'K14ISOM1':12}
                         
mcm_complex_rate_coefs = {f'KMT{i:02}':i for i in range(1,19)} # KMT01 to KMT18
mcm_complex_rate_coefs.update({'KFPAN': 21, 'KBPAN': 22, 'KBPPN': 23})

mcm_manual_lists={
        # ssh_genoa_spec
        '3.8E-13*EXP(780/TEMP)*(1-1/(1+498*EXP(-1160/TEMP)))':'EXTRA 99 1', #CH3O2 + HO2 -> CH3OOH
        '3.8E-13*EXP(780/TEMP)*(1/(1+498*EXP(-1160/TEMP)))': 'EXTRA 99 2', #CH3O2 + HO2 -> HCHO
        #CH3O2 + RO2 -> CH3OH + HCHO
        #1.03E-13*math.exp(365/TEMP)*(1-7.18*EXP(-885/TEMP))
        '2*KCH3O2*RO2*(1-7.18*EXP(-885/TEMP))':'RO2 1 EXTRA 99 3',
        # In case that lump CH3COCH3
        '8.8E-12*EXP(-1320/TEMP)+1.7E-14*EXP(423/TEMP)':'EXTRA 99 4',
        #HCOCO -> 1.0000 CO + 1.0000 OH
        #5.00E-12*O2*3.2*(1-EXP(-550/TEMP))
        '5.00E-12*O2*3.2*(1-EXP(-550/TEMP))':'TB O2 EXTRA 99 5',
        #CH3O2 + RO2 -> CH3OH + HCHO
        #1.03E-13*math.exp(365/TEMP)*(1-7.18*EXP(-885/TEMP))
        '2*KCH3O2*RO2*0.5*(1-7.18*EXP(-885/TEMP))':'RO2 1 EXTRA 99 6',
        # Add for LIMONENE species: INDO
        '2.20E+10*EXP(-8174/TEMP)*EXP(1.00E+8/TEMP@(3))':'EXTRA 99 7',
        '8.14E+9*EXP(-8591/TEMP)*EXP(1.00E+8/TEMP@(3))':'EXTRA 99 8',
        
        # Non-EXTRA ratio
        '2*(K298CH3O2*8.0E-12)@(0.5)*RO2':'RO2 1 ARR 3.34664e-12 0 0',
        '2*(K298CH3O2*2.9E-12*EXP(500/TEMP))@(0.5)*RO2':'RO2 1 ARR 2.01494e-12 0 -250',
        '2*(KCH3O2*7.8E-14*EXP(1000/TEMP))@(0.5)*RO2':'RO2 1 ARR 1.7926e-13 0 -682.5',
        '(KCH3O2*7.8E-14*EXP(1000/TEMP))@(0.5)' :'ARR 8.963e-14 0 -682.5',
        '2*KCH3O2*RO2*7.18*EXP(-885/TEMP)':'RO2 1 ARR 1.4791e-12 0 520',
         # TOL MCM
        '2*(KCH3O2*2.4E-14*EXP(1620/TEMP))@(0.5)*RO2':'RO2 1 ARR 9.9438e-14 0 -992.5',
        # NC12H26
        # 2*(1.03E-13*6.4E-14*math.exp(365/TEMP))**0.5 => 1.6238e-13 182.5
        '2*(KCH3O2*6.4E-14)@(0.5)*RO2':'RO2 1 ARR 1.6238e-13 0 -182.5',
        # 2*(3.5E-13*3E-13)**0.5   
        '2*(K298CH3O2*3E-13)@(0.5)*RO2':'RO2 1 ARR 6.4807e-13 0 0',
        # 2 * (1.03E-13*1.6E-12)**0.5 1.6238e-12 (365-2200)*0.5 = -917.5
        '2*(KCH3O2*1.6E-12*EXP(-2200/TEMP))@(0.5)*RO2':'RO2 1 ARR 1.6238e-12 0 917.5',
        # Add for LIMONENE species: INDO
        '1.80E+13*(TEMP/298)@(1.7)*EXP(-4079/TEMP)':'ARR 1.119708E+09 1.7 4079.',
        '1.80E+13*(TEMP/298)@(1.7)*EXP(-4733/TEMP)':'ARR 1.119708E+09 1.7 4733.',
        }

# Need to keep update with MCM website
photolysis_coefs = {1: [6.073E-05, 1.743, 0.474],
                    2: [4.775E-04, 0.298, 0.080],
                    3: [1.041E-05, 0.723, 0.279],
                    4: [1.165E-02, 0.244, 0.267],
                    5: [2.485E-02, 0.168, 0.108],
                    6: [1.747E-01, 0.155, 0.125],
                    7: [2.644E-03, 0.261, 0.288],
                    8: [9.312E-07, 1.230, 0.307],
                    11: [4.642E-05, 0.762, 0.353],
                    12: [6.853E-05, 0.477, 0.323],
                    13: [7.344E-06, 1.202, 0.417],
                    14: [2.879E-05, 1.067, 0.358],
                    15: [2.792E-05, 0.805, 0.338],
                    16: [1.675E-05, 0.805, 0.338],
                    17: [7.914E-05, 0.764, 0.364],
                    18: [1.140E-05, 0.396, 0.298],
                    19: [1.140E-05, 0.396, 0.298],
                    20: [7.600E-04, 0.396, 0.298],
                    21: [7.992E-07, 1.578, 0.271],
                    22: [5.804E-06, 1.092, 0.377],
                    23: [1.836E-05, 0.395, 0.296],
                    24: [1.836E-05, 0.395, 0.296],
                    31: [6.845E-05, 0.130, 0.201],
                    32: [1.032E-05, 0.130, 0.201],
                    33: [3.802E-05, 0.644, 0.312],
                    34: [1.537E-04, 0.170, 0.208],
                    35: [3.326E-04, 0.148, 0.215],
                    41: [7.649E-06, 0.682, 0.279],
                    51: [1.588E-06, 1.154, 0.318],
                    52: [1.907E-06, 1.244, 0.335],
                    53: [2.485E-06, 1.196, 0.328],
                    54: [4.095E-06, 1.111, 0.316],
                    55: [1.135E-05, 0.974, 0.309],
                    56: [7.549E-06, 1.015, 0.324],
                    57: [3.363E-06, 1.296, 0.322],
                    61: [7.537E-04, 0.499, 0.266]}

def mcm_manual_rates(kin):
    """Check and return if contains manually written kinetic rate found in the dictionary: mcm_manual_lists"""
        
    for key in mcm_manual_lists:
        if key in kin: # Found pair
            val = mcm_manual_lists[key]
            if key == kin:
                if 'EXTRA 99' in val: return f'KINETIC {val} 1.0' # Add a ratio
                else: return f'KINETIC {val}'
            else: # Check if there is a factor
                res = kin.replace(key, '').strip()
                try: # Update factor - 1st case
                    factor = float(eval(res[1:]))
                except:
                    try: # Update factor - 2nd case
                        factor = float(eval(res[:-1]))
                    except:
                        print(f'MCM maunal can not process rate sting: {kin} with key: {key} and res: {res}')
                        return None
            # Process factor
            operator = res[0]
            if operator == '*': ratio = factor
            elif operator == '/': ratio = 1/factor
            else: raise ValueError(f'Found unrecognized operator {operator} in {kin}')

            if 'EXTRA 99' in val: return f'KINETIC {val} {ratio:6.3E}' # Add a ratio
            elif 'ARR' in val: # Update C1
                val = val.split(' ')
                ind = val.index('ARR')
                val[ind+1] = f'{float(val[ind+1])*ratio:6.3E}'
                return f'KINETIC {" ".join(val)}'
            else: raise ValueError(f'Need to manually update the rate {val} with factor {factor} in {kin}')
            

    return None

def mcm_photolysis_coefs(n):
    """
    Convert MCM photolysis index to photolysis coefficients in the format [l, m, n].

    These coefficients (read from function mcm_manual_rates()) can be used in the formula:
    J = l * np.cos(X)**m * EXP(-n * (1 / np.cos(X)))
    """

    #Get coefficients
    if n in photolysis_coefs:
        return photolysis_coefs[n]
    else:
        print(f'Non-recognized MCM photolysis index: {n}')
        return None

def mcm_to_ssh_photolysis(Jin):
    """Convert photolysis kinetic rate from MCM format to a format readable by SSH-aerosol."""
    
    # Check format
    if not Jin.startswith('J<'): return None
    
    # Extract photolysis index from Jin in the format "J<1>*0.03*2*..."
    values = Jin.split('*')
    
    # Photolysis index
    ind = values[0].split('J<')[1].split('>')[0]
    
    # Factor to photo;ysis rate
    factor = 1.0
    if len(values) > 1:
        for i in values[1:]: factor *= float(i)

    if ind.isdigit():
        tmp = mcm_photolysis_coefs(int(ind))
        if tmp:
            return f'KINETIC EXTRA 91 {tmp[0]} {tmp[1]:5.3f} {tmp[2]:5.3f} {factor:6.3E}'
        else:
            return f'KINETIC photolysis {Jin} with unknown index {ind}: !!!'
    else:
        print(f'Found unrecognized mcm photolysis index {Jin}')
        return f'KINETIC photolysis {Jin} with invalid index: {ind}!!!'

def parse_kinetic_rate(kin):
    """Extract C1 from the input string if apply the arrihenius law in the format k = C1 * T**C2 * exp(-C3/TEMP)"""
    
    # C, return C
    if isfloat(kin):
        return {'C1':float(kin)}
    # Generic rate coefficients
    elif kin in mcm_simple_rate_coefs:
        return {'GRC':kin}
    # Complex rate coefficients
    elif kin in mcm_complex_rate_coefs:
        return {'CRC':kin}
    # Third body
    elif kin in ['H2O','M','O2','N2','H2','RO2']:
        return {'TB':kin}
    # EXP(-Ea/TEMP), return Ea
    elif kin.startswith('EXP(') and kin.endswith('/TEMP)'):
        val = kin.split('EXP(')[1].split('/TEMP')[0]
        if isfloat(val): return {'C3':float(val)*-1}
        else: return {'UK':kin}
    # (TEMP/300)@(-n), return c2 = n, c1 = 1/(300 ** -n)
    elif kin.startswith('(TEMP/300)@(') and kin.endswith(')'):
        val = kin.split('TEMP/300)@(')[1].split(')')[0]
        if isfloat(val):
            return {'C2':float(val),'C1':1/(300.**float(val))}
        else: return {'UK':kin}
    # TEMP@(n), return n
    elif kin.startswith('TEMP@(') and kin.endswith(')'):
        val = kin.split('TEMP@(')[1].split(')')[0]
        if isfloat(val): return {'C2':float(val)}
        else: return {'UK':kin}
    else:
        return {'UK':kin}
    
def mcm_to_ssh_kinetic_rate(kin):
    """Convert a kinetic rate string from MCM format to SSH format."""

    # Remove space
    kin = kin.replace(' ', '')
    
    # Photolysis
    tmp = mcm_to_ssh_photolysis(kin)
    if tmp: return tmp
    
    # Manual sets
    tmp = mcm_manual_rates(kin)
    if tmp: return tmp

    # Break down to pieces
    klist = kin.split('*')
    df = {}
    for val in klist:
        val = val.strip()
        if val == '': continue
        tmp = parse_kinetic_rate(val)
        for key,value in tmp.items():
            # Record new element
            if key not in df: df[key] = value
            # Update C1
            elif key == 'C1': df[key] *= tmp[key]
            # Update C2 in T**c
            elif key == 'C2': df[key] += tmp[key]
            # Update C3 in exp(-c/T)
            elif key == 'C3': df[key] += tmp[key]
            # Update unknown
            elif key == 'UK': df[key] = '\t'.join([df[key],tmp[key]])
            else: print(f'Found multple values for the same key {key} that can not be merged: {df}, {tmp}, from {kin}')
                
    # Check pieces
    rate_ssh = 'KINETIC '
    if 'UK' in df: 
        print(f'Found unknown elements in {kin}: {df["UK"]}')
        return rate_ssh + "!!!"

    # if exists GRC, CRC, should not have C2, C3
    tmp = set(df.keys()) - set({'C1', 'TB'})
    if len(tmp) > 1 and tmp != {'C2', 'C3'}:
        print(f'Found multiple keys in {kin}: {df}, {tmp}')
        return rate_ssh + "!!!"

    # Build output in the SSH format
    if 'C1' in df: C1 = f"{df['C1']:6.3E}"
    else: C1 = 1
    
    # Third body or RO2-RO2
    if 'TB' in df:  
        if df['TB'] == 'RO2': rate_ssh += f"RO2 1 "
        else: rate_ssh += f"TB {df['TB']} "
    
    # MCM specific rates
    if 'GRC' in df: return rate_ssh + f"EXTRA 92 {mcm_simple_rate_coefs[df['GRC']]} {C1}"
    if 'CRC' in df: return rate_ssh + f"EXTRA 93 {mcm_complex_rate_coefs[df['CRC']]} {C1}"
    
    # Arrihenis
    if 'C2' in df: C2 = f"{df['C2']:6.3E}"
    else: C2 = 0
    if 'C3' in df: C3 = f"{df['C3']:6.3E}"
    else: C3 = 0
    return rate_ssh + f"ARR {C1} {C2} {C3}"

def gecko_to_ssh_kinetic_rate(rate_str):
    """Convert a kinetic rate string from GECKO format to SSH format."""
            
    RO2s = ['PERO1', 'PERO2', 'PERO3', 'PERO4', 'PERO5', 'MEPERO', 'PERO7', 'PERO8', 'PERO9']
    # Update rate_type and rate_info
    if 'T:' in rate_str:
        tmp = rate_str.split('T:',1)[1].split(' ')
        rate_type = tmp[0].strip()
        rate_info = [i for i in tmp[1:] if i != '']
    else: rate_type, rate_info = None, []

    # Get arrhenius's law coefficients
    arrhenius = []
    for i in rate_str.split(' '):
        if not isfloat(i): continue
        arrhenius.append(i.strip())
        if len(arrhenius) == 3: break
    # Build ssh kinetic
    rate_ssh = 'KINETIC '

    # Check type
    if rate_type:
        if rate_type in gecko_to_ssh_keywords.keys():
            rate_ssh += gecko_to_ssh_keywords[rate_type]
        elif rate_type in RO2s: # Check RO2 pool index
            rate_ssh += f"{gecko_to_ssh_keywords['PERO']}{RO2s.index(rate_type)+1} "
        else:
            raise NameError(f'Not found kinetic rate type: {rate_type}')

    # Photolysis
    if rate_type == 'HV': # no arrhenius(3) is used
        # Check read len
        if len(rate_info) != 2: 
            raise ValueError(f"Invalid parameters for 'HV' type: {rate_info}", ind0,fin[ind0:ind0+2])

        rate_ssh += ' '.join(map(str, rate_info)) # index and ratio
   
    # FALL OFF with 11 paramters
    elif rate_type == 'FALLOFF': # arrhenius(3) + FALLOFF(7) + ratio = 11
        # Check read len
        if len(rate_info) != 7: 
            raise ValueError(f"Invalid parameters for 'FALLOFF' type: {rate_info}")
            
        # Write 11 params
        rate_ssh += ' '.join(map(str, arrhenius + rate_info + [1.0]))
        
    elif rate_type == 'EXTRA': # arrhenius(3) + EXTRA(?) + ratio = ?

        # Check type
        if rate_info[0] not in ['100','500','501','502','550']:
            raise NameError(f'READ EXTRA, type not found: {rate_info[0]}')
        else: # write type
            rate_ssh += f'{rate_info[0]} '
            
        # Write the rest of params
        rate_ssh += ' '.join(map(str, arrhenius + rate_info[1:] + [1.]))
        
    elif rate_type == 'ISOM': # as EXTRA 200
        
        # Check read len isom(5)
        if len(rate_info) != 5: 
            print(rate_info)
            raise ValueError("READ ISOM: ",ind0,fin[ind0:ind0+2],len(rate_info))
        
        # specific case: if rate info is [0,0,0,0,1], then it is a simple rate
        rate_val = [float(i) for i in rate_info]
        if rate_val == [0.,0.,0.,0.,1.]:
            rate_ssh = 'KINETIC ARR ' + ' '.join(map(str, arrhenius))
        else: 
            # Build rate to 11 params: '200' + arrhenius(3) + isom(5) + 1.0 = 10
            rate_ssh += ' '.join(map(str, arrhenius + rate_info+[1.0]))
            #print('ISOM find complex rate: ',rate_ssh)
    else: 
        # Write 3 params
        rate_ssh += ' '.join(map(str, ['ARR'] + arrhenius))
                
    return rate_ssh

