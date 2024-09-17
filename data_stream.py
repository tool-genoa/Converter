import os
import re
from collections import OrderedDict

import module as md
from utils import isfloat
from parameters import species_list_aer_init, n_soap_group, basicsps_dict, BaseSpecies
from kinetic_rate_to_ssh import mcm_to_ssh_kinetic_rate, gecko_to_ssh_kinetic_rate


def get_species_from_reactions(reactions):
    """Return a list of species name that appear in the reactions"""
    slist = set()
    for rcn in reactions:
        if rcn.status: slist.update(rcn.get_species())
    return slist
    
def read_chem_sets(reaction_file, species_file=None, reaction_type='SSH', species_type='SSH', aero_vfile=None):
    """Reads and processes chemical mechanisms from given reaction and species files."""
    # Process reaction list
    if reaction_type == 'SSH':
        reactions = process_ssh_reactions(reaction_file)
    elif reaction_type == 'GECKO':
        reactions, species = process_gecko_reactions_and_species(reaction_file)
    elif reaction_type == 'FACSMILE':
        reactions = process_facsmile_reactions(reaction_file)
    else:
        raise ValueError(f'Unknown reaction type {reaction_type}.')

    # Process species list
    if reaction_type != 'GECKO':
        rcn_list = get_species_from_reactions(reactions) # Get a list of species in reactions
        if species_type == 'MCM': 
            species = process_mcm_species(species_file, check_exist=rcn_list)
        elif species_type == 'SSH': 
            species = process_ssh_species(species_file, check_exist=rcn_list)
        else:
            raise ValueError(f'Unknown species type {species_type}.')

    ## SOAPStructure - currently read SOAP structure
    if aero_vfile:
        if os.path.exists(aero_vfile): update_SOAPStructure_from_file(species, aero_vfile)
        else: raise FileNotFoundError(f'Error: Not found aero_vfile: {aero_vfile}. Please check.')

    return reactions, species
    

def process_ssh_reactions(reaction_file):
    """Processes a reaction list in the SSH-aerosol format."""

    reactions = []
    separator = '%==='

    with open(reaction_file, 'r') as file: 
        parts = file.read().split(separator)[1:] # Remove basic reactions if any
        for part in parts:
            if '->' in part and 'KINETIC' in part: # Find a reaction with kinetic
                # Add a new reaction
                reactions.append(md.Reaction())
                rcn=reactions[-1]
                # Update reaction with ssh string
                rcn.from_ssh_rcn(separator + part)
        # Print to check
        #print(f'Read # {len(reactions)} reaction out of # {len(parts)} lines from ssh-aerosol reaction list: {reaction_file}.')
    return reactions

def process_gecko_reactions_and_species(geckopath, pvapfile='pvap.nan.dat', remove_empty_pero=True):
    """Processes reaction and species lists from a GECKO-A path."""

    # All kinetic types and codes
    rtypes = ['TBODY' , 'FALLOFF', 'EXTRA', 'HV', 'OXYGEN', 'PERO1', 'PERO2', 'PERO3', 'PERO4', 'PERO5', 'MEPERO', 'PERO7', 'PERO8', 'PERO9', 'ISOM']
    #notes = ['HABS', 'MJ19', 'RA07', 'OLDH', 'IUPA', 'LV09', 'RXEX']

    ## Check primary VOCs
    with open(os.path.join(geckopath, 'listprimary.dat'), 'r') as f: lines =  f.read().splitlines()
    new_vocs = [] # Save read precursors
    for line in lines:
        info = [i for i in line.strip().split(' ') if i != '']
        if len(info) == 2: new_vocs.append(info[0]) # Find primary VOCs

    ## Read dictionary to build the species list
    species, sps_info = [], {} # Save species objects and name for index
    ifile = os.path.join(geckopath, 'dictionary.out')
    print(f'\nReading species list from {ifile} ...')
    
    with open(ifile, 'r') as f: lines =  f.read().splitlines()
    # Organic species with 14 columns
    # Default species with 12 columns - should been recorded as default
    # Example:
    # 1N3009   CH2(OH)CH(O.)CH2(ONO2) 1.NO 3  136.1  1  3  6  1  5  0  0  0  0
    for line in lines:
        info = [i for i in line.split(' ') if i != '']
        if len(info) == 14:
            # [index]: Meaning
            # [0]: sname; [1]: GECKO completed structure; [2]: GECKO-A simplified structure
            # [3]: No.gen; [4]: mass; [5]: if 1: radical; 
            # [6]-[13]: atoms in dicts: C's, H's, N's, O's, S's, F's, Cl's, Br's, and .'s

            sname = info[0] # Species name
            if sname in sps_info: # Check if name is repeated
                raise NameError(f'Repeat species {sname} found in GECKO input file: {geckopath}/dictionary.out')

            # Build new species
            sps_info[sname] = len(species)
            species.append(md.Species(name=sname,strin=info[1]))

            # Update Species Properties
            isp = species[-1]
            isp.note = info[2] # GECKO-A simplified structure
            isp.mass = float(info[4]) # molar mass

            # Update generation - different definition for GECKO-A
            if info[3] != '-0': # default species
                i = int(info[3])
                if i < 0: raise ValueError(f'Find generation number < 0 for GECKO species {isp.name}. Check.')
                isp.generation = i+1 # +1 to separate from primaryVOCs
            isp.note += f'_{isp.generation:d}' # Update GECKO generation in the note 

            if info[2] == '-': continue # For some inorganics: CO, CO2
            
            # Only for organics
            isp.organic = True
            
            # Check if radical
            if 'O.' in isp.string: isp.Radical = True
            if 'OO.' in isp.string:
                # Not count for cerigee radicals
                if not any(i in isp.string for i in ['ZOO.', 'EOO.']): isp.RO2 = True
            if info[5] == '1' and not isp.Radical:
                raise ValueError('Check RO2 species',line)

            # SMILES
            if isp.Radical: isp.gecko_str_to_smiles(tag_canonical = False)
            else: isp.gecko_str_to_smiles() 
            
            # Structure
            tmp = {}
            for i,j in enumerate(['C','H','N','O']): tmp[j] = int(info[6+i])
            
            isp.formula = 'C{:d}H{:d}N{:d}O{:d}'.format(tmp['C'],tmp['H'],tmp['N'],tmp['O'])
            isp.ratios = {  'OM/OC':round(isp.mass/(tmp['C']*12.),3),
                            'H/C':round(tmp['H']*1./tmp['C'],3),
                            'O/C':round(tmp['O']*1./tmp['C'],3),
                            'N/C':round(tmp['N']*1./tmp['C'],3)}
            isp.DU = int(0.5*(2+2*tmp['C']-tmp['H']+tmp['N']))
            isp.functionalGroups = tmp
            #isp.precursors = primaryVOCs
    
    # Print to check
    print(f'Read # {len(species)} species from # {len(lines)} lines.')
    
    ## RO2 info from pero[i].dat files
  
    npero = 9 # Number of RO2 pool
    pero_empty = [] # Count for empty pero files
    print('Reading RO2 species lists ...')
    for i in range(npero):
        pero_file = f'pero{i+1}.dat'
        with open(os.path.join(geckopath, pero_file), 'r') as f: lines = f.read().splitlines()
        icount = 0  # count for RO2 species
        for line in lines:
            line = line.strip()
            if line and line.startswith("G"): # Find species
                sname = line[1:]  # Extract species name without "G"
                # Check if exists
                if sname not in sps_info:
                    raise NameError(f'RO2 species {sname} from {pero_file} not found in dictionary.out')
                # Update RO2 index
                isp = species[sps_info[sname]]
                if not (isp.Radical and isp.RO2):
                    raise ValueError(f'RO2 species {sname} is not recognized as RO2 in dictionary.out.')
                isp.RO2 = i+1
                icount += 1
        # Print to check
        if icount > 0: print(f'Read # {icount} RO2s from {pero_file}.')
        else: pero_empty.append(f'PERO{i+1}')
    
    # PERO6 -> MEPERO - Update index for CH3O2
    if 'CH3O2' in sps_info: species[sps_info['CH3O2']].RO2 = 6
    else: 
        print('CH3O2 not found in species list.')
        if 'PERO6' not in pero_empty: raise ValueError('CH3O2 not found in species list but PERO6 is not empty.')

    ## Read Psat
    ifile = os.path.join(geckopath, pvapfile)
    print(f'\nReading aerosol properties from {ifile} ...')
    with open(ifile, 'r') as f: lines = f.read().splitlines()
    icount = 0 # count for SVOCs
    for line in lines:
        info = [i for i in line.split(' ') if i != '']
        if len(info) == 3: # In format species name, Past, evap
            sname = info[0][1:] # Remove "G": GO03001   1.073E-01    38.2
            if sname not in sps_info:
                raise NameError(f'Condensable species {sname} from {pvapfile} is not found in dictionary.out.')
            # Get speices
            isp = species[sps_info[sname]]
            # Update aero info
            isp.condensable = True
            isp.psat_atm = float(info[1])
            isp.dHvap_KJ = float(info[2])
            isp.psat_torr = isp.psat_atm * 760.
            isp.update_condensable()
            icount += 1
    # Print to check
    print(f'Read # {icount} SVOCs.')

    ## Add inorganics
    for s in BaseSpecies:
        if s in sps_info: continue
        sps_info[s] = len(species)
        species.append(md.Species(name=s, mass=basicsps_dict[s]))

    ## Read reaction list
    reactions, npero = [], 0
    ifile = os.path.join(geckopath, 'reactions.dum')
    print(f'\nReading reactions from {ifile} ...')

    with open(ifile, 'r') as f: lines = f.read().splitlines()
    for iline, line in enumerate(lines):
        line = line.strip()
        # No process for certain conditions
        if line == '' or line[0] != 'G' or any([f'{i} ' in line for i in ['AIN','AOU','WIN','WOU']]) or '=>' not in line: continue # Not find gas-phase reaction
        if remove_empty_pero and any([f'{i} ' in line for i in pero_empty]): 
            npero += 1
            continue
        
        info = line.replace('=>',' => ') # In case no space -> to have '=>' in the list
        info = [i for i in info.split(' ') if (i != '' and i != '+')]  # Pre-treat line

        # Get treatable elements
        parts = []
        for i in info:
            if '+' in i and not isfloat(i): # Find 'GNO2+GNO3', not '4.670E+14' !
                for j in i.split('+'):
                    if j != '': parts.append(j)
            else: parts.append(i)
                
        ind = parts.index('=>') # Index to separate reactants and products

        # Build new reaction
        reactions.append(md.Reaction())
        rcn=reactions[-1]

        # Reactants
        rtype = []
        for i in parts[:ind]:
            if i[0] == 'G': # Find species name
                s = i[1:] # Species name 
                if s in sps_info: rcn.reactants[s] += 1 # Add reactants
                else: raise ValueError(f'Find species {s} in {line} but not in dictionary.') 
            else: # Add for further analysis
                rtype.append(i)
                if i not in rtypes:
                    raise NameError(f'Find new type {i} not in rtypes from reaction: {line}')

        # Check number of rtype
        if len(rtype) > 1: raise NameError(f'Find multiple rtype: {rtype} for reaction {line}')

        # Products
        rt = 1.0 # default ratio
        for i in parts[ind+1:-3]: # Remove arrhenius coefs
            if i == 'NOTHING': continue # No products
            elif isfloat(i): rt = float(i) # Get number
            elif i[0] == 'G': # Get species name
                s = i[1:]
                if s in sps_info:
                    rcn.products[s] += rt
                    rt = 1.0 # Reset ratio
                else: raise ValueError(f'Find species {sname} in {line} but not in dictionary.')
            else: raise ValueError(f'Find unknown product {i} started not with G in {line}.')

        ## Kinetic rate
        # Prepare string
        rate_str = ' '.join(parts[-3:])
        # Read more paramters to the string if need
        if rtype != []:
            itype = rtype[0]
            rate_str += f' T:{itype} '
            if itype in ['HV','FALLOFF','EXTRA','ISOM']:
                rate_str += lines[iline+1].replace('/','').replace(itype,'').strip()

        # Process string
        rate_ssh = gecko_to_ssh_kinetic_rate(rate_str)
        # Update
        rcn.rate.str = rate_str
        rcn.rate.ssh = rate_ssh
        rcn.rate.update_flags_from_ssh()
        rcn.update_info()
        
    # Print to check
    print(f'Read # {len(reactions)} reactions. Removed # {npero} reactions with empty pero files: {pero_empty}.')

    # Reorder species list based on generation before output
    species.sort(key=lambda s: (s.generation, s.name))
    
    return reactions, species

def process_facsmile_reactions(reaction_file):
    """Processes a reaction list in the FACSMILE format.""" 

    with open (reaction_file,'r') as f: lines = f.read().splitlines() # mcm_export.fac

    reactions = [] # Save reactions
    for line in lines:
        line = line.strip()
        if line.startswith('%') and '=' in line: # Find reaction e.g., % 3.31D-11 : C722OOH + OH = C722O2 ;
            parts = [i.strip() for i in line.replace('%','').replace(';','').split(':')]
            if len(parts) != 2:
                raise ValueError(f'Find reaction that can not be processed: {line}')

            # Build new reaction
            reactions.append(md.Reaction())
            rcn = reactions[-1]
            # Reactants and products
            rcn.from_line_wo_kinetic(parts[1], separator=' =')
            # Check if need to update notation
            if 'D' in parts[0]:
                # Replace D+ with E+, D- with E-
                parts[0] = parts[0].replace('D+', 'E+').replace('D-', 'E-')
                # Replace 'D' with 'E', but not if it's part of 'KDEC'
                if 'D' in parts[0]: # Still find D, change to E+ if not followed by EC
                    parts[0] = re.sub(r'D(?!EC)', 'E+', parts[0])
            # Record kinetic
            rcn.rate.str = parts[0]
            rcn.rate.ssh = mcm_to_ssh_kinetic_rate(rcn.rate.str)
            # Check ssh string
            if '!!!' in rcn.rate.ssh:
                print('Manual convertion is required for:', line)
                #raise ValueError(f'Find !!! in ssh string: {rcn.rate.ssh}')
            # Update kinetic flags
            rcn.rate.update_flags_from_ssh()

    # Print to check
    print(f'Read # {len(reactions)} reactions from the FACSMILE reaction list {reaction_file}.')

    return reactions


def process_ssh_species(species_file, check_exist=None):
    """Processes MCM species from a species file and returns a list of Species instances."""

    species = []

    with open (species_file, 'r') as f: lines = f.read().splitlines()
    inds = OrderedDict()
    for i, line in enumerate(lines):
        if 'name\t"' in line:  # Find s species
            # For last species
            end_idx = i-1
            if len(inds) != 0: inds[sname].append(end_idx)
            # For new species
            start_idx = i
            sname = line.split('"',1)[1].replace('"','').strip()
            if sname in inds: raise ValueError(f'Find repeat species {sname} in {species_file}')
            else: inds[sname] = [start_idx]
    # Add last index
    if inds: inds[sname].append(len(lines)-1)
    else: raise ValueError(f'Can not read species info from {species_file}.')

    # Check if contains all species need
    if check_exist is not None:
        sps_list = set(inds.keys())
        sps_diff = check_exist - sps_list - BaseSpecies
        if sps_diff:
            raise NameError(f'DS: species {sps_diff} not found in the list {species_file}')
        # remove unused species
        sps_diff = sps_list - check_exist - BaseSpecies
        for sp in sps_diff: del inds[sp]
    
    # Process species
    for sp in inds: # Find in species list
        # Get ranges for inpput properties
        irange = inds[sp]
        # Add new species
        isp = md.Species(name=sp)
        isp.update_from_text(lines[irange[0]:irange[1]])
        species.append(isp)
    # Add BaseSpecies
    for sp in BaseSpecies:
        if sp not in inds:
            species.append(md.Species(name=sp, mass=basicsps_dict[sp]))
            
    return species

def process_mcm_species(species_file, check_exist=None):
    """
        Processes MCM species from a species .tsv file and returns a list of Species instances.
        Current format: \t is the separator, with 9 columns
        Name	Smiles	Inchi	InchiKey	Formula	Mass	Excited	PeroxyRadical	Synonyms

    """

    print(f'\nReading species list from {species_file} ...')
    species=[] # md.Species() instances

    # Process species properties
    spLists = OrderedDict() # Species name & properties
    ncol = 0 # No.columns
    with open (species_file,'r') as f: lines = f.read().splitlines()
    for line in lines:
        line = line.strip()
        if line == '' or line[0] == '*': continue # No line or with comments
        if line.startswith('Name'): # Find the header
            if ncol == 0: ncol = len(line.split('\t'))
            else: raise ValueError(f'DS: Find multiple headers in {species_file}.')
            if ncol != 9: raise ValueError(f'DS: Find {ncol} columns in {line}. Default is 9, please check {species_file}.')
        else: # Read species info
            parts = [i.strip() for i in line.split('\t') if i != '']
            if len(parts) < 6:
                print(f'!!! Not read line with {len(parts)} columns < {ncol}: {line}')
                continue
            
            # Add new species, read SMILES, InChI, InChIKey (no), formula, mass, excited, peroxyRadical, synonyms (no)
            read_info = [parts[1], parts[4], float(parts[5])] # smiles, formula, mass
            spLists[parts[0]] = read_info

    # Remove unused species if needed
    if check_exist is not None:
        sps_list = set(spLists.keys())
        # Check if species list contains all species need
        sps_diff = check_exist - sps_list - BaseSpecies
        if sps_diff:
            raise NameError(f'DS: species {sps_diff} not found in the list {species_file}')
        # Remove unused species
        sps_diff = sps_list - check_exist - BaseSpecies
        for sp in sps_diff: del spLists[sp]

    # Record in species
    for sp in spLists:
        # Add new species
        species.append(md.Species(name=sp, mass=spLists[sp][2]))
        isp = species[-1]
        if 'C' not in spLists[sp][1]: 
            print(f'Read inorganic: {sp}')
            continue # inorganic

        isp.organic = True
        isp.SMILES = spLists[sp][0]
        isp.string = '' # No additional string
        isp.note = 'MCM'
        isp.formula = spLists[sp][1] #- will be updated later based on SMILES
        isp.mass = spLists[sp][2] #- will be updated later based on SMILES
        if 'C' in isp.SMILES and 'Cl' not in isp.SMILES:
            print(isp.SMILES)
            isp.update_based_on_smiles()

        # Check name
        #if sp.endswith('O2') and not isp.RO2:
        #    print(f'Warning: Find species {sp} is not RO2: {isp.SMILES}.')
    print(f'Read # {len(species)} species.')

    # Add BaseSpecies if needed
    for sp in BaseSpecies:
        if sp not in spLists:
            species.append(md.Species(name=sp, mass=basicsps_dict[sp]))

    return species

def update_SOAPStructure_from_file(species, filename, mode='SOAP', ngroup = n_soap_group):

    if not (filename and os.path.exists(filename)):
        raise ValueError(f'Can not read SOAPStructure from {filename}')
        return False

    # Get species name as a list
    sps_info = {s.name: i for i, s in enumerate(species) if (s.status and s.organic)}

    # Record No. species found
    find_sp = []

    print(f'\nReading SOAPStructure from {filename} ...')
    with open (filename, 'r') as f: lines = f.read().splitlines()
    
    if mode == 'SOAP':
        init = -1 # init
        identifer = 'is constructed from smiles:' # !!! May need to be updated
        for i, line in enumerate(lines): # Read from the next line
            if identifer in line: 
                init = i
                break
        if init == -1:
            raise ValueError(f'Not find identifer \'{identifer}\' in file {filename}')

        strucs = {} # species: SOAPStructure
        # Read SOAP structure
        for line in lines[init:]:
            if identifer in line: # Find new species name
                sp = line.split(identifer)[0].strip()
                if sp in strucs: raise ValueError(f'Repeat species {sp} found in {filename}')
                else: strucs[sp] = {}
            else: # Read SOAP structure
                parts = [i.strip() for i in line.split(' ') if i != '']
                if len(parts) >= 4 and parts[0].isdigit(): # !!! May need to be updated
                    ind, inum = parts[0], int(parts[-3])
                    if ind in strucs[sp]: raise ValueError(f'Repeat group {ind} found in {filename}')
                    else: strucs[sp][ind] = inum

        # Update SOAP structure
        for sp in strucs:
            if sp not in sps_info: continue
            #print(f'Update SOAPStructure for species {sp}: {strucs[sp]}')
            species[sps_info[sp]].SOAPStructure = strucs[sp]
            find_sp.append(sp)

    elif mode == 'SSH': # SSH-aerosol aerosol species list with vector
        for line in lines:
            if '&' in line:
                parts = [i for i in line.split('&') if isfloat(i)] 
                if len(parts) == ngroup: # Find vector groups
                    # Get species name
                    sp = line.strip().split('\t',1)[0][1:] # Remove "P"
                    if sp in sps_info[sp]:
                        vals = {}
                        for i, val_str in parts:
                            val = float(val_str)
                            if val > 0.: vals[str(i)] = val
                        # Update SOAPStructure
                        isp = species[sps_info[sp][sp]]
                        isp.SOAPStructure = vals
                        find_sp.append(sp)

    elif mode == 'TABLE':
        for line in lines:
            parts = [i for i in i.split(',') if i != '']
            n = 1
            if len(parts) ==  ngroup + n: # Name + groups
                sname = parts[0].strip()

                # Check name
                if sname not in sps_info: continue

                # Update SOAPStructure
                vals = {}
                for i, val_str in parts[n: ngroup+n]:
                    val = float(val_str)
                    if val > 0.: vals[str(i)] = val
                species[sps_info[sp]].SOAPStructure = vals
                find_sp.append(sp)
    else:
        raise ValueError(f'Unrecognized input mode: {mode}')

    # Prepare for checking
    n_find_sp = len(find_sp)
    # No.condensables not found with SOAPStructure
    miss_sp = [s.name for s in species if s.condensable and s.name not in find_sp]

    # Print to check
    print(f'Updated # {n_find_sp} species with SOAP structure.')
    if miss_sp != []:
        print(f'Find # {len(miss_sp)} condensables with no SOAP atructure recorded in {filename}: {miss_sp}')

    return species

def to_ssh_sets(path, chem, reactions, species, tag_folder=True, tag_fake=False):
    """Return mechanism files in the SSH format"""

    # Check species list and reaction list
    sps_info = {s.name: i for i, s in enumerate(species) if s.status}
    nsps = len(sps_info) # Count No. species
    
    ## Update reduction list - Add fake radicals in reaction list
    if tag_fake: # symbol = 'FA'
        fake_species = {}
        for rcn in reactions:
            if not rcn.status: continue
            for s in list(rcn.products.keys()):  # Allow change dict in loop
                if species[sps_info[s]].Radical: # Add Fake species
                    sname = 'FA'+ s
                    rt = rcn.products[s]
                    rcn.products[sname] += rt
                    if sname not in fake_species:
                        fake_species[sname] = species[sps_info[s]].mass

    # Create output folder
    if tag_folder: ipath = os.path.join(path, chem)
    else: ipath = path
    os.makedirs(ipath, exist_ok=True)

    # Save reaction list for all mode
    RO2s_used, nRO2_used = set(), 0 # Record RO2 species with RO2-RO2 reactions
    with open(os.path.join(ipath, f'{chem}.reactions'), 'w') as f:
        nrcn = 0
        separator = '=' * 10
        for rcn in reactions:
            if not rcn.status: continue
            nrcn += 1 # Count
            if 'RO2 ' in rcn.rate.ssh: RO2s_used.update(rcn.reactants) # Update used RO2 pool
            f.write(f'%{separator}{nrcn}{separator}\n{rcn.to_SSH_rcn()}\n')
        f.write('END\n')

    ## Save species files based on mode
    RO2s = set()
    # Write species files
    with open(os.path.join(ipath, f'{chem}.mol'), 'w') as fmol, \
     open(os.path.join(ipath, f'{chem}.species'), 'w') as fsps, \
     open(os.path.join(ipath, f'{chem}.aer.vec'), 'w') as faero:
        # Write headers for fsps and faero
        fsps.write(f'# species list of {chem}: {nsps} gas-phase species\n')
        with open (species_list_aer_init,'r') as f: faer0 = f.readlines()
        faero.write(''.join(faer0[:10]))

        with open(os.path.join(ipath, f'{chem}.aer.1st'), 'w') as faero_1st, \
             open(os.path.join(ipath, f'{chem}.aer'), 'w') as faero_sml:
            faero_1st.write(''.join(faer0[:10]))
            faero_sml.write(''.join(faer0[:10]))
            for isp in species:
                fsps.write(f'{isp.name}    {isp.mass:.2f}\n')
                fmol.write('\n' + isp.write_to_text() + '\n')
                if isp.condensable:
                    faero.write(isp.to_ssh_aerosol_list(mode='vectors'))
                    faero_1st.write(isp.to_ssh_aerosol_list(mode='smiles', firstRun=True))
                    faero_sml.write(isp.to_ssh_aerosol_list(mode='smiles'))
                elif isp.RO2: RO2s.add(isp.name)
            # Add water
            faero.write(faer0[-1])
            faero_1st.write(faer0[-1])
            faero_sml.write(faer0[-1])

        # Add fake species if need
        if tag_fake and fake_species:
            for s in fake_species: fsps.write(f'{s}    {fake_species[s]:.2f}\n')

    # Check and write RO2 species
    RO2s_used = RO2s & RO2s_used
    nRO2, nRO2_used = len(RO2s), len(RO2s_used)
    with open(os.path.join(ipath, f'{chem}.RO2'), 'w') as f:
        # Hearder
        f.write(f'# {nRO2_used} out of {nRO2} RO2 species with RO2 reactions.\n')
        # Content
        for s in RO2s_used: f.write(f'{s}    {species[sps_info[s]].RO2:d}\n')

    n = len(reactions)
    if nrcn != n: print(f'Write # {nrcn} out of # {n} reactions. Find # {n-nrcn} invalid reactions in the list.')
    n = len(sps_info)
    if nsps != n: print(f'Write # {nsps} out of # {n} species. Find # {n-nsps} invalid species in the list.')
    if nRO2_used: print(f'Write # {nRO2_used} out of # {nRO2} RO2 species with RO2-RO2 reactions.') #'RO2 without RO2-RO2: ', RO2s - RO2s_used)

    return None


