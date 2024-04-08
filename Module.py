# ================================================================================
#
#   GENOA v3.0: the GENerator of reduced Organic Aerosol mechanism
#
#    Copyright (C) 2023 CEREA (ENPC) - INERIS.
#    GENOA is distributed under GPL v3.
#
# ================================================================================

"""
  Module.py contains three classes: Species, Reaction, and Kinetic. 

  These classes are utilized to handle chemical species and reactions in GENOA.

"""

import json
import numpy as np
from collections import defaultdict

import MolProperty as mp
from Convert import n_soap_group, BaseSpecies, psat_svoc, psat_nvoc, vpType, isfloat

# For species output
species_attr = [ 'name', 'status', 'string', 'mass', 'organic', 'condensable', 'Radical', 'RO2', \
                 'formula', 'SMILES', 'functionalGroups', 'SOAPStructure', 'ratios', 'DU', \
                 'non_volatile', 'psat_atm', 'psat_torr', 'dHvap_KJ','gamma','henry','Kp', \
                 'generation', 'reductions', 'precursors', 'type', 'note', 'groupID']

class Species:
    """Class to store and manage information related to chemical species"""
    def __init__(self, name=None, strin=None, mass=0.0):
        """Initialize attibutes for a Species object"""
        self.status = True   # If it is activated 
        self.string = strin  # Input format (Can be MechGen or GECKO-A format)
        self.note = None     # Notes
        # Species name and mass
        self.name= name        # Species name
        self.mass= mass        # Molar mass in mol/g
        
        # Properties - for all species
        self.organic=False     # if it is organic
        # For organics
        self.RO2=False         # if it is a peroxy radical (consider in the RO2 pool)
        self.Radical=False     # if it is a radical
        self.condensable=False # if it is condensable

        # Structure info
        self.formula='-'    # Formula as C[i]H[i]N[i]O[i]
        self.SMILES ='-'    # SMILES structure
        self.functionalGroups={} # Functional group info (used to define lumpable species)
        self.SOAPStructure = None # No.non-zero function groups in SSH-aerosol vector format
        self.ratios=None    # Atomic ratios: OM/OC, N/C, H/C, O/C
        self.DU=None        # degree of unsaturation

        # Aerosol Properties - update during reduction
        self.non_volatile = 0 # if it is a non-volatile compound
        self.psat_torr = 0.0  # Saturation vapor pressure at Tref (torr)
        self.psat_atm = 0.0   # Saturation vapor pressure at Tref (atm)
        self.dHvap_KJ = 0.0   # Enthalpy of vaporization (kJ/mol)

        # Aerosol Properties - recompute in SSH-aerosol
        self.gamma = None    # Activity coefficient at infinite dilution in water
        self.henry = 0.0     # Henry's Law constant (M/atm, mol/atm*l)
        self.Kp = None       # partitioning coefficient

        # Info related to its location in the reaction list
        self.generation = -1   # No.generation
        self.precursors = None # Precursors where it derived from 

        # Reduction related
        self.reductions = {} # Records
        self.type = None # Only species within the same type can be merged

        # SSH output - see SSH-aerosol aerosol list for more info
        self.groupID = 2 # group ID: organic
        
    def update_attribute(self, key, value, check=False):
        """Adds a new attribute to the object or updates an existing one."""
        if check:
            if not hasattr(self, key): return
                #raise AttributeError(f'{self.name} does not have the attribute "{key}".')
        setattr(self, key, value)
        return True
            
    def write_to_text(self, attrs = species_attr):
        """Converts selected object attributes to text."""
        attr_vals = {s: getattr(self, s) for s in attrs}
        return '\n'.join([f'{s}\t{json.dumps(attr_vals[s])}' for s in attr_vals if attr_vals[s]])

    def update_with_dict(self, new_dict, check=True):
        """Check species properties based on input dictionary. It is used in lumping."""

        # Update from input dict
        for key, val in new_dict.items():
            self.update_attribute(key, val, check=check)

        # Update psat_torr
        if self.condensable: self.psat_torr = self.psat_atm * 760.

        # Update ratios and DU
        val = {i: self.functionalGroups[i] for i in ['C', 'H', 'N', 'O']}
        self.ratios['OM/OC'] = round(self.mass / (val['C'] * 12.), 3)
        self.ratios['H/C'] = round(val['H'] / val['C'], 3)
        self.ratios['O/C'] = round(val['O'] / val['C'], 3)
        self.ratios['N/C'] = round(val['N'] / val['C'], 3)
        # DU: 1/2 (2 + 2C + N - H - X)
        self.DU = 0.5 * (2 + 2 * val['C'] + val['N'] - val['H'])

        # Update flag/type info
        #self.update_flags() # Set Radical/SVOC - not used in lumping
        self.update_type()

    def update_from_text(self,lines):
        """ Updates object attributes from a list of text strings."""
        for line in lines:
            if not line.strip(): continue
            key, val_json = line.strip().split('\t', 1)
            val = json.loads(val_json)
            self.update_attribute(key, val, True)
        
    def update_based_on_smiles(self):
        """Updates various properties based on SMILES string. Only can be used before reduction."""
        
        # Check for active status
        if not self.status: return None

        # Check and clean SMILES
        if self.SMILES == '': raise ValueError(f'Read non SMILES for species {self.name}')
        elif self.SMILES == '-': raise ValueError(f'Species Update: No SMILES for the given species: {self.name}')
        elif '/' in self.SMILES or '\\' in self.SMILES:
            print(f'Find / or \\ in {self.name}\'s SMILES: {self.SMILES}')
            self.SMILES = self.SMILES.replace('/', '').replace('\\', '')

        # Update molecular properties if not already updated
        pymol = mp.get_pybelmol(self.SMILES)
        self.mass = pymol.molwt
        self.formula = pymol.formula
        self.ratios, self.DU = mp.get_organic_ratios_and_du(pymol)
        self.functionalGroups = mp.get_functional_groups(pymol)
        # Nof finished
        #self.SOAPStructure = mp.get_toSOAPStructure(pymol)

        # Update Psat related properties at 298 K
        self.psat_atm = mp.get_saturation_vapor_pressure(pymol, vpType)
        self.psat_torr = self.psat_atm * 760.  # unit conversion
        self.dHvap_KJ = mp.get_enthalpy_vaporization(pymol, vpType) / 1000.  # unit in KJ
        
        # Check for non-volatility
        if psat_nvoc and self.psat_atm < psat_nvoc:
            print(f'species {self.name} is set to non_volatile, with P_sat = {self.psat_atm} < {psat_nvoc}')
            self.non_volatile = 1

        # Not tested
        if 0:
            # activity coefficient
            self.gamma=mp.get_activity_coefficient_inf(pymol)
    
            # Additional calculations if gamma is available
            if hasattr(self, 'gamma') and self.gamma:
                # henry's law constant (M/atm)
                # H = 1000*760/(18.0*GAMMAinf*Psat)
                self.henry = 1000 / (18 * self.gamma * self.psat_atm)   # Henry's law constant
                self.Kp = mp.get_partitioning_coefficient(self, temp)  # Partitioning coefficient

        # Update type information
        self.update_flags()
        self.update_type()

    def update_flags(self):
        """Determines flags for targeted species."""
        # Check for radicals and RO2
        if '[O+]' in self.SMILES or '[O]' in self.SMILES or 'O.' in self.string:
            self.Radical = True

            # Check species name for radical types - !!! need to change !!!
            if any(i in self.SMILES for i in ['O[O]','[O]O']) or 'OO.' in self.string: self.RO2 = 1 # Default Pool 1

        # Check for condensable species
        self.update_condensable()
            
    def update_condensable(self, Psat_condense=psat_svoc):
        """Update the volatility of  a species based on SVOC setting. Caution!! It may change VOC to SVOC."""
        if not self.Radical and (self.psat_atm != 0.0 and self.psat_atm <= Psat_condense) : 
            self.condensable = True
        else: self.condensable = False

    def update_type(self, mode='none'):
        """update the type of species based on its characteristics."""
        if not self.organic:
            self.type = None
            return  # Skip inorganic species

        # Basic type
        itype = 'R' if self.Radical else 'A' if self.condensable else 'G'

        if self.Radical: # Radical
            if self.RO2: itype += 'O2'
            elif '[O+][O-]' in self.SMILES: itype += 'OO'
            elif '[O]' in self.SMILES: itype += 'O'
        else: # Non-radical species
            if mode == 'none':
                pass             
            elif mode == 'group':
                # Additional checks for functional groups
                itype += ' '.join(self.functionalGroups)
            elif mode == 'vector':
                if self.SOAPStructure != None:
                    itype += ' '.join([f'{self.SOAPStructure[i]}:.1f' for i in sorted(self.SOAPStructure.keys())])
            else:
                raise ValueError(f'Unknown mode {mode}')
        # Update
        self.type = itype
    
    ## Convert information to formats can be processed in the SSH-aerosol
    def SOAPStructure_to_ssh_str(self):
        """Converts the SOAPStructure attribute to a formatted string used for SSH-aerosol aerosol species list."""

        # Check if SOAPStructure is None or empty
        if not self.SOAPStructure: return "-"

        vals = [] # Init
        # Process SOAPStructure
        for i in range(n_soap_group):
            j = str(i)
            if j in self.SOAPStructure: num = self.SOAPStructure[j]
            else: num = 0.0
            vals.append(f'&{num:5.2E}')

        return ''.join(vals) # Reformat into a single string

    def to_ssh_aerosol_list(self, mode='smiles', firstRun=False):
        """Formats species data into a line for the SSH-aerosol aerosol species list."""

        # Check if it is an aerosol
        if not (self.organic and self.condensable):
            print(f'species {self.name} is not an aerosol.')
            return None
        
        # Structure
        if mode == 'vectors': struc = self.SOAPStructure_to_ssh_str()
        elif mode == 'smiles': struc = self.SMILES
        else: raise NameError(f'MD: Unknown mode: {mode}')
        
        # Define the content and formatting style
        content = [ f'P{self.name}',    # Aerosol name
                    '4',                # Organic in SSH-aerosol
                    f'{self.groupID:d}',# group ID: organic
                    f'{self.mass:.2f}', # MWs
                    '--' if firstRun else self.name, # Remove precursor if it is the 1st run
                    '687d0',            # coll_fac
                    '8.39d0',           # mole_diam
                    '30.D-03',          # surf_tens
                    '1.0',              # accomod
                    '1.30D-06',         # mass_dens
                    '0', #'1' if self.non_volatile else '0',  # non_volatile
                    'BOTH',             # partitioning
                    struc,              # smiles/vector
                    f'{self.psat_torr:5.2E}', # psat_torr
                    f'{self.dHvap_KJ:5.2E}',  # dHvap
                    '0.0', # Henry - compute in ssh-aerosol
                    '0.0'  # Tref - compute in ssh-aerosol
                    ]

        return '\t'.join(content) + '\n'
            
    ## For GECKO strings
    def gecko_str_to_smiles(self, tag_canonical=True):
        """Update SMILES structure if a GECKO structure is provided."""

        if self.string == None: 
            print(f'MD: GECKO structure is None for {self.name}, cannot be converted to SMILES.')
            return False

        self.SMILES = mp.gecko_to_smiles(self.string, tag_canonical)
        return True

class Reaction:
    """Class to store and manage information related to chemical reactions"""

    def __init__(self):
        """Initialize attibutes for a Reaction object"""
        self.status = True
        self.note = ''      # Input information

        self.reactants = defaultdict(int)    # List of reactants: ratio
        self.products =  defaultdict(float)  # List of products: ratio
        
        self.rt_fnc = None   # Dict of species that have a function for ratio

        # Kinetic rate
        self.rate = Kinetic() # Kinetic rate constant
        self.info = None      # For mergin

    def add_product(self,pd,rt):
        """Add a single product and its ratio."""
        self.products[pd] += rt

    def del_product(self,pd):
        """Delete a single product and its ratio."""
        self.products.pop(pd, None)

    def clean_products(self, lim=0.):
        """Delete products with a ratio lower than threshold"""
        for s in list(self.products.keys()):
            if self.products[s] <= lim: self.products.pop(s, None)
        
    def reset_species(self):
        """Empty reactants and products and their ratio."""
        self.reactants = defaultdict(float)
        self.reset_products()

    def reset_products(self):
        """Empty all product and its ratio."""
        self.products = defaultdict(float)
        self.rt_fnc = None

    def get_species(self, tag_sort=False):
        """Get reactant and product species."""
        species = set(self.reactants) | set(self.products)
        if tag_sort: species = sorted(species)
        return species
        
    def get_ratios(self, slist, update_rate=False, rate_inds=[0], get_fnc=False):
        """Get product ratios."""

        ns = len(slist) # No. products
        nr = len(rate_inds) # No. kinetic values
        
        ratios = np.zeros((ns, nr)) # Ouput ratios
        if get_fnc: fnc = {}

        if update_rate:
            rates = np.zeros(nr)
            for i, s in enumerate(self.rate.values):
                if i in rate_inds and s > 0.: rates[i] = s
        else: rates = np.ones(nr) # No kinetic rate as a value
        
        # Get basic ratios
        for i, s in enumerate(slist):
            if s in self.products:
                for j in range(nr): ratios[i][j] = self.products[s] * rates[j]

        # Update ratios based on function ratio, if applicable
        if self.rt_fnc is not None:
            for i, s in enumerate(slist):
                if s in self.rt_fnc:
                    rt_inds = self.rt_fnc[s]
                    rt_rate = self.rt_fnc['values']
                    if get_fnc: fnc[s] = [rt_rate[j] for j in rt_inds]
                    else:
                        for j in range(nr):
                            irt = 0.
                            for k in range(0, len(rt_inds), 2):
                                k0 = rt_rate[rt_inds[k]][rate_inds[j]]
                                k1 = rt_rate[rt_inds[k+1]][rate_inds[j]]
                                if k0 + k1 > 0.0: irt += k0/(k0+k1)
                            # Update ratio considering function ratio
                            ratios[i][j] *= irt
                            
        if nr == 1: ratios = ratios.reshape(-1) # reshape nr if need

        if get_fnc: return [ratios, fnc]
        else: return ratios

    def update_info(self):
        """Get info from this reaction to check."""
        # Reactants
        rcts = ' '.join(sorted(self.reactants))
        # Kinetic
        kinetics, ratios = self.rate.get_kinetic_info()
        rcn = [rcts] + kinetics
        self.info = '\t'.join(rcn + [f'{ratios:6.3E}'])

    def to_line_wo_kinetic(self, separator = '->'):
        """Output reactions without kinetics"""

        # Reactants
        if len(self.reactants) == 0:
            print(f'No reactant is found for reaction {self.note}')
            self.status = False
            return None
        parts = []
        for s, rt in self.reactants.items():
            for i in range(int(rt)): parts.append(s)
        outline = ' + '.join(sorted(parts)) + f' {separator} '
        
        # Products
        parts = []
        if self.rt_fnc is None:
            for s in sorted(self.products.keys()):
                if self.products[s] != 1.0: parts.append(f'{self.products[s]:6.3E} {s}')
                else: parts.append(s)

            return outline + ' + '.join(parts)
            
        else: # Some products are with function ratio
            knc = [] # list of kinetics use for ratio computation
            for s in sorted(self.products.keys()):
                rt = self.products[s]
                if s in self.rt_fnc and s not in BaseSpecies: # in format speies: k1, k2
                    tmp = [] # tmp list for k1, k2
                    for i in self.rt_fnc[s]:
                        iknc = self.rt_fnc['rates'][i]
                        if iknc not in knc: knc.append(iknc)
                        tmp.append(f'{knc.index(iknc)+1}')
                    
                    if len(tmp) % 2 == 0: # Check length
                        #if len(tmp) > 2: 
                            #print(f'!!! Find more than 1 function ratio for {s}: {tmp}', outline, parts)
                        # Write
                        for p in range(0, len(tmp), 2):
                            parts.append(f'{rt:6.3E} {tmp[p]} {tmp[p+1]} {s}')
                    else:
                        raise ValueError(f'Not recognize function ratio for {s}: {tmp}', outline)
                # No function
                elif rt == 1.: parts.append(s)
                else: parts.append(f'{rt:6.3E} {s}')

            return outline + ' + '.join(parts) + '\n' + '\n'.join(knc)

    def from_line_wo_kinetic(self, line, separator = '->'):
        """Reconstructs a reaction from a line without kinetics."""
        
        if separator not in line:
            print(f'!!! Invalid separator {separator} to parse a reaction: {line}')
            return

        rcline, pdline = line.split(separator, 1)
        # Reset reactants and products
        self.reset_species()

        # Record reactants
        parts = [s for s in rcline.strip().split(' ') if (s != '' and s != '+')]
        for s in parts: self.reactants[s] += 1

        if len(self.reactants) == 0: # No read reactants
            print(f'!!! No reactant can be read from line: {line}')
            self.status = False
            return
 
        # Products
        parts = [i for i in pdline.strip().split(' ') if i != '']
        n = len(parts)
        if n == 0: return True # No product
        elif n == 1: # Only one product
            # Get species name
            s = parts[0]
            self.products[s] += 1.0
            return 

        # Multiple products
        sps = []
        for i, s in enumerate(parts+['+']): # Add + to find end
            if s == '+': # Find end
                n1 = len(sps)
                if n1 not in [1,2,4]:
                    raise ValueError(f'Find {n1} elements for one product part: {sps} in {line}. n1 should be 1,2,4. parts: {parts}', i, s)

                # Get species name
                isp = sps[-1]
                # Get ratios
                if n1 == 1: rt = 1.0
                else: rt = float(sps[0])
                
                if n1 == 4: # Find function ratio
                    if self.rt_fnc is None: self.rt_fnc = {'rates':[], 'values':[]}
                    # with more than one fratio
                    if isp in self.products:
                        if rt != self.products[isp]: raise ValueError(f'Find different ratio for {isp}: {rt} and {self.products[isp]} in {line}')
                    else: self.products[isp] = rt # Add ratio
                    #print(f'Find function ratio for {isp}: {sps} in {line}')

                    if isp not in self.rt_fnc: self.rt_fnc[isp] = []
                    self.rt_fnc[isp] += [int(i)-1 for i in sps[1:3]] # Add index
                else: self.products[isp] += rt # Add ratio
                sps = []  # Reset
            else: sps.append(s) # Find one element: ratio or species
        return

    ## Write reaction in SSH-aerosol format
    def to_SSH_rcn(self, mode='all'):
        """Transfer to the SSH output."""

        # Reactants and products
        rline = self.to_line_wo_kinetic() + '\n'

        if mode == 'simple': return rline

        # Kinetic comment
        rline += f'%{self.rate.str}\n'
        # Kinetic SSH-aerosol excutable format
        rline += f'{self.rate.ssh}'

        return rline

    def from_ssh_rcn(self,strin):
        """Build from the SSH string info."""

        # Get string for reaction and kinetic
        rcn_str, rate_str, comments, n = '','',[], 0
        for line in strin.split('\n'):
            line = line.strip()
            # Skip comment or empty line
            if not line: continue
            if line.startswith('%'): # Find comments
                if '===' not in line: # Record comments for kinetic
                    comments.append(line[1:]) # Remove %
            elif '->' in line: # Find initial reaction line
                rcn_str = line
                self.from_line_wo_kinetic(line)
            elif line.startswith('KINETIC'): # Read kinetic line
                n += 1 # Count
                if n > 1: # Find for ratio
                    if self.rt_fnc is None:
                        raise ValueError(f'Find multiple kinetic lines but the no sign for function ratio.', strin)
                    else: # Update rate 
                        self.rt_fnc['rates'].append(self.rate.ssh)
                    self.rate = Kinetic() # Reset
                rate_str = line
                self.rate.ssh = line
                self.rate.update_flags_from_ssh()
            elif line not in ['END']:
                raise ValueError(f'Can not understand line for load reactions: {line}')

        # Check and update updates
        if rcn_str == '' or rate_str == '':
            raise ValueError('Can not update reaction with a string: {strin}')

        self.status = True
        self.note = strin
        self.rate.str = '\n%'.join(comments) # Update rate string to preserve comments

class Kinetic:
    """
    Class to store and manage kinetic information.
    Supports MCM and GECKO-A string formats. SSH-aerosol format as default format.
    """

    def __init__(self, strin=None, sshin=None):
        """Initialize attibutes for a Kinetic object"""
        self.str = strin        # Orignial input kinetic information
        self.ssh = sshin        # SSH-aerosol format        
        
        self.Photolysis = False # Flag for photolysis reaction
        self.ratios = []        # List of lumping ratios
        self.values = None      # Kinetic rate computed based on current conditions (env_params)
        
        if sshin != None: self.update_flags_from_ssh()
        
    def update_flags_from_ssh(self):
        """ Update properties based on the ssh format"""
        
        # Update string
        if self.str is None: self.str = self.ssh
                    
        # Check Photolysis
        if "PHOT" in self.ssh or "EXTRA 91" in self.ssh or 'J<' in self.str: self.Photolysis = True

        # Check lumping ratios in the format: "Rs:ratio1*ratio2*..."
        if "\tlRs:" in self.str: 
            self.ratios.extend([float(i) for i in self.str.split('\tlRs:')[1].split('\t')[0].split('*')])

        # Check values in the format:
        if "\tK:" in self.str:
            self.values = [float(i) for i in self.str.split('\tK:')[1].split('\t')[0].split(',')]

    def update_by_ratio(self, ratio, mode='add'):
        """Update rate with a reduction ratio"""

        # Primary check ratio range
        if ratio == 1 and mode == 'add': return None
        elif ratio == 0.0:
            self.status = False
            return None
        elif ratio < 0.0:
            raise ValueError(f'Ratio should not be less than 0.')
            return None

        # Fresh rate
        self.update_flags_from_ssh()

        # Update
        # if self.values is None: print(f'No kinetic value for adding ratio {ratio}: self.to_SSH_kinetic()')
        if mode == 'add':
            if self.values: self.values = [i * ratio for i in self.values]
            self.ratios.append(ratio)
        elif mode == 'new':
            if self.values: self.values = [ratio for i in self.values]
            self.ratios = [ratio]
        else:
            raise ValueError(f'Not recognize mode {mode}')
        
        # Update string
        self.update_ssh_by_ratio(ratio, mode) # ssh
        self.update_str()

    def update_str(self):
        """Update string."""
        ratios_str = 'lRs:'+'*'.join([f'{i:6.3E}' for i in self.ratios])
        if self.values: k_str = "K:{','.join([f'{i:6.3E}' for i in self.values])}"
        else: k_str = ''
        new_str = ''
        for val in self.str.split('\t'):
            if "lRs:" in val: new_str += ratios_str
            elif "K:" in val: new_str += k_str
            else: new_str += val
        self.str = new_str

    def update_ssh_by_ratio(self, ratio, mode='add'):
        """Update the SSH format of kinetic ratio with a ratio"""
        
        # Check ratio
        if ratio == 1.0: return # No need to update
        elif ratio < 0.0: raise ValueError(f'Ratio {ratio} < 0.0 for reaction with kientic: {self.ssh}')
        
        # Get parts and label
        parts = [i for i in self.ssh.split(' ') if i != '']
        label = parts[1] # KINETIC ARR c1 c2 c3
        
        # Update label if there are two keyword
        if label in ['TB','RO2']: # KINETIC RO2 1 ARR c1 c2 c3
            label = parts[3]
        
        # Arrhenius - ratio the first
        if label == 'ARR': ind = parts.index(label) + 1
        # Check label for those with ratio at the tail
        elif label in ['PHOT','PHOTOLYSIS','EXTRA','FALLOFF']: ind = -1
        
        if mode == 'add': parts[ind] = f'{float(parts[ind])*ratio:6.3E}'
        elif mode == 'new': parts[ind] = f'{ratio:6.3E}'
        else: raise ValueError(f'Not recognize mode {mode}')

        # Update ssh string
        self.ssh = ' '.join(parts)

    def get_kinetic_info(self):
        "Output two lists of keys that uniquely identifies the type of a kinetic reaction."

        detailed_keys = [] # For output
        parts = [i for i in self.ssh.split(' ') if i != '']
        
        ilabel = 1 # Get 1st keyword
        if parts[ilabel] in ['TB','RO2']: # KINETIC RO2 1 ARR c1 c2 c3
            detailed_keys.append(' '.join(parts[ilabel:ilabel+2])) # Add 1st keyword
            ilabel = ilabel + 2 # Update

        # Arrhenius - ratio the first
        if parts[ilabel] == 'ARR': # ARR + nums
            detailed_keys.append(' '.join(['ARR']+[format(float(i),'.2e') for i in parts[ilabel+2:]]))
            ratio = float(parts[ilabel+1]) # C1
        # Check label for those with ratio at the tail
        elif parts[ilabel] in ['PHOT','PHOTOLYSIS','EXTRA','FALLOFF']:
            detailed_keys.append(' '.join([parts[ilabel]]+[format(float(i),'.2e') for i in parts[ilabel+1:-1]])) # Remove ratio
            ratio = float(parts[-1]) # ratio
        else: raise ValueError(f'Unknown keyword: ', ilabel, ' '.join(parts))
        # Output
        return [detailed_keys, ratio]
