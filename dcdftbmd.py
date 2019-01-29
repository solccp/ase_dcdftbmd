"""This module defines an ASE interface to DC-DFTB-MD

http://www.chem.waseda.ac.jp/dcdftbmd/

sol.chou@gmail.com
v20190129


"""

import os

import numpy as np

from ase.calculators.calculator import FileIOCalculator
from ase.units import Hartree, Bohr

import io
import re


class DCDFTBMD(FileIOCalculator):
    """ A DC-DFTB-MD calculator with ase-FileIOCalculator

        Example:
            dcdftb_keywords = \"\"\"
            SCC=(DAMPXH=TRUE DAMPXHZETA=4.00 THIRDFULL=TRUE)
            DISP=(DISPTYPE=5)
            DC=FALSE
            MISC=(FORCE=TRUE)
            PBC=(STRESS=TRUE)
            \"\"\"

            dcdftb_skinfo = \"\"\"
            2
            O 2 -0.1575
                O-O.skf O-H.skf
            H 1 -0.1857
                H-O.skf H-H.skf
            \"\"\"

            atoms.set_calculator(DCDFTBMD(keywords=dcdftb_keywords, skinfo=dcdftb_skinfo))

    """
    if 'DCDFTBMD_COMMAND' in os.environ:
        command = os.environ['DCDFTBMD_COMMAND']
    else:
        command = 'dcdftbmd.x'

    implemented_properties = ['energy', 'forces', 'stress']

    def __init__(self, keywords ,skinfo, label='dcdftbmd', atoms=None,  **kwargs):
        """Construct a DC-DFTB-MD calculator.

        """

        if 'DCDFTBMD_SKPATH' in os.environ:
            self.slako_dir = os.environ['DCDFTBMD_SKPATH'].rstrip('/') + '/'
        else:
            self.slako_dir = './'

        
        self.lines = None
        self.atoms = None
        self.atoms_input = None
        self.outfilename = 'dftb.out'

        FileIOCalculator.__init__(self, label=label, atoms=atoms, **kwargs)

        spinpol = False
        self.nspin = 2 if spinpol else 1
        self.keywords = keywords
        self.skinfo = skinfo

        
    def write_dftb_in(self, filename, atoms):
        """ Write the input file for the dftb calculation.
        """

        with open(filename, 'w') as f:    
            print(self.keywords.strip(), file=f)

            print("", file=f)

            print("title", file=f)

            print("", file=f)

           
            sio = io.StringIO(self.skinfo.strip())
            nelem = int(next(sio))
            print(nelem, file=f)
            for i in range(nelem):
                line = next(sio)
                print(line, end='', file=f)
                skline = next(sio)
                sknames = skline.split()
                final_names = []
                for name in sknames:
                    if not name.startswith('/'):
                        final_names.append(self.slako_dir + name)
                    else:
                        final_names.append(name)
                print('    ' + ' '.join(final_names), file=f)

            print("", file=f)

            nat = len(atoms)
            charge = 0
            spin = 1
            print(f"{nat} {charge} {spin}", file=f)
            ispbc = atoms.get_pbc()
            
            box = atoms.get_cell()

            # print(box)
            for i in atoms:
                print('{:3s} {:20.12F} {:20.12F} {:20.12F}'.format(i.symbol, *i.position), file=f)


            if (any(ispbc)):
                for i in range(3):
                    print('TV {:20.12F} {:20.12F} {:20.12F}'.format(box[i][0], box[i][1], box[i][2]), file=f)

            print("", file=f)
 

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()
        return changed_parameters

    def check_state(self, atoms):
        system_changes = FileIOCalculator.check_state(self, atoms)
        # Ignore unit cell for molecules:
        if not atoms.pbc.any() and 'cell' in system_changes:
            system_changes.remove('cell')
        return system_changes

    def write_input(self, atoms, properties=None, system_changes=None):
        from ase.io import write
        
        FileIOCalculator.write_input(
            self, atoms, properties, system_changes)
        
        
        self.write_dftb_in(os.path.join(self.directory, 'dftb.inp'), atoms)
        
        self.atoms_input = atoms
        self.atoms = None

    def read_charges(self):
        """Get partial charges on atoms
            in case we cannot find charges they are set to None
        """
        infile = open(os.path.join(self.directory, 'dftb.dat'), 'r')
        lines = infile.readlines()
        infile.close()

        qm_charges = []

        lit = iter(lines)
        found_charge = False
        for line in lit:
            if line.strip().startswith('MULLIKEN POPULATIONS AND NET ATOMIC CHARGES'):
                found_charge = True
                next(lit)
                for _ in range(len(self.atoms_input)):
                    line = next(lit)
                    arr = line.split()
                    qm_charges.append(float(arr[2]))

        if (not found_charge):
            return None
       
        return np.array(qm_charges)

    def get_charges(self, atoms):
        """ Get the calculated charges
        this is inhereted to atoms object """
        if 'charges' in self.results:
            return self.results['charges']
        else:
            return None

    def read_results(self):
        myfile = open(os.path.join(self.directory, 'dftb.out'), 'r')
        self.lines = myfile.readlines()
        myfile.close()

        self.atoms = self.atoms_input

        #read energy and force        
        p = re.compile(r'\s*Final .*DFTB.*\sEnergy = (?P<energy>.*) Eh    after\s+[0-9]+\siterations\s*')
        
        lit = iter(self.lines)
        have_stress = False
        for line in lit:
            
            m = re.match(p, line)
            if (m):
                self.results['energy'] = float(m['energy'])
            
            if line.strip() == 'Atom           F(x)                F(y)                F(z)':
                next(lit)
                gradients = []
                for _ in range(len(self.atoms_input)):
                    line = next(lit)
                    word = line.split()
                    gradients.append([float(word[k]) for k in range(1, 4)])
                
                self.results['forces'] = np.array(gradients) * Hartree / Bohr
        
            if (line.strip().startswith('Stress tensor')):
                #stress
                

                line = next(lit)
                if (line.strip() != '==============='):
                    continue

                have_stress = True
                stress = list()
                next(lit)
                next(lit)
                next(lit)

                for _ in range(3):
                    line = next(lit)
                    arr = line.split()
                    print(arr)
                    stress.append(float(arr[1]))
                    stress.append(float(arr[2]))
                    stress.append(float(arr[3]))
                
                if have_stress:
                    # stress = -np.array(stress) * Hartree / Bohr**3
                    stress = np.array(stress)
                    self.results['stress'] = stress.flat[[0, 4, 8, 5, 2, 1]]
                #stress ends

        #charge
        charges = self.read_charges()
        if charges is not None:
            self.results['charges'] = charges

        
        os.rename(os.path.join(self.directory, 'dftb.out'), os.path.join(self.directory, 'dftb.out.bak'))
        os.rename(os.path.join(self.directory, 'dftb.dat'), os.path.join(self.directory, 'dftb.dat.bak'))

    


