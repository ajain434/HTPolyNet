"""
.. module:: cgenff
   :synopsis: Manages execution of CGenFF parameter generation and conversion to GROMACS format
   
.. moduleauthor: [Your Name]
"""
import os
import logging
import shutil
from parmed.gromacs import GromacsTopologyFile
import parmed
from HTPolyNet.command import Command
from HTPolyNet.coordinates import Coordinates
import HTPolyNet.software as sw

logger = logging.getLogger(__name__)

# Get the path to the cgenff_to_gmx.sh script, from $SILCSBIODIR
SICS_BIO_PATH = os.environ.get('SILCSBIODIR', None)

def setup_charmm36(ff_dir):
    """Setup CHARMM36 force field directory
    
    :param ff_dir: Directory where CHARMM36 force field files should be extracted
    :type ff_dir: str
    """
    if not os.path.exists(ff_dir):
        os.makedirs(ff_dir)
    
    # Extract CHARMM36 force field if not already present
    if not os.path.exists(os.path.join(ff_dir, 'charmm36.ff')):
        c = Command(f'tar xzf charmm36-jul2022.ff.tgz -C {ff_dir}')
        c.run(quiet=False)

def create_atom_mappings(original_mol2):
    """Create atom name mappings for initial molecules, ensuring unique numbered names.
    
    :param original_mol2: Molecule coordinates object from mol2 file
    :return: Tuple of dictionaries (atom_name_map, residue_name_map, residue_num_map)
    """
    atom_name_map = {}
    residue_name_map = {}
    residue_num_map = {}
    
    # Track highest number for each atom type
    atom_counters = {}
    
    # First pass: find highest numbers for existing numbered atoms
    for atom_name in original_mol2.A['atomName']:
        base = ''.join(c for c in atom_name if not c.isdigit())
        number = ''.join(c for c in atom_name if c.isdigit())
        
        if number:
            atom_counters[base] = max(atom_counters.get(base, 0), int(number))
    
    # Second pass: create the mapping
    for idx, row in original_mol2.A.iterrows():
        atom_idx = idx + 1
        atom_name = row['atomName']
        base = ''.join(c for c in atom_name if not c.isdigit())
        number = ''.join(c for c in atom_name if c.isdigit())
        
        # If atom already has a number, keep it
        if number:
            atom_name_map[atom_idx] = atom_name
        else:
            # If no number, assign next available number
            if base not in atom_counters:
                atom_counters[base] = 0
            atom_counters[base] += 1
            atom_name_map[atom_idx] = f"{base}{atom_counters[base]}"
        
        # Map residue information
        residue_name_map[atom_idx] = row['resName']
        residue_num_map[atom_idx] = row['resNum']
    
    return atom_name_map, residue_name_map, residue_num_map


def create_original_mappings(original_mol2):
    """Create atom mappings preserving original names for reaction products.
    
    :param original_mol2: Molecule coordinates object from mol2 file
    :return: Tuple of dictionaries (atom_name_map, residue_name_map, residue_num_map)
    """
    atom_name_map = {}
    residue_name_map = {}
    residue_num_map = {}
    
    for idx, row in original_mol2.A.iterrows():
        atom_idx = idx + 1
        atom_name_map[atom_idx] = row['atomName']
        residue_name_map[atom_idx] = row['resName']
        residue_num_map[atom_idx] = row['resNum']
    
    return atom_name_map, residue_name_map, residue_num_map

def CGenFFParameterize(
    inputPrefix, outputPrefix, input_structure_format='mol2', **kwargs
):
    """CGenFFParameterize manages execution of cgenff_to_gmx.sh to generate 
    CGenFF parameters

    :param inputPrefix: basename of input structure file
    :param outputPrefix: basename of output files
    :param input_structure_format: format of input structure file, defaults to 'mol2'
    """
    structin = f'{inputPrefix}.{input_structure_format}'
    
    # Backup original input file since some tools may overwrite it
    new_structin = f'{inputPrefix}-orig.{input_structure_format}'
    if os.path.exists(structin):
        logger.debug(
            f'Backing up input {structin} to {new_structin}'
        )
        shutil.copy(structin, new_structin)

    # Get atom mapping from original mol2
    # if input_structure_format == 'mol2':
    #     original_mol2 = Coordinates.read_mol2(new_structin)
    #     # Create mappings from index to atom properties
    #     h_count = 1
    #     atom_name_map = {}
    #     residue_name_map = {}
    #     residue_num_map = {}
        
    #     for idx, row in original_mol2.A.iterrows():
    #         atom_idx = idx + 1  # enumerate from 1
    #         atom_name = row['atomName']
    #         res_name = row['resName']
    #         res_num = row['resNum']
            
    #         # Map atom names (with sequential H numbering)
    #         if atom_name.startswith('H'):
    #             atom_name_map[atom_idx] = f'H{h_count}'
    #             h_count += 1
    #         else:
    #             atom_name_map[atom_idx] = atom_name
            
    #         # Map residue information
    #         residue_name_map[atom_idx] = res_name
    #         residue_num_map[atom_idx] = res_num
    
    # Read original mol2 to get atom name mapping
    if input_structure_format == 'mol2':
        original_mol2 = Coordinates.read_mol2(new_structin)
        
        # Choose mapping function based on whether this is a reaction product
        if '~' in inputPrefix:
            # For reaction products, preserve original names
            atom_name_map, residue_name_map, residue_num_map = create_original_mappings(original_mol2)
        else:
            # For initial molecules, ensure unique numbered names
            atom_name_map, residue_name_map, residue_num_map = create_atom_mappings(original_mol2)
            

        # CGENFF requires no "=" in the filename
    # Create a clean version of the input file without special characters
    if '=' in inputPrefix:
        clean_prefix = inputPrefix.replace('=', '_')
    else:
        clean_prefix = inputPrefix
    clean_structin = f'{clean_prefix}.{input_structure_format}'
    if '=' in inputPrefix:
        shutil.copy(structin, clean_structin)

    # Run cgenff_to_gmx.sh
    cmd = Command(
        f'{SICS_BIO_PATH}/cgenff/cgenff_to_gmx.sh mol={clean_structin}',
    )
    cmd.run()
    
    logger.debug(f'CGenFF generated {clean_prefix}_gmx.top')
    logger.debug(f'CGenFF generated {clean_prefix}_gmx.pdb')
    
    # Check for output files
    pdb_file = f'{inputPrefix}_gmx.pdb'
    top_file = f'{inputPrefix}_gmx.top'
    
    if not (os.path.exists(f'{clean_prefix}_gmx.top') and os.path.exists(f'{clean_prefix}_gmx.pdb')):
        raise Exception("CGenFF parameterization failed to generate output files")
    
    if '=' in inputPrefix:
        # Now return the pdb and top to the original name after cgenff
        shutil.move(f'{clean_prefix}_gmx.top', top_file)
        shutil.move(f'{clean_prefix}_gmx.pdb', pdb_file)
        
    # Fix the PDB file atom names, residue names, and residue numbers
    fixed_pdb = f'{inputPrefix}_gmx_fixed.pdb'
    
    with open(pdb_file, 'r') as f_in, open(fixed_pdb, 'w') as f_out:
        for line in f_in:
            if line.startswith('ATOM'):
                atom_num = int(line[6:11])
                if atom_num in atom_name_map:
                    # Replace atom name and residue info while maintaining PDB format
                    orig_name = atom_name_map[atom_num]
                    res_name = residue_name_map[atom_num]
                    res_num = residue_num_map[atom_num]
                    # PDB format: maintain proper spacing
                    new_line = (f"{line[:12]}{orig_name:<5}{res_name:<8}{res_num:<7}{line[32:]}")
                    f_out.write(new_line)
                else:
                    f_out.write(line)
            else:
                f_out.write(line)
   
    # Replace the original PDB with fixed version
    shutil.move(fixed_pdb, pdb_file)

    # Fix the TOP file atom names and residues
    fixed_top = f'{inputPrefix}_gmx_fixed.top'
    
    with open(top_file, 'r') as f_in, open(fixed_top, 'w') as f_out:
        in_atoms_section = False
        for line in f_in:
            if line.startswith('[ atoms ]'):
                in_atoms_section = True
                f_out.write(line)
            elif in_atoms_section and line.strip() and not line.startswith(';') and not line.startswith('['):
                # Parse atom line in topology file
                parts = line.split()
                atom_num = int(parts[0])
                if atom_num in atom_name_map:
                    # Replace atom name and residue info while maintaining alignment
                    parts[4] = atom_name_map[atom_num]  # atom name
                    parts[3] = residue_name_map[atom_num]  # residue name
                    parts[2] = str(residue_num_map[atom_num])  # residue number
                    f_out.write(f"{parts[0]:>6} {parts[1]:>6} {parts[2]:>6} {parts[3]:>6} {parts[4]:<4} {' '.join(parts[5:])}\n")
                else:
                    f_out.write(line)
            elif line.strip().startswith('['):
                in_atoms_section = False
                f_out.write(line)
            else:
                f_out.write(line)
    
    # Replace the original TOP with fixed version
    shutil.move(fixed_top, top_file)
    
    
    # Convert PDB to GRO using GROMACS
    cmd = Command(
        f'{sw.gmx} editconf -f {pdb_file} -o {outputPrefix}.gro',
    )
    cmd.run()
    
    logger.debug(f'GROMACS generated {outputPrefix}.gro')
    # Convert topology to GROMACS format
    groOut = f'{outputPrefix}.gro'
    topOut = f'{outputPrefix}.top'
    itpOut = f'{outputPrefix}.itp'
    Command(f'cp {top_file} {topOut}').run()
    
    # save the results of cgenff_to_gmx.sh as Gromacs gro and top files
    try:
        file = GromacsTopologyFile(top_file)
        logger.debug(f'Writing {topOut}')
        file.save(topOut, parameters=itpOut, overwrite=True)
    except Exception as m:
        logger.error('Unspecified parmed error')
        raise parmed.exceptions.GromacsError(m)
    
    logger.debug(f"good generation of {outputPrefix} parameters")

