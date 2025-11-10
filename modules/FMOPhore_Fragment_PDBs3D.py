import os
import glob
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import math
import argparse
from rdkit.Chem import Recap, AllChem

class Fragmentation3Dlig:
    def __init__(self, pdb_file, distance_cutoff, same_target=False, Binding_Energy=False):
        self.pdb_file = pdb_file
        self.same_target = same_target
        self.Binding_Energy = Binding_Energy
        self.chain_id_list = []
        self.ligands_list = []

    def Fragmentation3Dlig_processor(self):
        pdb_lines = self.load_pdb_file()
        pdb_lines = self.renumber_atoms(pdb_lines)
        pdb_lines = self.recap_fragmentation_and_save(pdb_lines)
        pdb_lines = self.process_multiple_ligands(pdb_lines)
        self.save_pdb_file(pdb_lines)

    ######################################################################################################################################
        if args.PDB:
            pdb_path = args.PDB
            lig_files = []
            if os.path.isdir(pdb_path):
                lig_files = glob.glob(os.path.join(pdb_path, "*lig*.pdb"))
            elif os.path.isfile(pdb_path) and pdb_path.endswith(".pdb"):
                lig_files.append(pdb_path)
            else:
                raise ValueError("Invalid PDB file or directory path.")
            if not lig_files:
                raise ValueError("No ligand files found in the specified path.")
            for pdb_file in lig_files:
                renumber_atoms(pdb_file)
                fragment_files, _ = recap_fragmentation_and_save(pdb_file)
                if fragment_files:
                    process_multiple_ligands(fragment_files)
    ################################################################
    def load_pdb_file(self):
        with open(self.pdb_file, 'r') as f:
            return f.readlines()
    ###################################################################            
    def renumber_atoms(pdb_file):
        updated_lines = []
        atom_counters = {}
        with open(pdb_file, 'r') as infile:
            for line in infile:
                if line.startswith("HETATM"):
                    residue_name = line[17:20].strip()
                    residue_number = line[22:26].strip()
                    chain_id = line[21].strip()
                    atom_name = line[12:16].strip()
                    residue_key = (residue_name, residue_number, chain_id)
                    if residue_key not in atom_counters:
                        atom_counters[residue_key] = {}
                    if atom_name.startswith("HX"):
                        base_name = "H"
                    elif atom_name[0].isdigit():
                        base_name = ''.join([char for char in atom_name[1:] if not char.isdigit()])
                    else:
                        base_name = ''.join([char for char in atom_name if not char.isdigit()])
                    if base_name not in atom_counters[residue_key]:
                        atom_counters[residue_key][base_name] = 0
                    atom_counters[residue_key][base_name] += 1
                    new_atom_name = f"{base_name}{atom_counters[residue_key][base_name]}"
                    new_atom_name = new_atom_name[:4].ljust(4)
                    updated_line = line[:12] + new_atom_name + line[16:]
                    updated_lines.append(updated_line)
                else:
                    updated_lines.append(line)
        with open(pdb_file, 'w') as outfile:
            outfile.writelines(updated_lines)
        print(f"Renamed and renumbered PDB file saved as: {pdb_file}")
    ###################################################################    
    def recap_fragmentation_and_save(pdb_file):
        try:
            mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
            if mol is None:
                raise ValueError(f"Could not parse PDB file: {pdb_file}")
            recap_tree = Recap.RecapDecompose(mol)
            if recap_tree is None or not recap_tree.children:
                raise ValueError("No RECAP fragmentation possible for this molecule.")
            fragment_files = []
            fragment_coords = []
            idx = 0
            for submol_node in recap_tree.GetLeaves().values():
                submol = submol_node.mol 
                if submol is None:
                    continue
                non_hydrogen_atoms = [atom for atom in submol.GetAtoms() if atom.GetSymbol() != "H"]
                if not non_hydrogen_atoms:
                    print(f"Fragment {idx} contains only hydrogen atoms. Skipping...")
                    continue
                idx += 1
                pdb_block = Chem.MolToPDBBlock(submol)
                filtered_pdb_lines = [
                    line for line in pdb_block.splitlines() if '*' not in line]
                if not any(
                    line.startswith("HETATM") and line[76:78].strip() != "H"
                    for line in filtered_pdb_lines):
                    print(f"Fragment {idx} contains only hydrogen atoms after filtering. Skipping...")
                    continue
                frag_file = f"{pdb_file[:-4]}_frag{idx}.pdb"
                with open(frag_file, 'w') as f:
                    f.write('\n'.join(filtered_pdb_lines))
                fragment_files.append(frag_file)
                conf = mol.GetConformer()
                frag_coords = [
                    (atom.GetIdx(), conf.GetAtomPosition(atom.GetIdx()))
                    for atom in submol.GetAtoms()
                ]
                fragment_coords.append(frag_coords)
            return fragment_files, fragment_coords
        except Exception as e:
            print(f"Error during Recap fragmentation: {e}")
            return None, None
    ###################################################################        
    def calculate_distance(coord1, coord2):
        return math.sqrt(sum((c1 - c2) ** 2 for c1, c2 in zip(coord1, coord2)))
    ###################################################################    
    def find_closest_atoms(coords1, coords2):
        min_distance = float('inf')
        closest_pair = None
        for idx1, coord1 in enumerate(coords1):
            for idx2, coord2 in enumerate(coords2):
                distance = calculate_distance(coord1[1], coord2[1])  
                if distance < min_distance:
                    min_distance = distance
                    closest_pair = (idx1, idx2)
        return closest_pair
    ###################################################################    
    def add_dummy_atom_or_replace(mol, conf, closest_idx, target_coord, closest_symbol):
        editable_mol = Chem.EditableMol(mol)
        closest_atom = mol.GetAtomWithIdx(closest_idx)
        residue_info = closest_atom.GetPDBResidueInfo()
        residue_name = residue_info.GetResidueName() if residue_info else "UNK"
        residue_number = residue_info.GetResidueNumber() if residue_info else 1
        chain_id = residue_info.GetChainId() if residue_info else "A"
        existing_dummy_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "*")
        dummy_suffix = existing_dummy_count + 1 

        if closest_symbol == "H":
            editable_mol.ReplaceAtom(closest_idx, Chem.Atom("*"))
            mol_updated = editable_mol.GetMol()
            conf_updated = mol_updated.GetConformer()
            coord = np.array([conf.GetAtomPosition(closest_idx).x,
                              conf.GetAtomPosition(closest_idx).y,
                              conf.GetAtomPosition(closest_idx).z])
            updated_atom = mol_updated.GetAtomWithIdx(closest_idx)
            pdb_info = Chem.AtomPDBResidueInfo()
            pdb_info.SetResidueName(residue_name)
            pdb_info.SetResidueNumber(residue_number)
            pdb_info.SetChainId(chain_id)
            pdb_info.SetName(f" *{dummy_suffix} ")
            pdb_info.SetIsHeteroAtom(True)
            updated_atom.SetPDBResidueInfo(pdb_info)
            return mol_updated, f"Replaced H with *{dummy_suffix} at: ({coord[0]:.3f}, {coord[1]:.3f}, {coord[2]:.3f})"
        else:
            atom_idx = editable_mol.AddAtom(Chem.Atom("*"))
            editable_mol.AddBond(closest_idx, atom_idx, Chem.BondType.SINGLE)
            mol_updated = editable_mol.GetMol()
            conf_updated = mol_updated.GetConformer()
            closest_coord = np.array([conf.GetAtomPosition(closest_idx).x,
                                       conf.GetAtomPosition(closest_idx).y,
                                       conf.GetAtomPosition(closest_idx).z])
            target_vector = target_coord - closest_coord
            normalized_vector = target_vector / np.linalg.norm(target_vector)
            arbitrary_vector = np.array([1.0, 0.0, 0.0])
            if np.allclose(normalized_vector, arbitrary_vector):  # If vectors are parallel
                arbitrary_vector = np.array([0.0, 1.0, 0.0])
            perpendicular_vector = np.cross(normalized_vector, arbitrary_vector)
            perpendicular_vector /= np.linalg.norm(perpendicular_vector)
            bond_length = 1.5
            new_coord = closest_coord + normalized_vector * bond_length + perpendicular_vector * 0.5
            conf_updated.SetAtomPosition(atom_idx, new_coord.tolist())
            new_atom = mol_updated.GetAtomWithIdx(atom_idx)
            pdb_info = Chem.AtomPDBResidueInfo()
            pdb_info.SetResidueName(residue_name)
            pdb_info.SetResidueNumber(residue_number)
            pdb_info.SetChainId(chain_id)
            pdb_info.SetName(f" *{dummy_suffix} ")  
            pdb_info.SetIsHeteroAtom(True)
            new_atom.SetPDBResidueInfo(pdb_info)
            return mol_updated, f"*{dummy_suffix} added at: ({new_coord[0]:.3f}, {new_coord[1]:.3f}, {new_coord[2]:.3f})"
    ###################################################################        
    def contains_nan_coordinates(pdb_lines):
        for line in pdb_lines:
            if line.startswith("HETATM") and (
                "nan" in line[30:38] or "nan" in line[38:46] or "nan" in line[46:54]):
                return True
        return False
    ###################################################################
    def process_multiple_ligands(pdb_files):
        ligand_data = []
        for pdb_file in pdb_files:
            mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
            if mol is None:
                print(f"Failed to load PDB file: {pdb_file}")
                continue
            try:
                Chem.SanitizeMol(mol)
            except Chem.SanitizeException as e:
                print(f"Sanitization failed for {pdb_file}: {e}")
                continue
            try:
                pdb_lines1 = Chem.MolToPDBBlock(mol, kekulize=False).splitlines()
            except Exception as e:
                print(f"Failed to generate PDB block for {pdb_file}: {e}")
                continue
            conf = mol.GetConformer()
            coords = [(atom.GetSymbol(), np.array([conf.GetAtomPosition(atom.GetIdx()).x, 
                                                   conf.GetAtomPosition(atom.GetIdx()).y, 
                                                   conf.GetAtomPosition(atom.GetIdx()).z])) 
                      for atom in mol.GetAtoms()]
            ligand_data.append((pdb_file, mol, coords))
        for i, (pdb_file1, mol1, coords1) in enumerate(ligand_data):
            for j, (pdb_file2, mol2, coords2) in enumerate(ligand_data):
                if i >= j:
                    continue
                closest_idx1, closest_idx2 = find_closest_atoms(coords1, coords2)
                closest_atom1 = coords1[closest_idx1]
                closest_atom2 = coords2[closest_idx2]
                mol1, dummy1 = add_dummy_atom_or_replace(mol1, mol1.GetConformer(), closest_idx1, closest_atom2[1], closest_atom1[0])
                mol2, dummy2 = add_dummy_atom_or_replace(mol2, mol2.GetConformer(), closest_idx2, closest_atom1[1], closest_atom2[0])
                output_file1 = os.path.join(f"{os.path.basename(pdb_file1).replace('.pdb', '_updated.pdb')}")
                output_file2 = os.path.join(f"{os.path.basename(pdb_file2).replace('.pdb', '_updated.pdb')}")
                pdb_lines1 = Chem.MolToPDBBlock(mol1).splitlines()
                pdb_lines2 = Chem.MolToPDBBlock(mol2).splitlines()
                if contains_nan_coordinates(pdb_lines1):
                    print(f"Skipping {output_file1} due to NaN coordinates.")
                else:
                    Chem.MolToPDBFile(mol1, output_file1)
                if contains_nan_coordinates(pdb_lines2):
                    print(f"Skipping {output_file2} due to NaN coordinates.")
                else:
                    Chem.MolToPDBFile(mol2, output_file2)
###################################################################################################
if __name__ == "__main__":
    """
    FMOPhore V.0.1 - Frag_processor - Copyright "©" 2024, Peter E.G.F. Ibrahim.
    """
    parser = argparse.ArgumentParser(description="Process PDB files with ligands")
    parser.add_argument('-pdb', '--PDB', type=str, help="Directory containing PDB files to be fragmented into pieces")
    args = parser.parse_args()
    Frag_processor = Fragmentation3Dlig(args.pdb_file)
    Frag_processor.Fragmentation3Dlig_processor()
