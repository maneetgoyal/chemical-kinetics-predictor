# Helps contruct feature vectors for indivdual reactions
import requests
import json
import numpy as np
import pandas as pd

class FeatureConstructor:

    def __init__(self, elements_csv, bonds_csv):
        """
        Initializer. Reads the csv which contains all the chemical elements
        :param elements_csv: Input CSV which contains all the chemical elements
        """
        self.elem = pd.read_csv(elements_csv, header=0)  # Datafame storing chemical elements
        self.periodic = [''] + list(self.elem.values)  # List Containing all the chemical elements
        self.bond_type = [' ', '-', '=', '#']  # List containing Bond types
        self.bond_vec = pd.read_csv(bonds_csv, header=0).values  # NP array Containing types of bonds.
        # It is the Feature Blueprint for the chemicals.

    @staticmethod
    def get_full(cid):
        """
        Gives the information on a chemical via its CID in a JSON format.
        :param cid: PubChem ID of a Chemical
        :return: Stringified JSON of Chemical's Data
        """
        r = requests.get('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/record/json'.format(cid))
        return r.text

    def json_str(self, stringified_json):
        """
        Get chemical bonds from a JSON format.
        :param stringified_json: Strigified JSON containing bonds info,
        fetched from Pubchem using Pubchem ID of chemical
        :return: a list containing all the bonds in text form, as in, C-H, S-H, etc.
        """

        # Loads the stringified data in a JSON format
        prop_json = json.loads(stringified_json)

        # Stacking up Atoms index along side their ID (Atomic Number)
        atoms_index = prop_json['PC_Compounds'][0]['atoms']['aid']
        atoms_identity = prop_json['PC_Compounds'][0]['atoms']['element']
        atoms = np.vstack((atoms_index, atoms_identity))

        # This 2-D array gives us a list of atoms and their index
        bonds = []

        # Obtain a list of bonds available from their index of position, element types and bond types
        if 'bonds' in prop_json['PC_Compounds'][0]:

            for i in range(len(prop_json['PC_Compounds'][0]['bonds']['order'])):

                atom_index_1 = prop_json['PC_Compounds'][0]['bonds']['aid1'][i]
                atom1 = atoms[1, atom_index_1 - 1]
                atom_index_2 = prop_json['PC_Compounds'][0]['bonds']['aid2'][i]
                atom2 = atoms[1, atom_index_2 - 1]

                # we will create a combination of atoms and bond to in a format like "C-C"
                bond_str = self.periodic[atom1] + \
                           self.bond_type[prop_json['PC_Compounds'][0]['bonds']['order'][i]] + \
                           self.periodic[atom2]

                bonds.append(bond_str)

        if 'radical' in prop_json['PC_Compounds'][0]['atoms']:

            radical_atom_index = prop_json['PC_Compounds'][0]['atoms']['radical'][0]['aid']
            radical_atom = atoms[1, radical_atom_index - 1]
            radical_str = self.periodic[radical_atom] + '.'
            bonds.extend([radical_str])

        return bonds

    def bonds_count_json(self, cid, stringified_json):
        """
        Creates feature vectors by comparing bonds of the chemical (given by cid) with our bond vec
        :param cid: Pubchem CID
        :param stringified_json: Strigified JSON containing bonds info,
        fetched from Pubchem using Pubchem ID of chemical
        :return: Feature vector of the input chemical
        """
        if cid == -1:  # Artificial CID for .CH Radical
            feature_vec = np.zeros(len(self.bond_vec))
            feature_vec[21] = 1
            feature_vec[35] = 1
        else:
            bonds = FeatureConstructor.json_str(self, stringified_json)
            feature_vec = np.zeros(len(self.bond_vec))
            bond_vec = np.array(self.bond_vec)

            # we construct a feature vector containing counting on different bonds
            for i, bond in enumerate(bonds):
                if np.any(bond_vec == bond):
                    idx = np.where(bond_vec == bond)
                    feature_vec[idx[0]] = feature_vec[idx[0]] + 1
        return feature_vec

    def bond_brk(self, reactant_cid, reactant_str_json, product_cid, product_str_json):

        bond_change = \
            sum([FeatureConstructor.bonds_count_json(self, prod_cid, product_str_json) for prod_cid in product_cid]) \
            - sum([FeatureConstructor.bonds_count_json(self, reac_cid, reactant_str_json) for reac_cid in reactant_cid])

        return list(bond_change)


# Code run check
# my_constructor = FeatureConstructor("FeatureLibrary/elements.csv", "FeatureLibrary/bonds.csv")
# print(my_constructor.bonds_count_json(1000))
