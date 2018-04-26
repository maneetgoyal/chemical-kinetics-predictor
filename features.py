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

    def bond_brk(self, reactant_cid_list, reactant_vec_list, product_cid_list, product_vec_list):
        """
        Calculates the feature vector of reactions.
        :param reactant_cid_list: List of CIDs of the reactants
        :param reactant_vec_list: List of feature vectors of the reactant species
        :param product_cid_list: List of CIDs of the products
        :param product_vec_list: list of Feature Vectors of the product species
        :return: Reaction Feature Vector
        """

        prod_sum = np.zeros(len(self.bond_vec.index))
        reac_sum = np.zeros(len(self.bond_vec.index))

        for idx, ele in enumerate(product_cid_list):
            if str(ele) == "-1":
                prod_sum = prod_sum + FeatureConstructor.bonds_count_json(self, -1, None)
            else:
                prod_sum = prod_sum + FeatureConstructor.bonds_count_json(self, None, product_vec_list[idx])

        for idx, ele in enumerate(reactant_cid_list):
            if str(ele) == "-1":
                reac_sum = reac_sum + FeatureConstructor.bonds_count_json(self, -1, None)
            else:
                reac_sum = reac_sum + FeatureConstructor.bonds_count_json(self, None, reactant_vec_list[idx])
        #
        # bond_change = \
        #     sum([FeatureConstructor.bonds_count_json(self, prod_cid, product_str_json_list) for prod_cid in product_cid_list]) \
        #     - sum([FeatureConstructor.bonds_count_json(self, reac_cid, reactant_str_json_list) for reac_cid in reactant_cid_list])

        bond_change = prod_sum - reac_sum

        return list(bond_change)

    def get_species_subset(self, output_hdf, species_df_key, subset_species_df_key):
        """
        Creates a subset dataframe of the species whose feature vector has been identified
        :param output_hdf: HDF5 file path
        :param species_df_key: Species Dataframe key
        :param subset_species_df_key: Subset-Species Dataframe Key
        :return: None
        """
        species_df = pd.read_hdf(output_hdf, species_df_key)
        subset_species_df = species_df.query('CID >= -1')
        subset_species_df['CID'] = subset_species_df['CID'].astype(int)
        subset_species_df = subset_species_df.reset_index()
        subset_species_df = subset_species_df.set_index(keys='CID', drop='False', verify_integrity=True)

        # Adding feature vectors to the species subset dataframe
        subset_species_df['FeatureVector'] = [""] * len(subset_species_df.index)

        # Appending feature vectors
        for idx, row in subset_species_df.iterrows():
            if idx == -1:
                subset_species_df.at[idx, 'FeatureVector'] = FeatureConstructor.bonds_count_json(self, -1, None)
            else:
                subset_species_df.at[idx, 'FeatureVector'] = FeatureConstructor.bonds_count_json(self, None, row['BondsInfo'])

        # Storing the subset dataframe
        subset_species_df = subset_species_df.reset_index()
        subset_species_df = subset_species_df.set_index(keys='SID', drop='False', verify_integrity=True)
        subset_species_df.to_hdf(output_hdf, subset_species_df_key)


# Code run check
# my_constructor = FeatureConstructor("FeatureLibrary/elements.csv", "FeatureLibrary/bonds.csv")
# print(my_constructor.bonds_count_json(1000))
