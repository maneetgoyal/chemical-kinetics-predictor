# Helps contruct feature vectors for indivdual reactions
import requests
import json
import numpy as np
import pandas as pd
import math

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
        feature_vec = np.zeros(len(self.bond_vec))
        if cid == -1:  # Artificial CID for .CH Radical
            feature_vec[21] = 1
            feature_vec[35] = 1
        elif cid == -2:  # Handling M (Catalyst)
            pass
        elif cid in [-3, -4]:  # Handling missing products/reactants cases
            return None
        else:  # Handling valid CIDs
            bonds = FeatureConstructor.json_str(self, stringified_json)
            bond_vec = np.array(self.bond_vec)
            # we construct a feature vector containing counting on different bonds
            for i, bond in enumerate(bonds):
                if np.any(bond_vec == bond):
                    idx = np.where(bond_vec == bond)
                    feature_vec[idx[0]] = feature_vec[idx[0]] + 1

        return feature_vec

    def bond_brk(self, input_hdf, species_df_key, input_rxn_df):
        """
        Calculates the feature vectors of all the reactions in the given dataframe
        :param input_hdf: As title
        :param species_df_key: Species Dataframe key
        :param input_rxn_df: Reaction Dataframe
        :return: Reaction Dataframe with feature vectors
        """

        # Species dataframe; it's indexed at SID
        species_df = pd.read_hdf(input_hdf, species_df_key)

        # Start populating feature vectors
        feat_vec = [""]*len(input_rxn_df.index)
        input_rxn_df['FeatureVector'] = feat_vec
        for idx, row in input_rxn_df.iterrows():

            # Initializing feature vector components
            prod_sum = np.zeros(len(self.bond_vec))
            reac_sum = np.zeros(len(self.bond_vec))

            # SID lists
            product_sid_list = row['Reactants_SIDs_List']
            reactant_sid_list = row['Products_SIDs_List']

            # Getting feature of products
            for ele in product_sid_list:
                prod_sum = prod_sum + species_df.at[ele, 'FeatureVector']

            # Getting feature vector of reactants
            for ele in reactant_sid_list:
                reac_sum = reac_sum + species_df.at[ele, 'FeatureVector']

            # Getting feature vector of reaction
            bond_change = prod_sum - reac_sum

            # Updating the reaction data frame
            input_rxn_df.at[idx, 'FeatureVector'] = bond_change

        return input_rxn_df

    def create_species_feat_vec(self, output_hdf, species_df_key):
        """
        Creates a subset dataframe of the species whose feature vector has been identified
        :param output_hdf: HDF5 file path
        :param species_df_key: Species Dataframe key
        :return: None
        """
        species_df = pd.read_hdf(output_hdf, species_df_key)

        # Adding feature vectors to the species dataframe
        species_df['FeatureVector'] = [""] * len(species_df.index)

        # Appending feature vectors
        for idx, row in species_df.iterrows():
            if row['CID'] != "" and not math.isnan(row['CID']):
                if row['CID'] <= 0:  # Handling artificial/invalid CIDs
                    artficial_bond_count = FeatureConstructor.bonds_count_json(self, row['CID'], None)
                    if artficial_bond_count is not None:
                        species_df.at[idx, 'FeatureVector'] = artficial_bond_count
                else:  # Handling valid CIDs
                    species_df.at[idx, 'FeatureVector'] = FeatureConstructor.bonds_count_json(self, None, row['BondsInfo'])

        # Updating the species dataframe in the HDF file
        species_df.to_hdf(output_hdf, species_df_key)


# Code run check
# my_constructor = FeatureConstructor("FeatureLibrary/elements.csv", "FeatureLibrary/bonds.csv")
# print(my_constructor.bonds_count_json(1000))

my_constructor = FeatureConstructor("FeatureLibrary/elements.csv", "FeatureLibrary/bonds.csv")
my_constructor.create_species_feat_vec('PreliminaryOutput/DemoGenerated/DataDF.h5', 'Species')
