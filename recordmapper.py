import pandas as pd
import math
import statistics
from random import *

class RecordMapper:

    @staticmethod
    def fill_rxn_order(output_hdf5, output_reaction_df, output_record_df):
        """
        Augments the given reaction dataframe with the reaction order
        :param output_hdf5: Output HDF5 file path
        :param output_reaction_df: Reaction df key
        :param output_record_df: Record df key
        :return: None
        """

        # Reading dataframes from the HDF5 file
        data_store = pd.HDFStore(output_hdf5)  # Opening the HDF5 file
        reaction_dataframe = data_store[output_reaction_df]  # Reading the dataframe
        record_dataframe = data_store[output_record_df]  # Reading the dataframe
        data_store.close()

        # Creating an list of lists: [list1, list2, ...., list-n]
        # Index in parent list correspond to reaction ID
        # Child list correspond to reaction order occurences
        parent_list = [[] for ele in range((max(reaction_dataframe.index)+1))]  # +1 because we want to include
        # max(reaction_dataframe.index) also
        for _, row in record_dataframe.iterrows():
            if not math.isnan(row['ReactionOrder']):
                parent_list[row['RID']].append(row['ReactionOrder'])

        # Extracting Modes of reaction orders
        for idx in range(max(reaction_dataframe.index)+1):
            try:
                if len(parent_list[idx]) == 0:
                    parent_list[idx] = math.nan
                else:
                    parent_list[idx] = statistics.mode(parent_list[idx])
            except Exception:
                if len(parent_list[idx]) > 0:
                    rand_index = randint(0, len(parent_list[idx])-1)
                    parent_list[idx] = parent_list[idx][rand_index]
                else:
                    print('Something wrong with Reaction ID {}'.format(idx))

        # Append Reaction Order Column to Reaction Dataframe
        rxn_order = [""]*len(reaction_dataframe.index)
        reaction_dataframe['ReactionOrder'] = rxn_order
        for idx, rows in reaction_dataframe.iterrows():
            reaction_dataframe.at[idx, 'ReactionOrder'] = parent_list[idx]

        # Updating dataframe in the HDF5 file
        reaction_dataframe.to_hdf(path_or_buf=output_hdf5, key=output_reaction_df, mode='a')

        print("--Records stored into HDF5 file--")

    @staticmethod
    def fill_activ_enrgy(output_hdf5, output_reaction_df, output_record_df):
        """
        Augments the given reaction dataframe with the reaction order
        :param output_hdf5: Output HDF5 file path
        :param output_reaction_df: Reaction df key
        :param output_record_df: Record df key
        :return: None
        """

        # Reading dataframes from the HDF5 file
        data_store = pd.HDFStore(output_hdf5)  # Opening the HDF5 file
        reaction_dataframe = data_store[output_reaction_df]  # Reading the dataframe
        record_dataframe = data_store[output_record_df]  # Reading the dataframe
        data_store.close()

        # Creating an list of lists: [list1, list2, ...., list-n]
        # Index in parent list correspond to reaction ID
        # Child list correspond to reaction order occurences
        parent_list = [[] for ele in range((max(reaction_dataframe.index) + 1))]  # +1 because we want to include
        # max(reaction_dataframe.index) also
        for _, row in record_dataframe.iterrows():
            if not math.isnan(row['ActivationEnergy']):
                parent_list[row['RID']].append(row['ActivationEnergy'])

        # Extracting Modes of reaction orders
        for idx in range(max(reaction_dataframe.index) + 1):
            if len(parent_list[idx]) == 0:
                parent_list[idx] = math.nan
            else:
                parent_list[idx] = statistics.mean(parent_list[idx])

        # Append Reaction Order Column to Reaction Dataframe
        act_enrg = [""] * len(reaction_dataframe.index)
        reaction_dataframe['ActivationEnergy'] = act_enrg
        for idx, rows in reaction_dataframe.iterrows():
            reaction_dataframe.at[idx, 'ActivationEnergy'] = parent_list[idx]

        # Updating dataframe in the HDF5 file
        reaction_dataframe.to_hdf(path_or_buf=output_hdf5, key=output_reaction_df, mode='a')

        print("--Records stored into HDF5 file--")

    @staticmethod
    def map_rid_to_cid(output_hdf5, output_reaction_df, input_species_df):
        """
        Map SID list to CID list
        :param output_hdf5: HDF5 data store path
        :param output_reaction_df: Reaction dataframe key
        :param input_species_df: Species dataframe key
        :return: None
        """

        # Reading dataframes from the HDF5 file
        data_store = pd.HDFStore(output_hdf5)  # Opening the HDF5 file
        reaction_dataframe = data_store[output_reaction_df]  # Reading the dataframe
        species_dataframe = data_store[input_species_df]  # Reading the dataframe
        data_store.close()

        # Creating a dict key:value ==> sid:cid
        sid_to_cid = {}
        for sid, species_row in species_dataframe.iterrows():
            if species_row['CID'] == "" or math.isnan(species_row['CID']):  # If empty or nan
                sid_to_cid[sid] = math.nan
            else:
                sid_to_cid[sid] = int(species_row['CID'])

        # Appending reactant and product cid list column
        reaction_dataframe['ReactantCID'] = [[] for _ in reaction_dataframe.index]
        reaction_dataframe['ProductCID'] = [[] for _ in reaction_dataframe.index]
        for idx, rxn_row in reaction_dataframe.iterrows():
            reaction_dataframe.at[idx, 'ReactantCID'] = [sid_to_cid[ele] for ele in rxn_row['Reactants_SIDs_List']]
            reaction_dataframe.at[idx, 'ProductCID'] = [sid_to_cid[ele] for ele in rxn_row['Products_SIDs_List']]

        # Store to HDF Data store
        reaction_dataframe.to_hdf(path_or_buf=output_hdf5, key=output_reaction_df, mode='a')

        print("--Records stored into HDF5 file--")


# Code Run Check
# RecordMapper.fill_rxn_order('NewGen2Output/NewGen.h5', 'Reactions', 'Records')
# RecordMapper.fill_activ_enrgy('NewGen2Output/NewGen.h5', 'Reactions', 'Records')
# RecordMapper.map_rid_to_cid('PreliminaryOutput/DemoGenerated/DataDF.h5', 'Reactions', 'Species')
