# Functions that help expand the current databse when new PubChem IDs have been fetched
import pandas as pd
import re
import math
import features as ft
import recordmapper as rm

class Extender:

    @staticmethod
    def export_all_to_excel(input_hdf5, out_directory_path):
        """
        Exports all the dataframe in HDF5 file to .xlsx format with their Keys as the output file names
        :param input_hdf5: As title
        :param out_directory_path: Path of the directory/folder where the exported xlsx are to be stored
        :return: None
        """
        data_store = pd.HDFStore(input_hdf5)  # Opening the HDF5 file
        for each_key in data_store.keys():
            data_store[each_key].to_excel(out_directory_path + each_key + ".xlsx")
            # '/' missing between folder name and
            # file name because file name already includes it.
        data_store.close()

        print("-- Dataframes written to Excel files (.xlsx) --")

    @staticmethod
    def get_rxn_subset(input_hdf5, rxn_df_key):
        """
        Fetch the subset of reactions that can be used for machine learning experiments
        :param input_hdf5: As title
        :param rxn_df_key: As title
        :return: None
        """
        rxn_df = pd.read_hdf(input_hdf5, rxn_df_key)
        # Subset filtering criteria; can be updated as required
        # For instance, when the Pubchem ids of all the chemicals with score >= 50 have been fetched,
        # Status75 below can be chnaged to Status 50. Currently, all the chemicals with scores >= 90
        # have been fetched
        rxn_df = rxn_df.query("Products_Available == True & Status_75 == True")

        # Now we will get rid of all the rows that have 'NaN' in their CID list or if they have any adduct
        to_be_dropped = []
        for idx, row in rxn_df.iterrows():
            cid_list = row['ReactantCID'] + row['ProductCID']
            for ele in cid_list:
                if math.isnan(ele) or re.search("adduct", row['Reactants']+row['Products'], re.IGNORECASE):
                    to_be_dropped.append(idx)
                    break

        rxn_df = rxn_df.drop(to_be_dropped)
        return rxn_df

    @staticmethod
    def expand_data(new_species_xlsx, output_hdf5, species_df_key, rxn_df_key, elements_csv, bonds_csv, new_xlsx_path):
        """
        Helps to inject new data to the species dataframe as more CIDs are fetched manually
        :param new_species_xlsx: New Species xlsx file which stores newly fetched CIDs
        :param output_hdf5: Output HDF% file that houses species df
        :param species_df_key: Specied df key in output_hdf5 file
        :param rxn_df_key:
        :param elements_csv:
        :param bonds_csv:
        :param new_xlsx_path:
        :return: New data for ML experiments
        """

        # Reading xlsx files which contains newly fetched PubChem IDs into pandas df
        new_df_from_xlsx = pd.read_excel(new_species_xlsx, header=0)

        # Reading old Species dataframe to which new PubChem ids have to be transfered
        old_df_from_hdf = pd.read_hdf(output_hdf5, species_df_key)

        # Setting 'Species' name as index for efficiency
        old_df_from_hdf = old_df_from_hdf.reset_index()
        old_df_from_hdf = old_df_from_hdf.set_index(keys="Species", verify_integrity=True)

        # Initializing FeatureConstructor
        my_constructor = ft.FeatureConstructor(elements_csv, bonds_csv)

        # Transfering CID, adding BondsInfo (stringified PubChem JSON), adding species feature vector
        new_species_count = 0
        for idx, row in new_df_from_xlsx.iterrows():
            if not math.isnan(row['CID']) and row['CID'] != "":
                if math.isnan(old_df_from_hdf.at[row['Species'], 'CID']) or old_df_from_hdf.at[row['Species'], 'CID'] == "":
                    old_df_from_hdf.at[row['Species'], 'CID'] = row['CID']
                    pubchem_str_json = my_constructor.get_full(row['CID'])
                    print("--Data fetched for CID {}--".format(int(row['CID'])))
                    old_df_from_hdf.at[row['Species'], 'BondsInfo'] = pubchem_str_json
                    old_df_from_hdf.at[row['Species'], 'FeatureVector'] = my_constructor.bonds_count_json(None, pubchem_str_json)
                    new_species_count = new_species_count + 1

        print('--Status--')
        print('--{} New Species Added--'.format(new_species_count))

        if new_species_count == 0:

            print('No new changes were made as there were no new species to add.')
            return

        else:

            # Updating HDF with updated species df
            old_df_from_hdf = old_df_from_hdf.reset_index()
            old_df_from_hdf = old_df_from_hdf.set_index(keys="SID", verify_integrity=True)
            old_df_from_hdf.to_hdf(output_hdf5, species_df_key)

            # Updating Reactions DF with new CID list
            rm.RecordMapper.map_rid_to_cid(output_hdf5, rxn_df_key, species_df_key)

            # Filetring out reactions whose feature vectors can be calculated
            reduced_rxn_df = Extender.get_rxn_subset(output_hdf5, rxn_df_key)

            # Creating feature vectors of the filtered out reactions
            reduced_rxn_df = my_constructor.bond_brk(output_hdf5, species_df_key, reduced_rxn_df)

            print('--Status--')
            print('--Reactions Feature Vectors Created--')

            # Creating the new reactions xlsx for ML Training
            reduced_rxn_df.to_excel(new_xlsx_path)

            print('--Status--')
            print('--Database Expansion Routine Complete--')


# Code run check
# ExpansionUtils.get_rxn_subset('PreliminaryOutput/DemoGenerated/DataDF.h5', 'Reactions')
