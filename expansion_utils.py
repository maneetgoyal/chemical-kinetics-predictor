# Functions that help expand the current databse when new PubChem IDs have been fetched
import pandas as pd
import re
import math

class ExpansionUtils:

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
                if math.isnan(ele) or re.search("adduct", row['Reactants']+row['Products'], re.IGNORECASE) is not None:
                    to_be_dropped.append(idx)
                    break

        rxn_df = rxn_df.drop(to_be_dropped)
        return rxn_df

    @staticmethod
    def inject_new_data_to_species(new_species_xlsx, output_hdf5, species_df_key):
        """
        Helps to inject new data to the species dataframe as more CIDs are fetched manually
        :param new_species_xlsx: New Species xlsx file which stores newly fetched CIDs
        :param output_hdf5: Output HDF% file that houses species df
        :param species_df_key: Specied df key in output_hdf5 file
        :return: None
        """
        new_species_df = pd.read_excel(new_species_xlsx, index_col=0, header=0)
        hdf5_fp = pd.HDFStore(output_hdf5)
        old_species_df = hdf5_fp[species_df_key]
        old_species_df['CID'] = new_species_df['CID']
        hdf5_fp[species_df_key] = old_species_df
        hdf5_fp.close()

        # Pushing new species data into df
        # Populator.get_pubchem_data(output_hdf5, species_df_key)
        return None

# Code run check
# ExpansionUtils.get_rxn_subset('PreliminaryOutput/DemoGenerated/DataDF.h5', 'Reactions')
