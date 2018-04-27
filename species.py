import pandas as pd
from chemspipy import ChemSpider
import logging
import requests
import math
logging.basicConfig(level=logging.DEBUG)

class Populator:
    """
    Populates the dataframe with ChemSpider/PubChem Results
    """
    def __init__(self):
        """
        Initializes all the object variables
        """

        # Reaction Dataframe
        self.reactions_dataframe = None

        # Reactant Dataframe
        self.species_df = None

        # Unique Reactants Dictionary
        self.unique_species_dict = None

        # Creating a transator for cleaning individual reactants off non-familiar characters
        self.translator = str.maketrans("Î", "α", "Â±€™")  # Argument style
        # (# intab,outtab,character string that should be mapped to None)

        # Autheticating ChemSpider API using the token
        self.security_token = "99c9f388-12be-4b22-8f83-00b6f1e2d7d0"  # Maneet's token
        self.cs = ChemSpider(self.security_token, user_agent="StudentResearcher, ChemSpiPy 1.0.5, Python 3.6")

        print('--Populator Initialized--')

    def reactions_and_species(self, reactions_tsv, output_hdf5, output_reaction_df, output_species_df):
        """
        Reads the reactions from a TSV file and creates a dataframe out of it with
        2 columns containiting reactants and products as a list for each reaction,
        another 2 columns containing species ID as a list corresponding to each reaction.
        Also stores a dataframe containing just the unique reactants and their corresponding reactant ID
        :param reactions_tsv: TSV file containing reactions
        :param output_hdf5: Output HDF5 where final Dataframe will be stored
        :param output_reaction_df: Name by which the final reactions dataframe will be stored inside the output HDF5 file
        :param output_species_df: Name by which the final species dataframe will be stored inside the output HDF5 file
        :return: None
        """

        # Reading the input tsv in a data frame
        self.reactions_dataframe = pd.read_csv(reactions_tsv, header=0, index_col=0, sep="\t")

        # Cleaning individual reactants
        reactant_as_names = [""]*len(self.reactions_dataframe.index)  # Column that will store cleaned reactants names
        product_as_names = [""]*len(self.reactions_dataframe.index)  # Column that will store cleaned products names
        idx = 0
        unique_species_set = set()
        for _, row in self.reactions_dataframe.iterrows():
            reactant_as_names[idx] = [ele.translate(self.translator).strip() for ele in row['Reactants']
                .replace(" + ", "$").replace("â‰¡", "#").split("$")]
            product_as_names[idx] = [ele.translate(self.translator).strip() for ele in str(row['Products'])
                .replace(" + ", "$").replace("â‰¡", "#").split("$")]
            unique_species_set.update(reactant_as_names[idx])
            unique_species_set.update(product_as_names[idx])
            idx = idx + 1

        # Appending the column containing cleaned reactants and products to the 'Reaction' dataframe
        self.reactions_dataframe['Reactants_List'] = reactant_as_names
        self.reactions_dataframe['Products_List'] = product_as_names

        # Set doesn't preserve the order; the order of element may differ from the order they were added
        # into the set. So, converting set to list, then sorting it so that we always get same order
        # and consequently, same Species ID
        unique_species_list = list(unique_species_set)
        unique_species_list.sort()

        # Converting Species List to Species Dict (list will have unique species by default)
        self.unique_species_dict = {}
        for idx, ele in enumerate(unique_species_list):
            self.unique_species_dict[ele] = idx

        # Converting individual reactants to reactant ID
        reactants_as_sids = [""]*len(self.reactions_dataframe.index)  # Column that will store cleaned reactants IDs
        products_as_sids = [""]*len(self.reactions_dataframe.index)  # Column that will store cleaned reactants IDs
        idx = 0
        for _, row in self.reactions_dataframe.iterrows():
            reactants_as_sids[idx] = [self.unique_species_dict[ele] for ele in row['Reactants_List']]
            products_as_sids[idx] = [self.unique_species_dict[ele] for ele in row['Products_List']]
            idx = idx + 1

        # Appending the column containing cleaned reactants' RIDs to the reaction dataframe
        self.reactions_dataframe['Reactants_SIDs_List'] = reactants_as_sids
        self.reactions_dataframe['Products_SIDs_List'] = products_as_sids

        # Writing the Reaction Dataframe to a HDF5 file
        self.reactions_dataframe.to_hdf(path_or_buf=output_hdf5, key=output_reaction_df, mode='a')

        # Writing unique reactants into a data frame
        just_the_keys = unique_species_list
        just_the_values = range(len(unique_species_list))
        input_to_reactant_df = {'Species': just_the_keys, 'SID': just_the_values}
        self.species_df = pd.DataFrame(data=input_to_reactant_df)
        self.species_df = self.species_df.set_index('SID')

        # Writing the Reactant Dataframe to a HDF5 file
        self.species_df.to_hdf(path_or_buf=output_hdf5, key=output_species_df, mode='a')

        print('-- DataFrames Created and Stored in {} --'.format(output_hdf5))

    @staticmethod
    def print_from_hdf5(hdf5_store, dataframe_key, lines=5):
        """
        Reads the first few lines of a dataframe stored inside an HDF5 file
        :param hdf5_store: HDF5 file storing the dataframe
        :param dataframe_key: Name of the datframe inside the HDf5 file
        :param lines: How many lines to read
        :return: None
        """
        data_store = pd.HDFStore(hdf5_store)  # Opening the HDF5 file
        read_dataframe = data_store[dataframe_key]  # Reading the dataframe
        data_store.close()
        print(read_dataframe.head(lines))

    def set_and_initialize_token(self, input_token):
        """
        Stores you ChemSpider security token as an object attribute and Associate your token to the ChemSpider api
        :param input_token: your security token (for ChemSpider)
        :return: None
        """
        self.security_token = input_token
        self.cs = ChemSpider(self.security_token)

    def fetch_csid_and_messages(self, output_hdf5, output_reactant_df):
        """
        Augments reac_df with ChemSpider CSID and query status results.
        :param output_hdf5: HDF5 file where Reactant DataFrame is stored
        :param output_reactant_df: Name of the Reactant DataFrame
        :return: None
        """

        # Read DataFrame from HDF5File
        data_store = pd.HDFStore(output_hdf5)  # Opening HDF5 File
        reactant_df = data_store[output_reactant_df]  # Reading the desired DF
        data_store.close()

        # Intitialize the columns that will be appended to the datframe
        num_results = [0]*len(reactant_df.index)
        csids = [""]*len(reactant_df.index)
        messages = [""]*len(reactant_df.index)

        # Populate the columns initialized above with the ChemSpider API results
        idx = 0
        for _, row in reactant_df.iterrows():

            out_result = self.cs.search(row['Species'])  # Requesting the ChemSpider API for info on the input reactant
            out_result.wait()  # Waiting until the API response is completely received
            result_length = len(list(out_result))  # Number of matches for a particular query
            num_results[idx] = result_length  # Storing the number of the matches obtained above

            csid_list = []  # Initializing a list that will containg the csid matches for a particular query
            if result_length > 0:
                for ele in out_result:
                    csid_list.append(ele.csid)
            csids[idx] = csid_list  # CSID Matches obtained against the input query
            messages[idx] = out_result.message  # Storing the messsage obtained
            print(idx)  # Just to check retrieval status
            idx = idx + 1

        # Augmenting to the dataframe with ChemSpider results
        reactant_df['NumResults'] = num_results  # Adding a new column storing number of matches
        reactant_df['CSIDs'] = csids  # Adding a new column storing CSID matches
        reactant_df['Message'] = messages  # Adding a new column storing query message

        # Store the appended dataframe back to the to the parent HDF5 file
        reactant_df.to_hdf(path_or_buf=output_hdf5, key=output_reactant_df, mode='a')

    def smile_it(self, output_hdf5, output_reactant_df):
        """
        Reads Pandas dataframe, augment it with SMILE strings and MOL2d data and store it
        in a given HDF5 file under the specified dataFrame
        :param output_reactant_df: Name of output dataframe; type: str
        :param output_hdf5: Path of HDF5 file; type: str
        :return: None
        """

        # Read DataFrame from HDF5File
        data_store = pd.HDFStore(output_hdf5)  # Opening HDF5 File
        reactant_df = data_store[output_reactant_df]  # Reading the desired DF
        data_store.close()

        # List storing SMILE representation
        extended_info = [""]*len(reactant_df.index)  # molecular mass, inchi key, smile string, etc.
        mol2d_data = [""]*len(reactant_df.index)  # Mol2D data string

        # Accepted categories for pulling SMILE strings
        accepted_categories = ["Found by approved synonym",
                               "Found by conversion query string to chemical structure (full match)"]

        # Aughmenting DF with Molecular Info
        idx = 0
        for _, row in reactant_df.iterrows():
            if row['Message'] in accepted_categories:
                under_radar = row['CSIDs']  # CSID list under radar
                length_under_radar = len(under_radar)
                if length_under_radar == 0:
                    pass
                elif length_under_radar > 0:
                    try:
                        extended_info[idx] = str(self.cs.get_extended_compound_info(under_radar[0]))
                        mol2d_data[idx] = self.cs.get_original_mol(under_radar[0])
                        print(idx)  # Status check
                    except Exception as e:
                        # Handling Connection Error
                        print(e)
                        print("Error seen at", idx, "with compound", under_radar[0])
                        # // Handling premature exit by saving whatever we have obtained
                        reactant_df['ExtendedInfo'] = extended_info
                        reactant_df['Mol2d'] = mol2d_data
                        # Store the appended dataframe back to the to the parent HDF5 file
                        reactant_df.to_hdf(path_or_buf=output_hdf5, key=output_reactant_df, mode='a')
                        return
            else:
                pass
            idx = idx + 1

        # If everything goes well, augmenting to the dataframe with ChemSpider results
        reactant_df['ExtendedInfo'] = extended_info
        reactant_df['Mol2d'] = mol2d_data

        # Store the appended dataframe back to the to the parent HDF5 file
        reactant_df.to_hdf(path_or_buf=output_hdf5, key=output_reactant_df, mode='a')

        return

    @staticmethod
    def status_check(output_hdf5, output_reaction_df, output_species_df):
        """
        Assign scores to each reactant
        :param output_hdf5: Output HDF5 file
        :param output_reaction_df: Reactions Dataframe
        :param output_species_df: Reactants Dataframe
        :return: None
        """

        # Reading dataframes from the HDF5 file
        data_store = pd.HDFStore(output_hdf5)  # Opening the HDF5 file
        reaction_dataframe = data_store[output_reaction_df]  # Reading the dataframe
        species_dataframe = data_store[output_species_df]  # Reading the dataframe
        data_store.close()

        # Creating and Appending Column which will contain the score of the reactants
        score = [0]*len(species_dataframe.index)
        species_dataframe['Scores'] = score

        # Assigning Scores to species
        for _, row in reaction_dataframe.iterrows():
            list_under_consider = row['Reactants_SIDs_List'] + row['Products_SIDs_List']
            for species in list_under_consider:
                species_dataframe.at[species, 'Scores'] = species_dataframe.at[species, 'Scores'] + 1

        # Updating dataframe in the HDF5 file
        species_dataframe.to_hdf(path_or_buf=output_hdf5, key=output_species_df, mode='a')

        print("-- Scores Assigned --\n")

    def fetch_more_smiles(self, output_hdf5, output_reactant_df):
        """
                Reads Pandas dataframe, augment it with SMILE strings and MOL2d data and store it
                in a given HDF5 file under the specified dataFrame.
                The reactants augmented are defined by the user via a custom criteria
                :param output_reactant_df: Name of output dataframe; type: str
                :param output_hdf5: Path of HDF5 file; type: str
                :return: None
                """

        # Read DataFrame from HDF5File
        data_store = pd.HDFStore(output_hdf5)  # Opening HDF5 File
        reactant_df = data_store[output_reactant_df]  # Reading the desired DF
        data_store.close()

        # Aughmenting DF with Molecular Info
        for idx, row in reactant_df.iterrows():
                if len(row['CSIDs']) == 1 and row['Mol2d'] == "":  # Custom Criteria
                        reactant_df.at[idx, 'ExtendedInfo'] = str(self.cs.get_extended_compound_info(row['CSIDs'][0]))
                        reactant_df.at[idx, 'Mol2d'] = self.cs.get_original_mol(row['CSIDs'][0])
                        print(idx)  # Status check

        # Store the appended dataframe back to the to the parent HDF5 file
        reactant_df.to_hdf(path_or_buf=output_hdf5, key=output_reactant_df, mode='a')

        return

    @staticmethod
    def reaction_status(output_hdf5, output_reaction_df, output_species_df):
        """
        Assign boolean flags to each reaction
        :param output_hdf5: Output HDF5 file
        :param output_reaction_df: Reactions Dataframe
        :param output_species_df: Reactants Dataframe
        :return: None
        """

        # Reading dataframes from the HDF5 file
        data_store = pd.HDFStore(output_hdf5)  # Opening the HDF5 file
        reaction_dataframe = data_store[output_reaction_df]  # Reading the dataframe
        species_dataframe = data_store[output_species_df]  # Reading the dataframe
        data_store.close()

        # Creating and Appending Column which will contain the score of the reactants
        flag_product_available = [True]*len(reaction_dataframe.index)
        status_50 = [False]*len(reaction_dataframe.index)
        status_75 = [False]*len(reaction_dataframe.index)
        status_100 = [False]*len(reaction_dataframe.index)
        reaction_dataframe['Products_Available'] = flag_product_available
        reaction_dataframe['Status_50'] = status_50
        reaction_dataframe['Status_75'] = status_75
        reaction_dataframe['Status_100'] = status_100

        # Assigning Scores to reactants
        for _, row in reaction_dataframe.iterrows():

            # Checking whether products are available for that reaction
            for prod in row['Products_List']:
                if prod in ['Products', 'Other Products']:
                    reaction_dataframe.at[_, 'Products_Available'] = False
                    break

            # Checking whether species occur in more than 50, 75, and 100 reactions.
            spec_id_len = len(row['Reactants_SIDs_List'] + row['Products_SIDs_List'])
            marker_50 = 0
            marker_75 = 0
            marker_100 = 0
            for spec_id in (row['Reactants_SIDs_List'] + row['Products_SIDs_List']):
                if species_dataframe.at[spec_id, 'Scores'] >= 50:
                    marker_50 = marker_50 + 1
                    if species_dataframe.at[spec_id, 'Scores'] >= 75:
                        marker_75 = marker_75 + 1
                        if species_dataframe.at[spec_id, 'Scores'] >= 100:
                            marker_100 = marker_100 + 1
                else:
                    break

            # Changing Markers based on whether the reaction qualified the specified criterion.
            if marker_50 == spec_id_len:
                reaction_dataframe.at[_, 'Status_50'] = True
                if marker_75 == spec_id_len:
                    reaction_dataframe.at[_, 'Status_75'] = True
                    if marker_100 == spec_id_len:
                        reaction_dataframe.at[_, 'Status_100'] = True

        # Updating dataframe in the HDF5 file
        reaction_dataframe.to_hdf(path_or_buf=output_hdf5, key=output_reaction_df, mode='a')

        print("-- Boolean Flags Assigned --\n")

    @staticmethod
    def get_pubchem_data(output_hdf, species_df_key):
        """
        Augments the species dataframe with pubchem data based on CID
        :param output_hdf: File where the older species df is read from and
        where the updated species df will be stored
        :param species_df_key: species df key in the output_hdf
        :return: None
        """
        species_df = pd.read_hdf(output_hdf, species_df_key)  # Reading the DF

        # Creating BondsInfo column if it doesnt exist already
        if 'BondsInfo' not in species_df.columns:
            bonds_info = [""]*(len(species_df.index))
            species_df['BondsInfo'] = bonds_info

        for idx, row in species_df.iterrows():
            if not math.isnan(row['CID']) and row['BondsInfo'] == "":
                cid = int(row['CID'])
                if cid > 0:  # Handling valid CIDs
                    r = requests.get('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/record/json'.format(cid))
                    species_df.at[idx, 'BondsInfo'] = r.text
                    print("{} done".format(idx))

        # Writing back to HDF5 file
        my_hdf = pd.HDFStore(output_hdf)
        my_hdf[species_df_key] = species_df
        my_hdf.close()


# Code Run Check
# my_populator = Populator()
# my_populator.reactions_and_reactants('DataFiles/kineticsDB_Parent/reactions.tsv', 'NewGenOutput/Fledged.h5', 'Reactions', 'Reactants')
# Populator.print_from_hdf5('NewGenOutput/Fledged.h5', 'Reactions')
# Populator.print_from_hdf5('NewGenOutput/Fledged.h5', 'Reactants')
# my_populator.fetch_csid_and_messages('NewGenOutput/Fledged.h5', 'Reactants')
# my_populator.smile_it('NewGenOutput/Fledged.h5', 'Reactants')
# Populator.status_check('NewGenOutput/Fledged.h5', 'Reactions', 'Reactants')
# my_populator.fetch_more_smiles('NewGenOutput/Fledged.h5', 'Reactants')
# Populator.status_100('NewGenOutput/Fledged.h5', 'Reactions', 'Reactants')

# Code Run Check
# my_populator = Populator()
# my_populator.reactions_and_species('DataFiles/kineticsDB_Parent/reactions.tsv', 'NewGen2Output/NewGen.h5', 'Reactions', 'Species')
# Populator.status_check('NewGen2Output/NewGen.h5', 'Reactions', 'Species')
# Populator.reaction_status('NewGen2Output/NewGen.h5', 'Reactions', 'Species')
# Populator.print_all_to_excel('NewGen2Output/NewGen.h5')

# Populator.get_pubchem_data("PreliminaryOutput/DemoGenerated/DataDF.h5", "Species")
