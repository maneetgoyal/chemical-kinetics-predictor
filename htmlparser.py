import bs4
import pandas as pd
from urllib.request import urlopen
from urllib.request import urlretrieve
import zipfile

class TableCreator:
    """
    Converts HTML to .txt
    Output Format: Reaction-ID, Link, #Records, Reactants, Products
    """
    def __init__(self, html_path):
        """
        Function to create a soup using BeautifulSoup out of the given HTML file path.
        :param html_path: Path of the HTML file which contains all the reactions. Format: str
        """
        try:
            pointer_to_html = open(html_path)  # Pointing to HTML file
            html_data = pointer_to_html.read()  # Reading all HTML data in memory
            pointer_to_html.close()
            print("----Status----")
            print("HTML File read in memory.")
        except IOError:
            print("Input HTML file couldn't be opened.")
        else:
            soup = bs4.BeautifulSoup(html_data, "html5lib")  # Creating soup out of the in-memory database
            print("----Status----")
            print("Soup created out of the HTML file")
            table = soup.find('table').find_all('table')
            self.reaction_table = table[1]  # Contains the table which has all the reactions stored.
            self.txt_created = False  # Because CSV file still not created
            print("----Status----")
            print("Done")

    def extrct_rxn_to_txt(self, tsv_path):
        """
        Function to write reactions from in-memory HTML table into a file pointed to by txt_path.
        :param tsv_path: Path of the .tsv file to which your reactions should be extracted. Format: str
        :return: None
        """
        try:
            pointer_to_tsv = open(tsv_path, "w", encoding='utf-8')
        except IOError:
            print("Output CSV file couldn't be opened.")
        else:
            pointer_to_tsv.write("RID" + "\t" + "Reaction Link" + "\t" + "Records" + "\t" + "Reactants" + "\t" + "Products" + "\n")
            idx = 1  # Reaction ID
            for each_row in self.reaction_table.find_all('tr'):
                three_columns = each_row.find_all('td')  # Extract all three columns
                if len(three_columns) == 3:
                    first_column = three_columns[0]  # First column has HREF and #Records
                    reaction_detail = [ele.strip() for ele in three_columns[2].text.split("â†’")]  # Third column has
                    # reactants and products
                    pointer_to_tsv.write((str(idx) + "\t" + first_column.find('a')['href'] + "\t" + first_column.text.split(" ")[0] +
                                         "\t" + reaction_detail[0] +
                                         "\t" + reaction_detail[1] + "\n"))
                    idx = idx + 1
            pointer_to_tsv.close()  # Closing the output file
            self.txt_created = True
            print("----Status----")
            print("All reactions written into the output .tsv file in TSV format.")

# Creating table in-memory from the input HTML file
# myTableCreator = TableCreator(r"G:\References\MS1\Spring2018\CHBE\Project\NIST Chemical Kinetics Database.html")
# Creating output .tsv file after reading all the reactions
# myTableCreator.extrct_rxn_to_txt("reactions.tsv")

class RxnDetailsExtractor:
    """
        Read HREFs in a tsv, then fetches reaction details from the corresponding HTML
        Output Format: Link, #Records, Reactants, Products
    """
    def __init__(self, tsv_file_path):
        """
        Function to read reactions from the given .tsv file into a pandas dataframe.
        :param tsv_file_path: Path of the .tsv file in which your reactions are stored. Format: str
        """
        self.reactions_df = pd.read_csv(tsv_file_path, usecols=[0, 1, 2], sep='\t', index_col=0, header=0)  # Reactions
        # dataframe, containing reaction links.
        print("----Status----")
        print("Reactions '.tsv' file read into a Pandas DataFrame.")
        self.tsv_created = False  # Reactions' record are not yet written into a TSV file.
        self.tsv_reference_created = False  # Reference Reactions are not yet written into a TSV file.

    def extrct_rec_to_tsv(self, tsv_file_path, ref_tsv_file_path):
        """
        Function to read all the records coresponding to all the reaction links into the given csv file.
        :param tsv_file_path: Path to .tsv file where all the reactions' records will be stored. Format: str
        :param ref_tsv_file_path: Path to .tsv file where all the reference reactions' will be stored. Format: str
        :return: None
        """
        try:  # Opening output tsv containing reaction records
            pointer_to_tsv = open(tsv_file_path, "a", encoding='utf-8')
        except IOError:
            raise IOError("Output reactions' record TSV file couldn't be opened.")
        else:
            # Writing Headers
            pointer_to_tsv.write("RecordID" + "\t" + "RID" + "\t" + "RecordType" + "\t" + "Squib" + "\t" +
                                 "PaperDetails" + "\t" + "Temperature(K)" + "\t" +
                                 "FrequencyFactor, A" + "\t" + "TemperatureRatioExponent, n" +
                                 "\t" + "ActivationEnergy, J/mol" + "\t" + "RateConstant, k(298 K)" +
                                 "\t" + "ReactionOrder" + "\n")
            idx = 1  # Reaction ID
            record_type = ""  # Not yet defined; will be defined during iterations as encouded in the table rows

        try:  # Opening output tsv containing reference reactions
            pointer_to_ref_tsv = open(ref_tsv_file_path, "a", encoding='utf-8')
        except IOError:
            raise IOError("Output reference reactions' TSV file couldn't be opened.")
        else:
            # Writing Headers
            pointer_to_ref_tsv.write("RecordID" + "\t" + "RID" + "\t" + "Squib" + "\t" + "ReferenceReaction" + "\n")

        for rid, row in self.reactions_df.iterrows():
            try:
                reaction_html = urlopen(row['Reaction Link'])  # Opening html in memory
                print(reaction_html)
            except ConnectionError:
                print("Could not go to link {}".format(row['Reaction Link']))
                return None
            reaction_soup = bs4.BeautifulSoup(reaction_html, "html5lib")  # Making soup out of it
            try:
                table_rows = reaction_soup.find('table').find_all('table')[5].find_all('tr')  # Finding all the rows
            except IndexError:
                print("Bad HTML | Something is wrong with this HTML {}. Proceeding with the next one.".format(row['Reaction Link']))
                continue
            for ele_rows in table_rows[1:]:  # Iterating through all the rows
                all_cols = ele_rows.find_all('td')  # Finding all the columns in the concerned rows
                if len(all_cols) == 1 and all_cols[0].text != "\xa0":  # Handling cases when record type
                    # may have been encountered
                    record_type = all_cols[0].text  # updating record type if it was encountered
                elif len(all_cols) == 3:  # for handling reference reaction data
                    td_list = ele_rows.find_all('td')  # Creating list of all the columns
                    reference_rxn = td_list[1].text[21:]  # Storing reference reaction
                    try:
                        pointer_to_ref_tsv.write(str(idx) + "\t" + str(rid) + "\t" + squib + "\t" +
                                                 reference_rxn + "\n")
                        # Writing all those values to the output tsv file
                    except IOError:
                        print("Could not write Record ID {} for RID {}".format(idx, rid))
                    continue
                elif len(all_cols) > 1:  # Handling case when data entry is encountered
                    td_list = ele_rows.find_all('td')  # Creating list of all the columns
                    squib = "http://kinetics.nist.gov" + td_list[2].find('a')["href"]  # Squib URL
                    paper_details = td_list[2].find('a')["onmouseover"][9:-10]  # Paper details
                    temperature = td_list[4].text
                    frequency_factor = td_list[6].text
                    temperature_exponent = td_list[8].text
                    activation_energy = td_list[10].text
                    rate_constant = td_list[12].text
                    reaction_order = td_list[14].text
                    try:
                        pointer_to_tsv.write(str(idx) + "\t" + str(rid) + "\t" + record_type + "\t" + squib + "\t" +
                                             paper_details + "\t" + temperature + "\t" + frequency_factor + "\t" +
                                             temperature_exponent + "\t" + str(activation_energy) + "\t" + rate_constant
                                             + "\t" + reaction_order + "\n")  # Writing all those values to
                        # the output tsv file
                    except IOError:
                        print("Could not write Record ID {} for RID {}".format(idx, rid))
                    idx = idx + 1  # incrementing record ID
            print("RID: {} parsed".format(rid))
        pointer_to_tsv.close()
        pointer_to_ref_tsv.close()
        self.tsv_created = True
        self.tsv_reference_created = True

    @staticmethod
    def send_records_to_hdf(records_file, dataframe_key, output_hdf):
        """
        This method transfer the records (xlsx) file to a pandas dataframe and then store that DF to a HDF5.
        :param records_file: Input xlsx file
        :param dataframe_key: Key (name) of the dataframe created from records_csv
        :param output_hdf: Output HDF path
        :return: None
        """
        # Reading records_csv as a DF
        records_df = pd.read_excel(records_file, index_col=0, header=0, na_values=["ï¿½", "�"])
        pointer_to_df = pd.HDFStore(output_hdf)  # Opening the output_hdf file
        pointer_to_df.put(dataframe_key, records_df)  # Putting DF into HDF5
        pointer_to_df.close()  # Closing HDF5 file

# # Read HREFs into the input tsv file path
# myRxnExtrator = RxnDetailsExtractor(r"PreliminaryOutput\reactions.tsv")
# print(myRxnExtrator.reactions_df)
# # Does all the scraping from the url suplied from tsv file path (argument 1)
# #  and writes the data to the inpput file path (argument 2)
# myRxnExtrator.extrct_rec_to_tsv("records.tsv", "ref_reaction.tsv")

class Urldownloader:
    """
    Downloads the network objected denoted by the url named as object variable, "pointing_to"
    """
    def __init__(self):
        """
        Initializes the url_dowloader object. No significant action is taking place.
        """
        self.pointing_to = None
        self.file_path_with_name = None
        self.info_tuple = None  # will contain tuple: (local file name under which object can be found,
        # meta-information of the page, such as headers)

    def set_url_and_path(self, input_url, output_path_file_name=None):
        """
        Sets the input url and output file name + path
        :param input_url: Url pointing towards the object you want to download
        :param output_path_file_name: Path along with the file name of output file
        :return: None
        """

        self.pointing_to = input_url  # Setting the input url

        # Setting the output path
        if output_path_file_name is None and self.file_path_with_name is None:  # If output path is neither pre-stored
            # in the object nor specified in the function call.
            raise ValueError("Output file path is not specified.")
        elif output_path_file_name is None and self.file_path_with_name is not None:  # If output path is not specifed
            # in the function call but a pre-stored value exists.
            pass
        else:  # If output file path is specifed
            self.file_path_with_name = output_path_file_name

    def retrieve_file(self):
        """
        Downloads the netowrk objects which your input url is pointing to at you putput file path
        :return: None
        """
        self.info_tuple = urlretrieve(self.pointing_to, self.file_path_with_name)

    @staticmethod
    def unzip_it(input_file_path, output_folder):
        """
        Unzips a given file
        :param input_file_path: input file path, type: str
        :param output_folder: output folder path, type: str
        :return: None
        """
        pointer_to_file_to_be_unzipped = zipfile.ZipFile(input_file_path)
        pointer_to_file_to_be_unzipped.extractall(output_folder)

