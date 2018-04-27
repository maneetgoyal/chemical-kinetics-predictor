# Script to tranfer the CIDs from old Species.xlsx to new Species.xlsx
import pandas as pd
import math

def transfer_cid(old_xlsx_file_path, hdf_file_path, species_df_key, new_xlsx_file_path):
    """
    Transfers the CIDs of old species.xlsx to new_species.xlsx and Species dataframe of
    the old HDF
    :param old_xlsx_file_path: As title
    :param hdf_file_path: As title
    :param species_df_key: As title
    :param new_xlsx_file_path: As title
    :return: None
    """

    # Reading xlsx files into pandas df
    old_df = pd.read_excel(old_xlsx_file_path, header=0)
    new_df = pd.read_hdf(hdf_file_path, species_df_key)

    # Adding CID column ot new df
    new_df['CID'] = [""]*len(new_df.index)

    # Setting Species name as index for efficiency
    new_df = new_df.reset_index()
    new_df = new_df.set_index(keys="Species", verify_integrity=True)

    # Transfer CID
    transfer_count = 0
    for idx, row in old_df.iterrows():
        if not math.isnan(row['CID']) and row['CID'] != "":
            new_df.at[row['Species'], 'CID'] = row['CID']
            transfer_count = transfer_count + 1
    print('--Status--')
    print('--{} Transfers Made--'.format(transfer_count))

    # Restoring the df index
    new_df = new_df.reset_index()
    new_df = new_df.set_index(keys="SID")

    # Creating the new_xlsx,
    new_df.to_excel(new_xlsx_file_path)

    #  Updating the HDF file also
    new_df.to_hdf(hdf_file_path, "Species")


transfer_cid("species.xlsx", '../PreliminaryOutput/DemoGenerated/DataDF.h5', 'Species', "new_species.xlsx")
