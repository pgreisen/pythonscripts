import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import argparse, re
from datetime import date
import datetime
import json

class ParserSynergyPlateReader:

    def get_well_label(self,key):
        alphabet = {
            0: 'A',
            1: 'B',
            2: 'C',
            3: 'D',
            4: 'E',
            5: 'F',
            6: 'G',
            7: 'H'
        }

        return alphabet[key]


    def __init__(self):
        self.json_dic = {}
        self.excelfile = ""
        self.plateformatfile = ""
        self.today = datetime.date.today().isoformat()
        self.debug = False
        # Meta-index
        self.metadata_index = 27
        self.numeric_columns = ['Read 2:485,515']

    def myconverter(self, o):
        if isinstance(o, datetime.datetime):
            return o.__str__()


    def assign_well_variant(self, list_of_plasmids, platetype=96 ):

        #ToDo make method here such that format of plate can change easily
        alphabet = {
            0: 'A',
            1: 'B',
            2: 'C',
            3: 'D',
            4: 'E',
            5: 'F',
            6: 'G',
            7: 'H'
        }
        dummy = 0
        cindex = 0
        platemap_ = {}
        lim_plate = len(list_of_plasmids) * len(alphabet.keys())
        while dummy < lim_plate:
            for i in list_of_plasmids:
                cindex += 1
                for j in alphabet.keys():
                    plateid = alphabet[j] + str(cindex)
                    if (i not in platemap_):
                        platemap_[i] = []
                    platemap_[i].append(plateid)
                    dummy += 1
        df_ = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in platemap_.items()]))  # B (Better)
        df_ = df_.melt(var_name='plasmidID', value_name='PlatemapID')
        return df_

    def clean_data_frame(self):
        '''

        :param excel file: with plate layout
        :return: pandas dataframe
        '''
        # First sheet contains the layout of the plates
        # It is explicit stated here to avoid confusion
        df = pd.read_excel(self.excelfile,sheet_name=0 )
        cols_ = df.columns
        df[ cols_[0] ].ffill(axis=0, inplace=True)
        # Remove all rows with all NaNs
        df.dropna(axis=0, how='all')
        cols = df.columns
        # initializing substring
        subs = 'StrainID'
        # using re + search()
        # to get string with substring
        res = [x for x in cols if re.search(subs, x)]

        plates_ = []
        df2 = df.drop(res, axis=1)
        df = df[res]

        # Remove all rows with all NaNs
        df = df.dropna(axis=0, how='all')
        for index, row in df.iterrows():
            if(self.debug):
                print("Index, row: ", index,row)
            tmpdf = self.assign_well_variant(list(row.dropna()))
            tmpdf["Plate_Number"] = df2.iloc[index, 0]
            tmpdf["Plate Info"] = df2.iloc[index, 2]
            plates_.append(tmpdf)

        df = pd.concat(plates_)
        df = df.dropna(axis=0)
        if(self.debug):
            df.to_excel(str(self.today)+"_test.xlsx")
            df2.to_excel(str(self.today)+"_df2_test.xlsx")

        return df,df2


    def pivot_table_from_cleaned_dataframe( self, df ):

        df = df.pivot(index='Unnamed: 0', columns='Unnamed: 1', values='Plate Info')
        df = df.reset_index()
        df.rename(columns={ 'Unnamed: 0': 'Plate_Number'}, inplace=True)

        params_ = {}
        for i, j in df.iterrows():
            tmplist = list(j)
            for k in tmplist:
                if (k not in params_.keys()):
                    key = k
                    params_[key] = []
                else:
                    params_[key].append(k)
        return df

    def get_merge_dataframe_objects(self, df, df2):
        merge_df1 = pd.merge(right=df, right_on='Plate_Number', left=df2, left_on='Plate_Number', suffixes=("", "_remove"))
        merge_df1.drop(axis=0, columns='Plate Info', inplace=True)
        return merge_df1

    def dilute_concentration(self, df, row="Row",conc="Conc (uM)",dilu="Dilution Factor"):
        return (df["Conc (uM)"].astype(float) / ((df["Row"].astype(int)-1) * df["Dilution Factor"].astype(float)))

    def dilute_columns(self, totdata):

        # ToDo
        dilution_factor_table = {
            'A': 0,
            'B': 1,
            'C': 2,
            'D': 3,
            'E': 4,
            'F': 5,
            'G': 6,
            'H': 7 }

        totdata["DL_factor"] = 1

        for i, j in zip(totdata.index, totdata["Column"]):
            totdata.iloc[i, totdata.columns.get_loc("DL_factor")] = dilution_factor_table[j]
        mask = (totdata['Row'] != "-1")
        # ToDo is this necessay?
        totdata_valid = totdata[mask]

        totdata['c'] = totdata['Conc (uM)']
        totdata.loc[mask, 'c'] = totdata["Conc (uM)"].astype(float) / (totdata["Dilution Factor"].astype(float) ** (totdata["DL_factor"].astype(int)))

        totdata.rename(columns={
            'Conc (uM)': 'StartConc (uM)',
            'c': 'Conc (uM)'}, inplace=True)
        totdata['Conc (uM)'] = totdata['Conc (uM)'].map('{:,.2E}'.format)
        totdata = totdata[['Plate_Number', 'Ligands', 'PlateID', 'plasmidID', 'PlatemapID', 'Read 1:600',
                           'Read 2:485,515', 'Conc (uM)']]
        return totdata

    def get_assure_numeric_values(self, df):
        '''
        :param df: pandas dataframe
        :return: dataframe where non float values have been set to nan
        '''
        for i in self.numeric_columns:
            df[i] = pd.to_numeric(df[i], errors='coerce')
        return df

    def get_plate_reader_results(self,excelfile, excel_sheet):
        df = pd.read_excel(excelfile,sheet_name=excel_sheet, header=None)
        df.dropna(how="all", inplace=True)
        df.reset_index(inplace=True,drop=True)

        for i, j in zip(df.index, df[0]):
            if (j == 'Actual Temperature:'):
                self.metadata_index = i + 2
                if(self.debug):
                    print("Metadata_index: ", self.metadata_index)

        # get the json information
        jsondic = {}
        for i, j in zip(df[:self.metadata_index-1][0], df[:self.metadata_index-1][1]):
            if (i not in jsondic.keys() and not pd.isna(i)):
                key = i
                jsondic[key] = []
            if(pd.isna(j)):
                j = None
            jsondic[key].append(j)
        self.json_dic[jsondic["Plate Number"][0]] = jsondic

        # ToDo
        plate = df[self.metadata_index:]
        plate = plate.dropna(axis=1, how='all')
        plate = plate.dropna(axis=0, how='all')
        plate = plate.reset_index(drop=True)
        plate_normalized_plate = plate.iloc[::2]
        plate_normalized_plate = plate_normalized_plate.set_index(1)
        plate_normalized_plate.index.name = 'Index'
        for i in plate_normalized_plate.columns:
            plate_normalized_plate.rename(columns={i: str(i - 1)}, inplace=True)
        plate_normalized_plate = plate_normalized_plate.reset_index()
        plate_not_normalized_plate = plate[1:].iloc[::2]
        plate_not_normalized_plate = plate_not_normalized_plate.set_index(1)
        plate_not_normalized_plate.index.name = 'Index'
        for i in plate_not_normalized_plate.columns:
            plate_not_normalized_plate.rename(columns={i: str(i - 1)}, inplace=True)
        plate_not_normalized_plate = plate_not_normalized_plate.reset_index()
        plate_not_normalized_plate['Index'] = plate_normalized_plate['Index']
        plate_normalized_plate = plate_normalized_plate.melt(id_vars=['Index'], value_name="Value")
        plate_normalized_plate["PlatemapID"] = plate_normalized_plate["Index"] + plate_normalized_plate["variable"].astype(str)
        plate_normalized_plate.rename(columns={'Value': 'Read 1:600'}, inplace=True)
        plate_not_normalized_plate = plate_not_normalized_plate.melt(id_vars=['Index'], value_name="Value")
        plate_not_normalized_plate["PlatemapID"] = plate_not_normalized_plate["Index"] + plate_not_normalized_plate["variable"].astype(str)
        plate_not_normalized_plate.rename(columns={'Value': 'Read 2:485,515'}, inplace=True)

        plate_not_normalized_plate["Plate_Number"] = jsondic["Plate Number"][0]
        plate_reformattet = pd.merge(right=plate_not_normalized_plate, right_on="PlatemapID",left=plate_normalized_plate, left_on="PlatemapID")
        plate_reformattet = plate_reformattet[['Read 1:600', 'PlatemapID', 'Read 2:485,515', 'Plate_Number']]
        # 96 wells plus header hence 96 and not 0:95
        return plate_reformattet.iloc[0:96,:]

    def main(self):

        # ToDo go through each sheet

        parser = argparse.ArgumentParser(description="")
        # get the initial rosetta design as input
        parser.add_argument("-f", "--file", dest="excelfile", help="Excel file with plate layout")
        parser.add_argument("-p", "--plateformatfile", dest="plateformatfile", help="Excel file")

        args_dict = vars(parser.parse_args())
        for item in args_dict:
            setattr(self, item, args_dict[item])
        if(self.plateformatfile == None):
            self.plateformatfile = self.excelfile

        df3s = []
        xls = pd.ExcelFile(self.plateformatfile)
        # Now you can list all sheets in the file
        for i in xls.sheet_names[1:]:
            df3s.append(self.get_plate_reader_results( self.plateformatfile, i))

        df3 = pd.concat(df3s)
        df3.reset_index(inplace=True)
        if(self.debug):
            df3.to_excel("debug.xlsx")
        # take self.excelfile as input which has been assign with the parser
        df, df2 = self.clean_data_frame()
        df2 = self.pivot_table_from_cleaned_dataframe(df2)
        if (self.debug):
            df2.to_excel("debug_excel.xlsx")

        tot_df = self.get_merge_dataframe_objects(df, df2)

        totdata = pd.merge(right=tot_df, right_on=["Plate_Number","PlatemapID"], left=df3, left_on=["Plate_Number","PlatemapID"])
        totdata["Row"] = totdata["PlatemapID"].str[1:]
        totdata["Column"] = totdata["PlatemapID"].str[0]
        totdata = self.dilute_columns( totdata )

        if (self.debug):
            totdata.to_excel("debug_excel2.xlsx")

        # set dilution to zero
        # ToDo vectorize it and add test that the value exists on the plate
        reformat_dilutions = dict(zip(df2['Plate_Number'], df2['Zero']))
        for z in reformat_dilutions:
            for i,j,k in zip(totdata.index, totdata['PlatemapID'], totdata['Plate_Number']):
                if( j[0] == reformat_dilutions[z] and k == z ):
                    totdata.iloc[i,totdata.columns.get_loc('Conc (uM)')] = '{:.2e}'.format(0)

        # type check columns problems could be overflow etc.
        totdata  = self.get_assure_numeric_values( totdata)
        totdata.to_excel(self.today + "_format_parser_input.xlsx", index=False, na_rep='Inf')

        with open(self.today+"_header_information.json", "w") as test:
            json.dump(self.json_dic, test, indent=4, sort_keys=True, default=self.myconverter, allow_nan=False)



if __name__ == "__main__":
   run = ParserSynergyPlateReader()
   run.main()