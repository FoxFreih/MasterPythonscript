import datetime
import time
import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None

#check if the key is in the dictionary
def checkKey(dict, key):
    if key in dict.keys():
        return 1
    else:
        return 0


def remove_last_characters(key):
    size = 0
    for i in key:
        if i != '.':
            size += 1
        else:
            key_size = len(key)
            remove_size = key_size - size
            # remove all the characters which came right after the point
            new_key_string = key[:key_size - remove_size]
            return new_key_string


# function that build a dictionary of minimal E-value
def insert(dict, key, i,csvFile):
    Evalue = csvFile['E-value'][i]
    if checkKey(dict, key):
        if Evalue < dict[key]:
            dict[key] = Evalue
    if checkKey(dict, key) == 0:
        dict[key] = Evalue
def insert_FPKM_TPM_values(RSEMfile,FPKMindex,TPMindex,key,i):
    FPTP = RSEMfile.loc[RSEMfile['transcript_id(s)'] == key, ['FPKM', 'TPM']]
    RSEM.at[i, FPKMindex] = FPTP['FPKM']
    RSEM.at[i, TPMindex] = FPTP['TPM']
#function to find middle of three distinct numbers
def middleOfThree(TempList):
    # Checking for number at the middle index
    count = TempList.count(0) #if we have 2 zeros return the max number as middle
    if count == 2:
        return max(TempList[0], TempList[1], TempList[2])
    if ((TempList[0] < TempList[1] and TempList[1] < TempList[2]) or (TempList[2] < TempList[1] and TempList[1] < TempList[0])):
        return TempList[1]

        # Checking for number at the first index
    if ((TempList[1] < TempList[0] and TempList[0] < TempList[2]) or (TempList[2] < TempList[0] and TempList[0] < TempList[1])):
        return TempList[0]
    else:
        return TempList[2]
def Avg(TempList,middle,scale):
    sum=0
    n=0
    for val in TempList:
        if middle!=0:
            if val/middle!=1: #find the other two values which aren't the middle
                if  middle/scale <= val <= middle*scale: #check if the val is in the range
                    sum+=val
                    n+=1
    if n>=1 and sum!=0:
        sum+=middle
        #if 0 in TempList: #n+1 in case if we have val!=0 in the range and val=0 (0 should be calculated in the average)
            #n+=1
        n+=1 #if we have a val that is in the range, enter the value of the middle in the calculating.n+1 for the middle
        avg=sum/n
        avg=round(avg, 2)
        return avg
    else:
        return middle/2



#adding leaves&galls columns to the RSEM file with calculating the average of FPKM&TPM in leaves&galls
def build_RSEM(ind):
        #print("index",ind)
        key = RSEM.iloc[ind, 0]
        insert_FPKM_TPM_values(RSEML4, "FPKM_L4", "TPM_L4", key, ind)
        insert_FPKM_TPM_values(RSEML5, "FPKM_L5", "TPM_L5", key, ind)
        insert_FPKM_TPM_values(RSEML6, "FPKM_L6", "TPM_L6", key, ind)
        insert_FPKM_TPM_values(RSEMG13, "FPKM_G13", "TPM_G13", key, ind)
        insert_FPKM_TPM_values(RSEMG14, "FPKM_G14", "TPM_G14", key, ind)
        insert_FPKM_TPM_values(RSEMG15, "FPKM_G15", "TPM_G15", key, ind)
        FPKMLeavesList = [RSEM.at[ind, "FPKM_L4"], RSEM.at[ind, "FPKM_L5"], RSEM.at[ind, "FPKM_L6"]]
        FPKMLeavesMiddle = middleOfThree(FPKMLeavesList)
        # print("middle= ", FPKMLeavesMiddle)
        RSEM.at[ind, "FPKM_L_AVG"] = Avg(FPKMLeavesList, FPKMLeavesMiddle, 20)
        FPKMGallsList = [RSEM.at[ind, "FPKM_G13"], RSEM.at[ind, "FPKM_G14"], RSEM.at[ind, "FPKM_G15"]]
        #print("AVG:",Avg(FPKMLeavesList, FPKMLeavesMiddle, 20))
        FPKMGallsMiddle = middleOfThree(FPKMGallsList)
        # print("middle= ", FPKMGallsMiddle)
        RSEM.at[ind, "FPKM_G_AVG"] = Avg(FPKMGallsList, FPKMGallsMiddle, 20)
        TPMLeavesList = [RSEM.at[ind, "TPM_L4"], RSEM.at[ind, "TPM_L5"], RSEM.at[ind, "TPM_L6"]]
        TPMLeavesMiddle = middleOfThree(TPMLeavesList)
        # print("middle= ", TPMLeavesMiddle)
        RSEM.at[ind, "TPM_L_AVG"] = Avg(TPMLeavesList, TPMLeavesMiddle, 20)
        TPMGallsList = [RSEM.at[ind, "TPM_G13"], RSEM.at[ind, "TPM_G14"], RSEM.at[ind, "TPM_G15"]]
        TPMGallsMiddle = middleOfThree(TPMGallsList)
        # print("middle= ", TPMGallsMiddle)
        RSEM.at[ind, "TPM_G_AVG"] = Avg(TPMGallsList, TPMGallsMiddle, 20)
        RSEM.at[ind, "Fold_increase_in_galls_FPKM"] = RSEM.at[ind, "FPKM_G_AVG"] / RSEM.at[ind, "FPKM_L_AVG"] if RSEM.at[ind, "FPKM_L_AVG"] > 0 else None
        RSEM.at[ind, "Fold_increase_in_galls_TPM"] = RSEM.at[ind, "TPM_G_AVG"] / RSEM.at[ind, "TPM_L_AVG"] if RSEM.at[ind, "TPM_L_AVG"] > 0 else None

# adding PFAM, blastp and blastx data to RSEM file
def add_PFAM_blastp_blastx_columns(ind,dataFile,dictFile,flag):
    similarEvalue = []
    key = RSEM.iloc[ind, 0]
    if key in dictFile.keys():
        print("flag",flag)
        print("key",key)
        df = dataFile[(dataFile['query name'] == key) & (dataFile['E-value'] == dictFile[key])]
        if dataFile.columns[0]=='target name':
            df = df.drop_duplicates(subset=['target name', 'accession', 'E-value'])
        if dataFile.columns[0]=='query name':
            df = df.drop_duplicates(subset=['accession', 'E-value'])
        if (len(df) > 1):
            print("2 evalue-s!!!")
            time.sleep(1) #test if everything is ok
            if flag=='f':
                RSEM.iloc[ind, RSEM.columns.get_loc('PFAM_accession')] = df.iloc[0]['accession']
                RSEM.iloc[ind, RSEM.columns.get_loc('PFAM_E-value')] = df.iloc[0]['E-value']
                RSEM.iloc[ind, RSEM.columns.get_loc('description of target')] = df.iloc[0]['description of target']
                RSEM.iloc[ind, RSEM.columns.get_loc('target name')] = df.iloc[0]['target name']
            if  flag=='x':
                RSEM.iloc[ind, RSEM.columns.get_loc('blastx_accession')] = df.iloc[0]['accession']
                RSEM.iloc[ind, RSEM.columns.get_loc('blastx_E-value')] = df.iloc[0]['E-value']
                RSEM.iloc[ind, RSEM.columns.get_loc('blastx_alignment length')] = df.iloc[0]['alignment length']
                RSEM.iloc[ind, RSEM.columns.get_loc('blastx_%identity')] = df.iloc[0]['%identity']
            if(flag=='p'):
                RSEM.iloc[ind, RSEM.columns.get_loc('blastp_accession')] = df.iloc[0]['accession']
                RSEM.iloc[ind, RSEM.columns.get_loc('blastp_E-value')] = df.iloc[0]['E-value']
                RSEM.iloc[ind, RSEM.columns.get_loc('blastp_alignment length')] = df.iloc[0]['alignment length']
                RSEM.iloc[ind, RSEM.columns.get_loc('blastp_%identity')] = df.iloc[0]['%identity']
            for i in range(len(df)):
                similarEvalue.append(key)
                similarEvalue.append(ind + 2)
                if flag == 'f':
                    similarEvalue.append(df.iloc[i]['target name'])
                    similarEvalue.append(df.iloc[i]['E-value'])
                if flag == 'x' or flag=='p':
                    similarEvalue.append(df.iloc[i]['accession'])
                    similarEvalue.append(df.iloc[i]['E-value'])
        if (len(df) == 1):
            if flag=='f':
                RSEM.iloc[ind, RSEM.columns.get_loc('PFAM_accession')] = df['accession']
                RSEM.iloc[ind, RSEM.columns.get_loc('PFAM_E-value')] = df['E-value']
                RSEM.iloc[ind, RSEM.columns.get_loc('description of target')] = df['description of target']
                RSEM.iloc[ind, RSEM.columns.get_loc('target name')] = df['target name']
            if  flag=='x':
                RSEM.iloc[ind, RSEM.columns.get_loc('blastx_accession')] = df['accession']
                RSEM.iloc[ind, RSEM.columns.get_loc('blastx_E-value')] = df['E-value']
                RSEM.iloc[ind, RSEM.columns.get_loc('blastx_alignment length')] = df['alignment length']
                RSEM.iloc[ind, RSEM.columns.get_loc('blastx_%identity')] = df['%identity']
            if(flag=='p'):
                RSEM.iloc[ind, RSEM.columns.get_loc('blastp_accession')] = df['accession']
                RSEM.iloc[ind, RSEM.columns.get_loc('blastp_E-value')] = df['E-value']
                RSEM.iloc[ind, RSEM.columns.get_loc('blastp_alignment length')] = df['alignment length']
                RSEM.iloc[ind, RSEM.columns.get_loc('blastp_%identity')] = df['%identity']
    return similarEvalue
# remove the last characters from query name in PFAM file & blastp file
def remove_sort(csvFile,flag,file_dict):
    #csvFile_dict={}
    for ind in range(len(csvFile)):
        #csvFile['query name'][ind]
        if flag=='f' or flag=='p':
            newkey = remove_last_characters(csvFile['query name'][ind])
            csvFile.at[ind, "query name"] = newkey
            print("index", ind)
            insert(file_dict, newkey, ind, csvFile)
        if flag=='x':
            insert(file_dict, csvFile['query name'][ind], ind, csvFile)
    return file_dict

def CSV_builder():
    SimilarEvalue_PFAM=[]
    SimilarEvalue_blastp=[]
    SimilarEvalue_blastx=[]
    for ind in range(len(RSEM)):
        print("RSEM index",ind)
        build_RSEM(ind) #adding 6 examples of leaves and galls
        print("add PFAM columns:\n")
        #building a dictionary of similar PFAM E-value
        SimilarEvalue_PFAM.extend(add_PFAM_blastp_blastx_columns(ind,PFAM,PFAM_dict,'f'))     #combine two lists
        print("add blastp columns:\n")
        #building a dictionary of similar blastp E-value
        SimilarEvalue_blastp.extend(add_PFAM_blastp_blastx_columns(ind,blastp, blastp_dict, 'p'))  # combine two lists
        print("add blastx columns:\n")
        #building a dictionary of similar blastx E-value
        SimilarEvalue_blastx.extend(add_PFAM_blastp_blastx_columns(ind,blastx, blastx_dict, 'x'))  # combine two lists
    RSEM.to_csv("RSEM.csv",index = False)
    print("list of similar Evalues in PFAM file:",SimilarEvalue_PFAM)
    print("list of similar Evalues in blastp file:", SimilarEvalue_blastp)
    print("list of similar Evalues in blastx file:", SimilarEvalue_blastx)



begin_time = datetime.datetime.now()
print("time now:",datetime.datetime.now())
##example for middle & average:
#LL=[0,5,9]
#middle=middleOfThree(LL)
#print("the middle of three numbers is: ",middle)
#print("Average= ",Avg(LL,middle,20))
#time.sleep(700) #Wait 700 seconds
RSEM = pd.read_csv(r'C:\master\Carmel\RSEM\Crmel RSEM all.csv',low_memory=False)
RSEML4=pd.read_csv(r'C:\master\Carmel\RSEM\leaves\RSEM_4all.genes.csv',low_memory=False)
RSEML5=pd.read_csv(r'C:\master\Carmel\RSEM\leaves\RSEM_5all.genes.csv',low_memory=False)
RSEML6=pd.read_csv(r'C:\master\Carmel\RSEM\leaves\RSEM_6.genes.csv',low_memory=False)
RSEMG13=pd.read_csv(r'C:\master\Carmel\RSEM\Galls\RSEM_13.genes.csv',low_memory=False)
RSEMG14=pd.read_csv(r'C:\master\Carmel\RSEM\Galls\RSEM_14.genes.csv',low_memory=False)
RSEMG15=pd.read_csv(r'C:\master\Carmel\RSEM\Galls\RSEM_CG_15all.genes.csv',low_memory=False)
PFAM = pd.read_csv(r'C:\master\PFAM Carmel all.csv',low_memory=False)
blastp=pd.read_csv(r'C:\master\Carmelall.swissprot.blastp.csv',low_memory=False)
blastx=pd.read_csv(r'C:\master\Carmelall.swissprot.blastx.csv',low_memory=False)
#print(PFAM.iloc[449476][6])
# print("PFAM header:",PFAM.columns[0])
# print("blastp header:",blastp.columns[0])
# print("blastx header:",blastx.columns[0])
# print("FPKM",RSEML4.at[0,'FPKM'])
#time.sleep(700) #Wait 700 seconds
print("end of loading!!\n")
#print("RSEM len:",len(RSEM))
print("time of loading:",datetime.datetime.now() - begin_time)
# adding new headers to RSEM file
heads = ["PFAM_E-value", "target name", "PFAM_accession", "description of target",
         "blastp_E-value", "blastp_accession", "blastp_%identity", "blastp_alignment length",
         "blastx_E-value", "blastx_accession", "blastx_%identity", "blastx_alignment length",
         "FPKM_L4", "FPKM_L5", "FPKM_L6", "FPKM_L_AVG", "Fold_increase_in_galls_FPKM",
         "FPKM_G_AVG", "FPKM_G13", "FPKM_G14", "FPKM_G15",
         "TPM_L4", "TPM_L5", "TPM_L6", "TPM_L_AVG", "Fold_increase_in_galls_TPM",
         "TPM_G_AVG", "TPM_G13", "TPM_G14", "TPM_G15"]
#adding new headers to RSEM file
for i in range(len(heads)):
    RSEM.insert(i + 1, heads[i], np.nan)
PFAM_dict = {}
blastp_dict={}
blastx_dict={}
#removing the last characters from query name in PFAM file & blastp file
print("the start time of insert:",datetime.datetime.now())
print(len(PFAM))
#removing the last characters from query name (".p1") in PFAM & blastp files + finding the min E-value for each query name in the 3 files(PFAM,blastp and blastx)
PFAM_dict=remove_sort(PFAM,'f',PFAM_dict)
blastp_dict=remove_sort(blastp,'p',blastp_dict)
blastx_dict=remove_sort(blastx,'x',blastx_dict)
print("time of insert PFAM:",datetime.datetime.now() - begin_time)
print("end of insert\n")
PFAM.to_csv("PFAM.csv",index = False)
print("CSV builder:\n")
#building the RSEM file with 6 examples for each of galls and leaves + adding the results of PFAM,blastp and blastx + calculating the average of leaves&galls in FPKM&TPM
CSV_builder()
print("time of running:",datetime.datetime.now() - begin_time)
