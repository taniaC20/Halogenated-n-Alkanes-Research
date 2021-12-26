#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
from io import StringIO
import pandas as pd
import math
import sys
import re


# # Get .com files

# In[ ]:


#This is for the compounds that you are using just one route (like B3LYP)
def get_com_files_B3LYP():
    # Change the path "Gaussian\AMAP_Test\AMAP_TEST_coordinates.xyz" below, for the path where your .xyz file
    # (file with coordinates from openBabel) is.
    with open("Gaussian\AMAP_Test\AMAP_TEST_coordinates.xyz") as F:
    Lines = F.readlines()[2:]
    columns='Object    X    Y    Z\n'
    dfs = [columns]
    for line in Lines:
        if line.strip(): dfs[-1] += line
        elif dfs[-1]!=columns: dfs.append(columns)
        else: continue

    dfs = [pd.read_table(StringIO(df)) for df in dfs]

    total_df=len(dfs)
    pd.set_option('display.max_colwidth', -1)
    bar = progressbar.ProgressBar(max_value=total_df)

    for i in range(0,total_df):
        if i !=total_df-1:
            dfs[i].drop(dfs[i].tail(1).index,inplace=True)
            #Create a Directory like "Com_Files" below, where you will save all the .com files. Do this before running the code.
            # Change "Gaussian\AMAP_Test\Com_Files\" for the path where you want to save the .com files.
            dfs[i].to_csv(f"Gaussian\AMAP_Test\Com_Files\Compound{i+1}.com", header=None, index=None)
        else:
            # Change "Gaussian\AMAP_Test\Com_Files\" for the path where you want to save the .com files.
            dfs[i].to_csv(f"Gaussian\AMAP_Test\Com_Files\Compound{i+1}.com", header=None, index=None)
        bar.update(i)

    bar = progressbar.ProgressBar(max_value=197)
    for i in range(1,197):
        try:
            # Change "Gaussian\AMAP_Test\Com_Files\" for the path where you want to save the .com files.
            f = open(f'Gaussian\AMAP_Test\Com_Files\Compound{i}.com','r+')
        except IOError as e:
            print("File not found",i)
            continue
        # Below we add the directions in the .com file. If you need to change the route, charge or multiplicity, you can do so below.
        new_content=f'%nosave\n%mem=4GB\n%chk=Compound{i}.chk\n%nproc=16\n# B3LYP/6-311+G**\n\n Compound{i}\n\n1 2\n'
        lines = f.readlines()
        f.seek(0)
        f.write(new_content)
        for line in lines:
            f.write(line)
        f.write("\n")
        f.close()
        bar.update(i)


# In[ ]:


# This is for compounds that you are using three different routes (air,water,octanol)
# "solvent" -- is the parameter of the solvent you are using (air,water or octanol)
# "route" -- is the parameter for the route you are using.

#   Example: for solvent= "Air" you are gonna use route = "# M062X/6-311+G** opt freq"
#            for solvent= "Water" you are gonna use route = "# M062X/6-311+G** SCRF=(SMD, Solvent=Water) opt freq"
#            for solvent= "Octanol" you are gonna use route = "# M062X/6-311+G** SCRF=(SMD, Solvent=n-Octanol) opt freq"

def get_com_files_Water_Air_Octanol(solvent,route):
    # Change the path "Gaussian\PXAs\PXAs_coordinates.xyz" below, for the path where your .xyz file
    # (file with coordinates from openBabel) is.
    with open("Gaussian\PXAs\PXAs_coordinates.xyz") as F:
    Lines = F.readlines()[2:]
    columns='Object    X    Y    Z\n'
    dfs = [columns]
    for line in Lines:
        if line.strip(): dfs[-1] += line
        elif dfs[-1]!=columns: dfs.append(columns)
        else: continue

    dfs = [pd.read_table(StringIO(df)) for df in dfs]
    
    total_df=len(dfs)
    pd.set_option('display.max_colwidth', -1)
    bar = progressbar.ProgressBar(max_value=total_df)

    for i in range(0,total_df):
        if i !=total_df-1:
            dfs[i].drop(dfs[i].tail(1).index,inplace=True)
            #Create a Directory like "Com_Files" below, where you will save all the .com files. Do this before running the code.
            # Change "Gaussian\PXAs\Com_Files_" for the path where you want to save the .com files
            dfs[i].to_csv(f"Gaussian\PXAs\Com_Files_{solvent}\Compound{i+1}_{solvent}.com", header=None, index=None)
        else:
            # Change "Gaussian\PXAs\Com_Files_" for the path where you want to save the .com files
            dfs[i].to_csv(f"Gaussian\PXAs\Com_Files_{solvent}\Compound{i+1}_{solvent}.com", header=None, index=None)
        bar.update(i)
        
    bar = progressbar.ProgressBar(max_value=197)
    for i in range(1,197):
        try:
            # Change "Gaussian\PXAs\Com_Files_" for the path where you want to save the .com files
            f = open(f'Gaussian\PXAs\Com_Files_{solvent}\Compound{i}_{solvent}.com','r+')
        except IOError as e:
            print("File not found",i)
            continue
        # Below we add the directions in the .com file. If you need to change the route, charge or multiplicity, you can do so below.
        new_content=f'%nosave\n%mem=4GB\n%chk=Compound{i}_{solvent}.chk\n%nproc=16\n{route}\n\n Compound{i}_{solvent}\n\n1 2\n'
        lines = f.readlines()
        f.seek(0)
        f.write(new_content)
        for line in lines:
            f.write(line)
        f.write("\n")
        f.close()
        bar.update(i)
    


# # Extract data from .log files

# In[ ]:


#Function used below, always run before reading log files
def check(file):
    blabla="NImag=0"
    with open(file) as f:
        datafile = f.readlines()
    for line in datafile:
        if blabla in line:
            return True
    return False 


# In[ ]:


#Function used below, always run before reading log files for getting HF
def check_HF(file):
    blabla="HF="
    with open(file) as f:
        datafile = f.readlines()
    for line in datafile:
        if blabla in line:
            return True
    return False 


# In[ ]:


#This is for the compounds that you are using just one route (like B3LYP)
# For getting HF results
def read_log_files_B3LYP():
    final_df=pd.DataFrame(columns=["Compound","Data","Imag"])
    octanol_array=[]
    line="NImag=0"
    row=0
    for i in range(1,105):
        
        results=pd.DataFrame()
        found=False
        found_sum=False
        # Change 'Gaussian\DT_Test\Log Files\Compound{i}.log' for the path where your log files are.
        found=check(f'Gaussian\DT_Test\Log Files\Compound{i}.log')
        found_sum=check_HF(f'Gaussian\DT_Test\Log Files\Compound{i}.log')
        
        if found_sum==False:
            continue
        # Change 'Gaussian\DT_Test\Log Files\Compound{i}.log' for the path where your log files are.    
        myfile =f'Gaussian\DT_Test\Log Files\Compound{i}.log'
        lines_arr=[]
        
        with open(myfile) as f:
            previous = None
            current = next(f).strip()
            for line in f:
                line=line.strip()
                if "HF=" in current:
                    lines_arr.append(previous+current+line)
                previous = current
                current = line

        hf_line_arr=lines_arr[0].split('\\') 
        for elem in hf_line_arr:
            if "HF=" in elem:
                final_df.loc[row,'Compound']=i
                final_df.loc[row,'Data']=elem.replace('HF=','')      
        if found==True:
            final_df.loc[row,f"Imag"]=0
        else:
            final_df.loc[row,f"Imag"]=1
        row+=1
    # Change 'Gaussian\DT_Test\DT_TEST_B3LYP_HF_ENERGIES(-1-2).csv' for the path where you want to save the csv file.
    final_df.to_csv('Gaussian\DT_Test\DT_TEST_B3LYP_HF_ENERGIES(-1-2).csv', index=False)
    return final_df


# In[ ]:


# This is for compounds that you are using three diferent routes (air,water,octanol)
# For getting "Sum of electronic and thermal Free Energies=" results
def read_log_water_octanol_air()
    state_arr=["Air","Water","Octanol"]
    final_df=pd.DataFrame(columns=["Compound_Air","Description_Air","Data_Air","Compound_Water","Description_Water","Data_Water","Compound_Octanol","Description_Octanol","Data_Octanol"])
    octanol_array=[]
    line="NImag=0"
    for i in range(1,30):
        results=pd.DataFrame()
        found=False
        for elem in state_arr:
            found=False
            df_1 =pd.DataFrame(columns=[f"Compound_{elem}",f"Description_{elem}",f"Data_{elem}",f"Imag_{elem}"])

            try:
                # Change ''Gaussian\PAHs\Literature\Log Files\' for the path where your log files are.
                df = pd.read_csv(f'Gaussian\PAHs\Literature\Log Files\Compound{i}_{elem}.log', sep="\n", header=None)
            except IOError as e:
                print(i)
                print("File not found",elem)
                continue
            # Change ''Gaussian\PAHs\Literature\Log Files\' for the path where your log files are.
            found=check(f'Gaussian\PAHs\Literature\Log Files\Compound{i}_{elem}.log')
            df = df[df.iloc[:,0].str.contains("Sum of electronic and thermal Free Energies=", regex=True)]


            if len(df)==0:
                continue
            df = df.iloc[:,0].replace(" ","")
            df = pd.DataFrame(df.values.tolist())
            total_rows = len(df.index)
            substring='(Hartree/Particle)'
            for row in range(0,total_rows):
                if substring in df.loc[row,0]:
                    df.loc[row,0]=df.loc[row,0].replace('(Hartree/Particle)','')
            df_1[[f'Description_{elem}',f'Data_{elem}']]=df[0].str.split('=', 1, expand=True)
            df_1[[f'Compound_{elem}']]=i
            if found==True:
                df_1[[f"Imag_{elem}"]]=0
            else:
                df_1[[f"Imag_{elem}"]]=1
            results=pd.concat([results, df_1], axis=1)   
        final_df = pd.concat([final_df,results])
    final_df=final_df.drop(columns=['Description_Water', 'Description_Octanol',"Compound_Water","Compound_Octanol"])
    result=final_df.rename(columns={'Description_Air': 'Description', 'Compound_Air': 'Compound'})
    # Change 'Gaussian\PAHs\Literature\PAHs_Literature_M06_ENERGIES.csv' for the path where you want to save the csv file.
    result.to_csv('Gaussian\PAHs\Literature\PAHs_Literature_M06_ENERGIES.csv', index=False)
    result


# # Check Optimized Geometries

# In[ ]:


def extract_coordinates_from_log():
    solvents=["Air","Water","Octanol"]
    for x in range(1,101):
        for solvent in solvents:
            start = 0
            end = 0
            # Change "Gaussian\PAHs\Hundred\Log Files\" below for the path where your log files are.
            filename = f"Gaussian\PAHs\Hundred\Log Files\Compound{x}_{solvent}.log"
            # Change "Gaussian\PAHs\Hundred\Coordinates\" below for the path where you want to save the optimized coordinates
            # Create the Directory "Coordinates" before running the code
            newfile = f"Gaussian\PAHs\Hundred\Coordinates\Compound{x}_{solvent}.txt"

            openold = open(filename,"r")
            try: 
                    openold = open(filename,"r")
            except IOError as e:
                    continue
            # Change "Gaussian\PAHs\Hundred\Log Files\" below for the path where your log files are.
            df = pd.read_csv(f"Gaussian\PAHs\Hundred\Log Files\Compound{x}_{solvent}.log", sep="\n", header=None)        
            df = df[df.iloc[:,0].str.contains("Sum of electronic and thermal Free Energies=", regex=True)]
            if len(df)==0:
                continue
            opennew = open(newfile,"w")
            
            rline = openold.readlines()

            for i in range (len(rline)):
                if "Standard orientation:" in rline[i]:
                    start = i

            for m in range (start + 5, len(rline)):
                if "---" in rline[m]:
                    end = m
                    break

            # Conversion section

            for line in rline[start+5 : end] :
                words = line.split()
                word1 = int(words[1])
                word3 = str(words[3])
                if word1 == 17 :
                    word1 = "Cl"
                elif word1 == 9 :
                    word1 = "F "
                elif word1 == 35 :
                    word1 = "Br"
                elif word1 == 5 :
                    word1 = "B "
                elif word1 == 46 :
                    word1 = "Pd"
                elif word1 == 6 :
                    word1 = "C "
                elif word1 == 1:
                    word1 = "H "
                elif word1 == 7:
                    word1 = "N "
                elif word1 == 8:
                    word1 = "O "
                opennew.write( "%s%s" % (word1,line[30:-1])+"\n")
            opennew.close()


# In[ ]:


#This is used in the check_geometries function to check if the gaussian calculations for the compounds was completed for all air,water and octanol.
def check_completed(x):
    solvents=["Air","Water","Octanol"]
    not_completed=False

    for sol in solvents:
        try: 
            # Change "Gaussian\PAHs\Hundred\Coordinates\" below for the path where your optimized coordinates are.
            file=open(f'Gaussian\PAHs\Hundred\Coordinates\Compound{x}_{sol}.txt',"r")
        except IOError as e:
            not_completed=True
            break
            
    return not_completed


# In[ ]:


def check_geometries():
    solvents=["Air","Water","Octanol"]
    not_completed=[]
    not_equal=[]
    for x in range(1,102):
        not_completed=check_completed(x)
        if not_completed:
            continue
        sol_1="Air"
        
        # Change "Gaussian\PAHs\Hundred\Coordinates\" below for the path where your optimized coordinates are.
        df_1=pd.read_table(f'Gaussian\PAHs\Hundred\Coordinates\Compound{x}_{sol_1}.txt', delim_whitespace=True,names=['Element','X','Y','Z'])
        for sol_2 in solvents:
            df_difference=pd.DataFrame(columns=['Element','X','Y','Z'])
            if sol_1==sol_2:
                continue
            # Change "Gaussian\PAHs\Hundred\Coordinates\" below for the path where your optimized coordinates are.
            df_2=pd.read_table(f'Gaussian\PAHs\Hundred\Coordinates\Compound{x}_{sol_2}.txt', delim_whitespace=True,names=['Element','X','Y','Z'])

            for row in range(0,len(df_1)):
                df_difference.loc[row,'Element']=df_1.loc[row,'Element']
                df_difference.loc[row,'X']=("Close" if abs(df_1.loc[row,'X']- df_2.loc[row,'X']) < 1 else "Not Close")
                df_difference.loc[row,'Y']=("Close" if abs(df_1.loc[row,'Y']- df_2.loc[row,'Y']) < 1 else "Not Close")
                df_difference.loc[row,'Z']=("Close" if abs(df_1.loc[row,'Z']- df_2.loc[row,'Z']) < 1 else "Not Close")

            if "Not Close" in df_difference["X"].values or "Not Close" in df_difference["Y"].values or "Not Close" in df_difference["X"].values:
                not_equal.append([x,sol_2])           
    data=pd.DataFrame(not_equal, columns=["Compound","Solvent to change"])
    # Change "Gaussian\Middle Batch\List_All_Outliners.csv" for the path where you want to save the csv file.
    data.to_csv("Gaussian\Middle Batch\List_All_Outliners.csv", index=False)
    not_equal

