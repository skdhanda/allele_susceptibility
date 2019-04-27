#!/usr/bin/python
from __future__ import division
import pandas as pd
import sys,os
import numpy as np
import scipy.stats as stats
from scipy.stats import chi2_contingency
import argparse

def OR_CI_abcd(a,b,c,d):
   print a,b,c,d
   chi2, pvalue_chi2, dof, ex= stats.chi2_contingency([[a,c],[b,d]],correction=False)
   chi2_corrected, pvalue_corrected, dof_corrected, ex_corrected= stats.chi2_contingency([[a,c],[b,d]])
   odd_ratio_fisher,pvalue_fisher=stats.fisher_exact([[a,c],[b,d]])
   if a==0 or b==0 or c==0 or d==0: 
     a+=0.5
     b+=0.5
     c+=0.5
     d+=0.5
   odd_ratio=(a*d)/(b*c)
   se=np.sqrt(1/a+1/b+1/c+1/d)
   ci_l=np.exp(np.log(odd_ratio)-1.96*se)
   ci_h=np.exp(np.log(odd_ratio)+1.96*se)
   return odd_ratio,se,ci_l,ci_h,pvalue_fisher,chi2,pvalue_chi2,chi2_corrected,pvalue_corrected

def OR_CI(row):
   a=row['Patient Count']
   b=row['Healthy Count']
   c=row['Remaining Patients']
   d=row['Remaining Healthy']
   ###print a,b,c,d
   odd_ratio_fisher,pvalue_fisher=stats.fisher_exact([[a,c],[b,d]])
   if a==0 or b==0 or c==0 or d==0: 
     a+=0.5
     b+=0.5
     c+=0.5
     d+=0.5
   odd_ratio=(a*d)/(b*c)
   chi2, pvalue_chi2, dof, ex= stats.chi2_contingency([[a,c],[b,d]],correction=False)
   chi2_corrected, pvalue_corrected, dof_corrected, ex_corrected= stats.chi2_contingency([[a,c],[b,d]])
   se=np.sqrt(1/a+1/b+1/c+1/d)
   ci_l=np.exp(np.log(odd_ratio)-1.96*se)
   ci_h=np.exp(np.log(odd_ratio)+1.96*se)
   return odd_ratio,se,ci_l,ci_h,pvalue_fisher,chi2,pvalue_chi2,chi2_corrected,pvalue_corrected


parser = argparse.ArgumentParser(description='Uses disease and control data in tsv files and the output file to print the resutl')
parser.add_argument('-p','--patient_data', dest='pfile',type=argparse.FileType('r'),help='The file containing the patient data along with their alleles')
parser.add_argument('-c', '--control_data', dest='cfile',type=argparse.FileType('r'),help='The file containing the data from control/healthy individuals along with their alleles')
parser.add_argument('-o','--output', dest='ofile',type=argparse.FileType('w'),default='output.tsv',help='The output file to save the results in a file, default is output.tsv')
args = parser.parse_args()

dis=pd.read_csv(args.pfile,sep="\t")
con=pd.read_csv(args.cfile,sep="\t")



def run_param(dis,con):
  parameters=['A_ag','B_ag','C_ag','DR_ag','DQ_ag']
  final_pd=pd.DataFrame()
  for param in parameters:
    dis_list=[]
    dis[param+"_1"].replace(" ","",inplace=True,regex=True)
    con[param+"_1"].replace(" ","",inplace=True,regex=True)
    dis[param+"_2"].replace(" ","",inplace=True,regex=True)
    con[param+"_2"].replace(" ","",inplace=True,regex=True)
    
    dis[param+"_2"].replace("-",dis[param+"_1"],inplace=True)
    con[param+"_2"].replace("-",con[param+"_1"],inplace=True)
    dis[param+"_1"].replace("-",dis[param+"_2"],inplace=True)
    con[param+"_1"].replace("-",con[param+"_2"],inplace=True)
    for i in ["_1","_2"]:
      dis[param+i].fillna('Blank',inplace=True)
      con[param+i].fillna('Blank',inplace=True)
      dis_list.extend(dis[param+i].unique().tolist())
      dis_list.extend(con[param+i].unique().tolist())
    total_list=list(set(dis_list))

    data_pd=pd.DataFrame()
    data_pd['Patient Count']=0
    data_pd['Healthy Count']=0
    for allele in total_list:
      data_pd.at[allele,'allele']=allele
      data_pd.at[allele,'param']=param
      data_pd.at[allele,'Patient Count'] = len(dis[(dis[param+"_1"]==allele) | (dis[param+"_2"]==allele)])
      data_pd.at[allele,'Healthy Count'] = len(con[(con[param+"_1"]==allele) | (con[param+"_2"]==allele)])
    dsum=data_pd['Patient Count'].sum()
    csum=data_pd['Healthy Count'].sum()
    data_pd['Remaining Patients']=dsum-data_pd['Patient Count']
    data_pd['Remaining Healthy']=csum-data_pd['Healthy Count']
    data_pd['Patient frequency (%)']=100*data_pd['Patient Count']/dsum
    data_pd['Healthy frequency (%)']=100*data_pd['Healthy Count']/csum
    data_pd=data_pd.reset_index()
    data_pd[['odds ratio','se','ci lower','ci higher','pvalue_fisher','chi2','pvalue_chi2','chi2_corrected','pvalue_corrected']]=data_pd.apply(OR_CI,axis=1,result_type="expand")
    final_pd=final_pd.append(data_pd,sort=True)
  #dis.to_csv('cml_mod.tsv',sep="\t")
  #con.to_csv('control_mod.tsv',sep="\t")
  return final_pd

df_output=run_param(dis,con)
df_output['Category']='Overall'
#### Gender wise analysis
dis['Sex'].replace(" ","",inplace=True,regex=True)
con['Sex'].replace(" ","",inplace=True,regex=True)
print "Different genders present in patient data  ", dis['Sex'].unique()
print "Different genders present in heathy individuals ", con['Sex'].unique()
#df_gender=pd.DataFrame()
print con['Age'].unique()
print con['Age'].unique()
for gender in ['M','F']:
  dis_gender=dis.loc[dis['Sex']==gender].copy()
  con_gender=con.loc[con['Sex']==gender].copy()
  print("number of patients with gender", gender,len(dis_gender))
  print("number of healthy individuals with gender", gender,len(con_gender))
  df_gen=run_param(dis_gender,con_gender)
  df_gen['Category']=gender
  df_output=df_output.append(df_gen,ignore_index=True,sort=True)
dis_younger=dis.loc[dis['Age']<=18].copy()
con_younger=con.loc[con['Age']<=18].copy()
print("number of patients with age group of  <=18 years",len(dis_younger))
print("number of healthy individual with age group of <=18years",len(con_younger))
#print con_younger,dis_younger
df_gen=run_param(dis_younger,con_younger)
df_gen['Category']='Younger <=18 years'
df_output=df_output.append(df_gen,ignore_index=True,sort=True)
dis_adults=dis.loc[dis['Age']>18].copy()
con_adults=con.loc[con['Age']>18].copy()
print("number of patients with age group of  >18 years",len(dis_adults))
print("number of healthy individual with age group of >18years",len(con_adults))
df_gen=run_param(dis_adults,con_adults)
df_gen['Category']='Adults >18 years'
df_output=df_output.append(df_gen,ignore_index=True,sort=True)
df_output.to_csv(args.ofile,index=False)
