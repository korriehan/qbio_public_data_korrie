#!/usr/bin/env python
# coding: utf-8


### Creation of BoxPlots for ESR1 Protien Expression ###

import cptac  # import cptac to download cptac protein and clinical data
import pandas as pd  # import pandas for
import matplotlib.pyplot as plt  # import matplotlib for boxplot
import seaborn as sns  # import seaborn for boxplot
from scipy import stats
from statannot import add_stat_annotation

cptac.download(dataset="Brca")   # download breast cancer dataset
br = cptac.Brca()  # save data in br variable

protein_data = br.get_proteomics()  # save proteomic data
protein_data = protein_data.droplevel(1, axis=1)     # remove multi index
clinical_data = br.get_clinical()  # save clinical data

esr1 = protein_data["ESR1"]  # save ESR1 protein expression column
clinical_data["ER.IHC.Score"] = clinical_data["ER.IHC.Score"].fillna("Not reported")   # fill in null values

er_mask = clinical_data["ER.IHC.Score"] == "3+"  # 3+ is ER-positive, create mask of ER-positive patients
patients = esr1[er_mask]  # apply mask to protein expression data
ages = clinical_data["Age.in.Month"][er_mask]/12  # calculate ages in years
er_positive_patients = pd.DataFrame(patients, columns=["ESR1"])  # create new dataframe
er_positive_patients["Age"] = ages  # apply ages to new column in dataframe
category = []  # set categories list
for age in er_positive_patients["Age"]:  # for each age in the age column
    if age < 50:  # if it's less than 50, append young/patient is considered young
        category.append("Young")
    else:  # if more than/equal to 50, append old/patient is considered old
        category.append("Old")
category = pd.array(category)  # make category a pandas array
er_positive_patients["Age_category"] = category   # create new column in df titled age category
plt.figure()
df = er_positive_patients  # set df to data frame
x = "Age_category"   # set x to age_category
y = "ESR1"  # set y to esr1 counts
order = ["Young", "Old"]  # order of x axis
ax = sns.boxplot(data=df, x=x, y=y, order=order)  # create boxplot
add_stat_annotation(ax, data=df, x=x, y=y, order=order,
                                   box_pairs=[("Young", "Old")],
                                   test='Mann-Whitney', text_format='star',
                                   loc='outside', verbose=2)     # calculate p-value
plt.savefig("/PATHWAY/esr1_protein_expression_young_old_boxplot.png")

# calculate regression
slope, intercept, r_value, p_value, std_err = stats.linregress(er_positive_patients["Age"], er_positive_patients["ESR1"])
plt.figure()
# create scatterplot w/regression line
ax = sns.regplot(x="Age", y="ESR1", data=er_positive_patients, color='b', line_kws={'label': "y={0:.1f}x+{1:.1f}".format(slope, intercept)})
# plot legend
ax.legend()
plt.savefig("/PATHWAY/esr1_protein_expression_young_old_linear.png")




### Creation of Spearman Plots ###
# In[17]:


get_ipython().system('jupyter nbconvert --to script Spearman_ESR1.ipynb')


# In[1]:


import cptac
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats


# In[2]:


cptac.download(dataset="Brca")
br = cptac.Brca()

protein_data = br.get_proteomics()

#The dataframes are MultIndex pandas dataframes.
#However, to teach the basics of pandas, we will remove the "multi" part of the dataframe.
protein_data = protein_data.droplevel(1, axis=1)

rna_data = br.get_transcriptomics()
clinical_data = br.get_clinical()


# In[3]:


clinical_data["Age_in_years"] = clinical_data["Age.in.Month"]/12


# In[4]:


assert list(rna_data.index) == list(protein_data.index)


# In[5]:


rna_esr1 = rna_data.loc[: , "ESR1"]
protein_esr1 = protein_data.loc[: , "ESR1"]

rho, spear_pvalue = stats.spearmanr( rna_esr1, protein_esr1 )

rho_check, spear_pvalue_check = stats.spearmanr( protein_esr1, rna_esr1 )

assert rho == rho_check


# In[16]:


plt.figure( figsize=(10,10) )

m, b = np.polyfit(protein_esr1, rna_esr1, 1)
print(m)
print(b)
#Replace x and y with appropriate variables
plt.scatter( protein_esr1, rna_esr1, c='black')
plt.plot(protein_esr1, m*protein_esr1 + b, 'g')

title = "rho: {} for ESR1 (all ages)".format(rho) #This is string formatting. The variable in the () will print in the {}
plt.title(title)

#Fill in informative x and y labels
plt.xlabel("Proteomic Data")
plt.ylabel("RNA Data")

#plt.show() #Comment out when running in script
plt.savefig( "/Users/Christopher/Desktop/Datanalysis/qbio_data_analysis_Chris/final_proj/SpearGraphALL.png", bbox_inches="tight" )


# In[11]:


#What column of clinical_data is referring to age?
young_mask = clinical_data["Age_in_years"] < 50.0
old_mask = clinical_data["Age_in_years"] >= 50.0

#Check for understanding: Why do the below lines work?
rna_esr1_young = rna_data["ESR1"][ young_mask ]
protein_esr1_young = protein_data["ESR1"][ young_mask ]

#We want all patients of the ESR1 column
rna_esr1_old = rna_data["ESR1"][ old_mask ]
protein_esr1_old = protein_data["ESR1"][ old_mask ]


# In[14]:


#YOUNG PLOT
rho_young, spear_pvalue_young = stats.spearmanr( rna_esr1_young, protein_esr1_young )

plt.figure( figsize=(10,10) )
m, b = np.polyfit(protein_esr1_young, rna_esr1_young, 1)
print(m)
print(b)
#Replace x and y with appropriate variables
plt.scatter( protein_esr1_young, rna_esr1_young, c='black' )
plt.plot(protein_esr1_young, m*protein_esr1_young + b, 'g')

title = "rho: {} for ESR1 (Patients < 50 years old)".format(rho_young) #This is string formatting. The variable in the () will print in the {}
plt.title(title)

plt.xlabel("Young Proteomic Data")
plt.ylabel("Young RNA Data")

#plt.show() #Comment out when running in script
plt.savefig( "/Users/Christopher/Desktop/Datanalysis/qbio_data_analysis_Chris/final_proj/SpearGraphYOUNG.png", bbox_inches="tight" )


# In[15]:


#OLD PLOT
rho_old, spear_pvalue_old = stats.spearmanr( rna_esr1_old, protein_esr1_old )

plt.figure( figsize=(10,10) )
m, b = np.polyfit(protein_esr1_old, rna_esr1_old, 1)
print(m)
print(b)
#Replace x and y with appropriate variables
plt.scatter( protein_esr1_old, rna_esr1_old, c='black' )
plt.plot(protein_esr1_old, m*protein_esr1_old + b, 'g')

title = "rho: {} for ESR1 (Patients >= 50 years old)".format(rho_old) #This is string formatting. The variable in the () will print in the {}
plt.title(title)

plt.xlabel("Old Proteomic Data")
plt.ylabel("Old RNA Data")

#plt.show() #Comment out when running in script
plt.savefig( "/Users/Christopher/Desktop/Datanalysis/qbio_data_analysis_Chris/final_proj/SpearGraphOLD.png", bbox_inches="tight" )


# In[ ]:


def spear_rho_plot(rna, protein, genename, figsz=10, pathout):
    rho = stats.spearmanr(rna, protein)
    m, b = np.polyfit(protein, rna, 1)
    plt.figure(figsize=(figsz, figsz))

    plt.scatter(protein, rna, c="black")
    plt.plot(protein, m*protein + b, 'g')

    title = "rho: {0} for {1}".format(rho_old, genename) #This is string formatting. The variable in the () will print in the {}
    plt.title(title)
    plt.xlabel("Proteomic Data")
    plt.ylabel("RNA Data")

    plt.savefig(pathout, bbox_inches="tight" )

### Spearman Plot Function ###
    #Spearman Rho Correlation Plotting Function for Specific Gene
def spear_rho_plot(rna, protein, genename, figsz=10, pathout):

    rho = stats.spearmanr(rna, protein)
    #Calculates spearman rho to use given rna and proteomic data for a specific gene
    m, b = np.polyfit(protein, rna, 1)
    #Calculates slope of data points in the correlation

    plt.figure(figsize=(figsz, figsz))

    #Plots the data points with a green trendline going through it
    plt.scatter(protein, rna, c="black")
    plt.plot(protein, m*protein + b, 'g')

    #Graph labels
    title = "rho: {0} for {1}".format(rho_old, genename) #Title of graph using string formatting
    plt.title(title)
    plt.xlabel("Proteomic Data")
    plt.ylabel("RNA Data")

    #Saves figure to the path given in function if running as python script
    plt.savefig(pathout, bbox_inches="tight" )
