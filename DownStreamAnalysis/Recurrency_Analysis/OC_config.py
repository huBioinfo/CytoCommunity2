import math
import os
import pandas as pd

def Overlap_Coefficient(set1, set2):
    intersection = set1.intersection(set2)  
    min_size = min(len(set1), len(set2))  
    if min_size == 0:
        return 0  
    else:
        return len(intersection) / min_size

threshold = -math.log10(0.05)
Input_folder = 'data/EnrichScoreMatrix/'
Region_filename = "SKCM_Input/ImageNameList.txt"
Output_path ='data/OverlapCoefficient'
os.makedirs(Output_path,exist_ok=True)
Output_path= os.path.join(Output_path,"config.csv")

region_name_list = pd.read_csv(
        Region_filename,
        sep="\t",  # tab-separated
        header=None,  # no heading row
        names=["Image"],  # set our own names for the columns
    )

results = []
for TCN in range(0,10):
    for graph_A in range(0, len(region_name_list)): 
        region_A = region_name_list.Image[graph_A]
        df_A = pd.read_csv(Input_folder + region_A +'_EnrichScoreMatrix.csv',header=None)
        
        for graph_B in range(graph_A+1, len(region_name_list)):
            region_B = region_name_list.Image[graph_B]
            df_B = pd.read_csv(Input_folder + region_B +'_EnrichScoreMatrix.csv',header=None)
            
            A_TCN_vector = df_A.iloc[:,TCN]
            B_TCN_vector = df_B.iloc[:,TCN]
            
            A_EnrichCT = df_A[A_TCN_vector > threshold].index.tolist()
            B_EnrichCT = df_B[B_TCN_vector > threshold].index.tolist()

            set_A = set(A_EnrichCT)  
            set_B = set(B_EnrichCT)
            OC = Overlap_Coefficient(set_A, set_B)
            
            results.append({
                'TCN': "TCN"+ str(TCN+1),
                'Region_A': region_A,
                'Region_B': region_B,
                'OverlapCoefficient': OC
            })

results_df = pd.DataFrame(results)

results_df.to_csv(Output_path, index=False)