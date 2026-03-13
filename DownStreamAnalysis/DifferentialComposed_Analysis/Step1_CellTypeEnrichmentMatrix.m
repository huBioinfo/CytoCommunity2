clear
%%% Hyperparameters
num_TCN=10;

%Path
result_folder='../../';
image_file=fullfile(result_folder,'TNBC_Input', 'ImageNameList.txt');
fileID = fopen(image_file);
image_list=textscan(fileID,'%s','delimiter','\t','Headerlines',0);
fclose(fileID);
image_list=image_list{1};

fid0_path=fullfile(result_folder,'Step1_Output/UniqueCellTypeList.txt');
fid0=fopen(fid0_path);
CellTypeVec_List=textscan(fid0,'%s','delimiter','\t','Headerlines',0);
fclose(fid0);
CellTypeVec_List=CellTypeVec_List{1};

for image = 1:length(image_list)

    image_name=image_list{image};

    filename_1=strcat(result_folder, 'Step4_Output/ResultTable_File/ResultTable_', image_name,'.csv');
    fid1=fopen(filename_1);
    fmt_1 = ['%f','%f','%s','%d'];
    TargetGraph=textscan(fid1,fmt_1,'delimiter',',','Headerlines',1);
    fclose(fid1);
    TargetGraph_CellType=TargetGraph{3};
    TargetGraph_TCN=TargetGraph{4};
    
    CellTypeByTCN_matrix=ones(length(CellTypeVec_List),num_TCN); %initialization.
    for j=1:num_TCN
        idx=find(TargetGraph_TCN==j);
    
        if ~isempty(idx)
            component_CellType=TargetGraph_CellType(idx);
    
            for i=1:length(CellTypeVec_List)
                overlap_CellType_Cyto=length(find(strcmp(CellTypeVec_List{i},component_CellType)));
    
                if overlap_CellType_Cyto~=0
                    overlap_CellType_Image=length(find(strcmp(CellTypeVec_List{i},TargetGraph_CellType)));
                    enrich_p=1-hygecdf(overlap_CellType_Cyto-1,length(TargetGraph{1}),overlap_CellType_Image,length(idx));
                    CellTypeByTCN_matrix(i,j)=enrich_p;
                end
    
            end %end of i.
    
        end
    
    end %end of j.
    %(1)conduct BH-based multi-test correction.
    [h, crit_p, adj_ci_cvrg, CellTypeByTCN_matrix_pBH]=fdr_bh(CellTypeByTCN_matrix,0.05,'pdep','no');
    %(2)conduct log transformation.
    CellTypeByTCN_matrix_pBH(CellTypeByTCN_matrix_pBH<1.0E-20)=1.0E-20; %deal with the p-value of zero.   
    CellTypeByTCN_EnrichScoreMatrix=-log10(CellTypeByTCN_matrix_pBH);
    %CellTypeByTCN_EnrichScoreMatrix
    
    save CellTypeEnrichment_Res CellTypeByTCN_EnrichScoreMatrix CellTypeVec_List
    %writematrix(CellTypeByTCN_matrix,'CellTypeByTCN_matrix.csv');
    
    % 加载 .mat 文件
    data = load('CellTypeEnrichment_Res.mat');
    %data.CellTypeByTCN_EnrichScoreMatrix
    EnrichScoreMatrix='data/EnrichScoreMatrix/'
    output=[result_folder EnrichScoreMatrix]
    % 创建文件夹
    if ~exist(output, 'dir')
        mkdir(output); 
    end
    % 将 CellTypeByTCN_EnrichScoreMatrix 保存为csv文件
    output_folder = fullfile(output,[image_name '_EnrichScoreMatrix.csv']);
    writematrix(data.CellTypeByTCN_EnrichScoreMatrix, output_folder);
    
    % 将 CellTypeVec_List 保存为csv文件
    cellTypeVecTable = cell2table(data.CellTypeVec_List);
    CT_folder = fullfile(output,'CellTypeVec_List.csv');
    writetable(cellTypeVecTable, CT_folder);
end