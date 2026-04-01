clear
close
clc

% MAGICAL_circuits=readtable('S16_S19/V2_Kidney_S16_S19_all_samples_trajectory_MAGICAL_circuits.txt');
% TF_binding=readtable('S16_S19/S16_S19_TF_freq_per_cluster.xlsx');

% MAGICAL_circuits=readtable('S9_S14/V2_Kidney_S9_S14_all_samples_trajectory_MAGICAL_circuits_fewer_S9.txt');
% TF_binding=readtable('S9_S14/S9_S14_TF_freq_per_cluster_fewer_S9.xlsx');

% MAGICAL_circuits=readtable('S9_S12/V2_Kidney_S9_S12_all_samples_trajectory_MAGICAL_circuits_fewer_S9.txt');
% no SNPs found

% MAGICAL_circuits=readtable('P15_P4/V2_Kidney_P15_P4_all_samples_trajectory_MAGICAL_circuits_V2.txt');
% TF_binding=readtable('P15_P4/P15_P4_TF_freq_per_cluster_V2.xlsx');

% MAGICAL_circuits=readtable('P15_P17/V2_Kidney_P15_P17_all_samples_trajectory_MAGICAL_circuits.txt');
% TF_binding=readtable('P15_P17/P15_P17_TF_freq_per_cluster.xlsx');

% MAGICAL_circuits=readtable('P15_P10/V2_Kidney_P15_P10_all_samples_trajectory_MAGICAL_circuits_V2.txt');
% TF_binding=readtable('P15_P10/P15_P10_TF_freq_per_cluster_V2.xlsx');

% MAGICAL_circuits=readtable('P15_P22/V2_Kidney_P15_P22_all_samples_trajectory_MAGICAL_circuits_V2.txt');
% TF_binding=readtable('P15_P22/P15_P22_TF_freq_per_cluster_V2.xlsx');


% MAGICAL_circuits=readtable('D17_D8/V2_Kidney_D17_D8_all_samples_trajectory_MAGICAL_circuits_V2.txt');
% TF_binding=readtable('D17_D8/D17_D8_TF_freq_per_cluster_V2.xlsx');

% MAGICAL_circuits=readtable('D17_D15/V2_Kidney_D17_D15_all_samples_trajectory_MAGICAL_circuits.txt');
% TF_binding=readtable('D17_D15/D17_D15_TF_freq_per_cluster.xlsx');

% MAGICAL_circuits=readtable('I22_I17/V2_Kidney_I22_I17_all_samples_trajectory_MAGICAL_circuits.txt');
% no SNPs found


% MAGICAL_circuits=readtable('I22_I18/V2_Kidney_I22_I18_all_samples_trajectory_MAGICAL_circuits.txt');
% TF_binding=readtable('I22_I18/I22_I18_TF_freq_per_cluster.xlsx');


MAGICAL_circuits=readtable('I14_I15/V2_Kidney_I14_I15_all_samples_trajectory_MAGICAL_circuits.txt');
TF_binding=readtable('I14_I15/I14_I15_TF_freq_per_cluster.xlsx');



%GWAS
[GWAS.chr, GWAS.point1, GWAS.point2, GWAS.SNP_ID, GWAS.gene]=textread('GWAS/GWAS_eGFR_hg38.txt', '%s %d %d %s %s');

GWAS_SNP_ID=cell(height(MAGICAL_circuits),1);
for g=1:424
    index1=find(strcmp(MAGICAL_circuits.Peak_chr, GWAS.chr(g))>0);
    index2=find(abs(GWAS.point1(g)+GWAS.point2(g)-MAGICAL_circuits.Peak_point1(index1)-MAGICAL_circuits.Peak_point2(index1))/2<1000);
    if ~isempty(index2)
        GWAS_SNP_ID(index1(index2))=GWAS.SNP_ID(g);
    end
end

GWAS_circuits_index=find(~cellfun(@isempty,GWAS_SNP_ID));

Motif_GWAS_count=zeros(length(TF_binding.Var1),1);
for m=1:length(TF_binding.Var1)
    index=find(contains(MAGICAL_circuits.TFs,[TF_binding.Var1(m),'_('])>0);
    index2=intersect(GWAS_circuits_index, index);
    if ~isempty(index2)
        Motif_GWAS_count(m,1)=length(unique(GWAS_SNP_ID(index2)));
    end
end



%GWAS
[GWAS.chr, GWAS.point1, GWAS.point2, GWAS.SNP_ID, GWAS.gene]=textread('GWAS/GWAS_eGFR_hg38.txt', '%s %d %d %s %s');

GWAS_SNP_ID=cell(height(MAGICAL_circuits),1);
for g=1:424
    index1=find(strcmp(MAGICAL_circuits.Peak_chr, GWAS.chr(g))>0);
    index2=find(abs(GWAS.point1(g)+GWAS.point2(g)-MAGICAL_circuits.Peak_point1(index1)-MAGICAL_circuits.Peak_point2(index1))/2<1000);
    if ~isempty(index2)
        GWAS_SNP_ID(index1(index2))=GWAS.SNP_ID(g);
    end
end

aa=~cellfun(@isempty,GWAS_SNP_ID);
Number_GWAS_SNPs=length(GWAS_SNP_ID(aa>0))


GWAS_locus_ID=cell(height(MAGICAL_circuits),1);
for g=425:length(GWAS.SNP_ID)
    index1=find(strcmp(MAGICAL_circuits.Peak_chr, GWAS.chr(g))>0);
    index2=find(GWAS.point1(g)<MAGICAL_circuits.Peak_point2(index1) & GWAS.point2(g)>MAGICAL_circuits.Peak_point1(index1));
    if ~isempty(index2)
        GWAS_locus_ID(index1(index2))=GWAS.SNP_ID(g);
    end
end
aa=~cellfun(@isempty,GWAS_locus_ID);
Number_GWAS_locus=length(unique(GWAS_locus_ID(aa>0)))


GWAS_genes_flag=zeros(height(MAGICAL_circuits),1);
for g=1:length(GWAS.SNP_ID)
    index1=find(strcmp(MAGICAL_circuits.Gene_symbol, GWAS.gene(g))>0);
    GWAS_genes_flag(index1)=1;
end
Number_GWAS_genes=length(unique(MAGICAL_circuits.Gene_symbol(GWAS_genes_flag>0)))



%validate TF-gene regulons
gene_list=unique(MAGICAL_circuits.Gene_symbol);
aa=~cellfun(@isempty,string(MAGICAL_circuits.CollecTRI_regulons));
mm=~strcmp(string(MAGICAL_circuits.CollecTRI_regulons), "Gene not included");
validated_genes=unique(MAGICAL_circuits.Gene_symbol(aa>0 & mm>0));
bb(1)=length(gene_list);
bb(2)=length(validated_genes);
[h,p_value]=fishertest([length(validated_genes),length(gene_list);6564,36588], 'tail', 'right')

length(validated_genes)/length(unique(MAGICAL_circuits.Gene_symbol(mm>0)))


%validated enhancer-gene looping
GeneHancer_elements=readtable('GeneHancer/GeneHancer_AnnotSV_elements_v5.19.txt', 'Delimiter','\t');
GeneHancer_elements.chr=string(GeneHancer_elements.chr);
GeneHancer_elements.chr='chr'+GeneHancer_elements.chr;
GeneHancer_geneassociations=readtable('GeneHancer/GeneHancer_AnnotSV_gene_association_scores_v5.19.txt', 'Delimiter','\t');

Looping_flag=zeros(height(MAGICAL_circuits),2);
for p=1:height(MAGICAL_circuits)
    index1=find(strcmp(GeneHancer_elements.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(GeneHancer_elements.element_start(index1)+GeneHancer_elements.element_end(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<500);
    if ~isempty(index2)
        Looping_flag(p,1)=1;
        index3=find(ismember(GeneHancer_geneassociations.GHid, unique(GeneHancer_elements.GHid(index1(index2))))>0);
        aa = intersect(GeneHancer_geneassociations.symbol(index3),MAGICAL_circuits.Gene_symbol(p));
        if ~isempty(aa)
            Looping_flag(p,2)=1;
        end
    end
end
sum(Looping_flag)
length(unique(MAGICAL_circuits.Peak_ID(Looping_flag(:,1)>0)))/length(unique(MAGICAL_circuits.Peak_ID(Looping_flag(:,1)>=0)))
length(unique(MAGICAL_circuits.Gene_symbol(Looping_flag(:,2)>0)))/length(intersect(MAGICAL_circuits.Gene_symbol,GeneHancer_geneassociations.symbol))


%CUT&RUN
[H3K27ac_peaks.chr, H3K27ac_peaks.point1, H3K27ac_peaks.point2]=textread('CUT_RUN-selected/Merge/Merge_biopsy/H3K27ac_broadPeak', '%s %d %d');
[H3K27me3_peaks.chr, H3K27me3_peaks.point1, H3K27me3_peaks.point2]=textread('CUT_RUN-selected/Merge/Merge_biopsy/H3K27me3_broadPeak', '%s %d %d');
[H3K4me1_peaks.chr, H3K4me1_peaks.point1, H3K4me1_peaks.point2]=textread('CUT_RUN-selected/Merge/Merge_biopsy/H3K4me1_broadPeak', '%s %d %d');
[H3K4me3_peaks.chr, H3K4me3_peaks.point1, H3K4me3_peaks.point2]=textread('CUT_RUN-selected/Merge/Merge_biopsy/H3K4me3_narrowPeak', '%s %d %d');

H3K27ac_peaks_flag=zeros(height(MAGICAL_circuits),1);
H3K27me3_peaks_flag=zeros(height(MAGICAL_circuits),1);
H3K4me1_peaks_flag=zeros(height(MAGICAL_circuits),1);
H3K4me3_peaks_flag=zeros(height(MAGICAL_circuits),1);

for p=1:height(MAGICAL_circuits)
    index1=find(strcmp(H3K27ac_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(H3K27ac_peaks.point1(index1)+H3K27ac_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        H3K27ac_peaks_flag(p)=1;
    end

    index1=find(strcmp(H3K27me3_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(H3K27me3_peaks.point1(index1)+H3K27me3_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        H3K27me3_peaks_flag(p)=1;
    end

    index1=find(strcmp(H3K4me1_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(H3K4me1_peaks.point1(index1)+H3K4me1_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        H3K4me1_peaks_flag(p)=1;
    end

    index1=find(strcmp(H3K4me3_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(H3K4me3_peaks.point1(index1)+H3K4me3_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        H3K4me3_peaks_flag(p)=1;
    end
end

CUT_Run_flag=[H3K27ac_peaks_flag, H3K4me1_peaks_flag, H3K4me3_peaks_flag, H3K27me3_peaks_flag];
sum(CUT_Run_flag)/size(CUT_Run_flag,1)

length(find(sum(CUT_Run_flag(:,1:3),2)>0))/size(CUT_Run_flag,1)
