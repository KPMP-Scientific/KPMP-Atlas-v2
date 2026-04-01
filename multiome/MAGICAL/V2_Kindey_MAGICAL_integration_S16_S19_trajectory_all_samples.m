% /Applications/MATLAB_R2019a.app/bin/matlab -batch Demo_main

clear
close all hidden
clc
warning('off','all')

%filter out samples with too few cells in this analysis
Sample_meta=readtable('Adventitial_Fibroblasts_pseudobulk/Samples_Conditions_and_cell_numbers_S16_S19.csv');
Sample_meta=Sample_meta(Sample_meta.Freq>50,:);

%load peaks called from all cell types
[Peak_ID, Peaks(:,1), Peaks(:,2), Peaks(:,3)]=textread('V2_peaks.csv', '%s %d %d %d', 'headerlines', 1, 'delimiter', ',');    
TF_binding=load('V2_Peak_Motifs_mapping.mat');


fprintf('load scATAC differential statistics from R output\n\n')
% select significant peaks and genes for all selected clusters
Peak_sc_diff_all = readtable('Adventitial_Fibroblasts_pseudobulk/All_samples_S16_S19_ATAC_markers.csv');
DA_peak_sc_index=find(abs(Peak_sc_diff_all.avg_log2FC)>0.3 & Peak_sc_diff_all.p_val<0.05 & (Peak_sc_diff_all.pct_1>0.1 | Peak_sc_diff_all.pct_2>0.1));
Peak_sc_diff_all = Peak_sc_diff_all(DA_peak_sc_index,:);


fprintf('load scRNA differential statistics from R output\n\n')
RNA_sc_diff_all = readtable('Adventitial_Fibroblasts_pseudobulk/All_samples_S16_S19_RNA_markers.csv');
DE_gene_sc_index = find(abs(RNA_sc_diff_all.avg_log2FC)>0.5 & RNA_sc_diff_all.p_val_adj<0.05 & max([RNA_sc_diff_all.pct_1, RNA_sc_diff_all.pct_2], [], 2)>0.1);
RNA_sc_diff_all = RNA_sc_diff_all(DE_gene_sc_index,:);

%RNA differential calling is wrong as there is no data normalization
%performed and not in log2 space

Candidate_clusters=intersect(unique(Peak_sc_diff_all.cluster), unique(RNA_sc_diff_all.cluster));
for p=1:length(Candidate_clusters)
    peak_index=find(strcmp(Peak_sc_diff_all.cluster, Candidate_clusters{p})>0);
    Peak_sc_diff{p}=Peak_sc_diff_all(peak_index,:);
    Peak_diff_ID{p}=Peak_sc_diff{p}.gene;
    fprintf(2, '%d DA peaks from single cell differential analysis of cluster %s\n\n', length(Peak_diff_ID{p}), Candidate_clusters{p})


    RNA_index=find(strcmp(RNA_sc_diff_all.cluster, Candidate_clusters{p})>0);
    RNA_sc_diff{p}=RNA_sc_diff_all(RNA_index,:);
    Gene_diff_symbol{p}=RNA_sc_diff{p}.gene;
    fprintf(2, '%d DE genes from single cell differential analysis of cluser %s\n\n', length(Gene_diff_symbol{p}), Candidate_clusters{p})
end



%% load ATAC and RNA pseudo bulk data from all clusters since we will regress the whole trajectory now
for p=1:length(Candidate_clusters) 
    ATAC = readtable(['Adventitial_Fibroblasts_pseudobulk/Pseudo_bulk_S16_S19_count_',Candidate_clusters{p},'_All_samples.txt'], 'ReadVariableNames',true);
    ATAC_samples = ATAC.Properties.VariableDescriptions(2:end);
    ATAC_count = sparse(ATAC{:,2:end});
    [ATAC_samples, sample_index] = intersect(ATAC_samples, Sample_meta.Var1);
    ATAC_count = ATAC_count(:,sample_index);


    RNA = readtable(['Adventitial_Fibroblasts_pseudobulk/Pseudo_bulk_S16_S19_RNA_count_',Candidate_clusters{p},'_All_samples.txt'], 'ReadVariableNames',true);
    RNA_symbols = RNA.Var1;
    RNA_samples = RNA.Properties.VariableDescriptions(2:end);
    RNA_count = sparse(RNA{:,2:end});
    [RNA_samples, sample_index] = intersect(RNA_samples, Sample_meta.Var1);
    RNA_count = RNA_count(:,sample_index);

    if p==1
        ATAC_samples_whole = insertBefore(ATAC_samples, 1, [Candidate_clusters{p},'_']);
        ATAC_count_whole = ATAC_count;

        RNA_samples_whole = insertBefore(RNA_samples, 1, [Candidate_clusters{p},'_']);
        RNA_count_whole = RNA_count;
    else
        ATAC_samples_whole = [ATAC_samples_whole;insertBefore(ATAC_samples, 1, [Candidate_clusters{p},'_'])];
        ATAC_count_whole = [ATAC_count_whole,ATAC_count];

        RNA_samples_whole = [RNA_samples_whole;insertBefore(RNA_samples, 1, [Candidate_clusters{p},'_'])];
        RNA_count_whole = [RNA_count_whole,RNA_count];
    end
end


total_ATAC_reads=5e6;
ATAC_raw_count_sum=sum(ATAC_count_whole)+1;
for s=1:size(ATAC_count_whole,2)
    ATAC_count_whole(:,s)=ATAC_count_whole(:,s)*(total_ATAC_reads/ATAC_raw_count_sum(s));
end

total_RNA_reads=5e6;
RNA_raw_count_sum=sum(RNA_count_whole)+1;
for s=1:size(RNA_count,2)
    RNA_count_whole(:,s)=RNA_count_whole(:,s)*(total_RNA_reads/RNA_raw_count_sum(s));
end

% RNA_count_whole_flag=RNA_count_whole~=0;
% RNA_sample_index=find(sum(RNA_count_whole_flag)/size(RNA_count_whole_flag,1)>0.05);
% RNA_samples_whole=RNA_samples_whole(RNA_sample_index);
% RNA_count_whole=RNA_count_whole(:,RNA_sample_index);
% 
% ATAC_count_whole_flag=ATAC_count_whole~=0;
% ATAC_sample_index=find(sum(ATAC_count_whole_flag)/size(ATAC_count_whole_flag,1)>0.05);
% ATAC_samples_whole=ATAC_samples_whole(ATAC_sample_index);
% ATAC_count_whole=ATAC_count_whole(:,ATAC_sample_index);

[Common_samples, bidex, cidex]=intersect(ATAC_samples_whole, RNA_samples_whole);
ATAC_count_whole=ATAC_count_whole(:,bidex);
RNA_count_whole=RNA_count_whole(:,cidex);

%% MAGICAl analysis
TF_flag=zeros(1,870);
for p=1:length(Candidate_clusters)
    p
    Candidate_Gene_Symbol = Gene_diff_symbol{p};
    Candidate_Peak_ID = Peak_diff_ID{p};
    

    %Enrichment analysis on differential peaks (with active peaks as background), to select candidate TFs regulating these peaks
    fprintf('Select candidate TFs that are enriched in differential peaks, with active peaks as background\n\n')
    [Candidate_Peak_ID, bidex]=intersect(TF_binding.Peak_ID, Candidate_Peak_ID);
    Motif_diff_peak_binding=TF_binding.Peak_motif_mapping(bidex,:);

    TF_pct(p,:)=sum(Motif_diff_peak_binding)/length(Candidate_Peak_ID);
    TF_background_pct=sum(TF_binding.Peak_motif_mapping)/length(TF_binding.Peak_ID);
    TF_enrichment_FC(p,:)=TF_pct(p,:)./TF_background_pct;
    TF_index=find(TF_pct(p,:)>0.10 & log2(TF_enrichment_FC(p,:))>1);
    TF_flag(p, TF_index)=1;

    if ~isempty(TF_index)
        Candidate_TFs=TF_binding.Motifs(TF_index);
        Candidate_TF_Peak_Binding=Motif_diff_peak_binding(:,TF_index);
        Binding_peak_index=find(sum(Candidate_TF_Peak_Binding, 2)>0);
        Candidate_Peak_ID=Candidate_Peak_ID(Binding_peak_index);
        Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding(Binding_peak_index,:);

        [Candidate_Peak_ID,bidex,cidex]=intersect(Candidate_Peak_ID, Peak_ID, 'stable');
        Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding(bidex,:);
        Candidate_Peaks=Peaks(cidex,:);
        fprintf(2, '%d enriched TFs\n\n', length(Candidate_TFs))
    else
        fprintf('Too few peaks with TF binding sites.\nMAGICAL not applicable to this cell type!\n\n\n\n')
        %         continue;
    end

    %% MAGICAL integration
    fprintf('MAGICAL integration starts!\n\n')
    %     MAGICAL_post=MAGICAL_V1(Candidate_Peak_ID, Candidate_Peaks, ATAC, Candidate_Gene_Symbol, RNA, Candidate_TFs, Candidate_TF_Peak_Binding, Sample_meta);
    %************************* Input data *********************
    % Candidate_Peak_ID: P x 1 string vector
    % Candidate_Peaks: P x 3 matrix with chr, point1, point2
    % ATAC: P x S matrix, raw pseudo bulk ATAC count for all peaks and S samples
    % Candidate_Gene_Symbol: G x 1 string vector with gene names
    % RNA: G x S matrix, raw pseudo bulk RNA count for all genes and S samples
    % Candidate_TFs: T x 1 string vector for enriched TFs
    % Candidate_TF_Peak_Binding: P x T binary matrix, binding state for P peaks and T TFs

    %********** Initial Integration to select TFs, peaks and genes ***********
    [Candidate_TFs, Candidate_TF_log2Count,...
        Candidate_Peak_ID, Candidate_Peaks, Candidate_Peak_log2Count,...
        Candidate_Gene_Symbol, Candidate_Gene_TSS, Candidate_Gene_log2Count,...
        Candidate_TF_Peak_Binding, Candidate_Peak_Gene_looping, Selected_samples]= Initial_integration(Candidate_Peak_ID, Candidate_Peaks, Peak_ID, ATAC_count_whole,...
        Candidate_Gene_Symbol, RNA_symbols,  RNA_count_whole, Candidate_TFs, Candidate_TF_Peak_Binding, Common_samples);

    S=length(Selected_samples);
    T=length(Candidate_TFs);
    F=length(Candidate_Peaks);
    G=length(Candidate_Gene_Symbol);

    fprintf(2, 'Initial integration associated %d TFs, %d DA peaks, and %d DE genes together\n\n', T, F, G)

    %********** Model variable prior calculation *****************************
    %Note: in MCMC, initial values may not be that important, but using more
    %importantive prior could speed up the sampling convergence

    fprintf('MAGICAL model initialization\n\n')

    [P_prior, P_mean, P_var, B_prior, B_mean, B_var, B_prob, L_prior, L_mean, L_var, L_prob]=...
        MAGICAL_prior(Candidate_TF_log2Count,Candidate_Peak_log2Count,Candidate_Gene_log2Count,Candidate_TF_Peak_Binding, Candidate_Peak_Gene_looping, S, T, F, G);

    MAGICAL_post(p).TFs=Candidate_TFs;
    MAGICAL_post(p).Peak_ID=Candidate_Peak_ID;
    MAGICAL_post(p).Peaks=Candidate_Peaks;
    MAGICAL_post(p).Genes=Candidate_Gene_Symbol;
    MAGICAL_post(p).Gene_TSS=Candidate_Gene_TSS;
    MAGICAL_post(p).TF_Peak_Binding_weight=B_prior;
    MAGICAL_post(p).TF_Peak_Binding_prob=B_prob;
    MAGICAL_post(p).Peak_Gene_Looping_weight=L_prior;
    MAGICAL_post(p).Peak_Gene_Looping_prob=L_prob;

    %********** Initial round of sampling, using their prior values **********

    A=Candidate_Peak_log2Count;
    R=Candidate_Gene_log2Count;

    P=P_prior;
    B=B_prior;
    L=L_prior;

    B_state=full(Candidate_TF_Peak_Binding);
    L_state=full(Candidate_Peak_Gene_looping);

    B_state_frq=full(Candidate_TF_Peak_Binding);
    L_state_frq=full(Candidate_Peak_Gene_looping);

    alpha_A=1;
    beta_A=1;
    ATAC_fitting_residue=A-B*P;
    sigma_A_noise=var(ATAC_fitting_residue(:));

    % RNA_fitting_residue=R-L'*(B*P);
    % sigma_R_noise=var(RNA_fitting_residue(:));
    alpha_R=1;
    beta_R=1;
    RNA_fitting_residue=R-L'*A;
    sigma_R_noise=var(RNA_fitting_residue(:));

    %************************* MAGICAL sampling ******************************
    iteration_num=10000;
    iteration_seg=iteration_num/10;
    for i=1:iteration_num

        %********** scATAC-seq fitting and variable sampling  ****************

        %Step 1: TF activity sampling
        P = TF_activity_P_sampling(A, B, P, P_mean, P_var, sigma_A_noise);

        %Step 2: TF-peak binding weight sampling
        B = TF_peak_binding_B_sampling(A, B, P, B_state, B_mean, B_var, sigma_A_noise);

        %Step 3: TF-peak binding state update
        [B_state, B] = TF_peak_binary_binding_B_state_sampling(A, B, P, B_state, B_mean, B_var, B_prob, sigma_A_noise);

        %Step 4: ATAC fitting residue variance control
        ATAC_fitting_residue=A-B*P;
        sigma_A_noise = 1/gamrnd(alpha_A+1/2,1/(beta_A+sum(sum(ATAC_fitting_residue.^2))/(2*F*S)));
        %         sigma_A_noise = (beta_A+sum(sum(ATAC_fitting_residue.^2))/(2*F*S))/chi2rnd(2*alpha_A+1);
        %         aa(i)=sigma_A_noise;

        %********** scRNA-seq fitting and variable sampling  *****************
        A_estimate=B*P;% or true A
        %Step 5: Peak-Gene looping weight sampling
        L = Peak_gene_looping_L_samping(R, L, A_estimate, L_state, L_mean, L_var, sigma_R_noise);

        %Step 6: Peak-Gene looping state update
        [L_state, L]=Peak_gene_binary_looping_L_state_samping(R, L, A_estimate, L_state, L_mean, L_var, L_prob, sigma_R_noise);

        %Step 7: RNA fitting residue variance control
        RNA_fitting_residue=R-L'*A_estimate;
        sigma_R_noise = 1/gamrnd(alpha_R+1/2,1/(beta_R+sum(sum(RNA_fitting_residue.^2))/(2*G*S)));
        %         sigma_R_noise = (beta_R+sum(sum(RNA_fitting_residue.^2))/(2*G*S))/chi2rnd(2*alpha_R+1);
        %         bb(i) = sigma_R_noise;

        %Step 8: Sample Summary
        B_state_frq=B_state_frq+B_state;
        L_state_frq=L_state_frq+L_state;

        if mod(i,iteration_seg)==0
            fprintf(2, 'MAGICAL finished %d percent\n\n', 10*i/iteration_seg)
        end
    end


    % output posterior TFs, Peaks, Genes, TF-Peak binding and Peak-Gene looping

    MAGICAL_post.TF_Peak_Binding_prob=B_state_frq/iteration_num;
    MAGICAL_post.Peak_Gene_Looping_prob=L_state_frq/iteration_num;
    fprintf('MAGICAL integrated %d peaks and %d genes together.\n\n', length(unique(xx)), length(unique(yy)))

end

save('S16_S19_trajectory_MAGICAL.mat')

clear

load S16_S19_trajectory_MAGICAL.mat

CollecTRI_regulons=readtable('CollecTRI_regulons.txt', 'delimiter', ',');

fid=fopen('V2_Kidney_S16_S19_all_samples_trajectory_MAGICAL_circuits.txt', 'w');
fprintf(fid, 'Cluster_ID\tGene_symbol\tGene_chr\tGene_TSS\tGene_pct_1\tGene_pct_2\tGene_sc_Log2FC\tGene_sc_p_val_adj\tPeak_ID\tPeak_pct_1\tPeak_pct_2\tPeak_sc_log2FC\tPeak_sc_p_val_adj\tATAC_RNA_looping_weight\tATAC_RNA_Bayesian_prob\tTFs\tCollecTRI_regulons\n');
for p=1:length(Candidate_clusters)
    p
    [xx,yy]=find(MAGICAL_post(p).Peak_Gene_Looping_prob>0.95);
    for i=1:length(xx)
        Gene_index_1=find(strcmp(RNA_sc_diff{p}.gene, MAGICAL_post(p).Genes(yy(i)))>0);
        Peak_index_1=find(strcmp(Peak_sc_diff{p}.gene, MAGICAL_post(p).Peak_ID(xx(i)))>0);

        fprintf(fid, '%s\t%s\tchr%d\t%d\t%f\t%f\t%f\t%G\t%s\t%f\t%f\t%f\t%G\t%f\t%f\t', ...
            Candidate_clusters{p},...
            MAGICAL_post(p).Genes{yy(i)}, MAGICAL_post(p).Gene_TSS(yy(i), :),...
            RNA_sc_diff{p}.pct_1(Gene_index_1(1)), RNA_sc_diff{p}.pct_2(Gene_index_1(1)),...
            RNA_sc_diff{p}.avg_log2FC(Gene_index_1(1)), RNA_sc_diff{p}.p_val_adj(Gene_index_1(1)),...
            MAGICAL_post(p).Peak_ID{xx(i)},...
            Peak_sc_diff{p}.pct_1(Peak_index_1(1)), Peak_sc_diff{p}.pct_2(Peak_index_1(1)),...
            Peak_sc_diff{p}.avg_log2FC(Peak_index_1(1)), Peak_sc_diff{p}.p_val(Peak_index_1(1)),...
            full(MAGICAL_post(p).Peak_Gene_Looping_weight(xx(i),yy(i))),...
            full(MAGICAL_post(p).Peak_Gene_Looping_prob(xx(i),yy(i))));

        [TF_prob, TF_index]=sort(full(MAGICAL_post(p).TF_Peak_Binding_prob(xx(i),:)), 'descend');
        for t=1:length(TF_index)
            if TF_prob(t)>0.8
                fprintf(fid, '%s (%f), ', MAGICAL_post(p).TFs{TF_index(t)}, MAGICAL_post(p).TF_Peak_Binding_weight(xx(i),TF_index(t)));
            end
        end
        fprintf(fid, '\t');


        target_index=find(strcmp(CollecTRI_regulons.target, MAGICAL_post(p).Genes{yy(i)})>0);
        if ~isempty(target_index)
            cn=0;
            for t=1:length(TF_index)
                if TF_prob(t)>0.8
                    source_index=find(strcmp(CollecTRI_regulons.source(target_index), MAGICAL_post(p).TFs{TF_index(t)})>0);
                    if ~isempty(source_index)
                        cn=cn+1;
                        fprintf(fid, '%s (%f), ', MAGICAL_post(p).TFs{TF_index(t)}, MAGICAL_post(p).TF_Peak_Binding_weight(xx(i),TF_index(t)));
                    end
                end
            end

        else
            fprintf(fid, 'Gene not included');
        end
        fprintf(fid, '\n');
    end
end


MAGICAL_circuits=readtable('V2_Kidney_S16_S19_all_samples_trajectory_MAGICAL_circuits.txt');

aa=~cellfun(@isempty,string(MAGICAL_circuits.TFs));
length(unique(MAGICAL_circuits.Gene_symbol(aa>0)))
length(unique(MAGICAL_circuits.Peak_ID(aa>0)))


for p=1:length(Candidate_clusters)
    [xx,yy]=find(MAGICAL_post(p).Peak_Gene_Looping_prob>0.95);
    [xx,yy]=find(MAGICAL_post(p).TF_Peak_Binding_weight(unique(xx),:)>0.8);
    if p==1
        selected_TFs=MAGICAL_post(p).TFs(unique(yy));
    else
        selected_TFs=union(selected_TFs,MAGICAL_post(p).TFs(unique(yy)));
    end
end
length(selected_TFs)

%validate TF-gene regulons
gene_list=unique(MAGICAL_circuits.Gene_symbol);
aa=~cellfun(@isempty,string(MAGICAL_circuits.CollecTRI_regulons));
mm=~strcmp(string(MAGICAL_circuits.CollecTRI_regulons), "Gene not included");
validated_genes=unique(MAGICAL_circuits.Gene_symbol(aa>0 & mm>0));
bb(1)=length(gene_list);
bb(2)=length(validated_genes);
[h,p_value]=fishertest([length(validated_genes),length(gene_list);6564,36588], 'tail', 'right')

length(validated_genes)/length(unique(MAGICAL_circuits.Gene_symbol(mm>0)))
[h,p_value]=fishertest([length(validated_genes),length(unique(MAGICAL_circuits.Gene_symbol(mm>0)));6564,36588], 'tail', 'right')


