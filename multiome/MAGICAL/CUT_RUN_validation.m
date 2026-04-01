close
clc
clear

MAGICAL_circuits=readtable('PT _TAL_Firbo_MAGICAL_circuits.xlsx', 'Sheet','V2_Kidney_P15_P4_all_samples_tr');

%CUT&RUN

[AKI_H3K27ac_peaks.chr, AKI_H3K27ac_peaks.point1, AKI_H3K27ac_peaks.point2]=textread('CUT&RUN_Nov_2024/AKI/AKI_biopsy_H3K27ac_S.2401.009918_S.2307.003357_110824_broad_RPKMnorm.bigwig_peaks.txt', '%s %d %d');
[AKI_H3K27me3_peaks.chr, AKI_H3K27me3_peaks.point1, AKI_H3K27me3_peaks.point2]=textread('CUT&RUN_Nov_2024/AKI/AKI_biopsy_H3K27me3_S.2401.009918_S.2307.003357_110824_broad_RPKMnorm.bigwig_peaks.txt', '%s %d %d');
[AKI_H3K4me1_peaks.chr, AKI_H3K4me1_peaks.point1, AKI_H3K4me1_peaks.point2]=textread('CUT&RUN_Nov_2024/AKI/AKI_biopsy_H3K4me1_S.2401.009918_S.2307.003357_110824_broad_RPKMnorm.bigwig_peaks.txt', '%s %d %d');
[AKI_H3K4me3_peaks.chr, AKI_H3K4me3_peaks.point1, AKI_H3K4me3_peaks.point2]=textread('CUT&RUN_Nov_2024/AKI/AKI_biopsy_H3K4me3_S.2401.009918_S.2307.003357_110824_narrow_RPKMnorm.bigwig_peaks.txt', '%s %d %d');

AKI_H3K27ac_peaks_flag=zeros(height(MAGICAL_circuits),1);
AKI_H3K27me3_peaks_flag=zeros(height(MAGICAL_circuits),1);
AKI_H3K4me1_peaks_flag=zeros(height(MAGICAL_circuits),1);
AKI_H3K4me3_peaks_flag=zeros(height(MAGICAL_circuits),1);

for p=1:height(MAGICAL_circuits)
    index1=find(strcmp(AKI_H3K27ac_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(AKI_H3K27ac_peaks.point1(index1)+AKI_H3K27ac_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        AKI_H3K27ac_peaks_flag(p)=1;
    end

    index1=find(strcmp(AKI_H3K27me3_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(AKI_H3K27me3_peaks.point1(index1)+AKI_H3K27me3_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        AKI_H3K27me3_peaks_flag(p)=1;
    end

    index1=find(strcmp(AKI_H3K4me1_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(AKI_H3K4me1_peaks.point1(index1)+AKI_H3K4me1_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        AKI_H3K4me1_peaks_flag(p)=1;
    end

    index1=find(strcmp(AKI_H3K4me3_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(AKI_H3K4me3_peaks.point1(index1)+AKI_H3K4me3_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        AKI_H3K4me3_peaks_flag(p)=1;
    end
end

AKI_CUT_Run_flag=[AKI_H3K27ac_peaks_flag, AKI_H3K4me1_peaks_flag, AKI_H3K4me3_peaks_flag, AKI_H3K27me3_peaks_flag];



[DKD_H3K27ac_peaks.chr, DKD_H3K27ac_peaks.point1, DKD_H3K27ac_peaks.point2]=textread('CUT&RUN_Nov_2024/DKD/DKD_biopsy_H3K27ac_S.2108.021702_S.2102.006723_S.2311.011848_S.2305.012652_110824_broad_peaks.txt', '%s %d %d');
[DKD_H3K27me3_peaks.chr, DKD_H3K27me3_peaks.point1, DKD_H3K27me3_peaks.point2]=textread('CUT&RUN_Nov_2024/DKD/DKD_biopsy_H3K27me3_S.2108.021702_S.2102.006723_S.2311.011848_S.2305.012652_110824_broad_peaks.txt', '%s %d %d');
[DKD_H3K4me1_peaks.chr, DKD_H3K4me1_peaks.point1, DKD_H3K4me1_peaks.point2]=textread('CUT&RUN_Nov_2024/DKD/DKD_biopsy_H3K4me1_S.2108.021702_S.2102.006723_S.2311.011848_S.2305.012652_110824_broad_peaks.txt', '%s %d %d');
[DKD_H3K4me3_peaks.chr, DKD_H3K4me3_peaks.point1, DKD_H3K4me3_peaks.point2]=textread('CUT&RUN_Nov_2024/DKD/DKD_biopsy_H3K4me3_S.2108.021702_S.2102.006723_S.2311.011848_S.2305.012652_110824_narrow_peaks.txt', '%s %d %d');
 
DKD_H3K27ac_peaks_flag=zeros(height(MAGICAL_circuits),1);
DKD_H3K27me3_peaks_flag=zeros(height(MAGICAL_circuits),1);
DKD_H3K4me1_peaks_flag=zeros(height(MAGICAL_circuits),1);
DKD_H3K4me3_peaks_flag=zeros(height(MAGICAL_circuits),1);

for p=1:height(MAGICAL_circuits)
    index1=find(strcmp(DKD_H3K27ac_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(DKD_H3K27ac_peaks.point1(index1)+DKD_H3K27ac_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        DKD_H3K27ac_peaks_flag(p)=1;
    end

    index1=find(strcmp(DKD_H3K27me3_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(DKD_H3K27me3_peaks.point1(index1)+DKD_H3K27me3_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        DKD_H3K27me3_peaks_flag(p)=1;
    end

    index1=find(strcmp(DKD_H3K4me1_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(DKD_H3K4me1_peaks.point1(index1)+DKD_H3K4me1_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        DKD_H3K4me1_peaks_flag(p)=1;
    end

    index1=find(strcmp(DKD_H3K4me3_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(DKD_H3K4me3_peaks.point1(index1)+DKD_H3K4me3_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        DKD_H3K4me3_peaks_flag(p)=1;
    end
end

DKD_CUT_Run_flag=[DKD_H3K27ac_peaks_flag, DKD_H3K4me1_peaks_flag, DKD_H3K4me3_peaks_flag, DKD_H3K27me3_peaks_flag];


[DMR_H3K27ac_peaks.chr, DMR_H3K27ac_peaks.point1, DMR_H3K27ac_peaks.point2]=textread('CUT&RUN_Nov_2024/DM-R/H3K27ac_S2203.008743_20_biopsy_peaks.txt', '%s %d %d');
[DMR_H3K27me3_peaks.chr, DMR_H3K27me3_peaks.point1, DMR_H3K27me3_peaks.point2]=textread('CUT&RUN_Nov_2024/DM-R/H3K27me3_biopsy_21_S2203008743_peaks.txt', '%s %d %d');
[DMR_H3K4me1_peaks.chr, DMR_H3K4me1_peaks.point1, DMR_H3K4me1_peaks.point2]=textread('CUT&RUN_Nov_2024/DM-R/H3K4me1_biopsy_23_S2203008743_peaks.txt', '%s %d %d');
[DMR_H3K4me3_peaks.chr, DMR_H3K4me3_peaks.point1, DMR_H3K4me3_peaks.point2]=textread('CUT&RUN_Nov_2024/DM-R/H3K4me3.S2203.008743_22_biopsy_peaks.txt', '%s %d %d');

DMR_H3K27ac_peaks_flag=zeros(height(MAGICAL_circuits),1);
DMR_H3K27me3_peaks_flag=zeros(height(MAGICAL_circuits),1);
DMR_H3K4me1_peaks_flag=zeros(height(MAGICAL_circuits),1);
DMR_H3K4me3_peaks_flag=zeros(height(MAGICAL_circuits),1);

for p=1:height(MAGICAL_circuits)
    index1=find(strcmp(DMR_H3K27ac_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(DMR_H3K27ac_peaks.point1(index1)+DMR_H3K27ac_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        DMR_H3K27ac_peaks_flag(p)=1;
    end

    index1=find(strcmp(DMR_H3K27me3_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(DMR_H3K27me3_peaks.point1(index1)+DMR_H3K27me3_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        DMR_H3K27me3_peaks_flag(p)=1;
    end

    index1=find(strcmp(DMR_H3K4me1_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(DMR_H3K4me1_peaks.point1(index1)+DMR_H3K4me1_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        DMR_H3K4me1_peaks_flag(p)=1;
    end

    index1=find(strcmp(DMR_H3K4me3_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(DMR_H3K4me3_peaks.point1(index1)+DMR_H3K4me3_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        DMR_H3K4me3_peaks_flag(p)=1;
    end
end

DMR_CUT_Run_flag=[DMR_H3K27ac_peaks_flag, DMR_H3K4me1_peaks_flag, DMR_H3K4me3_peaks_flag, DMR_H3K27me3_peaks_flag];


[HCKD_H3K27ac_peaks.chr, HCKD_H3K27ac_peaks.point1, HCKD_H3K27ac_peaks.point2]=textread('CUT&RUN_Nov_2024/H-CKD/HCKD_biopsy_H3K27ac_S.2108.021749_S.2305.012746_S.2310.018319_110824_broad_RPKMnorm.bigwig_peaks.txt', '%s %d %d');
[HCKD_H3K27me3_peaks.chr, HCKD_H3K27me3_peaks.point1, HCKD_H3K27me3_peaks.point2]=textread('CUT&RUN_Nov_2024/H-CKD/HCKD_biopsy_H3K27me3_S.2108.021749_S.2305.012746_S.2310.018319_110824_broad_RPKMnorm.bigwig_peaks.txt', '%s %d %d');
[HCKD_H3K4me1_peaks.chr, HCKD_H3K4me1_peaks.point1, HCKD_H3K4me1_peaks.point2]=textread('CUT&RUN_Nov_2024/H-CKD/HCKD_biopsy_H3K4me1_S.2108.021749_S.2305.012746_S.2310.018319_110824_broad_RPKMnorm.bigwig_peaks.txt', '%s %d %d');
[HCKD_H3K4me3_peaks.chr, HCKD_H3K4me3_peaks.point1, HCKD_H3K4me3_peaks.point2]=textread('CUT&RUN_Nov_2024/H-CKD/HCKD_biopsy_H3K4me3_S.2108.021749_S.2305.012746_S.2310.018319_110824_narrow_RPKMnorm.bigwig_peaks.txt', '%s %d %d');

HCKD_H3K27ac_peaks_flag=zeros(height(MAGICAL_circuits),1);
HCKD_H3K27me3_peaks_flag=zeros(height(MAGICAL_circuits),1);
HCKD_H3K4me1_peaks_flag=zeros(height(MAGICAL_circuits),1);
HCKD_H3K4me3_peaks_flag=zeros(height(MAGICAL_circuits),1);

for p=1:height(MAGICAL_circuits)
    index1=find(strcmp(HCKD_H3K27ac_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(HCKD_H3K27ac_peaks.point1(index1)+HCKD_H3K27ac_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        HCKD_H3K27ac_peaks_flag(p)=1;
    end

    index1=find(strcmp(HCKD_H3K27me3_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(HCKD_H3K27me3_peaks.point1(index1)+HCKD_H3K27me3_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        HCKD_H3K27me3_peaks_flag(p)=1;
    end

    index1=find(strcmp(HCKD_H3K4me1_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(HCKD_H3K4me1_peaks.point1(index1)+HCKD_H3K4me1_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        HCKD_H3K4me1_peaks_flag(p)=1;
    end

    index1=find(strcmp(HCKD_H3K4me3_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(HCKD_H3K4me3_peaks.point1(index1)+HCKD_H3K4me3_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        HCKD_H3K4me3_peaks_flag(p)=1;
    end
end

HCKD_CUT_Run_flag=[HCKD_H3K27ac_peaks_flag, HCKD_H3K4me1_peaks_flag, HCKD_H3K4me3_peaks_flag, HCKD_H3K27me3_peaks_flag];
 
[REF_H3K27ac_peaks.chr, REF_H3K27ac_peaks.point1, REF_H3K27ac_peaks.point2]=textread('CUT&RUN_Nov_2024/REF/H3K27ac_reference_21-020.Rep1_2_3.687_688_LN2_OCT.KTRC.A_B_D.22-078.Rep1_2_3_4_103024_broad_peaks.txt', '%s %d %d');
[REF_H3K27me3_peaks.chr, REF_H3K27me3_peaks.point1, REF_H3K27me3_peaks.point2]=textread('CUT&RUN_Nov_2024/REF/H3K27me3_reference_21-020.Rep1_2_3.687_688_LN2_OCT.KTRC.A_B_D.22-078.Rep1_2_3_4_103024_broad_peaks.txt', '%s %d %d');
[REF_H3K4me1_peaks.chr, REF_H3K4me1_peaks.point1, REF_H3K4me1_peaks.point2]=textread('CUT&RUN_Nov_2024/REF/H3K4me1_reference_21-020.Rep1_2_3.22-078.Rep1_2_3_4_103024_broad_peaks.txt', '%s %d %d');
[REF_H3K4me3_peaks.chr, REF_H3K4me3_peaks.point1, REF_H3K4me3_peaks.point2]=textread('CUT&RUN_Nov_2024/REF/H3K4me3_reference_21-020.Rep1_2_3.22-078.Rep1_2_3_4_103024_narrow_peaks.txt', '%s %d %d');

REF_H3K27ac_peaks_flag=zeros(height(MAGICAL_circuits),1);
REF_H3K27me3_peaks_flag=zeros(height(MAGICAL_circuits),1);
REF_H3K4me1_peaks_flag=zeros(height(MAGICAL_circuits),1);
REF_H3K4me3_peaks_flag=zeros(height(MAGICAL_circuits),1);

for p=1:height(MAGICAL_circuits)
    index1=find(strcmp(REF_H3K27ac_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(REF_H3K27ac_peaks.point1(index1)+REF_H3K27ac_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        REF_H3K27ac_peaks_flag(p)=1;
    end

    index1=find(strcmp(REF_H3K27me3_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(REF_H3K27me3_peaks.point1(index1)+REF_H3K27me3_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        REF_H3K27me3_peaks_flag(p)=1;
    end

    index1=find(strcmp(REF_H3K4me1_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(REF_H3K4me1_peaks.point1(index1)+REF_H3K4me1_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        REF_H3K4me1_peaks_flag(p)=1;
    end

    index1=find(strcmp(REF_H3K4me3_peaks.chr, MAGICAL_circuits.Peak_chr(p))>0);
    index2=find(abs(REF_H3K4me3_peaks.point1(index1)+REF_H3K4me3_peaks.point2(index1)-MAGICAL_circuits.Peak_point1(p)-MAGICAL_circuits.Peak_point2(p))/2<1000);
    if ~isempty(index2)
        REF_H3K4me3_peaks_flag(p)=1;
    end
end

REF_CUT_Run_flag=[REF_H3K27ac_peaks_flag, REF_H3K4me1_peaks_flag, REF_H3K4me3_peaks_flag, REF_H3K27me3_peaks_flag];

aa=[AKI_CUT_Run_flag, DKD_CUT_Run_flag, DMR_CUT_Run_flag, HCKD_CUT_Run_flag, REF_CUT_Run_flag];

