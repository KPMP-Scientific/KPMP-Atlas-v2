import { tsvParse } from 'd3-dsv';

const abbrevTsv = /*`abbreviation (see cluster_dm table)	cluster_name (CZI)*/ `clusterAbbrev	clusterName
PEC	Parietal Epithelial Cell
POD	Podocyte
dPOD	Podocyte (degenerative)
EC-GC	Glomerular Capillary Endothelial Cell
EC-GC-FILIP1+	"Glomerular Capillary Endothelial Cell, FILIP1+"
aEC-GC	Glomerular Capillary Endothelial Cell (adaptive/maladaptive)
dEC-GC	Glomerular Capillary Endothelial Cell (degenerative)
MC	Mesangial Cell
PT	Proximal Tubule Epithelial Cell
aPT	Proximal Tubule Epithelial Cell (adaptive/maladaptive)
aPT-KIM1+	Proximal Tubule Epithelial Cell, KIM1+ (adaptive/maladaptive)
aPT-KIM1+PROM1+	Proximal Tubule Epithelial Cell, KIM1+, PROM1+ (adaptive/maladaptive)
aPT-VCAM1+	Proximal Tubule Epithelial Cell, VCAM1+ (adaptive/maladaptive)
cycPT	Proximal Tubule Epithelial Cell (cycling)
dPT	Proximal Tubule Epithelial Cell (degenerative)
frPT	Proximal Tubule Epithelial Cell (failed repair)
frPT_(maPT)	Proximal Tubule Epithelial Cell (failed repair) (?)
PT-S1	Proximal Tubule Epithelial Cell Segment 1
aPT1	Proximal Tubule Epithelial Cell Type 1 (adaptive/maladaptive)
dPT-S1	Proximal Tubule Epithelial Cell Segment 1 (degenerative)
PT-S1/S2	Proximal Tubule Epithelial Cell Segment 1/Segment 2
aPT-S1/S2	Proximal Tubule Epithelial Cell Segment 1/Segment 2 (adaptive/maladaptive)
aPT-S1_S2	Proximal Tubule Epithelial Cell Segment 1/Segment 2 (adaptive/maladaptive)
frPT-S1/S2	Proximal Tubule Epithelial Cell Segment 1/Segment 2 (failed repair)
frPT-S1_S2	Proximal Tubule Epithelial Cell Segment 1/Segment 2 (failed repair)
PT-S2	Proximal Tubule Epithelial Cell Segment 2
aPT2	Proximal Tubule Epithelial Cell Type 2 (adaptive/maladaptive)
dPT-S2	Proximal Tubule Epithelial Cell Segment 2 (degenerative)
PT-S2/S3	Proximal Tubule Epithelial Cell Segment 2/Segment 3
PT-S3	Proximal Tubule Epithelial Cell Segment 3
aPT-S3	Proximal Tubule Epithelial Cell Segment 3 (adaptive/maladaptive)
dPT-S3	Proximal Tubule Epithelial Cell Segment 3 (degenerative)
frPT-S3	Proximal Tubule Epithelial Cell Segment 3 (failed repair)
DTL	Descending Thin Limb of Loop of Henle
aDTL	Descending Thin Limb of Loop of Henle (adaptive/maladaptive)
dDTL	Descending Thin Limb of Loop of Henle (degenerative)
DTL1	Descending Thin Limb of Loop of Henle Cell Type 1
DTL2	Descending Thin Limb of Loop of Henle Cell Type 2
aDTL2	Descending Thin Limb of Loop of Henle Cell Type 2 (adaptive/maladaptive)
DTL3	Descending Thin Limb of Loop of Henle Cell Type 3
ATL	Ascending Thin Limb of Loop of Henle
aATL	Ascending Thin Limb of Loop of Henle (adaptive/maladaptive)
dATL	Ascending Thin Limb of Loop of Henle (degenerative)
TAL	Thick Ascending Limb Cell of Loop of Henle
aTAL	Thick Ascending Limb of Loop of Henle Cell (adaptive/maladaptive)
cycTAL	Thick Ascending Limb of Loop of Henle Cell (cycling)
dTAL	Thick Ascending Limb of Loop of Henle Cell (degenerative)
frTAL	Thick Ascending Limb of Loop of Henle Cell (failed repair)
aTAL1	Thick Ascending Limb of Loop of Henle Cell Type 1 (adaptive/maladaptive)
aTAL2	Thick Ascending Limb of Loop of Henle Cell Type 2 (adaptive/maladaptive)
C-TAL	Cortical Thick Ascending Limb of Loop of Henle Cell
dC-TAL	Cortical Thick Ascending Limb of Loop of Henle Cell (degenerative)
C-TAL-A	Cortical Thick Ascending Limb of Loop of Henle Cell Type A
C-TAL-B	Cortical Thick Ascending Limb of Loop of Henle Cell Type B
C/M-TAL	Cortico-Medullary Thick Ascending Limb of Loop of Henle Cell
C/M-TAL-A	Cortico-Medullary Thick Ascending Limb of Loop of Henle Cell Type A
C_M-TAL-A	Cortico-Medullary Thick Ascending Limb of Loop of Henle Cell Type A
C/M-TAL-B	Cortico-Medullary Thick Ascending Limb of Loop of Henle Cell Type B
C_M-TAL-B	Cortico-Medullary Thick Ascending Limb of Loop of Henle Cell Type B
M-TAL	Medullary Thick Ascending Limb of Loop of Henle Cell
dM-TAL	Medullary Thick Ascending Limb of Loop of Henle Cell (degenerative)
MD	Macula Densa Cell
DCT	Distal Convoluted Tubule Cell
aDCT	Distal Convoluted Tubule Cell (adaptive/maladaptive)
dDCT	Distal Convoluted Tubule Cell (degenerative)
frDCT	Distal Convoluted Tubule Cell (failed repair)
DCT1	Distal Convoluted Tubule Cell Type 1
DCT2	Distal Convoluted Tubule Cell Type 2
DN	Distal Nephron
dDN	Distal Nephron (degenerative)
CNT	Connecting Tubule Cell
aCNT	Connecting Tubule Cell (adaptive/maladaptive)
dCNT	Connecting Tubule Cell (degenerative)
CNT-PC	Connecting Tubule Principal Cell
aCNT-PC	Connecting Tubule Principal Cell (adaptive/maladaptive)
dCNT-PC	Connecting Tubule Principal Cell (degenerative)
C-PC	Cortical Principal Cell
CCD-PC	Cortical Collecting Duct Principal Cell
aCCD-PC	Cortical Collecting Duct Principal Cell (adaptive/maladaptive)
dCCD-PC	Cortical Collecting Duct Principal Cell (degenerative)
IMCD	Inner Medullary Collecting Duct Cell
dIMCD	Inner Medullary Collecting Duct Cell (degenerative)
M-PC	Medullary Principal Cell
OMCD-PC	Outer Medullary Collecting Duct Principal Cell
aOMCD-PC	Outer Medullary Collecting Duct Principal Cell (adaptive/maladaptive)
dOMCD-PC	Outer Medullary Collecting Duct Principal Cell (degenerative)
tPC-IC	Principal-Intercalated Cell (transitional)
OMCD-IC-A	Outer Medullary Collecting Duct Intercalated Cell Type A
dOMCD-IC-A	Outer Medullary Collecting Duct Intercalated Cell Type A (degenerative)
tOMCD-PC-IC	Outer Medullary Collecting Duct Principal-Intercalated Cell (transitional)
IC	Intercalated Cell
IC-A	Intercalated Cell Type A
aIC-A	Collecting Duct Intercalated Cell Type A (adaptive/maladaptive)
CCD-IC-A	Cortical Collecting Duct Intercalated Cell Type A
dCCD-IC-A	Cortical Collecting Duct Intercalated Cell Type A (degenerative)
IC-B	Intercalated Cell Type B
PapE	Papillary Tip Epithelial Cell
FIB	Fibroblast
IM-FIB	Inner Medullary Fibroblast
dIM-FIB	Inner Medullary Fibroblast (degenerative)
OM-FIB	Outer Medullary Fibroblast
dOM-FIB	Outer Medullary Fibroblast (degenerative)
C-FIB	Cortical Interstitial Fibroblast
C/M-FIB	Corticomedullary Fibroblast
C_M-FIB	Corticomedullary Fibroblast
pvFIB	Perivascular Fibroblast
pvFIB-PI16+	"Perivascular Fibroblast, PI16+CD34+"
pvFIB-CD34+	"Perivascular Fibroblast, CD34+"
pvFIB-RSPO3+	"Perivascular Fibroblast, RSPO3+"
dFIB	Interstitial Fibroblast (degenerative)
C-FIB	Inflammatory Cortical Interstitial Fibroblast
C-intFIB	Inflammatory Cortical Interstitial Fibroblast
C-FIB-OSMRhi	OSMRhi Inflammatory Cortical Interstitial Fibroblast
C-intFIB-OSMRhi	OSMRhi Inflammatory Cortical Interstitial Fibroblast
C-FIB-OSMRlo	OSMRlo Inflammatory Cortical Interstitial Fibroblast
C-intFIB-OSMRlo	OSMRlo Inflammatory Cortical Interstitial Fibroblast
C-FIB-PATH	Pathology (Injury) Expanded Cortical Interstitial Fibroblast
C-intFIB-PATH	Pathology (Injury) Expanded Cortical Interstitial Fibroblast
MYOF	Myofibroblast
C-MYOF	Cortical Interstitial Myofibroblast
C-intMYOF	Cortical Interstitial Myofibroblast
IM-pvMYOF	Inner Medullary Perivascular Myofibroblast
pvMYOF	Perivascular Myofibroblast
MON	Classical Monocyte
cMON	Classical Monocyte
ncMON	Non-Classical Monocyte
MAC	Macrophage
cycMAC	Macrophage (cycling)
moMAC	Monocyte-Derived Inflammatory Macrophage
moMAC-C3+	"Monocyte-Derived Macrophage, C3+"
moMAC-CXCL10+	"Monocyte-Derived Inflammatory Macrophage, CXCL10+"
moMAC-LUCAT1+	"Monocyte-Derived Macrophage, LUCAT1+"
moMAC-HBEGF+	"Monocyte-Derived Inflammatory Macrophage, HBEGF+"
moMAC-FAM	"Monocyte-Derived Inflammatory Macrophage, FAM"
resMAC	Resident Macrophage
resMAC-HLAIIhi	"Resident Macrophage, LYVE1+ HLAIIhi"
resMAC-HLAIIhi_	"Resident Macrophage, LYVE1+ HLAIIhi_"
resMAC-LYVE1+	"Resident Macrophage, LYVE1+"
moFAM	Monocyte-Derived Fibrosis-Associated Macrophage
DC	
cDC1	Classical Type 1 Dendritic Cell
cDC2	Classical Type 2 Dendritic Cell
mDC	Migratory Dendritic Cell
pDC	Plasmacytoid Dendritic Cell
ERY	Erythrocyte
T	
cycT	T Cell (cycling)
T-CXCR4+	"T Cell, CXCR4+"
Naive_Th	"T Helper Cell, Naive"
T-REG	Regulatory T Cell
Th-memory	T Helper memory Cell
CD8+ TEM/TEMRA	"Cytotoxic T Cell, Tem/Temra (Terminally Differentiated) CD8+"
CD8+_TEM_TEMRA	"Cytotoxic T Cell, Tem/Temra (Terminally Differentiated) CD8+"
CD8+ TEM/TRM	"Cytotoxic T Cell, Tem/Trm CD8+"
CD8+_TEM_TRM	"Cytotoxic T Cell, Tem/Trm CD8+"
CD8+ TEMRA	"Cytotoxic T Cell, Temra (Terminally Differentiated) CD8+"
Th17	"T Helper Cell, Type 17"
CD8+ T-CYT	"Cytotoxic T Cell, CD8+"
CD8+ TEM	"Cytotoxic T Effector Memory (TEM) Cell, CD8+"
CD8+ TRM	"Cytotoxic T Resident Memory (TRM) Cell, CD8+"
CD8+ T-STR	"Cytotoxic Stressed T Cell, CD8+"
ILC3	Innate Lymphoid Cell Subpopulation III
MAIT	Mucosal-Associated Invariant T Cell
NK	Natural Killer Cell
CD16+ NK	"Natural Killer Cell, CD16+"
CD16- NK	"Natural Killer Cell, CD16-"
B	B Cell
B memory	B Memory Cell
B activated	B Intermediate Cell
B naive	"B Cell, Naive"
PL	Plasma Cell
MAST	Mast Cell
N	Neutrophil
SC/NEU	Schwann Cell / Neural
SC_NEU	Schwann Cell / Neural
Ad	Adipocyte
NA	ambiguous low quality
VSMC	Vascular Smooth Muscle Cell
dVSMC	Vascular Smooth Muscle Cell (degenerative)
VSMC/P	Vascular Smooth Muscle Cell / Pericyte
VSMC_P	Vascular Smooth Muscle Cell / Pericyte
M-VSMC/P	Medullary Vascular Smooth Muscle Cell / Pericyte
M-VSMC_P	Medullary Vascular Smooth Muscle Cell / Pericyte
REN	Renin+ Juxtaglomerular Granular Cell
EC	Endothelial Cell
cycEC	Endothelial Cell (cycling)
EC-AA	Afferent Artery / Arteriole Endothelial Cell
EC-EA	Efferent Arteriole Endothelial Cell
dEC-EA	Efferent Arteriole (degenerative)
iaEC-AVR	Interferon Activated Vasa Recta Endothelial Cell
EC-AVR	Ascending Vasa Recta Endothelial Cell
dEC-AVR	Ascending Vasa Recta Endothelial Cell (degenerative)
EC-DVR	Descending Vasa Recta Endothelial Cell
aEC-DVR	Descending Vasa Recta Endothelial Cell (adaptive/maladaptive)
EC-PTC	Peritubular Capillary Endothelial Cell
aEC-PTC	Peritubular Capillary Endothelial Cell (adaptive/maladaptive)
infEC-PTC	Inflammatory Peritubular Capillary Endothelial Cell
iaEC-PTC	Interferon Activated Peritubular Capillary Endothelial Cell
C-EC-PTC	Cortical Peritubular Capillary Endothelial Cell
dEC-PTC	Cortical Peritubular Capillary Endothelial Cell (degenerative)
angEC-PTC	Angiogenic Cortical Peritubular Capillary Endothelial Cell
M-EC-PTC	Medullary Peritubular Capillary Endothelial Cell
EC-V	Cortical Vein Endothelial Cell
EC-LYM	Lymphatic Endothelial Cell
EC-PCV	Postcapillary Venule Endothelial Cell`;

// Note: the rows after 'NA' are not sorted. There may be a better way to order these rows among the full set of table rows.
// For now, we just use the current ordering however.

export const abbrevRows = tsvParse(abbrevTsv);


export const abbrevToName = Object.fromEntries(
	abbrevRows.map(d => [d.clusterAbbrev, d.clusterName])
);

export const nameToAbbrev = Object.fromEntries(
	abbrevRows.map(d => [d.clusterName, d.clusterAbbrev])
);

export function abbrevToHumanReadable(abbrev) {
    if (abbrev.startsWith('C_M')) {
        abbrev = abbrev.replace('C_M', 'C/M');
    }
    if (abbrev.startsWith('C-int')) {
        // TODO: do we want to keep the -int here?
        abbrev = abbrev.replace('C-int', 'C-');
    }
    if(abbrev.includes('VSMC_P')) {
        return abbrev.replace('VSMC_P', 'VSMC/P');
    }
    if(abbrev === 'SC_NEU') {
        return 'SC/NEU';
    }
    if(abbrev.includes('_TEM_')) {
        return abbrev.replace('_TEM_', ' TEM/');
    }
    if(abbrev.includes('PT-S1_S2')) {
        return abbrev.replace('PT-S1_S2', 'PT-S1/S2');
    }
    return abbrev;
}