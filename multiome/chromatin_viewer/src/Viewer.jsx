/* eslint-disable */
import React, { Suspense, useRef, useState, useMemo, useCallback, useEffect } from 'react';

import { HiGlassComponent } from 'higlass';
import register from 'higlass-register';
import { BigwigDataFetcher } from 'higlass-bigwig-datafetcher';

import clsx from 'clsx';
import { makeStyles, useTheme } from '@material-ui/core/styles';
import Input from '@material-ui/core/Input';
import InputLabel from '@material-ui/core/InputLabel';
import MenuItem from '@material-ui/core/MenuItem';
import FormControl from '@material-ui/core/FormControl';
import ListItemText from '@material-ui/core/ListItemText';
import Select from '@material-ui/core/Select';
import Checkbox from '@material-ui/core/Checkbox';
import Chip from '@material-ui/core/Chip';
import Button from '@material-ui/core/Button';
import Grid from '@material-ui/core/Grid';
import TextField from '@material-ui/core/TextField';
import Autocomplete from '@material-ui/lab/Autocomplete';
import Typography from '@material-ui/core/Typography';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogContentText from '@material-ui/core/DialogContentText';
import DialogTitle from '@material-ui/core/DialogTitle';
import Snackbar from '@material-ui/core/Snackbar';
import IconButton from '@material-ui/core/IconButton';
import CloseIcon from '@material-ui/icons/Close';
import Fade from '@material-ui/core/Fade';

import { useQueryState, parseAsString } from 'nuqs';
import { abbrevRows, abbrevToName, nameToAbbrev, abbrevToHumanReadable } from './abbreviations.js';

register(
  { dataFetcher: BigwigDataFetcher, config: BigwigDataFetcher.config },
  { pluginType: "dataFetcher" }
);

// Lazy load the HiGlass React component,
// using a dynamic import.
// Note: the ref is not working as expected with React.lazy,
// so using non-lazy import instead.
/*
const LazyHiGlassComponent = React.lazy(async () => {
  const { HiGlassComponent } = await import('higlass');
  return { default: HiGlassComponent };
});
*/


const REFERENCE_TILESETS = {
  hg38: {
    chromosomes: 'NyITQvZsS_mOFNlz5C2LJg',
    genes: 'P0PLbQMwTYGy-5uPIQid7A',
  },
  hg19: {
    chromosomes: 'N12wVGG9SPiTkk03yUayUw',
    genes: 'OHJakQICQD6gTD7skx4EWA',
  },
  mm9: {
    chromosomes: 'WAVhNHYxQVueq6KulwgWiQ',
    genes: 'GUm5aBiLRCyz2PsBea7Yzg',
  },
  mm10: {
    chromosomes: 'EtrWT0VtScixmsmwFSd7zg',
    genes: 'QDutvmyiSrec5nX4pA5WGQ',
  },
};

// Reference: https://github.com/vitessce/vitessce/blob/b661a04380a7937fd8a3c576f0601e40922d9a0a/packages/view-types/genomic-profiles/src/GenomicProfilesSubscriber.js#L35
const REFERENCE_STATIC_FILES = {
  hg38: {
    chromosomes: 'https://raw.githubusercontent.com/vitessce/negspy/master/negspy/data/hg38/chromSizes.tsv',
  },
  hg19: {
    chromosomes: 'https://raw.githubusercontent.com/vitessce/negspy/master/negspy/data/hg19/chromSizes.tsv',
  },
  mm9: {
    chromosomes: 'https://raw.githubusercontent.com/vitessce/negspy/master/negspy/data/mm9/chromSizes.tsv',
  },
  mm10: {
    chromosomes: 'https://raw.githubusercontent.com/vitessce/negspy/master/negspy/data/mm10/chromSizes.tsv',
  },
};

const CHRS = [
	'chr1',
	'chr2',
	'chr3',
	'chr4',
	'chr5',
	'chr6',
	'chr7',
	'chr8',
	'chr9',
	'chr10',
	'chr11',
	'chr12',
	'chr13',
	'chr14',
	'chr15',
	'chr16',
	'chr17',
	'chr18',
	'chr19',
	'chr20',
	'chr21',
	'chr22',
	'chrX',
	'chrY',
];

const ITEM_HEIGHT = 48;
const ITEM_PADDING_TOP = 8;
const MenuProps = {
  PaperProps: {
    style: {
      //maxHeight: ITEM_HEIGHT * 4.5 + ITEM_PADDING_TOP,
      //width: 250,
    },
  },
};

const useStyles = makeStyles((theme) => ({
	formControl: {
		margin: theme.spacing(1),
		minWidth: 200,
		maxWidth: 'none',
	},
	buttonContainer: {
		display: 'inline-block'
	},
	button: {
		margin: theme.spacing(1)
	},
	chips: {
		display: 'flex',
		flexWrap: 'wrap',
	},
	chip: {
		margin: 2,
		height: 'auto',
	},
	container: {
		display: 'flex',
		flexDirection: 'column',
		height: 'calc(100%)'
	},
	gridContainer: {
		flexGrow: 0,
	},
	higlassContainer: {
		flexGrow: 1,
		flexBasis: 1,
		height: 'calc(100%)'
	},
}));


// From https://personal.sron.nl/~pault/#sec:qualitative
const PALETTE = [
  [68, 119, 170],
  [136, 204, 238],
  [68, 170, 153],
  [17, 119, 51],
  [153, 153, 51],
  [221, 204, 119],
  [204, 102, 119],
  [136, 34, 85],
  [170, 68, 153],
];


const baseUrlV3 = 'https://data-2.vitessce.io/kpmp-atlas-v2/transepinome/bigwigs_from_blue_20251110';
//const baseUrlV3 = 'http://localhost:3004';

const subclass_l3_from_table = abbrevRows.map(d => d.clusterAbbrev);

let subclass_l3 = ['ATL', 'B', 'C-EC-PTC', 'C-FIB', 'C-FIB-OSMRhi', 'C-FIB-OSMRlo', 'C-FIB-PATH', 'C-MYOF', 'C-TAL-A', 'C-TAL-B', 'CCD-IC-A', 'CCD-PC', 'CD8+_TEM_TEMRA', 'CD8+_TEM_TRM', 'CNT', 'CNT-PC', 'C_M-FIB', 'C_M-TAL-A', 'C_M-TAL-B', 'DCT1', 'DCT2', 'DTL1', 'DTL2', 'DTL3', 'EC-AA', 'EC-AVR', 'EC-DVR', 'EC-EA', 'EC-GC', 'EC-GC-FILIP1+', 'EC-LYM', 'EC-PCV', 'EC-V', 'ERY', 'IC-B', 'ILC3', 'IM-FIB', 'IM-pvMYOF', 'IMCD', 'M-EC-PTC', 'M-TAL', 'M-VSMC_P', 'MAIT', 'MAST', 'MC', 'MD', 'MON', 'NK', 'Naive_Th', 'OM-FIB', 'OMCD-IC-A', 'OMCD-PC', 'PEC', 'PL', 'POD', 'PT-S1', 'PT-S2', 'PT-S3', 'PapE', 'REN', 'SC_NEU', 'T-REG', 'VSMC', 'VSMC_P', 'aCNT', 'aDCT', 'aDTL2', 'aEC-GC', 'aPT-S1_S2', 'aPT-S3', 'aPT1', 'aPT2', 'aTAL1', 'aTAL2', 'angEC-PTC', 'cDC1', 'cDC2', 'cycEC', 'cycMAC', 'cycPT', 'cycT', 'cycTAL', 'dATL', 'dC-TAL', 'dCCD-IC-A', 'dCNT', 'dCNT-PC', 'dDCT', 'dDTL', 'dEC-AVR', 'dEC-GC', 'dEC-PTC', 'dFIB', 'dIM-FIB', 'dIMCD', 'dM-TAL', 'dOM-FIB', 'dOMCD-IC-A', 'dOMCD-PC', 'dPOD', 'dPT', 'dPT-S1', 'dPT-S3', 'dVSMC', 'frDCT', 'frPT-S1_S2', 'frPT-S3', 'frTAL', 'iaEC-AVR', 'iaEC-PTC', 'infEC-PTC', 'mDC', 'moFAM', 'moMAC-C3+', 'moMAC-CXCL10+', 'moMAC-HBEGF+', 'ncMON', 'pDC', 'pvFIB', 'pvFIB-PI16+', 'pvFIB-RSPO3+', 'pvMYOF', 'resMAC-HLAIIhi_', 'resMAC-LYVE1+', 'tOMCD-PC-IC']
	//.toSorted((a, b) => a.toLowerCase().localeCompare(b.toLowerCase()));
	.toSorted((a, b) => {
		const aIndex = subclass_l3_from_table.indexOf(a);
		const bIndex = subclass_l3_from_table.indexOf(b);
		return aIndex - bIndex;
	});

const SUBTYPE_TO_BIGWIG_URL = Object.fromEntries(
	subclass_l3.map(subType => {
		const url = `${baseUrlV3}/${encodeURIComponent(subType)}.RPGCnorm.bigwig`;
		return [subType, url];
	})
);


const technologyToBigWigUrl = {
	/* V1
	'Methylation levels (WGBS): microdissected glomeruli (GLOM)': `${baseUrl}/mlGLOM.bw`,
	'Methylation levels (WGBS): microdissected tubulointerstitium (TI)': `${baseUrl}/mlTI.bw`,
	'Histone modifications (CUT&RUN): H3K27Ac (active chromatin)': `${baseUrl}/H3K27Ac.20.Rep1_Rep2_Rep3.KTRC.A_B_D.merge.norm.bw`,
	'Histone modifications (CUT&RUN): H3K4me1 (active chromatin)': `${baseUrl}/H3K4me1.20.Rep1_2_3.merge.norm.bw`,
	'Histone modifications (CUT&RUN): H3K4me3 (active chromatin)': `${baseUrl}/H3K4me3.20.Rep1_2_3.merge.norm.bw`,
	'Histone modifications (CUT&RUN): H3K27me3 (repressive chromatin)': `${baseUrl}/H3K27me3.20.Rep1_Rep2_Rep3.687_688_LN2_OCT.KTRC.A_B_D.merge.norm.bw`,
	*/
};

const technologyToColor = {
	/* V1
	'Methylation levels (WGBS): microdissected glomeruli (GLOM)': [0x11, 0x38, 0x91],
	'Methylation levels (WGBS): microdissected tubulointerstitium (TI)': [0xbb, 0x00, 0x1e],
	'Histone modifications (CUT&RUN): H3K27Ac (active chromatin)': [0x00, 0x69, 0x50],
	'Histone modifications (CUT&RUN): H3K4me1 (active chromatin)': [0x00, 0xa6, 0xdd],
	'Histone modifications (CUT&RUN): H3K4me3 (active chromatin)': [0x00, 0xb1, 0x2a],
	'Histone modifications (CUT&RUN): H3K27me3 (repressive chromatin)': [0x80, 0x17, 0x37],
	*/
};

const ALL_TECHNOLOGIES = Object.keys(technologyToBigWigUrl); 

const cellTypeClassToBigWigUrl = {
	/* V1
	'aPT': `${baseUrl}/Multiome_paper/aPT.bw`,
	'aTAL12': `${baseUrl}/Multiome_paper/aTAL12.bw`,
	'C-TAL': `${baseUrl}/Multiome_paper/C-TAL.bw`,
	'POD': `${baseUrl}/Multiome_paper/POD.bw`,
	'PT-S1': `${baseUrl}/Multiome_paper/PT-S1.bw`,
	'PT-S12': `${baseUrl}/Multiome_paper/PT-S12.bw`,
	'PT-S2': `${baseUrl}/Multiome_paper/PT-S2.bw`,
	'PT-S3': `${baseUrl}/Multiome_paper/PT-S3.bw`,
	*/
	...SUBTYPE_TO_BIGWIG_URL,

};

const cellTypeSubclassToBigWigUrl = {
	/* V1
	'ATL (subclass)': `${baseUrl}/Multiome_subclassl1/ATL_subclass.l1.bw`,
	'CNT (subclass)': `${baseUrl}/Multiome_subclassl1/CNT_subclass.l1.bw`,
	'DCT (subclass)': `${baseUrl}/Multiome_subclassl1/DCT_subclass.l1.bw`,
	'DTL (subclass)': `${baseUrl}/Multiome_subclassl1/DTL_subclass.l1.bw`,
	'EC (subclass)': `${baseUrl}/Multiome_subclassl1/EC_subclass.l1.bw`,
	'FIB (subclass)': `${baseUrl}/Multiome_subclassl1/FIB_subclass.l1.bw`,
	'IC (subclass)': `${baseUrl}/Multiome_subclassl1/IC_subclass.l1.bw`,
	'IMM (subclass)': `${baseUrl}/Multiome_subclassl1/IMM_subclass.l1.bw`,
	'NEU (subclass)': `${baseUrl}/Multiome_subclassl1/NEU_subclass.l1.bw`,
	'PapE (subclass)': `${baseUrl}/Multiome_subclassl1/PapE_subclass.l1.bw`,
	'PC (subclass)': `${baseUrl}/Multiome_subclassl1/PC_subclass.l1.bw`,
	'PEC (subclass)': `${baseUrl}/Multiome_subclassl1/PEC_subclass.l1.bw`,
	'POD (subclass)': `${baseUrl}/Multiome_subclassl1/POD_subclass.l1.bw`,
	'PT (subclass)': `${baseUrl}/Multiome_subclassl1/PT_subclass.l1.bw`,
	'TAL (subclass)': `${baseUrl}/Multiome_subclassl1/TAL_subclass.l1.bw`,
	'VSM-P (subclass)': `${baseUrl}/Multiome_subclassl1/VSM-P_subclass.l1.bw`,
	*/
};

const cellTypeToBigWigUrl = {
	...cellTypeClassToBigWigUrl,
	...cellTypeSubclassToBigWigUrl,
};


const ALL_REFERENCE_TRACKS = Object.keys(technologyToBigWigUrl);
const ALL_CELLTYPES = Object.keys(cellTypeToBigWigUrl);
const ALL_COARSE = Object.keys(cellTypeClassToBigWigUrl);
const ALL_FINE = Object.keys(cellTypeSubclassToBigWigUrl);


const hgOptionsBase = {
	sizeMode: "bounded", // Stretch the height of HiGlass to its container <div/>
	pixelPreciseMarginPadding: false,
	bounded: true,
	containerPaddingX: 0,
	containerPaddingY: 0,
	viewMarginTop: 0,
	viewMarginBottom: 0,
	viewMarginLeft: 0,
	viewMarginRight: 0,
	viewPaddingTop: 0,
	viewPaddingBottom: 0,
	viewPaddingLeft: 0,
	viewPaddingRight: 0,
};

const initialXDomain = [
	0,
	300000000
];

// Map from a lock name (arbitrary) to a list of tracks (assumed to be in the "main" view).
const lockGroups = {
	"methylation": [
		"Methylation levels (WGBS): microdissected glomeruli (GLOM)",
		"Methylation levels (WGBS): microdissected tubulointerstitium (TI)"
	],
	/*"active-chromatin": [
		'Histone modifications (CUT&RUN): H3K27Ac (active chromatin)',
		'Histone modifications (CUT&RUN): H3K4me1 (active chromatin)',
		'Histone modifications (CUT&RUN): H3K4me3 (active chromatin)',
	],*/
	// only one track for repressive chromatin, so no need for lock.
	"coarse": ALL_COARSE,
	"fine": ALL_FINE,
};

export default function Viewer(props) {
	const {
		assembly = 'hg38',
		higlassServer = 'https://higlass.io/api/v1',
		theme = 'light',
	} = props;
	
	const hgRef = useRef();
	const containerRef = useRef(); // Used for passing anchorEl to MUI components.
	const [viewConfLoaded, setViewConfLoaded] = useState(false);
	const [gene, setGene] = useQueryState('gene', parseAsString);

	// TODO: store selected cell types in URL query params as well?
	// Or in browser local storage?

	const [toastOpen, setToastOpen] = useState(false);

	const handleToastClose = (event, reason) => {
		if (reason === 'clickaway') {
			return;
		}

		setToastOpen(false);
	};


	useEffect(() => {
		if(hgRef.current && gene) {
			const viewId = 'main';
			try {
				hgRef.current.api.zoomToGene(viewId, gene, 10_000);
				console.log('Zoomed to gene:', gene);
				setToastOpen(true);
			} catch (error) {
				console.error('Error zooming to gene:', error);
			}
		}
	}, [hgRef, gene, viewConfLoaded]);

	const [selectedCellTypes, setSelectedCellTypes] = useState(['aPT2', 'aPT1', 'aPT-S3', 'frPT-S3']);

	function handleChange(event) {
		setSelectedCellTypes(event.target.value)
	}

	const classes = useStyles();

	const hgViewConfig = useMemo(() => {
		const profileTrackHeight = 40;

		// Set up the colors to use in the HiGlass view config based on the current theme.
		const foregroundColor = (theme === 'dark' ? '#C0C0C0' : '#000000');
		const backgroundColor = (theme === 'dark' ? '#000000' : '#fff');

		// Define the "reference tracks" for chromosome labels and gene annotations.
		const referenceTracks = [
			{
				type: 'horizontal-chromosome-labels',
				server: higlassServer,
				chromInfoPath: REFERENCE_STATIC_FILES[assembly].chromosomes,
				uid: 'chromosome-labels',
				options: {
					color: foregroundColor,
					fontSize: 12,
					fontIsLeftAligned: false,
					showMousePosition: true,
					mousePositionColor: foregroundColor,
				},
				height: 30,
			},
			{
				type: 'horizontal-gene-annotations',
				server: higlassServer,
				tilesetUid: REFERENCE_TILESETS[assembly].genes,
				uid: 'gene-annotations',
				options: {
					name: 'Gene Annotations (hg38)',
					fontSize: 10,
					labelPosition: 'hidden',
					labelLeftMargin: 0,
					labelRightMargin: 0,
					labelTopMargin: 0,
					labelBottomMargin: 0,
					minHeight: 24,
					geneAnnotationHeight: 16,
					geneLabelPosition: 'outside',
					geneStrandSpacing: 4,
					showMousePosition: true,
					mousePositionColor: foregroundColor,
					plusStrandColor: foregroundColor,
					minusStrandColor: foregroundColor,
					labelColor: 'black',
					labelBackgroundColor: backgroundColor,
					trackBorderWidth: 0,
					trackBorderColor: 'black',
				},
				height: 70,
			},
		];

		const mainTracks = ALL_TECHNOLOGIES.map((technology) => {
			// Get the uid for the HiGlass track.
			const trackUid = technology;
			// Create the HiGlass track definition for this profile.

			const color = technologyToColor[technology];
			const track = {
				type: 'horizontal-bar',
				uid: trackUid,
				data: {
					type: 'bbi',
					url: technologyToBigWigUrl[technology],
					chromSizesUrl: REFERENCE_STATIC_FILES[assembly].chromosomes,
				},
				options: {
					name: technology,
					showMousePosition: true,
					mousePositionColor: foregroundColor,
					labelColor: (theme === 'dark' ? 'white' : 'black'),
					labelBackgroundColor: (theme === 'dark' ? 'black' : 'white'),
					labelShowAssembly: false,
					barFillColor: `rgb(${color[0]},${color[1]},${color[2]})`,
				},
				height: profileTrackHeight,
			};
			return track;
		});

		const profileTracks = selectedCellTypes.map((cellType) => {
			// Get the uid for the HiGlass track.
			const trackUid = cellType;
			const fullCellTypeName = abbrevToName[cellType];
			const color = PALETTE[ALL_CELLTYPES.indexOf(cellType) % PALETTE.length];
			// Create the HiGlass track definition for this profile.
			const track = {
			type: 'horizontal-bar',
			uid: trackUid,
			data: {
				type: 'bbi',
				url: cellTypeToBigWigUrl[cellType],
				chromSizesUrl: REFERENCE_STATIC_FILES[assembly].chromosomes,
			},
			options: {
				name: ` ${abbrevToHumanReadable(cellType)}: ${fullCellTypeName}`,
				showMousePosition: true,
				mousePositionColor: foregroundColor,
				labelColor: (theme === 'dark' ? 'white' : 'black'),
				labelBackgroundColor: (theme === 'dark' ? 'black' : 'white'),
				labelShowAssembly: false,
				barFillColor: `rgb(${color[0]},${color[1]},${color[2]})`,
			},
			height: profileTrackHeight,
			};
			return track;
		});

		const hgView = {
			tracks: {
				top: [
					...referenceTracks,
					...mainTracks,
					...profileTracks,
				],
				left: [],
				center: [],
				right: [],
				bottom: [],
				whole: [],
				gallery: [],
			},
			layout: {
				w: 12,
				h: 12,
				x: 0,
				y: 0,
				static: false,
			},
		};
		return hgView;
	}, [selectedCellTypes, theme, assembly]);

	const hgFullConfig = useMemo(() => {
		// Construct locks for value scales.
		const visibleTrackUids = [...ALL_REFERENCE_TRACKS, ...selectedCellTypes];
		const locksByViewUid = {};
		Object.entries(lockGroups).forEach(([lockName, trackUids]) => {
			trackUids.forEach(trackUid => {
				if(visibleTrackUids.includes(trackUid)) {
					locksByViewUid[`main.${trackUid}`] = lockName;
				}
			});
		});
		const locksDict = {};
		Object.entries(lockGroups).forEach(([lockName, trackUids]) => {
			locksDict[lockName] = { uid: lockName };
			trackUids.forEach(trackUid => {
				if(visibleTrackUids.includes(trackUid)) {
					locksDict[lockName][`main.${trackUid}`] = {
						view: "main",
						track: trackUid,
					};
				}
			});
		});

		// Return view config.
		return {
			editable: true,
			zoomFixed: false,
			trackSourceServers: [
				'//higlass.io/api/v1',
			],
			exportViewUrl: '//higlass.io/api/v1/viewconfs',
			views: [
				{
					uid: 'main',
					...hgViewConfig,
					initialXDomain,
					initialYDomain: initialXDomain,
					"autocompleteSource": `/api/v1/suggest/?d=${REFERENCE_TILESETS[assembly].genes}&`,
					"genomePositionSearchBox": {
						"autocompleteServer": "//higlass.io/api/v1",
						"autocompleteId": REFERENCE_TILESETS[assembly].genes,
						"chromInfoServer": "https://higlass.io/api/v1",
						"chromInfoId": assembly,
						"visible": true
					},
					"chromInfoPath": REFERENCE_STATIC_FILES[assembly].chromosomes,
				},
			],
			zoomLocks: {
				locksByViewUid: {},
				locksDict: {},
			},
			locationLocks: {
				locksByViewUid: {},
				locksDict: {},
			},
			// Reference: https://github.com/higlass/higlass/blob/b2ee5940c519982dc53685153ff863d64443d0bb/docs/examples/viewconfs/lots-of-locks.json#L158
			valueScaleLocks: {
				locksByViewUid,
				locksDict,
			},
		};
	}, [hgViewConfig, assembly]);

	const [acronymsOpen, setAcronymsOpen] = useState(false);

	function onAllCoarse() {
		setSelectedCellTypes(ALL_COARSE);
	}
	function onAllFine() {
		setSelectedCellTypes(ALL_FINE);
	}
	function onClearAll() {
		setSelectedCellTypes([]);
	}
	function onOpenAcronyms() {
		setAcronymsOpen(true);
	}
	function onCloseAcronyms() {
		setAcronymsOpen(false);
	}

	const hgOptions = useMemo(() => ({
		...hgOptionsBase,
		onViewConfLoaded: () => {
			setViewConfLoaded(true);
		},
	}), []);

	// Listen for higlass view config changes.
	useEffect(() => {
		hgRef.current.api.on('viewConfig', newViewConfigString => {
			// Note: the setViewConfLoaded in hgOptions is not being called as expected,
			// so using this callback instead.
			const hasMainView = JSON.parse(newViewConfigString).views.some(v => v.uid === 'main');
			if(hasMainView) {
				setViewConfLoaded(true);
			}
		});

		return () => hgRef.current.api.off('viewConfig');
	}, [hgRef]);

	return (
		<div className={classes.container}>
			<div className={classes.gridContainer} style={{ maxHeight: '130px', overflowY: 'scroll' }} ref={containerRef}>
				<Grid
					container
					direction="row"
					justifyContent="start"
				>
					<Grid item xs={6} style={{ padding: '4px'}}>
						<p>KPMP Atlas V2: Chromatin viewer&nbsp;
						<div style={{ display: 'inline'}}>
							{/* Acronym dialog */}
							<Button size="small" variant="outlined" onClick={onOpenAcronyms}>
								About
							</Button>
							<Dialog
								open={acronymsOpen}
								onClose={onCloseAcronyms}
								aria-labelledby="alert-dialog-title"
								aria-describedby="alert-dialog-description"
							>
								<DialogContent>
									<DialogContentText id="alert-dialog-description">
										<Typography variant="h5">About</Typography>
										This viewer displays chromatin accessibility profiles per cell type (at subclass level 3) that have been generated for the KPMP Atlas V2.
										The profiles are derived from MACS peak calling output from 10x Multiome experiments performed by the transepinome working group within KPMP.
										Associated fragment pileups were saved as bedGraph files, then converted to bigWig using the bedGraphToBigWig tool.
										<br/>
										<br/>
										<a href="https://doi.org/10.1101/2025.09.26.678707" target="_blank" style={{ fontWeight: 'bold' }}>KPMP Atlas V2 Preprint on bioRxiv</a><br/>
										<a href="https://github.com/KPMP-Scientific/KPMP-Atlas-v2" target="_blank">KPMP Atlas V2 GitHub repository</a><br/>
										{/*<Typography variant="h5">Publication</Typography>
										
										Gisch, D. L., Brennan M., Lake B. B., Basta J. et al. "The chromatin landscape of healthy and injured cell types in the human kidney." <b>bioRxiv</b> (2023).
										*/}
										{/*<Typography variant="h5">Cell type acronyms</Typography>
										TAL: thick ascending limb<br/>
										aTAL: adaptive TAL<br/>
										C-TAL: cortical TAL<br/>
										POD: podocyte<br/>
										PT: proximal tubule<br/>
										aPT: adaptive PT<br/>
										PT-S1: PT segment 1<br/>
										PT-S2: PT segment 2<br/>
										PT-S3: PT segment 3<br/>
										PT-S12: PT-S1 merged with PT-S2<br/>*/}
									</DialogContentText>
								</DialogContent>
								<DialogActions>
								<Button onClick={onCloseAcronyms} color="primary" autoFocus>
									Close
								</Button>
								</DialogActions>
							</Dialog>
						</div>	
						<br/><br/>
						Select cell types from the list to the right to visualize corresponding aggregate snATAC-seq tracks.</p>
					</Grid>
					<Grid item xs={4}>
						<FormControl className={classes.formControl}>
							<InputLabel id="demo-mutiple-chip-label">Cell Types</InputLabel>
							<Select
								labelId="demo-mutiple-chip-label"
								id="demo-mutiple-chip"
								multiple
								value={selectedCellTypes}
								onChange={handleChange}
								input={<Input id="select-multiple-chip" />}
								renderValue={(selected) => (
									<div className={classes.chips}>
										{selected.map((value) => (
											<Chip key={value} label={value} className={classes.chip} />
										))}
									</div>
								)}
								MenuProps={{
									...MenuProps,
									container: containerRef.current,
								}}
							>
								{ALL_CELLTYPES.map((abbrev) => (
									<MenuItem key={abbrev} value={abbrev} style={{ fontWeight: 'normal', flexDirection: 'column', alignItems: 'start', paddingTop: '3px', paddingBottom: '3px' }}>
										<span>{abbrevToHumanReadable(abbrev)}</span>
										<span style={{ fontSize: '11px'}}>{abbrevToName[abbrev]}</span>
									</MenuItem>
								))}
							</Select>
						</FormControl>
					</Grid>
					<Grid item xs={2}>
						<FormControl className={classes.buttonContainer}>
							{/*<Button size="small" variant="outlined" className={classes.button} onClick={onAllCoarse}>All classes</Button>*/}
							{/*<Button size="small" variant="outlined" className={classes.button} onClick={onAllFine}>All subclasses</Button>*/}
							<Button size="small" variant="outlined" className={classes.button} onClick={onClearAll}>Clear cell type tracks</Button>
						</FormControl>
					</Grid>
				</Grid>
			</div>
			<div className={classes.higlassContainer}>
				<HiGlassComponent
					ref={hgRef}
					zoomFixed={false}
					viewConfig={hgFullConfig}
					options={hgOptions}
				/>
				{/*<Suspense fallback={<div>Loading...</div>}>
					<LazyHiGlassComponent
						ref={hgRef}
						zoomFixed={false}
						viewConfig={hgFullConfig}
						options={hgOptions}
					/>
				</Suspense>*/}
			</div>
			<div>
				<Snackbar
					anchorOrigin={{
						vertical: 'bottom',
						horizontal: 'right',
					}}
					open={toastOpen}
					autoHideDuration={5000}
					onClose={handleToastClose}
					message={`Zoomed to gene ${gene}`}
					action={
						<IconButton size="small" aria-label="close" color="inherit" onClick={handleToastClose}>
							<CloseIcon fontSize="small" />
						</IconButton>
					}
					TransitionComponent={Fade}
				/>
			</div>
		</div>
	);
}

