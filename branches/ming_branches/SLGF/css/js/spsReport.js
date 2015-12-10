////////////////////////////////////////////////////////////////////////////////
// spsReport dynamic generation
//
// By jcanhita@eng.ucsd.edu on nov/2011
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Constants
////////////////////////////////////////////////////////////////////////////////

// report entry point & internal components
var INITIAL_PAGE = "index.html";

var PAGE_FIELD_CELLS_PER_LINE = "cellsPerLine";
var PAGE_FIELD_SERVER         = "serverLocation";
var PAGE_FIELD_PROJECT_DIR    = "projectDir";
var PAGE_FIELD_FASTA_FILENAME = "FastaFilename";


// table ID for reportServer
var TABLE_PROTEIN_ID  = 2;
var TABLE_CONTIG_ID   = 3;
var TABLE_CLUSTER_ID  = 4;
var TABLE_SPECTRA_ID  = 5;
var TABLE_COVERAGE_ID = 6;

// table in field separators.
var TABLE_SEP_L1 = "|";
var TABLE_SEP_L2 = "@";
var TABLE_SEP_L3 = "&";
var TABLE_SEP_L4 = "!";

// IDs for objects
var PEP_ELEM_PREFIX = "PEP";  // peptide element in protein coverage page
var INP_ELEM_PREFIX = "INP";  // peptide input element in protein coverage page

// CSV file separator
var CSV_SEP = ";";

// number of elements per page
var PAGE_LENGTH = 10;

// Cells per line in protein coverage reports pages
var CELLS_PER_LINE = 20;

// Amino acids per line in protein sequence, in proteins page
var AAS_PER_LINE = 50;

// Amino acids group size in protein sequence, in proteins page
var AAS_GROUP_SIZE = 10;

// Server location.
var SERVER = ""; //"http://usa.ucsd.edu:8080/cgi-bin/";

// define if the image cache is used
var USE_IMAGE_CACHE = true;

// Project directory on server.
var PROJECT_DIR = ".";

// Fasta filename in project directory subtree.
var GLOBAL_FASTA_FILENAME = "";

// Status file filename
var STATUS_FILENAME = "status.txt";

// relauch script name
var RELAUNCH_SCRIPT = "relauncher.sh";

// status to keep pooling
var STATUS_RETRIEVE = "Running";

// status to allow navigation
var STATUS_OK = "Finished";

// tag delimiters
var TAG_OPEN  = '<';
var TAG_CLOSE = '>';

////////////////////////////////////////////////////////////////////////////////
// Report table column headers

var REPORT_HEADER_PROTEINS_0  = "Protein";
var REPORT_HEADER_PROTEINS_1  = "Description";
var REPORT_HEADER_PROTEINS_2  = "Contigs";
var REPORT_HEADER_PROTEINS_3  = "Spectra";
var REPORT_HEADER_PROTEINS_4  = "Amino acids";
var REPORT_HEADER_PROTEINS_5  = "Coverage (%)";

var REPORT_HEADER_CONTIGS_0   = "Index";
var REPORT_HEADER_CONTIGS_1   = "Spectra";
var REPORT_HEADER_CONTIGS_2   = "Contig";
var REPORT_HEADER_CONTIGS_3   = "Contig sequence";
var REPORT_HEADER_CONTIGS_4   = "Protein";

var REPORT_HEADER_CLUSTERS_0  = "Index";
var REPORT_HEADER_CLUSTERS_1  = "Spectrum";
var REPORT_HEADER_CLUSTERS_2  = "Peptide";
var REPORT_HEADER_CLUSTERS_3  = "Mass (m)";
var REPORT_HEADER_CLUSTERS_4  = "Charge (z)";
var REPORT_HEADER_CLUSTERS_5  = "B (%)";
var REPORT_HEADER_CLUSTERS_6  = "Y (%)";
var REPORT_HEADER_CLUSTERS_7  = "BY Intensity (%)";

var REPORT_HEADER_SPECTRA_0  = "Index";
var REPORT_HEADER_SPECTRA_1  = "Scan";
var REPORT_HEADER_SPECTRA_2  = "Spectrum";
var REPORT_HEADER_SPECTRA_3  = "Sequences";
var REPORT_HEADER_SPECTRA_4  = "Protein";
var REPORT_HEADER_SPECTRA_5  = "Mass (m)";
var REPORT_HEADER_SPECTRA_6  = "Charge (z)";
var REPORT_HEADER_SPECTRA_7  = "B (%)";
var REPORT_HEADER_SPECTRA_8  = "Y (%)";
var REPORT_HEADER_SPECTRA_9  = "BY Intensity (%)";

// Sequence names
var REPORT_SEQ_NAME_REFERENCE = "Reference";
var REPORT_SEQ_NAME_HOMOLOG   = "Homolog";
var REPORT_SEQ_NAME_DENOVO    = "de Novo";
var REPORT_SEQ_NAME_USER      = "User";

// buttons
var REPORT_BUTTON_UPDATE      = "Update";
var REPORT_BUTTON_ALIGN       = "Re-align";

////////////////////////////////////////////////////////////////////////////////
// Table fields' contents

// Proteins table
var TABLE_PROTEINS_FIELD_ID       = 0;
var TABLE_PROTEINS_FIELD_NAME     = 1;
var TABLE_PROTEINS_FIELD_DESC     = 2;
var TABLE_PROTEINS_FIELD_CONTIGS  = 3;
var TABLE_PROTEINS_FIELD_SPECTRA  = 4;
var TABLE_PROTEINS_FIELD_AAS      = 5;
var TABLE_PROTEINS_FIELD_COVERAGE = 6;
var TABLE_PROTEINS_FIELD_SEQUENCE = 7;

var TAG_TABLE_PROTEINS_FIELD_ID       = makeTag(TABLE_PROTEINS_FIELD_ID);
var TAG_TABLE_PROTEINS_FIELD_NAME     = makeTag(TABLE_PROTEINS_FIELD_NAME);
var TAG_TABLE_PROTEINS_FIELD_DESC     = makeTag(TABLE_PROTEINS_FIELD_DESC);
var TAG_TABLE_PROTEINS_FIELD_CONTIGS  = makeTag(TABLE_PROTEINS_FIELD_CONTIGS);
var TAG_TABLE_PROTEINS_FIELD_SPECTRA  = makeTag(TABLE_PROTEINS_FIELD_SPECTRA);
var TAG_TABLE_PROTEINS_FIELD_AAS      = makeTag(TABLE_PROTEINS_FIELD_AAS);
var TAG_TABLE_PROTEINS_FIELD_COVERAGE = makeTag(TABLE_PROTEINS_FIELD_COVERAGE);


// Contigs table
var TABLE_CONTIGS_FIELD_ID            = 0;
var TABLE_CONTIGS_FIELD_PROTEIN       = 1;
var TABLE_CONTIGS_FIELD_SPECTRA       = 2;
var TABLE_CONTIGS_FIELD_SEQ_REF       = 3;
var TABLE_CONTIGS_FIELD_SEQ_HOM       = 4;
var TABLE_CONTIGS_FIELD_SEQ_NOVO      = 5;
var TABLE_CONTIGS_FIELD_SEQ_USER      = 6;
var TABLE_CONTIGS_FIELD_PROTEIN_NAME  = 7;
var TABLE_CONTIGS_FIELD_PROTEIN_DESC  = 8;
var TABLE_CONTIGS_FIELD_MASS_REF      = 10;
var TABLE_CONTIGS_FIELD_MASS_HOM      = 11;
var TABLE_CONTIGS_FIELD_OFF_REF       = 12;
var TABLE_CONTIGS_FIELD_OFF_HOM       = 13;
var TABLE_CONTIGS_FIELD_REVERSE       = 14;
var TABLE_CONTIGS_FIELD_FILE_ABRUIJN  = 15;
var TABLE_CONTIGS_FIELD_FILE_STARS    = 16;
var TABLE_CONTIGS_FIELD_FILE_SEQS     = 17;

var TAG_TABLE_CONTIGS_FIELD_ID            = makeTag(TABLE_CONTIGS_FIELD_ID);
var TAG_TABLE_CONTIGS_FIELD_PROTEIN       = makeTag(TABLE_CONTIGS_FIELD_PROTEIN);
var TAG_TABLE_CONTIGS_FIELD_SPECTRA       = makeTag(TABLE_CONTIGS_FIELD_SPECTRA);
var TAG_TABLE_CONTIGS_FIELD_SEQ_REF       = makeTag(TABLE_CONTIGS_FIELD_SEQ_REF);
var TAG_TABLE_CONTIGS_FIELD_SEQ_HOM       = makeTag(TABLE_CONTIGS_FIELD_SEQ_HOM);
var TAG_TABLE_CONTIGS_FIELD_SEQ_NOVO      = makeTag(TABLE_CONTIGS_FIELD_SEQ_NOVO);
var TAG_TABLE_CONTIGS_FIELD_SEQ_USER      = makeTag(TABLE_CONTIGS_FIELD_SEQ_USER);
var TAG_TABLE_CONTIGS_FIELD_PROTEIN_NAME  = makeTag(TABLE_CONTIGS_FIELD_PROTEIN_NAME);
var TAG_TABLE_CONTIGS_FIELD_PROTEIN_DESC  = makeTag(TABLE_CONTIGS_FIELD_PROTEIN_DESC);
var TAG_TABLE_CONTIGS_FIELD_MASS_REF      = makeTag(TABLE_CONTIGS_FIELD_MASS_REF);
var TAG_TABLE_CONTIGS_FIELD_MASS_HOM      = makeTag(TABLE_CONTIGS_FIELD_MASS_HOM);
var TAG_TABLE_CONTIGS_FIELD_OFF_REF       = makeTag(TABLE_CONTIGS_FIELD_OFF_REF);
var TAG_TABLE_CONTIGS_FIELD_OFF_HOM       = makeTag(TABLE_CONTIGS_FIELD_OFF_HOM);
var TAG_TABLE_CONTIGS_FIELD_REVERSE       = makeTag(TABLE_CONTIGS_FIELD_REVERSE);
var TAG_TABLE_CONTIGS_FIELD_FILE_ABRUIJN  = makeTag(TABLE_CONTIGS_FIELD_FILE_ABRUIJN);
var TAG_TABLE_CONTIGS_FIELD_FILE_STARS    = makeTag(TABLE_CONTIGS_FIELD_FILE_STARS);
var TAG_TABLE_CONTIGS_FIELD_FILE_SEQS     = makeTag(TABLE_CONTIGS_FIELD_FILE_SEQS);

// Cluster table
var TABLE_CLUSTER_FIELD_ID            = 0;
var TABLE_CLUSTER_FIELD_CONTIG        = 1;
var TABLE_CLUSTER_FIELD_PROTEIN       = 2;
var TABLE_CLUSTER_FIELD_REFERENCE     = 3;
var TABLE_CLUSTER_FIELD_HOMOLOG       = 4;
var TABLE_CLUSTER_FIELD_DENOVO        = 5;
var TABLE_CLUSTER_FIELD_USER          = 6;
var TABLE_CLUSTER_FIELD_MASS          = 7;
var TABLE_CLUSTER_FIELD_CHARGE        = 8;
var TABLE_CLUSTER_FIELD_B_PER         = 9;
var TABLE_CLUSTER_FIELD_Y_PER         = 10;
var TABLE_CLUSTER_FIELD_BY_INT        = 11;
var TABLE_CLUSTER_FIELD_PROTEIN_NAME  = 12;

var TAG_TABLE_CLUSTER_FIELD_ID        = makeTag(TABLE_CLUSTER_FIELD_ID);
var TAG_TABLE_CLUSTER_FIELD_CONTIG    = makeTag(TABLE_CLUSTER_FIELD_CONTIG);
var TAG_TABLE_CLUSTER_FIELD_REFERENCE = makeTag(TABLE_CLUSTER_FIELD_REFERENCE);
var TAG_TABLE_CLUSTER_FIELD_HOMOLOG   = makeTag(TABLE_CLUSTER_FIELD_HOMOLOG);
var TAG_TABLE_CLUSTER_FIELD_DENOVO    = makeTag(TABLE_CLUSTER_FIELD_DENOVO);
var TAG_TABLE_CLUSTER_FIELD_USER      = makeTag(TABLE_CLUSTER_FIELD_USER);
var TAG_TABLE_CLUSTER_FIELD_MASS      = makeTag(TABLE_CLUSTER_FIELD_MASS);
var TAG_TABLE_CLUSTER_FIELD_CHARGE    = makeTag(TABLE_CLUSTER_FIELD_CHARGE);
var TAG_TABLE_CLUSTER_FIELD_B_PER     = makeTag(TABLE_CLUSTER_FIELD_B_PER);
var TAG_TABLE_CLUSTER_FIELD_Y_PER     = makeTag(TABLE_CLUSTER_FIELD_Y_PER);
var TAG_TABLE_CLUSTER_FIELD_BY_INT    = makeTag(TABLE_CLUSTER_FIELD_BY_INT);

var TAG_CLUSTER_PEPTIDE_ALL  = makeTag(TABLE_CLUSTER_FIELD_REFERENCE + '|' + TABLE_CLUSTER_FIELD_HOMOLOG + '|' + TABLE_CLUSTER_FIELD_DENOVO);

// Spectra table
var TABLE_SPECTRA_FIELD_ID              = 0;
var TABLE_SPECTRA_FIELD_IDX             = 1;
var TABLE_SPECTRA_FIELD_SCAN            = 2;
var TABLE_SPECTRA_FIELD_CLUSTER         = 3;
var TABLE_SPECTRA_FIELD_PROTEIN_NAME    = 4;
var TABLE_SPECTRA_FIELD_PKLBIN_IDX      = 5;
var TABLE_SPECTRA_FIELD_PKLBIN_FILENAME = 6;
var TABLE_SPECTRA_FIELD_SEQ_REFERENCE   = 7;
var TABLE_SPECTRA_FIELD_SEQ_HOMOLOG     = 8;
var TABLE_SPECTRA_FIELD_SEQ_DENOVO      = 9;
var TABLE_SPECTRA_FIELD_SEQ_USER        = 10;
var TABLE_SPECTRA_FIELD_MASS            = 11;
var TABLE_SPECTRA_FIELD_CHARGE          = 12;
var TABLE_SPECTRA_FIELD_B_PER           = 13;
var TABLE_SPECTRA_FIELD_Y_PER           = 14;
var TABLE_SPECTRA_FIELD_BY_INTENSITY    = 15;
var TABLE_SPECTRA_FIELD_INPUT_FILENAME  = 16;

var TAG_TABLE_SPECTRA_FIELD_ID              = makeTag(TABLE_SPECTRA_FIELD_ID);
var TAG_TABLE_SPECTRA_FIELD_IDX             = makeTag(TABLE_SPECTRA_FIELD_IDX);
var TAG_TABLE_SPECTRA_FIELD_SCAN            = makeTag(TABLE_SPECTRA_FIELD_SCAN);
var TAG_TABLE_SPECTRA_FIELD_CLUSTER         = makeTag(TABLE_SPECTRA_FIELD_CLUSTER);
var TAG_TABLE_SPECTRA_FIELD_PROTEIN_NAME    = makeTag(TABLE_SPECTRA_FIELD_PROTEIN_NAME);
var TAG_TABLE_SPECTRA_FIELD_PKLBIN_IDX      = makeTag(TABLE_SPECTRA_FIELD_PKLBIN_IDX);
var TAG_TABLE_SPECTRA_FIELD_PKLBIN_FILENAME = makeTag(TABLE_SPECTRA_FIELD_PKLBIN_FILENAME);
var TAG_TABLE_SPECTRA_FIELD_SEQ_REFERENCE   = makeTag(TABLE_SPECTRA_FIELD_SEQ_REFERENCE);
var TAG_TABLE_SPECTRA_FIELD_SEQ_HOMOLOG     = makeTag(TABLE_SPECTRA_FIELD_SEQ_HOMOLOG);
var TAG_TABLE_SPECTRA_FIELD_SEQ_DENOVO      = makeTag(TABLE_SPECTRA_FIELD_SEQ_DENOVO);
var TAG_TABLE_SPECTRA_FIELD_SEQ_USER        = makeTag(TABLE_SPECTRA_FIELD_SEQ_USER);
var TAG_TABLE_SPECTRA_FIELD_MASS            = makeTag(TABLE_SPECTRA_FIELD_MASS);
var TAG_TABLE_SPECTRA_FIELD_CHARGE          = makeTag(TABLE_SPECTRA_FIELD_CHARGE);
var TAG_TABLE_SPECTRA_FIELD_B_PER           = makeTag(TABLE_SPECTRA_FIELD_B_PER);
var TAG_TABLE_SPECTRA_FIELD_Y_PER           = makeTag(TABLE_SPECTRA_FIELD_Y_PER);
var TAG_TABLE_SPECTRA_FIELD_BY_INTENSITY    = makeTag(TABLE_SPECTRA_FIELD_BY_INTENSITY);
var TAG_TABLE_SPECTRA_FIELD_INPUT_FILENAME  = makeTag(TABLE_SPECTRA_FIELD_INPUT_FILENAME);

var TAG_SPECTRA_PEPTIDE_ALL  = makeTag(TABLE_SPECTRA_FIELD_SEQ_REFERENCE + '|' + TABLE_SPECTRA_FIELD_SEQ_HOMOLOG + '|' + TABLE_SPECTRA_FIELD_SEQ_DENOVO);

// Proteins coverage table
var TABLE_COVERAGE_FIELD_ID              = 0;
var TABLE_COVERAGE_FIELD_NAME            = 1;
var TABLE_COVERAGE_FIELD_SEQ_REFERENCE   = 2;
var TABLE_COVERAGE_FIELD_PROT_SEQUENCE   = 3;
var TABLE_COVERAGE_CSPS_DATA             = 4;
var TABLE_COVERAGE_SPS_DATA              = 5;

////////////////////////////////////////////////////////////////////////////////
// Object types in report elements
var REPORT_CELL_TYPE_IOD              = "iod";
var REPORT_CELL_TYPE_STRING           = "str";
var REPORT_CELL_TYPE_STRING_MULTIPLE  = "strM";
var REPORT_CELL_TYPE_BOX              = "box";

////////////////////////////////////////////////////////////////////////////////
// Image id tags prefixes:
//
// image shown (src)       ->  im_ , imc_
// image on demand (href)  ->  io_ , ioc_
//
var IMAGE_ICON_ID_PREFIX        = "im_";
var IMAGE_ICON_CTRL_ID_PREFIX   = "imc_";
var IMAGE_LARGE_ID_PREFIX       = "io_";
var IMAGE_LARGE_CTRL_ID_PREFIX  = "ioc_";

////////////////////////////////////////////////////////////////////////////////
// System internal tags

var INTERNAL_ROW      = "row";
var INTERNAL_COL      = "col";
var INTERNAL_PROJDIR  = "projectdir";

var TAG_INTERNAL_ROW      = makeTag(INTERNAL_ROW);
var TAG_INTERNAL_COL      = makeTag(INTERNAL_COL);
var TAG_INTERNAL_PROJDIR  = makeTag(INTERNAL_PROJDIR);

////////////////////////////////////////////////////////////////////////////////
// renderers & server objects
var SPS_REPORTS       = "spsReports.cgi";
var CONTIG_RENDERER   = "contplot.cgi";
var SPECTRUM_RENDERER = "specplot.cgi";

// spsReports parameters names
var REPORT_SERVER_PAR_GET     = "--get";
var REPORT_SERVER_PAR_UPDATE  = "--update";
var REPORT_SERVER_PAR_LAUNCH  = "--launch";
var REPORT_SERVER_PAR_STATUS  = "--status";

var REPORT_SERVER_PAR_REQUEST_ID    = "--request-id";
var REPORT_SERVER_PAR_TABLE         = "--table";
var REPORT_SERVER_PAR_FILTER_FIELD  = "--filter-field";
var REPORT_SERVER_PAR_FILTER_DATA   = "--filter-data";
var REPORT_SERVER_PAR_UPDATE_FIELD  = "--update-field";
var REPORT_SERVER_PAR_UPDATE_DATA   = "--update-data";
var REPORT_SERVER_PAR_CLEAR_DATA    = "--clear-data";

var REPORT_SERVER_PAR_PROJECT_DIR   = "--project-dir";
var REPORT_SERVER_PAR_FILENAME      = "--filename";
var REPORT_SERVER_PAR_DESCRIPTION   = "--description";
var REPORT_SERVER_PAR_SEQUENCE      = "--sequence";
var REPORT_SERVER_PAR_ID            = "--ID";

// Renderers' parameters names (contplot)
var CONTPLOT_PAR_STAR           = "--star";
var CONTPLOT_PAR_ABRUIJN        = "--abruijn";
var CONTPLOT_PAR_SEQS           = "--seqs";
var CONTPLOT_PAR_TITLE          = "--title";
var CONTPLOT_PAR_NO_TITLE       = "--notitle";
var CONTPLOT_PAR_ZOOM           = "--zoom";
var CONTPLOT_PAR_CONTIG         = "--contig";
var CONTPLOT_PAR_MASS_REF       = "--mass-reference";
var CONTPLOT_PAR_MASS_HOM       = "--mass-homolog";
var CONTPLOT_PAR_OFF_REF        = "--offset-reference";
var CONTPLOT_PAR_OFF_HOM        = "--offset-homolog";
var CONTPLOT_PAR_REVERSE        = "--reverse";
var CONTPLOT_PAR_TARGET         = "--target";
var CONTPLOT_PAR_ENCODING       = "--encoding";
var CONTPLOT_PAR_SEQ_REFERENCE  = "--reference";
var CONTPLOT_PAR_SEQ_HOMOLOG    = "--homolog";
var CONTPLOT_PAR_SEQ_DENOVO     = "--denovo";
var CONTPLOT_PAR_SEQ_USER       = "--user";

// Renderers' parameters names (specplot)
var SPECPLOT_PAR_PKLBIN         = "--pklbin";
var SPECPLOT_PAR_SPECTURM       = "--spectrum";
var SPECPLOT_PAR_PEPTIDE        = "--peptide";
var SPECPLOT_PAR_TARGET         = "--target";
var SPECPLOT_PAR_ENCODING       = "--encoding";
var SPECPLOT_PAR_ZOOM           = "--zoom";
var SPECPLOT_PAR_TITLE          = "--title";
var SPECPLOT_PAR_NOTITLE        = "--notitle";



// renderers' parameters values (contplot) " +  + "
var CONTPLOT_VAL_STAR           = TAG_TABLE_CONTIGS_FIELD_FILE_STARS;
var CONTPLOT_VAL_ABRUIJN        = TAG_TABLE_CONTIGS_FIELD_FILE_ABRUIJN;
var CONTPLOT_VAL_SEQS           = TAG_TABLE_CONTIGS_FIELD_FILE_SEQS;
var CONTPLOT_VAL_TITLE          = "\"Contig " + TAG_TABLE_CONTIGS_FIELD_ID + "\"";
var CONTPLOT_VAL_CONTIG         = TAG_TABLE_CONTIGS_FIELD_ID;
var CONTPLOT_VAL_MASS_REF       = "\"" + TAG_TABLE_CONTIGS_FIELD_MASS_REF + "\"";
var CONTPLOT_VAL_MASS_HOM       = "\"" + TAG_TABLE_CONTIGS_FIELD_MASS_HOM + "\"";
var CONTPLOT_VAL_OFF_REF        = TAG_TABLE_CONTIGS_FIELD_OFF_REF;
var CONTPLOT_VAL_OFF_HOM        = TAG_TABLE_CONTIGS_FIELD_OFF_HOM;
var CONTPLOT_VAL_REVERSE        = "";
var CONTPLOT_VAL_TARGET         = "cout";
var CONTPLOT_VAL_ENCODING       = "uu64";
var CONTPLOT_VAL_SEQ_REFERENCE  = "\"" + TAG_TABLE_CONTIGS_FIELD_SEQ_REF + "\"";
var CONTPLOT_VAL_SEQ_HOMOLOG    = "\"" + TAG_TABLE_CONTIGS_FIELD_SEQ_HOM + "\"";
var CONTPLOT_VAL_SEQ_DENOVO     = "\"" + TAG_TABLE_CONTIGS_FIELD_SEQ_NOVO + "\"";
var CONTPLOT_VAL_SEQ_USER       = "\"" + TAG_TABLE_CONTIGS_FIELD_SEQ_USER + "\"";
var CONTPLOT_VAL_NO_TITLE       = "";
var CONTPLOT_VAL_ZOOM           = "0.4";


// renderers' parameters values (specplot)
var SPECPLOT_VAL_TARGET         = "cout";
var SPECPLOT_VAL_ENCODING       = "uu64";

// renderers' parameters values (specplot for cluster)
var SPECPLOT_VAL_CLUSTER_NOTITLE        = "";
var SPECPLOT_VAL_CLUSTER_TITLE          = "\"Consensus Spectrum " + TAG_TABLE_CLUSTER_FIELD_ID + "(contig " + TAG_TABLE_CLUSTER_FIELD_CONTIG + ")\"";
var SPECPLOT_VAL_CLUSTER_PKLBIN         = TAG_INTERNAL_PROJDIR + "/spectra/specs_ms.pklbin";
var SPECPLOT_VAL_CLUSTER_SPECTURM       = TAG_TABLE_CLUSTER_FIELD_ID;
var SPECPLOT_VAL_CLUSTER_ZOOM           = "0.35";
var SPECPLOT_VAL_CLUSTER_PEPTIDE_ALL    = "\"" + TAG_CLUSTER_PEPTIDE_ALL + "\"";  // <3|4|5>
var SPECPLOT_VAL_CLUSTER_PEPTIDE_REF    = "\"" + TAG_TABLE_CLUSTER_FIELD_REFERENCE + "\"";
var SPECPLOT_VAL_CLUSTER_PEPTIDE_HOM    = "\"" + TAG_TABLE_CLUSTER_FIELD_HOMOLOG + "\"";
var SPECPLOT_VAL_CLUSTER_PEPTIDE_NOVO   = "\"" + TAG_TABLE_CLUSTER_FIELD_DENOVO + "\"";
var SPECPLOT_VAL_CLUSTER_PEPTIDE_USER   = "\"" + TAG_TABLE_CLUSTER_FIELD_USER + "\"";


// renderers' parameters values (specplot for input spectra)
var SPECPLOT_VAL_SPECTRA_NOTITLE        = "";
var SPECPLOT_VAL_SPECTRA_TITLE          = "\"Spectrum Scan " + TAG_TABLE_SPECTRA_FIELD_SCAN + "\"";
var SPECPLOT_VAL_SPECTRA_PKLBIN         = TAG_INTERNAL_PROJDIR + "/" + TAG_TABLE_SPECTRA_FIELD_PKLBIN_FILENAME;
var SPECPLOT_VAL_SPECTRA_SPECTURM       = TAG_TABLE_SPECTRA_FIELD_IDX;
var SPECPLOT_VAL_SPECTRA_ZOOM           = "0.35";
var SPECPLOT_VAL_SPECTRA_PEPTIDE_ALL    = "\"" + TAG_SPECTRA_PEPTIDE_ALL + "\"";
var SPECPLOT_VAL_SPECTRA_PEPTIDE_REF    = "\"" + TAG_TABLE_SPECTRA_FIELD_SEQ_REFERENCE + "\"";
var SPECPLOT_VAL_SPECTRA_PEPTIDE_HOM    = "\"" + TAG_TABLE_SPECTRA_FIELD_SEQ_HOMOLOG + "\"";
var SPECPLOT_VAL_SPECTRA_PEPTIDE_NOVO   = "\"" + TAG_TABLE_SPECTRA_FIELD_SEQ_DENOVO + "\"";
var SPECPLOT_VAL_SPECTRA_PEPTIDE_USER   = "\"" + TAG_TABLE_SPECTRA_FIELD_SEQ_USER + "\"";


// renderers' parameters conditions (contplot)
var CONTPLOT_COND_STAR           = "";
var CONTPLOT_COND_ABRUIJN        = "";
var CONTPLOT_COND_SEQS           = "";
var CONTPLOT_COND_TITLE          = "";
var CONTPLOT_COND_CONTIG         = TAG_TABLE_CONTIGS_FIELD_ID;
var CONTPLOT_COND_MASS_REF       = TAG_TABLE_CONTIGS_FIELD_MASS_REF;
var CONTPLOT_COND_MASS_HOM       = TAG_TABLE_CONTIGS_FIELD_MASS_HOM;
var CONTPLOT_COND_OFF_REF        = TAG_TABLE_CONTIGS_FIELD_OFF_REF;
var CONTPLOT_COND_OFF_HOM        = TAG_TABLE_CONTIGS_FIELD_OFF_HOM;
var CONTPLOT_COND_REVERSE        = TAG_TABLE_CONTIGS_FIELD_REVERSE;
var CONTPLOT_COND_TARGET         = "";
var CONTPLOT_COND_ENCODING       = "";
var CONTPLOT_COND_SEQ_REFERENCE  = TAG_TABLE_CONTIGS_FIELD_SEQ_REF;
var CONTPLOT_COND_SEQ_HOMOLOG    = TAG_TABLE_CONTIGS_FIELD_SEQ_HOM;
var CONTPLOT_COND_SEQ_DENOVO     = TAG_TABLE_CONTIGS_FIELD_SEQ_NOVO;
var CONTPLOT_COND_SEQ_USER       = TAG_TABLE_CONTIGS_FIELD_SEQ_USER;
var CONTPLOT_COND_NO_TITLE       = "";
var CONTPLOT_COND_ZOOM           = "";


// renderers' parameters conditions (specplot)
var SPECPLOT_COND_TARGET         = "";
var SPECPLOT_COND_ENCODING       = "";

// renderers' parameters condition (specplot for cluster)
var SPECPLOT_COND_CLUSTER_NOTITLE        = "";
var SPECPLOT_COND_CLUSTER_TITLE          = "";
var SPECPLOT_COND_CLUSTER_PKLBIN         = "";
var SPECPLOT_COND_CLUSTER_SPECTURM       = "";
var SPECPLOT_COND_CLUSTER_ZOOM           = "";
var SPECPLOT_COND_CLUSTER_PEPTIDE_ALL    = "";
var SPECPLOT_COND_CLUSTER_PEPTIDE_REF    = "";
var SPECPLOT_COND_CLUSTER_PEPTIDE_HOM    = "";
var SPECPLOT_COND_CLUSTER_PEPTIDE_NOVO   = "";
var SPECPLOT_COND_CLUSTER_PEPTIDE_USER   = "";

// renderers' parameters condition (specplot for input spectra)
var SPECPLOT_COND_SPECTRA_NOTITLE        = "";
var SPECPLOT_COND_SPECTRA_TITLE          = "";
var SPECPLOT_COND_SPECTRA_PKLBIN         = "";
var SPECPLOT_COND_SPECTRA_SPECTURM       = "";
var SPECPLOT_COND_SPECTRA_ZOOM           = "";
var SPECPLOT_COND_SPECTRA_PEPTIDE_ALL    = "";
var SPECPLOT_COND_SPECTRA_PEPTIDE_REF    = "";
var SPECPLOT_COND_SPECTRA_PEPTIDE_HOM    = "";
var SPECPLOT_COND_SPECTRA_PEPTIDE_NOVO   = "";
var SPECPLOT_COND_SPECTRA_PEPTIDE_USER   = "";


////////////////////////////////////////////////////////////////////////////////
// page elements
var DIV_PAGE_MAIN  = "mainDiv2";
var DIV_PAGE_DATA  = "mainDiv";

// page type IDs
var PAGE_INITIAL          = 10;
var PAGE_PROTEINS         = 20;
var PAGE_PROTEIN          = 30;
var PAGE_PROTEIN_COVERAGE = 40;
var PAGE_CONTIGS          = 50;
var PAGE_CONTIG           = 60;
var PAGE_CLUSTER          = 70;
var PAGE_SPECTRA          = 80;


////////////////////////////////////////////////////////////////////////////////
// Message to be displayed on connection error.
var CONN_ERR = "XMLHTTP not available. Could not connect to server.";

////////////////////////////////////////////////////////////////////////////////
// Variables
////////////////////////////////////////////////////////////////////////////////

// keep track of what king of page is being displayed (page type)
var currentPageType;
// keep track of what page is being displayed (data)
var currentPageData;

// the current table
var tab;

// Global image ID. Used to generate unique IDs for images for image placement upon arrival from server
var globalImage = 0;

// Current protein coverage. Used to gather sequence in protein coverage pages to send to server for reprocessing.
var globalProteinLength = -1;

// status variables
var globalStatus = "";

// Image processing queue
var queue = new Array();

// Image cache
var iCache = new ImageCache();

// Global request ID, used to fool "smart" browsers that cache AJAX requests
var requestID = 0;

// Protein coverage buffer -- used to keep protien edition state while navigatin away
var proteinCoverageState = new Array();

////////////////////////////////////////////////////////////////////////////////
// Images
////////////////////////////////////////////////////////////////////////////////

var iconHome = "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";

var iconProteinList = "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";

var iconProtein = "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";

var iconProteinCoverage = "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";

var iconContig = "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";

var iconContigList = "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";

var loadingImage = "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";

var loadingImage2 = "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";

var workingImage = "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";

////////////////////////////////////////////////////////////////////////////////
// Initializers
////////////////////////////////////////////////////////////////////////////////

// Initialize page variables from HTML initial page
function init()
{
  var obj = document.getElementById(PAGE_FIELD_CELLS_PER_LINE);
  if(obj)
    CELLS_PER_LINE = parseInt(obj.value);

  var obj = document.getElementById(PAGE_FIELD_SERVER);
  if(obj)
    SERVER = obj.value;

  var obj = document.getElementById(PAGE_FIELD_PROJECT_DIR);
  if(obj)
    PROJECT_DIR = obj.value;

  var obj = document.getElementById(PAGE_FIELD_FASTA_FILENAME);
  if(obj)
    GLOBAL_FASTA_FILENAME = obj.value;

  getStatus();
}
////////////////////////////////////////////////////////////////////////////////
function showMain(ok)
{
	var main  = document.getElementById(DIV_PAGE_MAIN);
	var pages = document.getElementById(DIV_PAGE_DATA);

	if(ok) {
    main.style.display  = "block";
    pages.style.display = "none";
  }	else {
    main.style.display  = "none";
    pages.style.display = "block";
	}
}
////////////////////////////////////////////////////////////////////////////////
// AJAX main call
////////////////////////////////////////////////////////////////////////////////
function AjaxGlobal(sURL, sMethod, sVars, fnDone, sync)
{
  var myConn = new XHConn();
  if (!myConn) alert(CONN_ERR);
  myConn.connect(sURL, sMethod, sVars + '&' + REPORT_SERVER_PAR_REQUEST_ID + '=' + requestID++, fnDone, sync);
}
////////////////////////////////////////////////////////////////////////////////
function getStatus()
{
  var field = document.getElementById("status");

  fnDone7 = function (oXML) {
    if(field) {
      globalStatus = oXML.responseText;
      field.innerHTML = globalStatus;
    }
    setTimeout(function() {getStatus();}, 10000);
  };

  if(field) {
    var thePage = SERVER + SPS_REPORTS;
    var params = REPORT_SERVER_PAR_STATUS;
    params += '&' + REPORT_SERVER_PAR_PROJECT_DIR + '&' + PROJECT_DIR;
    params += '&' + REPORT_SERVER_PAR_FILENAME    + '&' + STATUS_FILENAME;
    AjaxGlobal(thePage, "GET", params, fnDone7, true);
  }
}
////////////////////////////////////////////////////////////////////////////////
// Image Cache
////////////////////////////////////////////////////////////////////////////////
function cacheElement()
{
  // renderer to be used
  this.renderer   = "";
  // parameters for renderer
  this.params     = "";
  // where to store the non changed parameters for the image (only the user sequence should change) <input> ID
  this.tag2params = "";
  // the unique ID used to identify the image
  this.id         = 0;
  // the user sequence, to check for change
  this.sequence   = "";
  // the image contents
  this.image      = "";
}
////////////////////////////////////////////////////////////////////////////////
// Translate string to number
function string2Number(str)
{
  var v = 0;
  for(var i = 0 ; i < str.length ; i++)
    v += str[i].charCodeAt() * Math.pow(10,i);
  return v;
}
////////////////////////////////////////////////////////////////////////////////
// Image cache object. Stores images do avoid duplicate requests to server
function ImageCache()
{
  // the actual cache
  this.cache = new BST();
  //this.cache = new Array();
}

ImageCache.prototype.add = function(e)
{
  var id = string2Number(e.id);
  this.cache.insert(id, e);
  //this.cache.push(e);
}

ImageCache.prototype.update = function(i, element)
{
  var obj = this.cache.search(i);
  obj.value = element;

  //this.cache[i] = element;
}

ImageCache.prototype.get = function(i)
{
  var aux = this.cache.search(i);
  return aux.value;
  
  //return this.cache[i];
}

ImageCache.prototype.find = function(key)
{
  var id = string2Number(key);
  node = this.cache.search(id);
  if(node != null)
    return node.key;

  //for(var i = 0 ; i < this.cache.length ; i++)
  //  if(this.cache[i].id == key)
  //    return i;
  
  return -1;
}
////////////////////////////////////////////////////////////////////////////////
// Image queue
////////////////////////////////////////////////////////////////////////////////
function queueElement()
{
  // renderer to be used
  this.renderer   = "";
  // parameters for renderer
  this.params     = "";
  // where to put the image <image> ID
  this.tag        = "";
  // where to store the params for later usage (image update) <input> ID
  this.tag2       = "";
  // where to store the non changed parameters for the image (only the user sequence should change) <input> ID
  this.tag2params = "";
  // user sequence
  this.sequence   = "";
  // href or src. Used to now which object tag should contain the image
  this.target     = "";
}
////////////////////////////////////////////////////////////////////////////////
// Send assyncronous image requests to server and put them in place
function processQueue()
{
  for(var i = 0 ; i < queue.length ; i++) {

    // get image
    loadImage(queue[i].renderer, queue[i].params, queue[i].tag, queue[i].target, false, true, queue[i].sequence);

    // store parameters for image update
    var obj = document.getElementById(queue[i].tag2);
    if(obj) {
      obj.value = queue[i].tag2params;
    }
  }

  // clear queue
  queue = [];
}
////////////////////////////////////////////////////////////////////////////////
// Reference Utilities
////////////////////////////////////////////////////////////////////////////////
// Create a variable reference to allow primitive data types to be passed-by-reference
function createReference(value)
{
   var newVar = new Array();
   newVar[0] = value;
   return newVar;
}
////////////////////////////////////////////////////////////////////////////////
// set reference variable's value
function setReference(val, value)
{
   val[0] = value;
}
////////////////////////////////////////////////////////////////////////////////
// add to variable
function addToReference(val, value)
{
   val[0] += value;
}
////////////////////////////////////////////////////////////////////////////////
// Auxiliary functions
////////////////////////////////////////////////////////////////////////////////
// Sort comparison function.
function orderSortingFunction(v1, v2)
{
  return (v1[1] == v2[1] ? v1[2] - v2[2] : v1[1] - v2[1]);
}
////////////////////////////////////////////////////////////////////////////////
// generate a click event for the object
function triggerClick(id)
{
  var clicky = document.createEvent('HTMLEvents');clicky.initEvent('click', true, true);
  var aux = document.getElementById(id);
  if(aux)
    aux.dispatchEvent(clicky);
}
////////////////////////////////////////////////////////////////////////////////
// Load entry point function and image load functions
////////////////////////////////////////////////////////////////////////////////
function loadImage(renderer, params, theTag, target, force, sync, sequence)
{
  // image chache index
  var index = -1;

  // cache search section

  if(!force && USE_IMAGE_CACHE) {
    // check for image in cache
    index = iCache.find(theTag);
    // ignore cache if not found
    if(index != -1) {
      // get element
      var elem = iCache.get(index);
      //check for sequence changes
      if(elem.sequence == sequence) {
        // set the image
        var obj = document.getElementById(theTag);
        if(obj) {
          if(target == "href")
            obj.href = "data:image/png;base64," + elem.image;
          else
            obj.src = "data:image/png;base64," + elem.image;
        }
        // and exit
        return;
      }
    }
  }

  // we need to get the image

  // create a local variable so that it can be used by callback function
  var tag = theTag;

  // callback function for assyncronous call to server
  fnDone = function (oXML) {

    var obj2 = document.getElementById(tag);
    if(obj2) {
      if(target == "href")
        obj2.href = "data:image/png;base64," + oXML.responseText;
      else
        obj2.src = "data:image/png;base64," + oXML.responseText;
    }

    // add element to cache
    var aux = new cacheElement();
    aux.renderer   = renderer;
    aux.params     = params;
    aux.tag2params = "";
    aux.id         = tag;
    aux.sequence   = sequence;
    aux.image      = oXML.responseText;

    if(USE_IMAGE_CACHE)
      if(index == -1)
        iCache.add(aux);
      else
        iCache.update(index, aux);
  };


  // what to get
  var thePage = SERVER + renderer;
  // connect to server
  AjaxGlobal(thePage, "GET", params, fnDone, sync);
}
////////////////////////////////////////////////////////////////////////////////
function validateStatus()
{
  if(globalStatus != STATUS_OK)
    return false;
  showMain(false);
  return true;
}
////////////////////////////////////////////////////////////////////////////////
function afterPageload(page, data)
{
  switch(page) {

  case PAGE_INITIAL:
    break;

  case PAGE_PROTEIN:
    break;

  case PAGE_PROTEIN_COVERAGE:
    restoreCoverage(data);
    break;

  case PAGE_CONTIGS:
    break;

  case PAGE_CONTIG:
    break;

  case PAGE_CLUSTER:
    break;

  case PAGE_SPECTRA:
    break;

  }
  // set current page type
  currentPageType = page;

  // set current page data;
  currentPageData = data;
}
////////////////////////////////////////////////////////////////////////////////
function onPageExit(page, data)
{
  switch(page) {

  case PAGE_INITIAL:
    break;

  case PAGE_PROTEINS:
    break;

  case PAGE_PROTEIN:
    break;

  case PAGE_PROTEIN_COVERAGE:
    saveCoverage(data);
    break;

  case PAGE_CONTIGS:
    break;

  case PAGE_CONTIG:
    break;

  case PAGE_CLUSTER:
    break;

  case PAGE_SPECTRA:
    break;

  }
}
////////////////////////////////////////////////////////////////////////////////
function loadProteinsPage()
{
  if(!validateStatus())
    return;

  // exit page procedure
  onPageExit(currentPageType, currentPageData);

  // set loading image
  setLoadingImage(DIV_PAGE_DATA);

  // build report object
  rep = new ReportProtein();
  rep.proteinsPage();

  // retrieve tables data
  rep.getData();

  setWorkingImage(DIV_PAGE_DATA);
  setTimeout(function() {

  // define renderer
  var r = new renderer();
  // build the report
  r.buildReport(DIV_PAGE_DATA, rep);
  // entry page procedure
  afterPageload(PAGE_PROTEINS);

  }, 10);
}
////////////////////////////////////////////////////////////////////////////////
function loadProteinPage(protein, start)
{
  if(!validateStatus())
    return;

  // exit page procedure
  onPageExit(currentPageType, currentPageData);

  // set loading image
  setLoadingImage(DIV_PAGE_DATA);

  // build report object
  rep = new ReportProtein();
  rep.proteinPage(protein);

  // retrieve tables data
  rep.getData();

  setWorkingImage(DIV_PAGE_DATA);
  setTimeout(function() {

  // define renderer
  var r = new renderer();
  // build entry point for pagination
  var fn = "loadProteinPage(" + protein + ",";
  // build the report
  r.buildReport(DIV_PAGE_DATA, rep, start, fn, 1);
  // entry page procedure
  afterPageload(PAGE_PROTEIN, protein);

  }, 10);
}
////////////////////////////////////////////////////////////////////////////////
function loadContigsPage(start)
{
  if(!validateStatus())
    return;

  // exit page procedure
  onPageExit(currentPageType, currentPageData);

  // set loading image
  setLoadingImage(DIV_PAGE_DATA);

  // build report object
  rep = new ReportContig();
  rep.contigsPage();

  // retrieve tables data
  rep.getData();

  setWorkingImage(DIV_PAGE_DATA);
  setTimeout(function() {

  // define renderer
  var r = new renderer();
  // build entry point for pagination
  var fn = "loadContigsPage(";
  // build the report
  r.buildReport(DIV_PAGE_DATA, rep, start, fn);
  // entry page procedure
  afterPageload(PAGE_CONTIGS);

  }, 10);
}
////////////////////////////////////////////////////////////////////////////////
function loadContigPage(contig, start)
{
  if(!validateStatus())
    return;

  // exit page procedure
  onPageExit(currentPageType, currentPageData);

  // set loading image
  setLoadingImage(DIV_PAGE_DATA);

  // build report object
  rep = new ReportContig();
  rep.contigPage(contig);

  // retrieve tables data
  rep.getData();

  // get protein ID
  protein = rep.getField(0, 0, 1);
  var navType = 2;
  if(protein[0] == -1) navType = 3;

  setWorkingImage(DIV_PAGE_DATA);
  setTimeout(function() {

  // define renderer
  var r = new renderer();
  // build entry point for pagination
  var fn = "loadContigPage(" + contig + ",";
  // build the report
  r.buildReport(DIV_PAGE_DATA, rep, start, fn, navType);
  // entry page procedure
  afterPageload(PAGE_CONTIG, contig);

  }, 10);
}
////////////////////////////////////////////////////////////////////////////////
function loadClusterPage(cluster, start)
{
  if(!validateStatus())
    return;

  // exit page procedure
  onPageExit(currentPageType, currentPageData);

  // set loading image
  setLoadingImage(DIV_PAGE_DATA);

  // build report object
  rep = new ReportCluster();
  rep.clusterPage(cluster);

  // retrieve tables data
  rep.getData();

  // get protein ID
  protein = rep.getField(0, 0, 2);
  var navType = 2;
  if(protein[0] == 0) navType = 4;

  setWorkingImage(DIV_PAGE_DATA);
  setTimeout(function() {

  // define renderer
  var r = new renderer();
  // build entry point for pagination
  var fn = "loadClusterPage(" + cluster + ",";
  // build the report
  r.buildReport(DIV_PAGE_DATA, rep, start, fn, navType);
  // entry page procedure
  afterPageload(PAGE_CLUSTER, cluster);

  }, 10);
}
////////////////////////////////////////////////////////////////////////////////
function loadInputSpectraPage(fileIndex, start)
{
  if(!validateStatus())
    return;

  // exit page procedure
  onPageExit(currentPageType, currentPageData);

  // set loading image
  setLoadingImage(DIV_PAGE_DATA);

  // build report object
  rep = new ReportInputSpectra();
  rep.inputSpectraPage(fileIndex);

  // retrieve tables data
  rep.getData();

  setWorkingImage(DIV_PAGE_DATA);
  setTimeout(function() {

  // define renderer
  var r = new renderer();
  // build entry point for pagination
  var fn = "loadInputSpectraPage(" + fileIndex + ",";
  // build the report
  r.buildReport(DIV_PAGE_DATA, rep, start, fn);
  // entry page procedure
  afterPageload(PAGE_SPECTRA, fileIndex);

  }, 10);
}
////////////////////////////////////////////////////////////////////////////////
function loadProteinDetails(protein)
{
  if(!validateStatus())
    return;

  // exit page procedure
  onPageExit(currentPageType, currentPageData);

  // set loading image
  setLoadingImage(DIV_PAGE_DATA);

  // build report object
  rep = new ReportProteinCoverage();
  rep.proteinCoveragePage(protein);

  // retrieve tables data
  rep.getData();

  setWorkingImage(DIV_PAGE_DATA);
  setTimeout(function() {

  // define renderer
  var r = new renderer();
  // build the report
  r.buildReport(DIV_PAGE_DATA, rep, null, "", 2);
  // entry page procedure
  afterPageload(PAGE_PROTEIN_COVERAGE, protein);

  }, 10);
}
////////////////////////////////////////////////////////////////////////////////
function loadProteinDetailsCSV(protein)
{
  if(!validateStatus())
    return;

  // build report object
  rep = new ReportProteinCoverage();
  rep.proteinCoverageCSVPage(protein);

  // retrieve tables data
  rep.getData();

  // define renderer
  var r = new renderer();
  // build the report. NEW means it will open a new window/tab
  r.buildReport('NEW', rep, null, "", 2);
}
////////////////////////////////////////////////////////////////////////////////
function setLoadingImage(div)
{
  setImage(div, loadingImage2);
}

function setWorkingImage(div)
{
  setImage(div, workingImage);
}


function setImage(div, image)
{
  // set loading image
  var aux = "<table align='center'><tr><td align='center'><img src='data:image/gif;base64," + image + "' /></td></tr></table>";
  // set data on page
  var target = document.getElementById(div);
  if(target)
    target.innerHTML = aux;
}
////////////////////////////////////////////////////////////////////////////////
function pausecomp(ms)
{
  ms += new Date().getTime();
  while (new Date() < ms) {}
}
////////////////////////////////////////////////////////////////////////////////
function makeTag(content)
{
  return TAG_OPEN + content + TAG_CLOSE;
}
////////////////////////////////////////////////////////////////////////////////
function updateContig(fieldId, textId, img, ctrl, target, row, col)
{
  var inData  = document.getElementById(fieldId);
  var outData = document.getElementById(textId);

  // test for inexistent input element
  if(!inData)
    return;

  // get data
  var str = inData.value.toUpperCase();

  // test for target and set string
  if(outData)
    outData.innerHTML = str;

  // update the tables thru the server
  updateField(TABLE_CONTIG_ID, TABLE_CONTIGS_FIELD_ID, row, col, str);

  // get params for image generation
  var paramsObj = document.getElementById(ctrl);
  if(!paramsObj) return;
  var params = paramsObj.value;

  if(str)
    params += '&' + CONTPLOT_PAR_SEQ_USER + '="' + str + '"';
  // load the image
  loadImage(CONTIG_RENDERER, params, img, target, false, true, str)
}
////////////////////////////////////////////////////////////////////////////////
function updateCluster(fieldId, textId, img, ctrl, target, row, col, sync)
{
  var inData  = document.getElementById(fieldId);
  var outData = document.getElementById(textId);

  var str = inData.value.toUpperCase();

  if(outData)
    outData.innerHTML = str;

  // update the tables thru the server
  updateField(TABLE_CLUSTER_ID, TABLE_CLUSTER_FIELD_ID, row, col, str);

  // get params for image generation
  var paramsObj = document.getElementById(ctrl);
  if(!paramsObj) return;
  var params = paramsObj.value;

  if(str)
    params += '&' + SPECPLOT_PAR_PEPTIDE+ '="' + str + '"';
  // load the image
  loadImage(SPECTRUM_RENDERER, params, img, target, false, sync, str)
}
////////////////////////////////////////////////////////////////////////////////
function updateSpectra(fieldId, textId, img, ctrl, target, row, col, sync)
{
  var inData  = document.getElementById(fieldId);
  var outData = document.getElementById(textId);

  var str = inData.value.toUpperCase();

  if(outData)
    outData.innerHTML = str;

  // update the tables thru the server
  updateField(TABLE_SPECTRA_ID, TABLE_SPECTRA_FIELD_ID, row, col, str);

 // get params for image generation
  var paramsObj = document.getElementById(ctrl);
  if(!paramsObj) return;
  var params = paramsObj.value;

  if(str)
    params += '&' + SPECPLOT_PAR_PEPTIDE + '="' + str + '"';
  // load the image
  loadImage(SPECTRUM_RENDERER, params, img, target, false, sync, str)
}
////////////////////////////////////////////////////////////////////////////////
function fillCoverageRow(obj, i)
{
  if(!obj.checked) return;
  for(var j = i ; j < i + CELLS_PER_LINE && j < globalProteinLength ; j++) {
    var aux1 = PEP_ELEM_PREFIX + j;
    var aux2 = INP_ELEM_PREFIX + j;
    var obj1 = document.getElementById(aux1);
    var obj2 = document.getElementById(aux2);
    if(obj1 && obj2)
      obj2.value = obj1.innerHTML;
  }
}
////////////////////////////////////////////////////////////////////////////////
function saveCoverage(data)
{
  // rows holder. Holds the rows data
  var rows = new Array();

  // row index counter
  var k = 0;

  // go thru all AA cells
  for(var i = 0 ; i < globalProteinLength ; i += CELLS_PER_LINE) {

    // where to save the sequence. 1st element is 'checked', subsequence are cell contents for that row
    var row = new Array();

    // get row state checkbox
    var ckId = 'ck' + i;
    var ckObj = document.getElementById(ckId);

    // save row state
    row[0] = ckObj.checked;

    // go thru all AA cells for the row
    for(var j = i ; j < i + CELLS_PER_LINE && j < globalProteinLength ; j++) {

      // get the cell ID
      var dataId = INP_ELEM_PREFIX + j;
      var dataVal = document.getElementById(dataId);
      // get the cell content
      row[j - i + 1] = dataVal.value;
    }
    // store the row
    rows[k++] = row;
  }
  proteinCoverageState[data] = rows;
}
////////////////////////////////////////////////////////////////////////////////
function restoreCoverage(data)
{
  // check if data is present
  if(typeof(data) == 'undefined' || typeof(proteinCoverageState[data]) == 'undefined')
    return;

  // row index counter
  var k = 0;

  // go thru all AA cells
  for(var i = 0 ; i < globalProteinLength && k < proteinCoverageState[data].length ; i += CELLS_PER_LINE, k++) {

    // get row state checkbox
    var ckId = 'ck' + i;
    var ckObj = document.getElementById(ckId);

    // get row state
    ckObj.checked = proteinCoverageState[data][k][0];

    // go thru all AA cells for the row
    for(var j = i ; j < i + CELLS_PER_LINE && j < globalProteinLength; j++) {

      // get the cell ID
      var dataId = INP_ELEM_PREFIX + j;
      var dataVal = document.getElementById(dataId);

      // put the cell content
      if(j - i + 1 < proteinCoverageState[data][k].length)
        dataVal.value = proteinCoverageState[data][k][j - i + 1];
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
function submitCoverage()
{
  var sequence = "";

  for(var i = 0 ; i < globalProteinLength ; i += CELLS_PER_LINE) {

    var ckId = 'ck' + i;
    var ckObj = document.getElementById(ckId);

    for(var j = i ; j < i + CELLS_PER_LINE && j < globalProteinLength ; j++) {

      var dataId = INP_ELEM_PREFIX + j;
      var dataVal = document.getElementById(dataId);
      var pepId = PEP_ELEM_PREFIX + j;
      var pepVal = document.getElementById(pepId);

      //if(dataVal.value.length && ckObj.checked)
      //  sequence += dataVal.value.toUpperCase();
      //else
      //  sequence += pepVal.innerHTML;

      if(ckObj.checked) {
        if(dataVal.value.length)
          sequence += dataVal.value.toUpperCase();
      } else
        sequence += pepVal.innerHTML;
    }
  }

  var ProtID   = document.getElementById('ProtID').value;
  var ProtDesc = document.getElementById('ProtDESC').value;

  if(ProtID.length == 0 || ProtDesc.length == 0 || sequence.length ==0) {
    var msg = "Attention!\n\nThe following fields are empty:\n\n";
    msg += (ProtID.length == 0 ? "Protein ID\n" : "");
    msg += (ProtDesc.length == 0 ? "Description\n" : "");
    msg += (sequence.length == 0 ? "Protein sequence\n" : "");
    alert(msg);
    return;
  }

  var params = REPORT_SERVER_PAR_LAUNCH + '=relauncher.sh';
  params += '&' + REPORT_SERVER_PAR_PROJECT_DIR + '=' + PROJECT_DIR;
  params += '&' + REPORT_SERVER_PAR_FILENAME    + '&' + GLOBAL_FASTA_FILENAME;
  params += '&' + REPORT_SERVER_PAR_ID          + '=' + ProtID;
  params += '&' + REPORT_SERVER_PAR_DESCRIPTION + '=' + ProtDesc;
  params += '&' + REPORT_SERVER_PAR_SEQUENCE    + '=' + sequence;

  // callback
  fnDone8 = function (oXML) {
    // goto main page
    document.location = INITIAL_PAGE;
  };
  // AJAX call
  var thePage = SERVER + SPS_REPORTS;
  AjaxGlobal(thePage, "POST", params, fnDone8, false);
  // reset global status
  globalStatus = "";
}
////////////////////////////////////////////////////////////////////////////////
function updateField(table, filterField, filterData, updateField, updateData)
{
  var params = REPORT_SERVER_PAR_UPDATE;
  params += '&' + REPORT_SERVER_PAR_TABLE         + '=' + table;
  params += '&' + REPORT_SERVER_PAR_PROJECT_DIR   + '=' + PROJECT_DIR;

  params += '&' + REPORT_SERVER_PAR_FILTER_FIELD  + '=' + filterField;
  params += '&' + REPORT_SERVER_PAR_FILTER_DATA   + '=' + filterData;
  params += '&' + REPORT_SERVER_PAR_UPDATE_FIELD  + '=' + updateField;
  if(updateData.length)
    params += '&' + REPORT_SERVER_PAR_UPDATE_DATA + '=' + updateData;
  else
    params += '&' + REPORT_SERVER_PAR_CLEAR_DATA;

  // callback
  fnDone = function (oXML) {
  };
  // SJAX call
  var thePage = SERVER + SPS_REPORTS;
  AjaxGlobal(thePage, "GET", params, fnDone, true);
}
////////////////////////////////////////////////////////////////////////////////
// Column types
////////////////////////////////////////////////////////////////////////////////
function ReportParamsOption(p, o, v, s)
{
  this.param      = p;
  this.option     = o;
  this.validator  = v;
  this.store      = typeof(s) != 'undefined' ? s : 1;
}
////////////////////////////////////////////////////////////////////////////////
function ReportColumnTypeBase()
{
  this.cssClass     = "";     // CSS class name for HTML formating
  this.dynamic      = "true"; // specifies a cell with content dynamically used (content sent to server)
  this.columnLabel  = "";     // label on the table column
  this.link         = "";     // URL template for the link
  this.onClick      = "";     // template for onClick HTML method
  this.id           = "";     // template for field ID, needed to read or write data to
  this.validator    = "";     // validator must be not null in order for the cell to be displayed
  this.rtti         = "base"; // RTTI emulation in JS
}
////////////////////////////////////////////////////////////////////////////////
function ReportColumnTypeImageOnDemand()
{
  this.rtti         = REPORT_CELL_TYPE_IOD;
  this.iconRenderer = "";             // renderer used to generate the icon. If empty, iconParams treated as image/URL
  this.iconParams   = new Array();    // Icon path/image/URL   vector<ReportParamsOption>
  this.iconSequence = "";

  this.label        = "";             // label to show for the link (defined by renderer and params)
  this.renderer     = "";             // Object name used for rendering the Image On Demand
  this.params       = new Array();    // parameters passed to the renderer object for  the Image On Demand  vector<ReportParamsOption>
  this.imgSequence  = "";

  // When using a CGI call, the command is constructed in the following way:
  // /cgi-bin/<renderer> <params>
  //
  // when rendering local static pages, the renderer name is used to generate/request a render object by name (using a object factory model)
  // and <params> are passed to build the image
}

ReportColumnTypeImageOnDemand.prototype = new ReportColumnTypeBase();
////////////////////////////////////////////////////////////////////////////////
function ReportColumnTypeString()
{
  this.rtti       = REPORT_CELL_TYPE_STRING;
  this.text       = "";     // Text template for cell contents, button, input box.
  this.isButton   = false;  // If true, a button is drawn with the text in the "text" field
  this.isInput    = false;  // if True, an input box is drawn
}

ReportColumnTypeString.prototype = new ReportColumnTypeBase();
////////////////////////////////////////////////////////////////////////////////
function ReportColumnTypeStringMultiple()
{
  this.rtti       = REPORT_CELL_TYPE_STRING_MULTIPLE;
  this.linkPrefix = ""; // link filename prefix
  this.linkSuffix = ""; // link filename suffix
  this.text       = ""; // Text template for cell contents.
}

ReportColumnTypeStringMultiple.prototype = new ReportColumnTypeBase();
////////////////////////////////////////////////////////////////////////////////
function ReportColumnTypeBox()
{
  this.rtti       = REPORT_CELL_TYPE_BOX;
  this.sequences = new Array(); // Vector of several column types. vector<ReportColumnTypeBase *>
}

ReportColumnTypeBox.prototype = new ReportColumnTypeBase();
////////////////////////////////////////////////////////////////////////////////
// Table base class
////////////////////////////////////////////////////////////////////////////////
function table()
{
  // the table rows
  this.theArray;

  // the column types
  this.colTypes;

  // table type
  this.tableType;

  this.drawBorders = 1;

  // field to filter data by
  this.filterField = -1;
  // data filter
  this.filterData = "";

  // drawing exceptions:
  this.exception = "";
}
////////////////////////////////////////////////////////////////////////////////
table.prototype.clearView = function()
{
  this.colTypes = [];
}
////////////////////////////////////////////////////////////////////////////////
table.prototype.getId = function()
{
  return [];
}
////////////////////////////////////////////////////////////////////////////////
table.prototype.getSingleId = function(row)
{
  if(row < this.theArray.length)
    return this.theArray[row][0];
  return -1;
}
////////////////////////////////////////////////////////////////////////////////
table.prototype.getField = function(row, col)
{
  var aux = [];
  if(row < this.theArray.length)
    if(col < this.theArray[row].length)
      aux.push(this.theArray[row][col]);

  return aux;
}
////////////////////////////////////////////////////////////////////////////////
table.prototype.loadData = function(data)
{
  var aux = data.split("\n");
  this.theArray = new Array(aux.length);
  for(var i = 0 ; i < aux.length ; i++) {

    var aux2 = aux[i].split(";");

    var arrayItem = new Array();

    for(var j = 0 ; j < aux2.length ; j++) {
      arrayItem[j] = aux2[j];
    }
    this.theArray[i] = arrayItem;
  }
}
////////////////////////////////////////////////////////////////////////////////
table.prototype.getData = function()
{
  var params = REPORT_SERVER_PAR_GET;
  params += '&' + REPORT_SERVER_PAR_TABLE       + '&' + this.tableType;
  params += '&' + REPORT_SERVER_PAR_PROJECT_DIR + '&' + PROJECT_DIR;

  if(this.filterField >= 0) {
    params += '&' + REPORT_SERVER_PAR_FILTER_FIELD + '&' + this.filterField;
    params += '&' + REPORT_SERVER_PAR_FILTER_DATA  + '&' + this.filterData;
  }

  // this
  var This = this;
  // callback
  fnDone = function (oXML) {
    This.loadData(oXML.responseText);
  };
  // SJAX call
  var thePage = SERVER + SPS_REPORTS;
  AjaxGlobal(thePage, "GET", params, fnDone, false);
}
////////////////////////////////////////////////////////////////////////////////
// Table Protien
////////////////////////////////////////////////////////////////////////////////
function tableProtein()
{
  this.tableType = TABLE_PROTEIN_ID;
}
////////////////////////////////////////////////////////////////////////////////
tableProtein.prototype = new table();

////////////////////////////////////////////////////////////////////////////////
tableProtein.prototype.getId = function()
{
  var aux = [];
  if(this.theArray.length) {
    aux.push(this.theArray[0][0]);
  }
  return aux;
}
////////////////////////////////////////////////////////////////////////////////
tableProtein.prototype.defineView = function()
{
  // clear view
  this.clearView();

  // colTypes[0] . (CTstring) Protein name in fasta file
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_PROTEINS_0;
  auxS.text         = TAG_TABLE_PROTEINS_FIELD_NAME;
  auxS.link         = "#";
  auxS.onClick      = "javascript:loadProteinPage(" + TAG_TABLE_PROTEINS_FIELD_ID + ");";
  this.colTypes.push(auxS);

  // colTypes[1] . (CTstring) Protein description
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_PROTEINS_1;
  auxS.text         = TAG_TABLE_PROTEINS_FIELD_DESC;
  this.colTypes.push(auxS);

  // colTypes[2] . (CTstring) contigs
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_PROTEINS_2;
  auxS.text         = TAG_TABLE_PROTEINS_FIELD_CONTIGS;
  this.colTypes.push(auxS);

  // colTypes[3] . (CTstring) spectra
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_PROTEINS_3;
  auxS.text         = TAG_TABLE_PROTEINS_FIELD_SPECTRA;
  this.colTypes.push(auxS);

  // colTypes[4] . (CTstring) Amino acids
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_PROTEINS_4;
  auxS.text         = TAG_TABLE_PROTEINS_FIELD_AAS;
  this.colTypes.push(auxS);

  // colTypes[5] . (CTstring) coverage %
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_PROTEINS_5;
  auxS.text         = TAG_TABLE_PROTEINS_FIELD_COVERAGE;
  this.colTypes.push(auxS);
}
////////////////////////////////////////////////////////////////////////////////
tableProtein.prototype.defineView2 = function()
{
  // clear view
  this.clearView();

  this.exception = "ProteinException";
}
////////////////////////////////////////////////////////////////////////////////
// Table Protien Coverage
////////////////////////////////////////////////////////////////////////////////
function tableProteinCoverage()
{
  this.tableType = TABLE_COVERAGE_ID;
}
////////////////////////////////////////////////////////////////////////////////
tableProteinCoverage.prototype = new table();

////////////////////////////////////////////////////////////////////////////////
tableProteinCoverage.prototype.getId = function()
{
  var aux = [];
  if(this.theArray.length) {
    aux.push(this.theArray[0][0]);
  }
  return aux;
}
////////////////////////////////////////////////////////////////////////////////
tableProteinCoverage.prototype.defineView = function()
{
  // clear view
  this.clearView();

  this.exception = "ProteinCoverageException";
}
////////////////////////////////////////////////////////////////////////////////
tableProteinCoverage.prototype.defineView2 = function()
{
  // clear view
  this.clearView();

  this.exception = "ProteinCoverageExceptionCSV";
}
////////////////////////////////////////////////////////////////////////////////
// Table contig
////////////////////////////////////////////////////////////////////////////////
function tableContig()
{
  this.tableType = TABLE_CONTIG_ID;
}
////////////////////////////////////////////////////////////////////////////////
tableContig.prototype = new table();

////////////////////////////////////////////////////////////////////////////////
tableContig.prototype.getId = function()
{
  var aux = [];
  if(this.theArray.length) {
    aux.push(this.theArray[0][1]);
    //aux.push(this.theArray[0][0]);
  }
  return aux;
}
////////////////////////////////////////////////////////////////////////////////
tableContig.prototype.buildUpdateCall = function(inputID, suffix)
{
  var aux = "javascript:updateContig(\"" + inputID + "\", \"" + suffix + "\", \"" + IMAGE_ICON_ID_PREFIX + suffix + "\", \"" + IMAGE_ICON_CTRL_ID_PREFIX + suffix + "\", \"src\", " + TAG_TABLE_CONTIGS_FIELD_ID + ", " + TABLE_CONTIGS_FIELD_SEQ_USER + ");";
  return aux;
}
////////////////////////////////////////////////////////////////////////////////
tableContig.prototype.defineView = function()
{
  // clear view
  this.clearView();

  var conts  = "conts" + TAG_TABLE_CONTIGS_FIELD_ID;
  var inputID = "input_contig_" + TAG_INTERNAL_ROW + "_" + TAG_INTERNAL_COL;

  ///////////////////////////////////////////////////////////////////////////////
  // View for contig list
  ////////////////////////////////////////
  // Table colTypes has the following structure:

  // colTypes[0] . (CTstring) Contig index
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_CONTIGS_0;
  auxS.text         = TAG_TABLE_CONTIGS_FIELD_ID; //"«0»";
  this.colTypes.push(auxS);

  // colTypes[1] . (CTstring) Number of spectra
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_CONTIGS_1;
  auxS.text         = TAG_TABLE_CONTIGS_FIELD_SPECTRA;
  this.colTypes.push(auxS);

  // colTypes[2] . (CTIOD) Contig image
  var auxI = new ReportColumnTypeImageOnDemand();
  auxI.link              = "#";
  auxI.onClick           = "javascript:loadContigPage(" + TAG_TABLE_CONTIGS_FIELD_ID + ");";
  auxI.columnLabel       = REPORT_HEADER_CONTIGS_2;
  auxI.iconRenderer      = CONTIG_RENDERER;
  auxI.iconSequence      = TAG_TABLE_CONTIGS_FIELD_SEQ_USER;
  auxI.id                = conts;

  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_STAR,           CONTPLOT_VAL_STAR,          CONTPLOT_COND_STAR            ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_ABRUIJN,        CONTPLOT_VAL_ABRUIJN,       CONTPLOT_COND_ABRUIJN         ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_SEQS,           CONTPLOT_VAL_SEQS,          CONTPLOT_COND_SEQS            ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_CONTIG,         CONTPLOT_VAL_CONTIG,        CONTPLOT_COND_CONTIG          ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_MASS_REF,       CONTPLOT_VAL_MASS_REF,      CONTPLOT_COND_MASS_REF        ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_MASS_HOM,       CONTPLOT_VAL_MASS_HOM,      CONTPLOT_COND_MASS_HOM        ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_OFF_REF,        CONTPLOT_VAL_OFF_REF,       CONTPLOT_COND_OFF_REF         ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_OFF_HOM,        CONTPLOT_VAL_OFF_HOM,       CONTPLOT_COND_OFF_HOM         ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_REVERSE,        CONTPLOT_VAL_REVERSE,       CONTPLOT_COND_REVERSE         ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_TARGET,         CONTPLOT_VAL_TARGET,        CONTPLOT_COND_TARGET          ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_ENCODING,       CONTPLOT_VAL_ENCODING,      CONTPLOT_COND_ENCODING        ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_NO_TITLE,       CONTPLOT_VAL_NO_TITLE,      CONTPLOT_COND_NO_TITLE        ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_ZOOM,           CONTPLOT_VAL_ZOOM,          CONTPLOT_COND_ZOOM            ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_SEQ_REFERENCE,  CONTPLOT_VAL_SEQ_REFERENCE, CONTPLOT_COND_SEQ_REFERENCE   ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_SEQ_HOMOLOG,    CONTPLOT_VAL_SEQ_HOMOLOG,   CONTPLOT_COND_SEQ_HOMOLOG     ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_SEQ_DENOVO,     CONTPLOT_VAL_SEQ_DENOVO,    CONTPLOT_COND_SEQ_DENOVO      ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_SEQ_USER,       CONTPLOT_VAL_SEQ_USER,      CONTPLOT_COND_SEQ_USER,     0 ) );

  this.colTypes.push(auxI);


  // colTypes[2] . (CTseqsBox)
  var auxB = new ReportColumnTypeBox();
  auxB.columnLabel  = REPORT_HEADER_CONTIGS_3;

  auxI = new ReportColumnTypeImageOnDemand();
  auxI.columnLabel  = REPORT_SEQ_NAME_REFERENCE;
  auxI.label        = TAG_TABLE_CONTIGS_FIELD_SEQ_REF;
  auxI.validator    = TAG_TABLE_CONTIGS_FIELD_SEQ_REF;
  auxB.sequences.push(auxI);

  var auxI = new ReportColumnTypeImageOnDemand();
  auxI.columnLabel  = REPORT_SEQ_NAME_HOMOLOG;
  auxI.label        = TAG_TABLE_CONTIGS_FIELD_SEQ_HOM;
  auxI.validator    = TAG_TABLE_CONTIGS_FIELD_SEQ_HOM;
  auxB.sequences.push(auxI);

  var auxI = new ReportColumnTypeImageOnDemand();
  auxI.columnLabel  = REPORT_SEQ_NAME_DENOVO;
  auxI.label        = TAG_TABLE_CONTIGS_FIELD_SEQ_NOVO;
  auxI.validator    = TAG_TABLE_CONTIGS_FIELD_SEQ_NOVO;
  auxB.sequences.push(auxI);

  var auxI = new ReportColumnTypeImageOnDemand();
  auxI.dynamic      = true;
  auxI.columnLabel  = REPORT_SEQ_NAME_USER;
  auxI.label        = TAG_TABLE_CONTIGS_FIELD_SEQ_USER;
  auxI.id           = conts;
  auxB.sequences.push(auxI);

  var auxS = new ReportColumnTypeString();
  auxS.isInput      = true;
  auxS.dynamic      = true;
  auxS.id           = inputID;
  auxB.sequences.push(auxS);

  var auxS = new ReportColumnTypeString();
  auxS.isButton     = true;
  auxS.dynamic      = true;
  auxS.text         = REPORT_BUTTON_UPDATE;
  auxS.onClick      = this.buildUpdateCall(inputID, conts);
  auxB.sequences.push(auxS);

  this.colTypes.push(auxB);


  // colTypes[3] . (CTstring) protein
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_CONTIGS_4;
  auxS.text         = TAG_TABLE_CONTIGS_FIELD_PROTEIN_NAME + "<br><br>" + TAG_TABLE_CONTIGS_FIELD_PROTEIN_DESC;
  auxS.link         = "#";
  auxS.onClick      = "javascript:loadProteinPage(" + TAG_TABLE_CONTIGS_FIELD_PROTEIN + ");";
  this.colTypes.push(auxS);
}
////////////////////////////////////////////////////////////////////////////////
tableContig.prototype.defineView2 = function()
{
  // clear view
  this.clearView();

  var cont  = "cont" + TAG_TABLE_CONTIGS_FIELD_ID;
  var inputID = "input_contig_" + TAG_INTERNAL_ROW + "_" + TAG_INTERNAL_COL;

  var auxB = new ReportColumnTypeBox();
  auxB.columnLabel  = "";

  // the image
  auxI = new ReportColumnTypeImageOnDemand();
  auxI.columnLabel       = "";
  auxI.iconRenderer      = CONTIG_RENDERER;
  auxI.iconSequence      = TAG_TABLE_CONTIGS_FIELD_SEQ_USER;
  auxI.id                = cont;

  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_STAR,           CONTPLOT_VAL_STAR,          CONTPLOT_COND_STAR            ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_ABRUIJN,        CONTPLOT_VAL_ABRUIJN,       CONTPLOT_COND_ABRUIJN         ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_SEQS,           CONTPLOT_VAL_SEQS,          CONTPLOT_COND_SEQS            ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_TITLE,          CONTPLOT_VAL_TITLE,         CONTPLOT_COND_TITLE           ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_CONTIG,         CONTPLOT_VAL_CONTIG,        CONTPLOT_COND_CONTIG          ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_MASS_REF,       CONTPLOT_VAL_MASS_REF,      CONTPLOT_COND_MASS_REF        ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_MASS_HOM,       CONTPLOT_VAL_MASS_HOM,      CONTPLOT_COND_MASS_HOM        ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_OFF_REF,        CONTPLOT_VAL_OFF_REF,       CONTPLOT_COND_OFF_REF         ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_OFF_HOM,        CONTPLOT_VAL_OFF_HOM,       CONTPLOT_COND_OFF_HOM         ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_REVERSE,        CONTPLOT_VAL_REVERSE,       CONTPLOT_COND_REVERSE         ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_TARGET,         CONTPLOT_VAL_TARGET,        CONTPLOT_COND_TARGET          ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_ENCODING,       CONTPLOT_VAL_ENCODING,      CONTPLOT_COND_ENCODING        ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_SEQ_REFERENCE,  CONTPLOT_VAL_SEQ_REFERENCE, CONTPLOT_COND_SEQ_REFERENCE   ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_SEQ_HOMOLOG,    CONTPLOT_VAL_SEQ_HOMOLOG,   CONTPLOT_COND_SEQ_HOMOLOG     ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_SEQ_DENOVO,     CONTPLOT_VAL_SEQ_DENOVO,    CONTPLOT_COND_SEQ_DENOVO      ) );
  auxI.iconParams.push( new ReportParamsOption(CONTPLOT_PAR_SEQ_USER,       CONTPLOT_VAL_SEQ_USER,      CONTPLOT_COND_SEQ_USER,     0 ) );

  auxB.sequences.push(auxI);


  // the input sequence
  var auxS = new ReportColumnTypeString();
  auxS.isInput      = true;
  auxS.dynamic      = true;
  auxS.id           = inputID;
  auxB.sequences.push(auxS);

  // the "update" button
  var auxS = new ReportColumnTypeString();
  auxS.isButton     = true;
  auxS.dynamic      = true;
  auxS.text         = REPORT_BUTTON_UPDATE;
  auxS.onClick      = this.buildUpdateCall(inputID, cont);
  auxB.sequences.push(auxS);

  this.colTypes.push(auxB);
}
////////////////////////////////////////////////////////////////////////////////
// Table Cluster consensus
////////////////////////////////////////////////////////////////////////////////
function tableCluster()
{
  this.tableType = TABLE_CLUSTER_ID;
}
////////////////////////////////////////////////////////////////////////////////
tableCluster.prototype = new table();

////////////////////////////////////////////////////////////////////////////////
tableCluster.prototype.getId = function()
{
  var aux = [];
  if(this.theArray.length) {
    aux.push(this.theArray[0][2]);
    aux.push(this.theArray[0][1]);
    //aux.push(this.theArray[0][0]);
  }
  return aux;
}
////////////////////////////////////////////////////////////////////////////////
tableCluster.prototype.buildUpdateCall = function(inputID, suffix, icon)
{
  var aux = "";

  if(typeof(icon) != 'undefined' && icon) {
    aux = "javascript:updateCluster(\"" + inputID + "\", \"" + suffix + "\", \"" + IMAGE_LARGE_ID_PREFIX + suffix + "\", \"" + IMAGE_LARGE_CTRL_ID_PREFIX + suffix + "\", \"href\", " + TAG_TABLE_CLUSTER_FIELD_ID + "," + TABLE_CLUSTER_FIELD_USER + ", false);";
    aux += "triggerClick(\"" + IMAGE_LARGE_ID_PREFIX + suffix + "\", false);";
  } else {
    aux = "javascript:updateCluster(\"" + inputID + "\", null, \"" + IMAGE_ICON_ID_PREFIX + suffix + "\", \"" + IMAGE_ICON_CTRL_ID_PREFIX + suffix + "\", \"src\", " + TAG_TABLE_CLUSTER_FIELD_ID + "," + TABLE_CLUSTER_FIELD_USER + ", true);";
  }

  return aux;
}
////////////////////////////////////////////////////////////////////////////////
tableCluster.prototype.defineView = function()
{
  // clear view
  this.clearView();

  var clusts  = "clusts" + TAG_TABLE_CLUSTER_FIELD_ID;
  var clustR  = "clustR" + TAG_TABLE_CLUSTER_FIELD_ID;
  var clustH  = "clustH" + TAG_TABLE_CLUSTER_FIELD_ID;
  var clustN  = "clustN" + TAG_TABLE_CLUSTER_FIELD_ID;
  var clustU  = "clustU" + TAG_TABLE_CLUSTER_FIELD_ID;
  var inputID = "input_cluster_" + TAG_INTERNAL_ROW + "_" + TAG_INTERNAL_COL;

  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_CLUSTERS_0;
  auxS.text         = TAG_TABLE_CLUSTER_FIELD_ID;
  this.colTypes.push(auxS);

  var auxI = new ReportColumnTypeImageOnDemand();
  auxI.link            = "#";
  auxI.onClick         = "javascript:loadClusterPage(" + TAG_TABLE_CLUSTER_FIELD_ID + ");";
  auxI.columnLabel     = REPORT_HEADER_CLUSTERS_1;
  auxI.label           = "";
  auxI.iconRenderer    = SPECTRUM_RENDERER;
  auxI.validator       = TAG_CLUSTER_PEPTIDE_ALL;
  auxI.id              = clusts;

  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_PKLBIN,   SPECPLOT_VAL_CLUSTER_PKLBIN,      SPECPLOT_COND_CLUSTER_PKLBIN      ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_SPECTURM, SPECPLOT_VAL_CLUSTER_SPECTURM,    SPECPLOT_COND_CLUSTER_SPECTURM    ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_PEPTIDE,  SPECPLOT_VAL_CLUSTER_PEPTIDE_ALL, SPECPLOT_COND_CLUSTER_PEPTIDE_ALL ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_TARGET,   SPECPLOT_VAL_TARGET,              SPECPLOT_COND_TARGET              ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_ENCODING, SPECPLOT_VAL_ENCODING,            SPECPLOT_COND_ENCODING            ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_ZOOM,     SPECPLOT_VAL_CLUSTER_ZOOM,        SPECPLOT_COND_CLUSTER_ZOOM        ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_NOTITLE,  SPECPLOT_VAL_CLUSTER_NOTITLE,     SPECPLOT_COND_CLUSTER_NOTITLE     ) );

  this.colTypes.push(auxI);


  // colTypes[2] . (CTseqsBox)
  var auxB = new ReportColumnTypeBox();
  auxB.columnLabel  = REPORT_HEADER_CLUSTERS_2;

  var auxI = new ReportColumnTypeImageOnDemand();
  auxI.columnLabel  = REPORT_SEQ_NAME_REFERENCE;
  auxI.label        = TAG_TABLE_CLUSTER_FIELD_REFERENCE;
  auxI.renderer     = SPECTRUM_RENDERER;
  auxI.validator    = TAG_TABLE_CLUSTER_FIELD_REFERENCE;
  auxI.id           = clustR;

  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TITLE,    SPECPLOT_VAL_CLUSTER_TITLE,       SPECPLOT_COND_CLUSTER_TITLE       ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PKLBIN,   SPECPLOT_VAL_CLUSTER_PKLBIN,      SPECPLOT_COND_CLUSTER_PKLBIN      ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_SPECTURM, SPECPLOT_VAL_CLUSTER_SPECTURM,    SPECPLOT_COND_CLUSTER_SPECTURM    ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PEPTIDE,  SPECPLOT_VAL_CLUSTER_PEPTIDE_REF, SPECPLOT_COND_CLUSTER_PEPTIDE_REF ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TARGET,   SPECPLOT_VAL_TARGET,              SPECPLOT_COND_TARGET              ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_ENCODING, SPECPLOT_VAL_ENCODING,            SPECPLOT_COND_ENCODING            ) );

  auxB.sequences.push(auxI);


  var auxI = new ReportColumnTypeImageOnDemand();
  auxI.columnLabel  = REPORT_SEQ_NAME_HOMOLOG;
  auxI.label        = TAG_TABLE_CLUSTER_FIELD_HOMOLOG;
  auxI.renderer     = SPECTRUM_RENDERER;
  auxI.validator    = TAG_TABLE_CLUSTER_FIELD_HOMOLOG;
  auxI.id           = clustH;

  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TITLE,    SPECPLOT_VAL_CLUSTER_TITLE,       SPECPLOT_COND_CLUSTER_TITLE       ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PKLBIN,   SPECPLOT_VAL_CLUSTER_PKLBIN,      SPECPLOT_COND_CLUSTER_PKLBIN      ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_SPECTURM, SPECPLOT_VAL_CLUSTER_SPECTURM,    SPECPLOT_COND_CLUSTER_SPECTURM    ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PEPTIDE,  SPECPLOT_VAL_CLUSTER_PEPTIDE_HOM, SPECPLOT_COND_CLUSTER_PEPTIDE_HOM ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TARGET,   SPECPLOT_VAL_TARGET,              SPECPLOT_COND_TARGET              ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_ENCODING, SPECPLOT_VAL_ENCODING,            SPECPLOT_COND_ENCODING            ) );

  auxB.sequences.push(auxI);

  var auxI = new ReportColumnTypeImageOnDemand();
  auxI.columnLabel  = REPORT_SEQ_NAME_DENOVO;
  auxI.label        = TAG_TABLE_CLUSTER_FIELD_DENOVO;
  auxI.renderer     = SPECTRUM_RENDERER;
  auxI.validator    = TAG_TABLE_CLUSTER_FIELD_DENOVO;
  auxI.id           = clustN;

  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TITLE,    SPECPLOT_VAL_CLUSTER_TITLE,         SPECPLOT_COND_CLUSTER_TITLE         ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PKLBIN,   SPECPLOT_VAL_CLUSTER_PKLBIN,        SPECPLOT_COND_CLUSTER_PKLBIN        ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_SPECTURM, SPECPLOT_VAL_CLUSTER_SPECTURM,      SPECPLOT_COND_CLUSTER_SPECTURM      ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PEPTIDE,  SPECPLOT_VAL_CLUSTER_PEPTIDE_NOVO,  SPECPLOT_COND_CLUSTER_PEPTIDE_NOVO  ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TARGET,   SPECPLOT_VAL_TARGET,                SPECPLOT_COND_TARGET                ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_ENCODING, SPECPLOT_VAL_ENCODING,              SPECPLOT_COND_ENCODING              ) );

  auxB.sequences.push(auxI);


  var auxI = new ReportColumnTypeImageOnDemand();
  auxI.dynamic      = true;
  auxI.columnLabel  = REPORT_SEQ_NAME_USER;
  auxI.label        = TAG_TABLE_CLUSTER_FIELD_USER;
  auxI.renderer     = SPECTRUM_RENDERER;
  auxI.imgSequence  = TAG_TABLE_CLUSTER_FIELD_USER;
  auxI.id           = clustU;
  //auxI.validator    = TAG_TABLE_CLUSTER_FIELD_USER;

  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TITLE,    SPECPLOT_VAL_CLUSTER_TITLE,         SPECPLOT_COND_CLUSTER_TITLE           ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PKLBIN,   SPECPLOT_VAL_CLUSTER_PKLBIN,        SPECPLOT_COND_CLUSTER_PKLBIN          ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_SPECTURM, SPECPLOT_VAL_CLUSTER_SPECTURM,      SPECPLOT_COND_CLUSTER_SPECTURM        ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PEPTIDE,  SPECPLOT_VAL_CLUSTER_PEPTIDE_USER,  SPECPLOT_COND_CLUSTER_PEPTIDE_USER, 0 ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TARGET,   SPECPLOT_VAL_TARGET,                SPECPLOT_COND_TARGET                  ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_ENCODING, SPECPLOT_VAL_ENCODING,              SPECPLOT_COND_ENCODING                ) );

  auxB.sequences.push(auxI);

  var auxS = new ReportColumnTypeString();
  auxS.isInput      = true;
  auxS.dynamic      = true;
  auxS.id           = inputID;
  auxB.sequences.push(auxS);

  var auxS = new ReportColumnTypeString();
  auxS.isButton     = true;
  auxS.dynamic      = true;
  auxS.text         = REPORT_BUTTON_UPDATE;
  auxS.onClick      = this.buildUpdateCall(inputID, clustU, true);
  auxB.sequences.push(auxS);

  this.colTypes.push(auxB);


// colTypes[2] . (CTstring) mass
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_CLUSTERS_3;
  auxS.text         = TAG_TABLE_CLUSTER_FIELD_MASS;
  this.colTypes.push(auxS);

// colTypes[3] . (CTstring) charge
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_CLUSTERS_4;
  auxS.text         = TAG_TABLE_CLUSTER_FIELD_CHARGE;
  this.colTypes.push(auxS);

// colTypes[4] . (CTstring) B%
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_CLUSTERS_5;
  auxS.text         = TAG_TABLE_CLUSTER_FIELD_B_PER;
  this.colTypes.push(auxS);

// colTypes[5] . (CTstring) Y%
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_CLUSTERS_6;
  auxS.text         = TAG_TABLE_CLUSTER_FIELD_Y_PER;
  this.colTypes.push(auxS);

// colTypes[6] . (CTstring) BY intensity %
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_CLUSTERS_7;
  auxS.text         = TAG_TABLE_CLUSTER_FIELD_BY_INT;
  this.colTypes.push(auxS);
}
////////////////////////////////////////////////////////////////////////////////
tableCluster.prototype.defineView2 = function()
{
  // clear view
  this.clearView();

  var clust   = "clust" + TAG_TABLE_CLUSTER_FIELD_ID;
  var inputID = "input_cluster_" + TAG_INTERNAL_ROW + "_" + TAG_INTERNAL_COL;

  var auxB = new ReportColumnTypeBox();
  auxB.columnLabel  = "Sequences";

  var auxI = new ReportColumnTypeImageOnDemand();
  auxI.columnLabel     = "";
  auxI.label           = "";
  auxI.iconRenderer    = SPECTRUM_RENDERER;
  auxI.id              = clust;

  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_TITLE,    SPECPLOT_VAL_CLUSTER_TITLE,       SPECPLOT_COND_CLUSTER_TITLE           ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_PKLBIN,   SPECPLOT_VAL_CLUSTER_PKLBIN,      SPECPLOT_COND_CLUSTER_PKLBIN          ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_SPECTURM, SPECPLOT_VAL_CLUSTER_SPECTURM,    SPECPLOT_COND_CLUSTER_SPECTURM        ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_TARGET,   SPECPLOT_VAL_TARGET,              SPECPLOT_COND_TARGET                  ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_ENCODING, SPECPLOT_VAL_ENCODING,            SPECPLOT_COND_ENCODING                ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_PEPTIDE,  SPECPLOT_VAL_CLUSTER_PEPTIDE_ALL, SPECPLOT_COND_CLUSTER_PEPTIDE_ALL,  0 ) );

  auxB.sequences.push(auxI);


  var auxS = new ReportColumnTypeString();
  auxS.isInput      = true;
  auxS.dynamic      = true;
  auxS.id           = inputID;
  auxB.sequences.push(auxS);

  var auxS = new ReportColumnTypeString();
  auxS.isButton     = true;
  auxS.dynamic      = true;
  auxS.text         = REPORT_BUTTON_UPDATE;
  auxS.onClick      = this.buildUpdateCall(inputID, clust, false);
  auxB.sequences.push(auxS);

  this.colTypes.push(auxB);
}
////////////////////////////////////////////////////////////////////////////////
// Table Input Spectra
////////////////////////////////////////////////////////////////////////////////
function tableInputSpectra()
{
  this.tableType = TABLE_SPECTRA_ID;
}
////////////////////////////////////////////////////////////////////////////////
tableInputSpectra.prototype = new table();

////////////////////////////////////////////////////////////////////////////////
tableInputSpectra.prototype.getId = function()
{
  var aux = [];
  if(this.theArray.length) {
    //aux.push(this.theArray[0][2]);
    //aux.push(this.theArray[0][1]);
    //aux.push(this.theArray[0][0]);
  }
  return aux;
}
////////////////////////////////////////////////////////////////////////////////
tableInputSpectra.prototype.getSingleId = function(row)
{
  if(row < this.theArray.length)
    return this.theArray[row][2];
  return -1;
}
////////////////////////////////////////////////////////////////////////////////
tableInputSpectra.prototype.buildUpdateCall = function(inputID, suffix)
{
  var aux = "javascript:updateSpectra(\"" + inputID + "\", \"" + suffix + "\", \"" + IMAGE_LARGE_ID_PREFIX + suffix + "\", \"" + IMAGE_LARGE_CTRL_ID_PREFIX + suffix + "\", \"href\", " + TAG_TABLE_SPECTRA_FIELD_ID + ", " + TABLE_SPECTRA_FIELD_SEQ_USER + ", false);";
  aux += "triggerClick(\"" + IMAGE_LARGE_ID_PREFIX + suffix + "\");";
  return aux;
}
////////////////////////////////////////////////////////////////////////////////
tableInputSpectra.prototype.defineView = function()
{
  // clear view
  this.clearView();

  var specs   = "specs" + TAG_TABLE_SPECTRA_FIELD_ID;
  var specR   = "specR" + TAG_TABLE_SPECTRA_FIELD_ID;
  var specH   = "specH" + TAG_TABLE_SPECTRA_FIELD_ID;
  var specN   = "specN" + TAG_TABLE_SPECTRA_FIELD_ID
  var specU   = "specU" + TAG_TABLE_SPECTRA_FIELD_ID;
  var inputID = "input_spectra_" + TAG_INTERNAL_ROW + "_" + TAG_INTERNAL_COL;

  // colTypes[0] . (CTstring) Spectrum index
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel   = REPORT_HEADER_SPECTRA_0;
  auxS.text          = TAG_TABLE_SPECTRA_FIELD_IDX;
  this.colTypes.push(auxS);

  var auxS = new ReportColumnTypeString();
  auxS.columnLabel   = REPORT_HEADER_SPECTRA_1;
  auxS.text          = TAG_TABLE_SPECTRA_FIELD_SCAN;
  this.colTypes.push(auxS);

  // colTypes[1] . (CTseqsBox)
  var auxB = new ReportColumnTypeBox();
  auxB.columnLabel   = REPORT_HEADER_SPECTRA_2;
  auxB.link          = "javascript:loadClusterPage(" + TAG_TABLE_SPECTRA_FIELD_CLUSTER + ");";

  var auxI = new ReportColumnTypeImageOnDemand();
  auxI.iconRenderer  = SPECTRUM_RENDERER;
  auxI.validator     = TAG_SPECTRA_PEPTIDE_ALL;
  auxI.id            = specs;

  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_PKLBIN,   SPECPLOT_VAL_SPECTRA_PKLBIN,      SPECPLOT_COND_SPECTRA_PKLBIN      ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_SPECTURM, SPECPLOT_VAL_SPECTRA_SPECTURM,    SPECPLOT_COND_SPECTRA_SPECTURM    ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_PEPTIDE,  SPECPLOT_VAL_SPECTRA_PEPTIDE_ALL, SPECPLOT_COND_SPECTRA_PEPTIDE_ALL ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_TARGET,   SPECPLOT_VAL_TARGET,              SPECPLOT_COND_TARGET              ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_ENCODING, SPECPLOT_VAL_ENCODING,            SPECPLOT_COND_ENCODING            ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_ZOOM,     SPECPLOT_VAL_SPECTRA_ZOOM,        SPECPLOT_COND_SPECTRA_ZOOM        ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_NOTITLE,  SPECPLOT_VAL_SPECTRA_NOTITLE,     SPECPLOT_COND_SPECTRA_NOTITLE     ) );

  auxB.sequences.push(auxI);


  var auxS = new ReportColumnTypeString();
  auxS.text          = TAG_TABLE_SPECTRA_FIELD_INPUT_FILENAME;
  auxB.sequences.push(auxS);

  this.colTypes.push(auxB);


  // colTypes[2] . (CTseqsBox)
  var auxB = new ReportColumnTypeBox();
  auxB.columnLabel   = REPORT_HEADER_SPECTRA_3;

  var auxI = new ReportColumnTypeImageOnDemand();
  auxI.columnLabel  = REPORT_SEQ_NAME_REFERENCE;
  auxI.label        = TAG_TABLE_SPECTRA_FIELD_SEQ_REFERENCE;
  auxI.renderer     = SPECTRUM_RENDERER;
  auxI.validator    = TAG_TABLE_SPECTRA_FIELD_SEQ_REFERENCE;
  auxI.id           = specR;

  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PKLBIN,   SPECPLOT_VAL_SPECTRA_PKLBIN,      SPECPLOT_COND_SPECTRA_PKLBIN      ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_SPECTURM, SPECPLOT_VAL_SPECTRA_SPECTURM,    SPECPLOT_COND_SPECTRA_SPECTURM    ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PEPTIDE,  SPECPLOT_VAL_SPECTRA_PEPTIDE_REF, SPECPLOT_COND_SPECTRA_PEPTIDE_REF ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TARGET,   SPECPLOT_VAL_TARGET,              SPECPLOT_COND_TARGET              ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_ENCODING, SPECPLOT_VAL_ENCODING,            SPECPLOT_COND_ENCODING            ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TITLE,    SPECPLOT_VAL_SPECTRA_TITLE,       SPECPLOT_COND_SPECTRA_TITLE       ) );

  auxB.sequences.push(auxI);


  var auxI = new ReportColumnTypeImageOnDemand();
  auxI.columnLabel  = REPORT_SEQ_NAME_HOMOLOG;
  auxI.label        = TAG_TABLE_SPECTRA_FIELD_SEQ_HOMOLOG;
  auxI.renderer     = SPECTRUM_RENDERER;
  auxI.validator    = TAG_TABLE_SPECTRA_FIELD_SEQ_HOMOLOG;
  auxI.id           = specH;

  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PKLBIN,   SPECPLOT_VAL_SPECTRA_PKLBIN,      SPECPLOT_COND_SPECTRA_PKLBIN      ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_SPECTURM, SPECPLOT_VAL_SPECTRA_SPECTURM,    SPECPLOT_COND_SPECTRA_SPECTURM    ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PEPTIDE,  SPECPLOT_VAL_SPECTRA_PEPTIDE_HOM, SPECPLOT_COND_SPECTRA_PEPTIDE_HOM ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TARGET,   SPECPLOT_VAL_TARGET,              SPECPLOT_COND_TARGET              ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_ENCODING, SPECPLOT_VAL_ENCODING,            SPECPLOT_COND_ENCODING            ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TITLE,    SPECPLOT_VAL_SPECTRA_TITLE,       SPECPLOT_COND_SPECTRA_TITLE       ) );

  auxB.sequences.push(auxI);


  var auxI = new ReportColumnTypeImageOnDemand();
  auxI.columnLabel  = REPORT_SEQ_NAME_DENOVO;
  auxI.label        = TAG_TABLE_SPECTRA_FIELD_SEQ_DENOVO;
  auxI.renderer     = SPECTRUM_RENDERER;
  auxI.validator    = TAG_TABLE_SPECTRA_FIELD_SEQ_DENOVO;
  auxI.id           = specN;

  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PKLBIN,   SPECPLOT_VAL_SPECTRA_PKLBIN,        SPECPLOT_COND_SPECTRA_PKLBIN        ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_SPECTURM, SPECPLOT_VAL_SPECTRA_SPECTURM,      SPECPLOT_COND_SPECTRA_SPECTURM      ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PEPTIDE,  SPECPLOT_VAL_SPECTRA_PEPTIDE_NOVO,  SPECPLOT_COND_SPECTRA_PEPTIDE_NOVO  ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TARGET,   SPECPLOT_VAL_TARGET,                SPECPLOT_COND_TARGET                ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_ENCODING, SPECPLOT_VAL_ENCODING,              SPECPLOT_COND_ENCODING              ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TITLE,    SPECPLOT_VAL_SPECTRA_TITLE,         SPECPLOT_COND_SPECTRA_TITLE         ) );

  auxB.sequences.push(auxI);


  var auxI = new ReportColumnTypeImageOnDemand();
  auxI.dynamic      = true;
  auxI.columnLabel  = REPORT_SEQ_NAME_USER;
  auxI.label        = TAG_TABLE_SPECTRA_FIELD_SEQ_USER;
  auxI.renderer     = SPECTRUM_RENDERER;
  auxI.imgSequence  = TAG_TABLE_SPECTRA_FIELD_SEQ_USER;
  auxI.id           = specU;
  //auxI.validator    = TAG_TABLE_SPECTRA_FIELD_SEQ_USER;

  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PKLBIN,   SPECPLOT_VAL_SPECTRA_PKLBIN,        SPECPLOT_COND_SPECTRA_PKLBIN          ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_SPECTURM, SPECPLOT_VAL_SPECTRA_SPECTURM,      SPECPLOT_COND_SPECTRA_SPECTURM        ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PEPTIDE,  SPECPLOT_VAL_SPECTRA_PEPTIDE_USER,  SPECPLOT_COND_SPECTRA_PEPTIDE_USER, 0 ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TARGET,   SPECPLOT_VAL_TARGET,                SPECPLOT_COND_TARGET                  ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_ENCODING, SPECPLOT_VAL_ENCODING,              SPECPLOT_COND_ENCODING                ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TITLE,    SPECPLOT_VAL_SPECTRA_TITLE,         SPECPLOT_COND_SPECTRA_TITLE           ) );

  auxB.sequences.push(auxI);


  var auxS = new ReportColumnTypeString();
  auxS.dynamic      = true;
  auxS.isInput      = true;
  auxS.id           = inputID;
  auxB.sequences.push(auxS);

  var auxS = new ReportColumnTypeString();
  auxS.dynamic      = true;
  auxS.isButton     = true;
  auxS.text         = REPORT_BUTTON_UPDATE;
  auxS.onClick      = this.buildUpdateCall(inputID, specU);
  auxB.sequences.push(auxS);

  this.colTypes.push(auxB);


  // colTypes[3] . (CTstring) mass
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_SPECTRA_4;
  auxS.text         = TAG_TABLE_SPECTRA_FIELD_PROTEIN_NAME;
  this.colTypes.push(auxS);

  // colTypes[3] . (CTstring) mass
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_SPECTRA_5;
  auxS.text         = TAG_TABLE_SPECTRA_FIELD_MASS;
  this.colTypes.push(auxS);

  // colTypes[4] . (CTstring) charge
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_SPECTRA_6;
  auxS.text         = TAG_TABLE_SPECTRA_FIELD_CHARGE;
  this.colTypes.push(auxS);

  // colTypes[5] . (CTstring) B%
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_SPECTRA_7;
  auxS.text         = TAG_TABLE_SPECTRA_FIELD_B_PER;
  this.colTypes.push(auxS);

  // colTypes[6] . (CTstring) Y%
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_SPECTRA_8;
  auxS.text         = TAG_TABLE_SPECTRA_FIELD_Y_PER;
  this.colTypes.push(auxS);

  // colTypes[7] . (CTstring) BY intensity %
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_SPECTRA_9;
  auxS.text         = TAG_TABLE_SPECTRA_FIELD_BY_INTENSITY;
  this.colTypes.push(auxS);
}
////////////////////////////////////////////////////////////////////////////////
tableInputSpectra.prototype.defineView2 = function()
{
  // clear view
  this.clearView();

  var specs   = "specs" + TAG_TABLE_SPECTRA_FIELD_ID;
  var specR   = "specR" + TAG_TABLE_SPECTRA_FIELD_ID;
  var specH   = "specH" + TAG_TABLE_SPECTRA_FIELD_ID;
  var specN   = "specN" + TAG_TABLE_SPECTRA_FIELD_ID
  var specU   = "specU" + TAG_TABLE_SPECTRA_FIELD_ID;
  var inputID = "input_spectra_" + TAG_INTERNAL_ROW + "_" + TAG_INTERNAL_COL;

  // colTypes[0] . (CTstring) Spectrum index
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel   = REPORT_HEADER_SPECTRA_0;
  auxS.text          = TAG_TABLE_SPECTRA_FIELD_IDX;
  this.colTypes.push(auxS);

  // colTypes[0] . (CTstring) Spectrum index
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel   = REPORT_HEADER_SPECTRA_1;
  auxS.text          = TAG_TABLE_SPECTRA_FIELD_SCAN;
  this.colTypes.push(auxS);

  // colTypes[1] . (CTseqsBox)
  var auxB = new ReportColumnTypeBox();
  auxB.columnLabel   = REPORT_HEADER_SPECTRA_2;
  auxB.link          = "javascript:loadClusterPage(" + TAG_TABLE_SPECTRA_FIELD_CLUSTER + ")";

  var auxI = new ReportColumnTypeImageOnDemand();
  auxI.iconRenderer  = SPECTRUM_RENDERER;
  auxI.validator     = TAG_SPECTRA_PEPTIDE_ALL;
  auxI.id            = specs;

  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_PKLBIN,   SPECPLOT_VAL_SPECTRA_PKLBIN,       SPECPLOT_COND_SPECTRA_PKLBIN       ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_SPECTURM, SPECPLOT_VAL_SPECTRA_SPECTURM,     SPECPLOT_COND_SPECTRA_SPECTURM     ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_PEPTIDE,  SPECPLOT_VAL_SPECTRA_PEPTIDE_ALL,  SPECPLOT_COND_SPECTRA_PEPTIDE_ALL  ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_TARGET,   SPECPLOT_VAL_TARGET,               SPECPLOT_COND_TARGET               ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_ENCODING, SPECPLOT_VAL_ENCODING,             SPECPLOT_COND_ENCODING             ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_ZOOM,     SPECPLOT_VAL_SPECTRA_ZOOM,         SPECPLOT_COND_SPECTRA_ZOOM         ) );
  auxI.iconParams.push( new ReportParamsOption(SPECPLOT_PAR_NOTITLE,  SPECPLOT_VAL_SPECTRA_NOTITLE,      SPECPLOT_COND_SPECTRA_NOTITLE      ) );

  auxB.sequences.push(auxI);


  var auxS = new ReportColumnTypeString();
  auxS.text          = TAG_TABLE_SPECTRA_FIELD_INPUT_FILENAME;
  auxB.sequences.push(auxS);

  this.colTypes.push(auxB);


  // colTypes[2] . (CTseqsBox)
  var auxB = new ReportColumnTypeBox();
  auxB.columnLabel   = REPORT_HEADER_SPECTRA_3;

  var auxI = new ReportColumnTypeImageOnDemand();
  auxI.columnLabel  = REPORT_SEQ_NAME_REFERENCE;
  auxI.label        = TAG_TABLE_SPECTRA_FIELD_SEQ_REFERENCE;
  auxI.renderer     = SPECTRUM_RENDERER;
  auxI.validator    = TAG_TABLE_SPECTRA_FIELD_SEQ_REFERENCE;
  auxI.id           = specR;

  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PKLBIN,   SPECPLOT_VAL_SPECTRA_PKLBIN,      SPECPLOT_COND_SPECTRA_PKLBIN      ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_SPECTURM, SPECPLOT_VAL_SPECTRA_SPECTURM,    SPECPLOT_COND_SPECTRA_SPECTURM    ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PEPTIDE,  SPECPLOT_VAL_SPECTRA_PEPTIDE_REF, SPECPLOT_COND_SPECTRA_PEPTIDE_REF ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TARGET,   SPECPLOT_VAL_TARGET,              SPECPLOT_COND_TARGET              ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_ENCODING, SPECPLOT_VAL_ENCODING,            SPECPLOT_COND_ENCODING            ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TITLE,    SPECPLOT_VAL_SPECTRA_TITLE,       SPECPLOT_COND_SPECTRA_TITLE       ) );

  auxB.sequences.push(auxI);


  var auxI = new ReportColumnTypeImageOnDemand();
  auxI.columnLabel  = REPORT_SEQ_NAME_HOMOLOG;
  auxI.label        = TAG_TABLE_SPECTRA_FIELD_SEQ_HOMOLOG;
  auxI.renderer     = SPECTRUM_RENDERER;
  auxI.validator    = TAG_TABLE_SPECTRA_FIELD_SEQ_HOMOLOG;
  auxI.id           = specH;

  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PKLBIN,   SPECPLOT_VAL_SPECTRA_PKLBIN,      SPECPLOT_COND_SPECTRA_PKLBIN      ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_SPECTURM, SPECPLOT_VAL_SPECTRA_SPECTURM,    SPECPLOT_COND_SPECTRA_SPECTURM    ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PEPTIDE,  SPECPLOT_VAL_SPECTRA_PEPTIDE_HOM, SPECPLOT_COND_SPECTRA_PEPTIDE_HOM ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TARGET,   SPECPLOT_VAL_TARGET,              SPECPLOT_COND_TARGET              ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_ENCODING, SPECPLOT_VAL_ENCODING,            SPECPLOT_COND_ENCODING            ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TITLE,    SPECPLOT_VAL_SPECTRA_TITLE,       SPECPLOT_COND_SPECTRA_TITLE       ) );

  auxB.sequences.push(auxI);


  var auxI = new ReportColumnTypeImageOnDemand();
  auxI.columnLabel  = REPORT_SEQ_NAME_DENOVO;
  auxI.label        = TAG_TABLE_SPECTRA_FIELD_SEQ_DENOVO;
  auxI.renderer     = SPECTRUM_RENDERER;
  auxI.validator    = TAG_TABLE_SPECTRA_FIELD_SEQ_DENOVO;
  auxI.id           = specN;

  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PKLBIN,   SPECPLOT_VAL_SPECTRA_PKLBIN,        SPECPLOT_COND_SPECTRA_PKLBIN        ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_SPECTURM, SPECPLOT_VAL_SPECTRA_SPECTURM,      SPECPLOT_COND_SPECTRA_SPECTURM      ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PEPTIDE,  SPECPLOT_VAL_SPECTRA_PEPTIDE_NOVO,  SPECPLOT_COND_SPECTRA_PEPTIDE_NOVO  ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TARGET,   SPECPLOT_VAL_TARGET,                SPECPLOT_COND_TARGET                ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_ENCODING, SPECPLOT_VAL_ENCODING,              SPECPLOT_COND_ENCODING              ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TITLE,    SPECPLOT_VAL_SPECTRA_TITLE,         SPECPLOT_COND_SPECTRA_TITLE         ) );

  auxB.sequences.push(auxI);


  var auxI = new ReportColumnTypeImageOnDemand();
  auxI.dynamic      = true;
  auxI.columnLabel  = REPORT_SEQ_NAME_USER;
  auxI.label        = TAG_TABLE_SPECTRA_FIELD_SEQ_USER;
  auxI.renderer     = SPECTRUM_RENDERER;
  auxI.imgSequence  = TAG_TABLE_SPECTRA_FIELD_SEQ_USER;
  auxI.id           = specU;
  //auxI.validator    = TAG_TABLE_SPECTRA_FIELD_SEQ_USER;

  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PKLBIN,   SPECPLOT_VAL_SPECTRA_PKLBIN,        SPECPLOT_COND_SPECTRA_PKLBIN          ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_SPECTURM, SPECPLOT_VAL_SPECTRA_SPECTURM,      SPECPLOT_COND_SPECTRA_SPECTURM        ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_PEPTIDE,  SPECPLOT_VAL_SPECTRA_PEPTIDE_USER,  SPECPLOT_COND_SPECTRA_PEPTIDE_USER, 0 ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TARGET,   SPECPLOT_VAL_TARGET,                SPECPLOT_COND_TARGET                  ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_ENCODING, SPECPLOT_VAL_ENCODING,              SPECPLOT_COND_ENCODING                ) );
  auxI.params.push( new ReportParamsOption(SPECPLOT_PAR_TITLE,    SPECPLOT_VAL_SPECTRA_TITLE,         SPECPLOT_COND_SPECTRA_TITLE           ) );

  auxB.sequences.push(auxI);


  var auxS = new ReportColumnTypeString();
  auxS.dynamic      = true;
  auxS.isInput      = true;
  auxS.id           = inputID;
  auxB.sequences.push(auxS);

  var auxS = new ReportColumnTypeString();
  auxS.dynamic      = true;
  auxS.isButton     = true;
  auxS.text         = REPORT_BUTTON_UPDATE;
  auxS.onClick      = this.buildUpdateCall(inputID, specU);
  auxB.sequences.push(auxS);

  this.colTypes.push(auxB);


  // colTypes[3] . (CTstring) protein name
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_SPECTRA_4;
  auxS.text         = TAG_TABLE_SPECTRA_FIELD_PROTEIN_NAME;
  this.colTypes.push(auxS);

  // colTypes[3] . (CTstring) mass
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_SPECTRA_5;
  auxS.text         = TAG_TABLE_SPECTRA_FIELD_MASS;
  this.colTypes.push(auxS);

  // colTypes[4] . (CTstring) charge
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_SPECTRA_6;
  auxS.text         = TAG_TABLE_SPECTRA_FIELD_CHARGE;
  this.colTypes.push(auxS);

  // colTypes[5] . (CTstring) B%
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_SPECTRA_7;
  auxS.text         = TAG_TABLE_SPECTRA_FIELD_B_PER;
  this.colTypes.push(auxS);

  // colTypes[6] . (CTstring) Y%
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_SPECTRA_8;
  auxS.text         = TAG_TABLE_SPECTRA_FIELD_Y_PER;
  this.colTypes.push(auxS);

  // colTypes[7] . (CTstring) BY intensity %
  var auxS = new ReportColumnTypeString();
  auxS.columnLabel  = REPORT_HEADER_SPECTRA_9;
  auxS.text         = TAG_TABLE_SPECTRA_FIELD_BY_INTENSITY;
  this.colTypes.push(auxS);
}
////////////////////////////////////////////////////////////////////////////////
// Document
////////////////////////////////////////////////////////////////////////////////
function Report()
{
  // the table vector
  this.tables = new Array();
}
////////////////////////////////////////////////////////////////////////////////
Report.prototype.getData = function()
{
  for(var i = 0 ; i < this.tables.length ; i++)
    this.tables[i].getData();
}
////////////////////////////////////////////////////////////////////////////////
Report.prototype.getId = function()
{
  if(this.tables.length)
    return this.tables[0].getId();
}
////////////////////////////////////////////////////////////////////////////////
Report.prototype.getField = function(tab, row, col)
{
  if(tab < this.tables.length)
    return this.tables[tab].getField(row, col);
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
function ReportProtein()
{
  this.tables = new Array();
}

ReportProtein.prototype = new Report();


ReportProtein.prototype.proteinsPage = function()
{
  // Define contig list view, not filtered
  var c = new tableProtein();
  // define view for this table
  c.defineView();
  // Add the table to table list
  this.tables.push(c);
}

ReportProtein.prototype.proteinPage = function(protein)
{
  // Define contig list view, not filtered
  var c = new tableProtein();
  // define view for this table
  c.defineView2();
  // set fitler field, value
  c.filterField = 0;
  c.filterData  = protein;
  // Add the table to table list
  this.tables.push(c);

  // Define contig list view, not filtered
  var c = new tableContig();
  // define view for this table
  c.defineView();
  // set fitler field, value
  c.filterField = 1;
  c.filterData  = protein;
  // Add the table to table list
  this.tables.push(c);
}
////////////////////////////////////////////////////////////////////////////////
function ReportProteinCoverage()
{
  this.tables = new Array();
}
////////////////////////////////////////////////////////////////////////////////
ReportProteinCoverage.prototype = new Report();

////////////////////////////////////////////////////////////////////////////////
ReportProteinCoverage.prototype.proteinCoveragePage = function(protein)
{
  // Define contig list view, not filtered
  var c = new tableProteinCoverage();
  // define view for this table
  c.defineView();
  // set fitler field, value
  c.filterField = 0;
  c.filterData  = protein;
  // Add the table to table list
  this.tables.push(c);
}
////////////////////////////////////////////////////////////////////////////////
ReportProteinCoverage.prototype.proteinCoverageCSVPage = function(protein)
{
  // Define contig list view, not filtered
  var c = new tableProteinCoverage();
  // define view for this table
  c.defineView2();
  // set fitler field, value
  c.filterField = 0;
  c.filterData  = protein;
  // Add the table to table list
  this.tables.push(c);
}
////////////////////////////////////////////////////////////////////////////////
function ReportContig()
{
  this.tables = new Array();
}

ReportContig.prototype = new Report();


ReportContig.prototype.contigsPage = function()
{
  // Define contig list view, not filtered
  var c = new tableContig();
  // define view for this table
  c.defineView();
  // Add the table to table list
  this.tables.push(c);
}

ReportContig.prototype.contigPage = function(contig)
{
  // Define contig list view, not filtered
  var c = new tableContig();
  // define view for this table
  c.defineView2();
  // set fitler field, value
  c.filterField = 0;
  c.filterData  = contig;
  //disable borders
  c.drawBorders = 0;
  // Add the table to table list
  this.tables.push(c);

  // Define contig list view, not filtered
  var c = new tableCluster();
  // define view for this table
  c.defineView();
  // set fitler field, value
  c.filterField = 1;
  c.filterData  = contig;
  // Add the table to table list
  this.tables.push(c);
}
////////////////////////////////////////////////////////////////////////////////
function ReportCluster()
{
  this.tables = new Array();
}

ReportCluster.prototype = new Report();


ReportCluster.prototype.clusterPage = function(cluster)
{
  // Define contig list view, not filtered
  var c = new tableCluster();
  // define view for this table
  c.defineView2();
  // set fitler field, value
  c.filterField = 0;
  c.filterData  = cluster;
  //disable borders
  c.drawBorders = 0;
  // Add the table to table list
  this.tables.push(c);

  // Define contig list view, not filtered
  var c2 = new tableInputSpectra();
  // define view for this table
  c2.defineView();
  // set fitler field, value
  c2.filterField = 3;
  c2.filterData  = cluster;
  // Add the table to table list
  this.tables.push(c2);
}
////////////////////////////////////////////////////////////////////////////////
function ReportInputSpectra()
{
  this.tables = new Array();
}

ReportInputSpectra.prototype = new Report();


ReportInputSpectra.prototype.inputSpectraPage = function(fileIndex)
{
  // Define contig list view, not filtered
  //var c = new tableInputSpectra();
  // define view for this table
  //c.defineView2();
  // set fitler field, value
  //c.filterField = 0;
  //c.filterData  = cluster;
  // Add the table to table list
  //this.tables.push(c);

  // Define contig list view, not filtered
  var c = new tableInputSpectra();
  // define view for this table
  c.defineView2();
  // set fitler field, value
  c.filterField = 5;
  c.filterData  = fileIndex;
  // Add the table to table list
  this.tables.push(c);
}
////////////////////////////////////////////////////////////////////////////////
// Renderer
////////////////////////////////////////////////////////////////////////////////
function renderer()
{
  this.currentRow = 0;
  this.currentCol = 0;
}
////////////////////////////////////////////////////////////////////////////////
// build the report:
// div      -> where to put the generate page code
// rep      -> the report object
// barType  -> navigation bar
renderer.prototype.buildReport = function(div, rep, start, functionName , barType)
{
  // global row counter for this report
  this.currentRow = 0;

  // create reference for storing return data
  var out = createReference("");

  // at this point, tables are supposed to already have all data
  for(var i = 0 ; i < rep.tables.length ; i++) {
    // build table
    this.buildTable(out, rep.tables[i], start, functionName);
    // clear
    //setReference(out, "");
  };

  // create reference for storing the navigation bar data
  var nav = createReference("");
  // build the navigation bar
  this.buildNavigationBar(nav, rep, barType);
  // horizontal line under

  nav += "<div><table width='100%'><tr></td><td class='ln'></td></tr><tr><td class='VHSep'></td></tr></table></div>";

  if(div == 'NEW') {

    // open in a new window
    var winRef = window.open( "","" )
    winRef.document.writeln("<pre>" + out + "</pre>");
    winRef.document.close();

  } else {

    // set data on page
    var target = document.getElementById(div);
    if(target)
      target.innerHTML = nav + out;
  }

  // call image queue processor
  processQueue();
}
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.addBarLink = function(out, func, img, text)
{
  var sep = "&nbsp;&nbsp;&nbsp;";

  addToReference(out, "<a href='#' onclick='javascript:" + func + "'>");
  addToReference(out, "<img height='24' width='24' src='data:image/png;base64," + img + "' />");
  addToReference(out, "<span style='font-family:Calibri;font-size:140%;color:blue'><u>" + text + "</u></span>");
  addToReference(out, "</a>");
  addToReference(out, sep);
}

renderer.prototype.buildNavigationBar = function(out, rep, barType)
{
  var sep = "&nbsp;&nbsp;&nbsp;";
  // initial page
  this.addBarLink(out, "showMain(true);", iconHome, "Initial page");

  // if no bartype defined, end here
  if(typeof(barType) == 'undefined')
    return;

  // protein --> add protein list
  if(barType == 1 || barType == 2) {
    this.addBarLink(out, "loadProteinsPage();", iconProteinList, "Protein list");
  }

  // contig --> add jump to protein
  if(barType == 2) {
    var aux = rep.getId();
    if(aux.length >= 1) {
      this.addBarLink(out, "loadProteinPage(" + aux[0] + ");", iconProtein, "Protein");
      this.addBarLink(out, "loadProteinDetails(" + aux[0] + ");", iconProteinCoverage, "Protein coverage");
    }

    if(aux.length >= 2) {
      this.addBarLink(out, "loadContigPage(" + aux[1] + ");", iconContig, "Contig");
    }
  }

  if(barType == 3) {
    this.addBarLink(out, "loadContigsPage();", iconContigList, "Contig list");
  }

  if(barType == 4) {

    var aux = rep.getId();
    this.addBarLink(out, "loadContigsPage();", iconContigList, "Contig list");
    this.addBarLink(out, "loadContigPage(" + aux[1] + ");", iconContig, "Contig");
  }
}
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.buildTable = function(out, tab, start, functionName)
{
  if(tab.exception == "ProteinException")
    return this.renderTableExceptionProteinHeader(out, tab);
  if(tab.exception == "ProteinCoverageException")
    return this.renderTableExceptionProteinCoverage(out, tab);
  if(tab.exception == "ProteinCoverageExceptionCSV")
    return this.renderTableExceptionProteinCoverageCSV(out, tab);

  // create reference for storing return data
  var paginationSequenceBefore = createReference("");
  var paginationSequenceAfter  = createReference("");

  // draw table pagination, in case there is need for it
  this.paginate(paginationSequenceBefore, paginationSequenceAfter, tab, functionName, start);

  // add pagination sequence at beggining of table
  addToReference(out, paginationSequenceBefore);


  // table initiator
  var aux;
  if(tab.drawBorders == 1)
    aux = "<table class='result sortable' align='center'>";
  else
    aux = "<table border='0' align='center'>";

  addToReference(out, aux);

  // column types vector for this tables
  var colTypes = tab.colTypes;

  // draw table header
  if(tab.drawBorders == 1)
    this.buildTableHeaderRow(out, colTypes);

  // check for empty tables
  if(tab.theArray.length < 1)
    return;

  if(tab.theArray.length == 1 && tab.theArray[0].length == 1)
    return;

  // set initial elemnt index
  var theStart = 0;
  if(typeof(start) != 'undefined')
    theStart = start;

  // cycle thru table rows (all rows)
  for(var i = theStart ; i < tab.theArray.length && i < theStart + PAGE_LENGTH; i++) {
    // the the data row
    var row = tab.theArray[i];
    // add row
    this.buildTableRow(out, row, colTypes);
    // increase row counter
    this.currentRow++;
  }
  addToReference(out, "</table>");

  // add pagination sequence at beggining of table
  addToReference(out, paginationSequenceAfter);
}
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.buildTableHeaderRow = function(out, colTypes)
{
  addToReference(out, "<tr>");
  for(var j = 0 ; j < colTypes.length ; j++) {
    this.buildTableHeaderCell(out, colTypes[j]);
  }
    addToReference(out, "</tr>");
}
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.buildTableHeaderCell = function(out, colType)
{
  addToReference(out, "<th>");
  //addToReference(out, "<a style='text-decoration: none; font-weight: bold; color: black; ' onmouseover='this.oldstyle=this.style.cssText;this.style.color=\"blue\"' onmouseout='this.style.cssText=this.oldstyle;' class='abc' href='#' title='Click here to Sort!' onclick='javascript:ml_tsort.resortTable(this.parentNode);return false'><span style='color:white'>");

  addToReference(out, "<span style='color: white; '>");
  addToReference(out, colType.columnLabel.replace(/ /g, '&nbsp;'));
  addToReference(out, "</span>");

  //addToReference(out, "</a>");
  addToReference(out, "</th>");
}
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.buildTableRow = function(out, row, colTypes)
{
  this.currentCol = 0;
  addToReference(out, "<tr>");
  for(var j = 0 ; j < colTypes.length ; j++) {
    this.buildTableCell(out, row, colTypes[j]);
    this.currentCol++;
  }
  addToReference(out, "</tr>");
}
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.buildTableCell = function(out, row, colType)
{
  // auxiliary variables
  var cls = "", slink = "";

  // gather needed attributes
  if(colType.link.length)
    slink = this.parseTemplates(colType.link, row, TAG_OPEN, TAG_CLOSE, 0);

  if(colType.cssClass.length)
    cls = " class='" + base.cssClass + "'";

  // cell begin HTML tag with class
  var aux = "<td align='center' valign='middle'" + cls + ">";
  addToReference(out, aux);
  // Link section
  if(colType.link.length) {
    aux = "<a href='" + slink + "'>";
    addToReference(out, aux);
  }
  // process base class cell renderer
  this.buildTableCellSpecific(out, row, colType);
  // Link section
  if(colType.link.length) {
    aux = "</a>";
    addToReference(out, aux);
  }
  // cell end HTML tag
  addToReference(out, "</td>");

  //addToReference(out, "<td>");
  //var content = this.parseTemplates(colType.text, row, TAG_OPEN, TAG_CLOSE, 0);
  //addToReference(out, content);
  //addToReference(out, "</td>");
}
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.buildTableCellSpecific = function(out, row, colType)
{
  // check for valid cells
  if(colType.validator.length) {
    var aux = this.parseTemplates(colType.validator, row, TAG_OPEN, TAG_CLOSE, 0);
    if(aux.length == 0)
      return;
  }

  // test for Image On Demand column type
  if(colType.rtti == REPORT_CELL_TYPE_IOD)
    return this.buildCellImageOnDemand(out, row, colType);

  // test for String column type
  if(colType.rtti == REPORT_CELL_TYPE_STRING)
    return this.buildCellString(out, row, colType);

  // test for String column type
  if(colType.rtti == REPORT_CELL_TYPE_STRING_MULTIPLE)
    return this.buildCellStringMultiple(out, row, colType);

  // test for Box column type
  if(colType.rtti == REPORT_CELL_TYPE_BOX)
    return this.buildCellBox(out, row, colType);

  // something's wrong...
  return -1;
}
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.buildCellImageOnDemand = function(out, row, colType)
{
  // auxiliary variables
  var icon = "", label = "", url = "", tag = "", tag2 = "", tag3 = "", tag4 = "", id = "";

  // set ids for images.
  if(colType.id.length) {
    id = this.parseTemplates(colType.id, row, TAG_OPEN, TAG_CLOSE, 0);
    tag  = IMAGE_ICON_ID_PREFIX       + id;
    tag2 = IMAGE_ICON_CTRL_ID_PREFIX  + id;
    tag3 = IMAGE_LARGE_ID_PREFIX      + id;
    tag4 = IMAGE_LARGE_CTRL_ID_PREFIX + id;
  } else {
    tag  = IMAGE_ICON_ID_PREFIX       + globalImage;
    tag2 = IMAGE_ICON_CTRL_ID_PREFIX  + globalImage++;
    tag3 = IMAGE_LARGE_ID_PREFIX      + globalImage;
    tag4 = IMAGE_LARGE_CTRL_ID_PREFIX + globalImage++;
  }

  // Icon path/image (src)
  if(colType.iconParams.length) {
    // get parsed parameters / URL
    var pars = new Array();
    var pars2 = new Array();
    this.parseParamsVector(pars, pars2, colType.iconParams, row);
    // if there is a renderer, use it to render the image
    if(colType.iconRenderer.length) {
      var seqi = this.parseTemplates(colType.iconSequence, row, TAG_OPEN, TAG_CLOSE, 0);
      this.getImage(colType.iconRenderer, pars, pars2, tag, "src", tag2, seqi);
    }
    // URL
    if(colType.link.length) {
      url = colType.link;
    }
    // the tag
    var aux = "<img id='" + tag + "'";
    if(colType.onClick.length) {
      url = this.parseTemplates(colType.onClick, row, TAG_OPEN, TAG_CLOSE, 0);
      aux += " onclick='" + url + "'";
    }
    aux += " src='data:image/gif;base64," + loadingImage + "' />";
    // the hidden command
    aux += "<input id='" + tag2 + "' type='hidden' text=''>";
    // store
    addToReference(out, aux);
  }


  // URL template to be used to get the image (href)
  if(colType.renderer.length) {
    // get parsed parameters / URL
    var pars = new Array();
    var pars2 = new Array();
    this.parseParamsVector(pars, pars2, colType.params, row);
    var seq = this.parseTemplates(colType.imgSequence, row, TAG_OPEN, TAG_CLOSE, 0);
    this.getImage(colType.renderer, pars, pars2, tag3, "href", tag4, seq);
    // url for IOD
    var image64 = "data:image/png;base64," + loadingImage;
    var aux =  "<a id='" + tag3 + "' href=\"" + image64 + "\" rel=\"lightbox\">";
    // the hidden command
    aux += "<input id='" + tag4 + "' type='hidden' text=''>";
    // store
    addToReference(out, aux);
  }

  // clickable label
  if(colType.label.length) {
    label = this.parseTemplates(colType.label, row, TAG_OPEN, TAG_CLOSE, 0);
    var labelId = "";
    if(colType.id.length)
      labelId = " ID='" + id + "'";
    var aux = "<div" + labelId + ">" + label + "</div>";
    addToReference(out, aux);
  }

  // url for IOD - terminator
  if(colType.renderer.length) {
    var aux = "</a>";
    addToReference(out, aux);
  }

  return 1;
}
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.buildCellString = function(out, row, colType)
{
  // auxiliary variables
  var onclick = "", id = "", text = "", aux= "";

  // gather needed attributes
  if(colType.text.length)
    text = this.parseTemplates(colType.text, row, TAG_OPEN, TAG_CLOSE, 0);

  if(colType.onClick.length)
    onclick = " onclick='" + this.parseTemplates(colType.onClick, row, TAG_OPEN, TAG_CLOSE, 0) + "'";

  if(colType.id.length)
    id = " ID='" + this.parseTemplates(colType.id, row, TAG_OPEN, TAG_CLOSE, 0) + "'";


  // edit box
  if(colType.isInput) {
    aux = "<br /><input class='iud' type='text'" + id + onclick + " />";
  }

  // button
  else if(colType.isButton) {
    aux = "<br /><input type='button' value='" + text + "'" + onclick + id + " />";

  // regular text
  } else {
    aux = "<p" + id + onclick + ">" + text + "</p>";
  }

  addToReference(out, aux);

  return 1;
}
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.buildCellStringMultiple = function(out, row, colType)
{
  // auxiliary variables
  var text = "", link = "";
  var links = new Array();
  var texts = new Array();

  // gather needed attributes
  if(colType.link.size())
    link = this.parseTemplates(colType.link, row, TAG_OPEN, TAG_CLOSE, 0);

  if(colType.text.size())
    text = this.parseTemplates(colType.text, row, TAG_OPEN, TAG_CLOSE, 0);

  links = link.split(TABLE_SEP_L1);
  texts = text.split(TABLE_SEP_L1);

  for(var i = 0 ; i < texts.length ; i++) {

    // the ith link
    if(links.length > i) {
      aux = "<a href='" + colType.linkPrefix + links[i] + colType.linkSuffix + "'>";
      addToReference(out, aux);
    }

    // regular text
    aux = "<p>" + texts[i] + "</p>";
    addToReference(out, aux);

    // the ith link
    if(links.size() > i) {
      aux = "</a>";
      addToReference(out, aux);
    }

    // line break
    if(i < texts.size()-1) {
      aux = "<br>";
      addToReference(out, aux);
    }
  }

  return 1;
}
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.buildCellBox = function(out, row, colType)
{
  // sequences box begin sequence
  var aux = "<table align='center' border='0' width='100%'><tr align='center'><td>";
  addToReference(out, aux);

  // sequences box cycle thru all cells
  var sequences = colType.sequences;

  var begin = true;

  for(var i = 0 ; i < sequences.length ; i++) {

    // validate cell
    if( sequences[i].validator.length) {
      aux = this.parseTemplates(sequences[i].validator, row, TAG_OPEN, TAG_CLOSE, 0);
      if(!aux.length) {
        continue;
      }
    }

    // if there is a column label
    if(sequences[i].columnLabel.length) {
      //new line, if it is not the first one
      if(!begin) {
        aux = "<br>"; //"<br><br>";
        addToReference(out, aux);
      }

      // the column label
      aux = "<b>" + sequences[i].columnLabel + "</b>";
      //new line
      aux += "<br>";
      // add to queue
      addToReference(out, aux);
      // subsequent new lines between items will be inserted
      begin = false;
    }

    var link = "";

    // gather needed attributes
    if(sequences[i].link.length)
      link = this.parseTemplates(sequences[i].link, row, TAG_OPEN, TAG_CLOSE, 0);

    // Link section
    if(sequences[i].link.length) {
      aux = "<a href='" + link + "'>";
      addToReference(out, aux);
    }

    // call base class to render cell
    this.buildTableCellSpecific(out, row, sequences[i]);

    // Link section
    if(sequences[i].link.length) {
      aux = "</a>";
      addToReference(out, aux);
    }
  }

  // sequences box end sequence
  aux = "</td></tr></table>";
  addToReference(out, aux);

  return 1;
}
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.parseTemplates = function(str, row, tagOpen, tagClose, add)
{
  // return string and auxiliary string to store tag contents
  var ret = "", tag;
  // store the position after each tag
  var lastPosition = 0;
  // get the first tag
  var first = str.indexOf(tagOpen, lastPosition);
  var second = str.indexOf(tagClose, first+1);
  // repeat until all the tags are processed
  while(first != -1 && second != -1) {
    // copy in beetwen tags
    ret += str.substr(lastPosition, first - lastPosition);
    // copy the tag contents to an auxiliary string
    tag = str.substr(first+1, second - first - 1);
    // get the tag contents
    ret += this.getTag(tag, row, tagOpen, tagClose, add);
    // update last position to past the tag
    lastPosition = second + 1;
    // get next tag
    first = str.indexOf(tagOpen, lastPosition);
    second = str.indexOf(tagClose, first+1);
  }

  // the remaining of the string
  ret += str.substr(lastPosition);
  // return the string
  return ret;
}
///////////////////////////////////////////////////////////////////////////////
renderer.prototype.getTag = function(tag, row, tagOpen, tagClose, add)
{
  // set default return value
  var ret = "";
  //ret = tagOpen;
  //ret += tag;
  //ret += tagClose;

  // find multi-tag contents
  var aux = tag.split('|');

  for(var i = 0 ; i < aux.length ; i++) {
    // translate the tag contents
    var content = this.translateTag(aux[i], row, add);
    // check for non-empty contents, and stop earching if not empty
    if(content.length)
      return content;
  }

  // return the translated attribute
  return ret;
}
///////////////////////////////////////////////////////////////////////////////
renderer.prototype.translateTag = function(tag, row, add)
{
  // set default return value
  var ret = "";

  // search for a number
  var myRegExp = /[a-z,A-Z]/;

  var matchPos1 = tag.search(myRegExp);

  if(matchPos1 == -1) {
    var aux = parseInt(tag);
    if(aux < row.length)
      ret = row[aux];

    if(add) {
      var aux2 = parseInt(ret);
      aux2 += add;
      ret = String(aux2);
    }
  }

  // seach for row
  if(tag == INTERNAL_ROW)
    ret = String(this.currentRow);
  // search for col
  else if(tag == INTERNAL_COL)
    ret = String(this.currentCol);
  // projectdir tag
  else if(tag == INTERNAL_PROJDIR)
    ret = PROJECT_DIR;

  // return the translated attribute
  return ret;
}
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.paginate = function(out, out2, tab, functionName, start)
{
  if(typeof(start) == 'undefined')
    start = 0;

  // check for the need of pagination
  if(tab.theArray.length <= PAGE_LENGTH )
    return;

  // auxiliary vars used
  var aux, aux2, aux3 = "";

  // initial sequence - table beggining
  aux = "<table border='0' align='center'><tr><td>";
  aux2  = aux + "&nbsp;</td></tr><tr><td>";

  for(var i = 0 ; i < tab.theArray.length ; i += PAGE_LENGTH) {
    var id1 = tab.getSingleId(i);
    var j = Math.min(i + PAGE_LENGTH - 1, tab.theArray.length - 1);
    var id2 = tab.getSingleId(j);

    var aux4 = "[" + id1;
    if(j > i)
      aux4 += "-" + id2;
    aux4 += "]";

    if(start >= i && start < i + PAGE_LENGTH) {
      aux3 += "<span style='color:#FF0000'>";
      aux3 += aux4;
      aux3 += "</span>";

    } else {
      aux3 += "<a href='#' onclick='javascript:" + functionName + i + ");'>";
      aux3 += aux4;
      aux3 += "</a>";
    }
    aux3 += "&nbsp;"
  }

  aux  += aux3 + "</td></tr><tr><td>&nbsp;</td></tr></table>";
  aux2 += aux3 + "</td></tr></table>";

  addToReference(out, aux);
  addToReference(out2, aux2);
}
///////////////////////////////////////////////////////////////////////////////
// generate params vector based on options object
renderer.prototype.parseParamsVector = function(params, paramsReload, options, row)
{
  // go thru all options
  for(var i = 0 ; i < options.length ; i++) {

    // validate option
    if(options[i].validator.length) {
      var aux = this.parseTemplates(options[i].validator, row, TAG_OPEN, TAG_CLOSE, 0);
      if(!aux.length) {
        continue;
      }
    }

    // store params
    params.push(options[i].param);
    if(options[i].option.length) {
      var aux = this.parseTemplates(options[i].option, row, TAG_OPEN, TAG_CLOSE, 0);
      params.push(aux);

      if(options[i].store == 1) {
        paramsReload.push(options[i].param);
        paramsReload.push(aux);
      }
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
// Store an image request in the image queue
renderer.prototype.getImage = function(renderer, params, paramsReload, tag, target, tag2, seq)
{
  // initialize parameters variables
  var allParams = "", tag2params = "";

  // compose parameters into a single string
  for(var i = 0 ; i < params.length ; i++, allParams += "&")
    allParams += params[i];

  for(var i = 0 ; i < paramsReload.length ; i++, tag2params += "&")
    tag2params += paramsReload[i];

  // build element
  var elem = new queueElement();
  elem.renderer   = renderer;
  elem.params     = allParams;
  elem.tag        = tag;
  elem.tag2       = tag2;
  elem.tag2params = tag2params;
  elem.sequence   = seq;
  elem.target   = target;

  // store it
  queue.push(elem);
}
///////////////////////////////////////////////////////////////////////////////
// Table rendering exceptions
///////////////////////////////////////////////////////////////////////////////
renderer.prototype.breakProteinIntoChunks = function(inData)
{
  var count = 0;
  for(var i = 0 ; i < inData.length ; i++) {
    // string container
    var aux = "";
    for(var j = 0 ; j < inData[i].length ; j++) {
      // line break;
      if(count % AAS_PER_LINE == 0) {
        aux += "<br>";
      // spacer
      } else if(count % AAS_GROUP_SIZE == 0) {
        aux += " ";
      }
      // copy element
      aux += inData[i][j];
      count++;
    }
    inData[i] = aux;
  }
  return count;
}
///////////////////////////////////////////////////////////////////////////////
renderer.prototype.colorProteinString = function(inData, out)
{
  var aux = "";

  for(var i = 0 ; i < inData.length ; i++) {
    if(inData[i].length) {
      aux += "<font color='";
      if(i % 2) {
        aux += "0";
      } else {
        aux += "#AAAAAA";
      }
      aux += "'>";
      aux += inData[i];
      aux += "</font>";
    }
  }
  addToReference(out, aux);
}
///////////////////////////////////////////////////////////////////////////////
renderer.prototype.renderTableExceptionProteinHeader = function(out, tab)
{
  var row = tab.theArray[0];

  // get the protein name, at indice 1
  var proteinId    = row[TABLE_PROTEINS_FIELD_ID];
  var proteinName  = row[TABLE_PROTEINS_FIELD_NAME];
  var contigs      = row[TABLE_PROTEINS_FIELD_CONTIGS];
  var spectra      = row[TABLE_PROTEINS_FIELD_SPECTRA];
  var aas          = row[TABLE_PROTEINS_FIELD_AAS];
  var coverage     = row[TABLE_PROTEINS_FIELD_COVERAGE];
  var sequence     = row[TABLE_PROTEINS_FIELD_SEQUENCE];

  // format protein sequence
  var count = 0;
  var coloredProtein = createReference("");
  var breaked = sequence.split(TABLE_SEP_L1);
  var count = this.breakProteinIntoChunks(breaked);
  this.colorProteinString(breaked, coloredProtein);

  // build legend for protein sequence
  var legend = "";
  var current = 1;
  while(current <= count) {
    legend += "<br>";
    legend += parseInt(current);
    current += AAS_PER_LINE;
  }

  // protein header HTML code
  var aux = "";
  aux += "<table class='result' width='100%' style='border-spacing: 0px;' align='center'>";
  aux += "<tr>";
  aux += "<td colspan='0'><h2><i>" + proteinName + "</i></h2>";
  aux += "<hr><b>" + contigs + " contigs, " + spectra + " spectra, " + aas + " amino acids, " + coverage + " coverage" + "</b></td>";
  aux += "<td></td>";
  aux += "</tr>";
  aux += "<tr> ";
  aux += "<td align='right'><tt>" + legend + "</tt></td>";
  aux += "<td><tt>" + coloredProtein + "</tt></td>";
  aux += "</tr>";
  aux += "</table>";

  // link to protein detains page
  aux += "<div align='center'>";
  aux += "<a href='#' onclick='javascript:loadProteinDetails(" + proteinId + ");'>Protein coverage</a>";
  aux += "</div>";
  aux += "<br>";

  addToReference(out, aux);

    // return status OK
  return 1;
}
////////////////////////////////////////////////////////////////////////////////
// Table Protein Coverage
////////////////////////////////////////////////////////////////////////////////
// Generated table, per row:
//
// cells[row][0] -> text          --> Protein ID
// cells[row][1] -> text          --> Protein name
// cells[row][2] -> text          --> Protein length (AAs)
// cells[row][3] -> text list     --> Protein sequence, separated by |
// cells[row][4] -> Contig data   --> CSPS Contigs, separated by |
// cells[row][5] -> Contig data   --> SPS Contigs, separated by |
//      Contig data: : items separated by &
//        0 -> Contig ID
//        1 -> Contig name
//        2 -> Contig start
//        3 -> Contig end
//        4 -> Contig items
//           Contig Item: items separated by @. Contents separated by !
//              0 -> Beginning
//              1 -> Span
//              0 -> Content
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.renderTableExceptionProteinCoverage = function(out, tab)
{
  var row = tab.theArray[0];

  // get the protein coverage data
  var proteinId       = parseInt(row[TABLE_COVERAGE_FIELD_ID]);
  var proteinName     = row[TABLE_COVERAGE_FIELD_NAME];
  var proteinLength   = parseInt(row[TABLE_COVERAGE_FIELD_SEQ_REFERENCE]);
  var proteinSequence = row[TABLE_COVERAGE_FIELD_PROT_SEQUENCE];
  var contigCsps      = row[TABLE_COVERAGE_CSPS_DATA];
  var contigSps       = row[TABLE_COVERAGE_SPS_DATA];

  // global protein length used to send protein sequence back to server
  globalProteinLength = proteinLength;

  // split the protein
  var proteinData = proteinSequence.split(TABLE_SEP_L1);

  // variables to hold contig data
  var contigDataCsps = new Array();
  var contigDataSps  = new Array();

  // get CSPS contig data
  this.buildContigDataStructure(contigDataCsps, contigCsps);
  // get SPS contig data
  this.buildContigDataStructure(contigDataSps, contigSps);

  // protein header HTML code
  var page = "";
  page += "<table class='result' width='100%' style='border-spacing: 0px;' align='center'>";
  page += "<tr>";
  page += "<td colspan='0'><h2><i>" + proteinName + "</i></h2><hr></td>";
  page += "</tr>";
  page += "<tr><td>&nbsp;</td></tr>";
  page += "</table>";

  // CSV protein coverage info
  page += "<table><tr><td><a href='#' onclick='javascript:loadProteinDetailsCSV(" + proteinId + ");'>Protein coverage as Excel-ready format (TXT file)</a></td></tr><tr><td>&nbsp;</td></tr></table>";


  addToReference(out, page);

  // general position indexer
  var i = 0;

  // Keep it under protein size
  while(i < proteinLength) {

    // Build a map key index. This is used to maintain the contig index order when outputing them under the protein sequence
    var spsID  = new Array();
    var cspsID = new Array();
    // get the csps contig indexes
    this.getOrder(contigDataCsps, i, CELLS_PER_LINE, cspsID);
    // get the sps contig indexes
    this.getOrder(contigDataSps, i, CELLS_PER_LINE, spsID);

    // if we are starting a new table, add table beggining
    page = "<table class=\"result2\" width=\"100%\" style=\"background-color: #CCCFFF\">\n";

    // store data so far
    addToReference(out, page);

    // output protein
    this.generateProteinSequence(out, i, proteinData, proteinLength);

    // generate input sequence slots
    this.generateInputSequence(out, i, proteinData, proteinLength);

    // Add CSPS contig information (if exists)
    this.generateOutputContig(i, proteinLength, cspsID, out, CELLS_PER_LINE, contigDataCsps, false);

    // Add SPS contig information (if exists)
    this.generateOutputContig(i, proteinLength, spsID, out, CELLS_PER_LINE, contigDataSps, true);

    // HTML table terminator
    page = "</table><br>\n";
    addToReference(out, page);

    i += CELLS_PER_LINE;
  }

  page += "<table width='100%'><tr><td>&nbsp;</td></tr>";
  page += "<table width='100%'><tr><td>By clicking '" + REPORT_BUTTON_ALIGN + "', the modified protein will be added to the protein database file with the ID and description provided below.</td></tr>";
  page += "<table width='100%'><tr><td>Protein ID: <input id='ProtID' type='text' size='40'></td></tr>";
  page += "<table width='100%'><tr><td>Description: <input id='ProtDESC' type='text' size='130'></td></tr>";
  page += "<table width='100%'><tr><td>&nbsp;</td></tr>";
  page += "<tr><td style='text-align:center;' width='100%'><input type='button' value='" + REPORT_BUTTON_ALIGN + "' onclick='javascript:submitCoverage();' width='300px' /></td></tr></table>";

  //page  = "</div></div>";
  //page += "</body>";
  addToReference(out, page);
}
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.renderTableExceptionProteinCoverageCSV = function(out, tab)
{
  var row = tab.theArray[0];

  // get the protein coverage data
  var proteinId       = parseInt(row[TABLE_COVERAGE_FIELD_ID]);
  var proteinName     = row[TABLE_COVERAGE_FIELD_NAME];
  var proteinLength   = parseInt(row[TABLE_COVERAGE_FIELD_SEQ_REFERENCE]);
  var proteinSequence = row[TABLE_COVERAGE_FIELD_PROT_SEQUENCE];
  var contigCsps      = row[TABLE_COVERAGE_CSPS_DATA];
  var contigSps       = row[TABLE_COVERAGE_SPS_DATA];

  // split the protein
  var proteinData = proteinSequence.split(TABLE_SEP_L1);

  // variables to hold contig data
  var contigDataCsps = new Array();
  var contigDataSps  = new Array();

  // get CSPS contig data
  this.buildContigDataStructure(contigDataCsps, contigCsps);
  // get SPS contig data
  this.buildContigDataStructure(contigDataSps, contigSps);

  // general position indexer
  var i = 0;

  // Keep it under protein size
  while(i < proteinLength) {

    // Build a map key index. This is used to maintain the contig index order when outputing them under the protein sequence
    var spsID  = new Array();
    var cspsID = new Array();
    // get the csps contig indexes
    this.getOrder(contigDataCsps, i, CELLS_PER_LINE, cspsID);
    // get the sps contig indexes
    this.getOrder(contigDataSps, i, CELLS_PER_LINE, spsID);

    // output protein
    this.generateProteinSequenceCSV(out, i, proteinData, proteinLength);

    // Add CSPS contig information (if exists)
    this.generateOutputContigCSV(i, proteinLength, cspsID, out, CELLS_PER_LINE, contigDataCsps, false);

    // Add SPS contig information (if exists)
    this.generateOutputContigCSV(i, proteinLength, spsID, out, CELLS_PER_LINE, contigDataSps, true);

    i += CELLS_PER_LINE;
  }
}
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.generateProteinSequence = function(out, i, proteinData, proteinLength)
{
  var page = "";

  page += "<tr>";
  page += "<td class='rc3' align='center'>";
  page += parseInt(i+1);
  page += "&nbsp;</td>";

  for(var j = i ; ( j < i + CELLS_PER_LINE ) && ( j < proteinLength ) ; j++) {

    page += "<td align='center' id='" + PEP_ELEM_PREFIX + j + "'";
    // if an empty cell, it's a separator column. It should be 1 pixel wide
    if( proteinData[j].length == 0 )
      page += "class='rh2' ";
    else {
      page += "class='rh1' ";
    }
    page += ">";
    // The AA from the protein sequence
    page +=  proteinData[j];
    // header cell terminator
    page += "</td>";
  }
  // Header row terminator
  page += "</tr>\n";

  // store data so far
  addToReference(out, page);
}
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.generateProteinSequenceCSV = function(out, i, proteinData, proteinLength)
{
  var page = "";

  // if we are starting a new table, add table beggining
  page += '\n';
  page += parseInt(i+1);
  page += CSV_SEP;

  for(var j = i ; ( j < i + CELLS_PER_LINE ) && ( j < proteinLength ) ; j++)
    page += CSV_SEP + proteinData[j] + CSV_SEP;

  // Header row terminator
  page += "\n";

  // store data so far
  addToReference(out, page);
}
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.generateInputSequence = function(out, i, proteinData, proteinLength)
{
  var page = "";

  page += "<tr>";
  page += "<td class='rc3' align='right'>";
  //page += parseInt(i+1);
  page += "<input type='checkbox' onchange='javascript:fillCoverageRow(this, " + i + ");' id='ck" + i + "' />";
  //page += "<input type='checkbox' id='ck" + i + "' />";
  page += "&nbsp;</td>";

  for(var j = i ; ( j < i + CELLS_PER_LINE ) && ( j < proteinLength ) ; j++) {

    page += "<td ";
    // if an empty cell, it's a separator column. It should be 1 pixel wide
    if( proteinData[j].length == 0 )
      page += "class='rh2' ";
    else {
      page += "class='rh1' ";
    }
    page += ">";
    // The AA from the protein sequence
    page += "<input class='iuc' type='text' id='" + INP_ELEM_PREFIX + j + "' />";
    // header cell terminator
    page += "</td>";
  }
  // Header row terminator
  page += "</tr>\n";

  // store data so far
  addToReference(out, page);
}
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.generateOutputContig = function(i, proteinSize, vectorID, out, cellPerLine, contig, link)
{
  var page = "";
  // Add contig information (if exists)
  for(var j = 0 ; j < vectorID.length ; j++)  {
    // get the contig sequence info
    var contigSequence = contig[vectorID[j]];
    // Check if sequence in the range we are outputing now
    if( (contigSequence[2] < i + cellPerLine) &&
        (contigSequence[3] > i              )    ) {

      // Write the contig id and link
      if(link) {
        page += "<tr><th align='right'><a href='#' onclick='javascript:loadContigPage(";
        page += this.getIntFromSeqName(contigSequence[1]); // 1 = name
        page += ");' style='color: white'>";
        page += contigSequence[1];
        page += "</a></th>";

        //page += "<tr><th align=\"right\">";
        //page += contigSequence.name;
        //page += "</th>";
      } else {
        page += "<tr>";
        page += "<td class='rc4'>";
        page += "CSPS ";
        //page += contigSequence.name;
        page += parseInt(contigSequence[0] + 1);
        page += "&nbsp;";
        page += "</TD>";
      }

      // find first cell to output
      var l = 0;
      while( (l < contigSequence[4].length) &&
             (i > contigSequence[4][l][0] + contigSequence[4][l][1]) )
        l++;

      // cycle thru
      for(var k = i ; (k < i + cellPerLine) && (k < proteinSize) ; k++) {

        // if start position is lower than current position, output an empty cell
        if(k < contigSequence[2])
          page += "<td class='rc2' />\n";
        else if( (l >= contigSequence[4].length) )
          page += "<td class='rc2' />\n";
        // otherwise, the content
        else {

          page += "<td ";
          // page += " class=\"rc1\" style=\"background-color: transparent; border: solid 0 #060; border-left-width:1px;border-right-width:1px;\" ";

          var border = 0;

          var outputString = contigSequence[4][l][2];

          // Calc colspan nr of occupied cells
          var colspan = contigSequence[4][l][1] + 1;
          // careful with split cells at the beggining or end -- end
          if(contigSequence[4][l][0] + contigSequence[4][l][1] >= i + cellPerLine) {
            colspan -= contigSequence[4][l][0] + contigSequence[4][l][1] - i - cellPerLine + 1;
            border += 2;
            if(colspan < (contigSequence[4][l][1] + 1) / 2 )
              outputString = "";
          }
          // beggining
          if(contigSequence[4][l][0] < i) {
            colspan -= i - contigSequence[4][l][0];
            border++;
            if(colspan <= (contigSequence[4][l][1] + 1) / 2 )
              outputString = "";
          }

          page += " class='rca";
          page += parseInt(border);
          page += "' ";

          if(colspan > 1) {
            page += " colspan='";
            page += parseInt(colspan);
            page += "'";
            k += colspan-1;
          }
          page += '>';
          page +=  outputString;
          page += "</td>";
          l++;
        }
      }
      page += "</tr>\n";
    }
  }
  addToReference(out, page);
}
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.generateOutputContigCSV = function(i, proteinSize, vectorID, out, cellPerLine, contig, link)
{
  var page = "";
  // Add contig information (if exists)
  for(var j = 0 ; j < vectorID.length ; j++)  {
    // get the contig sequence info
    var contigSequence = contig[vectorID[j]];
    // Check if sequence in the range we are outputing now
    if( (contigSequence[2] < i + cellPerLine) &&
        (contigSequence[3] > i              )    ) {

      // Write the contig id and link
      if(link) {
        page += contigSequence[1];
      } else {
        page += parseInt(contigSequence[0] + 1);
      }
      // separator
      page += CSV_SEP;

      // find first cell to output
      var l = 0;
      while( (l < contigSequence[4].length) &&
             (i > contigSequence[4][l][0] + contigSequence[4][l][1]) )
        l++;

      // cycle thru
      for(var k = i ; (k < i + cellPerLine) && (k < proteinSize) ; k++) {

        // if start position is lower than current position, output an empty cell
        if(k < contigSequence[2]) {
          page += CSV_SEP;
          page += CSV_SEP;
        } else if( (l >= contigSequence[4].length) ) {
          page += CSV_SEP;
          page += CSV_SEP;
        // otherwise, the content
        } else {

          var border = 0;

          var outputString = contigSequence[4][l][2];

          // Calc colspan nr of occupied cells
          var colspan = contigSequence[4][l][1] + 1;
          // careful with split cells at the beggining or end -- end
          if(contigSequence[4][l][0] + contigSequence[4][l][1] >= i + cellPerLine) {
            colspan -= contigSequence[4][l][0] + contigSequence[4][l][1] - i - cellPerLine + 1;
            border += 2;
            if(colspan < (contigSequence[4][l][1] + 1) / 2 )
              outputString = "";
          }
          // beggining
          if(contigSequence[4][l][0] < i) {
            colspan -= i - contigSequence[4][l][0];
            border++;
            if(colspan <= (contigSequence[4][l][1] + 1) / 2 )
              outputString = "";
          }

          if(!(border & 0x01))
            page += '|';

          if(colspan == 1) {
            // cell content
            page += CSV_SEP;
            page += outputString;
          } else {

            if((outputString.length > 0) && (outputString[0] == '(')) {
              // outputing (xx, yy), and empty cells following
              page += CSV_SEP;
              page += outputString;
              // until the end of cell
              while(--colspan && (k < i+cellPerLine)) {
                page += CSV_SEP;
                page += CSV_SEP;
                k++;
              }

            } else {
              // outputing A.B.C.D
              var cells = 0;
              while( (cells < colspan-1) && (k < i+cellPerLine)) {
                page += CSV_SEP;
                page += outputString[cells];
                page += CSV_SEP;
                page += '.';
                k++;
                cells++;
              }
              page += CSV_SEP;
              page += outputString[cells];

            }

          }
          // cell end
          page += CSV_SEP;
          // right border content (separator)
          if( (!(border & 0x02)) && (true) && (true) ) ;
//            page += '|';

          l++;
        }
      }
      page += '\n';
    }
  }
  addToReference(out, page);
}
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.getOrder = function(contig, i, size, order)
{
  // temporary vector to hold contig start position per contig
  var aux = new Array();
  // cycle thru all contigs
  for(var j = 0 ; j < contig.length ; j++ ) {
    // Check for empty contigs
    if(contig[j][4].length == 0) continue;
    // Check if sequence in the range we are outputing now
    if( (contig[j][2] < i + size) &&
        (contig[j][3] > i              )    ) {
      // create ordering cell with contig info
      var orderingCell = new Array();
      orderingCell[0] = j; //contig[j][0];
      orderingCell[1] = contig[j][2];
      orderingCell[2] = contig[j][3];
      aux.push(orderingCell);
    }
  }

  // Order
  aux.sort(orderSortingFunction);
  // set order in order vector
  for(var k = 0 ; k < aux.size() ; k++)
    order.push(aux[k][0]);
}
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.getIntFromSeqName = function(seq)
{
  var aux = seq.split(":");
  if(aux.length > 1)
    return aux[1];
  return "";
}
////////////////////////////////////////////////////////////////////////////////
renderer.prototype.buildContigDataStructure = function(contigData, contig)
{
  // split by contig
  var aux = contig.split(TABLE_SEP_L1);
  // cycle thru all contigs
  for(var i = 0 ; i < aux.length ; i++) {
    // array to hold data for one contig
    contigData[i] = new Array();
    // get contig header info + elements
    var aux2 = aux[i].split(TABLE_SEP_L2);
    // put header data
    // ID
    contigData[i][0] = parseInt(aux2[0]);
    // Name
    contigData[i][1] = aux2[1];
    // start
    contigData[i][2] = parseInt(aux2[2]);
    // end
    contigData[i][3] = parseInt(aux2[3]);
    // array for contig elements
    contigData[i][4] = new Array();
    // get elements
    if(typeof(aux2[4]) != 'undefined') {
      var elems = aux2[4].split(TABLE_SEP_L3);
      // cycle thru elems
      for(var j = 0 ; j < elems.length ; j++) {
        // elemet storage space
        contigData[i][4][j] = new Array();
        // split elements
        var elem = elems[j].split(TABLE_SEP_L4);
        // copy start position
        contigData[i][4][j][0] = parseInt(elem[0]);
        // copy colspan
        contigData[i][4][j][1] = parseInt(elem[1]);
        // copy element data
        contigData[i][4][j][2] = elem[2];
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
