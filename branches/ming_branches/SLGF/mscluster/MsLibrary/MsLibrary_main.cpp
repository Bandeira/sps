/*************************************************************************//**
@file MsCluster_main.cpp 

\brief This is the main file for MsCluster it manages the command line interface and 
calls the high-level functions that perfrom actions like dat creation, 
clustering, benchmarks etc.

*****************************************************************************/


#include "../Common/auxfun.h"
#include "../PepNovo/AllScoreModels.h"
#include "../MsCluster/MsParameterStruct.h"
#include "SpectrumLibrary.h"


const string build_name = "20090427"; 


/*!
Prints the help message for the user, specifing all command line flags.
*/
void print_help(const char *message) 
{
	if (message)
		printf("***************************************************************************\n\n%s\n",message);

	printf("\nMsLibrary v1.00 - Spectrum Library Identification.\n");
	printf("Release %s.\nAll rights reserved to the Regents of the University of California.\n\n",build_name.c_str());

	
	printf("Command line Arguments:\n");
	printf("-----------------------\n\n");
	printf("--model <model name> \n\n");
	printf("--load-library   <path> - path to dat list of library files.\n");
	printf("--window     <X> - compare with library spectra +-X Da from the spectrum's precursor mass.\n");
	printf("--out-dir <path> - where files get written (default ./out)\n");
	printf("--create-library  <name>  - name of file prefixes for library.\n");
	printf("--create-from-dat <list>  - list of dat files to be converted (otherwise use --library-list)\n");
	printf("--only-annotated		 - use only annotated spectra to create the library.\n");
	printf("--library-list  <path> - path to text file listing input library files (or a file path)\n");
	printf("--library-test  <path> - path to text file listing spectra files (or a file path)\n");
	printf("--min-mz <X> - only consider spectra with m/z >= X (default X=0)\n");
	printf("--max-mz <X> - only consider spectra with m/z <  X (default X=inf)\n");
	printf("--batch-idx <X> - adds the batch index to the results file names (for split jobs)\n");
	printf("--generation <X> - generation idx\n");
	printf("\nFor testing from mgfs and learning pvalues:\n");
	printf("-------------------------------------------\n");
	printf("--benchmark-library - computes ROC curves for different number of compared peptides (must have --library-list and --library-test too).\n");
	printf("--learn-pvalues - learns the CDFs for the random match p-values (must have --library-list and --library-test also).\n");
	
	printf("\n\nPlease send comments and bug reports to Ari Frank (arf@cs.ucsd.edu).\n\n");
#ifdef WIN32
	system("Pause");
#endif
	exit(1);
}


/*!
\brief The main function entry point.
*/
int main(int argc, char** argv)
{
	
	seedRandom(112233);		// use a fixed seed for "random" so behavior is the same in different runs

	if (argc <= 1)
		print_help(NULL);

	MsParameterStruct params;
	params.copyCommandLine(argc, argv);
	params.batchIdx = -1;

	// This portion of the code parses the command line arguments
	// There is probably a more elegant way to do this
	for (int i=1; i<argc; i++)
	{

		if (! strcmp(argv[i],"--library-list"))
		{
			if (++i == argc)
				print_help("missing library-list!");
			params.spectraListToLoad =argv[i];
			continue;
		}

		if (! strcmp(argv[i],"--library-test"))
		{
			if (++i == argc)
				print_help("missing library-list!");
			params.spectraListToTest =argv[i];
			continue;
		}

		if (! strcmp(argv[i],"--model"))
		{
			if (++i == argc)
				print_help("missing model name!");
			params.modelName=argv[i];
			continue;
		}

		if (! strcmp(argv[i],"--generation"))
		{
			if (++i == argc)
				print_help("missing generation idx!");
			params.generationIdx = atoi(argv[i]);
			continue;
		}

		if (! strcmp(argv[i],"--learn-pvalues"))
		{
			params.gotLearnLibraryPValues = true;
			continue;
		}

		if (! strcmp(argv[i],"--model-dir"))
		{
			if (++i == argc)
				print_help("missing model dir!");
			params.modelDir=argv[i];
			continue;
		}


		if (! strcmp(argv[i],"--out-dir"))
		{
			if (++i == argc)
				print_help("missing output directory!");
			params.outDir=argv[i];
			continue;
		}

		if (! strcmp(argv[i],"--tmp-dir"))
		{
			if (++i == argc)
				print_help("missing dat directory!");
			params.tmpDir=argv[i];
			continue;
		}

		if (! strcmp(argv[i],"--output-name"))
		{
			if (++i == argc)
				print_help("missing name to assign to clusters and files!");
			params.outputName=argv[i];
			continue;
		}

		if (! strcmp(argv[i],"--PTMs"))
		{
			if (++i == argc)
				print_help("missing PTM string!");
			params.ptmString=argv[i];
			continue;
		}

		if (! strcmp(argv[i],"--batch-idx"))
		{
			if (++i == argc)
				print_help("missing sqs threshold!");
			
			params.batchIdx = atoi(argv[i]);
			continue;
		}

		if (! strcmp(argv[i],"--sqs"))
		{
			if (++i == argc)
				print_help("missing sqs threshold!");
			
			params.sqsThreshold = MIN_FLOAT;
			params.sqsThreshold = atof(argv[i]);
			if (params.sqsThreshold<0.0 || params.sqsThreshold>=1.0)
				print_help("--sqs <X>, the treshold should be in the range 0<X<1");
			continue;
		}

		if (! strcmp(argv[i],"--window"))
		{
			if (++i == argc)
				print_help("missing window size!");
			
			params.mzWindow = MIN_FLOAT;
			params.mzWindow = atof(argv[i]);
			if (params.mzWindow<=0.0)
				print_help("--window <X>, select X>0");
			continue;
		}

		if (! strcmp(argv[i],"--min-mz"))
		{
			if (++i == argc)
				print_help("missing min-mz!");
			
			params.minMz = MIN_FLOAT;
			params.minMz = atof(argv[i]);
			if (params.minMz<0.0)
				print_help("--min-mz <X>, select X>=0");
			continue;
		}

		if (! strcmp(argv[i],"--max-mz"))
		{
			if (++i == argc)
				print_help("missing max-mz!");
			
			params.maxMz = MIN_FLOAT;
			params.maxMz = atof(argv[i]);
			if (params.maxMz<=0.0)
				print_help("--max-mz <X>, select X>0");
			continue;
		}

		if (! strcmp(argv[i],"--start-file-index"))
		{
			if (++i == argc)
				print_help("missing start file index!");
			
			params.startFileIdx = 0;
			params.startFileIdx = atoi(argv[i]);
			continue;
		}

		if (! strcmp(argv[i],"--verbose-level"))
		{
			if (++i == argc)
				print_help("missing verbose level!");
			
			params.verboseLevel = atoi(argv[i]);
			continue;
		}

	
		if (! strcmp(argv[i],"--filter-only"))
		{
			params.gotFilterOnly=true;
			continue;
		}

		if (! strcmp(argv[i],"--second-pass"))
		{
			params.gotSecondPass=true;
			continue;
		}


		if (! strcmp(argv[i],"--benchmark-library"))
		{
			params.gotBenchmarkLibrary = true;
			continue;
		}

		if (! strcmp(argv[i],"--library-stats"))
		{
			params.gotLibraryStats = true;
			continue;
		}

		if (! strcmp(argv[i],"--create-library"))
		{
			if (++i == argc)
				print_help("missing name to assign to clusters and files!");
			params.outputName=argv[i];
			params.gotCreateLibrary = true;
			continue;
		}

		if (! strcmp(argv[i],"--create-from-dat"))
		{
			if (++i == argc)
				print_help("missing list file");
			params.datList = argv[i];
			params.gotCreateLibraryFromDat = true;
			params.gotKeepDat = true;
			continue;
		}

		if (! strcmp(argv[i],"--only-annotated"))
		{
			params.gotMakeLibraryWithPeptidesOnly = true;
			continue;
		}

		if (! strcmp(argv[i],"--load-library"))
		{
			if (++i == argc)
				print_help("missing path to library dat files!");
			params.datList = argv[i];
			continue;
		}

		if (! strcmp(argv[i],"--help") || ! strcmp(argv[i],"-h"))
			print_help(NULL);
	

		printf("**********************************************************\n");
		printf("\nError: Unkown command line option: %s\n\n",argv[i]);
		printf("**********************************************************\n");
		exit(1); 
	}

	// Always need a model because it creates the config (that is always needed)
	if (params.modelName.length() == 0)
		params.modelName = "LTQ_TRYP";

	if (params.gotCreateLibrary && params.spectraListToLoad.length()== 0 && ! params.gotCreateLibraryFromDat)
		error("when creating library must supply --library-list!");


	AllScoreModels model;
	Config* config = model.get_config();
	
	if (params.modelDir.length()>0)
		config->set_resource_dir(params.modelDir);

	model.read_model(params.modelName.c_str());
	if (params.ptmString.length()>0)
	{
		config->apply_selected_PTMs(params.ptmString.c_str());
	}
	else
		config->apply_selected_PTMs("C+57");

	

	// This is a library job...
	assert( params.maxMz >= params.minMz );

	if (params.verboseLevel>0)
		cout << "MsLibrary v1.0 " << endl << "Started job on: " << getTimeString() << endl << endl;
	
	SpectrumLibrary library;
	if (params.gotLearnLibraryPValues)
	{
		library.learnCdfsForPValues(params.spectraListToLoad, params.spectraListToTest, config);
	}
	else if (params.gotLibraryStats)
	{
		library.libraryStats(&params, config);
	}
	else if (params.gotBenchmarkLibrary)
	{
		library.benchmarkLibrary(config, &params);
	}
	else if (params.gotCreateLibrary)
	{
		library.createLibrary(&model, &params);
	}
	else if (params.datList.length()>0)
	{
		if (params.spectraListToTest.length()==0)
			error("when identifying with a library, must supply --library-test!");

		createDirIfDoesNotExist(params.outDir.c_str());
		cout << "Reading library: " << params.datList << endl;
		library.readLibrary(&params, config);
		cout << "Identifying spectra: " << params.spectraListToTest << endl;
		library.identifyUsingLibrary(&params, config);
	}

	if (params.verboseLevel>0)
		cout << endl << "Ended job on: " << getTimeString() << endl;
	return 0;
}





