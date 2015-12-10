#include "SpectrumLibrary.h"
#include "../MsCluster/MsClusterAlgorithm.h"
#include "../MsCluster/MsParameterStruct.h"
#include "../MsCluster/DatFileWriter.h"
#include "../MsCluster/DatFileManager.h"


/*! \fn convertDatFileToLibraryEntries
	\brief reads the content of a datFile and converts all the spectra to libraryEntries.

	@return number of libraryEntries written.
*/
size_t convertDatFileToLibraryEntries(const Config* config, DatFile& datFile, LibraryEntry* libraryEntries, bool indOnlyWithPepties)
{
	datFile.openForReading();

	size_t numSpectraWritten = 0;
	size_t numBytes = 0;
	float  sqs;
	char*  specStart = 0;
	long   peaksFilePosition = 0;
	char*  peaksBufferPosition = 0;
	
	Cluster cluster; // this peak list is used for parsing only
	while( datFile.getNextDatSpectrum(specStart, numBytes, &sqs, &peaksFilePosition, &peaksBufferPosition) )
	{
		SingleSpectrumHeader& ssh = libraryEntries[numSpectraWritten].header;
		ssh.setFileType(IFT_DAT);
		ssh.scanSpectrumHeaderFromBuffer(specStart, config);

		if (indOnlyWithPepties && ssh.getPeptideStr().length() < 2)
			continue;

		ssh.setIndexInFile(numSpectraWritten);
		ssh.setPositionInFile(peaksFilePosition);

		cluster.createNewCluster(numSpectraWritten, &ssh, reinterpret_cast<Peak*>(peaksBufferPosition),
								 ssh.getOriginalNumPeaks());

		libraryEntries[numSpectraWritten].distancePeaks = *(cluster.getDistancePeaks());

		numSpectraWritten++;
	}

	datFile.closeAfterReading();
	return numSpectraWritten;
}


float SpectrumLibrary::computeMatchPValue(float similarity, int numPeptidesCompared) const
{
	const double cdf = cdf_[ computeCdfBin(similarity) ];
	return (static_cast<float>(1.0 - pow(cdf, numPeptidesCompared)));
}


size_t SpectrumLibrary::createLibrary(AllScoreModels* model, const MsParameterStruct* params)
{
	const Config* config = model->get_config();

	// first pass, write all spectra into files with large increment
	DatFileWriter datWriter(model);
	
	string tmpDatList;
	
	if (params->gotCreateLibraryFromDat)
	{
		tmpDatList = params->datList;
	}
	else if (! params->gotSecondPass)
	{
		tmpDatList = datWriter.convertDataToDatFirstPass(params);
	}
	else
		tmpDatList = params->spectraListToLoad;
	
	// second pass, load spectra of certain m/z chunk into memory
	DatFileManager datManager;
	datManager.init(tmpDatList, config);

	libraryEntries_.resize(datManager.getMaxNumSpectraInBatch());

	datWriter.init(1000.0);
	datWriter.setDatDir(params->outDir);
	datWriter.setDatName(params->outputName);
	datWriter.setMaxDatFileSize(1<<27);
	PeakList tmpPl;
	Peak* tmpPeaks = new Peak[MAX_NUM_PEAKS_FOR_DISTANCE];
	tmpPl.setPeaksPtr(tmpPeaks);

	size_t numEntriesWritten=0;
	DatBatch nextBatch;

	while ( datManager.getNextBatch(nextBatch) )
	{
		size_t numEntriesInBatch = 0;
		for (size_t datIdx = nextBatch.startIdx; datIdx<= nextBatch.endIdx; datIdx++)
		{
			DatFile& datFile = datManager.getDatFile(datIdx);
			numEntriesInBatch+=convertDatFileToLibraryEntries(config, datFile, &libraryEntries_[0], params->gotMakeLibraryWithPeptidesOnly);
		}

		if (params->verboseLevel>0)
			cout << "Converting batch: " << nextBatch.startIdx << "-" << nextBatch.endIdx << ", "  
			     << numEntriesInBatch << " spectra" << endl;
		

		// sort entries and write them to a library datFile
		sort(libraryEntries_.begin(), libraryEntries_.begin()+numEntriesInBatch);

		for (size_t i=0; i<numEntriesInBatch; i++)
		{
			if (params->datasetIdx>=0)
				libraryEntries_[i].header.setDatasetIndex(params->datasetIdx);
			if (libraryEntries_[i].header.getMOverZ()<=0.0)
				continue;

			tmpPl.setHeader(&libraryEntries_[i].header);
			assert( libraryEntries_[i].header.getMOverZ()>0.0 );
			assert( i == 0 || libraryEntries_[i].header.getMOverZ() >= libraryEntries_[i-1].header.getMOverZ());

			// create peaks
			const DistancePeakList& dpl = libraryEntries_[i].distancePeaks;
			for (int j=0; j<dpl.numPeaks; j++)
			{
				tmpPeaks[j].mass      = dpl.peaks[j].mass;
				tmpPeaks[j].intensity = dpl.peaks[j].intensity;
			}
			tmpPl.setNumPeaks(dpl.numPeaks);
			if (tmpPl.getNumPeaks()<5)
				continue;

			datWriter.addPeakListToDat(tmpPl, true);
		}
		numEntriesWritten += numEntriesInBatch;
	}

	cout << endl << "Wrote " << numEntriesWritten << " spectra to library." << endl;
	datWriter.closeAllOpenDats();
	datWriter.writeDatPaths(false);
	
	// remove first pass files if necessary
	if (! params->gotKeepDat)
		removeFilesInList(tmpDatList.c_str());

	return 0;
}


void SpectrumLibrary::readLibrary(const MsParameterStruct* params, const Config* config)
{
	// library is assumed to be in --dat-list
	assert(params->datList.length()>0);

	vector<string> datPaths;
	readListOfPaths(params->datList.c_str(), datPaths);

	const mass_t minMzToLoad = params->minMz - params->mzWindow;
	const mass_t maxMzToLoad = params->maxMz + params->mzWindow;
	vector<string> pathsToLoad;
	size_t maxLibrarySize =0;
	for (size_t i=0; i<datPaths.size(); i++)
	{
		DatFile dat;
		dat.peekAtStats(datPaths[i].c_str());
		if (dat.getMinMOverZ() < maxMzToLoad && dat.getMaxMOverZ() > minMzToLoad)
		{
			maxLibrarySize += dat.getNumSpectra();
			pathsToLoad.push_back(datPaths[i]);
		}
	}

	cout << "Found " << pathsToLoad.size() << " paths to load." << endl;
	if (maxLibrarySize>0)
	{
		libraryEntries_.resize(maxLibrarySize);
		size_t entryIdx=0;
		for (size_t i=0; i<pathsToLoad.size(); i++)
		{
			DatFile datFile;
			datFile.openForReading(pathsToLoad[i].c_str());

			size_t numBytes = 0;
			float  sqs;
			char* specStart = 0;
			long peaksFilePosition = 0;
			char* peaksBufferPosition = 0;
			size_t numSpectraRead = 0;
			mass_t previousMz = MIN_FLOAT;
			while( datFile.getNextDatSpectrum(specStart, numBytes, &sqs, &peaksFilePosition, &peaksBufferPosition) )
			{
				SingleSpectrumHeader& ssh = libraryEntries_[entryIdx].header;
				ssh.setFileType(IFT_DAT);
				ssh.scanSpectrumHeaderFromBuffer(specStart, config);
				ssh.setIndexInFile(numSpectraRead++);
				ssh.setPositionInFile(peaksFilePosition);

				if (ssh.getMOverZ() < minMzToLoad || ssh.getMOverZ() > maxMzToLoad)
					continue;

				if (ssh.getMOverZ() < previousMz)
					continue; // HACK

				assert(ssh.getMOverZ() >= previousMz);
				previousMz = ssh.getMOverZ();

				const size_t numPeaks = ssh.getOriginalNumPeaks();
				assert(numPeaks <= MAX_NUM_PEAKS_FOR_DISTANCE);
				DistancePeakList& dpl = libraryEntries_[entryIdx].distancePeaks;
				dpl.numPeaks = numPeaks;

				const Peak* peaks = reinterpret_cast<const Peak*>(peaksBufferPosition);
				dpl.sumSqrAdjustedIntensity = 0.0;
				for (size_t j=0; j<numPeaks; j++)
				{
					dpl.peaks[j].mass = peaks[j].mass;
					dpl.peaks[j].intensity = peaks[j].intensity;
#ifndef SIMPLE_DISTANCE
					dpl.sumSqrAdjustedIntensity += (peaks[j].intensity*peaks[j].intensity);
#endif
				}
#ifdef SIMPLE_DISTANCE
				dpl.sumSqrAdjustedIntensity = static_cast<float>(numPeaks);
#endif
				entryIdx++;
			}

			datFile.closeAfterReading();
		}

		const size_t numSpectraLoaded = entryIdx;
		libraryEntries_.resize(numSpectraLoaded);

		
		// check if needs to be sorted
		size_t i;
		for (i=1; i<numSpectraLoaded; i++)
			if (libraryEntries_[i].header.getMOverZ() < libraryEntries_[i-1].header.getMOverZ())
				break;

		if (i<numSpectraLoaded)
		{
			cout << "Sorting entries... " << i << " : " << libraryEntries_[i].header.getMOverZ() 
				<< "    " << i-1 << " : " << libraryEntries_[i-1].header.getMOverZ() << endl;
			// for some reason the sort messes up the library entries

			sort(libraryEntries_.begin(), libraryEntries_.end());
		}

		if (params->verboseLevel > 0)
			cout << "Read " << pathsToLoad.size() << " library dat files with " << libraryEntries_.size() << " spectra." << endl;
	}
	else
	{
		cout << "No spectra loaded from " << params->datList << endl;
	}


	// load the cdf
	loadCdfs(config);
}


void SpectrumLibrary::libraryStats(const MsParameterStruct* params, const Config* config)
{
	assert(params->datList.length()>0);

	vector<string> datPaths;
	readListOfPaths(params->datList.c_str(), datPaths);

	const mass_t minMzToLoad = params->minMz - params->mzWindow;
	const mass_t maxMzToLoad = params->maxMz + params->mzWindow;
	vector<string> pathsToLoad;
	size_t maxLibrarySize =0;
	for (size_t i=0; i<datPaths.size(); i++)
	{
		DatFile dat;
		dat.peekAtStats(datPaths[i].c_str());
		if (dat.getMinMOverZ() < maxMzToLoad && dat.getMaxMOverZ() > minMzToLoad)
		{
			maxLibrarySize += dat.getNumSpectra();
			pathsToLoad.push_back(datPaths[i]);
		}
	}

	map<string, size_t> peptideCounts, peptideChargeCounts;
	size_t totalSpectraCount   =0;
	size_t totalAnnotatedCount =0;

	if (maxLibrarySize>0)
	{
		size_t entryIdx=0;
		for (size_t i=0; i<pathsToLoad.size(); i++)
		{
			DatFile datFile;
			datFile.openForReading(pathsToLoad[i].c_str());

			size_t numBytes = 0;
			float  sqs;
			char* specStart = 0;
			long peaksFilePosition = 0;
			char* peaksBufferPosition = 0;
			size_t numSpectraRead = 0;
			mass_t previousMz = MIN_FLOAT;
			while( datFile.getNextDatSpectrum(specStart, numBytes, &sqs, &peaksFilePosition, &peaksBufferPosition) )
			{
				SingleSpectrumHeader ssh;
				ssh.setFileType(IFT_DAT);
				ssh.scanSpectrumHeaderFromBuffer(specStart, config);

				if (ssh.getMOverZ() < minMzToLoad && ssh.getMOverZ() > maxMzToLoad)
					continue;

				assert(ssh.getMOverZ() >= previousMz);
				previousMz = ssh.getMOverZ();

				const size_t numPeaks = ssh.getOriginalNumPeaks();
				assert(numPeaks <= MAX_NUM_PEAKS_FOR_DISTANCE);

				totalSpectraCount++;
				if (ssh.getPeptideStr().length()>0)
				{
					assert(ssh.getCharge()>0);
					totalAnnotatedCount++;
					peptideCounts[ssh.getPeptideStr()]++;
					ostringstream oss;
					oss << ssh.getPeptideStr() << ssh.getCharge();
					peptideChargeCounts[oss.str()]++;
				}
				
			}

			datFile.closeAfterReading();
		}

		cout << "Read " << pathsToLoad.size() << " library dat files with " << totalAnnotatedCount << "/" 
			 << totalSpectraCount << " annotated spectra." << endl;
		cout << "Unique peptide/charge: " << peptideChargeCounts.size() << endl;
		cout << "Unique peptide       : " << peptideCounts.size() << endl;
	}
	else
	{
		cout << "No spectra loaded from " << params->datList << endl;
	}
}




struct LibrarySearchResult  {
	LibrarySearchResult() : testEntry(0), matchEntry(0), similarity(0.0), pvalue(1.0) {}
	bool operator< (const LibrarySearchResult& rhs) const
	{
		assert(testEntry && rhs.testEntry);
		return (testEntry->position < rhs.testEntry->position);
	}

	const LibraryEntry* testEntry;
	const LibraryEntry* matchEntry;
	float similarity;
	float pvalue;
};



void SpectrumLibrary::identifyUsingLibrary(const MsParameterStruct* params, 
										   const Config* config,
										   float maxPValue) const
{
	Cluster::setTolerances(config->getTolerance());
	assert(params->spectraListToTest.length()>0);
	vector<string> testPaths;
	readListOfPaths(params->spectraListToTest.c_str(), testPaths);

	size_t  peakBufferSize = 5000;
	Peak*   peakBuffer = new Peak[peakBufferSize];

	for (size_t i=0; i<testPaths.size(); i++)
	{
		SpectraAggregator sa;
		sa.initializeFromSpectraFilePath(testPaths[i].c_str(), config);
		SpectraList sl(sa);
		sl.selectHeaders(params->minMz, params->maxMz, 0, 9999999, params->sqsThreshold);

		if (sl.getNumHeaders() == 0)
			continue;

		string fname;
		getFileNameWithoutExtension(testPaths[i].c_str(), fname);
		ostringstream oss;
		oss << params->outDir + "/" + fname;
		if (params->batchIdx>=0)
			oss << "_" << params->batchIdx;
		oss << "_libres.txt";
		const string resPath = oss.str();

		ofstream ofs(resPath.c_str());
		if (! ofs.good())
			error("Could not open library search results file: ",resPath.c_str());

		if (params->verboseLevel>0)
			cout << i << "\r" << testPaths[i] << "\t";

		ofs << "#\tScan\tTitle\tm/z\tSQS\tSim\tPValue\tDelta\tTitle\tPeptide\tCharge" << endl;
		
		Cluster cluster;
		vector<LibraryEntry> testEntries(sl.getNumHeaders());
		vector<LibrarySearchResult> results(sl.getNumHeaders());
		size_t numToTest=0;

		// read all test spectra in the file and convert them into library entries
		for (size_t j=0; j<sl.getNumHeaders(); j+= 7)
		{
			const SingleSpectrumHeader* header = sl.getSpectrumHeader(j);

			if (header->getOriginalNumPeaks()>peakBufferSize)
			{
				peakBufferSize = 2*header->getOriginalNumPeaks();
				delete [] peakBuffer;
				peakBuffer = new Peak[peakBufferSize];
				cluster.setPeaksPtr(peakBuffer);
			}

			const int numPeaks = cluster.readPeaksToBuffer(sa, header, peakBuffer);
			if (numPeaks == 0)
				continue;

			cluster.initializePeakList(config, true);

			cluster.createNewCluster(j, header, peakBuffer, cluster.getNumPeaks());
			testEntries[numToTest].position = j;
			testEntries[numToTest].header   = *header;
			testEntries[numToTest].distancePeaks = *(cluster.getDistancePeaks());
			numToTest++;
		}

		// sort according to precursor m/z
		sort(testEntries.begin(), testEntries.begin() + numToTest);


		// find best match amongst library spectra
		size_t lowerIdx=0;
		size_t higherIdx=0;
		for (size_t j=0; j<numToTest; j++)
		{
			const mass_t spectrumMz   = testEntries[j].header.getMOverZ();
			const mass_t minLibraryMz = spectrumMz - params->mzWindow;
			const mass_t maxLibraryMz = spectrumMz + params->mzWindow;

//			cout << j << "\t" << fixed << setprecision(3) << spectrumMz << "\t[" << minLibraryMz << " - " << maxLibraryMz << "]" << endl;

			while (lowerIdx<libraryEntries_.size() && libraryEntries_[lowerIdx].header.getMOverZ()<minLibraryMz)
			{
//				cout << "   > > " << lowerIdx << "\t" << libraryEntries_[lowerIdx].header.getMOverZ() << endl;
				lowerIdx++;
			}

			while (higherIdx<libraryEntries_.size() && libraryEntries_[higherIdx].header.getMOverZ()<maxLibraryMz)
			{
//				cout << "   ] ] " << higherIdx << "\t" << libraryEntries_[higherIdx].header.getMOverZ() << endl;
				higherIdx++;
			}

			results[j].testEntry = &testEntries[j];
			results[j].matchEntry = NULL;
			results[j].pvalue = 1.0;
			results[j].similarity = 0.0;

//			cout << j << "\t" << spectrumMz << "\t" << lowerIdx << " - " << higherIdx << endl;
	
			if (lowerIdx == libraryEntries_.size())
				continue;

			const DistancePeakList* testDistancePeaks = &testEntries[j].distancePeaks;
			for (size_t idx = lowerIdx; idx<higherIdx; idx++)
			{
				const float sim = computeSimilarity(testDistancePeaks, 
													&libraryEntries_[idx].distancePeaks,
													Cluster::getPeakIndexTolerance());
													
				if (sim > results[j].similarity)
				{
					results[j].matchEntry = &libraryEntries_[idx];
					results[j].similarity = sim;
				}
			}

			if (results[j].similarity>0.0)
				results[j].pvalue = computeMatchPValue(results[j].similarity, (higherIdx - lowerIdx + 1));
		}
		
		// sort results according to position (matches the order of scans in the files)
		sort(results.begin(), results.begin()+numToTest);

		// ouput results
		for (size_t j=0; j<numToTest; j++)
		{
			if (! results[j].matchEntry)
				continue;

			const SingleSpectrumHeader& testHeader = results[j].testEntry->header;
		
			const int scan = testHeader.getScanNumber();
			const int pos  = testHeader.getIndexInFile();
			const mass_t mz = testHeader.getMOverZ();

			ofs << j << "\t" << (scan>0 ? scan : pos) << "\t" << testHeader.getTitle() << "\t" << mz << "\t" 
				<< testHeader.getSqs() << "\t" << results[j].similarity << "\t" << results[j].pvalue << "\t";

			if (results[j].matchEntry)
			{
				const SingleSpectrumHeader& matchHeader = results[j].matchEntry->header;
				ofs << mz - matchHeader.getMOverZ() << "\t" << matchHeader.getTitle();
				if (matchHeader.getPeptideStr().length()>0)
				{
					ofs << "\t" << matchHeader.getPeptideStr() << "\t" << matchHeader.getCharge();
				}
				else
					ofs << matchHeader.getDatasetIndex() << "\t" << matchHeader.getSpectraFileIndexInList() 
						<< "\t" << matchHeader.getScanNumber();
			}
			else
			{
				ofs << "-";
			}
			ofs << endl;
		}

		ofs.close();

		if (params->verboseLevel>0)
			cout << numToTest << endl;
	}
}



void SpectrumLibrary::learnCdfsForPValues(const string& backgroundSpectra, const string& testSpectraList,
								   const Config* config)
{
	Cluster::setTolerances(config->getTolerance());
	SpectraAggregator loadSa, testSa;
	vector<Cluster> loadClusters, testClusters;
	map<Annotation, vector<size_t> > loadAnns, testAnns;
	
	readAnnotatedSpectraIntoClusters(backgroundSpectra, config, loadSa, loadClusters, loadAnns);
	readAnnotatedSpectraIntoClusters(testSpectraList, config, testSa, testClusters, testAnns);
	sort(loadClusters.begin(), loadClusters.end());
	sort(testClusters.begin(), testClusters.end());

	// fix annotations after sorting
	loadAnns.clear();
	for (size_t i=0; i<loadClusters.size(); i++)
	{
		Annotation ann;
		ann.peptideStr = loadClusters[i].getHeader()->getPeptideStr();
		ann.charge	   = loadClusters[i].getHeader()->getCharge();
		assert(ann.peptideStr.length()>0 && ann.charge>0);
		loadAnns[ann].push_back(i);
	}

	// check if there are spectra that appear in both sets!
	size_t loadIdx=0, testIdx=0;
	while (loadIdx<loadClusters.size() && testIdx<testClusters.size())
	{
		if (loadClusters[loadIdx].getClusterMOverZ() == testClusters[testIdx].getClusterMOverZ())
		{
			if (loadClusters[loadIdx].getHeader()->getTitle() == testClusters[testIdx].getHeader()->getTitle())
			{
				cout << "Spectrum " << loadClusters[loadIdx].getHeader()->getTitle() << endl;
				error("Appears in both sets! Make sure the library and test spectra are not the same!");
			}
			loadIdx++;
			testIdx++;
			continue;
		}
		if (loadClusters[loadIdx].getClusterMOverZ() < testClusters[testIdx].getClusterMOverZ())
		{
			loadIdx++;
		}
		else
			testIdx++;
	}

	const int lastLoadIdx = loadClusters.size()-1;
	int midIdx=0;

	vector<double> simBinCounts(NUM_CDF_BINS,0);
	
	int numTestCases=0;
	for (size_t testIdx=0; testIdx<testClusters.size(); testIdx++)
	{
		const Cluster& testCluster = testClusters[testIdx];
		const mass_t currentMz = testCluster.getClusterMOverZ();
		numTestCases++;
	
		while (midIdx < lastLoadIdx && loadClusters[midIdx].getClusterMOverZ() < currentMz)
			midIdx++;	
	
		// compute similarity with one of the good spectra
		int lowIdx = midIdx - 10000;
		if (lowIdx<0)
			lowIdx = 0;
		int highIdx = lowIdx + 10000;
		if (highIdx>lastLoadIdx)
		{
			highIdx = lastLoadIdx;
			lowIdx = lastLoadIdx - 20000;
		}
		

			
		for (size_t j = lowIdx; j<highIdx; j++)
		{
			const float similarity = computeSimilarity(testCluster.getDistancePeaks(),
													 loadClusters[j].getDistancePeaks(),
													 Cluster::getPeakIndexTolerance());

			const size_t binIdx = computeCdfBin(similarity);
			simBinCounts[binIdx]++;	
		}
	}

	// compute cdf
	cdf_.clear();
	cdf_.resize(NUM_CDF_BINS,0);
	cdf_[0]=simBinCounts[0];
	for (size_t i=1; i<NUM_CDF_BINS; i++)
		cdf_[i]=cdf_[i-1]+simBinCounts[i];
	const double total = cdf_[NUM_CDF_BINS-1];
	for (size_t i=0; i<NUM_CDF_BINS; i++)
		cdf_[i] /= total;

	// if cdf is 1 , make it increase to that value with a fixed slope
	// (to avoid having p-values that are 0)
	size_t i;
	for (i=NUM_CDF_BINS-1; i>0; i--)
		if (cdf_[i]<1.0)
			break;

	if (i<NUM_CDF_BINS-2)
	{
		const double delta = (1.0 - cdf_[i])/(NUM_CDF_BINS - i + 1);
		i++;
		for ( ; i<NUM_CDF_BINS; i++)
			cdf_[i] = cdf_[i-1] + delta;
	}


	writeCdfs(config);

	// print report
	cout << "Pvalues learned from " << 20000*numTestCases << " pairs:" << endl;
	cout << "Sim";
	for (int r=1; r<16; r+=2)
		cout << "\t" << (1<<r) << "\t";
	cout << endl;
	for (int i=0; i<NUM_CDF_BINS; i+=2)
	{
		cout << fixed << setprecision(2) << i*0.005;
		for (int r=1; r<16; r+=2) 
			cout << "\t" << scientific << setprecision(3) << computeMatchPValue(i*0.005,(1<<r));
		cout << endl;
	}
}

void SpectrumLibrary::writeCdfs(const Config* config) const
{
	string pValueFile = config->get_resource_dir() + "/" + config->getModelName() + "_lib_cdfs.txt";
	ofstream ofs(pValueFile.c_str());
	if (! ofs.good())
		error("Couldn't open library pvalue file for writing: ",pValueFile.c_str());
	assert(cdf_.size() == NUM_CDF_BINS);
	ofs << cdf_.size() << endl;
	ofs << scientific << setprecision(6);
	for (size_t i=0; i<cdf_.size(); i++)
		ofs << 1.0 - cdf_[i] << endl;
	ofs.close();
}

void SpectrumLibrary::loadCdfs(const Config *config)
{
	string pValueFile = config->get_resource_dir() + "/" + config->getModelName() + "_lib_cdfs.txt";
	ifstream ifs(pValueFile.c_str());
	if (! ifs.good())
		error("Couldn't open library pvalue file for reading: ",pValueFile.c_str());

	size_t size;
	ifs >> size;
	cdf_.resize(size,1.0);
	assert(size == NUM_CDF_BINS);
	for (size_t i=0; i<size; i++)
	{
		ifs >> cdf_[i];
		assert(cdf_[i]<=1.0);
		cdf_[i] = 1.0 - cdf_[i];
	}
	ifs.close();
}


void SpectrumLibrary::benchmarkLibrary(const Config* config,
									   const MsParameterStruct* params)
{
	Cluster::setTolerances(config->getTolerance());

	if (params->spectraListToLoad.length() == 0 || params->spectraListToTest.length() == 0)
		error("Must supply list of files to load and test!");

	// the benchmark tests the identifications at different peptide density levels
	const int peptideDensities[]={10,100,1000,10000};
	const size_t numPeptideDensities=sizeof(peptideDensities)/sizeof(int);

	SpectraAggregator loadSa, testSa;
	vector<Cluster> loadClusters, testClusters;
	map<Annotation, vector<size_t> > loadAnns, testAnns;
	
	readAnnotatedSpectraIntoClusters(params->spectraListToLoad, config, loadSa, loadClusters, loadAnns);
	readAnnotatedSpectraIntoClusters(params->spectraListToTest, config, testSa, testClusters, testAnns);
	sort(loadClusters.begin(), loadClusters.end());
	sort(testClusters.begin(), testClusters.end());

	// fix annotations after sorting
	loadAnns.clear();
	for (size_t i=0; i<loadClusters.size(); i++)
	{
		Annotation ann;
		ann.peptideStr = loadClusters[i].getHeader()->getPeptideStr();
		ann.charge	   = loadClusters[i].getHeader()->getCharge();
		assert(ann.peptideStr.length()>0 && ann.charge>0);
		loadAnns[ann].push_back(i);
	}

	// check if there are spectra that appear in both sets!
	size_t loadIdx=0, testIdx=0;
	while (loadIdx<loadClusters.size() && testIdx<testClusters.size())
	{
		if (loadClusters[loadIdx].getClusterMOverZ() == testClusters[testIdx].getClusterMOverZ())
		{
			if (loadClusters[loadIdx].getHeader()->getTitle() == testClusters[testIdx].getHeader()->getTitle())
			{
				cout << "Spectrum " << loadClusters[loadIdx].getHeader()->getTitle() << endl;
				error("Appears in both sets! Make sure the library and test spectra are not the same!");
			}
			loadIdx++;
			testIdx++;
			continue;
		}
		if (loadClusters[loadIdx].getClusterMOverZ() < testClusters[testIdx].getClusterMOverZ())
		{
			loadIdx++;
		}
		else
			testIdx++;
	}

	const int lastLoadIdx = loadClusters.size()-1;
	int midIdx=0, lowerWindowIdx=0, higherWindowIdx=0;
	size_t numMissingInLib =0;

	vector< vector<double> > recall(numPeptideDensities+1, vector<double>(101,0.0));
	vector< vector<double> > precision(numPeptideDensities+1, vector<double>(101,0.0));
	vector<string> titles(numPeptideDensities+1);
	ostringstream oss;
	oss << "Window +/- " << setprecision(2) << params->mzWindow;
	titles[0]=oss.str();

	for (size_t i=0; i<numPeptideDensities; i++)
	{
		ostringstream oss;
		oss << "Density " << peptideDensities[i];
		titles[1+i] = oss.str();
	}
	
	int numTestCases=0;
	for (size_t testIdx=0; testIdx<testClusters.size(); testIdx++)
	{
		const Cluster& testCluster = testClusters[testIdx];
		const mass_t currentMz = testCluster.getClusterMOverZ();
		const mass_t lowerMz = currentMz - params->mzWindow;
		const mass_t higherMz = currentMz + params->mzWindow;
	
		if (! testCluster.getIndInPlay())
			continue;

		// make sure that at least one of the library spectra is a correct one
		Annotation ann;
		ann.charge = testCluster.getHeader()->getCharge();
		ann.peptideStr = testCluster.getHeader()->getPeptideStr();

		map<Annotation, vector<size_t> >::const_iterator it = loadAnns.find(ann);
		if (it == loadAnns.end())
		{
			numMissingInLib++;
			continue;
		}

		numTestCases++;
		// set idxs for each range
		vector<int> lowerIdxs(numPeptideDensities+1);
		vector<int> higherIdxs(numPeptideDensities+1);

		while (midIdx < lastLoadIdx && loadClusters[midIdx].getClusterMOverZ() < currentMz)
			midIdx++;	
		while (lowerWindowIdx < lastLoadIdx && loadClusters[lowerWindowIdx].getClusterMOverZ()< lowerMz)
			lowerWindowIdx++;
		while (higherWindowIdx < lastLoadIdx && loadClusters[higherWindowIdx].getClusterMOverZ()< higherMz)
			higherWindowIdx++;
	
		lowerIdxs[0] = lowerWindowIdx; 
		higherIdxs[0] = higherWindowIdx;

		for (size_t r=0; r<numPeptideDensities; r++)
		{
			int low = testIdx - peptideDensities[r]/2;
			if (low<0)
				low=0;
			int high = low + peptideDensities[r];
			if (high>lastLoadIdx)
				high=lastLoadIdx;
			if (high - low < peptideDensities[r])
				low = high - peptideDensities[r];
			if (low<0)
				low=0;
			lowerIdxs[r+1]=low;
			higherIdxs[r+1]=high;
		}

		// compute similarity with one of the good spectra
		for (size_t r=0; r<=numPeptideDensities; r++)
		{
			vector<int> goodBelow(101,0);
			vector<int> badBelow(101,0);
			int numTestedForSpectrum = 0;
		
			const size_t otherIdx = it->second[it->second.size()/2];
			assert(testCluster.getHeader()->getPeptideStr() == loadClusters[otherIdx].getHeader()->getPeptideStr());

			// test a good spectrum so we don't only test bad ones
			float maxCorrectSimilarity=0.0;
			if (otherIdx<lowerIdxs[r] || otherIdx>=higherIdxs[r])
			{
				float similarity =  computeSimilarity(testCluster.getDistancePeaks(),
												loadClusters[otherIdx].getDistancePeaks(),
												Cluster::getPeakIndexTolerance());
				if (similarity>maxCorrectSimilarity)
					maxCorrectSimilarity = similarity;
				const size_t binIdx = static_cast<size_t>(100.0*similarity+0.499);
				numTestedForSpectrum++;
				for (size_t i=0; i<=binIdx; i++)
					goodBelow[i]++;
			}

			
			for (size_t j = lowerIdxs[r]; j<=higherIdxs[r]; j++)
			{
				if (! loadClusters[j].getIndInPlay())
					continue;

				const float similarity = computeSimilarity(testCluster.getDistancePeaks(),
													 loadClusters[j].getDistancePeaks(),
													 Cluster::getPeakIndexTolerance());
				
				bool isCorrect = (loadClusters[j].getHeader()->getPeptideStr() == 
					    		  testCluster.getHeader()->getPeptideStr());
				const size_t binIdx = static_cast<size_t>(100.0*similarity+0.499);
				numTestedForSpectrum++;
				if (isCorrect)
				{
					if (similarity>maxCorrectSimilarity)
						maxCorrectSimilarity = similarity;
					for (size_t i=0; i<=binIdx; i++)
						goodBelow[i]++;
				}
				else
					for (size_t i=0; i<=binIdx; i++)
						badBelow[i]++;
			}

			for (size_t i=0; i<101; i++)
			{
				int sum = goodBelow[i]+badBelow[i];
				if (sum>0)
				{
					precision[r][i] += (static_cast<double>(goodBelow[i])/static_cast<double>(sum));
				}
				else
					precision[r][i]++;
			}

			const size_t binIdx = static_cast<size_t>(100.0*maxCorrectSimilarity+0.499);
			for (size_t i=0; i<binIdx; i++)
				recall[r][i]++;
		}
	}

	if (numMissingInLib>0)
	{
		cout << "Skipped " << numMissingInLib << "/" << testClusters.size() << endl;
		if (numTestCases == 0)
		{
			cout << "No test cases!" << endl;
			return;
		}
	}

	//produce reports
	// recall = numGood/numTestCases
	// precision = numGood/(numGood+numBad)
	for (size_t r=0; r<=numPeptideDensities; r++)
	{
		for (size_t i=0; i<101; i++)
		{
			precision[r][i] /= static_cast<double>(numTestCases);
			recall[r][i]   /= static_cast<double>(numTestCases);
		}

		// make sure the vectors are monotonic
		for (size_t i=1; i<101; i++)
		{
			if (precision[r][i]<precision[r][i-1])
				precision[r][i]=precision[r][i-1];
			if (recall[r][i]>recall[r][i-1])
				recall[r][i]=recall[r][i-1];
		}
	}

	// print roc curves
	cout << setprecision(4) << fixed;
	cout << endl;
	for (size_t i=0; i<=numPeptideDensities; i++)
		cout << "\t" << titles[i];
	cout << endl;
	cout << "Sim";
	for (size_t i=0; i<=numPeptideDensities; i++)
		cout << "\tPrec\tRecall";
	cout << endl;
	for (int i=0; i<101; i++)
	{
		cout << 0.01*i;
		for (int r=0; r<=numPeptideDensities; r++)
			cout << "\t" << precision[r][i] << "\t" << recall[r][i];
		cout << endl;
	}

	cout << "AUC:";
	for (size_t i=0; i<=numPeptideDensities; i++)
		cout << "\t" << computeRocAuc(precision[i],recall[i]) << "\t";
	cout << endl;
}


