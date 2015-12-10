Clustering Mass Spectra (MS-Clustering) 
Copyright 2007, The Regents of the University of California. All Rights Reserved. 
Created by Ari Frank (arf@cs.ucsd.edu).

MS-Clustering is designed to rapidly cluster large MS/MS datasets. The program merges similar spectra (having similar m/z values – within a given tolerance), and creates a single consensus spectrum as a representative. The input formats accepted are: dta, mgf, mzXML. The output format is mgf.

Setting up:
-----------
In the same directory as the MS-Clustering executable you should place the “Models” directory. You should also make directories “tmp” and “out” in that location (though these directories can reside elsewhere if their path is specified using the –tmp_dir and –out_dir). If running on Unix/Linux, run the command "dos2unix Models/*" to fix any Windows/Unix issues.

How it works (very short version):
----------------------------------
The clustering is performed in two stages. First all files are read and stored in binary format (dat files) in the “tmp” directory (these files might require manual deletion in the end of the run “rm tmp/*”). In the second stage, the program reads the files according to m/z and clusters “slices” of the dataset (typically slices of 2.5 Da). An additional filtration stage can be used to remove spectra of very low quality (by using the flag –filter_model_name CID_IT_TRYP). This stage uses a quality filtration program designed for LTQ ion-trap data (models for other instruments can be created however they will require training data). Each “slice” is clustered independently (though unassigned spectra in slice x can be added to clusters in slice x+1 to reduce cluster fragmentation). The clusters are outputted into the out directory in the form of mgf files, with 10000 spectra per file. Each file also has a .clust.txt file which serves as an indexing file that states who the original members of the cluster were (each such spectrum has a line with three numbers, the file index, the scan/location in file, and the m/z).The default setting keeps clusters of all sizes (to control the outputted sizes of clusters use the "-min_size" and "-max_size" flags.  


Filtering spectra (and avoiding over-filtering):
------------------------------------------------

Removing low quality spectra can greatly reduce the number of clusters outputted by the algorithm. Paradoxically, it can also help increase the number of identifications (since there are less spurious hits to the decoy database). You can invoke quality filtration by adding the flag “-filter_model_name LTQ_LOW_TRYP” (currently there is only a model for LTQ ion-trap data). This flag will discard all spectra of low quality that are not likely to yield peptide identifications. The default quality threshold is 0.1, however there might be times when this threshold should be decreased. Sometimes the filtering models may not be appropriate for certain datasets that were generated on different instruments than the ones used to generate the training data. In such cases the quality scores might be lower than they should and too many spectra get filtered. With typical runs ~ 40%-60% of the spectra should fall into clusters (due to internal settings, the value should be approximately 60% for datasets of less than 2 million and ~50% for larger datasets). If at the end of the clustering you notice that the ratio of spectra that fell into clusters is below 40%, the program might be over-filtering (which can lead to a loss of identifications). In such cases you can decrease the quality threshold with the flag "-min_filter_prob 0.05" (or some other value below 0.1).


Clustering large datasets:
--------------------------

The default settings are appropriate for datasets of up to 5 million spectra. For larger datasets you might need to increase the number of spectra that are simultaneously processed with the flag "-specs_per_slice". The default value is 25000, and should grow linearly with the total number of spectra processed (so for 10 million spectra, you should use 50000, etc.). Note that increasing this parameter greatly increases the memory demands (setting the parameter to 50000 requires approximately 1 GB of RAM). Note that if this parameter is too low the clustering program will still work however there might be more incidents of cluster fragmentation.

If you have several processing nodes to work with (and ample memory for each node), large clustering jobs can be parallelized. This can be done by first creating dat files for the whole dataset ("-make_only_dat" flag), and then splitting the dat files for separate jobs (in such a case the lists should be given with "-dat_list" flag and not the "-list" flag). To having multiple output files with the same name, each job can be given a unique batch index with the flag "-batch <xx>" (where xx is a unique integer for each job).


Output:
-------

The program writes to types of files to the output directory. The first are the clustered spectra files (in MGF format). For each MGF file, the program also writes a "clust.txt" file which contains a mapping of the cluster members (original members)for each spectrum in the MGF file. For each scan in the MGF, the "clust.txt" file has the following entries: <Cluster name> <cluster size> <cluster m/z>. Following this line, there is a line for each member of the cluster with the following fields: <file index> <scan number> <charge - if known, otherwise 0>.


Important parameters: 
---------------------

(all possible parameters can be seen when MS-Clustering is run without arguments).

-list <path_to_list_file> .  The list file is a text file that holds the paths to the input files (full path to each mgf/mzXML/dta), one file per line.

-dat_list <path_to_list_file> . If you have previously created dat files, you can skip their recreation and use them directly.

Optional parameters:
-------------------
-no_normalize        - do not normalize the peaks intensities (otherwise intensities are normalized to 1000).
-make_dat_only       - only makes dat files (from the input files)
-use_spectrum_charge - only cluster spectra that have the same charge (as assigned in the file), also outputs a "CHARGE=" field in the mgf
-tolerance    <xx>  - peak tolerance in Da. (default xx=0.4 Da. – peaks will be joined if they are approximately within this distance)
-slice_width  <xx>  - the width of the m/z clustering window (default xx=2.5 Da.)
-similarity   <xx>  - the minimal similarity for clustering (default xx=0.55, values that can be used 0.1-0.99)
-sim_peaks    <xx>  - the number of peaks to use for similarity calculations per 1000 Da. of mass (default xx=15)
-min_m_over_z <xx>  - minimal m/z to cluster (default xx=0)
-max_m_over_z <xx>  - maximal m/z to cluster (default xx=10000)
-min_size     <xx>  - the minimal size for an outputted cluster (default xx=2)
-max_size     <xx>  - the maximal size for an outputted cluster (default xx=10000)
-tmp_dir     <xx>  - directory where the temporary files are written (default xx=".\tmp")
-out_dir      <xx>  - directory where the cluster files are written (default xx=".\out")
-name         <xx>  - name to be given cluster files (default xx="Cluster")
-batch        <xx>  - the index given to this batch of clusters (default xx=0)
-file_idx_start <xx> - the number to add to the file indices (only use if dat files are created from split lists).
-specs_per_slice <xx> - number of expected spectra per slice (default xx = 25000). This number is suited for about 4M spectra. If a larger number of spectra is clustered, you should increase this figure (linearly). This parameter determines how much memory is used by the program.
-filter_model_name <xx> - name of model PMC and SQS model files (should reside in ".\Models")
-model_dir       <path> - directory where model files are kept (default ./Models)
-min_filter_prob <xx> - the minimal probability score to keep a spectrum ([0-1], default xx=0.075)
-assign_charges         - if filtering is performed, then a charge is typically asigned to each cluster, using this flag adds the CHARGE field to the outputted MGF.
-output_peak_density <xx> - the number of peaks per 100 Da. that will be in the outputted clusters' spectra, default xx=8.
-mass_exclude_range <xx yy> - mass ranges to exclude for similarity computation (e.g.,for iTraq) use  "-mass_exclude_range 113.5 118.5" (note: if normalization is used, the iTraq peaks of clusters may not represent the correct iTraq intensities).

