#include "QuickClustering.h"

void QCOutputter::init(string _name , string _dir, int _batch_idx,
			  mass_t min_m_over_z, mass_t max_m_over_z,
			  float min_similarity, int min_cluster_size)
{
	batch_idx = _batch_idx;
	

	if (batch_idx>=0)
	{

		ostringstream oss;
		oss << batch_idx;
		batch_str = oss.str();

		dir = _dir;
		name = _name;
		string sum_name = dir + "/" + name + "_" + batch_str + "_sum.txt";
		string list_name = dir + "/" + name + "_" + batch_str + "_list.txt";
		string param_name = dir + "/" + name + "_" + batch_str + "_params.txt";

		summary_stream.open(sum_name.c_str(),ios::out);               
		file_list_stream.open(list_name.c_str(),ios::out);
		fstream param_stream(param_name.c_str(),ios::out);

		
		if (! summary_stream.is_open() || 
			! file_list_stream.is_open() ||
			! param_stream.is_open() )
		{
			cout << "Error: couldn't open outputter file streams!" << endl;
			exit(1);
		}

		param_stream << "batch:     " << batch_idx << endl;
		param_stream << "out dir:   " << dir << endl;
		param_stream << "min m/z:   " << min_m_over_z << endl;
		param_stream << "max m/z:   " << max_m_over_z << endl;
		param_stream << "min similarity: " << min_similarity << endl;
		param_stream << "min cluster size: " << min_cluster_size << endl;

		param_stream.close();
	}
	else
	{
		dir = _dir;
		name = _name;
		cout << "Sending output MGF files called " << name << " to dir: " << dir << endl;
	}
}

QCOutputter::~QCOutputter()
{
	if (mgf_stream.is_open())
		mgf_stream.close();

	if (cluster_file_stream.is_open())
		cluster_file_stream.close();

	if (summary_stream.is_open())
		summary_stream.close();

	if (file_list_stream.is_open())
		file_list_stream.close();
}


void QCOutputter::output_basic_spectrum_to_mgf(BasicSpectrum &bs, Config *config)
{
	// check if a new file should be opened
	if (total_spectra_counter == 0 ||
		spectra_counter == NUM_CLUSTERS_PER_FILE)
	{
		if (mgf_stream.is_open())
		{
			mgf_stream.close();
			cluster_file_stream.close();
		}

		file_counter++;
		ostringstream oss;
		oss << file_counter;
		mgf_name = name +  "_" + oss.str() + ".mgf";
		string mgf_path = dir + "/" + mgf_name;
		mgf_stream.open(mgf_path.c_str(),ios::out);

		string list_name = dir + "/" + name +  "_" + oss.str() + "_list.txt";
		cluster_file_stream.open(list_name.c_str(),ios::out);

		spectra_counter=0;
	}

	const SingleSpectrumFile * ssf = bs.ssf;
	if (ssf->type == MZXML)
	{
		MZXML_single* mzxml_ssf = (MZXML_single *)ssf;

		cluster_file_stream << mzxml_ssf->file_idx << "\t" << mzxml_ssf->scan_number << "\t" <<
			mzxml_ssf->m_over_z << "\t" << mzxml_ssf->charge ;
	}
	else if (ssf->type == DAT)
	{
		DAT_single* dat_ssf = (DAT_single *)ssf;
		cluster_file_stream << dat_ssf->mzxml_file_idx << "\t" << dat_ssf->scan_number << "\t" <<
				dat_ssf->m_over_z << "\t" << dat_ssf->charge;
	}
	else if (ssf->type == MGF)
	{
		MGF_single* mgf_ssf = (MGF_single *)ssf;
		cluster_file_stream << mgf_ssf->file_idx << "\t" << mgf_ssf->idx_in_file << "\t" <<
			mgf_ssf->m_over_z << "\t" << mgf_ssf->charge;
	}
	else if (ssf->type == DTA)
	{
		cluster_file_stream << ssf->single_name << "\t" << ssf->m_over_z << "\t" << ssf->charge;
	}

	cluster_file_stream << "\t" << ssf->sqs << endl;
	bs.output_to_mgf(mgf_stream,config);

	spectra_counter++;
	total_spectra_counter++;

}



/****************************************************************
	Writes the cluster spectrum to the output file, and adds the 
	relevant info to the summary and list files
*****************************************************************/
void QCOutputter::output_cluster_spectrum(ClusterSpectrum& cs, bool ind_write_charge)
{
	int i;

	// check if a new file should be opened
	if (total_spectra_counter == 0 ||
		spectra_counter == NUM_CLUSTERS_PER_FILE)
	{
		if (mgf_stream.is_open())
		{
			mgf_stream.close();
			cluster_file_stream.close();
		}

		file_counter++;
		ostringstream oss;
		oss << file_counter;
		mgf_name = name + "_" + batch_str + "_" + oss.str() + ".mgf";
		string mgf_path = dir + "/" + mgf_name;
		mgf_stream.open(mgf_path.c_str(),ios::out);

		file_list_stream << mgf_path << endl;

		cluster_file_name = dir + "/" + name + "_" + batch_str + "_" + oss.str() + ".clust.txt";
		cluster_file_stream.open(cluster_file_name.c_str(),ios::out);

		spectra_counter=0;
	}

	cs.make_title(name,batch_idx,total_spectra_counter);

	// write spectrum to mgf
	cs.write_spectrum_to_mgf(mgf_stream,false,ind_write_charge);


	Config *config = cs.get_config();

	// write cluster info to cluster file
	cluster_file_stream << cs.get_title() << " " << cs.basic_spectra.size() << " " << cs.m_over_z << endl;
	for (i=0; i<cs.basic_spectra.size(); i++)
	{
		const SingleSpectrumFile * ssf = cs.basic_spectra[i].ssf;
		if (ssf->type == MZXML)
		{
			MZXML_single* mzxml_ssf = (MZXML_single *)ssf;


			cluster_file_stream << mzxml_ssf->file_idx << "\t" << mzxml_ssf->scan_number << "\t" <<
				mzxml_ssf->m_over_z << "\t" << mzxml_ssf->charge;
		}
		else if (ssf->type == DAT)
		{
			DAT_single* dat_ssf = (DAT_single *)ssf;
			cluster_file_stream << dat_ssf->mzxml_file_idx << "\t" << dat_ssf->scan_number << "\t" <<
				dat_ssf->m_over_z << "\t" << dat_ssf->charge;
		}
		else if (ssf->type == MGF)
		{
			MGF_single* mgf_ssf = (MGF_single *)ssf;
			cluster_file_stream << mgf_ssf->file_idx << "\t" << (mgf_ssf->scan_number>0 ? mgf_ssf->scan_number : mgf_ssf->idx_in_file) << 
				"\t" << mgf_ssf->m_over_z << "\t" << mgf_ssf->charge << "\t" << mgf_ssf->single_name;
		}
		else if (ssf->type == DTA)
		{
			cluster_file_stream << ssf->single_name << "\t" << ssf->m_over_z << "\t" << ssf->charge;
		}

		
		cluster_file_stream << endl;
	}
	cluster_file_stream << endl;

	// update summary file
	summary_stream << mgf_name << " " << spectra_counter << " " << cs.basic_spectra.size() << " " << cs.get_title()  << " " << cs.m_over_z << endl;

	spectra_counter++;
	total_spectra_counter++;

}



/****************************************************************
	Writes the cluster spectrum to the output file, and adds the 
	relevant info to the summary and list files
*****************************************************************/
void QCOutputter::output_cluster_spectrum_as_single_pkl(ClusterSpectrum& cs)
{
	int i;

	// check if a new file should be opened
	if (total_spectra_counter == 0)
	{
		if (cluster_file_stream.is_open())
		{
			cluster_file_stream.close();
		}

		file_counter++;
		ostringstream oss;
		oss << file_counter;
		
		cluster_file_name = dir + "/" + name + "_" + batch_str + "_" + oss.str() + ".clust.txt";
		cluster_file_stream.open(cluster_file_name.c_str(),ios::out);

		spectra_counter=0;
	}

	cs.make_title(name,batch_idx,total_spectra_counter);

	char scan_buff1[16];
	char scan_buff2[16];
	char charge_buff[16];

	const int end_scan_num = total_spectra_counter+cs.get_num_basic_spectra();
	sprintf(scan_buff1,"%d",total_spectra_counter);
	sprintf(scan_buff2,"%d",end_scan_num);
	sprintf(charge_buff,"%d",cs.get_charge());

	// write spectrum to mgf
	string pkl_name =  name + ".";
	if (total_spectra_counter<10)
	{
		pkl_name += "000";
	}
	else if (total_spectra_counter<100)
	{
		pkl_name += "00";
	}
	else if (total_spectra_counter<1000)
	{
		pkl_name += "0";
	}

	pkl_name += scan_buff1;
	pkl_name += ".";

	if (end_scan_num<10)
	{
		pkl_name += "000";
	}
	else if (end_scan_num<100)
	{
		pkl_name += "00";
	}
	else if (end_scan_num<1000)
	{
		pkl_name += "0";
	}

	pkl_name += scan_buff2;
	pkl_name += ".";
	pkl_name += charge_buff;
	pkl_name += ".pkl";

	
	string pkl_path = dir + "\\" + pkl_name;

	cs.set_title(pkl_name);

	cs.write_spectrum_to_pkl_single(pkl_path);

	Config *config = cs.get_config();

	// write cluster info to cluster file
	cluster_file_stream << cs.get_title() << " " << cs.basic_spectra.size() << " " << cs.m_over_z << endl;
	for (i=0; i<cs.basic_spectra.size(); i++)
	{
		const SingleSpectrumFile * ssf = cs.basic_spectra[i].ssf;
		if (ssf->type == MZXML)
		{
			MZXML_single* mzxml_ssf = (MZXML_single *)ssf;


			cluster_file_stream << mzxml_ssf->file_idx << "\t" << mzxml_ssf->scan_number << "\t" <<
				mzxml_ssf->m_over_z << "\t" << mzxml_ssf->charge << endl;
		}
		else if (ssf->type == DAT)
		{
			DAT_single* dat_ssf = (DAT_single *)ssf;
			cluster_file_stream << dat_ssf->mzxml_file_idx << "\t" << dat_ssf->scan_number << "\t" <<
				dat_ssf->m_over_z << "\t" << dat_ssf->charge << endl;
		}
		else if (ssf->type == MGF)
		{
			MGF_single* mgf_ssf = (MGF_single *)ssf;
			cluster_file_stream << mgf_ssf->file_idx << "\t" << mgf_ssf->idx_in_file << "\t" <<
				mgf_ssf->m_over_z << "\t" << mgf_ssf->charge << endl;
		}
		else if (ssf->type == PKL)
		{
			PKL_single * pkl_ssf = (PKL_single *)ssf;
		
			cluster_file_stream << pkl_ssf->file_idx << "\t" << pkl_ssf->scan_number << "\t" <<
				pkl_ssf->m_over_z << "\t" << pkl_ssf->charge << endl;
		}
		else if (ssf->type == DTA)
		{
			cluster_file_stream << ssf->single_name << "\t" << ssf->m_over_z << "\t" << 
				ssf->charge << endl;
		}
	}
	cluster_file_stream << endl;

		// update summary file
	summary_stream << pkl_name << spectra_counter << "\t" << cs.basic_spectra.size() << "\t" << cs.m_over_z << "\t" << cs.charge << endl;

	spectra_counter++;
	total_spectra_counter++;
}



/****************************************************************
	Writes the annotations of the cluster to an output.
	Bases the cluster's annotation on the majority of the spectra
	that have an annotation. If there is no peptide that has 75%
	majority, the annotations are not written. Otherwise, the annotaitons
	to all the spectra are set like the majority.

*****************************************************************
void QCOutputter::output_cluster_anns(ClusterSpectrum& cs)
{
	int i;

	if (cs.basic_spectra.size()==0)
		return;

	// check if spectra have annotations
	string pep_str;
	mass_t pep_mass;
	if (! cs.has_majority_annotation(pep_str,pep_mass))
		return;


	int charge=-1;
	mass_t m_over_z = cs.get_m_over_z();
	for (charge=1; charge<20; charge++)
		if (fabs(m_over_z*charge - charge +1 - pep_mass)<10)
			break;
	
	if (charge == 20)
	{
		cout << "Warning: couldn't find charge for " << pep_str << " m/z: " << m_over_z << endl;
		cout << "The peptides mass is: " << pep_mass << endl;
		return;
	}



	// write cluster info to cluster file
	cluster_file_stream << cs.get_title() << " " << cs.basic_spectra.size() << " " << cs.m_over_z << endl;
	for (i=0; i<cs.basic_spectra.size(); i++)
	{
		int file_idx,scan;
		const SingleSpectrumFile * ssf = cs.basic_spectra[i].ssf;
		if (ssf->type == MZXML)
		{
			MZXML_single* mzxml_ssf = (MZXML_single *)ssf;
			file_idx = mzxml_ssf->file_idx;
			scan = mzxml_ssf->scan_number;
		}
		else if (ssf->type == DAT)
		{
			DAT_single* dat_ssf = (DAT_single *)ssf;
			file_idx = dat_ssf->mzxml_file_idx;
			scan     = dat_ssf->scan_number;
		}
		else if (ssf->type == MGF)
		{
			MGF_single* mgf_ssf = (MGF_single *)ssf;
			file_idx = mgf_ssf->file_idx;
			scan     = mgf_ssf->idx_in_file;
		}
		else if (ssf->type == DTA)
		{
			cout << "Error: outputting anns for DTA!" << endl;
			exit(1);
		}


			// check if a new file should be opened
		if (total_spectra_counter == 0)
		{
			string anns_name = name + "_" + batch_str + "_anns_new.txt";
			string anns_path = dir + "/" + anns_name;
			anns_stream.open(anns_path.c_str(),ios::out);

			if (! anns_stream.is_open()  )
			{
				cout << "Error: couldn't open anns file for wirting: " << anns_path << endl;
				exit(1);
			}

			cout << "Opened: " << anns_path << endl;
		}

		if (! anns_stream.good())
		{
			cout << "Error: bad annoation_file_stream!!! " << endl;
			exit(1);

		}


		total_spectra_counter++;

		anns_stream << file_idx << " -1 " << scan << " " << charge << " " << pep_str << endl;
	//	cout << file_idx << " -1 " << scan << " " << charge << " " << pep_str << endl;
	}	
}

*/


