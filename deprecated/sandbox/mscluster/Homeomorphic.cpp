#include "Homeomorphic.h"
#include "AnnotatedSpectrum.h"
#include "DeNovoDP.h"
#include "auxfun.h"





void hgraph::create_graph(Config *_config, Peptide& p)
{
	config = _config;
	org_pep = p;
	int i,j;

	const int num_aa = p.get_num_aas();
	const vector<int>& amino_acids = p.get_amino_acids();
	const vector<int>& session_aas = config->get_session_aas();
	const vector<string>& aa2label = config->get_aa2label();
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	const mass_t peptide_mass = p.get_mass();
	const mass_t p18_mass = peptide_mass + MASS_H2O;
	const mass_t tolerance = config->get_tolerance();


	vector<mass_t> break_masses;
	p.calc_expected_breakage_masses(config,break_masses);

	nodes.resize(2*num_aa+2);
	nodes[0].mass=0;
	nodes[0].type = 0;
	nodes[0].idx =0;

	nodes[1].mass = p18_mass;
	nodes[1].type = 1;
	nodes[1].idx = 0;

	for (i=1; i<=num_aa; i++)
	{
		int n_idx=2*i;

		nodes[n_idx].mass = nodes[n_idx-2].mass + aa2mass[amino_acids[i-1]];
		nodes[n_idx].type = 0;
		nodes[n_idx].idx = i;

		n_idx++;

		nodes[n_idx].mass = p18_mass - nodes[n_idx-1].mass;
		nodes[n_idx].type = 1;
		nodes[n_idx].idx = i;
	}
	sort(nodes.begin(),nodes.end());

	// connect edges (single only)

	i=0;
	j=1;

	while (i<nodes.size()-2)
	{
		int k=i+1;
		while (nodes[k].type != 0)
			k++;

		hedge e;
		e.n1=i;
		e.n2=k;
	
		nodes[i].add_edge(e);
		i=k;
	}

	while (j<nodes.size()-2)
	{
		int k=j+1;
		while (nodes[k].type != 1)
			k++;

		hedge e;
		e.n1=j;
		e.n2=k;
	
		nodes[j].add_edge(e);
		j=k;	
	}



	// connect crossover edges
	for (i=0; i<nodes.size()-2; i++)
	{
		int j;
		for (j=0; j<session_aas.size(); j++)
		{
			const int aa= session_aas[j];
			if (aa==Ile)
				continue;

			int k;
			for (k=i+1; k<nodes.size(); k++)
			{
				if (nodes[i].type != nodes[k].type)
				{
					hedge e;
					e.n1=i;
					e.n2=k;
					e.num_aa=1;
					e.score=1;
					e.is_crossover = (nodes[i].type != nodes[k].type);
					e.aa[0] = aa;

					nodes[i].add_edge(e);
					break;
				}
			}
		}
	}
}


bool hgraph::has_two_paths()
{
	const int num_aa = org_pep.get_num_aas();
	const int goal_node_idx = nodes.size()-2;
	vector<int> out_idxs, prev_idxs;

	bool print = true;

	// init different vectors
	out_idxs.clear();
	out_idxs.resize(nodes.size(),0);
	prev_idxs.clear();
	prev_idxs.resize(nodes.size(),-1);

	// do a DFS
	int p=0;
	int count =0;
    
	while (1)
	{
		while (out_idxs[p]<nodes[p].edges.size())
		{

			// take this edge
			const hedge& e = nodes[p].edges[out_idxs[p]];
			const int node_idx = nodes[e.n2].idx ;
			const hnode& node = nodes[e.n2];

			int old_p=p;
			p=  e.n2;

			prev_idxs[p]=old_p;;
			out_idxs[old_p]++;

			if (p == goal_node_idx )
			{
				// check for forbidden
				count++;

				if (count>1)
					return true;
				
				break;
			}

		}

		// back track
		
		out_idxs[p]=0;
		p=prev_idxs[p];

		if (p<0)
			break;
	}

	return false;
}


void hgraph::get_delta_path_nums(vector<int>& counts, int max_print)
{
	const int num_aa = org_pep.get_num_aas();
	const int goal_node_idx = nodes.size()-2;
	vector<bool> used_idxs;
	vector<int> out_idxs, prev_idxs;

	// init different vectors
	counts.clear();
	counts.resize(num_aa*2,0);
	used_idxs.resize(num_aa+1,false);
	used_idxs[0]=true;
	out_idxs.clear();
	out_idxs.resize(nodes.size(),0);
	prev_idxs.clear();
	prev_idxs.resize(nodes.size(),-1);

	// do a DFS
	int score=0;
	int p=0;

//	int n;
//	for (n=0; n<nodes.size(); n++)
//		cout << n << " " << nodes[n].idx << (nodes[n].type==0? "  " : "' ") << 
//		        nodes[n].edges.size() << endl;

	while (1)
	{
		while (out_idxs[p]<nodes[p].edges.size())
		{

			// take this edge
			const hedge& e = nodes[p].edges[out_idxs[p]];
			const int node_idx = nodes[e.n2].idx ;
			const hnode& node = nodes[e.n2];

			if (used_idxs[node_idx]) // forbidden pair check
			{
				out_idxs[p]++;
				continue;
			}

			int old_p=p;
			p=  e.n2;
			score+= e.score;
			prev_idxs[p]=old_p;
			used_idxs[node_idx]=true;
			out_idxs[old_p]++;

			if (p == goal_node_idx )
			{
				// check for forbidden
				int bin = num_aa-score;
				counts[bin]++;
				if (bin<=max_print)
				{
					int q=p;
					vector<int> ns;
					while (q>=0)
					{
						ns.push_back(q);
						q=prev_idxs[q];
					}

					sort(ns.begin(),ns.end());

					// count alternating
					int type=0;
					int alt_count =0;
					int i;
					for (i=0; i<ns.size(); i++)
					{
						if (nodes[ns[i]].type != type)
						{
							alt_count++;
							type = nodes[ns[i]].type;
						}
					}

					if (alt_count>4 && bin == 0)
					{
						cout << alt_count << "  PEP: ";
						for (q=0; q<ns.size(); q++)
						{
							cout << nodes[ns[q]].idx;
							if (nodes[ns[q]].type == 1)
								cout << "'";
							cout << " ";
						}
						cout << endl;
					}
				}
				break;
			}

		}

		// back track
		
		out_idxs[p]=0;
		used_idxs[nodes[p].idx]=false;
		p=prev_idxs[p];

		if (p<0)
			break;

		score -= nodes[p].edges[out_idxs[p]-1].score;
	}


}




void hgraph::print() const
{
	int i;
	int type;
	for (type =0; type<=1; type++)
	{
		for (i=0; i<nodes.size(); i++)
		{
			if (nodes[i].type == type)
			{
				const vector<hedge>& edges = nodes[i].edges;
				int j;
				for (j=0; j<edges.size(); j++)
					if (edges[j].is_crossover)
					{
						if (type == 0)
						{
							cout << nodes[edges[j].n1].idx << " -> " << nodes[edges[j].n2].idx << "' *" << endl; 
						}
						else
							cout << nodes[edges[j].n1].idx << "' ->" << nodes[edges[j].n2].idx << " *" <<endl;
					}
				/*	else
					{
						if (nodes[edges[j].n1].type == 0)
						{
							cout << nodes[edges[j].n1].idx << " -> " << nodes[edges[j].n2].idx << endl;
						}
						else
						{
							cout << nodes[edges[j].n1].idx << "' -> " << nodes[edges[j].n2].idx << "'" << endl;
						}
					}*/
			}
			
		}
	}
}


void hexp1(Config *config)
{
	int i,size;
	int total = 5000;

	for (size =5; size<=20; size++)
	{
		int count =0;
		for (i=0; i<total; i++)
		{
			Peptide p;
			vector<int> counts;
			p.generate_random_peptide(config,size);

			hgraph hg;
			hg.create_graph(config,p);
		
			if (hg.has_two_paths())
			{
				count++;
				hg.print();
				cout << endl;
			}

		
		}

		cout << size << " " << count << "/" << total << "=" << (double)count/(double)total << endl;
	}
}





//-------------------------------------------


struct cross_jump {
	cross_jump() : idx1(-1), idx2(-1), aa(-1), aa2(-1) {}
	int idx1,idx2;
	int aa,aa2;
};

struct node_mass {
	bool operator< ( const node_mass& other) const
	{
		return mass<other.mass;
	}
		
	int source;
	mass_t mass;
	int idx;
};

int get_num_paths(Config *config, Peptide& peptide, 
						 mass_t tolerance, bool check_forbidden, bool allow_double,
						 bool print)
{
	const vector<int>& session_aas = config->get_session_aas();
	const vector<string>& aa2label = config->get_aa2label();
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	const mass_t peptide_mass = peptide.get_mass();

	vector<mass_t> break_masses;
	vector<mass_t> other_masses;
	peptide.calc_expected_breakage_masses(config,break_masses);

	int i;
	for (i=0; i<break_masses.size(); i++)
	{
		other_masses.push_back(peptide_mass - break_masses[i]+18);
		if (print)
			cout << i<< " " << setw(8) << left <<break_masses[i] << " " << setw(8) << left << other_masses[i] << endl;
	}

	vector<cross_jump> over,back;
	over.clear();
	back.clear();
	int length = break_masses.size();

	for (i=0; i<break_masses.size(); i++)
	{
		int j;
		for (j=0; j<other_masses.size(); j++)
		{
			int a;
			for (a=0; a<session_aas.size(); a++)
			{
				const int& aa = session_aas[a];
				if (aa == Ile || aa == Lys)
					continue;

				if (fabs(break_masses[i] + aa2mass[aa] - other_masses[j])<=tolerance)
				{
					cross_jump cj;
					cj.aa = aa;
					cj.idx1 = i;
					cj.idx2 = j;

					over.push_back(cj);					
				}
				else
				{
					if (allow_double)
					{
						int b;
						for (b=0; b<session_aas.size(); b++)
						{
							const int& bb = session_aas[b];
							if (bb == Ile || bb == Lys)
								continue;

							if (! config->is_allowed_double_edge(aa,bb))
								continue;

							if (fabs(break_masses[i] + aa2mass[aa] + aa2mass[bb]  - other_masses[j])<=tolerance)
							{
								cross_jump cj;
								cj.aa = aa;
								cj.aa2 = bb;
								cj.idx1 = i;
								cj.idx2 = j;
								over.push_back(cj);
							}				
						}
					}
				}
			}
		}
	}


	for (i=0; i<break_masses.size(); i++)
	{
		int j;
		for (j=other_masses.size()-1; j>=0; j--)
		{
			int a;
			for (a=0; a<session_aas.size(); a++)
			{
				const int& aa = session_aas[a];
				if (aa == Ile || aa == Lys)
					continue;

				if (fabs(break_masses[i] - aa2mass[aa] - other_masses[j])<=tolerance)
				{
					cross_jump cj;
					cj.aa = aa;
					cj.idx1 = i;
					cj.idx2 = j;

					back.push_back(cj);
				}
				else
				{
					if (allow_double)
					{
						int b;
						for (b=0; b<session_aas.size(); b++)
						{
							const int& bb = session_aas[b];
							if (bb == Ile || bb == Lys)
								continue;

							if (fabs(break_masses[i] - aa2mass[aa] - aa2mass[bb]  - other_masses[j])<=tolerance)
							{
								cross_jump cj;
								cj.aa = aa;
								cj.aa2 = bb;
								cj.idx1 = i;
								cj.idx2 = j;
								back.push_back(cj);
							}				
						}
					}
				}
			}
		}
	}

//	if (over.size()==0 || back.size()==0)
//		return;

	if (print && over.size()>0 && back.size() >0)
	{
		int i;
		for (i=0; i<over.size(); i++)
		{
			string aa_label = aa2label[over[i].aa] ;
			mass_t aa_mass = aa2mass[over[i].aa];
			if (over[i].aa2>=0)
			{
				aa_label += aa2label[over[i].aa2];
				aa_mass += aa2mass[over[i].aa2];
			}
			cout << break_masses[over[i].idx1] << " -> " << aa_label << " (" << aa_mass
						 << ") " << " -> " << other_masses[over[i].idx2] << endl;
		}

		for (i=0; i<back.size(); i++)
		{
			string aa_label = aa2label[back[i].aa] ;
			mass_t aa_mass = aa2mass[back[i].aa];
			if (back[i].aa2>=0)
			{
				aa_label += aa2label[back[i].aa2];
				aa_mass += aa2mass[back[i].aa2];
			}
			cout << break_masses[back[i].idx1] << " <- " << aa_label << " (" <<
						aa_mass << ") " << " <- " << other_masses[back[i].idx2] << endl;
		}
	}

	// create graph
	vector<node_mass> nodes;
	vector< vector<bool> > edges;
	vector< int > break_idxs, other_idxs;

	// fill in nodes

	
	for (i=0; i<break_masses.size(); i++)
	{
		node_mass nm;
		nm.idx=i;
		nm.mass = break_masses[i];
		nm.source=0;
		nodes.push_back(nm);
	}

	for (i=0; i<other_masses.size(); i++)
	{
		node_mass nm;
		nm.idx=i;
		nm.mass = other_masses[i];
		nm.source=1;
		nodes.push_back(nm);
	}

	sort(nodes.begin(),nodes.end());
	int num_nodes = nodes.size();

	// fill in edges
	edges.resize(nodes.size());
	for (i=0; i<nodes.size(); i++)
		edges[i].resize(nodes.size(),false);

	// make translation vectors
	break_idxs.resize(length,-1);
	other_idxs.resize(length,-1);
	for (i=0; i<nodes.size(); i++)
		if (nodes[i].source == 0)
		{
			break_idxs[nodes[i].idx]=i;
		}
		else
			other_idxs[nodes[i].idx]=i;

	
	// fill in regular edges
	for (i=0; i<length-1; i++)
	{
		edges[break_idxs[i]][break_idxs[i+1]]=true;
		edges[other_idxs[i+1]][other_idxs[i]]=true;
	}

	// fill in crossover edges
	for (i=0; i<over.size(); i++)
		edges[break_idxs[over[i].idx1]][other_idxs[over[i].idx2]]=true;

	for (i=0; i<back.size(); i++)
		edges[other_idxs[back[i].idx2]][break_idxs[back[i].idx1]]=true;

	// fill forbidden pairs
	vector<int> forbidden_pairs;
	forbidden_pairs.resize(nodes.size(),-1);
	for (i=0; i<break_idxs.size(); i++)
	{
		int n1=break_idxs[i];
		int n2=other_idxs[i];
		forbidden_pairs[n1]=n2;
		forbidden_pairs[n2]=n1;
	}


	// print graph

	if (print)
	{
		for (i=0; i<nodes.size(); i++)
			cout << i << " " << nodes[i].mass << " " << nodes[i].idx <<" (" << nodes[i].source << ")" << endl;

		cout << endl;
		cout << "   ";
		for (i=0; i<nodes.size(); i++)
			cout << i % 10;
		cout << endl;

		for (i=0; i<nodes.size(); i++)
		{
			int j;
			cout << setw(2) << left << i << " ";
			for (j=0; j<nodes.size(); j++)
				if (edges[i][j])
				{
					cout << "+";
				}
				else
					cout << " ";
			cout << endl;
		}
		cout << endl;
	}

	// find number of paths
	int num_paths_found=0;
	vector<bool> used;
	vector<int> out_idxs, prev_idxs;
	out_idxs.resize(num_nodes);
	for (i=0; i<num_nodes; i++)
		out_idxs[i]=i;

	used.resize(num_nodes,false);
	prev_idxs.resize(num_nodes,-1);

	int p=0;
	used[0]=true;
	while (1)
	{
		while (out_idxs[p]<num_nodes)
		{
			if (p == num_nodes -2 )
			{
				// check for forbidden
				
				int i=-1;
				if (check_forbidden)
					for (i=0; i<num_nodes; i++)
						if (used[i] && used[forbidden_pairs[i]])
							break;

				if (!check_forbidden || i == num_nodes)
				{
					if (print)
					{
						int j;
						for (j=0; j<num_nodes; j++)
							if (used[j])
								cout << j << " ";
						cout << endl;
					}

					num_paths_found++;

				//	if (num_paths_found>1)
				//		return num_paths_found;
				}

				break;
			}

			if (! edges[p][out_idxs[p]])
			{
				out_idxs[p]++;
				continue;
			}


			// take this edge
			int old_p=p;
			p=out_idxs[p];
			prev_idxs[p]=old_p;
			used[p]=true;
			out_idxs[old_p]++;
		}

		// back track
		out_idxs[p]=p;
		used[p]=false;
		p=prev_idxs[p];

		if (p<0)
			break;
	}
	
	if (print)
		cout << "NUM PATHS = " << num_paths_found << endl;

	return num_paths_found;
}


void random_check_homemorphic(Config *config, int num_peptides, int peptide_length)
{
	const vector<mass_t>& aa2mass = config->get_aa2mass();

	int i,l;

	for (l=5; l<=peptide_length; l++)
	{
		int sum=0;

		for (i=0; i<num_peptides; i++)
		{
			Peptide p;
			p.generate_random_peptide(config,l);
			
		//	cout << p.as_string(config) << endl;
			
		
			if (get_num_paths(config,p,config->get_tolerance(),false,false,false)>1)
				sum++;
		//	cout << endl;
		}

		cout << l << " " << setprecision(4) << (double) sum/(double)num_peptides << endl;
	}
}


void homeomorphic_exp1(Config *config, int num_peptides)
{
	int i,peptide_length = 10;
	mass_t tol;
	vector<mass_t> tolerances;
	for (tol=0; tol<=0.101; tol += 0.01)
		tolerances.push_back(tol);
	tolerances.push_back(0.25);
	tolerances.push_back(0.5);
	tolerances.push_back(0.75);

	int t;
	for (t=0; t<tolerances.size(); t++)
	{
		tol = tolerances[t];
		double sum_tf=0, sum_tt=0, sum_ff=0, sum_ft=0;
		double num_tf=0, num_tt=0, num_ff=0, num_ft=0;
		for (i=0; i<num_peptides; i++)
		{
			Peptide p;
			p.generate_random_peptide(config,peptide_length);
		
			double tf = get_num_paths(config,p,tol,true,false,false);
			double tt = get_num_paths(config,p,tol,true,true,false);
		//	double ff = get_num_paths(config,p,tol,false,false,false);
		//	double ft = get_num_paths(config,p,tol,false,true,false);

			if (tf>1)
			{
				num_tf++;
				sum_tf+=tf;
			}

			if (tt>1)
			{
				num_tt++;
				sum_tt+=tt;
			}

		//	if (ff>1)
		//	{
		//		num_ff++;
		//		sum_ff+=ff;
		//	}

		//	if (ft>1)
		//	{
		//		num_ft++;
		//		sum_ft+=ft;
		//	}
		}
		sum_tf /= num_tf;
		sum_tt /= num_tt;
	//	sum_ff /= num_ff;
	//	sum_ft /= num_ft;

		num_tf /= num_peptides;
		num_tt /= num_peptides;
	//	num_ff /= num_peptides;
	//	num_ft /= num_peptides;

	//	cout << setw(5) << left << setprecision(4) << tol << " " << 
	//		 setw(5) << left << num_tf << " " << setw(5) << left << sum_tf <<
	//		" " << setw(5) << left << num_ff << " " << setw(5) << left << sum_ff << " | ";

		cout << setw(5) << left << setprecision(4) << tol << " " << 
			 setw(5) << left << num_tf << " " << setw(5) << left << sum_tf <<
			" " << setw(5) << left << num_tt << " " << setw(5) << left << sum_tt << endl;
	}
}



void homeomorphic_exp2(Config *config, int num_peptides)
{
	int i,a;
	mass_t tol=0.5;
	vector<int> peptide_lengths;
	for (i=5; i<=32; i+=3)
		peptide_lengths.push_back(i);

//	peptide_lengths.push_back(8);
//	peptide_lengths.push_back(17);
//	peptide_lengths.push_back(26);
	
	tol = 0.5;
	for (a=0; a<peptide_lengths.size(); a++)
	{
		double sum_tf=0, sum_tt=0, sum_ff=0, sum_ft=0;
		double num_tf=0, num_tt=0, num_ff=0, num_ft=0;
		for (i=0; i<num_peptides; i++)
		{
			Peptide p;
			p.generate_random_peptide(config,peptide_lengths[a]);
		
		//	double tf=0;
		//	double ff=0;
			double tf = get_num_paths(config,p,tol,true,false,false);
			double tt = get_num_paths(config,p,tol,true,true,false);
		//	double ff = get_num_paths(config,p,tol,false,false,false);
		//	double ft = get_num_paths(config,p,tol,false,true,false);

		
			if (tf>1)
			{
				num_tf++;
				sum_tf+=tf;
			}

			if (tt>1)
			{
				num_tt++;
				sum_tt+=tt;
			}

		//	if (ff>1)
		//	{
		//		num_ff++;
		//		sum_ff+=ff;
		//	}

		//	if (ft>1)
		//	{
		//		num_ft++;
		//		sum_ft+=ft;
		//	}
		}
		sum_tf /= num_tf;
		sum_tt /= num_tt;
	//	sum_ff /= num_ff;
	//	sum_ft /= num_ft;

		num_tf /= num_peptides;
		num_tt /= num_peptides;
	//	num_ff /= num_peptides;
	//	num_ft /= num_peptides;

	//	cout << setw(5) << left << peptide_lengths[a] << " " << 
	//		 setw(5) << left << num_tf << " " << setw(5) << left << sum_tf <<
	//		" " << setw(5) << left << num_ff << " " << setw(5) << left << sum_ff << " | ";

		cout << setw(5) << left << peptide_lengths[a] << " " << 
			 setw(5) << left << num_tf << " " << setw(5) << left << sum_tf <<
			" " << setw(5) << left << num_tt << " " << setw(5) << left << sum_tt << endl;

	//	break;
	}
}



void homeomorphic_exp3(Config *config, int num_peptides)
{
	int i,peptide_length = 10;
	mass_t tol=0.1;
	
	for (i=0; i<num_peptides; i++)
	{
		Peptide p;
		p.generate_random_peptide(config,6);
		

		if (get_num_paths(config,p,tol,true,false,false)>1)
		{
			cout << p.as_string(config) << endl;
			get_num_paths(config,p,tol,true,false,true);
			cout << endl;
		}
		
	}
}




void read_fasta(char *file_name, int **int_array, int *total_size, Config *con)
{
	char buff[1024];
	int file_size;
	int seq_p=0;
	int *arr;

	const vector<int>& char2aa = con->get_char2aa();


	ifstream file (file_name, ios::in|ios::ate);
	if (file.is_open())
	{
		file_size = file.tellg();
		file.seekg (0, ios::beg);
	}
	else
	{
		cout << "Error: reading!"<< file_name << endl;
		exit(1);
	}



	// allocate sequence memory
	arr = new int[file_size];
	if (int_array == NULL)
	{
		printf("Couldn't allocate memory for sequences!\n");
		exit(1);
	}

	// add sequence terminatng symbol -1
	// before first sequence

	file.getline(buff,1024);
	while(1)
	{
		if (buff[0] != '>')
		{
			file.getline(buff,1024);
			if (! file.good())
				break;

			continue;
		}

		// push the protein name
	

		// read protein sequence
		while ( file.getline(buff,1024) )
		{
			if (buff[0] == '>' || ! file.good())
				break;

			int i,len=strlen(buff);

			for (i=0; i<len; i++)
			{
				if (buff[i]>= 'A' && buff[i]<'Z')
					arr[seq_p]=char2aa[buff[i]];
				if (arr[seq_p]==Gln)
					arr[seq_p]=Lys;
				if (arr[seq_p]==Ile)
					arr[seq_p]=Leu;

				if (arr[seq_p]<Ala || arr[seq_p]>Val)
					seq_p--;

				seq_p++;
			}
		}
	}

	*total_size = seq_p;
	*int_array = arr;

//	int i;
//	for (i=0; i<100; i++)
//		cout << arr[i] << " ";
//	cout << endl;
	cout << "Read fasta file:   " << file_name << endl;
	cout << "Total amino acids: " << seq_p <<endl;
}



struct hpeptide {

	void calc_pep_masses_with_missing_cleavages(int *_addr, int _len, 
								Config *config, int num_miss);

	void calc_pep_masses(int *_addr, int _len, Config *config)
	{
		const vector<mass_t>& aa2mass = config->get_aa2mass();
		addr=_addr;
		len=_len;
		
		const int max_len = (len*2)+2;
		masses.resize(max_len);
		
		int i;

		masses[0]=0;
		for (i=1; i<=len; i++)
			masses[i]=masses[i-1]+aa2mass[*_addr++];
		
		masses[i++]=MASS_H2O;
		for ( ; i< max_len; i++)
			masses[i]=masses[i-1]+aa2mass[*(--_addr)];

		sort(masses.begin(),masses.end());

	/*	cout << len << " : " ;
		for (i=0; i<masses.size(); i++)
			cout << masses[i] << " ";
		cout<<endl;*/
	}


	int find_delta(const hpeptide& other, mass_t tolerance) const
	{
		int count =0;
		int i,j=0;
	
		for (i=0; i<masses.size(); i++)
		{
			mass_t min_mass = masses[i]-tolerance;
			mass_t max_mass = masses[i]+tolerance;

			while (j<other.masses.size() && other.masses[j]<min_mass)
				j++;

			if (j<other.masses.size() && other.masses[j] <= max_mass)
				count++;
		}

		return (masses.size() - count)/2;
	}

	void print(Config *config)
	{
		const vector<string>& aa2label = config->get_aa2label();
		int i;
		for (i=0; i<len; i++)
			cout << aa2label[addr[i]];

		cout << "  " << masses.size() << ": " ;
		for (i=0; i<masses.size(); i++)
			cout << " " << masses[i];
		cout << endl;
	}

	int *addr;
	int len;
	vector<mass_t> masses;
};


// removes a certain number of cleavages at random
void hpeptide::calc_pep_masses_with_missing_cleavages(int *_addr, int _len, 
								Config *config, int num_miss)
{
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	addr=_addr;
	len=_len;
		
	const int max_len = (len*2)+2;
	masses.clear();
	masses.reserve(max_len-2*num_miss);

	if (len-num_miss<2)
	{
		calc_pep_masses_with_missing_cleavages(_addr,_len,config,num_miss-1);
		return;
	}

	vector<int> skip_flags; // mask which cleavages should be skipped
	skip_flags.resize(len+1,0);

	int i;

	if (num_miss>0)
	{
	//	cout << "miss ";
		for (i=0; i<num_miss; i++)
		{
			int idx=(int)(my_random()*(len-1))+1;

			if (!skip_flags[idx])
			{
				skip_flags[idx]=1;
		//		cout << " " << idx ;
			}
			else
				idx--;
		}
	//	cout << endl;
	}

	masses.push_back(0);
	mass_t previous=0;

	for (i=1; i<=len; i++)
	{
		previous+=aa2mass[*_addr++];
		if (skip_flags[i])
			continue;
		masses.push_back(previous);
	}
		
	previous+=MASS_H2O;

	const int size=masses.size();
	for (i=0; i<size; i++)
		masses.push_back(previous-masses[i]);

	sort(masses.begin(),masses.end());
}



int homeomorphic_levels(Config *config, int *int_array, int max_length, mass_t min_mass,
						 mass_t max_mass, int num_missing_cleavages,  char *res_file, vector<int>& total_counts,
						 vector< vector<int> >& md_counts)
{
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	vector<hpeptide> peps;
	int max_idx = max_length - 40;
	int i;

//	ofstream fs(res_file,ios::out|ios::app);

	peps.reserve(300000);
	peps.clear();

	total_counts.resize(100,0);
	md_counts.resize(100);
	for (i=0; i<md_counts.size(); i++)
		md_counts[i].resize(100,0);

	// simple fill all peptides

	for (i=0; i<max_idx; i++)
	{
		int *org_addr = int_array+i;
		int *addr= org_addr;
		mass_t m=0;

		while (m<min_mass)
			m+=aa2mass[*addr++];

		if (m>max_mass)
			continue;

		hpeptide pep;
		pep.calc_pep_masses(org_addr,addr-org_addr,config);
		peps.push_back(pep);
	}

	if (peps.size()<2)
		return 0;

	for (i=0; i<10000; i++)
	{
		int min_dis=999999;
		int j;

		int p_idx= (int)(my_random()*peps.size());

		hpeptide p = peps[p_idx];
		p.calc_pep_masses_with_missing_cleavages(p.addr,p.len,config,num_missing_cleavages);

	//	p.print(config);
	//	peps[p_idx].print(config);

		vector<int> counts;
		counts.resize(peps[p_idx].len,0);


		for (j=0; j<peps.size(); j++)
		{
			if (p_idx==j)
				continue;

			int dis = p.find_delta(peps[j],config->get_tolerance());
			counts[dis]++;

			if (dis < min_dis)
			{
				// check if this is the same peptide
				if (dis == 0)
				{
					if (peps[p_idx].len == peps[j].len)
					{
						int k;
						for (k=0; k<peps[p_idx].len; k++)
							if (peps[p_idx].addr[k] != peps[j].addr[k])
								break;

						if (k==peps[p_idx].len)
						{
							counts[0]--;
							continue;
						}
					}
				}

				min_dis = dis;

				if (dis<1)
				{

				//	if (peps[j].find_delta(peps[p_idx],0.5)<1)
				//	{
						cout << "d("<<p_idx<<","<<j<<")= "<<dis << endl;
						p.print(config);
						peps[j].print(config);
						cout << endl;
				//	}
				}
			}

			total_counts[peps[p_idx].len]++;
			md_counts[peps[p_idx].len][min_dis]++;
		}

	//	fs << peps[p_idx].len << " " << peps.size() << " ";
	//	for (j=0; j<counts.size(); j++)
	//		fs << " " << counts[j];
	//	fs << endl;


	//	cout << endl;
	}

//	fs.close();
	return 1;

}


void full_exp(Config *config, int num_missing, int *int_array, int max_length, char *res_file)
{
	mass_t mass = 700;
	const int max_aa = 50;
	int i;

	vector<int> total_counts; // aa
	vector< vector<int> > min_dis_counts; // aa, dis

	total_counts.resize(max_aa,0);
	min_dis_counts.resize(max_aa);
	for (i=0; i<min_dis_counts.size(); i++)
		min_dis_counts[i].resize(max_aa,0);

	mass_t tol = config->get_tolerance() * 2;
	while (mass<2200)
	{
		vector<int> counts;
		vector< vector<int> > md_counts;
		
		mass += (my_random()*99);
		if (! homeomorphic_levels(config, int_array, max_length , mass-tol, mass+tol, num_missing, res_file,
							counts, md_counts) )
			continue;

		int j;
		for (i=0; i<max_aa; i++)
		{
			total_counts[i]+=counts[i];
			for (j=0; j<max_aa; j++)
				min_dis_counts[i][j]+=md_counts[i][j];
		}
	}

	// print results
	ofstream fs(res_file,ios::out|ios::app);
	for (i=4; i<30; i++)
	{
		if (i != 7 && i != 14 && i != 21)
			continue;

		if (i> 21 && total_counts[i]==0)
			break;

		fs << i << " ";
		int j;
		for (j=0; j<=i; j++)
			fs << setw(5) << (double)min_dis_counts[i][j]/(double)total_counts[i] << "\t";
		fs << endl;
	}

	fs.close();
}
