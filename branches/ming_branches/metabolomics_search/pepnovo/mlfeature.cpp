#include "mlfeature.h"
#include "mlauxfun.h"
#include "auxfun.h"





void MlFeatureSet::initializeWithDefaultFeatures(size_t numFeatures)
{
	features_.clear();
	features_.resize(numFeatures);

	size_t i;
	for (i=0; i<numFeatures; i++)
	{
		char buff[16];
		sprintf(buff,"F%d",i);
		features_[i].name_ = static_cast<string>(buff);
	}
}



void MlFeatureSet::readFeatureFile(const char* path, bool verbose)
{
	features_.clear();
	ifstream ifs(path);
	if (! ifs.is_open())
		error("couldn't open feature file for reading: ",path);

	size_t counter=0;
	char buffer[256];
	while (! ifs.eof())
	{
		ifs.getline(buffer,256);
		if (ifs.gcount()<=1)
			continue;

		if (buffer[0] == '#')
		{
			size_t n=0;
			if (sscanf(buffer,"#NUM FEAUTRES %d",&n) == 1)
				features_.resize(n);
			continue;
		}

		MlFeature mlf;
		mlf.index_ = counter++;
		
		istringstream iss(buffer);
		if (sscanf(buffer,"%d",&mlf.index_)==1) // line contains feature index
			iss >> mlf.index_;

		string unaryFlag="-";
		string conditionalFlag="+";
		iss >> mlf.name_ >> unaryFlag >> conditionalFlag;

		if (mlf.name_.length()<1)
			error("bad feature line: ",buffer);

		if (unaryFlag[0]=='-')
		{
			mlf.flagConsiderForUnary_ = false;
		}
		else if (unaryFlag[0] != '+')
			error("Illegal unary flag (should be +/-): ",unaryFlag.c_str());
		
		if (conditionalFlag[0]=='-')
		{
			mlf.flagConsiderForConditional_ = false;
		}
		else if (conditionalFlag[0] != '+')
			error("Illegal conditional flag (should be +/-): ",conditionalFlag.c_str());


		// read groups
		while (true)
		{
			string g="";
			iss >> g;
			if (g.length()>0)
			{
				// get group index
				map<string,size_t>::const_iterator it=this->groupMapping_.find(g);
				if (it != groupMapping_.end())
				{
					mlf.groups_.push_back(it->second);
					continue;
				}

				groupMapping_[g]=groupNames_.size();
				groupNames_.push_back(g);
				continue;
			}
			break;
		}
		features_.push_back(mlf);
	}
	ifs.close();

	// check features order
	size_t i;
	for (i=1; i<features_.size(); i++)
		if (features_[i].index_<= features_[i-1].index_)
			break;
	if (i<features_.size())
		sort(features_.begin(), features_.end());

	// make sure there are no duplicate indices
	for (i=1; i<features_.size(); i++)
		if (features_[i].index_ == features_[i-1].index_)
			break;
	if (i<features_.size())
		error("duplicate feature indexes detected, index ", features_[i].index_);

	// create lists of group members
	groupMembers_.resize(groupNames_.size());
	for (i=0; i<features_.size(); i++)
	{
		size_t j;
		for (j=0; j<features_[i].groups_.size(); j++)
			groupMembers_[features_[i].groups_[j]].push_back(i);
	}

	if (verbose)
		cout << "Read " << features_.size() << " features from " << path << endl;
}

void MlFeatureSet::writeFeatureFile(const char* path, bool verbose) const
{
	ofstream ofs(path);
	if (! ofs.is_open())
		error("couldn't open feature file for writing: ",path);

	ofs << "#NUM FEATURES " << features_.size() << endl;
	size_t i;
	for (i=0; i<features_.size(); i++)
	{
		char buffer[256];
		size_t len=sprintf(buffer,"%d\t%s\t %c %c\n",features_[i].getIndex(),features_[i].getName().c_str(),
			(features_[i].getFlagConsiderUnary() ? '+' : '-'), 
			(features_[i].getFlagConsiderConditional() ? '+' : '-'));

		const vector<size_t>& groups = features_[i].getGroups();
		char* pos = buffer + len;
		size_t j;
		for (j=0; j<groups.size(); j++)
			pos+=sprintf(pos,"\t%s",this->groupNames_[groups[j]].c_str());

		ofs << buffer << endl;
	}
	ofs.close();
	if (verbose)
		cout << "Wrote " << features_.size() << " features to " << path << endl;
}


