#include "projectionutils.h"

float spectrum_similarity(psmPtr psm1, string annotation1,  psmPtr psm2, string annotation2,
                          int peptideLength, MS2ScoringModel &model,
                          vector<string> &ionsToExtract, string allIons){
    
    psmPtr psm1temp(new PeptideSpectrumMatch);
    psmPtr psm2temp(new PeptideSpectrumMatch);
    
    Spectrum temp_spectrum1 = *(psm1->m_spectrum);
    Spectrum temp_spectrum2 = *(psm2->m_spectrum);
    
    Spectrum * old_spectrum1 = psm1->m_spectrum;
    Spectrum * old_spectrum2 = psm2->m_spectrum;
    
    psm1temp->m_spectrum = &temp_spectrum1;
    psm2temp->m_spectrum = &temp_spectrum2;
    psm1temp->m_annotation = psm1->m_annotation;
    psm2temp->m_annotation = psm2->m_annotation;
    
    //Making max 1000
    preprocess_spectrum_intensities_max_intensity(psm1temp->m_spectrum, 1000.f);
    preprocess_spectrum_intensities_max_intensity(psm2temp->m_spectrum, 1000.f);
    
    //Applying SQRT
    preprocess_spectrum_intensities(psm1temp->m_spectrum, 0, 1);
    preprocess_spectrum_intensities(psm2temp->m_spectrum, 0, 1);
    
    psm1temp->annotate(psm1temp->m_annotation,allIons,model,0,0,.45);
    psm2temp->annotate(psm2temp->m_annotation,allIons,model,0,0,.45);
    
    vector<float> ion_intensities1;
    vector<float> ion_intensities2;
    
    extractIons(psm1temp,peptideLength,model,ionsToExtract,ion_intensities1, 0, 0);
    extractIons(psm2temp,peptideLength,model,ionsToExtract,ion_intensities2, 0, 0);

    normalize_extracted_ions(ion_intensities1);
    normalize_extracted_ions(ion_intensities2);
    
    float cosine_val = cosine(ion_intensities1, ion_intensities2);
    
    psm1->m_spectrum = old_spectrum1;
    psm2->m_spectrum = old_spectrum2;
    
    return cosine_val;
}



float spectrum_similarity(psmPtr psm1, psmPtr psm2, int peptideLength, MS2ScoringModel &model,
        vector<string> &ionsToExtract, string allIons){

    psmPtr psm1temp(new PeptideSpectrumMatch);
    psmPtr psm2temp(new PeptideSpectrumMatch);
    
    Spectrum temp_spectrum1 = *(psm1->m_spectrum);
    Spectrum temp_spectrum2 = *(psm2->m_spectrum);
    
    Spectrum * old_spectrum1 = psm1->m_spectrum;
    Spectrum * old_spectrum2 = psm2->m_spectrum;
    
    psm1temp->m_spectrum = &temp_spectrum1;
    psm2temp->m_spectrum = &temp_spectrum2;
    psm1temp->m_annotation = psm1->m_annotation;
    psm2temp->m_annotation = psm2->m_annotation;
    
    //Making max 1000
    preprocess_spectrum_intensities_max_intensity(psm1temp->m_spectrum, 1000.f);
    preprocess_spectrum_intensities_max_intensity(psm2temp->m_spectrum, 1000.f);
    
    //Applying SQRT
    preprocess_spectrum_intensities(psm1temp->m_spectrum, 0, 1);
    preprocess_spectrum_intensities(psm2temp->m_spectrum, 0, 1);
    
    psm1temp->annotate(psm1temp->m_annotation,allIons,model,0,0,.45);
    psm2temp->annotate(psm2temp->m_annotation,allIons,model,0,0,.45);
    
    vector<float> ion_intensities1;
    vector<float> ion_intensities2;
    
    extractIons(psm1temp,peptideLength,model,ionsToExtract,ion_intensities1, 0, 0);
    extractIons(psm2temp,peptideLength,model,ionsToExtract,ion_intensities2, 0, 0);

    normalize_extracted_ions(ion_intensities1);
    normalize_extracted_ions(ion_intensities2);
    
    float cosine_val = cosine(ion_intensities1, ion_intensities2);
    
    psm1->m_spectrum = old_spectrum1;
    psm2->m_spectrum = old_spectrum2;
    
    return cosine_val;
}


float spectrum_similarity(psmPtr psm1, psmPtr psm2, int peptideLength, MS2ScoringModel &model,
        vector<string> &ionsToExtract, string allIons, vector<vector<float> > average_ions){
    
    psmPtr psm1temp(new PeptideSpectrumMatch);
    psmPtr psm2temp(new PeptideSpectrumMatch);
    
    Spectrum temp_spectrum1 = *(psm1->m_spectrum);
    Spectrum temp_spectrum2 = *(psm2->m_spectrum);
    
    Spectrum * old_spectrum1 = psm1->m_spectrum;
    Spectrum * old_spectrum2 = psm2->m_spectrum;
    
    psm1temp->m_spectrum = &temp_spectrum1;
    psm2temp->m_spectrum = &temp_spectrum2;
    psm1temp->m_annotation = psm1->m_annotation;
    psm2temp->m_annotation = psm2->m_annotation;
    
    //Making max 1000
    preprocess_spectrum_intensities_max_intensity(psm1temp->m_spectrum, 1000.f);
    preprocess_spectrum_intensities_max_intensity(psm2temp->m_spectrum, 1000.f);
    
    //Applying SQRT
    preprocess_spectrum_intensities(psm1temp->m_spectrum, 0, 1);
    preprocess_spectrum_intensities(psm2temp->m_spectrum, 0, 1);
    
    psm1temp->annotate(psm1temp->m_annotation,allIons,model,0,0,.45);
    psm2temp->annotate(psm2temp->m_annotation,allIons,model,0,0,.45);
    
    vector<float> ion_intensities1;
    vector<float> ion_intensities2;
    
    extractIons(psm1temp,peptideLength,model,ionsToExtract,ion_intensities1, 0, 0);
    extractIons(psm2temp,peptideLength,model,ionsToExtract,ion_intensities2, 0, 0);

    normalize_extracted_ions(ion_intensities1);
    normalize_extracted_ions(ion_intensities2);
    
    //Subtracting off average
    for(int i = 0; i < ion_intensities1.size(); i++){
        ion_intensities1[i] = ion_intensities1[i] - average_ions[peptideLength][i];
        ion_intensities2[i] = ion_intensities2[i] - average_ions[peptideLength][i];
    }
    
    normalize_extracted_ions(ion_intensities1);
    normalize_extracted_ions(ion_intensities2);
    
    float cosine_val = cosine(ion_intensities1, ion_intensities2);
    
    psm1->m_spectrum = old_spectrum1;
    psm2->m_spectrum = old_spectrum2;
    
    return cosine_val;
}


float full_spectrum_similarity(Spectrum spec1, Spectrum spec2){
    float max_mass = 0.f;
    
    Spectrum spec1_temp = spec1;
    Spectrum spec2_temp = spec2;
    spec1_temp.setResolution(1.0005, 1);
    spec2_temp.setResolution(1.0005, 1);

    preprocess_spectrum_intensities(&spec1_temp, 0,1);
    preprocess_spectrum_intensities(&spec2_temp, 0,1);
    
    //Finding the maximum mass in a spectrum
    for(int spec1_idx = 0; spec1_idx < spec1_temp.size(); spec1_idx++){
        if(max_mass < spec1_temp[spec1_idx][0])
            max_mass = spec1_temp[spec1_idx][0];
    }
    
    for(int spec2_idx = 0; spec2_idx < spec2_temp.size(); spec2_idx++){
        if(max_mass < spec2_temp[spec2_idx][0])
            max_mass = spec2_temp[spec2_idx][0];
    }
    
    
    vector<float> spec1_peak_vector( (int)max_mass + 1);
    vector<float> spec2_peak_vector( (int)max_mass + 1);
    for(int i = 0; i < spec1_temp.size(); i++)
        spec1_peak_vector[(int)spec1_temp[i][0]] = spec1_temp[i][1];
    for(int i = 0; i < spec2_temp.size(); i++)
        spec2_peak_vector[(int)spec2_temp[i][0]] = spec2_temp[i][1];
    
    normalize_extracted_ions(spec1_peak_vector);
    normalize_extracted_ions(spec2_peak_vector);
    
    //float cumu = 0.f;
    //for(int i = 0; i < spec1_peak_vector.size(); i++){
    //    if(spec1_peak_vector[i] > 0 || spec2_peak_vector[i] > 0){
    //        cumu += spec1_peak_vector[i] * spec2_peak_vector[i];
    //        cout<<i<<"\t"<<spec1_peak_vector[i]<<"\t"<<spec2_peak_vector[i]<<"\t"<<cumu<<endl;
    //    }
    //}
    
    float cosine_val = cosine(spec1_peak_vector, spec2_peak_vector);
    
    return cosine_val;
}

//Calculating amino acid mass
float getMass(string peptide, map<char, float> amino_acid_map){
    float total_mass = 0.f;
    for(int i = 0; i < peptide.length(); i++){
        char amino_acid = tolower(peptide[i]);
        map<char,float>::iterator it;
        it = amino_acid_map.find(amino_acid);
        if(it != amino_acid_map.end()){
            total_mass += amino_acid_map[amino_acid];
        }
    }
    return total_mass;
}

float getSubstitutionCosineDepression(char a, char b, map<string, float> amino_acid_transformation_map){
	a = toupper(a);
	b = toupper(b);
	
	string lookup1 = "";
	lookup1 += a;
	lookup1 += b;
	string lookup2 = "";
	lookup2 += b;
	lookup2 += a;
	
	float lookup1_score = 0.f;
	float lookup2_score = 0.f;
	
	if(amino_acid_transformation_map.find(lookup1) != amino_acid_transformation_map.end()){
		lookup1_score = amino_acid_transformation_map[lookup1];
	}
	
	if(amino_acid_transformation_map.find(lookup2) != amino_acid_transformation_map.end()){
		lookup2_score = amino_acid_transformation_map[lookup2];
	}
	
	return max(lookup2_score, lookup1_score);
}

//Generate lookup table for amino acid substitutions and cosine depression
map<string, float> getAminoAcidSubstitutionLookup(){
    map<string, float> amino_acid_transform_map;
	amino_acid_transform_map["WT"] = 0.912276;
	amino_acid_transform_map["AS"] = 0.908659;
	amino_acid_transform_map["YA"] = 0.900532;
	amino_acid_transform_map["DI"] = 0.920199;
	amino_acid_transform_map["DL"] = 0.919106;
	amino_acid_transform_map["FE"] = 0.936028;
	amino_acid_transform_map["EV"] = 0.922741;
	amino_acid_transform_map["SG"] = 0.886207;
	amino_acid_transform_map["LN"] = 0.90118;
	amino_acid_transform_map["QV"] = 0.91915;
	return amino_acid_transform_map;
}

//Generate lookup table for amino acid masses
map<char, float> getAminoAcidLookup(){
    map<char, float> amino_acid_map;
    amino_acid_map['g'] = 57.02146;
    amino_acid_map['a'] = 71.03711f;
    amino_acid_map['s'] = 87.03203f;
    amino_acid_map['p'] = 97.05276f;
    amino_acid_map['v'] = 99.06841f;
    amino_acid_map['t'] = 101.04768f;
    amino_acid_map['c'] = 103.00919f;
    amino_acid_map['l'] = 113.08406f;
    amino_acid_map['i'] = 113.08406f;
    amino_acid_map['n'] = 114.04293f;
    amino_acid_map['d'] = 115.02694f;
    amino_acid_map['q'] = 128.05858f;
    amino_acid_map['k'] = 128.09496f;
    amino_acid_map['e'] = 129.04259f;
    amino_acid_map['m'] = 131.04049f;
    amino_acid_map['h'] = 137.05891f;
    amino_acid_map['f'] = 147.06841f;
    amino_acid_map['r'] = 156.10111f;
    amino_acid_map['y'] = 163.06333f;
    amino_acid_map['w'] = 186.07931f;
    return amino_acid_map;
}

string stripString(string peptide){
	string stripped_peptide = peptide;
	//Removing * and . at ends
	if(stripped_peptide[0] == '*'){
		stripped_peptide.erase(stripped_peptide.begin() + 0);
	}
	
	if(stripped_peptide[stripped_peptide.length()-1] == '*'){
		stripped_peptide.erase(stripped_peptide.begin() + stripped_peptide.length()-1);
	}
	
	if(stripped_peptide[0] == '.'){
		stripped_peptide.erase(stripped_peptide.begin() + 0);
	}
	
	if(stripped_peptide[stripped_peptide.length()-1] == '.'){
		stripped_peptide.erase(stripped_peptide.begin() + stripped_peptide.length()-1);
	}
	
	return stripped_peptide;
}

string getAdditionalMassString(int parent_mass_difference){
	//inserting the additional mass
	stringstream mass_addition_stream;
	if(parent_mass_difference > 0){
		mass_addition_stream<<"["<<parent_mass_difference<<"]";
	}
	else if(parent_mass_difference < 0){
		mass_addition_stream<<"["<<parent_mass_difference<<"]";
	}
	//otherwise equal, so do nothing
	
	return mass_addition_stream.str();
}


string cleanAnnotationEnds(string annotation){
    if(annotation[0] != '*'){
            annotation[0] = '*';
            //cout<<new_peptides[i]<<"\t"<<annotation<<endl;
    }
    if(annotation[annotation.length() - 1] != '*'){
        annotation[annotation.length()-1] = '*';
        //cout<<new_peptides[i]<<"\t"<<annotation<<endl;
    }
    return annotation;
}


string create_annotation_ends(string annotation){
    int firstdotloc = annotation.find_first_of(".");
    int lastdotloc = annotation.find_last_of(".");
    if(lastdotloc != (annotation.length() - 2)) annotation.insert(annotation.length(), ".*");
    if(firstdotloc != 1) annotation.insert(0, "*.");
    
    return annotation;    
}

string remove_annotation_ends(string annotation){
    int firstdotloc = annotation.find_first_of(".");
    int lastdotloc = annotation.find_last_of(".");
    if(lastdotloc == (annotation.length() - 2)) annotation = annotation.substr(0, annotation.length() - 2);
    if(firstdotloc == 1) annotation = annotation.substr(2, annotation.length() - 2);
    
    return annotation;
}



int getDifferenceIndex(string in1, string in2){
    if(in1.length() != in2.length())
        return -1;
    int difference = 0;
    for(int i = 0; i < in1.length(); i++){
        if(in1[i] != in2[i]){
            return i;
        }
    }
    return -1;
}

int getStringDifference(string in1, string in2){
    if(in1.length() != in2.length())
        return -1;
    int difference = 0;
    for(int i = 0; i < in1.length(); i++){
        if(in1[i] != in2[i]){
            difference++;
        }
    }
    return difference;
}

// Computes to cosine between two normalized vectors of the same size
float cosine(vector<float> &u, vector<float> &v) {
    if(u.size()!=v.size()) return 0;
    float c=0;
    for(unsigned int i=0; i<u.size(); i++) c+=u[i]*v[i];
    return c;
}

void preprocess_spectrum_intensities(Spectrum * s, int do_early_normalize, int do_filter){
    if(do_filter == 1){
        for(int peakIdx = 0; peakIdx < s->size(); peakIdx++){
            (*s)[peakIdx][1] =  sqrt((*s)[peakIdx][1]);
        }
    }
    if(do_filter == 2){
        for(int peakIdx = 0; peakIdx < s->size(); peakIdx++){
            (*s)[peakIdx][1] =  log((*s)[peakIdx][1]);
        }
    }
}


void preprocess_spectrum_intensities_max_intensity(Spectrum * s, float max){
    float current_max_intesity = 0.f;
    for(int i = 0; i < s->size(); i++){
        if(current_max_intesity < (*s)[i][1]) current_max_intesity = (*s)[i][1];
    }
    
    for(int i = 0; i < s->size(); i++){
        (*s)[i][1] = (*s)[i][1] * max / current_max_intesity;
    }
}

void normalize_extracted_ions(vector<float> &v) {
    //if(do_early_normalize == 1)
    //    return;
    
    float n=0;
    for(unsigned int i=0; i<v.size(); i++) n+=v[i]*v[i];
        if(n>0) {
        n=sqrt(n);
        for(unsigned int i=0; i<v.size(); i++) v[i]/=n;
    }
}

void extractIons(psmPtr psm, int peptideLength, MS2ScoringModel &model,
        vector<string> &ionsToExtract, vector<float> &ions, int do_early_normalize, int do_filter) {
	ions.resize(peptideLength*ionsToExtract.size());
	for(unsigned int i=0; i<ions.size(); i++) ions[i]=0;
	
	if(psm->m_peakAnnotations.size()!=psm->m_spectrum->size()) return;
	
	//Normalize before extraction, so that we take into account noise
    //normalize_spectrum(s, do_early_normalize, do_filter);

	// Select target ions
	unsigned int baseIdx=0;
	for(unsigned int ionType=0; ionType<ionsToExtract.size(); ionType++) {
        for(unsigned int peakIdx=0; peakIdx<psm->m_peakAnnotations.size(); peakIdx++)
            if(psm->m_peakAnnotations[peakIdx].first and psm->m_peakAnnotations[peakIdx].first->name==ionsToExtract[ionType]){
				ions[baseIdx+psm->m_peakAnnotations[peakIdx].second-1]=(*psm->m_spectrum)[peakIdx][1];
			}
		baseIdx+=peptideLength;
	}
}

void extractIons(psmPtr psm, int peptideLength, MS2ScoringModel &model,
        vector<string> &ionsToExtract, vector<pair <float, float> > &ions, int do_early_normalize, int do_filter) {
	ions.resize(peptideLength*ionsToExtract.size());
	for(unsigned int i=0; i<ions.size(); i++){
		ions[i].first=0;
		ions[i].second=0;
	}
	if(psm->m_peakAnnotations.size()!=psm->m_spectrum->size()) return;

	//Normalize before extraction, so that we take into account noise
    //normalize_spectrum(s, do_early_normalize, do_filter);

	// Select target ions
	unsigned int baseIdx=0;
	for(unsigned int ionType=0; ionType<ionsToExtract.size(); ionType++) {
		for(unsigned int peakIdx=0; peakIdx<psm->m_peakAnnotations.size(); peakIdx++)
			if(psm->m_peakAnnotations[peakIdx].first and psm->m_peakAnnotations[peakIdx].first->name==ionsToExtract[ionType]){
				ions[baseIdx+psm->m_peakAnnotations[peakIdx].second-1].second=(*psm->m_spectrum)[peakIdx][1];
				ions[baseIdx+psm->m_peakAnnotations[peakIdx].second-1].first=(*psm->m_spectrum)[peakIdx][0];
			}
		baseIdx+=peptideLength;
	}
}


void norm_vector(vector<float> & input){
    float sum = 0.f;
    for(int i = 0; i < input.size(); i++){
        sum += input[i]  * input[i];
    }
    
    if(sum == 0.f) return;
    
    float euclidian_norm = sqrt(sum);
    for(int i = 0; i < input.size(); i++){
        input[i] = input[i] / euclidian_norm;
    }
}

//Comparater function for sorting
bool search_results_comparator (score_results_tuple i, score_results_tuple j){
    return tr1::get<1>(i) > tr1::get<1>(j);
}


string reverse_string(string input){
    string output = "";
    for(int i = 0; i < input.length(); i++){
        output += input[input.length() - 1 - i];
    }
    return output;
}

float full_spectrum_dotbias(Spectrum spec1, Spectrum spec2, float spec_sim){
    float max_mass = 0.f;
    
    Spectrum spec1_temp = spec1;
    Spectrum spec2_temp = spec2;
    spec1_temp.setResolution(1.0005, 1);
    spec2_temp.setResolution(1.0005, 1);

    preprocess_spectrum_intensities(&spec1_temp, 0,1);
    preprocess_spectrum_intensities(&spec2_temp, 0,1);
    
    //Finding the maximum mass in a spectrum
    for(int spec1_idx = 0; spec1_idx < spec1_temp.size(); spec1_idx++){
        if(max_mass < spec1_temp[spec1_idx][0])
            max_mass = spec1_temp[spec1_idx][0];
    }
    
    for(int spec2_idx = 0; spec2_idx < spec2_temp.size(); spec2_idx++){
        if(max_mass < spec2_temp[spec2_idx][0])
            max_mass = spec2_temp[spec2_idx][0];
    }
    
    for(int i = 0; i < spec1_temp.size(); i++){
        //cout<<spec1_temp.peakList[i][0]<<"\t"<<spec1_temp.peakList[i][1]<<"\t"<<spec1.peakList[i][0]<<"\t"<<spec1.peakList[i][1]<<endl;
    }
    
    vector<float> spec1_peak_vector( (int)max_mass + 1);
    vector<float> spec2_peak_vector( (int)max_mass + 1);
    for(int i = 0; i < spec1_temp.size(); i++)
        spec1_peak_vector[(int)spec1_temp[i][0]] = spec1_temp[i][1]*spec1_temp[i][1];
    for(int i = 0; i < spec2_temp.size(); i++)
        spec2_peak_vector[(int)spec2_temp[i][0]] = spec2_temp[i][1]*spec2_temp[i][1];
    
    
    
    normalize_extracted_ions(spec1_peak_vector);
    normalize_extracted_ions(spec2_peak_vector);
    
    float cosine_val = cosine(spec1_peak_vector, spec2_peak_vector);
    
    float sqrt_cosine_val = sqrt(cosine_val);
    
    float DB = sqrt_cosine_val/spec_sim;
    
    if(DB < 0.1)
        return 0.12;
    if(DB <= 0.4 && DB > 0.35)
        return .12;
    if(DB > 0.4 && DB <= 0.45)
        return 0.18;
    if(DB > 0.45)
        return 0.24;
    
    return 0;
}

void sorted_vector_library(vector<Spectrum *> &library_ptr, SpectralLibrary & real_library){
    for(int i = 0; i < real_library.size(); i++){
        library_ptr.push_back(&real_library[i]);
    }
    
    sort(library_ptr.begin(), library_ptr.end(), spectrum_ptr_compare);
    
}

bool spectrum_ptr_compare (Spectrum* first, Spectrum* second){ 
    return (first->parentMZ<second->parentMZ);
}

int spectrum_ptr_startend(vector<Spectrum *> &library_ptr, float parentMZ, float tolerance, int &start_idx, int &end_idx){
    if(library_ptr.size() == 0){
        start_idx = 0;
        end_idx = -1;
        return 0;
    }
    
    float looking_mass = parentMZ - tolerance;
    
    int start_range = 0;
    int end_range = library_ptr.size()-1;
    while(end_range - start_range > 1){
        int middle = ( start_range + end_range )/2;
        if(looking_mass > library_ptr[middle]->parentMZ){
            start_range = middle + 1;
            continue;
        }
        if(looking_mass < library_ptr[middle]->parentMZ){
            end_range = middle - 1;
            continue;
        }
        if(looking_mass == library_ptr[middle]->parentMZ){
            start_range = middle;
            end_range = middle;
            break;
        }
    }
    
    start_idx = start_range;
    
    looking_mass= parentMZ + tolerance;
    start_range = 0;
    end_range = library_ptr.size()-1;
    while(end_range - start_range > 1){
        int middle = ( start_range + end_range )/2;
        if(looking_mass > library_ptr[middle]->parentMZ){
            start_range = middle + 1;
            continue;
        }
        if(looking_mass < library_ptr[middle]->parentMZ){
            end_range = middle - 1;
            continue;
        }
        if(looking_mass == library_ptr[middle]->parentMZ){
            start_range = middle;
            end_range = middle;
            break;
        }
    }
    
    end_idx = start_range;
    
    return 0;
}

int getpeptideLength(string annotation){
    int length = 0;
    for(int i = 0; i < annotation.length(); i++){
        if(annotation[i] == 'A' || 
            annotation[i] == 'R' ||
            annotation[i] == 'N' ||
            annotation[i] == 'D' ||
            annotation[i] == 'C' ||
            annotation[i] == 'E' ||
            annotation[i] == 'Q' ||
            annotation[i] == 'G' ||
            annotation[i] == 'H' ||
            annotation[i] == 'I' ||
            annotation[i] == 'L' ||
            annotation[i] == 'K' ||
            annotation[i] == 'M' ||
            annotation[i] == 'F' ||
            annotation[i] == 'P' ||
            annotation[i] == 'S' ||
            annotation[i] == 'T' ||
            annotation[i] == 'W' ||
            annotation[i] == 'Y' ||
            annotation[i] == 'V' ||
            annotation[i] == '[')
            length++;    
    }
    return length;
}


