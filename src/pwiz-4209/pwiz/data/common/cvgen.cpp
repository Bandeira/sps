//
// $Id: cvgen.cpp 3540 2012-04-16 20:31:12Z pcbrefugee $
//
//
// Original author: Darren Kessner <darren@proteowizard.org>
//
// Copyright 2007 Spielberg Family Center for Applied Proteomics
//   Cedars-Sinai Medical Center, Los Angeles, California  90048
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//


#include "obo.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/filesystem/fstream.hpp"
#include "boost/regex.hpp"
#include "pwiz/utility/misc/Std.hpp"

using namespace pwiz::data;
namespace bfs = boost::filesystem;


//
// This program selectively parses OBO format controlled vocabulary files
// and generates C++ code (one hpp file and one cpp file).
//


void writeCopyright(ostream& os, const string& filename)
{
    os << "//\n"
       << "// $Id: cvgen.cpp 3540 2012-04-16 20:31:12Z pcbrefugee $" << endl
       << "//\n"
          "//\n"
          "// Darren Kessner <darren@proteowizard.org>\n"
          "//\n"
          "// Copyright 2007 Spielberg Family Center for Applied Proteomics\n"
          "//   Cedars-Sinai Medical Center, Los Angeles, California  90048\n"
          "//\n"
          "// Licensed under the Apache License, Version 2.0 (the \"License\");\n"
          "// you may not use this file except in compliance with the License.\n"
          "// You may obtain a copy of the License at\n"
          "//\n"
          "// http://www.apache.org/licenses/LICENSE-2.0\n"
          "//\n"
          "// Unless required by applicable law or agreed to in writing, software\n"
          "// distributed under the License is distributed on an \"AS IS\" BASIS,\n"
          "// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.\n"
          "// See the License for the specific language governing permissions and\n"
          "// limitations under the License.\n"
          "//\n"
          "// This file was generated by cvgen.\n"
          "//\n"
          "// Do not edit this file! Your changes will be lost next time cvgen is run -\n"
          "// see pwiz/scripts/misc/update_cv.bat for info on how that works.\n"
          "// Instead, edit cvgen.cpp itself, or the cv.inl include file if adding static\n"
          "// code or data.\n"
          "//\n\n\n";
}


string includeGuardString(const string& basename)
{
    string includeGuard = basename;
    transform(includeGuard.begin(), includeGuard.end(), includeGuard.begin(), (int(*)(int))toupper);
    return "_" + includeGuard + "_HPP_";
}


void namespaceBegin(ostream& os, const string& name)
{
    os << "namespace pwiz {\nnamespace cv {\n\n\n";
}


void namespaceEnd(ostream& os, const string& name)
{
    os << "} // namespace cv\n} // namespace pwiz\n\n\n";
}


inline char toAllowableChar(char a)
{
    return isalnum(a) ? a : '_';
}


string enumName(const string& prefix, const string& name, bool isObsolete)
{
    string result = name;
    transform(result.begin(), result.end(), result.begin(), toAllowableChar);
    result = prefix + "_" + result + (isObsolete ? "_OBSOLETE" : "");
    return result;
}


string enumName(const Term& term)
{
    return enumName(term.prefix, term.name, term.isObsolete);
}


const size_t enumBlockSize_ = 100000000;


size_t enumValue(const Term& term, size_t index)
{
    return term.id + (enumBlockSize_ * index);
}


vector< map<Term::id_type, const Term*> > termMaps;
vector< map<Term::id_type, string> > correctedEnumNameMaps;

void writeHpp(const vector<OBO>& obos, const string& basename, const bfs::path& outputDir)
{
    string filename = basename + ".hpp";
    bfs::path filenameFullPath = outputDir / filename;
    bfs::ofstream os(filenameFullPath, ios::binary);

    writeCopyright(os, filename);

    string includeGuard = includeGuardString(basename);
    os << "#ifndef " << includeGuard << endl
       << "#define " << includeGuard << "\n\n\n"
       << "#include <string>\n"
       << "#include <vector>\n"
       << "#include <map>\n"
       << "#include \"pwiz/utility/misc/Export.hpp\"\n"
       << "\n\n";

    for (vector<OBO>::const_iterator obo=obos.begin(); obo!=obos.end(); ++obo)
    {
		string filename = bfs::path(obo->filename).filename();
        os << "// [" << filename << "]\n";
        string fileDefine = bal::to_upper_copy(filename);
        transform(fileDefine.begin(), fileDefine.end(), fileDefine.begin(), toAllowableChar);
        os << "#define _" << fileDefine << "_\n";

        for (vector<string>::const_iterator it=obo->header.begin(); it!=obo->header.end(); ++it)
            os << "//   " << *it << endl;

        os << "//\n";
    }
    os << "\n\n";

    namespaceBegin(os, basename);

    os << "/// enumeration of controlled vocabulary (CV) terms, generated from OBO file(s)\n"
          "enum PWIZ_API_DECL CVID\n{\n"
          "    CVID_Unknown = -1";
    for (vector<OBO>::const_iterator obo=obos.begin(); obo!=obos.end(); ++obo)
    BOOST_FOREACH(const Term& term, obo->terms)
    {
        string& eName = correctedEnumNameMaps[obo-obos.begin()][term.id];

        os << ",\n\n"
           << "    /// " << term.name << ": " << term.def << "\n"
           << "    " << eName << " = " << enumValue(term, obo-obos.begin());

        if (obo->prefix == "MS") // add synonyms for PSI-MS only
        {
            BOOST_FOREACH(const string& synonym, term.exactSynonyms)
            {
                os << ",\n\n"
                   << "    /// " << synonym << " (" << term.name << "): " << term.def << "\n"
                   << "    " << enumName(term.prefix, synonym, term.isObsolete) << " = " << eName;
            }
        }
    }
    os << "\n}; // enum CVID\n\n\n";

    os << "/// Information about an ontology or CV source and a short 'lookup' tag to refer to.\n"
          "struct PWIZ_API_DECL CV\n"
          "{\n"
          "    /// the short label to be used as a reference tag with which to refer to this particular Controlled Vocabulary source description (e.g., from the cvLabel attribute, in CVParamType elements).\n"
          "    std::string id;\n"
          "\n"
          "    /// the URI for the resource.\n"
          "    std::string URI;\n"
          "\n"
          "    /// the usual name for the resource (e.g. The PSI-MS Controlled Vocabulary).\n"
          "    std::string fullName;\n"
          "\n"
          "    /// the version of the CV from which the referred-to terms are drawn.\n"
          "    std::string version;\n"
          "\n"
          "    /// returns true iff id, URI, fullName, and version are all pairwise equal\n"
          "    bool operator==(const CV& that) const;\n"
          "\n"
          "    /// returns ture iff id, URI, fullName, and version are all empty\n"
          "    bool empty() const;\n"
          "};\n\n\n";

    os << "/// returns a CV object for the specified namespace (prefix);\n"
          "/// currently supported namespaces are: MS UO\n"
          "PWIZ_API_DECL const CV& cv(const std::string& prefix);\n\n\n";

    os << "/// structure for holding CV term info\n"
          "struct PWIZ_API_DECL CVTermInfo\n"
          "{\n"
          "    CVID cvid;\n"
          "    std::string id;\n"
          "    std::string name;\n"
          "    std::string def;\n"
          "    bool isObsolete;\n"
          "\n"
          "    typedef std::vector<CVID> id_list;\n"
          "    id_list parentsIsA;\n"
          "    id_list parentsPartOf;\n"
          "    std::multimap<std::string, CVID> otherRelations;\n"
          "    std::vector<std::string> exactSynonyms;\n"
          "    std::multimap<std::string, std::string> propertyValues;\n"
          "\n"
          "    CVTermInfo() : cvid((CVID)-1) {}\n"
          "    const std::string& shortName() const;\n"
          "    std::string prefix() const;\n"
          "};\n\n\n";

    os << "/// returns CV term info for the specified CVID\n"
          "PWIZ_API_DECL const CVTermInfo& cvTermInfo(CVID cvid);\n\n\n";

    os << "/// returns CV term info for the specified id (accession number)\n"
          "PWIZ_API_DECL const CVTermInfo& cvTermInfo(const char* id);\n"
          "PWIZ_API_DECL const CVTermInfo& cvTermInfo(const std::string& id);\n\n\n";

    os << "/// returns true iff child IsA parent in the CV\n"
          "PWIZ_API_DECL bool cvIsA(CVID child, CVID parent);\n\n\n";

    os << "/// returns vector of all valid CVIDs\n"
          "PWIZ_API_DECL const std::vector<CVID>& cvids();\n\n\n";

    namespaceEnd(os, basename);

    os << "#endif // " << includeGuard << "\n\n\n";
}


void writeCpp(const vector<OBO>& obos, const string& basename, const bfs::path& outputDir)
{
    string filename = basename + ".cpp";
    bfs::path filenameFullPath = outputDir / filename;
    bfs::ofstream os(filenameFullPath, ios::binary);

    writeCopyright(os, filename);

    os << "#define PWIZ_SOURCE\n\n"
       << "#include \"" << basename << ".hpp\"\n"
       << "#include \"pwiz/utility/misc/Std.hpp\"\n"
       << "#include \"pwiz/utility/misc/Singleton.hpp\"\n"
       << "\n\n";

    namespaceBegin(os, basename);

    os << "namespace {\n\n\n";

    os << "struct TermInfo\n"
          "{\n"
          "    CVID cvid;\n"
          "    const char* id;\n"
          "    const char* name;\n"
          "    const char* def;\n"
          "    bool isObsolete;\n"
          "};\n\n\n";

    os << "const TermInfo termInfos_[] =\n{\n";
    os << "    {CVID_Unknown, \"??:0000000\", \"CVID_Unknown\", \"CVID_Unknown\", false},\n";
    for (vector<OBO>::const_iterator obo=obos.begin(); obo!=obos.end(); ++obo)
    BOOST_FOREACH(const Term& term, obo->terms)
    {
        string& eName = correctedEnumNameMaps[obo-obos.begin()][term.id];

        string correctName = term.name;
        if (eName != enumName(term))
            correctName += " (" + term.prefix + ":" + lexical_cast<string>(enumValue(term, obo-obos.begin())) + ")";

        os << "    {" << eName << ", "
           << "\"" << term.prefix << ":" << (term.prefix != "UNIMOD" ? setw(7) : setw(1) )  << setfill('0') << term.id << "\", "
           << "\"" << correctName << "\", "
           << "\"" << term.def << "\", "
           << (term.isObsolete ? "true" : "false") // setw(7) screws up direct output
           << "},\n";
    }
    os << "}; // termInfos_\n\n\n";

    os << "const size_t termInfosSize_ = sizeof(termInfos_)/sizeof(TermInfo);\n\n\n";

    os << "struct CVIDPair\n"
          "{\n"
          "    CVID first;\n"
          "    CVID second;\n"
          "};\n\n\n";

    os << "CVIDPair relationsIsA_[] =\n{\n";
    for (vector<OBO>::const_iterator obo=obos.begin(); obo!=obos.end(); ++obo)
    BOOST_FOREACH(const Term& term, obo->terms)
    BOOST_FOREACH(const Term::id_type& id, term.parentsIsA)
        os << "    {" << correctedEnumNameMaps[obo-obos.begin()][term.id] << ", "
           << correctedEnumNameMaps[obo-obos.begin()][id] << "},\n";
    os << "}; // relationsIsA_\n\n\n";

    os << "const size_t relationsIsASize_ = sizeof(relationsIsA_)/sizeof(CVIDPair);\n\n\n";

    os << "CVIDPair relationsPartOf_[] =\n{\n";
    for (vector<OBO>::const_iterator obo=obos.begin(); obo!=obos.end(); ++obo)
    BOOST_FOREACH(const Term& term, obo->terms)
    BOOST_FOREACH(const Term::id_type& id, term.parentsPartOf)
        os << "    {" << correctedEnumNameMaps[obo-obos.begin()][term.id] << ", "
           << correctedEnumNameMaps[obo-obos.begin()][id] << "},\n";
    os << "}; // relationsPartOf_\n\n\n";

    os << "const size_t relationsPartOfSize_ = sizeof(relationsPartOf_)/sizeof(CVIDPair);\n\n\n";


    os << "struct OtherRelationPair\n"
          "{\n"
          "    CVID subject;\n"
          "    const char* relation;\n"
          "    CVID object;\n"
          "};\n\n\n";

    os << "OtherRelationPair relationsOther_[] =\n"
       << "{\n"
       << "    {CVID_Unknown, \"Unknown\", CVID_Unknown},\n";
    for (vector<OBO>::const_iterator obo=obos.begin(); obo!=obos.end(); ++obo)
    BOOST_FOREACH(const Term& term, obo->terms)
    for (Term::relation_map::const_iterator jt=term.relations.begin(); jt!=term.relations.end(); ++jt)
    {
        const Term* relatedTerm = NULL;
        vector<OBO>::const_iterator obo2;
        for (obo2=obos.begin(); obo2!=obos.end(); ++obo2)
        {
            if (jt->second.first != obo2->prefix)
                continue;

            set<Term>::const_iterator relatedTermItr = obo2->terms.find(Term(jt->second.second));
            if (relatedTermItr != obo2->terms.end())
            {
                relatedTerm = &*relatedTermItr;
                break;
            }
        }

        if (!relatedTerm)
             cerr << "[writeCpp] Warning: unable to find object of term relationship." << endl;
        else
            os << "    {" << correctedEnumNameMaps[obo-obos.begin()][term.id] << ", "
               << "\"" << jt->first << "\", " << correctedEnumNameMaps[obo2-obos.begin()][relatedTerm->id] << "},\n";
    }
    os << "}; // relationsOther_\n\n\n";

    os << "const size_t relationsOtherSize_ = sizeof(relationsOther_)/sizeof(OtherRelationPair);\n\n\n";


    os << "struct CVIDStringPair\n"
          "{\n"
          "    CVID first;\n"
          "    const char* second;\n"
          "};\n\n\n";

    os << "CVIDStringPair relationsExactSynonym_[] =\n"
          "{\n"
          "    {CVID_Unknown, \"Unknown\"},\n";
    for (vector<OBO>::const_iterator obo=obos.begin(); obo!=obos.end(); ++obo)
    BOOST_FOREACH(const Term& term, obo->terms)
    BOOST_FOREACH(const string& synonym, term.exactSynonyms)
        os << "    {" << correctedEnumNameMaps[obo-obos.begin()][term.id] << ", \"" << synonym << "\"},\n";
    os << "}; // relationsExactSynonym_\n\n\n";

    os << "const size_t relationsExactSynonymSize_ = sizeof(relationsExactSynonym_)/sizeof(CVIDStringPair);\n\n\n";


    os << "struct PropertyValuePair\n"
          "{\n"
          "    CVID term;\n"
          "    const char* name;\n"
          "    const char* value;\n"
          "};\n\n\n";


    typedef pair<string, string> NameValuePair;

    os << "PropertyValuePair propertyValue_[] =\n"
          "{\n"
          "    {CVID_Unknown, \"Unknown\", \"Unknown\"},\n";
    for (vector<OBO>::const_iterator obo=obos.begin(); obo!=obos.end(); ++obo)
    {
        if (obo->prefix != "UNIMOD") // we currently only use UNIMOD properties
            continue;

        BOOST_FOREACH(const Term& term, obo->terms)
        BOOST_FOREACH(const NameValuePair& nameValuePair, term.propertyValues)
        {
            if (!(bal::ends_with(nameValuePair.first, "_classification") ||
                  bal::ends_with(nameValuePair.first, "_position") ||
                  bal::ends_with(nameValuePair.first, "_hidden") ||
                  bal::ends_with(nameValuePair.first, "_site") ||
                  nameValuePair.first == "delta_composition" ||
                  nameValuePair.first == "approved"))
                  continue;
            os << "    {" << correctedEnumNameMaps[obo-obos.begin()][term.id]
                          << ", \"" << nameValuePair.first
                          << "\", \"" << nameValuePair.second
                          << "\"},\n";
        }
    }
    os << "}; // propertyValue_\n\n\n";
    os << "const size_t propertyValueSize_ = sizeof(propertyValue_)/sizeof(PropertyValuePair);\n\n\n";

    os << "class CVTermData : public boost::singleton<CVTermData>\n"
          "{\n"
          "    public:\n"
          "    CVTermData(boost::restricted)\n"
          "    {\n"
          "        for (const TermInfo* it=termInfos_; it!=termInfos_+termInfosSize_; ++it)\n"
          "        {\n"
          "            CVTermInfo temp;\n"
          "            temp.cvid = it->cvid;\n"
          "            temp.id = it->id;\n"
          "            temp.name = it->name;\n"
          "            temp.def = it->def;\n"
          "            temp.isObsolete = it->isObsolete;\n"
          "            infoMap_[temp.cvid] = temp;\n"
          "            cvids_.push_back(it->cvid);\n"
          "        }\n"
          "\n"
          "        for (const CVIDPair* it=relationsIsA_; it!=relationsIsA_+relationsIsASize_; ++it)\n"
          "            infoMap_[it->first].parentsIsA.push_back(it->second);\n"
          "\n"
          "        for (const CVIDPair* it=relationsPartOf_; it!=relationsPartOf_+relationsPartOfSize_; ++it)\n"
          "            infoMap_[it->first].parentsPartOf.push_back(it->second);\n"
          "\n"
          "        for (const OtherRelationPair* it=relationsOther_; it!=relationsOther_+relationsOtherSize_; ++it)\n"
          "            infoMap_[it->subject].otherRelations.insert(make_pair(it->relation, it->object));\n"
          "\n"
          "        for (const CVIDStringPair* it=relationsExactSynonym_; it!=relationsExactSynonym_+relationsExactSynonymSize_; ++it)\n"
          "            infoMap_[it->first].exactSynonyms.push_back(it->second);\n"
          "\n"
          "        for (const PropertyValuePair* it=propertyValue_; it!=propertyValue_+propertyValueSize_; ++it)\n"
          "            infoMap_[it->term].propertyValues.insert(make_pair(it->name, it->value));\n"
          "\n";

    // TODO: is there a way to get these from the OBOs?
    os << "        cvMap_[\"MS\"].fullName = \"Proteomics Standards Initiative Mass Spectrometry Ontology\";\n"
          "        cvMap_[\"MS\"].URI = \"http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo\";\n"
          "\n"
          "        cvMap_[\"UO\"].fullName = \"Unit Ontology\";\n"
          "        cvMap_[\"UO\"].URI = \"http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo\";\n"
          "\n"
          "        cvMap_[\"UNIMOD\"].fullName = \"UNIMOD\";\n"
          "        cvMap_[\"UNIMOD\"].URI = \"http://www.unimod.org/obo/unimod.obo\";\n"
          "\n";

    // populate CV ids and versions from OBO headers
    for (vector<OBO>::const_iterator obo=obos.begin(); obo!=obos.end(); ++obo)
    {
        os << "        cvMap_[\"" << obo->prefix << "\"].id = \"" << obo->prefix << "\";\n";

        string version;
        for (size_t i=0; i < obo->header.size(); ++i)
        {
            boost::regex e(".*?[^-]version: (\\S+)");
            boost::smatch what;
            if (regex_match(obo->header[i], what, e))
            {
                version = what[1];
                break;
            }
        }

        if (version.empty())
        {
            // since UNIMOD doesn't update its 'date' field,
            // we use the maximum "date_time_modified" property_value
            string max_date_time_modified;
            boost::regex e("(\\d+-\\d+-\\d+).*");
            boost::smatch what;
            BOOST_FOREACH(const Term& term, obo->terms)
            BOOST_FOREACH(const NameValuePair& nameValuePair, term.propertyValues)
                if (nameValuePair.first == "date_time_modified" &&
                    regex_match(nameValuePair.second, what, e))
                {
                    if (max_date_time_modified.empty() || what[1] > max_date_time_modified)
                        max_date_time_modified = what[1];
                    continue; // to the next term
                }
            version = max_date_time_modified;
        }

        if (version.empty())
        {
            for (size_t i=0; i < obo->header.size(); ++i)
            {
                boost::regex e("\\s*date: (\\S+).*");
                boost::smatch what;
                if (regex_match(obo->header[i], what, e))
                {
                    version = what[1];
                    break;
                }
            }
        }

        if (version.empty())
            version = "unknown";

        os << "        cvMap_[\"" << obo->prefix << "\"].version = \"" << version << "\";\n\n";
    }

    os << "    }\n"
          "\n"
          "    inline const map<CVID,CVTermInfo>& infoMap() const {return infoMap_;}\n"
          "    inline const map<string,CV>& cvMap() const {return cvMap_;}\n"
          "    inline const vector<CVID>& cvids() const {return cvids_;}\n"
          "\n"
          "    private:\n"
          "    map<CVID,CVTermInfo> infoMap_;\n"
          "    map<string,CV> cvMap_;\n"
          "    vector<CVID> cvids_;\n"
          "};\n\n\n";

    os << "const char* oboPrefixes_[] =\n"
          "{\n";
    for (vector<OBO>::const_iterator obo=obos.begin(); obo!=obos.end(); ++obo)
        os << "    \"" << obo->prefix << "\",\n";
    os << "};\n\n\n";

    os << "const size_t oboPrefixesSize_ = sizeof(oboPrefixes_)/sizeof(const char*);\n\n\n"

          "const size_t enumBlockSize_ = " << enumBlockSize_ << ";\n\n\n"

          "struct StringEquals\n"
          "{\n"
          "    bool operator()(const string& yours) {return mine==yours;}\n"
          "    string mine;\n"
          "    StringEquals(const string& _mine) : mine(_mine) {}\n"
          "};\n\n\n";

    os << "} // namespace\n\n\n";

    os << "#include \"cv.inl\" // code and data that doesn't change with CV releases\n\n";

    namespaceEnd(os, basename);
}


void generateFiles(const vector<OBO>& obos, const string& basename, const bfs::path& outputDir)
{
    // populate term maps for each OBO
    termMaps.resize(obos.size());
    correctedEnumNameMaps.resize(obos.size());
    for (vector<OBO>::const_iterator obo=obos.begin(); obo!=obos.end(); ++obo)
    {
        multiset<string> enumNames;
        BOOST_FOREACH(const Term& term, obo->terms)
            enumNames.insert(enumName(term));

        BOOST_FOREACH(const Term& term, obo->terms)
        {
            termMaps[obo-obos.begin()][term.id] = &term;

            string& eName = correctedEnumNameMaps[obo-obos.begin()][term.id] = enumName(term);
            if (enumNames.count(eName) > 1)
                eName += "_" + lexical_cast<string>(enumValue(term, obo-obos.begin()));
        }
    }

    writeHpp(obos, basename, outputDir);
    writeCpp(obos, basename, outputDir);
}


int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cout << "Usage: cvgen file.obo [...]\n";
        cout << "Parse input file(s) and output cv.hpp and cv.cpp.\n";
        return 1;
    }

    try
    {
        bfs::path exeDir(bfs::path(argv[0]).branch_path());

        vector<OBO> obos;
        for (int i=1; i<argc; i++)
            obos.push_back(OBO(argv[i]));

        generateFiles(obos, "cv", exeDir);

        return 0;
    }
    catch (exception& e)
    {
        cerr << "Caught exception: " << e.what() << endl;
    }
    catch (...)
    {
        cerr << "Caught unknown exception.\n";
    }

    return 1;
}

