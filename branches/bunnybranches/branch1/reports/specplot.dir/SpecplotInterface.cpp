///////////////////////////////////////////////////////////////////////////////
#include "SpecplotInterface.h"


#include "CommandLineParser.h"
#include "ParameterList.h"

#include "aminoacid.h"
#include "Logger.h"
#include "mzxml.h"

#include "Timer.h"

#include "copyright.h"

///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
SpecplotInterface::SpecplotInterface() :
  specSet(NULL), specSet_own(true)
{
  // Sets draw object in plotSpectrum as PlotGnu
  plotSpectrum.setDrawMethod(RENDERER_TYPE_GNU);
}

SpecplotInterface::~SpecplotInterface()
{
  if(specSet && specSet_own)
    delete specSet;
}
///////////////////////////////////////////////////////////////////////////////
void SpecplotInterface::setData(int type, void *data)
{
  specSet = (SpecSet *)data;
  specSet_own = false;
}
///////////////////////////////////////////////////////////////////////////////
int SpecplotInterface::plot(void)
{
  if(m_spectrumIndex > 0 && m_spectrumIndex <= specSet->size()) {
    // set the spectrum to plot
    plotSpectrum.setSpectrum(&(*specSet)[m_spectrumIndex - 1]);
    // Draws the image
    plotSpectrum.draw();
  }

  if(m_spectrumIndex == -2) {
    for(int i = 0 ; i < specSet->size() ; i++) {

      plotSpectrum.clear();
      plotSpectrum.setSpectrum(&(*specSet)[i]);
      plotSpectrum.setSpectrumIndex(i + 1);
      // set fn suffix
      plotSpectrum.setFnSuffix(parseInt(i+1));

      // Set the title
      string title = "Spectrum ";
      title += parseInt(i + 1);
      plotSpectrum.setTitle(title);

      // Draws the image
      plotSpectrum.draw();
    }
  }

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int SpecplotInterface::processOptions(int argc, char **argv)
{
  //--------------------------------------------------------------------------------------------
  // Initialize directories used
  //--------------------------------------------------------------------------------------------
  // Get the execultable directory. Unless specified, renderer directory and font directory is expected
  string str = argv[0];

  //convert path
  size_t found;
  found = str.find_last_of("/\\");
  string aux;
  aux = '.';
  if(found != string::npos)
    aux = str.substr(0, found);
  string exePath = aux;

//#if  defined(__MINGW32__) || defined(__CYGWIN__)
  exePath = replaceAll(exePath, "\\", "/");
//#endif

  // set renderer location
  plotSpectrum.setRendererLocation(aux);
  // Font location also defaults to executale directory
  plotSpectrum.setFontLocation(exePath);
  // Annotation file directory also defaults to executable directory
  plotSpectrum.setAnnotationModelDirectory(aux);
  // Amino acid file directory also defaults to executable directory
  plotSpectrum.setAminoAcidsFileDirectory(aux);

  //--------------------------------------------------------------------------------------------
  // Parse the command line parameters
  //--------------------------------------------------------------------------------------------
  vector<CommandLineParser::Option> listOptions;

  //  c++ instruction                            cmd line option    parameter name    takes value?

  listOptions.push_back(CommandLineParser::Option("-help",            "help",             false));
  listOptions.push_back(CommandLineParser::Option("-version",         "VERSION",          false));
  listOptions.push_back(CommandLineParser::Option("-annot",           "ANNOTATIONS",      true));
  listOptions.push_back(CommandLineParser::Option("-annot-by-intensity","ANNOTATIONS_BY_INTENSITY",      false));
//  listOptions.push_back(CommandLineParser::Option("-verbose",         "VERBOSE",          false));

  listOptions.push_back(CommandLineParser::Option("-notitle",         "NOTITLE",          false));
  listOptions.push_back(CommandLineParser::Option("-title",           "TITLE",            true));
//  listOptions.push_back(CommandLineParser::Option("-shift",           "SHIFT",            false));
  listOptions.push_back(CommandLineParser::Option("-shift-value",     "SHIFT_VALUE",      true));

  listOptions.push_back(CommandLineParser::Option("-aa",              "AMINO_ACID_LIST",  true));
  listOptions.push_back(CommandLineParser::Option("-annotation-model","ANNOTATION_MODEL", true));
  listOptions.push_back(CommandLineParser::Option("-model-default",   "MODEL_DEFAULT",    false));
  listOptions.push_back(CommandLineParser::Option("-renderer-dir",    "RENDERER_LOCATION",true));
  listOptions.push_back(CommandLineParser::Option("-font-dir",        "FONT_LOCATION",    true));

  listOptions.push_back(CommandLineParser::Option("-pklbin",          "PKLBIN",           true));
  listOptions.push_back(CommandLineParser::Option("-mgf",             "MGF",              true));
  listOptions.push_back(CommandLineParser::Option("-pkl",             "PKL",              true));
  listOptions.push_back(CommandLineParser::Option("-mzxml",           "MZXML",            true));
  listOptions.push_back(CommandLineParser::Option("-infile",          "INFILE",           true));


  listOptions.push_back(CommandLineParser::Option("-spectrum",        "SPECTRUM",         true));
  listOptions.push_back(CommandLineParser::Option("-spectrumscan",    "SPECTRUMSCAN",     true));
  listOptions.push_back(CommandLineParser::Option("-peptide",         "PEPTIDE",          true));
  listOptions.push_back(CommandLineParser::Option("-annotation-style-inspect","ANNOT_STYLE_INSPECT",       false));

  listOptions.push_back(CommandLineParser::Option("-peakmasstol",     "TOLERANCE_PEAK",   true));
  listOptions.push_back(CommandLineParser::Option("-parentmasstol",   "TOLERANCE_PM",     true));

  listOptions.push_back(CommandLineParser::Option("-outdir",          "OUTDIR",           true));
  listOptions.push_back(CommandLineParser::Option("-prefix",          "PREFIX",           true));
  listOptions.push_back(CommandLineParser::Option("-outfile",         "OUTFILE",          true));
  listOptions.push_back(CommandLineParser::Option("-format",          "FORMAT",           true));
  listOptions.push_back(CommandLineParser::Option("-encoding",        "ENCODING",         true));
  listOptions.push_back(CommandLineParser::Option("-target",          "TARGET",           true));

  listOptions.push_back(CommandLineParser::Option("-zoom",            "ZOOM",             true));
  listOptions.push_back(CommandLineParser::Option("-image-width",     "IMAGE_WIDTH",      true));
  listOptions.push_back(CommandLineParser::Option("-image-height",    "IMAGE_HEIGHT",     true));
  listOptions.push_back(CommandLineParser::Option("-font-size",       "FONT_SIZE",        true));


  listOptions.push_back(CommandLineParser::Option("-range",           "RANGE",            true));

  listOptions.push_back(CommandLineParser::Option("-gnuplot-file",    "GNUPLOT_FILE",     false));

  listOptions.push_back(CommandLineParser::Option("-label",           "LABEL",            true));
//  listOptions.push_back(CommandLineParser::Option("-spectruminfo",    "spectruminfo",     true));
//  listOptions.push_back(CommandLineParser::Option("-inspect",         "inspect",          true));
//  listOptions.push_back(CommandLineParser::Option("-inspect",         "inspect",          true));

  listOptions.push_back(CommandLineParser::Option("-request-id",      "REQ_ID",           true));

  // parameter file
  listOptions.push_back(CommandLineParser::Option("-p",               "PARAMETER_FILE",   true));

  ////////////////////////////////////////////////////////////////////////////////
  // Execute the command line parser
  CommandLineParser clp(argc, argv, 0, listOptions);

  ////////////////////////////////////////////////////////////////////////////////
  // Checking for errors
  string parser_error;
  if (!clp.validate(parser_error)) {
    ERROR_MSG(parser_error);
    return -1;
  }


  ////////////////////////////////////////////////////////////////////////////////
  // The parameters' values
  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);

  return processOptions(commandLineParams);
}
////////////////////////////////////////////////////////////////////////////////
int SpecplotInterface::processOptions(ParameterList & ip)
{
  ParameterList commandLineParams;

  ////////////////////////////////////////////////////////////////////////////////
  // Parameter file
  ////////////////////////////////////////////////////////////////////////////////

  if (ip.exists("PARAMETER_FILE")) {

    string parameterFilename = ip.getValue("PARAMETER_FILE");

    commandLineParams.readFromFile(parameterFilename);
    // Combine the command line parameters to the file ones
    //   Command line parameters take precedence (hence the overwrite flag set)
  }

  commandLineParams.addList(ip, true);

  ////////////////////////////////////////////////////////////////////////////////
  // "help" prints help and exits
  if (commandLineParams.exists("VERSION"))
    return version(cout);

  ////////////////////////////////////////////////////////////////////////////////
  // the same for "version"
  if (commandLineParams.exists("help"))
    return help(cout);

  ////////////////////////////////////////////////////////////////////////////////
  // Boolean parameters

//  if (commandLineParams.exists("SHIFT"))
//    plotSpectrum.setMassShift(0.0);

  if (commandLineParams.exists("SHIFT_VALUE"))
    plotSpectrum.setMassShift(getFloat(commandLineParams.getValue("SHIFT_VALUE").c_str()));

  if (commandLineParams.exists("NOTITLE"))
    plotSpectrum.setTitlePresence(0);

  if (commandLineParams.exists("ANNOT_STYLE_INSPECT"))
    plotSpectrum.setAnnotatinStyle(ANNOTATION_STYLE_INSPECT);




  ////////////////////////////////////////////////////////////////////////////////
  // File locations

  // Renderer location
  if (commandLineParams.exists("RENDERER_LOCATION")) {
    string aux = commandLineParams.getValue("RENDERER_LOCATION");
    aux = replaceAll(aux, "\\", "/");
    plotSpectrum.setRendererLocation(aux);
  }

  // Font file location
  if (commandLineParams.exists("FONT_LOCATION")) {
    string aux = commandLineParams.getValue("FONT_LOCATION");
    aux = replaceAll(aux, "\\", "/");
    plotSpectrum.setFontLocation(aux);
  }


  // Amino-acid masses file
  if (commandLineParams.exists("AMINO_ACID_LIST")) {
    string aux = commandLineParams.getValue("AMINO_ACID_LIST");
    aux = replaceAll(aux, "\\", "/");
    int found = aux.find_last_of("/\\");
    if(found != string::npos) {
      plotSpectrum.setAminoAcidsFile(aux.substr(found + 1));
      plotSpectrum.setAminoAcidsFileDirectory(aux.substr(0, found));
    } else {
      plotSpectrum.setAminoAcidsFile(aux);
      //plotSpectrum.setAminoAcidsFileDirectory(".");
    }
  }

  // Annotation model
  if (commandLineParams.exists("ANNOTATION_MODEL")) {
    plotSpectrum.setModelFromFile();
    string aux = commandLineParams.getValue("ANNOTATION_MODEL");
    aux = replaceAll(aux, "\\", "/");
    int found = aux.find_last_of("/\\");
    if(found != string::npos) {
      plotSpectrum.setAnnotationModel(aux.substr(found + 1));
      plotSpectrum.setAnnotationModelDirectory(aux.substr(0, found));
    } else {
      plotSpectrum.setAnnotationModel(aux);
      //plotSpectrum.setAnnotationModelDirectory(".");
    }
  }

  if (commandLineParams.exists("MODEL_DEFAULT"))
    plotSpectrum.setModelFromFile();


  ////////////////////////////////////////////////////////////////////////////////
  // Peptide and amino acid parameters

  // peptide(s)
  if (commandLineParams.exists("PEPTIDE")) {
    string peptide = commandLineParams.getValue("PEPTIDE");
    transform(peptide.begin(), peptide.end(), peptide.begin(), ::toupper);
    vector<string> aux;
    stringSplit(peptide, aux, "|");
    for(int i = 0 ; i < aux.size() ; i++)
      plotSpectrum.addPeptide(aux[i]);
  }


  ////////////////////////////////////////////////////////////////////////////////
  // Spectrum related parameters: part 1

  bool spectrumScanPresent  = false;
  bool spectrumPresent      = false;

  // Spectrumscan
  if (commandLineParams.exists("SPECTRUMSCAN")) {
    m_spectrumScan = commandLineParams.getValue("SPECTRUMSCAN");
    plotSpectrum.setSpectrumScan(m_spectrumScan);
    string title = "Scan ";
    title += m_spectrumScan;
    plotSpectrum.setTitle(title);
    spectrumScanPresent = true;
  }


  // spectrum index used
  if (commandLineParams.exists("SPECTRUM")) {

    spectrumPresent = true;

    if(spectrumScanPresent) {
      //stringstream err; err << "ERROR: Spectrum index and spectrum scan specified.";
      //return error(err.str());
    }

    string spectrumValue = commandLineParams.getValue("SPECTRUM");
    m_spectrumIndex = -1;
    if(!spectrumValue.compare("all"))
      m_spectrumIndex = -2;
    else
      m_spectrumIndex = getInt(spectrumValue.c_str());

    plotSpectrum.setSpectrumIndex(m_spectrumIndex);
    // Set the title
    string title = "Spectrum ";
    title += parseInt(m_spectrumIndex);
    plotSpectrum.setTitle(title);

  } else if(!spectrumScanPresent)
    // A spectrum index must be specified
    return error("No spectrum index or spectrum scan specified.");


  ////////////////////////////////////////////////////////////////////////////////
  // File load section

  // A file must be loaded
  bool fileLoaded = false;

  // file name processor class
  FilenameManager fm;

  if(commandLineParams.exists("PKLBIN")) {
    fm.filenameFull = commandLineParams.getValue("PKLBIN");
    fm.splitFilename();
    fm.extension = "pklbin";
  }

  if(commandLineParams.exists("MGF")) {
    fm.filenameFull = commandLineParams.getValue("MGF");
    fm.splitFilename();
    fm.extension = "mgf";
  }

  if(commandLineParams.exists("PKL")) {
    fm.filenameFull = commandLineParams.getValue("PKL");
    fm.splitFilename();
    fm.extension = "pkl";
  }

  if(commandLineParams.exists("MZXML")) {
    fm.filenameFull = commandLineParams.getValue("MZXML");
    fm.splitFilename();
    fm.extension = "mzxml";
  }

  if(commandLineParams.exists("INFILE")) {
    fm.filenameFull = commandLineParams.getValue("INFILE");
    fm.splitFilename();
    fm.lowerCaseExtension();
  }

  // load te file
  if(!specSet) {
    // create object
    specSet = new SpecSet();
    // translate path separators
    fm.filenameFull = replaceAll(fm.filenameFull, "\\", "/");

    // load mzxml file
    if(fm.extension.compare("mzxml") == 0) {
      // auxilizary vector needed
      vector<short> msLevel;
      // hold the return value
      int ret;
      // load method depends on the presence of spectrumscan specification
      try {
        if(m_spectrumScan.empty()) {
          fileLoaded = LoadMzxml( (char * const)(fm.filenameFull.c_str()), *specSet, & msLevel, 2);
        } else {
          fileLoaded = LoadMzxml(fm.filenameFull.c_str(), *specSet, m_spectrumScan.c_str(), & msLevel, 2);
        }
      } catch (...) {
        fileLoaded = false;
      }
    }

    if(!fileLoaded) {

      // load other formats

      fileLoaded = specSet->Load(fm.filenameFull.c_str());
      if(!fileLoaded) {
        stringstream err;
        err << "Error loading file: " << fm.filenameFull;
        return error(err.str());
      } else if(!spectrumPresent) {
          // get spectrumscan from index
          bool found = false;
          for(int i = 0 ; i < specSet->size() ; i++)
            if((*specSet)[i].scan == getInt(m_spectrumScan.c_str())) {
              m_spectrumIndex = i;
              plotSpectrum.setSpectrumIndex(m_spectrumIndex);
              found = true;
            }

          if(!found) {
            stringstream err; err << "ERROR: Spectrumscan not found.";
            return error(err.str());
          }
          //stringstream err; err << "ERROR: Spectrum index not specified.";
          //return error(err.str());
        }
    }
  }


  // Test spectrum index to draw
  if(spectrumPresent) {
    if( ((m_spectrumIndex < 1) || (m_spectrumIndex > specSet->size())) && (m_spectrumIndex != -2)) {
      stringstream err; err << "ERROR: Spectrum index out of range: " << specSet->size() << " < " << m_spectrumIndex;
      return error(err.str());
    }
  }

  // Check if a file was loaded
  if(!fileLoaded && !specSet)
    return error("No file loaded.");

  // Annotations input file
  if (commandLineParams.exists("LABEL"))
    plotSpectrum.setAnnotationInputFile(commandLineParams.getValue("LABEL"));


  ////////////////////////////////////////////////////////////////////////////////
  //  File save section

  // Output directory
  if (commandLineParams.exists("OUTDIR")) {
    string aux = commandLineParams.getValue("OUTDIR");
    aux = replaceAll(aux, "\\", "/");
    plotSpectrum.setFileOutDir(aux);
  }

  // Output file name.
  if (commandLineParams.exists("OUTFILE")) {
    string aux = commandLineParams.getValue("OUTFILE");
    aux = replaceAll(aux, "\\", "/");
    plotSpectrum.setFileName(aux);
  }

  // Filename prefix
  if (commandLineParams.exists("PREFIX"))
   plotSpectrum.setFilePrefix(commandLineParams.getValue("PREFIX"));

  // Output file format: png, ...
  if (commandLineParams.exists("FORMAT"))
    plotSpectrum.setFileFormat(commandLineParams.getValue("FORMAT"));

  // Output file encoding: uu64, ...
  if (commandLineParams.exists("ENCODING"))
    plotSpectrum.setEncoding(commandLineParams.getValue("ENCODING"));

  // Output target: file, cout, cout, internal
  if (commandLineParams.exists("TARGET"))
    plotSpectrum.setTarget(commandLineParams.getValue("TARGET"));

  // Annotations output file
  if (commandLineParams.exists("ANNOTATIONS"))
    plotSpectrum.setAnnotationOutputFile(commandLineParams.getValue("ANNOTATIONS"));

  if (commandLineParams.exists("ANNOTATIONS_BY_INTENSITY"))
    plotSpectrum.setAnnotationByIntensity();


  ////////////////////////////////////////////////////////////////////////////////
  //  Spectrum related parameters


  // Spectrum scan processing
  if(spectrumScanPresent) {
    m_spectrumIndex = 1;
    plotSpectrum.setSpectrumIndex(1);
//    unsigned scan = getInt(plotSpectrum.m_spectrumScan.c_str());
//    for(int i = 0 ; i < specSet->size() ; i++)
//      if( (*specSet)[i].scan == scan) {
//        plotSpectrum.m_spectrumIndex = i;
//        break;
//      }
  }



  // spectrumInfo
  if (commandLineParams.exists("spectruminfo"))
    plotSpectrum.setSpectrumInfo(commandLineParams.getValue("spectruminfo"));

  // peak mass tolerance
  if (commandLineParams.exists("TOLERANCE_PEAK"))
    plotSpectrum.setPeakMassTol(getFloat(commandLineParams.getValue("TOLERANCE_PEAK").c_str()));


  ////////////////////////////////////////////////////////////////////////////////
  //  Image related parameters

  // Title
  if (commandLineParams.exists("TITLE"))
    plotSpectrum.setTitle(commandLineParams.getValue("TITLE"));

  // Output image zoom factor
  if (commandLineParams.exists("ZOOM"))
    plotSpectrum.setZoom(getFloat(commandLineParams.getValue("ZOOM").c_str()));

  // Output image zoom factor
  if (commandLineParams.exists("IMAGE_WIDTH"))
    plotSpectrum.setImageWidth(getInt(commandLineParams.getValue("IMAGE_WIDTH").c_str()));

  // Output image zoom factor
  if (commandLineParams.exists("IMAGE_HEIGHT"))
    plotSpectrum.setImageHeight(getInt(commandLineParams.getValue("IMAGE_HEIGHT").c_str()));

  // Output image zoom factor
  if (commandLineParams.exists("FONT_SIZE"))
    plotSpectrum.setFontSize(getFloat(commandLineParams.getValue("FONT_SIZE").c_str()));



  // m/z  axis range.
  if (commandLineParams.exists("RANGE")) {
    string range = commandLineParams.getValue("RANGE");
    vector<string> aux;
    stringSplit(range, aux, ":");
    float faux;
    switch(aux.size()) {
    case 1:
      faux = getFloat(aux[0].c_str());
      if(range[0] ==':')
        plotSpectrum.setRangeMax(faux);
      else
        plotSpectrum.setRangeMin(faux);
      break;
    case 2:
      float rangeMin = getFloat(aux[0].c_str());
      float rangeMax = getFloat(aux[1].c_str());
      plotSpectrum.setRangeMin(rangeMin);
      plotSpectrum.setRangeMax(rangeMax);
      if(rangeMin >= rangeMax)
        return error("range: minimun value greater or equal than maximun value.");
      break;
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Other

  // keep gnuplot file
  if (commandLineParams.exists("GNUPLOT_FILE"))
    plotSpectrum.setDebug(1);


  // return status ok
  return 0;
}
///////////////////////////////////////////////////////////////////////////////
int SpecplotInterface::help(ostream & sout)
{
  version(sout);

  sout << "Usage: specplot [OPTION]\n";
  sout << "Options:\n";
  sout << "  --infile FILE               Spectrum data\n";
  sout << "  --pklbin FILE               Spectrum data (read in pklbin format)\n";
  sout << "  --mgf FILE                  Spectrum data (read in mgf format)\n";
  sout << "  --pkl FILE                  Spectrum data (read in pkl format)\n";
  sout << "  --mzxml FILE                Spectrum data (read in mzxml format)\n";
  sout << '\n';
  sout << "  --aa FILE                   Amino acids file (txt format)\n";
  sout << "  --annotation-model FILE     Annotation model file (txt format)\n";
  sout << "  --model-default             Use default annotation model\n";
  sout << "  --font-dir DIRECTORY        Font files location in the file system.\n";
  sout << "  --renderer-dir DIRECTORY    Renderer location in the file system.\n";
  sout << "  --label FILE                Explicit labels [mass, height, type, order, charge, peptide, prob]\n";
  sout << '\n';
  sout << "  --spectrum INDEX            Spectrum (X) 1-based\n";
  sout << "  --spectrumscan INDEX        Spectrum scan (X)\n";
//  sout << "  --spectruminfo STRING       Spectrum info\n";
  sout << "  --peakmasstol NUMBER        Amino acid peak mass tolerance\n";
  sout << "  --parentmasstol NUMBER      Parent mass tolerance\n";
  sout << "  --peptide STRING            Contiguous amino acid letters. Several peptides may be specified, separated by the '|' character\n";
  sout << "  --annotation-style-inspect  Set annotation style to INSPECT. Default is SPECNETS\n";
  sout << "  --annot FILE                Write annotations [mass, height, type, order, charge, peptide, prob] to file\n";
  sout << "  --annot-by-intensity        Order annotations by intensity instead of by m/z\n";
//  sout << "  --shift                     Apply mass shift\n";
  sout << "  --shift-value               Specify and apply mass shift (e.g., use --shift-value -1 for CID PRM spectra)\n";
  sout << "  --notitle                   Disables title generation\n";
  sout << "  --title STRING              Specifies title\n";
  sout << "  --zoom NUMBER               Zoom factor\n";
  sout << "  --image-width NUMBER        Image width, in pixels\n";
  sout << "  --image-height NUMBER       Image height, in pixels\n";
  sout << "  --font-size NUMBER          Font size\n";
  sout << "  --range X:Y                 Zoom spectrum to m/z range between m/z=X and m/z=Y\n";
  sout << '\n';
  sout << "  --format STRING             Image output format (png or eps)\n";
  sout << "  --outdir PATH               Output directory\n";
  sout << "  --prefix STRING             Output filename prefix\n";
  sout << "  --outfile FILE              Output filename (overrides prefix)\n";
  sout << "  --target STRING             Output target (file, cout, internal)\n";
  sout << "  --encoding STRING           Output encoding (uu64, default is no encoding)\n";
  sout << '\n';
  sout << "  --p FILE                    Read parameters from file\n";
  sout << '\n';
  sout << "  --help                      Display this help and exit\n";
  sout << "  --version                   Output version information and exit\n";
  sout << endl;

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int SpecplotInterface::version(ostream & sout)
{
  sout << "specplot 2.0." << XSTR(SPS_VERSION) << '\n';
  sout << COPYRIGHT1 << '\n';
  sout << COPYRIGHT2;
  sout << "\n\n";

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int SpecplotInterface::error(const string & a)
{
  cerr << "specplot: invalid or missing option: " << a << endl;
  cerr << "Type 'specplot --help' for more information." << endl << endl;

  return -1;
}
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
