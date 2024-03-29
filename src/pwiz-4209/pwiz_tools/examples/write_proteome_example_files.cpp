//
// $Id: write_proteome_example_files.cpp 2911 2011-08-05 17:29:12Z chambm $
//
//
// Original author: Matt Chambers <matt.chambers .@. vanderbilt.edu>
//
// Copyright 2010 Vanderbilt University - Nashville, TN 37232
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


#include "pwiz/data/proteome/ProteomeDataFile.hpp"
#include "pwiz/data/proteome/examples.hpp"
#include "pwiz/utility/misc/Std.hpp"


using namespace pwiz::proteome;


void writeTiny()
{
    ProteomeData pd;
    examples::initializeTiny(pd);

    // write out FASTA 
    string filename = "tiny.fasta";
    cout << "Writing file " << filename << endl;
    // call after writer creation
    ProteomeDataFile::write(pd, filename);
}


int main()
{
    try
    {
        writeTiny();

        cout << "\nhttp://proteowizard.sourceforge.net\n"
             << "support@proteowizard.org\n";

        return 0;
    }
    catch (exception& e)
    {
        cerr << e.what() << endl;
    }
    catch (...)
    {
        cerr << "Caught unknown exception.\n";
    }

    return 1;
}
