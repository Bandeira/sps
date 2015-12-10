//
// $Id: Reader_Bruker_Test.cpp 4129 2012-11-20 00:05:37Z chambm $
//
//
// Original author: Matt Chambers <matt.chambers .@. vanderbilt.edu>
//
// Copyright 2009 Vanderbilt University - Nashville, TN 37232
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


#include "pwiz/utility/misc/unit.hpp"
#include "Reader_Bruker.hpp"
#include "pwiz/utility/misc/VendorReaderTestHarness.hpp"
#include "pwiz/utility/misc/Filesystem.hpp"
#include "pwiz/utility/misc/Std.hpp"

struct IsDirectory : public pwiz::util::TestPathPredicate
{
    bool operator() (const string& rawpath) const
    {
        return bfs::is_directory(rawpath);
    }
};

int main(int argc, char* argv[])
{
    TEST_PROLOG(argc, argv)

    #ifdef PWIZ_READER_BRUKER
    const bool testAcceptOnly = false;
    #else
    const bool testAcceptOnly = true;
    #endif

    try
    {
        pwiz::util::testReader(pwiz::msdata::Reader_Bruker(), testArgs, testAcceptOnly, IsDirectory());
    }
    catch (exception& e)
    {
        TEST_FAILED(e.what())
    }
    catch (...)
    {
        TEST_FAILED("Caught unknown exception.")
    }

    TEST_EPILOG
}
