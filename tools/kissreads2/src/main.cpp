
/*****************************************************************************
 *   discoSnp++: discovering polymorphism from raw unassembled NGS reads
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2014  INRIA
 *   Authors: P.Peterlongo, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

/********************************************************************************/

#include <Kissreads2.h>


using namespace std;

using namespace gatb::core;
using namespace gatb::core::tools;
using namespace gatb::core::bank;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

/********************************************************************************/


//// We define a functor that will be cloned by the dispatcher
//struct Functor { void operator() (int i)
//    {
//        // In this instruction block, we are executing in one of the nbCores threads
//        // created by the dispatcher. Note that 'i' is one value of our range
//        cout<<i<<endl;
//    }};

int main (int argc, char* argv[])
{
    u_int64_t a=17;
    u_int64_t b=32;
    u_int64_t nbcreated ;
    
    
//    Hash16<u_int64_t, u_int64_t> * tested_starters_and_pwis = new Hash16<u_int64_t, u_int64_t>( 65 , &nbcreated );
//    cout<<"nbcreated = "<<nbcreated;
////    
//    if (tested_starters_and_pwis->get(a, &b) == 0){
//        cout<<"yé pas"<<endl;
//        tested_starters_and_pwis->insert(a, b);
//    }
//    cout<<tested_starters_and_pwis->get(a, &b)<<endl;
//    return 1;
    
    cout<<"TODO: / TESTS FASTQ / SEUIL PAR JEUX DE DONNEES"<<endl;
    
    
    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        /** We execute the tool. */
        Kissreads2().run (argc, argv);
    }
    
    catch (OptionFailure& e)
    {
        return e.displayErrors (cout);
    }
    
    catch (Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}

