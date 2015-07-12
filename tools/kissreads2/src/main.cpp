
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
    
//    // We get the number of cores to be used.  If we don't give any number,
//    // we set to 0 which implies the usage of all available cores
//    size_t nbCores = (0);
//    // We create an iterator over an integer range
//    Range<int>::Iterator it (1,1000);
//    // We create a dispatcher configured for 'nbCores' cores.
//    Dispatcher dispatcher (nbCores);
//    // We dispatch the range iteration with the dispatcher.
//    // This will create nbCores threads and each thread will be fed with
//    // one value of the defined range
//    // NOTE: we could also use lambda expression (easing the code readability)
//    IDispatcher::Status status = dispatcher.iterate (it, Functor());
//    // We dump some information about the dispatching
//    cout << "nbCores=" << status.nbCores << "  time=" << status.time << endl;
//    // IMPORTANT: usage of Dispatcher has sense only if the iterated items
//    // can be processed independently from each other.
//    // The point to understand with the Dispatcher is that it can
//    // iterate any instance of Iterator class. If you have any set of items
//    // that can be enumerated through an Iterator implementation, then you
//    // can parallelize the iteration with a Dispatcher instance
    
    cout<<"TODO: / TESTS FASTQ / Standart fasta option"<<endl;    
    
    
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

