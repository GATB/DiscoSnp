#include <gatb/gatb_core.hpp>
#include <iostream>
#define STR_OUT_COVERAGE_FILE_NAME              "-coverage_file"

/********************************************************************************/
int main (int argc, char* argv[])
{
    
    
    
    
    /** We create a command line parser. */
    OptionsParser parser ("CoverageDump");
    parser.push_back (new OptionOneParam (STR_URI_INPUT, " set of .h5 files from which one need to grab the coverage cutoff ",   true));
    parser.push_back (new OptionOneParam (STR_OUT_COVERAGE_FILE_NAME, " name of the created .h5 coverage file ",   true));
    
    try
    {
        IProperties* options = parser.parse (argc, argv);
        
        std::string input_h5_file_names = options->getStr(STR_URI_INPUT); // h5 files, separated by comma
        std::string delimiter = ",";
        
        // PREPARE FOR WRITING
        Storage* storage_out = StorageFactory(STORAGE_HDF5).create (options->getStr(STR_OUT_COVERAGE_FILE_NAME), true, false);
        LOCAL (storage_out);
        Group& root_out = storage_out->root();
        Collection<NativeInt64>& myIntegers_out = root_out.getCollection<NativeInt64> ("cutoffs");
        
        
        
       
        // PARSE ALL INPUT H5 files, Dump their content in the output h5.
        size_t pos = 0;
        std::string token;
        while ((pos = input_h5_file_names.find(delimiter)) != std::string::npos) {
            token = input_h5_file_names.substr(0, pos);
            cout<<"reading "<<token<<endl;
            // READING
            Storage* storage = StorageFactory(STORAGE_HDF5).load (token);
            LOCAL (storage);
            Group& root = storage->root();
    
            Group& configGroup = storage->getGroup("configuration");
            stringstream ss; ss << configGroup.getProperty ("xml");
            Properties props; props.readXML (ss);
            int asked_ab_min = props.getInt("abundance_min");
            
            
            if (asked_ab_min != -1) // if asked_ab_min != auto
            {
                myIntegers_out.insert (asked_ab_min);
            }
            else //ecover the automatic value computed by dsk in the /histogram/cutoff field.
            {
                Group& histo = root.getGroup("histogram");
                
                
                Collection<NativeInt64>& myIntegers = histo.getCollection<NativeInt64> ("cutoff");
                // We create an iterator for our collection.
                Iterator<NativeInt64>* iterInt = myIntegers.iterator();
                LOCAL (iterInt);
                // WRITING the content of the input h5.
                for (iterInt->first(); !iterInt->isDone(); iterInt->next())  {
                    int n=iterInt->item().toInt();
                    myIntegers_out.insert (n);
                    break;                          // only one value
                }
                myIntegers.flush();
            }
            
            

            
            
            // Get the next asked coverage valoue
            input_h5_file_names.erase(0, pos + delimiter.length());
        }
        // no need to get the last value as entry finishes with a ',' from run_discoSnp++_ML.sh
        
 

        
        myIntegers_out.flush();

        
    }
    catch (OptionFailure& e)
    {
        return e.displayErrors (std::cerr);
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }
}