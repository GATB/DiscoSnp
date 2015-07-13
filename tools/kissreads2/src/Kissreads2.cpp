//
//  Kissreads2.cpp
//  discoSnp_GATB
//
//  Created by Pierre Peterlongo on 03/07/15.
//  Copyright (c) 2015 Pierre Peterlongo. All rights reserved.
//

#include "kissreads2.h"



using namespace std;

/********************************************************************************
 *
 * A QUICK OVERVIEW OF THE CLASS...
 *
 * This class implements the detection of SNP in a provided de Bruijn graph.
 *
 * It is implemented as a subclass of Tool, so we get all the facilities of the
 * Tool class. We can therefore see two main parts here:
 *
 * 1) constructor: we define all the command line parameters available
 *
 * 2) 'execute' method: this is the main method of the class where the main loop
 *   (ie. iteration over nodes of the graph) is done.
 *
 ********************************************************************************/

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
Kissreads2::Kissreads2 () : Tool ("Kissreads2")
{
    /** We add options known by kissnp2. */
    
    getParser()->push_front (new OptionNoParam (STR_KISSREADS_GENOTYPE,             "Compute genotypes", false));
    getParser()->push_front (new OptionNoParam (STR_KISSREADS_OUTPUT_FASTA,         "Output stnadart Fasta. By default the output is formatted especially for the discoSnp++ pipeline", false));
    
    getParser()->push_front (new OptionOneParam (STR_KISSREADS_SIZE_SEEDS,          "Size of the used seeds (distinct from the size of k)",  false, "25"));
    getParser()->push_front (new OptionOneParam (STR_KISSREADS_INDEX_STRIDE,        "Index Stride", false, "2"));
    
    getParser()->push_front (new OptionOneParam (STR_KISSREADS_SIZE_K,              "Size of k, used as minial overlap and kmer spanning read coherence",  false, "31"));
    getParser()->push_front (new OptionOneParam (STR_KISSREADS_MIN_COVERAGE,        "Minimal coverage", false, "2"));
    getParser()->push_front (new OptionOneParam (STR_KISSREADS_MAX_HAMMING,         "Maximal hamming distance authorized while maping",     false, "1"));
    
    getParser()->push_front (new OptionOneParam (STR_URI_OUTPUT_COHERENT,           "Output coherent file name",                      true));
    getParser()->push_front (new OptionOneParam (STR_URI_OUTPUT_UNCOHERENT,         "Output uncoherent file name",                      true));
    getParser()->push_front (new OptionOneParam (STR_URI_READS_INPUT,               "Input reads",  true));
    getParser()->push_front (new OptionOneParam (STR_URI_PREDICTION_INPUT,          "Input predictions",  true));
    
}




/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/

void Kissreads2::execute ()
{
    
    IProperties* props= getInput();
    
    
   
    
    
   
    BankFasta predictions_bank = BankFasta(props->getStr(STR_URI_PREDICTION_INPUT));
    
    
    
    // We declare a Bank instance.
    BankAlbum banks (props->getStr(STR_URI_READS_INPUT));
    const std::vector<IBank*>& banks_of_queries = banks.getBanks();
    u_int64_t nbReads = banks.estimateNbItems();
    
    GlobalValues gv;
    gv.silent=                  false; //TODO add an option
    gv.size_seeds=              props->getInt (STR_KISSREADS_SIZE_SEEDS);
    gv.index_stride=            props->getInt (STR_KISSREADS_INDEX_STRIDE);
    gv.minimal_read_overlap=    props->getInt (STR_KISSREADS_SIZE_K);
    
    
    gv.number_of_read_sets=     banks_of_queries.size();
    gv.subst_allowed=           props->getInt (STR_KISSREADS_MAX_HAMMING);
    gv.min_coverage=            props->getInt (STR_KISSREADS_MIN_COVERAGE);
    
    gv.compute_genotypes=       props->get    (STR_KISSREADS_GENOTYPE) != 0;
    gv.standard_fasta=          props->get    (STR_KISSREADS_OUTPUT_FASTA) != 0;
    
    gv.set_mask_code_seed();
    
    
    
    
    
    ofstream coherent_out;
    coherent_out.open(props->getStr(STR_URI_OUTPUT_COHERENT));
    
    ofstream uncoherent_out;
    uncoherent_out.open(props->getStr(STR_URI_OUTPUT_UNCOHERENT));
      
     
    getTimeInfo().start ("indexing");
    
    FragmentIndex index(predictions_bank.estimateNbItems());
    
    cout<<"Indexing bank "<<props->getStr(STR_URI_PREDICTION_INPUT)<<" -- "<< gv.number_of_read_sets<< " read set(s) "<<nbReads<<" reads"<<endl;
    index.index_predictions (predictions_bank, gv);       // read and store all starters presents in the pointed file. Index by seeds of length k all these starters.
    
    getTimeInfo().stop ("indexing");
    
    
    
    
    
    

    
    vector<ReadMapper> RMvector;
    size_t  nb_cores = getDispatcher()->getExecutionUnitsNumber ();
    for (int read_set_id=0;read_set_id<gv.number_of_read_sets;read_set_id++) {
        RMvector.push_back(ReadMapper(banks_of_queries[read_set_id],read_set_id,nb_cores));
    }
    
    
    Range<u_int64_t> range (
                            0,gv.number_of_read_sets-1
                            );
    
    // We create an iterator over our integer range.
    // Note how we use the Tool::createIterator method. According to the value of the "-verbose" argument,
    // this method will add some progression bar if needed.
    Iterator<u_int64_t>* iter = createIterator<u_int64_t> (range, "");
    LOCAL (iter);
    

    // Total number of mapped reads
    u_int64_t totalNumberOfMappedReads = 0;
    // We want to get execution time. We use the Tool::getTimeInfo() method for this.
    getTimeInfo().start ("mapping reads");
    // We iterate the range through the Dispatcher we got from our Tool parent class.
    // The dispatcher is configured with the number of cores provided by the "-nb-cores" command line argument.
    for (int read_set_id=0;read_set_id<gv.number_of_read_sets;read_set_id++)
                                                           {
                                                               // MAP ALL READS OF THE READ SET read_set_id
                                                               totalNumberOfMappedReads+= RMvector[read_set_id].map_all_reads_from_a_file(gv,index);
                                                               // SET THE READ COHERENCY OF THIS READ SET.
                                                               RMvector[read_set_id].set_read_coherency(gv,index);
                                                           }
    
    
    
    
    
    
    getTimeInfo().stop ("mapping reads");
    
    
   
    
    
    getTimeInfo().start ("print results");
    print_results_2_paths_per_event(coherent_out, uncoherent_out, index, gv);
    getTimeInfo().stop ("print results");
    
    coherent_out.close();
    uncoherent_out.close();
    
    // We gather some statistics. Note how we use the getInfo() method from the parent class Tool
    // If the verbosity is not 0, all this information will be dumped in the console at the end
    // of the tool execution
    getInfo()->add (1, "Stats");
    getInfo()->add (2, "Total Number of Mapped reads",     "%ld",  totalNumberOfMappedReads);
    getInfo()->add (1, "Outputs");
    getInfo()->add (2, "Number of read coherent predictions",     "%ld",  index.nb_coherent);
    getInfo()->add (2, "Number of read uncoherent predictions",     "%ld",  index.nb_uncoherent);
    
    getInfo()->add (1, getTimeInfo().getProperties("Time"));
}
