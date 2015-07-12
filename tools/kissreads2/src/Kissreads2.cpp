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


// We define a functor that will be cloned by the dispatcher
struct Functor { void operator() (ReadMapper rm)
    {
//        rm.map_all_reads_from_a_file(gv,index);
//        rm.set_read_coherency(gv,index);
//        map_all_reads_from_a_file (banks_of_queries[read_set_id],
//                                   read_set_id,
//                                   gv,
//                                   index
//                                   );
        // In this instruction block, we are executing in one of the nbCores threads
        // created by the dispatcher. Note that 'i' is one value of our range
    }};


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
      
     
    
    
    FragmentIndex index(predictions_bank.estimateNbItems());
    
    cout<<"Indexing bank "<<props->getStr(STR_URI_PREDICTION_INPUT)<<" -- "<< gv.number_of_read_sets<< " read set(s)"<<endl;
    index.index_predictions (predictions_bank, gv);       // read and store all starters presents in the pointed file. Index by seeds of length k all these starters.
    
    
    cout<<"Mapping "<< gv.number_of_read_sets<< " read set(s) on prediction(s)"<<endl;
    
#ifdef OMP
#pragma omp parallel for if(nbthreads>1 && sam_out==NULL) num_threads(nbthreads) private(read_set_id)
#endif
    // We create an iterator over an integer range
    vector<ReadMapper> vector;
    for (int read_set_id=0;read_set_id<gv.number_of_read_sets;read_set_id++) vector.push_back(ReadMapper(banks_of_queries[read_set_id],read_set_id));
    
    
//    Dispatcher::Status status = getDispatcher()->iterate (*(vector.begin()), Functor());
    
    for(std::vector<ReadMapper>::iterator it = vector.begin(); it != vector.end(); ++it) {
        it->map_all_reads_from_a_file(gv,index);
        it->set_read_coherency(gv,index);
    }
    
    
//    for (int read_set_id=0;read_set_id<gv.number_of_read_sets;read_set_id++){
//        
//        cout<<"map read set "<<read_set_id<<endl;
//        ReadMapper rm (banks_of_queries[read_set_id],read_set_id);
//        rm.map_all_reads_from_a_file(gv,index);
//        rm.set_read_coherency(index,gv);
//        
//    }
    
    
    
    printf("Results analysis done, Print results\n");
    print_results_2_paths_per_event(coherent_out, uncoherent_out, index, gv);
    
    coherent_out.close();
    uncoherent_out.close();
    
    printf("printing done, free memory, and finish\n");
   
    
//    free(seed_table);
//    FreeHashTable((struct HashTable*)seeds_count);
//    for (i=0;i<nb_events_per_set;i++)
//    {
//        free(results_against_set[i]->read_coherent);
//        free(results_against_set[i]->number_mapped_reads);
//        int read_set_id;
//        for (read_set_id=0;read_set_id<number_of_read_sets;read_set_id++)
//        {
//            free(results_against_set[i]->read_coherent_positions[read_set_id]);
//        }
//        free(results_against_set[i]->read_coherent_positions);
//        free(results_against_set[i]->sum_qualities);
//        free(results_against_set[i]->nb_mapped_qualities);
//        free(results_against_set[i]);
//    }
//    free(results_against_set);
//    
//    
//    fclose(coherent_out);
//    if (!only_print) fclose(uncoherent_out);
//    if(sam_out) fclose(sam_out);
//    if(!silent) {
//        if(only_print) printf("Results are in %s", coherent_file_name);
//        else printf("Results are in %s and %s", coherent_file_name, uncoherent_file_name);
//        if(sam_out) printf(" (mapped reads are in %s)", samout_file_name);
//        printf(".\n");
//        after_all = time(NULL);
//        printf("Total time: %.0lf secs\n",  difftime(after_all, before_all));
//    }
//    return 0;
}
