//! [snippet1]

#include <Kissnp2.hpp>
#include <BubbleFinder.hpp>

using namespace std;

/********************************************************************************/

static const char* STR_LOW_COMPLEXITY       = "-l";
static const char* STR_AUTHORISED_BRANCHING = "-b";

/********************************************************************************/

template<size_t span>
static void executeAlgorithm (Kissnp2& tool, const Graph& graph, IProperties* input)
{
    bool extend_snps        = false;
    bool print_extensions   = true;
    int  min_size_extension = -1;
    bool strict_extension   = true;

    /** We create a bubble instance. */
    BubbleFinder<span> bubble (
        graph,
        (input->getStr(STR_URI_OUTPUT)+string(".fa")).c_str(),
        input->getInt(STR_LOW_COMPLEXITY),
        input->getInt(STR_AUTHORISED_BRANCHING),
        extend_snps,
        min_size_extension,
        print_extensions,
        strict_extension
    );

    // We launch the bubbles lookup.
    bubble.execute ();

    // We aggregate information
    tool.getInfo()->add (1, bubble.getInfo());

    ////    printf("\nFound %lu bubbles (not branching if demanded and with no restriction on closing). Among these, we select the closing bubbles with non-branching kmers (if demanded), from which %lu are high complexity bubbles and %lu are low complexity bubbles (by default, this is 0 if low bubbles do not have to be printed)\n",
    ////        nb_bubbles, nb_bubbles_high, nb_bubbles_low);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Kissnp2::Kissnp2 ()  : Tool ("Kissnp2")
{
    /** We add options specific to Minia (most important at the end). */
    getParser()->push_front (new OptionOneParam (STR_URI_INPUT,             "input file (likely a hdf5 file)",  true));
    getParser()->push_front (new OptionOneParam (STR_URI_OUTPUT,            "output name",                      true));
    getParser()->push_front (new OptionOneParam (STR_LOW_COMPLEXITY,        "conserve low complexity SNPs",     false, "0"));
    getParser()->push_front (new OptionOneParam (STR_AUTHORISED_BRANCHING,  "input file (likely a hdf5 file)",  false, "1"));

    /* authorised_branching =
    *   0: branching forbidden in any path
    *   1: same branching on both path forbidden (i.e. 2 distinct nucleotides may be used in both paths for extension)
    *   2: no restriction on branching
    */
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Kissnp2::execute ()
{
    /** We load the graph from the provided uri. */
    Graph graph = Graph::load (getInput()->getStr(STR_URI_INPUT));

    /** We retrieve the kmer size. */
    size_t kmerSize = graph.getKmerSize();

    /** According to the kmer size, we instantiate one DSKAlgorithm class and delegate the actual job to it. */
         if (kmerSize < KSIZE_1)  { executeAlgorithm <KSIZE_1>  (*this, graph, getInput());  }
    else if (kmerSize < KSIZE_2)  { executeAlgorithm <KSIZE_2>  (*this, graph, getInput());  }
    else if (kmerSize < KSIZE_3)  { executeAlgorithm <KSIZE_3>  (*this, graph, getInput());  }
    else if (kmerSize < KSIZE_4)  { executeAlgorithm <KSIZE_4>  (*this, graph, getInput());  }
    else  { throw Exception ("unsupported kmer size %d", kmerSize);  }
}
