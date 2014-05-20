/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
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

#include <Kissnp2.hpp>

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
 * 2) 'execute' method: 'main' method of the class where the main loop (ie. iteration
 *    over nodes of the graph) is done. Each node (and its revcomp) is used as the first
 *    branch of a bubble, the second branch being a mutated node of the first node.
 *    This bubble [node1,node2] is then expanded kmerSize if possible. If still ok, a few
 *    checks are done to avoid redundant computations.
 *
 * A 'configure' method is called at the beginning of 'execute' in order to configure
 * all the attributes (with command line parameters information)
 *
 * There are several 'checkXXX' methods used to avoid unwanted (or duplicate) bubbles.
 *
********************************************************************************/

/** We define string constants for command line options. */
static const char* STR_DISCOSNP_LOW_COMPLEXITY       = "-l";
static const char* STR_DISCOSNP_AUTHORISED_BRANCHING = "-b";
static const char* STR_DISCOSNP_TRAVERSAL_UNITIG     = "-t";
static const char* STR_DISCOSNP_TRAVERSAL_CONTIG     = "-T";
static const char* STR_DISCOSNP_EXTENSION_SIZE       = "-e";

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Kissnp2::Kissnp2 ()
    : Tool ("Kissnp2"), _outputBank(0), _synchronizer(0), nb_bubbles(0), nb_bubbles_high(0), nb_bubbles_low(0),
      low(false), authorised_branching(0),
      sizeKmer (0), min_size_extension(-1), threshold(0),
      traversalKind(Traversal::NONE)
{
    /** We add options known by kissnp2. */
    getParser()->push_front (new OptionNoParam  (STR_DISCOSNP_LOW_COMPLEXITY,       "conserve low complexity SNPs",     false));
    getParser()->push_front (new OptionOneParam (STR_DISCOSNP_AUTHORISED_BRANCHING, "branching mode\n"
            "\t0: forbid SNPs for wich any of the two paths is branching (high precision, low recall)\n"
            "\t1: forbid SNPs for wich the two paths are branching (e.g. the two paths can be created either with a 'A' or a 'C' at the same position (default value)\n"
            "\t2: No limitation on branching (low precision, high recall)",  false, "1"));
    getParser()->push_front (new OptionOneParam (STR_DISCOSNP_EXTENSION_SIZE,       "extend found SNPs  and conserve only those whose min(left and right extension) is bigger or equal to length",  false, "-1"));
    getParser()->push_front (new OptionNoParam  (STR_DISCOSNP_TRAVERSAL_UNITIG,     "extend found and stop at first polymorphism (strict extension=unitigs) SNPs. Uncompatible with -T",  false));
    getParser()->push_front (new OptionNoParam  (STR_DISCOSNP_TRAVERSAL_CONTIG,     "extend found and stop at large polymorphism (extension=contigs) SNPs. Uncompatible with -t",  false));
    getParser()->push_front (new OptionOneParam (STR_URI_OUTPUT,                    "output name",                      true));
    getParser()->push_front (new OptionOneParam (STR_URI_INPUT,                     "input file (likely a hdf5 file)",  true));

    /** Attribute initialization. */
    memset (nb_where_to_extend, 0, sizeof(nb_where_to_extend));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Kissnp2::~Kissnp2 ()
{
    setOutputBank    (0);
    setSynchronizer  (0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Kissnp2::configure ()
{
#if 1
BankFasta::setDataLineSize (100000);
#endif

    /** We load the graph from the provided uri. */
    graph = Graph::load (getInput()->getStr(STR_URI_INPUT));

    /** We retrieve the kmer size. */
    sizeKmer = graph.getKmerSize();

    /** We set the threshold. */
    threshold  = (sizeKmer/2-2)*(sizeKmer/2-3);

    /** We set attributes according to user choice. */
    low                  = getInput()->get    (STR_DISCOSNP_LOW_COMPLEXITY) != 0;
    authorised_branching = getInput()->getInt (STR_DISCOSNP_AUTHORISED_BRANCHING);
    min_size_extension   = getInput()->getInt (STR_DISCOSNP_EXTENSION_SIZE);

    /** We set the traversal kind. */
    if (getInput()->get(STR_DISCOSNP_TRAVERSAL_UNITIG) != 0)  { traversalKind = Traversal::UNITIG; }
    if (getInput()->get(STR_DISCOSNP_TRAVERSAL_CONTIG) != 0)  { traversalKind = Traversal::CONTIG; }

    /** We set the name of the output file. */
    stringstream ss;
    ss << getInput()->getStr(STR_URI_OUTPUT)  << "_k_" << sizeKmer  << "_c_" << graph.getInfo().getInt("abundance");
    if (min_size_extension > -1)  { ss << "_e_" << min_size_extension; }
    ss << ".fa";

    /** We set the output file. So far, we force FASTA output usage, but we could make it configurable. */
    setOutputBank (new BankFasta (ss.str()));

    /** We need a synchronizer for dumping high/low sequences into the output bank in an atomic way
     * (ie avoids potential issues with interleaved high/low sequences in multithread execution). */
    setSynchronizer (System::thread().newSynchronizer());
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
    /** We configure the object with command line arguments. */
    configure ();

    /** We get an iterator over the nodes of the graph. */
    ProgressGraphIterator<Node,ProgressTimer> it (graph.iterator<Node>(), "nodes");

    /** THIS IS THE MAIN ITERATION LOOP... We launch the iteration over all the nodes of the graph.
     *
     * NOTE: we provide an instance of BubbleFinder as functor of the iteration.
     * This instance will be cloned N times, one per thread created by the dispatcher.
     */
    IDispatcher::Status status = getDispatcher()->iterate (it, BubbleFinder(*this, graph));

    /** We aggregate information for user. */
    getInfo()->add (1, "config",   "");
    getInfo()->add (2, "kmer_size",      "%d", sizeKmer);
    getInfo()->add (2, "threshold",      "%d", threshold);
    getInfo()->add (2, "auth_branch",    "%d", authorised_branching);
    getInfo()->add (2, "low",            "%d", low);
    getInfo()->add (2, "traversal_kind", "%d", traversalKind);

    getInfo()->add (1, "bubbles",   "");
    getInfo()->add (2, "nb",      "%lu", nb_bubbles);
    getInfo()->add (2, "nb_high", "%lu", nb_bubbles_high);
    getInfo()->add (2, "nb_low",  "%lu", nb_bubbles_low);
    getInfo()->add (2, "extensions",  "");
    getInfo()->add (3, "none",       "%d", nb_where_to_extend[0]);
    getInfo()->add (3, "left",       "%d", nb_where_to_extend[1]);
    getInfo()->add (3, "right",      "%d", nb_where_to_extend[2]);
    getInfo()->add (3, "left|right", "%d", nb_where_to_extend[3]);

    getInfo()->add (1, "time", "");
    getInfo()->add (2, "find", "%d", status.time);
}

