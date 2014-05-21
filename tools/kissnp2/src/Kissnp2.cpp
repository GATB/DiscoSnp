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
#include <Bubble.hpp>

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
Kissnp2::Kissnp2 () : Tool ("Kissnp2")
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
    /** We define which kind of nodes we want to use as starter for bubble building.
     *  It could be either BranchingNode or Node */
    typedef BranchingNode StartingNode;

    /** We load the graph from the provided uri. */
    Graph graph = Graph::load (getInput()->getStr(STR_URI_INPUT));

    /** We get an iterator over the nodes of the graph. */
    ProgressGraphIterator<StartingNode,ProgressTimer> it (graph.iterator<StartingNode>(), "nodes");

    /** We want to get some statistics about the execution. */
    BubbleFinder::Stats stats;

    /** We create an instance of BubbleFinder, used as a functor by the dispatcher.
     * This instance will be cloned N times, one per thread created by the dispatcher.
     */
    BubbleFinder bubbleFinder (getInput(), graph, stats);

    /** THIS IS THE MAIN ITERATION LOOP... We launch the iteration over all the nodes of the graph.
     * Each iterated node is sent in one of N threads where it is provided to the operator() method
     * of one of the N BubbleFinder instance.
     */
    IDispatcher::Status status = getDispatcher()->iterate (it, bubbleFinder);

    /** We aggregate information for user. */
    getInfo()->add (1, bubbleFinder.getConfig());
    getInfo()->add (1, "nodes",    "%lu", it.size());
    getInfo()->add (1, "bubbles",  "");
    getInfo()->add (2, "nb",      "%lu", stats.nb_bubbles);
    getInfo()->add (2, "nb_high", "%lu", stats.nb_bubbles_high);
    getInfo()->add (2, "nb_low",  "%lu", stats.nb_bubbles_low);
    getInfo()->add (2, "extensions",  "");
    getInfo()->add (3, "none",       "%d", stats.nb_where_to_extend[0]);
    getInfo()->add (3, "left",       "%d", stats.nb_where_to_extend[1]);
    getInfo()->add (3, "right",      "%d", stats.nb_where_to_extend[2]);
    getInfo()->add (3, "left|right", "%d", stats.nb_where_to_extend[3]);
    getInfo()->add (1, "time", "");
    getInfo()->add (2, "find", "%d", status.time);
}

