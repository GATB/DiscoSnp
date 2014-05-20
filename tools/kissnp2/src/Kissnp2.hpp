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

#ifndef _TOOL_KISSNP2_HPP_
#define _TOOL_KISSNP2_HPP_

/********************************************************************************/
#include <gatb/gatb_core.hpp>
#include <Bubble.hpp>
#include <queue>
/********************************************************************************/


/********************************************************************************/
/** \brief Tool class that looks for SNP
 *
 * The Kissnp2 is the front class for SNP detection in a provided de Bruijn graph.
 * The output is a bank with pairs of sequences defining a bubble.
 *
 * This class may be inherited by refining some methods; for instance, the pairs
 * of sequences in the output bank may have some left/right extensions from the
 * found bubbles.
 */
class Kissnp2 : public Tool
{
public:

    /** Constructor. */
    Kissnp2 ();

    /** Destructor. */
    ~Kissnp2 ();

    /** Implementation of Tool::execute method. */
    void execute ();

protected:

    /** */
    void configure ();

    /** De Bruijn graph, likely built by the 'dbgh5' binary from GATB-CORE. */
    Graph graph;

    /** Shortcut attribute for the kmer size of the de Bruijn graph. */
    size_t sizeKmer;

    /** Threshold (computed from the kmer size). */
    int threshold;

    /** Output bank of the bubbles (as a pair of sequences). Note here: we use the IBank
     * interface here, and not a specific implementation (like BankFasta), so we could
     * deal with different kinds of banks. */
    IBank* _outputBank;
    void setOutputBank (IBank* outputBank)  { SP_SETATTR(outputBank); }

    bool low;

    /* authorised_branching =
    *   0: branching forbidden in any path
    *   1: same branching on both path forbidden (i.e. 2 distinct nucleotides may be used in both paths for extension)
    *   2: no restriction on branching */
    int authorised_branching;

    int  min_size_extension;

    /** We need a synchronizer for dumping the sequences into the output bank. */
    ISynchronizer* _synchronizer;
    void setSynchronizer (ISynchronizer* synchronizer)  { SP_SETATTR(synchronizer); }

    /** Statistics about the bubbles lookup. */
    size_t nb_bubbles;
    size_t nb_bubbles_high;
    size_t nb_bubbles_low;
    size_t nb_where_to_extend[4];

    /** */
    Traversal::Kind traversalKind;

    friend class BubbleFinder;
};

/********************************************************************************/

#endif /* _TOOL_KISSNP2_HPP_ */

