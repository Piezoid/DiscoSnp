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

#ifndef _TOOL_BUBBLE_HPP_
#define _TOOL_BUBBLE_HPP_

/********************************************************************************/
#include <cassert>
#include <string>

#include <Filter.hpp> //TODO: could be removed in the header


// #include <gatb/gatb_core.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankFasta.hpp>
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/debruijn/impl/Graph.hpp>
#include <gatb/debruijn/impl/Terminator.hpp>
#include <gatb/debruijn/impl/Traversal.hpp>


/********************************************************************************/

/** We define string constants for command line options. */
#define STR_MAX_AMBIGOUS_INDELS             "-max_ambigous_indel"
#define STR_DISCOSNP_LOW_COMPLEXITY         "-l"
#define STR_DISCOSNP_AUTHORISED_BRANCHING   "-b"
#define STR_DISCOSNP_TRAVERSAL_UNITIG       "-t"
#define STR_DISCOSNP_TRAVERSAL_CONTIG       "-T"
#define STR_KISSNP2_COVERAGE_FILE_NAME      "-coverage_file"
#define STR_KISSNP2_DONT_OUTPUT_FIRST_COV   "-dont_output_first_coverage"

#define STR_MAX_INDEL_SIZE                  "-D"
#define STR_MAX_POLYMORPHISM                "-P"
#define STR_MAX_SYMMETRICAL_CROSSROADS      "-max_symmetrical_crossroads"
#define STR_BFS_MAX_DEPTH                   "-bfs-max-depth"
#define STR_BFS_MAX_BREADTH                 "-bfs-max-breadth"


/********************************************************************************/


namespace dbg = gatb::core::debruijn::impl;
using gatb::core::kmer::Nucleotide;
using gatb::core::kmer::Strand;
using gatb::core::bank::Sequence;
template <size_t span>
class SpecializedGraphClient {
public:

    using Kmer = typename gatb::core::kmer::impl::Kmer<span>::Type;
    using Node = typename dbg::NodeFast<span>;
    using BranchingNode = dbg::BranchingNode_t<Node>;
    using Edge = dbg::Edge_t<Node>;
    using Path = dbg::Path_t<Node>;
    using Graph = dbg::GraphTemplate<Node, Edge, dbg::GraphDataVariantFast<span>>;
    using NodesVec = dbg::GraphVector<Node>;
    using EdgesVec = dbg::GraphVector<Edge>;
    using Traversal = dbg::TraversalTemplate<Node, Edge, Graph>;
    using Terminator = dbg::BranchingTerminatorTemplate<Node, Edge, Graph>;
};






template<typename Node>
inline bool eqNode(Node n1, Node n2) {
    return (n1.kmer == n2.kmer && n1.strand == n2.strand);
}

template<typename Edge>
inline bool eqEdge(Edge e1, Edge e2) {
    if(eqNode(e1.from, e2.from) && e1.nt == e2.nt) {
        assert(eqNode(e1.to, e2.to));
        return true;
    }
    return false;
}


/** \brief Extension state on one path of a bubble (either higher or lower)
 * It doesn't own the extension string but remember the extension lenght at this point
 */

// FIXME: remove
#define FIXED_sizeKmer 31

template<typename Node>
inline Nucleotide lastNT(size_t kmerSize, Node node) {
    //NOTE: Integer.operator[] is faster than graph.getNT because the visitation can be inlined
    if (node.strand == Strand::STRAND_FORWARD)  { return (Nucleotide) (node.kmer[0]); }
    else                                { return gatb::core::kmer::reverse ((Nucleotide) (node.kmer[kmerSize-1])); }
}

template<size_t span>
class PathStack {
public:
    using Graph = typename SpecializedGraphClient<span>::Graph;
    using Edge = typename SpecializedGraphClient<span>::Edge;
    using EdgesVec = typename SpecializedGraphClient<span>::EdgesVec;
    using Node = typename SpecializedGraphClient<span>::Node;
    using NodesVec = typename SpecializedGraphClient<span>::NodesVec;

    class State {
        friend class PathStack;
        Node old, head;
        size_t local_extension_length;

        State(const Node& node) : old(node), head(node), local_extension_length(0) {}

    public:
        inline Nucleotide lastNT() const {
            return ::lastNT(FIXED_sizeKmer, head);
        }

        const Node& getHead() const { return head; }
        Node& getHead() { return head; }

        size_t size() const { return local_extension_length; }
    };


    PathStack(Node& from) : origin(from), path() {}
    PathStack() = default;

    State getInitialState() {
        path.resize(0);
        return State(origin);
    }

    /** \brief Advance the path by one nucleotide
    * Mutate the extension_stack to match the length of the previous extension of this state, the add the nucleotide.
    * \returns true if the extension occurs (the next node differs from the last two previous.)
    */
    inline bool extend(const Node& next, State& state) {
        // Checks if state is in sync:
        assert(state.local_extension_length == path.size());
        assert(state.local_extension_length == 0 || path.back() == ascii(state.lastNT()));

        // Checks if there is no 1-repeats or 2-repeats
        if(eqNode(next, state.head) || eqNode(next, state.old))
            return false;

        state.old = state.head;
        state.head = next;
        state.local_extension_length++;
        path += ascii(state.lastNT());
        return true;
    }

    inline void backtrack(const State& state) {
        assert(path.length() >= state.local_extension_length);
        path.resize(state.local_extension_length);
        assert(state.local_extension_length == 0 || path.back() == ascii(state.lastNT()));
    }

    inline void append(const std::string& extension, const Node& head, const State& state) {
        path += extension;
        state.head = head;
        state.local_extension_length = path.size();
    }

    const std::string& getPath() const { return path; }

    const Node& getOrigin() const { return origin; }
    Node& getOrigin() { return origin; }

    inline string getSequence(const Graph& graph) {
        return graph.toString(origin) + path;
    }

    template<typename Iter>
    void dump(Iter& iter, const Graph& graph) const {
        for(const char c: graph.toString(origin)) *(iter++) = c;
        for(const char c: path) *(iter++) = c;
    }

    size_t extenssion_lenght() const { return path.size(); }

    size_t size(const Graph& graph) const {
        return graph.getKmerSize() + extenssion_lenght();
    }

protected:
    Node origin; //TODO: make me const
    string path;
};




/** We define a structure holding all the information about a bubble. */
template <size_t span>
class BubbleTemplate
{
public:
    using Edge = typename SpecializedGraphClient<span>::Edge;
    using Node = typename SpecializedGraphClient<span>::Node;
    using Path = typename SpecializedGraphClient<span>::Path;
    using PathStack = PathStack<span>;
    using PathState = typename PathStack::State;

    struct State {
        friend class BubbleTemplate;
        PathState higher, lower;
        size_t nb_polymorphism,
               sym_branches, // number of symmetrically branchings traversed (used in b 2 mode)
               stack_size;

        State(PathState&& higher, PathState&& lower) : higher(higher), lower(lower), nb_polymorphism(0), sym_branches(0), stack_size(0)
        {}
    };


//    PathCoupleStack(Node& higher, Node& lower) : higher_path(higher), lower_path(lower) {}

    State getInitialState() {
        return State(higher_path.getInitialState(), lower_path.getInitialState());
    }

    inline bool extend(const pair<Node, Node>& next_nodes, State& state) {
        assert(!eqNode(next_nodes.first, next_nodes.second));
        return higher_path.extend(next_nodes.first, state.higher)
            && lower_path.extend(next_nodes.first, state.lower);
    }

    inline bool extend(const pair<Edge, Edge>& next_edges, State& state) {
        return extend(std::make_pair(next_edges.first.to, next_edges.second.to), state);
    }

    inline void backtrack(const State& state) {
        higher_path.backtrack(state.higher);
        lower_path.backtrack(state.lower);
    }

    PathStack& getHigherPath() { return higher_path; }
    PathStack& getLowerPath() { return lower_path; }

//protected:
    PathStack higher_path, lower_path;

};


template<size_t span>
struct Closure {
    using Path = typename SpecializedGraphClient<span>::Path;
    using Node = typename SpecializedGraphClient<span>::Node;
    using NodesVec = typename SpecializedGraphClient<span>::NodesVec;
    using Graph = typename SpecializedGraphClient<span>::Graph;
    using Traversal = typename SpecializedGraphClient<span>::Traversal;


    Path extenssion;

    // The closure only exists if the bubble have a unique kmer on this side of the bubble
    Nucleotide unique_closure_nucleotide;
    /** Length of the unitig :
    * if the traversal is a simple traversal, then this length is equal to the length of the extension
    *  - else (if the traversal is a monument traversal), then this length is equal to the starting position of the first bubble (if exist))
    */
    size_t divergence;

    Closure(Node origin, const Graph& graph, Traversal* traversal) :
        unique_closure_nucleotide(Nucleotide::NUCL_UNKNOWN), divergence(0)
    {
        // Fist nucleotide is done by hand : this allows to check for the unicity of the bubble closure
        NodesVec successors = graph.successors(origin);

        if (successors.size()==1)
        {
            /** We compute right extension of the node. */
            unique_closure_nucleotide = lastNT(graph.getKmerSize(), successors[0]);
            traversal->traverse (successors[0], dbg::DIR_OUTCOMING, extenssion);
            divergence = traversal->getBubbles().empty() ? extenssion.size() : traversal->getBubbles()[0].first;
        }
    }

    size_t getDivergence() const { return divergence; }

    bool haveExtenssion() const
    { return unique_closure_nucleotide != Nucleotide::NUCL_UNKNOWN; }

    size_t size() const { return haveExtenssion() ? extenssion.size() + 1 : 0 ; }

    Nucleotide operator[](size_t pos) const {
        if(pos == 0) return unique_closure_nucleotide;
        else return extenssion[pos-1];
    }

    template<typename Iter>
    void dump(Iter& iter) const {
        for (size_t i=0; i < size(); i++)
            *(iter++) = tolower(ascii((*this)[i]));
    }

    template<typename Iter>
    void dump_reverse(Iter& iter) const {
        for (size_t i=size(); i-- > 0; )  {
            Nucleotide n = gatb::core::kmer::reverse((*this)[i]);
            *(iter++) = tolower(ascii(n));
        }
    }

};

template<size_t span>
struct FinishedBubble {
    using Path = typename SpecializedGraphClient<span>::Path;
    using Node = typename SpecializedGraphClient<span>::Node;
    using Graph = typename SpecializedGraphClient<span>::Graph;
    using Traversal = typename SpecializedGraphClient<span>::Traversal;

    using Bubble = BubbleTemplate<span>;
    using State = typename Bubble::State;
    using PathStack = typename Bubble::PathStack;
    using Closure = Closure<span>;



    FinishedBubble(Bubble& bubble, const State& state, const Graph& graph, Traversal* traversal) :
        bubble(bubble),
        lower_end(state.higher.getHead()), higher_end(state.lower.getHead()),
        closure_left(graph.reverse(bubble.getHigherPath().getOrigin()), graph, traversal),
        closure_right(higher_end, graph, traversal)

    {
        bubble.backtrack(state);
    }

    void setIndex(size_t _index) { index = _index; }

    PathStack& getHigherPath() { return bubble.getHigherPath(); }
    PathStack& getLowerPath() { return bubble.getLowerPath(); }


    bool haveIndel() const
    { return bubble.getHigherPath().extenssion_lenght() != bubble.getLowerPath().extenssion_lenght(); }


    /** Check whether the first kmer of the first path is smaller than the first kmer
     * of the revcomp(first path), this avoids repeated SNPs
     * \return set isCanonical to true if first path is lower than last reverse path.
     */
    bool isCanonical() {
        // Trick: we test the reverse complement of the second $k$-mer. Thus the last nuc. of the first kmer and the first nuc. of the second kmer does not influence the
        // comparison. Thus the choice of the path used to make this comparison (higher or lower) does not change the results.
        return bubble.higher_path.getOrigin()  < gatb::core::tools::math::revcomp(higher_end.kmer, FIXED_sizeKmer);
    }

    /** Check complexity for a bubble.
     */
    bool isLowComplexity(const Graph& graph) {
        size_t kmer_size = graph.getKmerSize();
        string path1 = graph.toString (getHigherPath().getOrigin()).substr(0, kmer_size-1)
                     + graph.toString (higher_end);
        string path2 = graph.toString (getHigherPath().getOrigin()).substr(0, kmer_size-1)
                     + graph.toString (lower_end);

        /** We compute the low complexity score of the two paths. */
        return filterLowComplexity2Paths (path1, path2);
    }


    size_t dump(bool higher_path, const Graph& graph, gatb::core::tools::misc::Data& data) {
        PathStack& path = higher_path ? getHigherPath() : getLowerPath();
        const size_t len = closure_left.size() + path.size(graph) + closure_right.size();
        if(data.size() <= len) data.resize(len + 1);
        char* ptr = data.getBuffer();

        closure_left.dump_reverse(ptr);
        assert(ptr - data.getBuffer() == closure_left.size());

        path.dump(ptr, graph);
        assert(ptr - data.getBuffer() == closure_left.size() + path.size(graph));

        closure_right.dump(ptr);
        assert(ptr - data.getBuffer() == closure_left.size() + path.size(graph) + closure_right.size());


        *(ptr++) = '\0';
        assert(ptr - data.getBuffer() == len + 1);
        return len;
    }

    string comment(bool higher_path) {
        stringstream comment;
        // PathStack& path = higher_path ? getHigherPath() : getLowerPath();

        comment << (haveIndel() ? "InDel_" : "SNP_") << index
                << (higher_path ? "_higher" : "_lower")
                << "|left_unitig_length=" << closure_left.getDivergence()
                << "|left_contig_length=" << closure_left.size()
                << "|right_unitig_length=" << closure_right.getDivergence()
                << "|right_contig_length=" << closure_left.size();

        return comment.str();
    }

    gatb::core::bank::Sequence sequence(bool higher_path, const Graph& graph) {
        gatb::core::bank::Sequence seq;
        dump(higher_path, graph, seq.getData());
        seq.setComment(comment(higher_path));
        return seq;
    }


protected:
    BubbleTemplate<span>& bubble;
    Node lower_end, higher_end;


    const Closure closure_left, closure_right;
    int final_nb_polymorphism; // number of SNPs in a bubble, could be one (isolated SNP or an insertion) or more (an indel+n SNPs (n+1)) or n SNPs (n)

    // Index of the bubble
    size_t index;
};


class BubbleFinderBase {
public:
    /** We define a structure gathering information during bubble detection. */
    struct Stats
    {
        Stats ()  : nb_bubbles(0), nb_bubbles_snp(0), nb_bubbles_snp_high(0), nb_bubbles_snp_low(0), nb_bubbles_del(0), nb_bubbles_del_high(0), nb_bubbles_del_low(0)  { memset (nb_where_to_extend_snp, 0, sizeof(nb_where_to_extend_snp)); memset (nb_where_to_extend_del, 0, sizeof(nb_where_to_extend_del)); }

        size_t nb_bubbles;

        size_t nb_bubbles_snp;
        size_t nb_bubbles_snp_high;
        size_t nb_bubbles_snp_low;
        size_t nb_where_to_extend_snp[4];



        size_t nb_bubbles_del;
        size_t nb_bubbles_del_high;
        size_t nb_bubbles_del_low;
        size_t nb_where_to_extend_del[4];
    };


    BubbleFinderBase (gatb::core::tools::misc::IProperties* props, Stats& stats);
    BubbleFinderBase (const BubbleFinderBase& bf);
    ~BubbleFinderBase ();


protected:

    /** Statistics about the bubbles lookup. */
    Stats& stats;

    /** Current Bubble instance built by this BubbleFinderTemplate instance. */
//     Bubble<span> bubble;

    /** Shortcut attribute for the kmer size of the de Bruijn graph. */
    size_t sizeKmer;

    /** Maximal number of polymorphism per bubble. Isolated = zero **/
    size_t max_polymorphism;

    /** Max deletion size **/
    int max_indel_size;

    /** Max indel size **/
    int max_indel_ambiguity;

    bool accept_low; // Option set: do we accept low complexity bubbles




    unsigned int max_recursion_depth;
    unsigned int current_recursion_depth;

    unsigned int max_depth;   // for unitigs/contigs extensions
    unsigned int max_breadth; // for unitigs/contigs extensions

    /* authorised_branching =
    *   0: branching forbidden in any path
    *   1: same branching on both path forbidden (i.e. 2 distinct nucleotides may be used in both paths for extension)
    *   2: no restriction on branching */
    int authorised_branching;


    /** In b 2: maximaml number of symetrically branches traversed while walking the bubble**/
    int max_sym_branches;


    /** Gives the kind of traversal to be done at left/right of the bubble. */
    gatb::core::tools::misc::TraversalKind traversalKind;

    /** Output bank of the bubbles (as a pair of sequences). Note here: we use the gatb::core::bank::IBank
     * interface here, and not a specific implementation (like BankFasta), so we could
     * deal with different kinds of banks. */
    gatb::core::bank::IBank* _outputBank;
    void setOutputBank (gatb::core::bank::IBank* outputBank)  { SP_SETATTR(outputBank); }

    /** We need a synchronizer for dumping the sequences into the output bank. */
    gatb::core::system::ISynchronizer* _synchronizer;
    void setSynchronizer (gatb::core::system::ISynchronizer* synchronizer)  { SP_SETATTR(synchronizer); }
};


/********************************************************************************/

/** \brief class that tries to build a bubble from a starting node
 *
 * This class does all the job for expanding a bubble from a starting node if possible and
 * potentially extend it with right and left unitigs/contigs.
 *
 * This class is intended to be instantiated N times, one per thread.
 *
 * The starting node is provided to the operator() (so, one can see this class as a functor).
 */
template <size_t span>
class BubbleFinderTemplate : public BubbleFinderBase, public SpecializedGraphClient<span>
{
public:
    using Bubble = BubbleTemplate<span>;
    using FinishedBubble = FinishedBubble<span>;
    using State = typename Bubble::State;


    using Graph = typename SpecializedGraphClient<span>::Graph;
    using Edge = typename SpecializedGraphClient<span>::Edge;
    using EdgesVec = typename SpecializedGraphClient<span>::EdgesVec;
    using Node = typename SpecializedGraphClient<span>::Node;
    using NodesVec = typename SpecializedGraphClient<span>::NodesVec;
    using Terminator = typename SpecializedGraphClient<span>::Terminator;
    using Traversal = typename SpecializedGraphClient<span>::Traversal;
    using BranchingNode = typename SpecializedGraphClient<span>::BranchingNode;
    using EdgePairsVec = dbg::GraphVector<pair<Edge, Edge>>;
    using NodePairsVec = dbg::GraphVector<pair<Node, Node>>;

    /** Constructor. */
    BubbleFinderTemplate (gatb::core::tools::misc::IProperties* props, const Graph& graph, Stats& stats);

    /** Copy constructor. It will be used when cloning N times the instance by the dispatcher
     * \param[in] bf : instance to be copied.*/
    BubbleFinderTemplate (const BubbleFinderTemplate& bf);

    /** Destructor. */
    ~BubbleFinderTemplate ();

    /** Starting method that gets a node as argument and tries to build a bubble from it.
     * NOTE: defined as a template, because we can use either Node or BranchingNode instances
     * as starting nodes. See also 'start' method, which is template too, with too possible
     * specializations: one for Node, one for BranchingNode
     * \param[in] node : the starting node. */
    template<typename T>
    void operator() (const T& node)
    {
        /** We start the SNP in both directions (forward and reverse). */
        start (node);
        start (graph.reverse(node));
    }

    /** Get a properties object with the configuration of the finder. */
    gatb::core::tools::misc::IProperties* getConfig () const;
    


protected:
    /** */
    const Graph& graph;

    /** Terminator for marking branching nodes (used by the Traversal instance) */
    Terminator* _terminator;
    void setTerminator (Terminator* terminator)  { SP_SETATTR(terminator); }

    /** Used for computing unitigs or contigs (according to traversal kind choice) at the left
     * and right of the bubble. */
    Traversal* _traversal;
    void setTraversal (Traversal* traversal) { SP_SETATTR(traversal); }

    /** Start a bubble detection from a given node. This node is mutated (its last nucleotide) in a
     * second node, and so this couple of nodes is set as the starting branch of a potential bubble.
     * NOTE: defined as template, with 2 specializations: one for Node, one for BranchingNode (see cpp file)
     * \param[in] node : the starting node. */
//     template<typename T>
    void start (const Node& node);
    void start (const BranchingNode& node);
    
    /** Extension of a bubble given two nextNodes to be tested.
     *
     */
    bool expand_heart(Bubble& bubble, State state,
                                 const pair<Edge, Edge>& nextEdge);
    
    /** Extension of a bubble by testing extensions from both branches of the bubble.
     *
     */
    inline bool expand (Bubble& bubble, State& state);
    
    /** Extend the bubble to the left/right with a small assembly part of the de Bruijn graph.
     * \return true if the bubble has been extended, false otherwise. */
    void extend (FinishedBubble &bubble);
    
    /** Finish the bubble, ie output the pair of sequences in the output bank.
     * \param[in] bubble: bubble to be dumped in the output bank
     */
    void finish (FinishedBubble &bubble);

    /** Check bubble according to user choice for branching.
     * \param[in] node 1 : bubble branch last node
     * \param[in] node 2 : bubble branch last node
     * \return true if bubble is ok */
    bool checkBranching (Bubble& bubble, State& state, const EdgePairsVec& successors) const;
    
    /** Check that indel bubbles respect the maximal size of the position ambiguity */
    bool checkRepeatSize (const string &extension1, const string &extension2) const;

    /** Fill a Sequence object for a given branch of a bubble.
     * \param[in] path : branch of the bubble.
     * \param[in] type : used to set the comment part of the sequence (likely 'higher' or 'lower')
     * \param[in] score : score for the bubble.
     * \param[in] where_to_extend : got from 'extend' method.
     * \param[in] seqIndex : index of the sequence (more exactly index for the pair of sequences)
     * \param[out] seq : sequence to be filled
     */
    void buildSequence (Bubble& bubble, size_t pathIdx, const char* type, Sequence& seq, std::string polymorphism_comments);

    /** */
    inline bool isAuthorizedSymmetricBranching(size_t& sym_branches) const;
    bool close(Bubble& bubble, const State& state);
private:
    bool recursive_indel_prediction(
                                    int extended_path_id,
                                    std::string tried_extension,
                                    Node current,
                                    size_t insert_size,
                                    const char end_insertion);
    void start_snp_prediction(Bubble& bubble);
    void start_indel_prediction(Bubble& bubble, State state);
};


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
//#include <Bubble.hpp>
#include <Filter.hpp>
#include <cassert>

#include <string>
using namespace std;

#define DEBUG(a) // a



/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
BubbleFinderBase::BubbleFinderBase (gatb::core::tools::misc::IProperties* props, Stats& stats)
: stats(stats), _outputBank(0), _synchronizer(0)
{
    assert (props != 0);

    /** We set attributes according to user choice. */
    accept_low                  = props->get    (STR_DISCOSNP_LOW_COMPLEXITY) != 0;
    authorised_branching = props->getInt (STR_DISCOSNP_AUTHORISED_BRANCHING);

    max_indel_size       = props->getInt (STR_MAX_INDEL_SIZE);
    max_indel_ambiguity  = props->getInt (STR_MAX_AMBIGOUS_INDELS);
    max_polymorphism     = props->getInt (STR_MAX_POLYMORPHISM);
    max_sym_branches     = props->getInt (STR_MAX_SYMMETRICAL_CROSSROADS);



    max_depth   = props->getInt (STR_BFS_MAX_DEPTH);
    max_recursion_depth=1000; // TODO: parameter?

    max_breadth = props->getInt (STR_BFS_MAX_BREADTH);
    /** We set the traversal kind. */
    traversalKind = gatb::core::tools::misc::TraversalKind::TRAVERSAL_NONE;
    if (props->get(STR_DISCOSNP_TRAVERSAL_UNITIG) != 0)  { traversalKind = gatb::core::tools::misc::TraversalKind::TRAVERSAL_UNITIG; }
    if (props->get(STR_DISCOSNP_TRAVERSAL_CONTIG) != 0)  { traversalKind = gatb::core::tools::misc::TraversalKind::TRAVERSAL_CONTIG; }

    /** We set the name of the output file. */
    stringstream ss;
    ss << props->getStr(STR_URI_OUTPUT);//  << "_k_" << sizeKmer  << "_c_" << graph.getInfo().getInt("abundance");
    //    ss << "_D_"<<max_indel_size;
    ss << ".fa";

    /** We set the output file. So far, we force FASTA output usage, but we could make it configurable. */
    setOutputBank (new gatb::core::bank::impl::BankFasta (ss.str()));

    /** We need a synchronizer for dumping high/low sequences into the output bank in an atomic way
     * (ie avoids potential issues with interleaved high/low sequences in multithread execution). */
    setSynchronizer (gatb::core::system::impl::System::thread().newSynchronizer());
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
template<size_t span>
BubbleFinderTemplate<span>::BubbleFinderTemplate (gatb::core::tools::misc::IProperties* props, const Graph& graph, Stats& stats)
: BubbleFinderBase(props, stats), graph(graph), _terminator(0), _traversal(0)
{
    /** We retrieve the kmer size. */
    sizeKmer = graph.getKmerSize();

    /** We set a terminator here. Note that the construction will set up its inner map
     * with the branching nodes as keys. Then, these keys can be shared with other instances
     * of BranchingTerminator, which leads to less memory usage since this branching nodes keys
     * is not supposed to changed. */
    setTerminator (new Terminator(graph));

    /** We need a synchronizer for dumping high/low sequences into the output bank in an atomic way
     * (ie avoids potential issues with interleaved high/low sequences in multithread execution). */
    setSynchronizer (gatb::core::system::impl::System::thread().newSynchronizer());
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
BubbleFinderBase::BubbleFinderBase (const BubbleFinderBase& bf)
:  stats(bf.stats), _outputBank(0), _synchronizer(0)
{
    sizeKmer             = bf.sizeKmer;
    accept_low           = bf.accept_low;
    authorised_branching = bf.authorised_branching;
    max_indel_size       = bf.max_indel_size;
    max_indel_ambiguity  = bf.max_indel_ambiguity;
    max_polymorphism     = bf.max_polymorphism;
    max_sym_branches     = bf.max_sym_branches;
    max_depth            = bf.max_depth;
    max_recursion_depth  = bf.max_recursion_depth;
    max_breadth          = bf.max_breadth;
    traversalKind        = bf.traversalKind;


    /** Copy by reference (not by value). */
    setOutputBank   (bf._outputBank);
    setSynchronizer (bf._synchronizer);
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
template<size_t span>
BubbleFinderTemplate<span>::BubbleFinderTemplate (const BubbleFinderTemplate& bf)
:  BubbleFinderBase(bf), graph(bf.graph), _terminator(0), _traversal(0)
{
    /** NOT A TRUE COPY: each instance created by this constructor will have its own
     *  Traversal/Terminator instances. */
    setTerminator (new Terminator(*(bf._terminator)));
    setTraversal  (Traversal::create (traversalKind, graph, *_terminator, 0, bf.max_depth, bf.max_breadth));
}


/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
BubbleFinderBase::~BubbleFinderBase ()
{
    setOutputBank   (0);
    setSynchronizer (0);
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
template<size_t span>
BubbleFinderTemplate<span>::~BubbleFinderTemplate ()
{
    setTerminator   (0);
    setTraversal    (0);
}




template<size_t span>
void BubbleFinderTemplate<span>::start_snp_prediction(Bubble& bubble){
    if (max_polymorphism<1) { // if the parameter P is set to 0, do not output any SNP
        return;
    }


    State state = bubble.getInitialState();
    state.nb_polymorphism = 1;
    expand (bubble, state);
}

/** Transform a nucleotide in ASCII form into an integer form as:
 *     - A=0
 *     - C=1
 *     - T=2
 *     - G=3
 * \param[in] nt : the nucleotide in ASCII
 * \return the translated nucleotide */
static int NT2int(char nt)  {  return (nt>>1)&3;  }

template<typename T>
void clear_queue( std::queue<T> &q )
{
    std::queue<T> empty;
    std::swap( q, empty );
}


/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS : breadth first search (queue already stored in the bubbleFinder class.
 Then limit the queue size. Thus early possible bubbles are tested even if longer insers are too complex to be treated.
 *********************************************************************/

/** BFS without have seen support */
template <typename OnExtend, typename OnDepthInc, typename Payload, size_t span>
void bfs(typename SpecializedGraphClient<span>::Graph& graph,
         typename SpecializedGraphClient<span>::Node& start,
         OnExtend on_extend,
         OnDepthInc on_depth_inc
        ) {
    using Node = typename SpecializedGraphClient<span>::Node;
    using NodesVec = typename SpecializedGraphClient<span>::NodesVec;

    // An extension paused on a out-branching path
    struct paused_extension {
        Node head;
        std::string extension;
        Payload payload;
    };

    NodesVec successors;

    // Double buffered stacks used as a queue :
    std::vector<paused_extension> current_stack, next_stack;
    current_stack.emplace_back(start, 0, Payload());
    size_t depth = 0;

    while(!current_stack.empty()) {
        for(paused_extension& paused_ext : current_stack) {
            while(true) { // While the path is single out-branching
                successors = graph.successors(paused_ext.head);
                if (successors.size() == 0 ) break; // dead-end
                else if( successors.size() == 1 ) {
                    // Continue the simple path while the callback returns true
                    paused_ext.head = successors[0];
                    paused_ext.extension.push_back(ascii(lastNT(FIXED_sizeKmer, paused_ext.head)));
                    if( !on_extend(paused_ext.extension, paused_ext.head, depth) )
                        break;
                } else {
                    // If branching, each path is placed in the queue when the callback returns true
                    paused_ext.push_back('\0'); // Will edit
                    for (size_t i=0; i<successors.size(); i++) {
                        paused_ext.head = successors[i];
                        paused_ext.back() = ascii(lastNT(FIXED_sizeKmer, paused_ext.head ));
                        if( on_extend(paused_ext.extension, paused_ext.head, depth) )
                        {
                            next_stack.push_back(paused_ext); // Push a copy of the current state
                        }
                    }
                    break;
                }
            }
        }

        // next depth:
        if(!on_depth_inc(++depth, next_stack.size()))
            break;
        current_stack.clear();
        swap(next_stack, current_stack);
    }
}





template<size_t span>
void BubbleFinderTemplate<span>::start_indel_prediction(Bubble& bubble, State state){
    if (max_indel_size==0)
        return; // no need to try to find indels
    bubble.type=1;
    // Consider a deletion in the upper path (avance on the lower) and then try the opposite



    EdgesVec successors;

    DEBUG(cout << "Try indel from " << graph.toString(bubble.bubble.higher_path.getOrigin()) << " " <<  graph.toString(bubble.bubble.lower_path.getOrigin()) << endl);


    const Edge begin_edges[2] = {
        Edge(Node(), bubble.getHigherPath().getOrigin(), graph.getNT(bubble.higher_path.getOrigin(), sizeKmer-1), dbg::DIR_OUTCOMING),
        Edge(Node(), bubble.getLowerPath().getOrigin(), graph.getNT(bubble.lower_path.getOrigin(), sizeKmer-1), dbg::DIR_OUTCOMING)
    };

    for(int extended_path_id=0;extended_path_id<2;extended_path_id++){
//         if(extended_path_id)
//         cout << "Invert low/high" << endl;
        //        nb_snp_start++;
        //        cout<<"nb_snp_start "<<nb_snp_start<<endl;


        const Edge& ins_begin = begin_edges[extended_path_id];
        const Edge& del_begin = begin_edges[!extended_path_id];

        std::queue<std::pair<Edge,std::string> > breadth_first_queue;
        breadth_first_queue.emplace(ins_begin, string(""));

        DEBUG((cout<<"start indel finding with  "<<ascii(del_begin.nt)<<" extending path "<<extended_path_id<<endl));

        int ps = -1;
        while(!breadth_first_queue.empty()){
            // TODO maybe we could stop the whole breadth first search in case a bubble was found at a found_del_size and the current insert_size is <= found_del_size-2
            DEBUG((cout<<"queue size  "<<breadth_first_queue.size()<<" max_breadth "<<max_recursion_depth<<endl));
            if(breadth_first_queue.size()>max_recursion_depth){
                break; // This bubble is too complex, we stop.
            }
            pair<Edge,string> element=breadth_first_queue.front();
            breadth_first_queue.pop();
            DEBUG((cout<<"and now queue size  "<<breadth_first_queue.size()<<endl));

            DEBUG((cout<<"insert size   "<<element.second.length()<<endl));

            /** if we already found an indel bubble: we need to check to other possible bubbles of the same size (indels of the same size) */
            /** however, if we reach a node at depth lower than the succesfull bubble, as no other bubbles are stored in the queue with */
            /** a higher length (property of the queue), then we can safelly stop the breadth first search */
//             if (element.second.length() < found_del_size){
//                 clear_queue(breadth_first_queue);
//                 break; // ...and stop
//             }
//             /** checks if an indel was already found at a lower depth */
//             // TODO: maybe impossible, to check.
//             if (element.second.length() > found_del_size) {
//                 (cout << "Skipping insert of size " << element.second.length() << " due to previous found of size " << found_del_size << endl);
//                 continue;
//             }


//             if (/*element.second.length() <= ps && */element.second.length() > 0) {
//                 cout << element.second<<endl; //WARNING: DEBUG
//             }
//             ps = element.second.length();



            bool abort_branch = false;
//            while(true) {
//                if (del_begin.nt  == element.first.nt){
//                    DEBUG((cout<<"start an INDEL detection  "<<endl));
//                    bubble.polymorphism_type="INDEL";//+(insert_size);

//                    /** try to close the bubble from the two initial (with one extended) node */


//                    State state;
//                    state.nb_polymorphism = 1;
//                    state.sym_branches = 0;
//                    state.stack_size = 0;
//                    if(extended_path_id == 0) {
//                        state.higher.last_edge = element.first;
//                        state.higher.local_extension_length = element.second.length();
//                        swap(path_stack.first, element.second);

//                        state.lower.last_edge = del_begin;
//                        state.lower.local_extension_length = 0;
//                        path_stack.second.clear();
//                    } else {
//                        state.higher.last_edge = del_begin;
//                        state.higher.local_extension_length = 0;
//                        path_stack.first.clear();

//                        state.lower.last_edge = element.first;
//                        state.lower.local_extension_length = element.second.length();
//                        swap(path_stack.second, element.second);
//                    }

//                    if (expand(bubble, state)) {
//                        found_del_size=element.second.length();
//                        abort_branch = true;
////                         cout << "Found bubble !" << endl;
//                        break;
//                    }

//                }

//                if ( element.second.length() >= found_del_size // no need to try longer extension than the one already found.
//                  || element.second.length() >= max_indel_size      // no need to try longer extension than the maximal length
//                ) {
//                    abort_branch = true;
//                    break;
//                }



//                successors = graph.successorsEdge (element.first.to);
//                if(successors.size() == 1) {
//                    element.first = successors[0];
//                    element.second += ascii(element.first.nt);
//                } else {
//                    if (successors.size() == 0) abort_branch = true;
//                    break;
//                }
//            }

            if(abort_branch) continue;

            if(authorised_branching > 0) {
                for (size_t successor_id=0; successor_id<successors.size() ; successor_id++) {
                    breadth_first_queue.emplace(successors[successor_id],
                                                element.second+ascii(successors[successor_id].nt));
                }
            }
            else {
                break; // Not authorized branching
            }
        }
    }
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
template<size_t span>
void BubbleFinderTemplate<span>::start (const BranchingNode& node)
{
    DEBUG((cout << "[BubbleSNPFinder::start] BRANCHING NODE " << graph.toString(node) << endl));
    DEBUG ((cout << "[BubbleSNPFinder::start] bubble.isCanonical " << bubble.isCanonical << endl));
    /** We compute the successors of the node. */
    NodesVec successors = graph.successors((Node&)node);
    DEBUG((cout << "successor size"<<successors.size()<<endl));
    if(successors.size()<2) return; // false branching (no extention in one or the other direction).

    BubbleTemplate<span> bubble;

    for (size_t i=0; i<successors.size(); i++)
    {
        bubble.lower_path = PathStack<span>(successors[i]);

        // In case two or more branching nodes lead to the same extensions : (b1 -> s1 and b1 -> s2 and b2 -> s1 and b2 -> s2), then
        // we need to construct the bubble s1... s2... from only one of the two branching nodes b1 or b2.
        // We chose the smallest from all the possible starting nodes to start a bubble.
        const NodesVec predecessors = graph.predecessors (bubble.lower_path.getOrigin());

        if (predecessors.size()>1)
        {
            for (size_t k=0; k<predecessors.size(); k++) { if (predecessors[k].kmer < node.kmer) { return; } }
        }


        bubble.higher_path = successors[i];
        for (size_t j=i+1; j<successors.size(); j++)
        {
            bubble.lower_path = successors[j];

            /*************************************************/
            /** Try a SNP                         **/
            /*************************************************/
            DEBUG ((cout << " start SNP detection with " << graph.toString(bubble.bubble.higher_path.getOrigin()) <<" and "<< graph.toString(bubble.bubble.lower_path.getOrigin()) << endl));
            start_snp_prediction(bubble);

            /*************************************************/
            /** Try an isolated insertion                   **/
            /*************************************************/
            //DEBUG ((cout << " start indel detection with " << graph.toString(bubble.bubble.higher_path.getOrigin()) <<" and "<< graph.toString(bubble.bubble.lower_path.getOrigin()) << endl));
            //start_indel_prediction(bubble);
        }
    }
}


template<size_t span>
bool BubbleFinderTemplate<span>::close(Bubble& bubble, const State &state) {
    DEBUG((cout<<"last  node1.value "<<graph.toString(state.higher.getHead())<<" node2.value "<<graph.toString(state.lower.getHead())<<endl));
    /** We check the branching properties of the next kmers. */ //TODO: replicate that



    /** We finish the bubble with last distinct nodes. */
    FinishedBubble finished_bubble(bubble, state, graph, _traversal);

    /** We check several conditions (the first path vs. its revcomp and low complexity). */
    if (finished_bubble.isCanonical() && (accept_low || finished_bubble.isLowComplexity(graph)) && checkRepeatSize(bubble.higher_path.getPath(), bubble.lower_path.getPath()))
    {
        finish (finished_bubble);
        return true;
    } else {
        return false;
    }
}

/*********************************************************************
 ** METHOD  : Heart of the expand method. Factorisation of code, used twice in the expand method.
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  : True if expended, else return false
 ** REMARKS :
 *********************************************************************/
template<size_t span>
bool BubbleFinderTemplate<span>::expand_heart(Bubble& bubble, State state, // NOTE: Here the state is copied
                                 const pair<Edge, Edge>& next_edge){
    /** We check whether the new nodes are different from previous ones. */

    /************************************************************/
    /**                   RECURSION FINISHED                   **/
    /************************************************************/
    if(eqNode(next_edge.first.to, next_edge.second.to))
    { return close(bubble, state);
    }

    /************************************************************/
    /**                   RECURSION CONTINUES                  **/
    /************************************************************/
    else
    {


        DEBUG((cout<<"continue with nextNode1.value "<<graph.toString(state.higher.last_edge.to)<<" nextNode2.value "<<graph.toString(state.lower.last_edge.to)<<endl));
        /** We call recursively the method (recursion on 'pos'). */
        state.stack_size++;
        return bubble.extend(next_edge, state) && expand(bubble, state);

        //            /** There's only one branch to expand if we keep non branching SNPs only, therefore we can safely stop the for loop */
        //            if ( authorised_branching==0 || authorised_branching==1 )   {  break; }

    }

}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  : True if expended, else return false
 ** REMARKS :
 *********************************************************************/
template<size_t span>
inline bool BubbleFinderTemplate<span>::expand (Bubble& bubble, State& state) // NOTE: Here we mutate the parent state
{
    DEBUG((cout<<" *"<<state.stack_size<<"+"<<endl));
    DEBUG((cout<<"expand with node1.value "<<graph.toString(state.higher.last_edge.to)<<" node2.value "<<graph.toString(state.lower.last_edge.to)<<endl));
    DEBUG((cout<<"expand with local_extended_string1 "<<path_stack.first<<" local_extended_string2 "<<path_stack.second<<endl));


    /****************************************************************************/
    /**************** OPTIMIZATION : AVOID RECURSIONS ***************************/


    EdgePairsVec successors;
    while (true){

        successors = graph.successorsEdge (state.higher.getHead(), state.lower.getHead()); // get next two nodes
        if (successors.size() != 1) break; // no more iterative programming

        if (checkBranching(bubble, state, successors) == false) return false; // no possibility to continue

        const pair<Edge,Edge>& successor = successors[0];
        if(eqNode(successor.first.to, successor.second.to))
            return close(bubble, state);

        if(!bubble.extend(successor, state)) { return false ; }
    }

    /**************** END OPTIMIZATION  *****************************************/
    /****************************************************************************/

    /** We may have to stop the extension according to the branching mode. */
    if (checkBranching(bubble, state, successors) == false)  {
        return false;
    }


    /** We get the common successors of node1 and node2. */
    /** Returns the successors of two nodes, ie with the same transition nucleotide from both nodes. */
//    GraphVector < pair<Node,Node> > successors = graph.successors (node1, node2);
    DEBUG((cout<<"successors size "<<successors.size()<<endl));



    bool dumped_bubble=false;
    /** We loop over the successors of the two nodes. */
    size_t i;
    for (i=0; i<successors.size(); i++)
    {
        const pair<Edge,Edge>& successor = successors[i];
//        assert(state.higher.last_edge.to == successor.first.from);
//        assert(state.lower.last_edge.to == successor.second.from);

        /** extend the bubble with the couple of nodes */
        bubble.backtrack(state);
        dumped_bubble |= expand_heart(bubble, state, successor);
        //assert(!dumped_bubble || (dumped_bubble && bubble.closed_bubble));
        /** B 2 special case: if two or more symmetrical branching close a bubble, the output is redundant. **/
        /** Thus, if successor.first = successor.second and if the bubble is dumped, we stop **/
        if (successor.first.to == successor.second.to && dumped_bubble) break;


        // /** Stop as soon as a bubble is dumped */
        // VERSION 2.2.5: commented this break line. Enable to explore all possible symmetrical paths, even in case of success on one of the paths.
        //if(dumped_bubble) break;
    }

    DEBUG((cout<<"stop try"<<endl));

    /** Maybe we can search for a close SNP */
    if (!dumped_bubble && state.nb_polymorphism < max_polymorphism) {
        DEBUG((cout<<"try with a new polymorphism ("<<state.nb_polymorphism<<") with node1.value "<<graph.toString(state.higher.last_edge.to)<<" node2.value "<<graph.toString(state.lower.last_edge.to)<<endl));
        EdgesVec successors1 = graph.successorsEdge (state.higher.getHead());
        EdgesVec successors2 = graph.successorsEdge (state.lower.getHead());

        /** We loop over the successors of the two nodes found with distinct extending nucleotides. */
        for (size_t i1=0; i1<successors1.size(); i1++){
            Edge& edge1 = successors1[i1];
//            assert(edge1.from == state.higher.last_edge.to);
            for (size_t i2=0; i2<successors2.size(); i2++){
                Edge& edge2 = successors2[i2];
//                assert(edge2.from == state.lower.last_edge.to);
                if ( edge1.nt == edge2.nt) continue; // This has already been tested in previous loop

                DEBUG((cout<<"TRYING"<<endl));

                state.nb_polymorphism++;
                bubble.backtrack(state);
                dumped_bubble |= expand_heart(bubble, state, make_pair(edge1, edge2));
                //assert(!dumped_bubble || (dumped_bubble && bubble.closed_bubble));
                /******************************************************************************************* **/
                /** Un-understood Sept 2015 (Pierre). Removed and replaced by the next "break"               **/
                /** if the bubble is finished with THIS couple of tested successors, we stop here.**/
                /** if we don't check this, in b 2 mode we may close a bubble with several distinct couple of node and thus create redondant bubbles **/
                /** if(dumped_bubble && successors1[i1]==successors2[i2]) break; **/
                /******************************************************************************************* **/
                if(dumped_bubble) break;
            }
            if(dumped_bubble) break;
        }
    }
    return dumped_bubble;
}



/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS : In case of indel, the extended nodes are those from the full path of length 2k-1
 *********************************************************************/
template<size_t span>
void BubbleFinderTemplate<span>::extend (FinishedBubble& bubble)
{
    Nucleotide closureLeft  = Nucleotide::NUCL_UNKNOWN;
    Nucleotide closureRight = Nucleotide::NUCL_UNKNOWN;

    /** We may have to extend the bubble according to the user choice. */
    if (traversalKind != gatb::core::tools::misc::TraversalKind::TRAVERSAL_NONE)
    {
        /** We ask for the predecessors of the first node and successors of the last node. */
        NodesVec successors   = graph.successors   (bubble.end[0]);
        NodesVec predecessors = graph.predecessors (bubble.bubble.higher_path.getOrigin());

        /** We need to reset branching nodes between extensions in case of overlapping extensions. */
        _terminator->reset ();

        /** If unique, we keep the left/right extensions. */
        if (successors.size()==1)
        {
            /** We compute right extension of the node. */
            closureRight = graph.getNT (successors  [0], sizeKmer-1);
            _traversal->traverse (successors[0], dbg::DIR_OUTCOMING, bubble.extensionRight);
            bubble.divergenceRight = _traversal->getBubbles().empty() ? bubble.extensionRight.size() : _traversal->getBubbles()[0].first;
        }

        if (predecessors.size()==1)
        {
            /** We compute left extension of the node. */
            closureLeft  = graph.getNT (predecessors[0], 0);
            Node rev_pred = graph.reverse(predecessors[0]);
            _traversal->traverse (rev_pred, dbg::DIR_OUTCOMING, bubble.extensionLeft);
            bubble.divergenceLeft = _traversal->getBubbles().empty() ? bubble.extensionLeft.size() : _traversal->getBubbles()[0].first;
        }
    }

    /** We return a code value according to left/right extensions status. */
    if (closureLeft==Nucleotide::NUCL_UNKNOWN && closureRight==Nucleotide::NUCL_UNKNOWN)  { bubble.where_to_extend = 0; }
    else if (closureLeft!=Nucleotide::NUCL_UNKNOWN && closureRight==Nucleotide::NUCL_UNKNOWN)  { bubble.where_to_extend = 1; }
    else if (closureLeft==Nucleotide::NUCL_UNKNOWN && closureRight!=Nucleotide::NUCL_UNKNOWN)  { bubble.where_to_extend = 2; }
    else if (closureLeft!=Nucleotide::NUCL_UNKNOWN && closureRight!=Nucleotide::NUCL_UNKNOWN)  { bubble.where_to_extend = 3; }

    bubble.closure_left.unique_closure_nucleotide  = closureLeft;
    bubble.closure_right.unique_closure_nucleotide = closureRight;

}
/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
template<size_t span>
void BubbleFinderTemplate<span>::finish (FinishedBubble& bubble)
{


    /** We build two Sequence objects from the information of the bubble. */

    /** We set the bubble index. NOTE: We have to make sure it is protected against concurrent
     * accesses since we may be called here from different threads. */
    bubble.setIndex(__sync_add_and_fetch (&(stats.nb_bubbles), 1));


    /** We have to protect the sequences dump wrt concurrent accesses. We use a {} block with
     * a gatb::core::system::LocalSynchronizer instance with the shared ISynchonizer of the Kissnp2 class. */
    {
        gatb::core::system::LocalSynchronizer sync (_synchronizer);

        /** We insert the two sequences into the output bank. */
        _outputBank->insert (bubble.sequence(false, graph));
        _outputBank->insert (bubble.sequence(true, graph));

        /** Stats update (in concurrent access protection block). */

//        stats.nb_bubbles++;

        if(bubble.haveIndel())
            __sync_add_and_fetch (&(stats.nb_bubbles_del), 1);
        else
            __sync_add_and_fetch (&(stats.nb_bubbles_snp), 1);

        if (!bubble.haveIndel()){
            //stats.nb_where_to_extend_snp[bubble.where_to_extend] ++;
            if (bubble.isLowComplexity(graph))  { stats.nb_bubbles_snp_low++; }
            else                           { stats.nb_bubbles_snp_high++;  }
        }
        else {
            //stats.nb_where_to_extend_del[bubble.where_to_extend] ++;
            if (bubble.isLowComplexity(graph))  { stats.nb_bubbles_del_low++; }
            else                           { stats.nb_bubbles_del_high++;  }
        }
    }


}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
template<size_t span>
void BubbleFinderTemplate<span>::buildSequence (Bubble& bubble, size_t pathIdx, const char* type, Sequence& seq, string polymorphism_comments)
{
    stringstream commentStream;

    /** We build the comment for the sequence. */
    commentStream << bubble.polymorphism_type << "_" << type << "_path_" << bubble.index << "|" << polymorphism_comments << "|" << (bubble.high_complexity ? "high" : "low")<< "|nb_pol_" <<bubble.final_nb_polymorphism;



    /** We may have extra information for the comment. */
    if (traversalKind == gatb::core::tools::misc::TraversalKind::TRAVERSAL_UNITIG)
    {
        commentStream << "|left_unitig_length_";
        commentStream << (bubble.where_to_extend%2==1 ? (bubble.extensionLeft.size() +1) : 0);
        commentStream << "|right_unitig_length_";
        commentStream << (bubble.where_to_extend>1    ? (bubble.extensionRight.size()+1) : 0);
    }
    if (traversalKind == gatb::core::tools::misc::TraversalKind::TRAVERSAL_CONTIG)
    {
        commentStream << "|left_unitig_length_";
        commentStream << (bubble.where_to_extend%2==1 ? (bubble.divergenceLeft +1) : 0);
        commentStream << "|right_unitig_length_";
        commentStream << (bubble.where_to_extend>1    ? (bubble.divergenceRight+1) : 0);

        commentStream << "|left_contig_length_";
        commentStream << (bubble.where_to_extend%2==1 ? (bubble.extensionLeft.size() +1) : 0);
        commentStream << "|right_contig_length_";
        commentStream << (bubble.where_to_extend>1    ? (bubble.extensionRight.size()+1) : 0);
    }

    /** We assign the comment of the sequence. */
    seq.setComment (commentStream.str());

    size_t lenLeft  = bubble.extensionLeft.size ();
    size_t lenRight = bubble.extensionRight.size ();
    size_t len      = sizeKmer + bubble.extended_string[pathIdx].length();

    if (bubble.closureLeft  != Nucleotide::NUCL_UNKNOWN)  { len += 1 + lenLeft;  }
    if (bubble.closureRight != Nucleotide::NUCL_UNKNOWN)  { len += 1 + lenRight; }

    /** We resize the sequence data if needed. Note: +1 for ending '\0'
     * NOTE: we use resize if we need more space, setSize otherwise. */
    if (seq.getData().size() < len+1)  {  seq.getData().resize  (len+1); }
    else                               {  seq.getData().setSize (len+1); }

    char* output = seq.getDataBuffer();

    /** We add the left extension if any. Note that we use lower case for extensions. */
    if (bubble.closureLeft != Nucleotide::NUCL_UNKNOWN)
    {
        for (size_t i=0; i<lenLeft; i++)  {  *(output++) = tolower(ascii (reverse(bubble.extensionLeft [lenLeft-i-1])));  }
        *(output++) = tolower(ascii(bubble.closureLeft));
    }

    /** We add the bubble path. */
    string begin = graph.toString (bubble.begin[pathIdx]);

    /** note that if path overlap > 0 central string is empty **/
    for (size_t i=0; i<sizeKmer; i++)  {  *(output++) = begin[i];  }
    for (size_t i=0; i<bubble.extended_string[pathIdx].length(); i++) { *(output++) = bubble.extended_string[pathIdx][i]; }

    /** We add the right extension if any. Note that we use lower case for extensions. */
    if (bubble.closureRight != Nucleotide::NUCL_UNKNOWN)
    {
        *(output++) =  tolower(ascii(bubble.closureRight));
        for (size_t i=0; i<lenRight; i++)  {  *(output++) = tolower(ascii (bubble.extensionRight[i]));  }
    }

    /** We add a null terminator for the strings. */
    *(output++) = '\0';
}


/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
template<size_t span>
bool BubbleFinderTemplate<span>::checkRepeatSize (const string &extension1, const string &extension2) const
{

    if (extension1.length() == extension2.length()) return true; // This is a SNP
    /** compute the bubble paths */
    /** Expected size of a path = 2k-2 (without first and last kmer common to both alleles
     * This can be smaller un case of repeat position ambiguity
     * k-1-size_of_smallest_extension (without the first kmer thus provides the size of the ambiguity
     * In this code we only compute the length of the extension, expected to be k-1 (see bellow)
     * ACCTGGGA
     * ACCTXXGGGA
     * ACCT -> CCTG -> CTGG -> TGGG -> GGGA  ---------------------> extended string is GGG (size k-1)
     * ACCT -> CCTX -> CTXX -> TXXG -> XXGG -> XGGG -> GGGA ------> extended string is XXGGG (size k-1+size insert
     *
     * Extreme case:
     * ACCTGGGA
     * ACCT|GGGA|GGGA
     *
     * ACCT->CCTG->CTGG->TGGG->GGGA->GGAX ------------------------------> extended string is empty (size 0 = k-1-ambiguity => ambiguity=k-1)
     * ACCT->CCTG->CTGG->TGGG->GGGA->GGAG->GAGG->AGGG->GGGA->GAAX ------> extended string is AGGG (size k = k-1-ambiguity+size_ins = k-1-(k-1)+k => insertion of length k
     **/


    const int size_repeat = FIXED_sizeKmer-2-min(extension1.length(), extension2.length());
    if (size_repeat>max_indel_ambiguity) {
        return false;
    }
    return true;

}


template<size_t span>
inline bool BubbleFinderTemplate<span>::isAuthorizedSymmetricBranching(size_t& sym_branches) const {
    // if authorised_branching==1 no symmetric branching is allowed.
    // Otherwise it depends on the number of symmetric branching seen so far.
    if(authorised_branching < 2 || sym_branches >= max_sym_branches ) {
        return false;
    } else { // The symmetric branching is tolerated beacause authorised_branching == 2 and the limit is not exceded.
        sym_branches++;
        return true;
    }
}

/*********************************************************************
 ** METHOD  : BubbleFinderTemplate<span>::checkBranching
 ** PURPOSE : Checks the branching properties of a couple of nodes. If authorised_branching==0: no branching is authorized in any of the two nodes.
 **           If authorised_branching==1: no symetrical branching authorized from the two nodes (ie. N1-->A and N1-->B and N2-->A and N2-->B not authorized)
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
template<size_t span>
bool BubbleFinderTemplate<span>::checkBranching (Bubble& bubble, State& state, const EdgePairsVec& successors) const
{
    if(authorised_branching==0) {
        return !(successors.size() > 1 // Shortcut (we already know the symmetric successors)
              || graph.isBranching(state.higher.getHead())
              || graph.isBranching(state.lower.getHead()));
    } else {
        // Analyze succesors :
        size_t symmetric_successors = 0, ends = 0;
        for (size_t i=0; i < successors.size(); i++) {
            const pair<Edge,Edge>& succ = successors[i];
            if (eqNode(succ.first.to, succ.second.to)) {
                ends++; // TODO: dump bubble
            } else {
                symmetric_successors++;
            }
        }

        // TODO: dispatch on every cases (ie. trigger bubble dump)
        if(symmetric_successors > 1) return isAuthorizedSymmetricBranching(state.sym_branches);

        // Analyze predecessors :
        // WARNING: graph.predecessors(Node,Node) is not implemented :(
        NodePairsVec predecessors = graph.successors(graph.reverse(state.higher.getHead()), graph.reverse(state.lower.getHead()));
        size_t symmetric_predecessors = 0, starts = 0;
        for (size_t i=0; i < predecessors.size(); i++) {
            const pair<Node,Node>& pred = predecessors[i];
            if (eqNode(pred.first, pred.second)) {
                if(state.higher.size() != 0 && state.lower.size() != 0)
                {
                    cout << "Found alternative begin " << graph.toString(graph.reverse(pred.first))
                         << " past original begin for bubble : " << endl
                         << "\t'" << bubble.higher_path.getPath() << "' ("<< graph.toString(state.higher.getHead()) <<  ")" << endl
                         << "\t'" << bubble.lower_path.getPath() << "' ("<< graph.toString(state.lower.getHead()) <<  ")" << endl;                }
                starts++; // TODO: abort if not cannonical (ie. will be visited with another start node)
            } else {
                symmetric_predecessors++;
            }
        }

        // TODO: dispatch on every cases (see above)
        if(symmetric_predecessors > 1)
            return isAuthorizedSymmetricBranching(state.sym_branches);
        else
            return true;
    }
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
template<size_t span>
gatb::core::tools::misc::IProperties* BubbleFinderTemplate<span>::getConfig () const
{
    gatb::core::tools::misc::IProperties* props = new gatb::core::tools::misc::impl::Properties();

    /** We aggregate information for user. */
    props->add (0, "config",   "");
    props->add (1, "kmer_size",        "%d", sizeKmer);
    props->add (1, "auth_branch",      "%d", authorised_branching);
    props->add (1, "max_indel_size",     "%d", max_indel_size);
    props->add (1, "max_polymorphism", "%d", max_polymorphism);
    props->add (1, "low",              "%d", accept_low);
    props->add (1, "traversal",        "%s", toString (traversalKind).c_str());

    return props;
}


/********************************************************************************/
/******************************** SNPs ******************************************/
/********************************************************************************/


#endif /* _TOOL_BUBBLE_HPP_ */


