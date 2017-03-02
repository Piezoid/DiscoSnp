/*****************************************************************************
x *   discoSnp++: discovering polymorphism from raw unassembled NGS reads
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

#include <Bubble.hpp>
#include <Filter.hpp>
#include <cassert>

#include <string>
using namespace std;

#define DEBUG(a) // a
const char* BubbleFinder::STR_BFS_MAX_DEPTH   = "-bfs-max-depth";
const char* BubbleFinder::STR_BFS_MAX_BREADTH = "-bfs-max-breadth";



/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
BubbleFinder::BubbleFinder (IProperties* props, const Graph& graph, Stats& stats)
: graph(graph), stats(stats), _terminator(0), _traversal(0), _outputBank(0), _synchronizer(0)
{
    assert (props != 0);
    
    /** We retrieve the kmer size. */
    sizeKmer = graph.getKmerSize();
    
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
    traversalKind = TRAVERSAL_NONE;
    if (props->get(STR_DISCOSNP_TRAVERSAL_UNITIG) != 0)  { traversalKind = TRAVERSAL_UNITIG; }
    if (props->get(STR_DISCOSNP_TRAVERSAL_CONTIG) != 0)  { traversalKind = TRAVERSAL_CONTIG; }
    
    /** We set the name of the output file. */
    stringstream ss;
    ss << props->getStr(STR_URI_OUTPUT);//  << "_k_" << sizeKmer  << "_c_" << graph.getInfo().getInt("abundance");
    //    ss << "_D_"<<max_indel_size;
    ss << ".fa";
    
    /** We set the output file. So far, we force FASTA output usage, but we could make it configurable. */
    setOutputBank (new BankFasta (ss.str()));
    
    /** We need a synchronizer for dumping high/low sequences into the output bank in an atomic way
     * (ie avoids potential issues with interleaved high/low sequences in multithread execution). */
    setSynchronizer (System::thread().newSynchronizer());

    /** We set a terminator here. Note that the construction will set up its inner map
     * with the branching nodes as keys. Then, these keys can be shared with other instances
     * of BranchingTerminator, which leads to less memory usage since this branching nodes keys
     * is not supposed to changed. */
    setTerminator (new BranchingTerminator(graph));
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
BubbleFinder::BubbleFinder (const BubbleFinder& bf)
:  graph(bf.graph), stats(bf.stats), _terminator(0), _traversal(0), _outputBank(0), _synchronizer(0)
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
    breadth_first_queue  = bf.breadth_first_queue;

    
    /** Copy by reference (not by value). */
    setOutputBank   (bf._outputBank);
    setSynchronizer (bf._synchronizer);
    
    /** NOT A TRUE COPY: each instance created by this constructor will have its own
     *  Traversal/Terminator instances. */
    setTerminator (new BranchingTerminator(*(bf._terminator)));
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
BubbleFinder::~BubbleFinder ()
{
    setOutputBank   (0);
    setSynchronizer (0);
    setTerminator   (0);
    setTraversal    (0);
}





void BubbleFinder::start_snp_prediction(){
    if (max_polymorphism<1) { // if the parameter P is set to 0, do not output any SNP
        return;
    }
    bubble.polymorphism_type="SNP";
    bubble.type=0;
    bubble.extended_string[0]="";
    bubble.extended_string[1]="";
    Node notzero = Node(~0);
    Node notzero2 = Node(~0);


    SNP_State state;
    state.nb_polymorphism = 1;
    state.sym_branches = 0;
    state.stack_size = 0;

    state.higher.last_edge.to = bubble.begin[0];
    state.higher.last_edge.from = Node(~0);
    state.higher.last_edge.nt = graph.getNT(bubble.begin[0], sizeKmer-1);
    state.higher.local_extension_length = 0;
    path_stack.first = "";

    state.lower.last_edge.to = bubble.begin[1];
    state.lower.last_edge.from = Node(~0);
    state.lower.last_edge.nt = graph.getNT(bubble.begin[1], sizeKmer-1);
    state.lower.local_extension_length = 0;
    path_stack.second = "";

    expand (state);
}

/** Transform a nucleotide in ASCII form into an integer form as:
 *     - A=0
 *     - C=1
 *     - T=2
 *     - G=3
 * \param[in] nt : the nucleotide in ASCII
 * \return the translated nucleotide */
static int NT2int(char nt)  {  return (nt>>1)&3;  }


void clear_queue_pair( std::queue<pair<Node, string> > &q )
{
    std::queue<pair<Node, string> > empty;
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
void BubbleFinder::start_indel_prediction(){
    if (max_indel_size==0)
        return; // no need to try to find indels
    bubble.type=1;
    // Consider a deletion in the upper path (avance on the lower) and then try the opposite
    
    Node current;
    int found_del_size=max_indel_size+10; // No indel found for now (we could also use MAXINT)
    

    for(int extended_path_id=0;extended_path_id<2;extended_path_id++){
        //        nb_snp_start++;
        //        cout<<"nb_snp_start "<<nb_snp_start<<endl;
        current= bubble.begin[extended_path_id]; // 0 or 1
        Nucleotide end_insertion=graph.getNT(bubble.begin[!extended_path_id], sizeKmer-1);
        
        DEBUG((cout<<"start indel finding with  "<<ascii(end_insertion)<<" extending path "<<extended_path_id<<endl));
        string tried_extension;
        breadth_first_queue.emplace(current, string(""));
        while(!breadth_first_queue.empty()){
            // TODO maybe we could stop the whole breadth first search in case a bubble was found at a found_del_size and the current insert_size is <= found_del_size-2
            DEBUG((cout<<"queue size  "<<breadth_first_queue.size()<<" max_breadth "<<max_recursion_depth<<endl));
            if(breadth_first_queue.size()>max_recursion_depth){
                clear_queue_pair(breadth_first_queue);
                break; // This bubble is too complex, we stop.
            }
            pair<Node,string> element=breadth_first_queue.front();
            breadth_first_queue.pop();
            DEBUG((cout<<"and now queue size  "<<breadth_first_queue.size()<<endl));
            current = element.first;
            tried_extension=element.second;
            int insert_size = tried_extension.length();
            DEBUG((cout<<"insert size   "<<insert_size<<endl));
            /** if we already found an indel bubble: we need to check to other possible bubbles of the same size (indels of the same size) */
            /** however, if we reach a node at depth lower than the succesfull bubble, as no other bubbles are stored in the queue with */
            /** a higher length (property of the queue), then we can safelly stop the breadth first search */
            if (insert_size == found_del_size-1){
                clear_queue_pair(breadth_first_queue);
                break; // ...and stop
            }
            /** checks if an indel was already found at a lower depth */
            // TODO: maybe impossible, to check.
            if (insert_size > found_del_size) {
                DEBUG(cout << "Skipping insert of size " << insert_size << " due to previous found of size " << found_del_size << endl);
                continue;
            }
            if (end_insertion  == graph.getNT(current, sizeKmer-1)){
                DEBUG((cout<<"start an INDEL detection  "<<endl));
                bubble.polymorphism_type="INDEL";//+(insert_size);
                
                /** try to close the bubble from the two initial (with one extended) node */


                SNP_State state;
                state.nb_polymorphism = 1;
                state.sym_branches = 0;
                state.stack_size = 0;
                if(extended_path_id == 0) {
                    state.higher.last_edge.to = current;
                    state.higher.last_edge.from = Node(~0);
                    state.higher.last_edge.nt = graph.getNT(current, sizeKmer-1);
                    state.higher.local_extension_length = tried_extension.length();
                    path_stack.first = tried_extension;

                    state.lower.last_edge.to = bubble.begin[1];
                    state.lower.last_edge.from = Node(~0);
                    state.lower.last_edge.nt = graph.getNT(bubble.begin[1], sizeKmer-1);
                    state.lower.local_extension_length = 0;
                    path_stack.second = "";
                } else {
                    state.higher.last_edge.to = bubble.begin[1];
                    state.higher.last_edge.from = Node(~0);
                    state.higher.last_edge.nt = graph.getNT(bubble.begin[1], sizeKmer-1);
                    state.higher.local_extension_length = 0;
                    path_stack.first = "";

                    state.lower.last_edge.to = current;
                    state.lower.last_edge.from = Node(~0);
                    state.lower.last_edge.nt = graph.getNT(current, sizeKmer-1);
                    state.lower.local_extension_length = tried_extension.length();
                    path_stack.second = tried_extension;
                }

                if (expand(state)) {
                    found_del_size=insert_size;
                    continue;
                }

            }
            
            if(
               insert_size >= found_del_size || // no need to try longer extension than the one already found.
               insert_size >= max_indel_size      // no need to try longer extension than the maximal length
               ) continue;
            Graph::Vector<Node> successors = graph.successors (current);
            
            /** No branching authorized in the insertion mode. */
            if (successors.size()>1 && authorised_branching==0) {
                clear_queue_pair(breadth_first_queue);
                break; // ...and stop
            }

            /** checks if a successor with the good starting letter (the one potentially closing the indel) exists */
            bool exists;
            Node successor = graph.successor(current,end_insertion,exists);
            if(exists)
                breadth_first_queue.emplace(successor,tried_extension+ascii(end_insertion));
            
            /** then checks for the other possible extensions */
            for (size_t successor_id=0; successor_id<successors.size() ; successor_id++) {
                Nucleotide last_extension = graph.getNT(successors[successor_id], sizeKmer-1);
                if(last_extension != end_insertion)
                    breadth_first_queue.emplace(successors[successor_id], tried_extension+ascii(last_extension));
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
template<>
void BubbleFinder::start (Bubble& bubble, const BranchingNode& node)
{
    this->bubble=bubble;
    DEBUG ((cout << "[BubbleSNPFinder::start] BRANCHING NODE " << graph.toString(node) << endl));
    DEBUG ((cout << "[BubbleSNPFinder::start] bubble.isCanonical " << bubble.isCanonical << endl));
    /** We compute the successors of the node. */
    Graph::Vector<Node> successors = graph.successors((Node&)node);
    DEBUG((cout << "successor size"<<successors.size()<<endl));
    if(successors.size()<2) return; // false branching (no extention in one or the other direction).
    for (size_t i=0; i<successors.size(); i++)
    {
        bubble.begin[0] = successors[i];
        
        // In case two or more branching nodes lead to the same extensions : (b1 -> s1 and b1 -> s2 and b2 -> s1 and b2 -> s2), then
        // we need to construct the bubble s1... s2... from only one of the two branching nodes b1 or b2.
        // We chose the smallest from all the possible starting nodes to start a bubble.
        Graph::Vector<Node> predecessors = graph.predecessors (successors[i]);
        
        if (predecessors.size()>1)
        {
            for (size_t k=0; k<predecessors.size(); k++) { if (predecessors[k].kmer < node.kmer) { return; } }
        }
        
        
        
        for (size_t j=i+1; j<successors.size(); j++)
        {
            bubble.begin[1] = successors[j];
            bubble.isCanonical=false;
            bubble.closed_bubble=false;
            /*************************************************/
            /** Try a SNP                         **/
            /*************************************************/
            DEBUG ((cout << " start SNP detection with " << graph.toString(bubble.begin[0]) <<" and "<< graph.toString(bubble.begin[1]) << endl));
            start_snp_prediction();
            
            /*************************************************/
            /** Try an isolated insertion                   **/
            /*************************************************/
            bubble.isCanonical=false;
            bubble.closed_bubble=false;
            DEBUG ((cout << " start indel detection with " << graph.toString(bubble.begin[0]) <<" and "<< graph.toString(bubble.begin[1]) << endl));
            start_indel_prediction();
        }
    }
}


bool BubbleFinder::close(const SNP_State& state) {
    bubble.closed_bubble=true;
    DEBUG((cout<<"last  node1.value "<<graph.toString(state.higher.last_edge.to)<<" node2.value "<<graph.toString(state.lower.last_edge.to)<<endl));
    /** We check the branching properties of the next kmers. */



    /** We finish the bubble with last distinct nodes. */
    bubble.end[0] = state.higher.last_edge.to;
    bubble.end[1] = state.lower.last_edge.to;

    checkPath();
    checkLowComplexity();
    /** We check several conditions (the first path vs. its revcomp and low complexity). */
    if (bubble.isCanonical && bubble.acceptable_complexity && checkRepeatSize(path_stack.first, path_stack.second))
    {

        /** We extend the bubble on the left and right (unitigs or contigs). */
        extend ();
        /** We got all the information about the bubble, we finish it. */
        bubble.extended_string[0] = path_stack.first;
        bubble.extended_string[1] = path_stack.second;
        bubble.final_nb_polymorphism=state.nb_polymorphism;
        finish ();
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
bool BubbleFinder::expand_heart( SNP_State state, // NOTE: Here the state is copied
                                 const pair<Edge, Edge>& next_edge){
    /** We check whether the new nodes are different from previous ones. */

    /************************************************************/
    /**                   RECURSION FINISHED                   **/
    /************************************************************/
    if(eqNode(next_edge.first.to, next_edge.second.to))
    { return close(state);
    }
    
    /************************************************************/
    /**                   RECURSION CONTINUES                  **/
    /************************************************************/
    else
    {
        

        DEBUG((cout<<"continue with nextNode1.value "<<graph.toString(state.higher.last_edge.to)<<" nextNode2.value "<<graph.toString(state.lower.last_edge.to)<<endl));
        /** We call recursively the method (recursion on 'pos'). */
        state.stack_size++;
        return state.extend(next_edge, path_stack) && expand(state);
        
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







bool BubbleFinder::expand (SNP_State& state) // NOTE: Here we mutate the parent state
{
    DEBUG((cout<<" *"<<state.stack_size<<"+"<<endl));
    DEBUG((cout<<"expand with node1.value "<<graph.toString(state.higher.last_edge.to)<<" node2.value "<<graph.toString(state.lower.last_edge.to)<<endl));
    DEBUG((cout<<"expand with local_extended_string1 "<<path_stack.first<<" local_extended_string2 "<<path_stack.second<<endl));


    /****************************************************************************/
    /**************** OPTIMIZATION : AVOID RECURSIONS ***************************/


    Graph::Vector < pair<Edge,Edge> > successors;
    while (true){
        if (checkBranching(state) == false) return false;       // no possibility to continue
        successors = graph.successorsEdge (state.higher.last_edge.to, state.lower.last_edge.to); // get next two nodes
        if (successors.size() != 1 || eqNode(successors[0].first.to, successors[0].second.to)) break;                                          // no more iterative programming
        assert(eqNode(successors[0].first.from, state.higher.last_edge.to));
        assert(eqNode(successors[0].second.from, state.lower.last_edge.to));
        if(!state.extend(successors[0], path_stack)) { return false ; }
    }

    /**************** END OPTIMIZATION  *****************************************/
    /****************************************************************************/

    /** We may have to stop the extension according to the branching mode. */
    if (checkBranching(state) == false)  {
        return false;
    }
    
    
    /** We get the common successors of node1 and node2. */
    /** Returns the successors of two nodes, ie with the same transition nucleotide from both nodes. */
//    Graph::Vector < pair<Node,Node> > successors = graph.successors (node1, node2);
    DEBUG((cout<<"successors size "<<successors.size()<<endl));


    
    bool dumped_bubble=false;
    /** We loop over the successors of the two nodes. */
    size_t i;
    for (i=0; i<successors.size(); i++)
    {

        assert(state.higher.last_edge.to == successors[i].first.from);
        assert(state.lower.last_edge.to == successors[i].second.from);





        /** extend the bubble with the couple of nodes */
        state.backtrack_stacks(path_stack);
        dumped_bubble = expand_heart(state, successors[i]);
        assert(!dumped_bubble || (dumped_bubble && bubble.closed_bubble));
        /** B 2 special case: if two or more symmetrical branching close a bubble, the output is redundant. **/
        /** Thus, if successors[i].first = successors[i].second and if the bubble is dumped, we stop **/
        if (successors[i].first.to == successors[i].second.to && dumped_bubble) break;
        

        // /** Stop as soon as a bubble is dumped */
        // VERSION 2.2.5: commented this break line. Enable to explore all possible symmetrical paths, even in case of success on one of the paths.
        //if(dumped_bubble) break;
    }
    
    DEBUG((cout<<"stop try"<<endl));
    if(dumped_bubble || bubble.closed_bubble) {
        return true;
    }
    
    /** NON DUMPED BUBBLE */
    /** if the bubble was closed and was not dumped, it means that it was not canonical. It will be find latter thus we should not find close SNPs from this bubble. */
    if(bubble.closed_bubble){
        return true;
    }
    
    /** Maybe we can search for a close SNP */
    if (state.nb_polymorphism < max_polymorphism && bubble.type==0) {
        DEBUG((cout<<"try with a new polymorphism ("<<state.nb_polymorphism<<") with node1.value "<<graph.toString(state.higher.last_edge.to)<<" node2.value "<<graph.toString(state.lower.last_edge.to)<<endl));
        Graph::Vector < Edge > successors1 = graph.successorsEdge (state.higher.last_edge.to);
        Graph::Vector < Edge > successors2 = graph.successorsEdge (state.lower.last_edge.to);
        
        /** We loop over the successors of the two nodes found with distinct extending nucleotides. */
        for (size_t i1=0; i1<successors1.size(); i1++){
            Edge& edge1 = successors1[i1];
            assert(edge1.from == state.higher.last_edge.to);
            for (size_t i2=0; i2<successors2.size(); i2++){
                Edge& edge2 = successors2[i2];
                assert(edge2.from == state.lower.last_edge.to);
                if ( edge1.nt == edge2.nt) continue; // This has already been tested in previous loop

                DEBUG((cout<<"TRYING"<<endl));

                state.nb_polymorphism++;
                state.backtrack_stacks(path_stack);
                dumped_bubble |= expand_heart(state, make_pair(edge1, edge2));
                assert(!dumped_bubble || (dumped_bubble && bubble.closed_bubble));
                /******************************************************************************************* **/
                /** Un-understood Sept 2015 (Pierre). Removed and replaced by the next "break"               **/
                /** if the bubble is finished with THIS couple of tested successors, we stop here.**/
                /** if we don't check this, in b 2 mode we may close a bubble with several distinct couple of node and thus create redondant bubbles **/
                /** if(dumped_bubble && successors1[i1]==successors2[i2]) break; **/
                /******************************************************************************************* **/
                if(dumped_bubble || bubble.closed_bubble) break;
            }
            if(dumped_bubble || bubble.closed_bubble) break;
        }
    }
    return dumped_bubble || bubble.closed_bubble;
}



/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS : In case of indel, the extended nodes are those from the full path of length 2k-1
 *********************************************************************/
void BubbleFinder::extend ()
{
    Nucleotide closureLeft  = NUCL_UNKNOWN;
    Nucleotide closureRight = NUCL_UNKNOWN;
    
    /** We may have to extend the bubble according to the user choice. */
    if (traversalKind != TRAVERSAL_NONE)
    {
        /** We ask for the predecessors of the first node and successors of the last node. */
        Graph::Vector<Node> successors   = graph.successors   (bubble.end[0]);
        Graph::Vector<Node> predecessors = graph.predecessors (bubble.begin[0]);
        
        /** We need to reset branching nodes between extensions in case of overlapping extensions. */
        _terminator->reset ();
        
        /** If unique, we keep the left/right extensions. */
        if (successors.size()==1)
        {
            /** We compute right extension of the node. */
            closureRight = graph.getNT (successors  [0], sizeKmer-1);
            _traversal->traverse (successors[0], DIR_OUTCOMING, bubble.extensionRight);
            bubble.divergenceRight = _traversal->getBubbles().empty() ? bubble.extensionRight.size() : _traversal->getBubbles()[0].first;
        }
        
        if (predecessors.size()==1)
        {
            /** We compute left extension of the node. */
            closureLeft  = graph.getNT (predecessors[0], 0);
            Node rev_pred = graph.reverse(predecessors[0]);
            _traversal->traverse (rev_pred, DIR_OUTCOMING, bubble.extensionLeft);
            bubble.divergenceLeft = _traversal->getBubbles().empty() ? bubble.extensionLeft.size() : _traversal->getBubbles()[0].first;
        }
    }
    
    /** We return a code value according to left/right extensions status. */
    if (closureLeft==NUCL_UNKNOWN && closureRight==NUCL_UNKNOWN)  { bubble.where_to_extend = 0; }
    else if (closureLeft!=NUCL_UNKNOWN && closureRight==NUCL_UNKNOWN)  { bubble.where_to_extend = 1; }
    else if (closureLeft==NUCL_UNKNOWN && closureRight!=NUCL_UNKNOWN)  { bubble.where_to_extend = 2; }
    else if (closureLeft!=NUCL_UNKNOWN && closureRight!=NUCL_UNKNOWN)  { bubble.where_to_extend = 3; }
    
    bubble.closureLeft  = closureLeft;
    bubble.closureRight = closureRight;
    
}
/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
void BubbleFinder::finish ()
{
    
    
    /** We build two Sequence objects from the information of the bubble. */
    
    /** We set the bubble index. NOTE: We have to make sure it is protected against concurrent
     * accesses since we may be called here from different threads. */
    if (bubble.type==0) __sync_add_and_fetch (&(stats.nb_bubbles_snp), 1);
    if (bubble.type==1) __sync_add_and_fetch (&(stats.nb_bubbles_del), 1);
    bubble.index = __sync_add_and_fetch (&(stats.nb_bubbles), 1);

    /** compute the bubble paths */
    string path_0 = graph.toString (bubble.begin[0])+bubble.extended_string[0];
    string path_1 = graph.toString (bubble.begin[1])+bubble.extended_string[1];
    stringstream comment;
    if ( bubble.polymorphism_type=="SNP" ){
        int polymorphism_id=1;
        for (unsigned int i=0;i<path_0.length();i++){
            if (path_0[i]!=path_1[i]) {
                if (polymorphism_id>1) {
                    comment << ",";
                }
                comment<<"P_" << polymorphism_id << ":" << i << "_" << path_0[i] << "/" << path_1[i];
                polymorphism_id++;
            }
        }
    }
    if ( bubble.polymorphism_type=="INDEL" ){
        const int insert_size = path_0.length()<path_1.length()?path_1.length()-path_0.length():path_0.length()-path_1.length();
        const int size_repeat = sizeKmer-2-min(bubble.extended_string[0].length(),bubble.extended_string[1].length()); // SEE checkRepeatSize function for explanations
        
            
        comment << "P_1:" << (sizeKmer-1) << "_" << (insert_size) << "_" << (size_repeat);
    }
    
    
    if (bubble.extended_string[0].length()<=bubble.extended_string[1].length()){
        buildSequence ( 0, "higher", bubble.seq1, comment.str());
        buildSequence ( 1, "lower",  bubble.seq2, comment.str());
    }
    else{ // put the smaller overlap as the first sequence.
        buildSequence ( 1, "higher", bubble.seq1, comment.str());
        buildSequence ( 0, "lower",  bubble.seq2, comment.str());
    }
    
    /** We have to protect the sequences dump wrt concurrent accesses. We use a {} block with
     * a LocalSynchronizer instance with the shared ISynchonizer of the Kissnp2 class. */
    {
        LocalSynchronizer sync (_synchronizer);
        
        /** We insert the two sequences into the output bank. */
        _outputBank->insert (bubble.seq1);
        _outputBank->insert (bubble.seq2);
        
        /** Stats update (in concurrent access protection block). */
        
//        stats.nb_bubbles++;
        if (bubble.type==0){
            stats.nb_where_to_extend_snp[bubble.where_to_extend] ++;
            if (bubble.high_complexity)  { stats.nb_bubbles_snp_high++; }
            else                         { stats.nb_bubbles_snp_low++;  }
        }
        if (bubble.type==1){
            stats.nb_where_to_extend_del[bubble.where_to_extend] ++;
            if (bubble.high_complexity)  { stats.nb_bubbles_del_high++; }
            else                         { stats.nb_bubbles_del_low++;  }
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
void BubbleFinder::buildSequence ( size_t pathIdx, const char* type, Sequence& seq, string polymorphism_comments)
{
    stringstream commentStream;
    
    /** We build the comment for the sequence. */
    commentStream << bubble.polymorphism_type << "_" << type << "_path_" << bubble.index << "|" << polymorphism_comments << "|" << (bubble.high_complexity ? "high" : "low")<< "|nb_pol_" <<bubble.final_nb_polymorphism;
    
    
    
    /** We may have extra information for the comment. */
    if (traversalKind == TRAVERSAL_UNITIG)
    {
        commentStream << "|left_unitig_length_";
        commentStream << (bubble.where_to_extend%2==1 ? (bubble.extensionLeft.size() +1) : 0);
        commentStream << "|right_unitig_length_";
        commentStream << (bubble.where_to_extend>1    ? (bubble.extensionRight.size()+1) : 0);
    }
    if (traversalKind == TRAVERSAL_CONTIG)
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
    
    if (bubble.closureLeft  != NUCL_UNKNOWN)  { len += 1 + lenLeft;  }
    if (bubble.closureRight != NUCL_UNKNOWN)  { len += 1 + lenRight; }
    
    /** We resize the sequence data if needed. Note: +1 for ending '\0'
     * NOTE: we use resize if we need more space, setSize otherwise. */
    if (seq.getData().size() < len+1)  {  seq.getData().resize  (len+1); }
    else                               {  seq.getData().setSize (len+1); }
    
    char* output = seq.getDataBuffer();
    
    /** We add the left extension if any. Note that we use lower case for extensions. */
    if (bubble.closureLeft != NUCL_UNKNOWN)
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
    if (bubble.closureRight != NUCL_UNKNOWN)
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
 ** REMARKS : Trick: we test the reverse complement of the second $k$-mer. Thus the last nuc. of the first kmer and the first nuc. of the second kmer does not influence the
 **           comparison. Thus the choice of the path used to make this comparison (higher or lower) does not change the results.
 *********************************************************************/
void BubbleFinder::checkPath ()
{
    /** We test whether the first kmer of the first path is smaller than
     * the first kmer of the revcomp(first path), this should avoid repeated SNPs */
    DEBUG((cout<<"check path "<<graph.toString (bubble.begin[0])  <<"<"<<  graph.toString (graph.reverse(bubble.end[0]))<<endl));
    if(graph.toString (bubble.begin[0])  <  graph.toString (graph.reverse(bubble.end[0])))
        bubble.isCanonical=true;
    else
        bubble.isCanonical=false;
    
    return ;
}



/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
bool BubbleFinder::checkRepeatSize (string &extension1, string &extension2) const
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

    
    const int size_repeat = sizeKmer-2-min(extension1.length(), extension2.length());
    if (size_repeat>max_indel_ambiguity) {
        return false;
    }
    return true;
    
}


inline bool BubbleFinder::isAuthorizedSymmetricBranching(size_t& sym_branches) const {
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
 ** METHOD  : BubbleFinder::checkBranching
 ** PURPOSE : Checks the branching properties of a couple of nodes. If authorised_branching==0: no branching is authorized in any of the two nodes.
 **           If authorised_branching==1: no symetrical branching authorized from the two nodes (ie. N1-->A and N1-->B and N2-->A and N2-->B not authorized)
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
bool BubbleFinder::checkBranching (SNP_State& state) const
{
    if(authorised_branching==0) {
        return !graph.isBranching(state.higher.last_edge.to) && !graph.isBranching(state.lower.last_edge.to);
    } else {
        // Analyze succesors :
        Graph::Vector< pair<Edge, Edge> > successors = graph.successorsEdge(state.higher.last_edge.to, state.lower.last_edge.to);
        size_t symmetric_successors = 0, ends = 0;
        for (size_t i=0; i < successors.size(); i++) {
            pair<Edge,Edge> succ = successors[i];
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
        Graph::Vector< pair<Node, Node> > predecessors = graph.successors(graph.reverse(state.higher.last_edge.to), graph.reverse(state.lower.last_edge.to));
        size_t symmetric_predecessors = 0, starts = 0;
        for (size_t i=0; i < predecessors.size(); i++) {
            pair<Node,Node> pred = predecessors[i];
            if (eqNode(pred.first, pred.second)) {
                if(state.higher.local_extension_length != 0 && state.lower.local_extension_length != 0)
                {
                    cout << "Found alternative begin past original begin " << graph.toString(graph.reverse(pred.first)) << " for bubble : " << endl
                         << "\t'" << path_stack.first << "' ("<< graph.toString(state.higher.last_edge.to) <<  ")" << endl
                         << "\t'" << path_stack.second << "' ("<< graph.toString(state.lower.last_edge.to) <<  ")" << endl;
                }
                starts++; // TODO: abbort if not cannonical (ie. will be visited with another start node)
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
void BubbleFinder::checkLowComplexity ()
{
    bubble.acceptable_complexity=true;
    string path1 = graph.toString (bubble.begin[0]).substr(0, sizeKmer-1) + graph.toString (bubble.end[0]);
    string path2 = graph.toString (bubble.begin[1]).substr(0, sizeKmer-1) + graph.toString (bubble.end[1]);
    
    /** We compute the low complexity score of the two paths. */
    bubble.high_complexity = filterLowComplexity2Paths (path1, path2);
    
    if (accept_low) return; // the complexity is acceptable for this bubble anyway.
    
    bubble.acceptable_complexity=bubble.high_complexity;
    DEBUG((cout<<"check low "<<accept_low<<" "<<bubble.high_complexity<<endl));
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
IProperties* BubbleFinder::getConfig () const
{
    IProperties* props = new Properties();
    
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

