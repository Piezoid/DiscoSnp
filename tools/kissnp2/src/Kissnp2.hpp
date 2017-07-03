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

#ifndef _TOOL_KISSNP2_HPP_
#define _TOOL_KISSNP2_HPP_

/********************************************************************************/
// #include <gatb/gatb_core.hpp>
#include <gatb/tools/misc/impl/Tool.hpp>
#include <string>


// using namespace gatb::core::tools::misc::impl;
// using namespace gatb::core::tools::misc;


using namespace std;



/********************************************************************************/

/** \brief Tool class that looks for SNP
 *
 * The Kissnp2 is the front class for SNP detection in a provided de Bruijn graph.
 * The output is a bank with pairs of sequences defining a bubble.
 */
class Kissnp2 : public gatb::core::tools::misc::impl::Tool
{
public:

    /** Constructor. */
    Kissnp2 ();
    
    /** Implementation of Tool::execute method. */
    void execute ();

    size_t getKmerSize() const;

    const string& getGraphUri() const { return graph_uri; }
private:
    string graph_uri;
};

/********************************************************************************/

#endif /* _TOOL_KISSNP2_HPP_ */

