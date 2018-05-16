//
//  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met: 
//
//     * Redistributions of source code must retain the above copyright 
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following 
//       disclaimer in the documentation and/or other materials provided 
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
//       nor the names of its contributors may be used to endorse or promote 
//       products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
/************************************************************************/
/*                                                                      */
/*    File:           layout.h                                          */
/*                                                                      */
/*    Purpose:        This file declares the funtions provided for      */
/*                    layouting connection tables. The functions can    */
/*                    layout molecules as well as the components of     */
/*                    reactions.                                        */
/*                                                                      */
/************************************************************************/

                    /* used to mark bonds not to be used for layout */
#define RUBBER_BOND     0x40
#define DONT_FLIP_BOND  0x20
	/* used to prevent atom from being initialized to random coordinates */
#define KEEP_POSITION   0x4000

struct reaccs_molecule_t *LayoutMolecule(struct reaccs_molecule_t *mp);
/*
 * Imposes a new layout onto the molecule *mp. It assumes that the
 * atoms and bonds of prelayouted fragments have been colored with
 * the same color.
 * The function returns a copy of the molecule with new coordinates.
 * The original molecule is not changed.
 */

void RecolorMolecule(struct reaccs_molecule_t *mp);
/*
 * Recolors the atoms and bonds of the molecule *mp, each one
 * with a different color.
 */

void LayoutBondStereo(struct reaccs_molecule_t *mp,
		      neighbourhood_t *nbp,
		      int ring_size[]);
/*
 * Uses the stereoparity information of the bonds of *mp
 * to flip double bonds if needed.
 */

int FloodColor(struct reaccs_molecule_t *mp,
               neighbourhood_t *nbp,
	           int aindex, int color);
/*
 * Recursively colors the fragment of *mp containing the atom with
 * index aindex with color. It uses a flood fill algorithm.
 * It returns the number of atoms recolored during the process.
 */

