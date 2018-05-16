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
/**********************************************************************/
/*                                                                    */
/*    File:           canonizer.c                                     */
/*                                                                    */
/*    Purpose:        This file implements the functions needed to    */
/*                    generate canonical SMILES from an MDL molecule  */
/*                    data structure.                                 */
/*                                                                    */
/*    History:        07-Feb-2002     Start of development.           */
/*                                                                    */
/**********************************************************************/

#ifdef __WIN32__
#include <windows.h>
#endif

#include <stdlib.h>

#include <string.h>
#include <ctype.h>

#include "reaccs.h"
#include "utilities.h"
#include "perceive.h"
#include "reaccsio.h"
#include "stereo.h"
#include "pattern.h"

#include "smi2mol.h"
#include "canonizer.h"

static int CheckHeap(void)
{
   int result = TRUE;

/* non-dummy only for 32-bit Windows */
#ifdef __WIN32__
   HGLOBAL heap;

   heap = GetProcessHeap();
   result = HeapValidate(heap,0,NULL);
   if (!result)
      fprintf(stderr, "HeapValidate(%d, %d, %d) returns %d\n", heap, 0, NULL, result);
#endif

   return (result);
}

static int PropertySplit(int n_atoms,
                         int atom_index[],
                         int pre_num[],
                         long split_prop[])
/*
 * atom_index[i] is the index of the atom that has rank i in the emerging
 * canonical ordering. pre_num[j] is the canonical class number for atom
 * j. split_prop[j] is the property of atom j on which the next class split
 * is going to happen.
 *
 * The function returns TRUE if pre_num[] has been refined and FALSE otherwise.
 */
{
   int *new_num;
   int i,j, tmp, last_class;
   int result = FALSE;

   /* Insertion sort of indexes */
   for (i=1; i<n_atoms; i++)
      for (j=i-1; j>=0; j--)
         if (pre_num[atom_index[j]] == pre_num[atom_index[j+1]]  &&
             split_prop[atom_index[j]] > split_prop[atom_index[j+1]])
         {
            tmp = atom_index[j];
            atom_index[j] = atom_index[j+1];
            atom_index[j+1] = tmp;
         }
         else
            break;
   /* update pre_numbering */
   last_class = 0;
   new_num = TypeAlloc(n_atoms, int);
   new_num[atom_index[0]] = last_class;
   for (i=1; i<n_atoms; i++)
      if (pre_num[atom_index[i-1]] != pre_num[atom_index[i]]  ||
          split_prop[atom_index[i-1]] != split_prop[atom_index[i]])
      {
         new_num[atom_index[i]] = i;
         last_class = i;
      }
      else
      {
         new_num[atom_index[i]] = last_class;
      }
   for (i=0; i<n_atoms; i++)
   {
      if ( pre_num[i] != new_num[i]) result = TRUE;
      pre_num[i] = new_num[i];
   }
   MyFree((char *)new_num);

   return (result);
}

static void ComputeNewGraphProperties(struct reaccs_molecule_t *mp,
                                      neighbourhood_t           nbp[],
                                      long                      split_prop[],
                                      int                       pre_num[])
/*
 * Refine properties using neighbourhood.
 */
{
   struct reaccs_atom_t *ap;
   int i, j, bt;
   unsigned int sum, prod;

   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      split_prop[i] = 0;
      sum = 0; prod = 1;
      for (j=0; j<nbp[i].n_ligands; j++)
      {
         bt = mp->bond_array[nbp[i].bonds[j]].bond_type;
         sum  += pre_num[nbp[i].atoms[j]]*bt;
         prod *= pre_num[nbp[i].atoms[j]]+bt;
         prod &= 0xFFFF;
      }
      split_prop[i] = 0xFFFF&(sum*0x100 + prod);
   }
}

typedef double point[2];

static void print_prenum(char *header, int *pre_num, struct reaccs_molecule_t *mp)
/*
 * Debug output
 */
{
   int i;
   fprintf(stderr, "%s\t", header);
   for (i=0; i<mp->n_atoms; i++)
      if (mp->atom_array[i].query_H_count)
         fprintf(stderr, " %s%d|%d(%d),", mp->atom_array[i].atom_symbol, i+1, pre_num[i], mp->atom_array[i].query_H_count);
      else
         fprintf(stderr, " %s%d|%d,", mp->atom_array[i].atom_symbol, i+1, pre_num[i]);
   fprintf(stderr, "\n");
}

static int *CanonicalNumbering(struct reaccs_molecule_t *mp,
                               neighbourhood_t *nbp,
                               int flags)
/*
 * Returns a pointer to an array of canonical numbers[0..mp->n_atoms-1].
 * The smallest number being 0. It uses refinement on atom and bond types
 * and also on stereochemistry.
 *
 * Superflous explicit hydrogens need to be removed before this function
 * is called.
 *
 * Note: As a side effect, this method orders the neighbourhood description
 * nbp by canonical rank.
 */
{
   int *atom_index;  /* atom indices sorted by rank */
   int *result;      /* resulting numbering. inverse of atom_index */
   long *split_prop; /* array of current property values */
   int *pre_num;     /* intermediate steps towards canonical numbering */
   int i, j, jj, tmp, ai, ai1, ai2;
   int l1, l2;
   point p1, p2, r12;
   double q;
   int changed;
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bp;
   int nequal;
   struct npoint_t tetra[4];
   double vol;
   int split_in_ring;
   int *atom_status, *bond_status;
   int min_tie, min_tie_index;
   int rgnum;

   if (mp->n_atoms <= 0) return TypeAlloc(0, int);      /* safe-guard */
   result = TypeAlloc(mp->n_atoms, int);

   atom_index = TypeAlloc(mp->n_atoms, int);
   for (i=0; i<mp->n_atoms; i++) atom_index[i] = i;
   pre_num = TypeAlloc(mp->n_atoms, int);
   for (i=0; i<mp->n_atoms; i++) pre_num[i] = 0;


   atom_status = TypeAlloc(mp->n_atoms, int);
   bond_status = TypeAlloc(mp->n_bonds, int);
   RingState(mp, atom_status, bond_status);

   split_prop = TypeAlloc(mp->n_atoms, long);

   /* atomic number and degree are first splitting criterion */
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      /* Prefer low ligance => to simplify DB processing */
      split_prop[i] = nbp[i].n_ligands;
      split_prop[i] *= 0x100;
      tmp = StringToInt(periodic_table, ap->atom_symbol);
      /* Move special atom types to front */
      if (tmp <= 114) split_prop[i] += tmp;
      else            split_prop[i] -= tmp;
      /* monomers need special treatment to make sure that */
      /* strings don't start with attachment point */
      if (ap->atom_symbol[0] == 'R'  &&
          ap->atom_symbol[1] == '#'  &&
          (flags & TO_MONOMER) != 0)
         split_prop[i] += 0x4000;
   }
   changed = PropertySplit(mp->n_atoms, atom_index, pre_num, split_prop);

   SetRingSizeFlags(mp, 12, nbp);
   /* Now, we split on ring sizes. */
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      split_prop[i] = ap->rsize_flags&0x3FFF;
   }
   changed = PropertySplit(mp->n_atoms, atom_index, pre_num, split_prop);

   /* Use other atom-bound properties, i.e. query_H_count, charges, radicals, and mass */
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      split_prop[i] = ap->mass_difference;
      /* hydrogen isotopes are special */
      if (0 == strcmp("D", ap->atom_symbol)) split_prop[i] = 1;
      if (0 == strcmp("T", ap->atom_symbol)) split_prop[i] = 2;
      split_prop[i] *= 0x04;
      split_prop[i] += ap->charge;
      split_prop[i] *= 0x04;
      split_prop[i] += ap->radical;
      split_prop[i] *= 0x04;
      split_prop[i] += ap->query_H_count;
      /* Special split for R-Group attachments */
      if (ap->atom_symbol[0] == 'R'  && ap->atom_symbol[1] == '#')
      {
         if (!GetNumProperty(mp->prop_lines, "M  RGP", i+1, &rgnum))
            if (ap->mass_difference > 0)
               rgnum = ap->mass_difference;
            else
               rgnum = 0;
         split_prop[i] += 0x1000*rgnum;
      }
   }
   changed = PropertySplit(mp->n_atoms, atom_index, pre_num, split_prop);

next_cycle:
   do
   {
      ComputeNewGraphProperties(mp, nbp, split_prop, pre_num);

      changed = PropertySplit(mp->n_atoms, atom_index, pre_num, split_prop);
   } while (changed);

   /* Now, we user the neighbourhood to refine the splits. */

   /* compute inverted atom_index permutation as result */
   for (i=0; i<mp->n_atoms; i++)
      result[atom_index[i]] = i;

   /* sorting the neighbourhood */
   for (i=0; i<mp->n_atoms; i++)
   {
      for (j=1; j<nbp[i].n_ligands; j++)
         for (jj=j-1; jj>=0; jj--)
            if (result[nbp[i].atoms[jj]] > result[nbp[i].atoms[jj+1]])
            {
               tmp=nbp[i].atoms[jj];
               nbp[i].atoms[jj] = nbp[i].atoms[jj+1];
               nbp[i].atoms[jj+1] = tmp;
               tmp=nbp[i].bonds[jj];
               nbp[i].bonds[jj] = nbp[i].bonds[jj+1];
               nbp[i].bonds[jj+1] = tmp;
            }
            else
               break;
   }

   /* refine by stereo parity */
   for (i=0; i<mp->n_atoms; i++)
   {
      split_prop[i] = 0;
      if (nbp[i].n_ligands < 3  ||  nbp[i].n_ligands > 4) continue;
      /* Check if two ligands have same renumbering */
      nequal = 0;
      for (j=1; j<nbp[i].n_ligands; j++)
         if (pre_num[nbp[i].atoms[j-1]] == pre_num[nbp[i].atoms[j]])
            nequal++;
      if (nequal > 0) continue; /* Cannot be TH stereo center */
      for (j=0; j<nbp[i].n_ligands; j++)
      {
         tetra[j].x = mp->atom_array[nbp[i].atoms[j]].x;
         tetra[j].y = mp->atom_array[nbp[i].atoms[j]].y;
         tetra[j].z = 0;
         bp = &mp->bond_array[nbp[i].bonds[j]];
         if (bp->atoms[0]-1 != i) continue;  /* not the cusp of that bond */
         if (bp->stereo_symbol == UP)   tetra[j].z =  1.0;
         if (bp->stereo_symbol == DOWN) tetra[j].z = -1.0;
      }
      if (nbp[i].n_ligands == 3)
      {
         tetra[3].x = mp->atom_array[i].x;
         tetra[3].y = mp->atom_array[i].y;
         tetra[3].z = 0;
      }
      vol = Volume(tetra);
      if (vol >  0.0001) split_prop[i] =  1;
      if (vol < -0.0001) split_prop[i] = -1;
   }
   changed = PropertySplit(mp->n_atoms, atom_index, pre_num, split_prop);
   if (changed) goto next_cycle;

   /* sorting the neighbourhood */
   for (i=0; i<mp->n_atoms; i++)
   {
      for (j=1; j<nbp[i].n_ligands; j++)
         for (jj=j-1; jj>=0; jj--)
            if (result[nbp[i].atoms[jj]] > result[nbp[i].atoms[jj+1]])
            {
               tmp=nbp[i].atoms[jj];
               nbp[i].atoms[jj] = nbp[i].atoms[jj+1];
               nbp[i].atoms[jj+1] = tmp;
               tmp=nbp[i].bonds[jj];
               nbp[i].bonds[jj] = nbp[i].bonds[jj+1];
               nbp[i].bonds[jj+1] = tmp;
            }
            else
               break;
   }

   /* refine ring centers by stereo ligand classes */
   for (i=0; i<mp->n_atoms; i++) split_prop[i] = 0;
   for (i=0; i<mp->n_atoms; i++)
   {
      if (nbp[i].n_ligands < 3  ||  nbp[i].n_ligands > 4) continue;
      /* Check if two ligands have same renumbering */
      nequal = 0;
      for (j=1; j<nbp[i].n_ligands; j++)
         if (pre_num[nbp[i].atoms[j-1]] == pre_num[nbp[i].atoms[j]])
            nequal++;
      if (nequal != 1) continue; /* Only one tie can be split */
      for (j=0; j<nbp[i].n_ligands; j++)
      {
         tetra[j].x = mp->atom_array[nbp[i].atoms[j]].x;
         tetra[j].y = mp->atom_array[nbp[i].atoms[j]].y;
         tetra[j].z = 0;
         bp = &mp->bond_array[nbp[i].bonds[j]];
         if (bp->atoms[0]-1 != i) continue;  /* not the cusp of that bond */
         if (bp->stereo_symbol == UP)   tetra[j].z =  1.0;
         if (bp->stereo_symbol == DOWN) tetra[j].z = -1.0;
      }
      if (nbp[i].n_ligands == 3)
      {
         tetra[3].x = mp->atom_array[i].x;
         tetra[3].y = mp->atom_array[i].y;
         tetra[3].z = 0;
      }
      vol = Volume(tetra);
      if (vol >  0.0001)      vol =  1.0;
      else if (vol < -0.0001) vol = -1.0;
      else                    continue;
      /* break tie depending on vol */
      split_in_ring = FALSE;
      for (j=1; j<nbp[i].n_ligands; j++)
         if (pre_num[nbp[i].atoms[j-1]] == pre_num[nbp[i].atoms[j]])
         {
            /* Only split non-ring pairs */
            if (atom_status[nbp[i].atoms[j]] != 0)
               split_in_ring = TRUE;
            if (atom_status[nbp[i].atoms[j-1]] != 0)
               split_in_ring = TRUE;
            if (!split_in_ring) break;
            split_prop[nbp[i].atoms[j-1]] +=  vol;
            split_prop[nbp[i].atoms[j]]   += -vol;
            break;
         }
   }
   changed = PropertySplit(mp->n_atoms, atom_index, pre_num, split_prop);
   if (changed) goto next_cycle;

   /* Clear stereo parity for remaining unresolved centers */
   for (i=0; i<mp->n_atoms; i++)
   {
      if (nbp[i].n_ligands < 3  ||  nbp[i].n_ligands > 4) continue;
      /* Check if two ligands have same renumbering */
      nequal = 0;
      for (j=1; j<nbp[i].n_ligands; j++)
         if (pre_num[nbp[i].atoms[j-1]] == pre_num[nbp[i].atoms[j]])
            nequal++;
      if (nequal == 0) continue;          /* resolved => no change */
      for (j=0; j<nbp[i].n_ligands; j++)
      {
         bp = &mp->bond_array[nbp[i].bonds[j]];
         if (bp->atoms[0]-1 != i) continue;  /* not the cusp of that bond */
         if (bp->stereo_symbol == UP)   bp->stereo_symbol = NONE;
         if (bp->stereo_symbol == DOWN) bp->stereo_symbol = NONE;
      }
   }

   /* Clear DB stereo with symmetrical substitution on at least one end */
   /* nbp[] is assumed to be sorted before */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      if (bp->bond_type != DOUBLE) continue;
      nequal = 0;
      ai = bp->atoms[0]-1;
      for (j=1; j<nbp[ai].n_ligands; j++)
         if (pre_num[nbp[ai].atoms[j-1]] ==
             pre_num[nbp[ai].atoms[j]])
            nequal++;
      ai = bp->atoms[1]-1;
      for (j=1; j<nbp[ai].n_ligands; j++)
         if (pre_num[nbp[ai].atoms[j-1]] ==
             pre_num[nbp[ai].atoms[j]])
            nequal++;
      if (nequal > 0  &&  bp->stereo_symbol == NONE)
         bp->stereo_symbol = CIS_TRANS_EITHER;
   }

   /* Split on CIS/TRANS properties */
   for (i=0; i<mp->n_atoms; i++) split_prop[i] = 0;
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      if (bp->bond_type != DOUBLE) continue;
      if (bp->stereo_symbol == CIS_TRANS_EITHER) continue;
      ai1 = bp->atoms[0]-1; ai2 = bp->atoms[1]-1;
      /* find reference ligands, i.e. the ones with smallest index */
      l1 = mp->n_atoms;
      for (j=0; j<nbp[ai1].n_ligands; j++)
	 if (nbp[ai1].atoms[j] != ai2)
         {
	    l1 = nbp[ai1].atoms[j];
            break;
         }
      if (l1 == mp->n_atoms) continue;
      l2 = mp->n_atoms;
      for (j=0; j<nbp[ai2].n_ligands; j++)
	 if (nbp[ai2].atoms[j] != ai1)
         {
	    l2 = nbp[ai2].atoms[j];
            break;
         }
      if (l2 == mp->n_atoms) continue;
      p1[0] = mp->atom_array[l1].x - mp->atom_array[ai1].x;
      p1[1] = mp->atom_array[l1].y - mp->atom_array[ai1].y;
      p2[0] = mp->atom_array[l2].x - mp->atom_array[ai2].x;
      p2[1] = mp->atom_array[l2].y - mp->atom_array[ai2].y;
      r12[0]= mp->atom_array[ai2].x - mp->atom_array[ai1].x;
      r12[1]= mp->atom_array[ai2].y - mp->atom_array[ai1].y;
      q = (p1[0]*r12[1]-p1[1]*r12[0]) * (p2[0]*r12[1]-p2[1]*r12[0]);
      if (bp->stereo_symbol != CIS_TRANS_SWAPPED) q *= -1;
      if (q < 0.00001)
      {
         split_prop[ai1] -= 1; split_prop[ai2] -= 1;
      }
      else if (q > 0.00001)
      {
         split_prop[ai1] += 1; split_prop[ai2] += 1;
      }
   }
   changed = PropertySplit(mp->n_atoms, atom_index, pre_num, split_prop);
   if (changed) goto next_cycle;

   /* Find a remaining tie if any */
   min_tie = mp->n_atoms;
   min_tie_index = -1;
   for (i=0; i<mp->n_atoms; i++)
   {
      if (min_tie <= pre_num[i]) continue;  /* cannot be a smaller tie */
      for (j=i+1; j<mp->n_atoms; j++)
         if (pre_num[i] == pre_num[j])
            break;
      if (j == mp->n_atoms) continue;       /* no tie */
      min_tie = pre_num[i];
      min_tie_index = i;
   }
   if (min_tie_index >= 0)                  /* do the split */
   {
      for (i=0; i<mp->n_atoms; i++)
      {
         if (pre_num[i] == min_tie  &&  i != min_tie_index)
            split_prop[i] = 1;
         else
            split_prop[i] = 0;
      }
      changed = PropertySplit(mp->n_atoms, atom_index, pre_num, split_prop);
      if (changed) goto next_cycle;
   }

   MyFree((char *)pre_num);
   MyFree((char *)split_prop);
   MyFree((char *)atom_index);
   MyFree((char *)atom_status);
   MyFree((char *)bond_status);

   return result;
}

static unsigned int seed=13;

void ScrambleMolecule(struct reaccs_molecule_t *mp)
/*
 * Renumber the molecule mp and consequently re-sort its atoms. This function is used to provide
 * test cases to show that canonical coding is independant of input order.
 */
{
   int *numbering;
   int i, itmp, i1, i2;
   struct reaccs_atom_t atom;

   srand(seed);
   seed = seed*4+3;
   numbering = TypeAlloc(mp->n_atoms, int);
   for (i=0; i<mp->n_atoms; i++)
      numbering[i] = i;
   for (i=0; i<2*mp->n_atoms; i++)
   {
      i1 = (int)(mp->n_atoms*(rand()/(double)RAND_MAX));
      i2 = (int)(mp->n_atoms*(rand()/(double)RAND_MAX));
      itmp = numbering[i1];
      numbering[i1] = numbering[i2];
      numbering[i2] = itmp;
      atom = mp->atom_array[numbering[i1]];
      mp->atom_array[numbering[i1]] = mp->atom_array[numbering[i2]];
      mp->atom_array[numbering[i2]] = atom;
   }
   for (i=0; i<mp->n_bonds; i++)
   {
      mp->bond_array[i].atoms[0] =
         numbering[mp->bond_array[i].atoms[0]-1]+1;
      mp->bond_array[i].atoms[1] =
         numbering[mp->bond_array[i].atoms[1]-1]+1;
   }
   MyFree((char *)numbering);
}

int SplitIonicBonds(struct reaccs_molecule_t *mp)
/*
 * This function converts hetero to main group metal bonds to the charge
 * separated no-bond state.
 * This is done to get a more effective conversion to the standard fragments.
 *
 * The function returns TRUE when it deleted a bond and FALSE otherwise.
 */
{
   int i, ai1, ai2;
   struct reaccs_bond_t *bp;
   int bond_deleted;
   int *good_atoms, *good_bonds;

   bond_deleted = FALSE;
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      if (bp->bond_type != SINGLE) continue;
      ai1 = bp->atoms[0]-1; ai2 = bp->atoms[1]-1;
      if (AtomSymbolMatch(mp->atom_array[ai1].atom_symbol, "alk,gr2")  &&
          AtomSymbolMatch(mp->atom_array[ai2].atom_symbol, "Q"))
      {
         bp->bond_type = -1;    // Mark as to be deleted
         mp->atom_array[ai1].charge++;
         mp->atom_array[ai2].charge--;
         bond_deleted = TRUE;
      }
      else if (AtomSymbolMatch(mp->atom_array[ai2].atom_symbol, "alk,gr2")  &&
               AtomSymbolMatch(mp->atom_array[ai1].atom_symbol, "Q"))
      {
         bp->bond_type = -1;    // Mark as to be deleted
         mp->atom_array[ai2].charge++;
         mp->atom_array[ai1].charge--;
         bond_deleted = TRUE;
      }
   }
   if (!bond_deleted) return (FALSE);
   good_atoms = TypeAlloc(mp->n_atoms+1, int);
   good_bonds = TypeAlloc(mp->n_bonds, int);
   for (i=0; i<mp->n_atoms; i++)
      good_atoms[i+1] = TRUE;
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      good_bonds[i] = bp->bond_type >= 0;
   StripMolecule(mp, good_atoms, good_bonds);
   MyFree((char *)good_atoms);
   MyFree((char *)good_bonds);
   return (TRUE);
}

static char *tfm_strings[][2] =
{
   {"N(-C,N,O)(=C,N,O)(=O:1)",    "N+1(-C,N,O)(=C,N,O)(-O-1)"},
   {"N(=C,N,O)(#C,N:1)",    "N+1(=C,N,O)(=C,N-1)"},
   {"B(-A)(-A)(-A)(-A)",    "B-1(-A)(-A)(-A)(-A)"},
};

/*
 * This array is filled with the parsed augmented atom transformation rules
 * built from the string pairs above.
 */
static aa_pair *trans_pairs = (aa_pair *)NULL;

static void ParseTfmTable(void)
{
   int i, npairs;

   if (trans_pairs) return;
   npairs = sizeof(tfm_strings)/sizeof(char *)/2;
   trans_pairs = TypeAlloc(npairs, aa_pair);
   for (i=0; i<sizeof(tfm_strings)/sizeof(char *)/2; i++)
   {
      StringToAugmentedAtom(tfm_strings[i][0], &trans_pairs[i][0]);
      StringToAugmentedAtom(tfm_strings[i][1], &trans_pairs[i][1]);
   }
}

/**
 * This function makes sure that all atom environments of first row elements
 * are converted such that they obey the octet rule. The most obvious case
 * is the nitro group. Non-first row elements have their bipolar bonds
 * converted to double bonds.
 * The function also protonates any unprotonated negatively charged atoms that
 * are not a member of a bipolar bond.
 */
void ForceOctet(struct reaccs_molecule_t *mp)
{
   int i, ai1, ai2, tmp;
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bp;
   int *H_count;

   /* prepare some cross references */
   H_count = TypeAlloc(mp->n_atoms+1, int);
   ComputeImplicitH(mp, H_count);
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      if (ap->query_H_count != NONE)
      {
         H_count[i+1] = ap->query_H_count-ZERO_COUNT;
      }
   }

   /* First, we convert all potentially bipolar bonds into double bonds */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      ai1 = bp->atoms[0]-1; ai2 = bp->atoms[1]-1;
      if (mp->atom_array[ai1].charge <= 0  &&  mp->atom_array[ai2].charge <= 0)
         continue;
      /* Now, we have at least one positive charge */
      if (mp->atom_array[ai2].charge > 0)
      {
         tmp = ai1; ai1 = ai2; ai2 = tmp;
      }
      /* Now, ai1 points to the positive atom and ai2 to the other one */
      if (mp->atom_array[ai2].charge < 0)  /* charges can be used to make DB */
      {
         if (bp->bond_type == SINGLE)
            bp->bond_type = DOUBLE;
         else if (bp->bond_type == DOUBLE)
            bp->bond_type = TRIPLE;
         else
            continue;
         mp->atom_array[ai1].charge--;
         mp->atom_array[ai2].charge++;
      }
      else if (H_count[ai2+1] > 0  &&
               AtomSymbolMatch(mp->atom_array[ai2].atom_symbol, "Q"))
      {
         mp->atom_array[ai1].charge--;
         H_count[ai2+1]--;
         bp->bond_type = DOUBLE;
      }
   }

   /* Uncharge any further charged non-metals */
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      if (ap->charge == 0) continue;
      if (!AtomSymbolMatch(ap->atom_symbol, "Q")) continue;

      /* Remove negative charge by protonation, if neighbour isn't positive */
      if (ap->charge < 0   &&  AtomSymbolMatch(ap->atom_symbol, "Q") )
      {
         if (ap->query_H_count == NONE)
         {
             H_count[i+1] -= ap->charge;
             ap->charge = 0;
         }
      }
      /* Remove positive charge by deprotonation as much as possible */
      else
      {
         tmp = ap->charge;
         if (H_count[i+1] < ap->charge) tmp = H_count[i+1];
         H_count[i+1] -= tmp;
         ap->charge -= tmp;
      }
      if (ap->query_H_count != NONE)
      {
         ap->query_H_count = H_count[i+1]+ZERO_COUNT;
      }
   }

   /* Now, we go through the atoms and transform the ones that need */
   /* bipolar bonds to obey the Octet-Rule. */
   ParseTfmTable();
   TransformAugmentedAtoms(mp, trans_pairs, sizeof(tfm_strings)/sizeof(char *)/2);

   MyFree((char *)H_count);
}

/**
 * This function removes all connected fully inorganic components from *mp.
 * It retains either the organic part(s) or the largest inorganic part.
 */
void RemoveInorganicFragments(struct reaccs_molecule_t *mp)
{
   int i;
   int ncarbon;
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bp;
   int *good_atoms, *good_bonds;

   /* Reset all colors to 0 */
   ncarbon = 0;
   ResetColors(mp);
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (strcmp(ap->atom_symbol, "C") == 0)    /* Carbon atom => label it */
      {
         ncarbon++;
         ap->color = 1;
      }
   if (ncarbon > 0)    /* Only delete inorgancs if there is a carbon atom */
   {
      FloodFillFragments(mp);
      good_atoms = TypeAlloc(mp->n_atoms+1, int);
      good_bonds = TypeAlloc(mp->n_bonds, int);
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
         if (ap->color == 1) good_atoms[i+1] = TRUE;
         else                good_atoms[i+1] = FALSE;
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
         if (mp->atom_array[bp->atoms[0]-1].color != NONE  ||
             mp->atom_array[bp->atoms[1]-1].color != NONE)
            good_bonds[i] = TRUE;
         else
            good_bonds[i] = FALSE;
      StripMolecule(mp, good_atoms, good_bonds);
      MyFree((char *)good_atoms);
      MyFree((char *)good_bonds);
   }
   ResetColors(mp);
}

/**
 * This function removes the SMILES of any remaining standard counter ions
 * or solvent molecules from smiles[] and returnes the modified string.
 * The table of 'salt_array' is provides as a parameter and assumed to
 * be already canonicalized.
 * Note: The function modifies smiles in situ.
 */
char *RemoveStandardComponents(char *smiles, char **salt_array)
{
   int i, j, ncopy, nsalt;
   int changed;
   char *match;
   char salt[100];

   /* catch some wierd input */
   if (!smiles  ||
       !smiles[0]  ||
       !salt_array ||
       !salt_array[0])
      return (smiles);

   do
   {
      changed = FALSE;
      for (i=0; salt_array[i]; i++)
      {
         /* short circuit obvious non-matches */
         if (!strstr(smiles, salt_array[i])) continue;
         /* squeeze out salt_array at all possible places */
         /* try beginning of smiles */
         strcpy(salt, salt_array[i]); strcat(salt, ".");
         nsalt = strlen(salt);
         match = strstr(smiles, salt);
         if (smiles == match)   /* smiles starts with salt */
         {
            ncopy = strlen(match)-nsalt;
            for (j=0; j<=ncopy; j++) match[j] = match[nsalt+j];
            changed = TRUE;
            continue;
         }
         /* try middle of smiles */
         strcpy(salt, "."); strcat(salt, salt_array[i]); strcat(salt, ".");
         nsalt = strlen(salt);
         match = strstr(smiles, salt);
         if (match)   /* smiles contains salt */
         {
            ncopy = strlen(match)-nsalt;
            match++;    /* keep first '.' */
            for (j=0; j<=ncopy; j++) match[j] = match[nsalt+j-1];
            changed = TRUE;
            continue;
         }
         /* try end of smiles */
         strcpy(salt, "."); strcat(salt, salt_array[i]);
         nsalt = strlen(salt);
         match = strstr(smiles, salt);
         if (match  &&  strlen(match) == nsalt)   /* smiles ends with salt */
         {
            match[0] = '\0';    /* truncate at match */
            changed = TRUE;
            continue;
         }
      }
   } while (changed);

   return (smiles);
}

/*
 * Canonization primitive function. It returns a pointer to a newly
 * allocated character string that contains the canonicalized form of
 * insmiles[]. It is used once or twice by CanSmiles(). The second call
 * is only done in special case such as the existence of stereo double bonds.
 *
 * It will return NULL if anything goes wrong.
 */
char *CanSmilesStep(char *insmiles, int flags)
{
   char *smiles, *cp, *cph;
   char *coordp;
   int *numbering;
   struct reaccs_molecule_t *mp;
   neighbourhood_t *nbp;
   int i;
   struct reaccs_atom_t *ap;

   char from[20], to[20];

   mp = SMIToMOL(insmiles, TRUE_DB_STEREO | DO_LAYOUT | DROP_TRIVIAL_HYDROGENS);
   if (!mp) return NULL;
   if (!(flags&CENTER_STEREO))  /* remove hashes and wedges */
   {
      for (i=0; i<mp->n_bonds; i++)
         if (mp->bond_array[i].bond_type == SINGLE)
            mp->bond_array[i].stereo_symbol = NONE;
   }

   if (flags&REMOVE_ISOTOPE)
      for (i=0; i<mp->n_atoms; i++)
         mp->atom_array[i].mass_difference = 0;

   if (flags&REMOVE_MAPPING)
      for (i=0; i<mp->n_atoms; i++)
         mp->atom_array[i].mapping = NONE;

   if (flags & (COMPONENT_SET|MAIN_PART))
   {
      SplitIonicBonds(mp);
   }

   if (flags & FORCE_OCTET) ForceOctet(mp);
   if (flags & MAIN_PART) RemoveInorganicFragments(mp);

   if (flags&SCRAMBLE) ScrambleMolecule(mp);

   if (!(flags&TO_MONOMER))
      MakeSomeHydrogensImplicit(mp,
                                NON_STEREO_HYDROGENS |
                                NO_QUERY_H_COUNT	|
                                ANCILLARY_STEREO);

   nbp   = TypeAlloc(mp->n_atoms, neighbourhood_t);
   if (!SetupNeighbourhood(mp,nbp,mp->n_atoms))
   {
      MyFree((char *)nbp);
      return NULL;
   }

   if (flags&DY_AROMATICITY)
      PerceiveDYAromaticity(mp, nbp);
   else
      PerceiveAromaticBonds(mp);

   numbering = CanonicalNumbering(mp, nbp, flags);
   if (flags & TO_MONOMER)
   {
      /* Convert simple R atoms to ones with a mass */
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
         if (0 == strcmp(ap->atom_symbol,"R")) ap->charge = 1;
      smiles = MOLToSMIExt(mp,
                           TO_MONOMER |
                           ISOMERIC_SMILES |
                           AROMATIC_SMILES,
                           numbering,
                           &coordp);
      ReplaceOnce(smiles,"([*+])", "&1");
      ReplaceOnce(smiles,"[*+]", "&1");
      for (i=2; i<20; i++)
      {
         sprintf(to, "&%d", i);
         sprintf(from, "([*+%d])", i); ReplaceOnce(smiles, from, to);
         sprintf(from, "[*+%d]", i);   ReplaceOnce(smiles, from, to);
      }
   }
   else
      smiles = MOLToSMIExt(mp,
                           ISOMERIC_SMILES |
                           AROMATIC_SMILES,
                           numbering,
                           &coordp);
   FreeMolecule(mp);
   if (!(flags&DB_STEREO))      /* prune away DB stereo designations */
   {
      cp=cph=smiles;
      while (*cp)
         if ((*cp)=='/'  ||  (*cp)=='\\')
            cp++;
         else
         {
            (*cph) = (*cp);
            cph++; cp++;
         }
      (*cph) = '\0';
   }
   if (coordp) MyFree((char *)coordp);
   if (numbering) MyFree((char *)numbering);
   MyFree((char *)nbp);
   return (smiles);
}

#define START_OF_FRAGMENT       1
#define SLASH_FOUND             2
#define BACKSLASH_FOUND         3

void UnifyFragmentDBStereo(char *smiles)
/*
 * Converts the '/' and '\'characters in smiles such that each fragment
 * (separted by '.') has '/' first. It uses a finite state machine approach.
 */
{
   int state = START_OF_FRAGMENT;
   char *cp;

   if (!smiles) return;
   for (cp=smiles; (*cp) != '\0'; cp++)
   {
      switch (*cp)
      {
         case '\\':
            if (state == START_OF_FRAGMENT)
            {
               state = BACKSLASH_FOUND;
               (*cp) = '/';
            }
            else if (state == BACKSLASH_FOUND)
               (*cp) = '/';
            else if (state == SLASH_FOUND)
            {
               /* NOP since fragment started with '/' */;
            }
            else
               fprintf(stderr, "1: illegal state %d\n", state);
            break;
         case '/':
            if (state == START_OF_FRAGMENT)
               state = SLASH_FOUND;
            else if (state == BACKSLASH_FOUND)
               (*cp) = '\\';
            else if (state == SLASH_FOUND)
            {
               /* NOP since fragment started with '/' */
            }
            else
               fprintf(stderr, "1: illegal state %d\n", state);
            break;
         case '.':
            state = START_OF_FRAGMENT;
            break;
         default:
            break;
      }
   }
}

#define MAXLINE 20000

/**
 * Parses smiles into unique '.'-separated components. The function returns
 * the modified smiles which has repeated fragments removed. Note, that this
 * function assumes that the input smiles is canonical and thus has identical
 * fragments listed next to each other.
 */
char *makeComponentSet(char *smiles)
{
   char *tok;
   char last_fragment[MAXLINE];
   char result[MAXLINE];

   /* safeguard against dubious input */
   if (!smiles  ||  smiles[0] == '\0') return (smiles);

   result[0] = '\0';
   last_fragment[0] = '\0';
   for (tok = strtok(smiles, ".");
        tok != (char *)NULL;
        tok = strtok((char*) NULL, "."))
   {
      if (last_fragment[0])     /* subsequent fragments */
      {
         if (0 != strcmp(last_fragment, tok))   /* new fragment? */
         {
            strcpy(last_fragment, tok);
            strcat(result, ".");
            strcat(result, tok);
         }
      }
      else                      /* first fragment */
      {
         strcpy(result, tok);
         strcpy(last_fragment, tok);
      }
   }
   strcpy(smiles, result);      /* can only get shorter than original smiles */
   return (smiles);
}

/*
 * Return a pointer to a newly allocated character string that contains
 * the canonicalized form of insmiles[]. Delegates actual processing to
 * CanSmilesStep().
 *
 * It will return NULL if anything goes wrong.
 */
char *CanSmiles(char *insmiles, int flags)
{
   char *smiles, *tmp;

if (insmiles)
   smiles = CanSmilesStep(insmiles, flags);
   if (smiles == NULL) return (NULL);
   UnifyFragmentDBStereo(smiles);
   /* Rerun canonicalization if there is a stereo double bond */
   if (!(flags&TO_MONOMER)  &&
       (strchr(smiles, '/')  ||  strchr(smiles, '\\')))
   {
      tmp = smiles;
      smiles = CanSmilesStep(tmp, flags & (~MAIN_PART));
      UnifyFragmentDBStereo(smiles);
      MyFree((char *)tmp);
   }
   if (flags&COMPONENT_SET) smiles = makeComponentSet(smiles);
   return (smiles);
}

#ifdef MAIN

static void explain(char * prognam, int argc)
{
    if (strstr(prognam, "canonizer"))
   {
      fprintf(stderr,"canonizer used with %d arguments\n",argc-1);
      fprintf(stderr,"usage: canonizer {flags} <input.tbl >output.tbl\n");
      fprintf(stderr,"The following flags are known: \n");
      fprintf(stderr,"-h      \tShows this list \n");
      fprintf(stderr,"-hasid  \tfirst token is an id and appended to output\n");
      fprintf(stderr,"-add    \tappends the result to the input line\n");
      fprintf(stderr,"-unique \talso outputs the SMILES w/o stereo\n");
      fprintf(stderr,"-monomer\tconvert the SMILES to a monomer definition\n");
      fprintf(stderr,"-keep_mapping\tdoes not clear away atom-mapping information\n");
      fprintf(stderr,"-mwmf   \tcalculate MW and MF\n");
      fprintf(stderr,"-d \tuse Daylight aromaticity perception\n");
      fprintf(stderr,"-do_test \tperform a regression test \n");
      exit (EXIT_FAILURE);
   }
}

#include "didepict.h"

/* used for 16-bit Borland Turbo C compilation */
#ifdef __TURBOC__
#include <process.h>
unsigned _stklen = 0xFEEE;
#endif

#define MXBUF	MAXLINE

FILE *log_file;

main(int argc, char *argv[])
{
   struct reaccs_molecule_t *mp;
   char *smiles, *usmiles;
   char *coordp;
   char line[MAXLINE];
   char buffer[MXBUF];
   char check_smiles[MXBUF];
   char *cp;
   char prognam[256];
   int nsmi, i;
   double mw;
   char mf[100];

   int show_help    = FALSE;
   int add_result   = FALSE;
   int unique_smiles= FALSE;
   int first_run, nruns;
   int do_test      = FALSE;
   int mwmf         = FALSE;
   int hasid        = FALSE;
   char *smistart;
   char id[32];
   int flags        = 0;

   log_file = stderr;

   flags |= REMOVE_MAPPING;
   /* Read options */
   while (argc > 1  &&  argv[1][0] == '-' && argv[1][1] != '\0' )
   {
      if (0 == strcmp(argv[1], "-add"))      add_result = TRUE;
      if (0 == strcmp(argv[1], "-unique"))   unique_smiles = TRUE;
      if (0 == strcmp(argv[1], "-do_test"))  do_test   = TRUE;
      if (0 == strcmp(argv[1], "-hasid"))    hasid   = TRUE;
      if (0 == strcmp(argv[1], "-mwmf"))     mwmf   = TRUE;
      if (0 == strcmp(argv[1], "-d"))       flags   |= DY_AROMATICITY;
      if (0 == strcmp(argv[1], "-monomer"))  flags   |= TO_MONOMER;
      if (0 == strcmp(argv[1], "-keep_mapping"))  flags   &= ~REMOVE_MAPPING;
      if (0 == strcmp(argv[1], "-h"))        show_help  = TRUE;
      for (i=2; i<argc; i++) argv[i-1] = argv[i];
      argc--;
   }

   coordp = NULL;
   strncpy(prognam, argv[0], 255);
   for (cp=prognam; (*cp); cp++)
      (*cp) = tolower(*cp);

   if( show_help )
   {
      explain(prognam, argc);
      return (EXIT_FAILURE);
   }

   /* Stand-alone driver for canonicalization */
   if (strstr(prognam, "canonizer"))
   {
      nsmi = 0;
      SetStripZeroes(TRUE);
      while (!feof(stdin))
      {
         fgets(line, MAXLINE-1, stdin);
if (strstr(line, "EXIT")) break;
         /* remove any remaining '\n' */
         line[MAXLINE-1] = '\0';
         if (strlen(line) > 0 && line[strlen(line)-1] == '\n')
            line[strlen(line)-1] = '\0';
	 if (hasid)
	 {
	    for (cp=line; (*cp) && !isspace(*cp); cp++)
	       ;
	    strncpy(id, line, cp-line);
	    id[cp-line] = '\0';
	    while (isspace(*cp)) cp++;
	    smistart = cp;
	 }
	 else
	    smistart = line;
         /* get SMILES as first word in line */
         for (cp=smistart; (*cp) && !isspace(*cp); cp++)
            ;
         strncpy(buffer, smistart, cp-smistart);
         buffer[cp-smistart] = '\0';
         if (feof(stdin)) break;
         if (strlen(buffer) == 0) continue;
         first_run = TRUE;
         nruns = 0;
         // if (do_test) fprintf(stderr, "%s\n", buffer);

         seed = 13;

repeat_do_test:
         if (do_test)
            smiles = CanSmiles(buffer,
                               SCRAMBLE |
                               flags |
                               DB_STEREO |
                               CENTER_STEREO);
         else
            smiles = CanSmiles(buffer,
                               flags |
                               DB_STEREO |
                               CENTER_STEREO);
         if (unique_smiles)
            usmiles = CanSmiles(buffer, flags);
         else
            usmiles=(char *)NULL;
         if (!smiles)
         {
            /* No result => just print line */
            if (first_run  &&  add_result) fprintf(stdout,"%s\n", line);
            continue;
         }
         nruns++;
         if (first_run)
         {
            nsmi++;
            strcpy(check_smiles, smiles);
         }

         if (add_result)
         {
            if (first_run)
            {
               fprintf(stdout, "%s\t%s", smistart, smiles);
               if (unique_smiles) fprintf(stdout, "\t%s", usmiles);
               fprintf(stdout, "\n");
            }
            else
            {
               if (0 != strcmp(check_smiles, smiles))
               {
                  fprintf(stderr, "'%s' <> '%s'\n", check_smiles, smiles);
                  fprintf(stderr, "original line = %s\n", line);
                  // exit(1);
               }
            }
         }
         else
         {
            if (first_run)
            {
               if (mwmf)
               {
                  mp = SMIToMOL(smiles, DROP_TRIVIAL_HYDROGENS);
                  if (!mp) fprintf(stdout,"%s\n", smistart);
                  else
                  {
                     fprintf(stdout,"%s", smistart);
                     mw = MolecularWeight(mp);
                     MolecularFormula(mp, mf);
                     fprintf(stdout,"\t%g\t%s",mw,mf);
                     FreeMolecule(mp);
                     if (unique_smiles) fprintf(stdout, "\t%s", usmiles);
		     if (hasid)
			fprintf(stdout, "\t%s\n", id);
		     else
			fprintf(stdout, "\n");
                  }
               }
               else
               {
                  fprintf(stdout, "%s", smiles);
                  if (unique_smiles) fprintf(stdout, "\t%s", usmiles);
		  if (hasid)
		     fprintf(stdout, "\t%s\n", id);
		  else
		     fprintf(stdout, "\n");
               }
            }
            else
            {
               if (0 != strcmp(check_smiles, smiles))
               {
                  fprintf(stderr, "'%s' <> '%s'\n", check_smiles, smiles);
                  fprintf(stderr, "original line = %s\n", line);
                  // exit(1);
               }
            }
         }
         first_run=FALSE;
         MyFree(smiles);
         if (usmiles) MyFree(usmiles);
	 if (coordp)
	 {
	    MyFree((char *)coordp);
	    coordp = NULL;
	 }
         if (nruns == 2)   /* check for oscillations */
         {
            strcpy(buffer, check_smiles);
         }
         if (do_test && nruns < 3) goto repeat_do_test;
      }
      return (EXIT_SUCCESS);
   }

   else
   {
      explain(prognam, argc);
   }

   return (EXIT_SUCCESS);
}

#endif
