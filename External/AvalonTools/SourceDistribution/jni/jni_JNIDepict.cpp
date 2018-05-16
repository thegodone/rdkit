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
/*
 * File: JNIDepict.cpp
 *
 * Implements the wrappers for the JNIDepict class.
 */

#include "jni_JNIDepict.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "local.h"
#include "reaccs.h"
#include "reaccsio.h"
#include "utilities.h"
#include "smi2mol.h"
#include "depictutil.h"
#include "ssmatch.h"
#include "patclean.h"

/* Referred to in local.h but not used */
FILE *log_file = (FILE *)NULL;

/* Values between 0 and 999 enable heap debugging completely */
/* Values above 1000 enable the segment that has the same thousends digit */
#define DEBUG_HEAP 0000

/**
 * Returns a full path name for a file fname in the directory pointed
 * to by the environment variables TMP or TEMP in that order.
 */
static char *GetTempName(const char *fname)
{
   static char buffer[1000];
   char *varp;
   varp = getenv("TMP");
   if (!varp) varp = getenv("TEMP");
   if (varp)
   {
      strcpy(buffer, varp);
      strcat(buffer, "\\");
   }
   else      strcpy(buffer, "C:\\TEMP\\");
   strcat(buffer, fname);
   return (buffer);
}

void logString(char *str)
{
   FILE *fp;
   fp = fopen(GetTempName("tmp.log"), "a");
   fprintf(fp, "%s\n", str);
   fclose(fp);
}


/*
 * Class:     JNIDepict
 * Method:    smilesMatchesQueryCTNative
 * Signature: (Ljava/lang/String;Ljava/lang/String;)Z
 */
JNIEXPORT
jboolean JNICALL Java_jni_JNIDepict_smilesMatchesQueryCTNative(JNIEnv *env,
                                                               jobject thisObj,
                                                               jstring smilesObj,
                                                               jstring queryCT)
{
   struct reaccs_molecule_t *mp;
   struct reaccs_molecule_t *qp;
   const char *const_smiles;
   const char *const_query;
   jboolean isCopy = JNI_FALSE;
   jboolean result = JNI_FALSE;
   Fortran_FILE *ffp;
   ssmatch_t * matches;
   struct reaccs_bond_t *bp;
   unsigned int i;
   int *H_count;


   const_smiles = env->GetStringUTFChars(smilesObj, &isCopy);
   mp = SMIToMOL(const_smiles, DO_LAYOUT | DROP_TRIVIAL_HYDROGENS);
   if (isCopy == JNI_TRUE) env->ReleaseStringUTFChars(smilesObj, const_smiles);
   /* Illegal molecules don't match at all */
   if (mp == (struct reaccs_molecule_t *)NULL)
      return (result);

   // Read query
   const_query = env->GetStringUTFChars(queryCT, &isCopy);
   ffp = FortranStringOpen((char *)const_query);
   qp = TypeAlloc(1,struct reaccs_molecule_t);
   if (FORTRAN_NORMAL != ReadREACCSMolecule(ffp,qp,""))
   {
      fprintf(stderr, "failed to read template\n");
      FreeMolecule(qp);
      qp = (struct reaccs_molecule_t *)NULL;
   }
   FortranClose(ffp);
   if (isCopy == JNI_TRUE) env->ReleaseStringUTFChars(queryCT, const_query);


   if (qp == (struct reaccs_molecule_t *)NULL)  /* illegal query => no match */
   {
      FreeMolecule(mp);
      return (result);
   }

   /* Compute query_H_count to match the required explicit hydrogens */
   MakeHydrogensImplicit(qp);

   /* Set up hydrogen count fields in structure for matching */
   H_count = TypeAlloc(mp->n_atoms+1, int);
   ComputeImplicitH(mp, H_count);
   /* Add the explicit hydrogens to the implicit counts */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      if (0 == strcmp("H", mp->atom_array[bp->atoms[0]-1].atom_symbol))
         H_count[bp->atoms[1]]++;
      else if (0 == strcmp("H", mp->atom_array[bp->atoms[1]-1].atom_symbol))
         H_count[bp->atoms[0]]++;
   }
   /* set the 'query_H_count' field to the correct value */
   for (i=0; i<mp->n_atoms; i++)
      if (H_count[i+1] >= 0) 
         mp->atom_array[i].query_H_count = ZERO_COUNT+H_count[i+1];




   matches = SSMatch(mp, qp, TRUE);
   if (matches != (ssmatch_t *)NULL)
   {
      result = JNI_TRUE;
      FreeSSMatch(matches);
   }
   FreeMolecule(mp); FreeMolecule(qp);


   MyFree((char *)H_count);
   return (result);
}


/*
 * Class:     JNIDepict
 * Method:    smilesToMOLFile
 * Signature: (Ljava/lang/String;Ljava/lang/String;)I
 *
 * Converts the SMILES string *smiles to a MOL-file named *fname.
 *
 * The file is cleared in case of error.
 */
JNIEXPORT
jint JNICALL Java_jni_JNIDepict_smilesToMOLFile(JNIEnv *env,
                                                jobject thisObj,
                                                jstring smilesObj,
                                                jstring fnameObj)
{
   int result = TRUE;
   FILE *fp;
   struct reaccs_molecule_t *mp;
   unsigned int i;
   const char *const_smiles;
   const char *const_fname;
   jboolean isCopy = JNI_FALSE;


   const_smiles = env->GetStringUTFChars(smilesObj, &isCopy);
   mp = SMIToMOL(const_smiles, DO_LAYOUT | DROP_TRIVIAL_HYDROGENS);
   if (isCopy == JNI_TRUE) env->ReleaseStringUTFChars(smilesObj, const_smiles);

   const_fname = env->GetStringUTFChars(fnameObj, &isCopy);
   fp = fopen(const_fname,"w");
   if (isCopy == JNI_TRUE) env->ReleaseStringUTFChars(fnameObj, const_fname);

   /* The following is a patch to get the correct sizes into ISIS */
   if (mp)
   {
      for (i=0; i<mp->n_atoms; i++)
      {
         mp->atom_array[i].x *= 0.5;
         mp->atom_array[i].y *= 0.5;
      }

      PrintREACCSMolecule(fp,mp,"");
      FreeMolecule(mp);
   }
   else
      result = FALSE;
   fclose(fp);
   return (result);
}

/*
 * Class:     JNIDepict
 * Method:    smilesToMOLFileWithTemplate
 * Signature: (Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;)I
 *
 * This function writes the MOL file corresponding to smiles to the
 * file fname using the connection table in tplCT to rotate it.
 *
 * It returns 0 if the operation was successful and non-0 otherwise.
 */
JNIEXPORT
jint JNICALL Java_jni_JNIDepict_smilesToMOLFileWithTemplate(JNIEnv *env,
                                                            jobject thisObj,
                                                            jstring smilesObj,
                                                            jstring fnameObj,
                                                            jstring tplCTObj)
{
   int result = TRUE;
   FILE *fp;
   struct reaccs_molecule_t *mp;
   unsigned int i, itmp;
   const char *const_smiles;
   const char *const_fname;
   jboolean isCopy = JNI_FALSE;

   const char *const_tplCT;
   Fortran_FILE *ffp;
   struct reaccs_molecule_t *tpl;


   // Read template
   const_tplCT = env->GetStringUTFChars(tplCTObj, &isCopy);
   ffp = FortranStringOpen((char *)const_tplCT);
   tpl = TypeAlloc(1,struct reaccs_molecule_t);
   if (FORTRAN_NORMAL != ReadREACCSMolecule(ffp,tpl,""))
   {
      fprintf(stderr, "failed to read template\n");
      FreeMolecule(tpl);
      tpl = (struct reaccs_molecule_t *)NULL;
   }
   else if (tpl->n_atoms <= 0)
   {
      fprintf(stderr, "template without atoms\n");
      FreeMolecule(tpl);
      tpl = (struct reaccs_molecule_t *)NULL;
   }
   else
   {
      MakeHydrogensImplicit(tpl); /* templates don't recognize hydrogens */
      for (i=0; i<tpl->n_atoms; i++)
         tpl->atom_array[i].query_H_count = NONE;
   }
   FortranClose(ffp);
   if (isCopy == JNI_TRUE) env->ReleaseStringUTFChars(tplCTObj, const_tplCT);

   const_smiles = env->GetStringUTFChars(smilesObj, &isCopy);
   mp = SMIToMOL(const_smiles, DO_LAYOUT | DROP_TRIVIAL_HYDROGENS);
   if (isCopy == JNI_TRUE) env->ReleaseStringUTFChars(smilesObj, const_smiles);

   if (tpl != (struct reaccs_molecule_t *)NULL)
   {
      itmp = TemplateRotate(mp, tpl);
   }

   const_fname = env->GetStringUTFChars(fnameObj, &isCopy);
   fp = fopen(const_fname,"w");
   if (isCopy == JNI_TRUE) env->ReleaseStringUTFChars(fnameObj, const_fname);

   /* The following is a patch to get the correct sizes into ISIS */
   if (mp)
   {
      for (i=0; i<mp->n_atoms; i++)
      {
         mp->atom_array[i].x *= 0.5;
         mp->atom_array[i].y *= 0.5;
      }

      PrintREACCSMolecule(fp,mp,"");
      FreeMolecule(mp);
   }
   else
      result = FALSE;
   fclose(fp);
   if (tpl != (struct reaccs_molecule_t *)NULL) FreeMolecule(tpl);
   return (result);
}

/*
 * Class:     JNIDepict
 * Method:    mwFromSmiles
 * Signature: (Ljava/lang/String;)D
 */
JNIEXPORT
jdouble JNICALL Java_jni_JNIDepict_mwFromSmiles(JNIEnv  *env,
                                                jobject  thisObj,
                                                jstring smilesObj)
{
   char buffer[100];
   double mw;
   const char *const_smiles;
   char smiles[1000];
   jboolean isCopy = JNI_FALSE;
   

   const_smiles = env->GetStringUTFChars(smilesObj, &isCopy);
   strcpy(smiles, const_smiles);
   SmilesToMWMF(smiles, &mw, buffer);
   if (isCopy == JNI_TRUE) env->ReleaseStringUTFChars(smilesObj, const_smiles);

   return (mw);
}

/*
 * Class:     JNIDepict
 * Method:    MOLFileToSmilesBytes
 * Signature: ([BLjava/lang/String;)I
 */
JNIEXPORT jint JNICALL Java_jni_JNIDepict_MOLFileToSmilesBytes
  (JNIEnv *env, jobject thisObj, jbyteArray buffer, jstring fnameObj)
{
   const char *const_fname;
   int size;
   char *smiles;
   jbyte *bytes;
   jboolean isCopy = JNI_FALSE;


   const_fname = env->GetStringUTFChars(fnameObj, &isCopy);
   size = 999;
   smiles = MOLFileToSMILESString(&size, (char*)const_fname);
   if (isCopy == JNI_TRUE) env->ReleaseStringUTFChars(fnameObj, const_fname);

   if (size <= 0) return (0);
   size = strlen(smiles);
   if (size+1 >= env->GetArrayLength(buffer)) return (-size);

   bytes = env->GetByteArrayElements(buffer, &isCopy);
   strncpy((char *)bytes, smiles, size);
   bytes[size] = '\0';
   if (isCopy == JNI_TRUE) env->ReleaseByteArrayElements(buffer, bytes, 0);
   return (size);
}

/*
 * Class:     JNIDepict
 * Method:    MOLFileToSMARTSBytes
 * Signature: ([BLjava/lang/String;)I
 */
JNIEXPORT jint JNICALL Java_jni_JNIDepict_MOLFileToSMARTSBytes
  (JNIEnv *env, jobject thisObj, jbyteArray buffer, jstring fnameObj)
{
   const char *const_fname;
   int size;
   char *smarts;
   jbyte *bytes;
   jboolean isCopy = JNI_FALSE;


   const_fname = env->GetStringUTFChars(fnameObj, &isCopy);
   size = 999;
   smarts = MOLFileToSMARTSString(&size, (char*)const_fname);
   if (isCopy == JNI_TRUE) env->ReleaseStringUTFChars(fnameObj, const_fname);

   if (smarts == 0  ||  size <= 0) return (0);
   size = strlen(smarts);
   if (size+1 >= env->GetArrayLength(buffer)) return (-size);

   bytes = env->GetByteArrayElements(buffer, &isCopy);
   strncpy((char *)bytes, smarts, size);
   bytes[size] = '\0';
   if (isCopy == JNI_TRUE) env->ReleaseByteArrayElements(buffer, bytes, 0);
   return (size);
}

#ifdef __cplusplus
}
#endif

