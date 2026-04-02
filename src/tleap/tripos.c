/*
 *  File:   tripos.c
 *
 ************************************************************************
 *                            LEAP                                      *
 *                                                                      *
 *                   Copyright (c) 1992, 1995                           *
 *           Regents of the University of California                    *
 *                     All Rights Reserved.                             *
 *                                                                      *
 *  This software provided pursuant to a license agreement containing   *
 *  restrictions on its disclosure, duplication, and use. This software *
 *  contains confidential and proprietary information, and may not be   *
 *  extracted or distributed, in whole or in part, for any purpose      *
 *  whatsoever, without the express written permission of the authors.  *
 *  This notice, and the associated author list, must be attached to    *
 *  all copies, or extracts, of this software. Any additional           *
 *  restrictions set forth in the license agreement also apply to this  *
 *  software.                                                           *
 ************************************************************************
 *                                                                      *
 *     Designed by:    Christian Schafmeister                           *
 *     Author:         Christian Schafmeister                           *
 *                                                                      *
 *     VERSION: 1.0                                                     *
 *     Programmers:                                                     *
 *             Christian Schafmeister                                   *
 *             David Rivkin                                             *
 *                                                                      *
 *     Principal Investigator: Peter A. Kollman                         *
 *                                                                      *
 ************************************************************************
 *
 *  Description:
 *      Code to read TRIPOS files from SYBL.
 */


#include    "basics.h"
#include    "classes.h"
#include    "tripos.h"
#include    "parmLib.h"

/*
 *  uTriposReadUnit
 *
 *  Author: Christian Schafmeister (1991)
 *
 *  Read a single UNIT from the TRIPOS file.
 *  Return the UNIT, or NULL if nothing was read.
 */
UNIT uTriposReadUnit(fIn)
FILE *fIn;
{
    STRING sLine;
    int iAtoms;
    VARARRAY vaAtoms;
    VARARRAY vaResidues;
    UNIT uUnit;
    STRING sUnitName;
    int iRet;
    double dX, dY, dZ, dDummy;
    STRING sType;
    int iResidue;
    STRING sTemp;
    VECTOR vPos;
    int iIndex, iA, iB;
    ATOM aA, aB;
    int iOrder;
    ATOM aAtom;
    BOOL bGotIt;
    int iFileAtoms;
    int iFileBonds;
    int iFileSubstructures;
    int iTemp, i;
    double dCharge;
    RESIDUE rRes;
    STRING sName;
    STRING sOrder;
    int iRootAtom;
    STRING sSize;
    STRING sDescriptor, sDesc;
    STRING sCmd;
    int iElement, iTag, iDummy;
    PARMSET psSet;


#define TRIPOS_MOLECULE     "@<TRIPOS>MOLECULE"
#define TRIPOS_ATOM     "@<TRIPOS>ATOM"
#define TRIPOS_BOND     "@<TRIPOS>BOND"
#define TRIPOS_SUBSTRUCTURE "@<TRIPOS>SUBSTRUCTURE"

#define T_FSCANF( iret,f,s,ss,fail) { fgets( s, sizeof(s), fIn );\
                      iret=sscanf ss;\
                      if (feof(f)) goto fail;}

    vaAtoms = NULL;
    vaResidues = NULL;

    /* Search for the TRIPOS_MOLECULE string */

    bGotIt = FALSE;
    while (TRUE) {
        T_FSCANF(iRet, fIn, sLine, (sLine, "%s", sCmd), FAIL);
        if (strcmp(sCmd, TRIPOS_MOLECULE) == 0) {
            bGotIt = TRUE;
            break;
        }
    }
    if (!bGotIt)
        goto FAIL;

    /* Read the UNITs name */

    uUnit = (UNIT) oCreate(UNITid);
    T_FSCANF(iRet, fIn, sLine, (sLine, "%s", sUnitName), FAIL);
    VP0 ( ( "Reading MOLECULE named %s", sLine ) ) ;
    ContainerSetName(uUnit, sUnitName);
    MESSAGE(("Reading unit: %s\n", sUnitName));
    T_FSCANF(iRet, fIn, sLine, (sLine, "%d %d %d %d %d", &iFileAtoms,
                                &iFileBonds, &iFileSubstructures,
                                &iTemp, &iTemp), FAIL);

    /* Search for the TRIPOS_ATOM string */

    bGotIt = FALSE;
    while (TRUE) {
        T_FSCANF(iRet, fIn, sLine, (sLine, "%s", sCmd), FAIL);
        if (strcmp(sCmd, TRIPOS_ATOM) == 0) {
            bGotIt = TRUE;
            break;
        }
    }
    if (!bGotIt)
        goto FAIL;

#if 0
    T_FSCANF(iRet, fIn, sLine, (sLine, "%s", sSize), FAIL);
    T_FSCANF(iRet, fIn, sLine, (sLine, "%s", sDescriptor), FAIL);
    T_FSCANF(iRet, fIn, sLine, (sLine, "%s", sTemp), FAIL);
    T_FSCANF(iRet, fIn, sLine, (sLine, "%s", sTemp), FAIL);

    T_FSCANF(iRet, fIn, sLine, (sLine, "%s", sCmd), FAIL);
    if (strcmp(sCmd, TRIPOS_ATOM) != 0)
        goto FAIL;
#endif

    /* Read in the ATOM records */

    vaAtoms = vaVarArrayCreate(sizeof(ATOM));
    for (i = 0; i < iFileAtoms; i++) {
        dCharge = 0.0;
        T_FSCANF(iRet, fIn, sLine,
                 (sLine, "%d %s %lf %lf %lf %s %d %s %lf",
                  &iIndex, sName,
                  &dX, &dY, &dZ, sType, &iResidue, sTemp, &dCharge), FAIL);
        if (iRet != 9)
            break;

        MESSAGE((" Atom: %s\n", sName));
        aAtom = (ATOM) oCreate(ATOMid);
        ContainerSetName(aAtom, sName);
        AtomSetTempInt(aAtom, iResidue);
        VectorDef(&vPos, dX, dY, dZ);
        AtomSetPosition(aAtom, vPos);

/*  iElement = iElementNumberFromAmber(sType);   */
/*    dac change: just use the first letter: how else to easily tell
      CA from calcium???                                              */

        sTemp[0] = cUpper(sType[0]);
        sTemp[1] = '\0';
        iElement = iElementNumberFromAmber(sTemp);

#if 0
/*  initial attempt to use the atom map, but how do we get that read in?  */
        if (bParmLibDefaultExists()) {
            fprintf(stderr, "have the library\n");
            PARMLIB_DEFAULT_LOOP(psSet,
                                 (iTag = iParmSetFindAtom(psSet, sType)));
            if (iTag != PARM_NOT_FOUND) {
                ParmSetAtom(psSet, iTag, sTemp,
                            &dDummy, &dDummy, &dDummy, &dDummy, &dDummy,
                            &dDummy, &iElement, &iDummy, sDesc);
            }
        }
#endif

        AtomSetElement(aAtom, iElement);
        AtomSetType(aAtom, sType);
        AtomSetCharge(aAtom, dCharge);
        VarArrayAdd(vaAtoms, (GENP) & aAtom);
    }

    T_FSCANF(iRet, fIn, sLine, (sLine, "%s", sCmd), FAIL);
    if (strcmp(sCmd, TRIPOS_BOND) != 0)
        goto FAIL;

    iAtoms = iVarArrayElementCount(vaAtoms);
    for (i = 0; i < iFileBonds; i++) {
        T_FSCANF(iRet, fIn, sLine,
                 (sLine, "%d %d %d %s", &iIndex, &iA, &iB, sOrder), FAIL);
        if (iRet != 4)
            break;
        MESSAGE((" Bond %d - %d\n", iA, iB));
        if (iA > iAtoms || iB > iAtoms) {
            printf("Cannot form bond between atoms %d and %d\n", iA, iB);
        } else {
            aA = *PVAI(vaAtoms, ATOM, iA - 1);
            aB = *PVAI(vaAtoms, ATOM, iB - 1);
            StringLower(sOrder);
            iOrder = BONDSINGLE;
            if (strcmp(sOrder, "1") == 0) {
                iOrder = BONDSINGLE;
            } else if (strcmp(sOrder, "2") == 0) {
                iOrder = BONDDOUBLE;
            } else if (strcmp(sOrder, "3") == 0) {
                iOrder = BONDTRIPLE;
            } else if (strcmp(sOrder, "ar") == 0) {
                iOrder = BONDAROMATIC;
            } else if (strcmp(sOrder, "am") == 0) {
                iOrder = BONDAROMATIC;
            }
            AtomBondToOrder(aA, aB, iOrder);
        }
    }

    /* Read the RESIDUE stuff */

    if( iFileSubstructures ){
        T_FSCANF(iRet, fIn, sLine, (sLine, "%s", sCmd), FAIL);
        if (strcmp(sCmd, TRIPOS_SUBSTRUCTURE) != 0) goto FAIL;

        vaResidues = vaVarArrayCreate(sizeof(RESIDUE));
        for (i = 0; i < iFileSubstructures; i++) {
            T_FSCANF(iRet, fIn, sLine,
                 (sLine, "%d %s %d %s",
                  &iIndex, sName, &iRootAtom, sType), FAIL);
            if (iRet != 4) break;

            MESSAGE((" Substructure: %s\n", sName));
            rRes = (RESIDUE) oCreate(RESIDUEid);
            ContainerSetName(rRes, sName);
            VarArrayAdd(vaResidues, (GENP) & rRes);
        }
    } else {

        /* create a single residue with the same name as the unit name:  */

        vaResidues = vaVarArrayCreate(sizeof(RESIDUE));
        rRes = (RESIDUE) oCreate(RESIDUEid);
        ContainerSetName(rRes, sUnitName);
        VarArrayAdd(vaResidues, (GENP) & rRes);
	}

    /* Connect everything together */

    /* Put the ATOMs within the RESIDUES */

    for (i = 0; i < iVarArrayElementCount(vaAtoms); i++) {
        aAtom = *PVAI(vaAtoms, ATOM, i);
        iResidue = iAtomTempInt(aAtom);
        rRes = *PVAI(vaResidues, RESIDUE, iResidue - 1);
        ContainerAdd((CONTAINER) rRes, (OBJEKT) aAtom);
    }

    /* Put the RESIDUES within the UNIT */

    for (i = 0; i < iVarArrayElementCount(vaResidues); i++) {
        rRes = *PVAI(vaResidues, RESIDUE, i);
        ContainerAdd((CONTAINER) uUnit, (OBJEKT) rRes);
    }
    VarArrayDestroy(&vaAtoms);
    VarArrayDestroy(&vaResidues);

    return (uUnit);

  FAIL:
    /* Destroy everything */
    if (vaAtoms) {
        for (i = 0; i < iVarArrayElementCount(vaAtoms); i++) {
            aAtom = *PVAI(vaAtoms, ATOM, i);
            Destroy((OBJEKT *) & aAtom);
        }
        VarArrayDestroy(&vaAtoms);
    }
    if (vaResidues) {
        for (i = 0; i < iVarArrayElementCount(vaResidues); i++) {
            rRes = *PVAI(vaResidues, RESIDUE, i);
            Destroy((OBJEKT *) & rRes);
        }
        VarArrayDestroy(&vaResidues);
    }
    return (NULL);
}
