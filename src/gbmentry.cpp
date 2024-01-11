// GBM by Greg Ridgeway  Copyright (C) 2003

#include "gbm.h"

extern "C" {

#include <R.h>
#include <Rinternals.h>

SEXP gbm_fit
(
    SEXP radY,       // outcome or response
    SEXP radOffset,  // offset for f(x), NA for no offset
    SEXP radX,        
    SEXP raiXOrder,        
    SEXP radWeight,
    SEXP radMisc,   // other row specific data (eg failure time), NA=no Misc
    SEXP rcRows,
    SEXP rcCols,
    SEXP racVarClasses,
    SEXP ralMonotoneVar,
    SEXP rszFamily, 
    SEXP rcTrees,
    SEXP rcDepth,       // interaction depth
    SEXP rcMinObsInNode,
    SEXP rcNumClasses,
    SEXP rdShrinkage,
    SEXP rdBagFraction,
    SEXP rcTrain,
    SEXP radFOld,
    SEXP rcCatSplitsOld,
    SEXP rcTreesOld,
    SEXP rfVerbose
)
{
    unsigned long hr = 0;

    SEXP rAns = NULL;
    SEXP rNewTree = NULL;
    SEXP riSplitVar = NULL;
    SEXP rdSplitPoint = NULL;
    SEXP riLeftNode = NULL;
    SEXP riRightNode = NULL;
    SEXP riMissingNode = NULL;
    SEXP rdErrorReduction = NULL;
    SEXP rdWeight = NULL;
    SEXP rdPred = NULL;

    SEXP rdInitF = NULL;
    SEXP radF = NULL;
    SEXP radTrainError = NULL;
    SEXP radValidError = NULL;
    SEXP radOOBagImprove = NULL;

    SEXP rSetOfTrees = NULL;
    SEXP rSetSplitCodes = NULL;
    SEXP rSplitCode = NULL;

    VEC_VEC_CATEGORIES vecSplitCodes;

    int i = 0;
    int iT = 0;
    int iK = 0;
    int cTrees = INTEGER(rcTrees)[0];
    const int cResultComponents = 7;
    // rdInitF, radF, radTrainError, radValidError, radOOBagImprove
    // rSetOfTrees, rSetSplitCodes
    const int cTreeComponents = 8;
    // riSplitVar, rdSplitPoint, riLeftNode,
    // riRightNode, riMissingNode, rdErrorReduction, rdWeight, rdPred
    int cNodes = 0;
    int cTrain = INTEGER(rcTrain)[0];
    int cNumClasses = INTEGER(rcNumClasses)[0];

    double dTrainError = 0.0;
    double dValidError = 0.0;
    double dOOBagImprove = 0.0;

    CGBM *pGBM = NULL;
    CDataset *pData = NULL;
    CDistribution *pDist = NULL;
    int cGroups = -1;

    // set up the dataset
    pData = new CDataset();
    if(pData==NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }

    // initialize R's random number generator
    GetRNGstate();

    // initialize some things
    hr = gbm_setup(REAL(radY),
                   REAL(radOffset),
                   REAL(radX),
                   INTEGER(raiXOrder),
                   REAL(radWeight),
                   REAL(radMisc),
                   INTEGER(rcRows)[0],
                   INTEGER(rcCols)[0],
                   INTEGER(racVarClasses),
                   INTEGER(ralMonotoneVar),
                   CHAR(STRING_ELT(rszFamily,0)),
                   INTEGER(rcTrees)[0],
                   INTEGER(rcDepth)[0],
                   INTEGER(rcMinObsInNode)[0],
                   INTEGER(rcNumClasses)[0],
                   REAL(rdShrinkage)[0],
                   REAL(rdBagFraction)[0],
                   INTEGER(rcTrain)[0],
                   pData,
                   pDist,
                   cGroups);
    if(GBM_FAILED(hr))
    {
        goto Error;
    }
       
    // allocate the GBM
    pGBM = new CGBM();
    if(pGBM==NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }

    // initialize the GBM
    hr = pGBM->Initialize(pData,
                          pDist,
                          REAL(rdShrinkage)[0], 
                          cTrain, 
                          REAL(rdBagFraction)[0],
                          INTEGER(rcDepth)[0],
                          INTEGER(rcMinObsInNode)[0],
                          INTEGER(rcNumClasses)[0],
                          cGroups);

    if(GBM_FAILED(hr))
    {
        goto Error;
    }

    // allocate the main return object
    PROTECT(rAns = allocVector(VECSXP, cResultComponents));

    // allocate the initial value
    PROTECT(rdInitF = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(rAns,0,rdInitF);
    UNPROTECT(1); // rdInitF

    // allocate the predictions
    PROTECT(radF = allocVector(REALSXP, (pData->cRows) * cNumClasses));
    SET_VECTOR_ELT(rAns,1,radF);
    UNPROTECT(1); // radF

    hr = pDist->Initialize(pData->adY,
                           pData->adMisc,
                           pData->adOffset,
                           pData->adWeight,
                           pData->cRows);
    if(GBM_FAILED(hr))
    {
      goto Error;
    }
    
    if(ISNA(REAL(radFOld)[0])) // check for old predictions
    {
        // set the initial value of F as a constant
        hr = pDist->InitF(pData->adY,
                          pData->adMisc,
                          pData->adOffset,
                          pData->adWeight,
                          REAL(rdInitF)[0], 
                          cTrain);
        if(GBM_FAILED(hr))
        {
          goto Error;
        }
      
        for(i=0; i < (pData->cRows) * cNumClasses; i++)
        {
            REAL(radF)[i] = REAL(rdInitF)[0];
        }
    }
    else
    {
        for(i=0; i < (pData->cRows) * cNumClasses; i++)
        {
            REAL(radF)[i] = REAL(radFOld)[i];
        }
    }

    // allocate space for the performance measures
    PROTECT(radTrainError = allocVector(REALSXP, cTrees));
    PROTECT(radValidError = allocVector(REALSXP, cTrees));
    PROTECT(radOOBagImprove = allocVector(REALSXP, cTrees));
    SET_VECTOR_ELT(rAns,2,radTrainError);
    SET_VECTOR_ELT(rAns,3,radValidError);
    SET_VECTOR_ELT(rAns,4,radOOBagImprove);
    UNPROTECT(3); // radTrainError , radValidError, radOOBagImprove

    // allocate the component for the tree structures
    PROTECT(rSetOfTrees = allocVector(VECSXP, cTrees * cNumClasses));
    SET_VECTOR_ELT(rAns,5,rSetOfTrees);
    UNPROTECT(1); // rSetOfTrees

    if(INTEGER(rfVerbose)[0])
    {
       Rprintf("Iter   TrainDeviance   ValidDeviance   StepSize   Improve\n");
    }
    for(iT=0; iT<cTrees; iT++)
    {
        // Update the parameters
        hr = pDist->UpdateParams(REAL(radF), 
                                 pData->adOffset, 
                                 pData->adWeight, 
                                 cTrain);
        if(GBM_FAILED(hr))
        {
           goto Error;
        }
        
        REAL(radTrainError)[iT] = 0.0;
        REAL(radValidError)[iT] = 0.0;
        REAL(radOOBagImprove)[iT] = 0.0;
        for (iK = 0; iK < cNumClasses; iK++)
        {
            hr = pGBM->iterate(REAL(radF),
                               dTrainError,dValidError,dOOBagImprove,
                               cNodes, cNumClasses, iK);

            if(GBM_FAILED(hr))
            {
                goto Error;
            }

            // store the performance measures
            REAL(radTrainError)[iT] += dTrainError;
            REAL(radValidError)[iT] += dValidError;
            REAL(radOOBagImprove)[iT] += dOOBagImprove;

            // allocate the new tree component for the R list structure
            PROTECT(rNewTree = allocVector(VECSXP, cTreeComponents));
            // riNodeID,riSplitVar,rdSplitPoint,riLeftNode,
            // riRightNode,riMissingNode,rdErrorReduction,rdWeight
            PROTECT(riSplitVar = allocVector(INTSXP, cNodes));
            PROTECT(rdSplitPoint = allocVector(REALSXP, cNodes));
            PROTECT(riLeftNode = allocVector(INTSXP, cNodes));
            PROTECT(riRightNode = allocVector(INTSXP, cNodes));
            PROTECT(riMissingNode = allocVector(INTSXP, cNodes));
            PROTECT(rdErrorReduction = allocVector(REALSXP, cNodes));
            PROTECT(rdWeight = allocVector(REALSXP, cNodes));
            PROTECT(rdPred = allocVector(REALSXP, cNodes));
            SET_VECTOR_ELT(rNewTree,0,riSplitVar);
            SET_VECTOR_ELT(rNewTree,1,rdSplitPoint);
            SET_VECTOR_ELT(rNewTree,2,riLeftNode);
            SET_VECTOR_ELT(rNewTree,3,riRightNode);
            SET_VECTOR_ELT(rNewTree,4,riMissingNode);
            SET_VECTOR_ELT(rNewTree,5,rdErrorReduction);
            SET_VECTOR_ELT(rNewTree,6,rdWeight);
            SET_VECTOR_ELT(rNewTree,7,rdPred);
            UNPROTECT(cTreeComponents); 
            SET_VECTOR_ELT(rSetOfTrees,(iK + iT * cNumClasses),rNewTree);
            UNPROTECT(1); // rNewTree
        
            hr = gbm_transfer_to_R(pGBM,
                                   vecSplitCodes,
                                   INTEGER(riSplitVar),
                                   REAL(rdSplitPoint),
                                   INTEGER(riLeftNode),
                                   INTEGER(riRightNode),
                                   INTEGER(riMissingNode),
                                   REAL(rdErrorReduction),
                                   REAL(rdWeight),
                                   REAL(rdPred),
                                   INTEGER(rcCatSplitsOld)[0]);
        } // Close for iK

        // print the information
        if((iT <= 9) ||
           ((iT+1+INTEGER(rcTreesOld)[0])/20 ==
            (iT+1+INTEGER(rcTreesOld)[0])/20.0) ||
            (iT==cTrees-1))
        {
            R_CheckUserInterrupt();
            if(INTEGER(rfVerbose)[0])
            {
               Rprintf("%6d %13.4f %15.4f %10.4f %9.4f\n",
                       iT+1+INTEGER(rcTreesOld)[0],
                       REAL(radTrainError)[iT],
                       REAL(radValidError)[iT],
                       REAL(rdShrinkage)[0],
                       REAL(radOOBagImprove)[iT]);
            }
        }
    }

    if(INTEGER(rfVerbose)[0]) Rprintf("\n");

    // transfer categorical splits to R
    PROTECT(rSetSplitCodes = allocVector(VECSXP, vecSplitCodes.size()));
    SET_VECTOR_ELT(rAns,6,rSetSplitCodes);
    UNPROTECT(1); // rSetSplitCodes

    for(i=0; i<(int)vecSplitCodes.size(); i++)
    {
        PROTECT(rSplitCode = 
                    allocVector(INTSXP, size_of_vector(vecSplitCodes,i)));
        SET_VECTOR_ELT(rSetSplitCodes,i,rSplitCode);
        UNPROTECT(1); // rSplitCode

        hr = gbm_transfer_catsplits_to_R(i,
                                         vecSplitCodes,
                                         INTEGER(rSplitCode));
    }

    // dump random number generator seed
    #ifdef NOISY_DEBUG
    Rprintf("PutRNGstate\n");
    #endif
    PutRNGstate();

Cleanup:
    UNPROTECT(1); // rAns
    #ifdef NOISY_DEBUG
    Rprintf("destructing\n");
    #endif

    if(pGBM != NULL)
    {
        delete pGBM;
        pGBM = NULL;
    }
    if(pDist != NULL)
    {
        delete pDist;
        pDist = NULL;
    }
    if(pData != NULL)
    {
        delete pData;
        pData = NULL;
    }

    return rAns;
Error:
    goto Cleanup;
}

SEXP gbm_pred
(
   SEXP radX,         // the data matrix
   SEXP rcRows,       // number of rows
   SEXP rcCols,       // number of columns
   SEXP rcNumClasses, // number of classes
   SEXP rcTrees,      // number of trees, may be a vector
   SEXP rdInitF,      // the initial value
   SEXP rTrees,       // the list of trees
   SEXP rCSplits,     // the list of categorical splits
   SEXP raiVarType,   // indicator of continuous/nominal
   SEXP riSingleTree  // boolean whether to return only results for one tree
)
{
   int iTree = 0;
   int iObs = 0;
   int cRows = INTEGER(rcRows)[0];
   int cPredIterations = LENGTH(rcTrees);
   int iPredIteration = 0;
   int cTrees = 0;
   int iClass = 0;
   int cNumClasses = INTEGER(rcNumClasses)[0];

   SEXP rThisTree = NULL;
   int *aiSplitVar = NULL;
   double *adSplitCode = NULL;
   int *aiLeftNode = NULL;
   int *aiRightNode = NULL;
   int *aiMissingNode = NULL;
   int iCurrentNode = 0;
   double dX = 0.0;
   int iCatSplitIndicator = 0;
   bool fSingleTree = (INTEGER(riSingleTree)[0]==1);

   SEXP radPredF = NULL;

   // allocate the predictions to return
   PROTECT(radPredF = allocVector(REALSXP, cRows*cNumClasses*cPredIterations));
   if(radPredF == NULL)
   {
      goto Error;
   }

   // initialize the predicted values
   if(!fSingleTree)
   {
      // initialize with the intercept for only the smallest rcTrees
      for(iObs=0; iObs<cRows*cNumClasses; iObs++)
      {
         REAL(radPredF)[iObs] = REAL(rdInitF)[0];
      }
   }
   else
   {
      for(iObs=0; iObs<cRows*cNumClasses*cPredIterations; iObs++)
      {
         REAL(radPredF)[iObs] = 0.0;
      }
   }

   iTree = 0;
   for(iPredIteration=0; iPredIteration<LENGTH(rcTrees); iPredIteration++)
   {
      cTrees = INTEGER(rcTrees)[iPredIteration];
      if(fSingleTree) iTree=cTrees-1;
      if(!fSingleTree && (iPredIteration>0))
      {
          // copy over from the last rcTrees
          for(iObs=0; iObs<cRows*cNumClasses; iObs++)
          {
             REAL(radPredF)[cRows*cNumClasses*iPredIteration+iObs] =
                REAL(radPredF)[cRows*cNumClasses*(iPredIteration-1)+iObs];
          }
      }
      while(iTree<cTrees*cNumClasses)
      {
         for (iClass = 0; iClass < cNumClasses; iClass++)
         {
            rThisTree   = VECTOR_ELT(rTrees,iTree);
            // these relate to columns returned by pretty.gbm.tree()
            aiSplitVar    = INTEGER(VECTOR_ELT(rThisTree,0));
            adSplitCode   = REAL   (VECTOR_ELT(rThisTree,1));
            aiLeftNode    = INTEGER(VECTOR_ELT(rThisTree,2));
            aiRightNode   = INTEGER(VECTOR_ELT(rThisTree,3));
            aiMissingNode = INTEGER(VECTOR_ELT(rThisTree,4));

            for(iObs=0; iObs<cRows; iObs++)
            {
               iCurrentNode = 0;
               while(aiSplitVar[iCurrentNode] != -1)
               {
                  dX = REAL(radX)[aiSplitVar[iCurrentNode]*cRows + iObs];
                  // missing?
                  if(ISNA(dX))
                  {
                     iCurrentNode = aiMissingNode[iCurrentNode];
                  }
                  // continuous?
                  else if(INTEGER(raiVarType)[aiSplitVar[iCurrentNode]] == 0)
                  {
                     if(dX < adSplitCode[iCurrentNode])
                     {
                        iCurrentNode = aiLeftNode[iCurrentNode];
                     }
                     else
                     {
                        iCurrentNode = aiRightNode[iCurrentNode];
                     }
                  }
                  else // categorical
                  {
                     iCatSplitIndicator = INTEGER(
                                 VECTOR_ELT(rCSplits,
                                            (int)adSplitCode[iCurrentNode]))[(int)dX];
                     if(iCatSplitIndicator==-1)
                     {
                        iCurrentNode = aiLeftNode[iCurrentNode];
                     }
                     else if(iCatSplitIndicator==1)
                     {
                        iCurrentNode = aiRightNode[iCurrentNode]; 
                     }  
                     else // categorical level not present in training
                     {
                        iCurrentNode = aiMissingNode[iCurrentNode];
                     }
                  }
               }
               REAL(radPredF)[cRows*cNumClasses*iPredIteration+cRows*iClass+iObs] += 
                  adSplitCode[iCurrentNode]; // add the prediction
            } // iObs
            iTree++;
         } // iClass
      } // iTree
   }  // iPredIteration
    
Cleanup:
    UNPROTECT(1); // radPredF
    return radPredF;
Error:
    goto Cleanup;
}


SEXP gbm_plot
(
    SEXP radX,          // vector or matrix of points to make predictions
    SEXP rcRows,        // number of rows in X
    SEXP rcCols,        // number of columns in X
    SEXP rcNumClasses,  // number of classes
    SEXP raiWhichVar,   // length=cCols, index of which var cols of X are
    SEXP rcTrees,       // number of trees to use
    SEXP rdInitF,       // initial value
    SEXP rTrees,        // tree list object
    SEXP rCSplits,      // categorical split list object
    SEXP raiVarType     // vector of variable types
)
{
    int i = 0;
    int iTree = 0;
    int iObs = 0;
    int iClass = 0;
    int cRows = INTEGER(rcRows)[0];
    int cCols = INTEGER(rcCols)[0];
    int cTrees = INTEGER(rcTrees)[0];
    int cNumClasses = INTEGER(rcNumClasses)[0];

    SEXP rThisTree = NULL;
    int *aiSplitVar = NULL;
    double *adSplitCode = NULL;
    int *aiLeftNode = NULL;
    int *aiRightNode = NULL;
    int *aiMissingNode = NULL;
    double *adW = NULL;
    int iCurrentNode = 0;
    double dCurrentW = 0.0;
    double dX = 0.0;
    int iCatSplitIndicator = 0;

    SEXP radPredF = NULL;
    int aiNodeStack[40];
    double adWeightStack[40];
    int cStackNodes = 0;
    int iPredVar = 0;

    // allocate the predictions to return
    PROTECT(radPredF = allocVector(REALSXP, cRows*cNumClasses));
    if(radPredF == NULL)
    {
        goto Error;
    }
    for(iObs=0; iObs<cRows*cNumClasses; iObs++)
    {
        REAL(radPredF)[iObs] = REAL(rdInitF)[0];
    }
    for(iTree=0; iTree<cTrees; iTree++)
    {
        for (iClass = 0; iClass < cNumClasses; iClass++)
        {
            rThisTree     = VECTOR_ELT(rTrees,iClass + iTree*cNumClasses);
            aiSplitVar    = INTEGER(VECTOR_ELT(rThisTree,0));
            adSplitCode   = REAL   (VECTOR_ELT(rThisTree,1));
            aiLeftNode    = INTEGER(VECTOR_ELT(rThisTree,2));
            aiRightNode   = INTEGER(VECTOR_ELT(rThisTree,3));
            aiMissingNode = INTEGER(VECTOR_ELT(rThisTree,4));
            adW           = REAL   (VECTOR_ELT(rThisTree,6));
            for(iObs=0; iObs<cRows; iObs++)
            {
                aiNodeStack[0] = 0;
                adWeightStack[0] = 1.0;
                cStackNodes = 1;
                while(cStackNodes > 0)
                {
                    cStackNodes--;
                    iCurrentNode = aiNodeStack[cStackNodes];

                    if(aiSplitVar[iCurrentNode] == -1) // terminal node
                    {
                        REAL(radPredF)[iClass*cRows + iObs] += 
                            adWeightStack[cStackNodes]*adSplitCode[iCurrentNode];
                    }
                    else // non-terminal node
                    {
                        // is this a split variable that interests me?
                        iPredVar = -1;
                        for(i=0; (iPredVar == -1) && (i < cCols); i++)
                        {
                            if(INTEGER(raiWhichVar)[i] == aiSplitVar[iCurrentNode])
                            {
                                iPredVar = i; // split is on one that interests me
                            }
                        }

                        if(iPredVar != -1) // this split is among raiWhichVar
                        {    
                            dX = REAL(radX)[iPredVar*cRows + iObs];
                            // missing?
                            if(ISNA(dX))
                            {
                                aiNodeStack[cStackNodes] = aiMissingNode[iCurrentNode];
                                cStackNodes++;                            
                            }
                            // continuous?
                            else if(INTEGER(raiVarType)[aiSplitVar[iCurrentNode]] == 0)
                                {
                                if(dX < adSplitCode[iCurrentNode])
                                {
                                    aiNodeStack[cStackNodes] = aiLeftNode[iCurrentNode];
                                        cStackNodes++;
                                }
                                else
                                {
                                    aiNodeStack[cStackNodes] = aiRightNode[iCurrentNode];
                                    cStackNodes++;
                                }
                            }
                            else // categorical
                            {
                                iCatSplitIndicator = INTEGER(
                                    VECTOR_ELT(rCSplits,
                                               (int)adSplitCode[iCurrentNode]))[(int)dX];
                                if(iCatSplitIndicator==-1)
                                {
                                    aiNodeStack[cStackNodes] = aiLeftNode[iCurrentNode];
                                    cStackNodes++;
                                }
                                else if(iCatSplitIndicator==1)
                                {
                                    aiNodeStack[cStackNodes] = aiRightNode[iCurrentNode];
                                    cStackNodes++;
                                }
                                else // handle unused level
                                {
                                    iCurrentNode = aiMissingNode[iCurrentNode];
                                }
                            }
                        } // iPredVar != -1
                        else // not interested in this split, average left and right 
                        {
                            aiNodeStack[cStackNodes] = aiRightNode[iCurrentNode];
                            dCurrentW = adWeightStack[cStackNodes];
                            adWeightStack[cStackNodes] = dCurrentW *
                                adW[aiRightNode[iCurrentNode]]/
                                (adW[aiLeftNode[iCurrentNode]]+
                                 adW[aiRightNode[iCurrentNode]]);
                            cStackNodes++;
                            aiNodeStack[cStackNodes] = aiLeftNode[iCurrentNode];
                            adWeightStack[cStackNodes] = 
                                    dCurrentW-adWeightStack[cStackNodes-1];
                            cStackNodes++;
                        }
                    } // non-terminal node
                } // while(cStackNodes > 0)
            } // iObs
        } // iClass
    } // iTree

Cleanup:
    UNPROTECT(1); // radPredF
    return radPredF;
Error:
    goto Cleanup;
} // gbm_plot

} // end extern "C"

