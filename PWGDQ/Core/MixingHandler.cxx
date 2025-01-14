// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/VarManager.h"

#include <iostream>
#include <fstream>
using namespace std;

#include <TMath.h>
#include <TTimeStamp.h>
#include <TRandom.h>

ClassImp(MixingHandler);

//_________________________________________________________________________
MixingHandler::MixingHandler() : TNamed(),
                                 fIsInitialized(kFALSE),
                                 fVariableLimits(),
                                 fVariables()
{
  //
  // default constructor
  //
}

//_________________________________________________________________________
MixingHandler::MixingHandler(const char* name, const char* title) : TNamed(name, title),
                                                                    fIsInitialized(kFALSE),
                                                                    fVariableLimits(),
                                                                    fVariables()
{
  //
  // Named constructor
  //
}

//_________________________________________________________________________
MixingHandler::~MixingHandler()
{
  //
  // destructor
  //
}

//_________________________________________________________________________
void MixingHandler::AddMixingVariable(int var, int nBins, float* binLims)
{
  //
  // add a mixing variable
  //
  fVariables.push_back(var);
  TArrayF varBins;
  varBins.Set(nBins, binLims);
  fVariableLimits.push_back(varBins);
  VarManager::SetUseVariable(var);
}

//_________________________________________________________________________
void MixingHandler::AddMixingVariable(int var, int nBins, std::vector<float> binLims)
{

  float* bins = new float[nBins];
  for (int i = 0; i < nBins; ++i) {
    bins[i] = binLims[i];
  }
  AddMixingVariable(var, nBins, bins);
}

//_________________________________________________________________________
int MixingHandler::GetMixingVariable(VarManager::Variables var)
{
  int varNum = -1;
  for (int iVar = 0; iVar < fVariables.size(); ++iVar) {
    if (fVariables[iVar] == var) {
      varNum = iVar;
    }
  }
  return varNum;
}

//_________________________________________________________________________
std::vector<float> MixingHandler::GetMixingVariableLimits(VarManager::Variables var)
{
  std::vector<float> binLimits;
  for (int iVar = 0; iVar < fVariables.size(); ++iVar) {
    if (fVariables[iVar] == var) {
      for (int iBin = 0; iBin < fVariableLimits[iVar].GetSize(); ++iBin) {
        binLimits.push_back(fVariableLimits[iVar].At(iBin));
      }
    }
  }
  return binLimits;
}

//_________________________________________________________________________
void MixingHandler::Init()
{
  //
  // Initialization of pools
  //       The correct event category will be retrieved using the function FindEventCategory()
  //
  int size = 1;
  for (int iVar = 0; iVar < fVariables.size(); ++iVar) {
    size *= (fVariableLimits[iVar].GetSize() - 1);
  }

  fIsInitialized = kTRUE;
}

//_________________________________________________________________________
int MixingHandler::FindEventCategory(float* values)
{
  //
  // Find the event category corresponding to the added mixing variables
  //
  if (fVariables.size() == 0) {
    return -1;
  }
  if (!fIsInitialized) {
    Init();
  }

  std::vector<int> bin;
  for (int i = 0; i < fVariables.size(); ++i) {
    int binValue = TMath::BinarySearch(fVariableLimits[i].GetSize(), fVariableLimits[i].GetArray(), values[fVariables[i]]);
    bin.push_back(binValue);
    if (bin[i] == -1 || bin[i] == fVariableLimits[i].GetSize() - 1) {
      return -1; // all variables must be inside limits
    }
  }

  int category = 0;
  for (int iVar = 0; iVar < fVariables.size(); ++iVar) {
    int tempCategory = 1;
    for (int iVar2 = iVar; iVar2 < fVariables.size(); ++iVar2) {
      if (iVar2 == iVar) {
        tempCategory *= bin[iVar2];
      } else {
        tempCategory *= (fVariableLimits[iVar2].GetSize() - 1);
      }
    }
    category += tempCategory;
  }
  return category;
}

//_________________________________________________________________________
int MixingHandler::GetBinFromCategory(VarManager::Variables var, int category) const
{
  //
  // find the bin in variable var for the n-dimensional "category"
  //
  if (fVariables.size() == 0) {
    return -1;
  }

  // Search for the position of the variable "var" in the internal variable list of the handler
  int tempVar;
  for (int i = 0; i < fVariables.size(); ++i) {
    if (fVariables[i] == var) {
      tempVar = i;
    }
  }

  // extract the bin position in variable "var" from the category
  int norm = 1;
  for (int i = fVariables.size() - 1; i > tempVar; --i) {
    norm *= (fVariableLimits[i].GetSize() - 1);
  }
  int truncatedCategory = category - (category % norm);
  truncatedCategory /= norm;
  return truncatedCategory % (fVariableLimits[tempVar].GetSize() - 1);
}
