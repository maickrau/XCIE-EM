#ifndef SNPCounter_h
#define SNPCounter_h

#include "AlleleSpecificExpression.h"

std::vector<CellMatch> countSNPsFromBamVcf(const std::string& vcfFile, const std::string& bamFile);

#endif
