#ifndef SNPCounter_h
#define SNPCounter_h

#include "AlleleSpecificExpression.h"

struct SNPMatch
{
public:
	std::string chromosome;
	size_t position;
	char ref;
	char alt;
	std::string barcode;
	size_t altFwCount;
	size_t altBwCount;
	size_t refFwCount;
	size_t refBwCount;
};

std::vector<SNPMatch> countSNPsFromBamVcf(const std::string& vcfFile, const std::string& bamFile);
std::vector<CellMatch> convertSNPMatchesToCellMatches(const std::vector<SNPMatch>& snpMatches);

#endif
