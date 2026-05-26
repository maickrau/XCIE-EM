#ifndef FileHandler_h
#define FileHandler_h

#include <vector>
#include <string>
#include <unordered_map>
#include "EM.h"

std::vector<CellMatch> readScReadCountsMatchCounts(const std::string& scReadCountsFile);
std::vector<CellMatch> readMatchCounts(const std::string& matchTableFile);
void writeResultCellOnlyNonescapeVariants(const EMOutput& result, const bool phasesAreMatPat, std::ostream& stream);
void writeResultCells(const EMOutput& result, const bool phasesAreMatPat, std::ostream& stream);
void writeResultVariants(const EMOutput& result, const bool phasesAreMatPat, std::ostream& stream);
void writePseudobulkVariantResults(const std::vector<PseudobulkVariantInfo>& variantPseudobulk, const bool phasesAreMatPat, const std::string filename);
std::unordered_map<std::string, bool> readForcedVariantPhases(const std::string& filename);
void writeCellMatchCounts(const std::vector<CellMatch>& cellMatches, const std::string& filename);

#endif
