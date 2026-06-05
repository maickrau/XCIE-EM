#ifndef FileHandler_h
#define FileHandler_h

#include <vector>
#include <string>
#include <unordered_map>
#include "EM.h"
#include "PseudobulkAnalysis.h"

std::vector<CellMatch> readScReadCountsMatchCounts(const std::string& scReadCountsFile);
std::vector<CellMatch> readMatchCounts(const std::string& matchTableFile);
void writeResultCellOnlyNonescapeVariants(const EMOutput& result, const bool phasesAreMatPat, std::ostream& stream);
void writeResultCells(const EMOutput& result, const bool phasesAreMatPat, std::ostream& stream);
void writeResultVariants(const EMOutput& result, const bool phasesAreMatPat, std::ostream& stream);
void writePseudobulkResults(const std::vector<PseudobulkInfo>& pseudobulk, const bool phasesAreMatPat, const std::string& firstColumnName, const std::string filename);
std::pair<std::unordered_map<std::string, bool>, bool> readForcedVariantPhases(const std::string& filename);
void writeCellMatchCounts(const std::vector<CellMatch>& cellMatches, const std::string& filename);
std::vector<std::tuple<size_t, size_t, std::string, std::string>> getGeneInfo(const std::string gff3Path, const bool onlyProteinCoding);
void writeGenesPerVariant(const EMOutput& output, const std::vector<std::vector<std::pair<std::string, std::string>>>& variantToGeneMatch, const std::string& filename);

#endif
