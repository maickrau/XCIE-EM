#ifndef Pseudobulkanalysis_h
#define Pseudobulkanalysis_h

#include <string>
#include <vector>
#include "AlleleSpecificExpression.h"
#include "EM.h"

struct PseudobulkInfo
{
public:
	size_t variantIndex;
	std::string name;
	size_t matXa;
	size_t matXi;
	size_t patXa;
	size_t patXi;
};

std::vector<PseudobulkInfo> getVariantPseudobulk(const EMOutput& output, const std::vector<CellMatch>& cellMatches, const double minConfidence);
std::vector<PseudobulkInfo> getVariantGroupPseudobulk(const EMOutput& output, const std::vector<CellMatch>& cellMatches, const double minConfidence, const std::unordered_map<std::string, std::string>& cellGrouping);
std::vector<std::vector<std::pair<std::string, std::string>>> getVariantGeneContainment(const EMOutput& output, const std::vector<std::tuple<size_t, size_t, std::string, std::string>>& geneInfo);
std::vector<PseudobulkInfo> getGenePseudobulk(const std::vector<PseudobulkInfo>& variantPseudobulk, const std::vector<std::vector<std::pair<std::string, std::string>>>& variantGeneMatch);
std::vector<PseudobulkInfo> getGeneGroupPseudobulk(const EMOutput& output, const std::vector<CellMatch>& cellMatches, const std::vector<std::vector<std::pair<std::string, std::string>>>& variantGeneMatch, const double minConfidence, const std::unordered_map<std::string, std::string>& cellGrouping);

#endif
