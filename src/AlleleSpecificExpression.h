#ifndef AlleleSpecificExpression_h
#define AlleleSpecificExpression_h

#include <string>
#include <vector>

struct CellMatch
{
public:
	std::string cell;
	std::string variant;
	bool alt;
	size_t count;
};

struct PseudobulkInfo
{
public:
	std::string name;
	size_t matXa;
	size_t matXi;
	size_t patXa;
	size_t patXi;
};

struct EMOutput;

std::vector<CellMatch> filterOutHomozygousSites(const std::vector<CellMatch>& raw);
std::vector<CellMatch> excludeRegions(const std::vector<CellMatch>& raw, const std::vector<std::pair<size_t, size_t>>& excludedRegions);
std::vector<PseudobulkInfo> getVariantPseudobulk(const EMOutput& output, const std::vector<CellMatch>& cellMatches, const double minConfidence);
std::vector<std::tuple<size_t, size_t, std::string, std::string>> getGeneInfo(const std::string gff3Path, const bool onlyProteinCoding);
std::vector<PseudobulkInfo> getGenePseudobulk(const std::vector<PseudobulkInfo>& variantPseudobulk, const std::vector<std::tuple<size_t, size_t, std::string, std::string>>& geneInfo);

#endif
