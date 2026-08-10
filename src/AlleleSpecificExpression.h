#ifndef AlleleSpecificExpression_h
#define AlleleSpecificExpression_h

#include <string>
#include <vector>
#include <unordered_set>

struct CellMatch
{
public:
	std::string cell;
	std::string variant;
	bool alt;
	size_t count;
};


std::vector<CellMatch> filterToValidBarcodes(const std::vector<CellMatch>& raw, const std::unordered_set<std::string>& validBarcodes);
std::vector<CellMatch> filterOutHomozygousSites(const std::vector<CellMatch>& raw);
std::vector<CellMatch> excludeRegions(const std::vector<CellMatch>& raw, const std::vector<std::pair<size_t, size_t>>& excludedRegions);

#endif
