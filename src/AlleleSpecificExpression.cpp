#include <fstream>
#include <regex>
#include <map>
#include <cassert>
#include <algorithm>
#include <unordered_set>
#include <zstr.hpp> //https://github.com/mateidavid/zstr
#include "Logger.h"
#include "EM.h"
#include "AlleleSpecificExpression.h"
#include "Common.h"

std::vector<CellMatch> filterOutHomozygousSites(const std::vector<CellMatch>& raw)
{
	const double minMinorAlleleRatio = 0.2;
	std::unordered_map<std::string, size_t> variantRefCount;
	std::unordered_map<std::string, size_t> variantAltCount;
	for (const auto& match : raw)
	{
		if (match.alt)
		{
			variantAltCount[match.variant] += match.count;
		}
		else
		{
			variantRefCount[match.variant] += match.count;
		}
	}
	std::vector<CellMatch> result;
	for (const auto& match : raw)
	{
		if ((double)variantRefCount[match.variant] < (double)variantAltCount[match.variant] * minMinorAlleleRatio) continue;
		if ((double)variantAltCount[match.variant] < (double)variantRefCount[match.variant] * minMinorAlleleRatio) continue;
		result.emplace_back(match);
	}
	return result;
}

std::vector<CellMatch> excludeRegions(const std::vector<CellMatch>& raw, const std::vector<std::pair<size_t, size_t>>& excludedRegions)
{
	std::vector<CellMatch> result;
	for (const auto& t : raw)
	{
		size_t variantPos = parseVariantPosition(t.variant);
		bool excluded = false;
		for (const auto& region : excludedRegions)
		{
			if (region.first > variantPos) break;
			if (region.second < variantPos) continue;
			excluded = true;
			break;
		}
		if (excluded) continue;
		result.emplace_back(t);
	}
	return result;
}

