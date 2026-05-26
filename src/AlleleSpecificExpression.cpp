#include <fstream>
#include <regex>
#include <map>
#include <cassert>
#include <algorithm>
#include <unordered_set>
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

std::vector<PseudobulkInfo> getVariantPseudobulk(const EMOutput& output, const std::vector<CellMatch>& cellMatches, const double minConfidence)
{
	std::unordered_map<size_t, size_t> variantCoverageMatXa;
	std::unordered_map<size_t, size_t> variantCoverageMatXi;
	std::unordered_map<size_t, size_t> variantCoveragePatXa;
	std::unordered_map<size_t, size_t> variantCoveragePatXi;
	std::unordered_set<size_t> includedVariants;
	std::vector<bool> ignoreNothing;
	ignoreNothing.resize(output.helpers.numVariants(), false);
	for (size_t variantIndex = 0; variantIndex < output.helpers.numVariants(); variantIndex++)
	{
		double phaseScoreDifference = output.additions.variantPhaseConfidence[variantIndex];
		if (phaseScoreDifference < minConfidence) continue;
		includedVariants.insert(variantIndex);
	}
	std::unordered_set<size_t> includedCells;
	for (size_t cellIndex = 0; cellIndex < output.helpers.numCells(); cellIndex++)
	{
		double scoreDifference = output.additions.cellActiveConfidence[cellIndex];
		if (scoreDifference < minConfidence) continue;
		includedCells.insert(cellIndex);
	}
	for (const auto& t : cellMatches)
	{
		const size_t variantIndex = output.helpers.variantNameToIndex.at(t.variant);
		const size_t cellIndex = output.helpers.cellNameToIndex.at(t.cell);
		if (includedVariants.count(variantIndex) == 0) continue;
		if (includedCells.count(cellIndex) == 0) continue;
		bool activeCoverage = (output.result.variantIsMatRef[variantIndex] == output.result.cellIsMatActive[cellIndex]) == (!t.alt);
		if (activeCoverage && output.result.cellIsMatActive[cellIndex])
		{
			variantCoverageMatXa[variantIndex] += t.count;
		}
		if (!activeCoverage && output.result.cellIsMatActive[cellIndex])
		{
			variantCoverageMatXi[variantIndex] += t.count;
		}
		if (activeCoverage && !output.result.cellIsMatActive[cellIndex])
		{
			variantCoveragePatXa[variantIndex] += t.count;
		}
		if (!activeCoverage && !output.result.cellIsMatActive[cellIndex])
		{
			variantCoveragePatXi[variantIndex] += t.count;
		}
	}
	std::vector<std::string> variantOrder = getVariantOrder(output.helpers.variantNameToIndex);
	std::vector<PseudobulkInfo> result;
	result.reserve(output.helpers.numVariants());
	for (const std::string& name : variantOrder)
	{
		size_t index = output.helpers.variantNameToIndex.at(name);
		result.emplace_back();
		result.back().name = name;
		result.back().matXa = variantCoverageMatXa[index];
		result.back().matXi = variantCoverageMatXi[index];
		result.back().patXa = variantCoveragePatXa[index];
		result.back().patXi = variantCoveragePatXi[index];
	}
	return result;
}

std::string parseTag(const std::string& tags, const std::string& tagname)
{
	std::regex regex { "[;^]" + tagname + "=([^;]+)[;$]" };
	std::smatch matches;
	if (std::regex_match(tags, matches, regex))
	{
		return matches[1];
	}
	return "";
}

std::vector<std::tuple<size_t, size_t, std::string, std::string>> getGeneInfo(const std::string gff3Path, const bool onlyProteinCoding)
{
	std::vector<std::tuple<size_t, size_t, std::string, std::string>> result;
	std::ifstream file { gff3Path };
	while (file.good())
	{
		std::string line;
		std::getline(file, line);
		if (!file.good()) break;
		if (line.size() < 10) continue;
		if (line[0] == '#') continue;
		auto parts = split(line, '\t');
		std::string chromosome = parts[0];
		if (chromosome != "23" && lowercase(chromosome) != "x" && lowercase(chromosome) != "chrx") continue;
		std::string type = parts[2];
		if (type != "gene") continue;
		std::string geneType = parseTag(parts[8], "gene_type");
		if (onlyProteinCoding && geneType != "protein_coding") continue;
		std::string geneId = parseTag(parts[8], "gene_id");
		std::string geneName = parseTag(parts[8], "gene_name");
		size_t startPos = std::stoull(parts[3]);
		size_t endPos = std::stoull(parts[4]);
		if (geneName == "") geneName = geneId;
		result.emplace_back(startPos, endPos, geneId, geneName);
	}
	std::sort(result.begin(), result.end());
	return result;
}

std::vector<PseudobulkInfo> getGenePseudobulk(const std::vector<PseudobulkInfo>& variantPseudobulk, const std::vector<std::tuple<size_t, size_t, std::string, std::string>>& geneInfo)
{
	Logger::Log.log(Logger::LogLevel::DebugInfo) << geneInfo.size() << " genes included" << std::endl;
	if (geneInfo.size() == 0) return std::vector<PseudobulkInfo> {};
	std::map<std::vector<std::string>, std::tuple<size_t, size_t, size_t, size_t>> counts;
	for (const auto& t : variantPseudobulk)
	{
		size_t position = parseVariantPosition(t.name);
		std::vector<std::string> key;
		for (const auto& t2 : geneInfo)
		{
			if (std::get<0>(t2) > position) break;
			if (std::get<1>(t2) < position) continue;
			key.push_back(std::get<3>(t2));
		}
		if (key.size() == 0) continue;
		std::sort(key.begin(), key.end());
		std::tuple<size_t, size_t, size_t, size_t>& values = counts[key];
		std::get<0>(values) += t.matXa;
		std::get<1>(values) += t.matXi;
		std::get<2>(values) += t.patXa;
		std::get<3>(values) += t.patXi;
	}
	std::vector<PseudobulkInfo> result;
	for (const auto& pair : counts)
	{
		assert(pair.first.size() != 0);
		result.emplace_back();
		for (const std::string& gene : pair.first)
		{
			result.back().name += gene + ",";
		}
		assert(result.back().name.size() >= 2);
		result.back().name.pop_back();
		result.back().matXa = std::get<0>(pair.second);
		result.back().matXi = std::get<1>(pair.second);
		result.back().patXa = std::get<2>(pair.second);
		result.back().patXi = std::get<3>(pair.second);
	}
	std::sort(result.begin(), result.end(), [](const auto& left, const auto& right) { return left.name < right.name; });
	return result;
}
