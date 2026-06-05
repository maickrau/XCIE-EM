#include <map>
#include <cassert>
#include <tuple>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <iostream>
#include <zstr.hpp>
#include "Common.h"
#include "PseudobulkAnalysis.h"

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
		result.back().variantIndex = index;
		result.back().name = name;
		result.back().matXa = variantCoverageMatXa[index];
		result.back().matXi = variantCoverageMatXi[index];
		result.back().patXa = variantCoveragePatXa[index];
		result.back().patXi = variantCoveragePatXi[index];
	}
	return result;
}

std::vector<PseudobulkInfo> getGenePseudobulk(const std::vector<PseudobulkInfo>& variantPseudobulk, const std::vector<std::vector<std::pair<std::string, std::string>>>& variantGeneMatch)
{
	assert(variantGeneMatch.size() >= variantPseudobulk.size());
	std::map<std::vector<std::string>, std::tuple<size_t, size_t, size_t, size_t>> counts;
	for (const auto& t : variantPseudobulk)
	{
		std::vector<std::string> key;
		for (const std::pair<std::string, std::string>& gene : variantGeneMatch[t.variantIndex])
		{
			key.push_back(gene.second);
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
		if (std::get<0>(pair.second) == 0 && std::get<1>(pair.second) == 0 && std::get<2>(pair.second) == 0 && std::get<3>(pair.second) == 0) continue;
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

std::vector<std::vector<std::pair<std::string, std::string>>> getVariantGeneContainment(const EMOutput& output, const std::vector<std::tuple<size_t, size_t, std::string, std::string>>& geneInfo)
{
	std::vector<std::vector<std::pair<std::string, std::string>>> result;
	result.resize(output.helpers.numVariants());
	for (const auto& pair : output.helpers.variantNameToIndex)
	{
		const size_t position = parseVariantPosition(pair.first);
		const size_t i = pair.second;
		assert(i < result.size());
		for (const auto& t2 : geneInfo)
		{
			if (std::get<0>(t2) > position) break;
			if (std::get<1>(t2) < position) continue;
			result[i].emplace_back(std::get<2>(t2), std::get<3>(t2));
		}
		std::sort(result[i].begin(), result[i].end());
	}
	return result;
}
