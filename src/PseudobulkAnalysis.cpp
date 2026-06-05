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

std::vector<bool> getIncludedVariants(const EMOutput& output, const double minConfidence)
{
	std::vector<bool> result;
	result.resize(output.helpers.numVariants(), false);
	for (size_t variantIndex = 0; variantIndex < output.helpers.numVariants(); variantIndex++)
	{
		double phaseScoreDifference = output.additions.variantPhaseConfidence[variantIndex];
		if (phaseScoreDifference < minConfidence) continue;
		result[variantIndex] = true;
	}
	return result;
}

std::vector<bool> getIncludedCells(const EMOutput& output, const double minConfidence)
{
	std::vector<bool> result;
	result.resize(output.helpers.numCells(), false);
	for (size_t cellIndex = 0; cellIndex < output.helpers.numCells(); cellIndex++)
	{
		double scoreDifference = output.additions.cellActiveConfidence[cellIndex];
		if (scoreDifference < minConfidence) continue;
		result[cellIndex] = true;
	}
	return result;
}

std::vector<std::unordered_map<std::string, std::tuple<size_t, size_t, size_t, size_t>>> getCellGroupedVariantPseudobulk(const EMOutput& output, const std::vector<CellMatch>& cellMatches, const double minConfidence, const std::unordered_map<std::string, std::string>& cellGrouping)
{
	std::vector<std::unordered_map<std::string, std::tuple<size_t, size_t, size_t, size_t>>> coveragesPerGroupPerVariant;
	coveragesPerGroupPerVariant.resize(output.helpers.numVariants());
	const std::vector<bool> includedVariants = getIncludedVariants(output, minConfidence);
	const std::vector<bool> includedCells = getIncludedCells(output, minConfidence);
	std::vector<std::string> cellGroup;
	cellGroup.resize(output.helpers.numCells(), "");
	for (const auto& pair : cellGrouping)
	{
		if (output.helpers.cellNameToIndex.count(pair.first) == 0) continue;
		assert(cellGroup[output.helpers.cellNameToIndex.at(pair.first)] == "");
		cellGroup[output.helpers.cellNameToIndex.at(pair.first)] = pair.second;
	}
	for (const auto& t : cellMatches)
	{
		const size_t variantIndex = output.helpers.variantNameToIndex.at(t.variant);
		const size_t cellIndex = output.helpers.cellNameToIndex.at(t.cell);
		if (!includedVariants[variantIndex]) continue;
		if (!includedCells[cellIndex]) continue;
		if (cellGroup[cellIndex] == "") continue;
		std::tuple<size_t, size_t, size_t, size_t>& values = coveragesPerGroupPerVariant[variantIndex][cellGroup[cellIndex]];
		bool activeCoverage = (output.result.variantIsMatRef[variantIndex] == output.result.cellIsMatActive[cellIndex]) == (!t.alt);
		if (activeCoverage && output.result.cellIsMatActive[cellIndex])
		{
			std::get<0>(values) += t.count;
		}
		if (!activeCoverage && output.result.cellIsMatActive[cellIndex])
		{
			std::get<1>(values) += t.count;
		}
		if (activeCoverage && !output.result.cellIsMatActive[cellIndex])
		{
			std::get<2>(values) += t.count;
		}
		if (!activeCoverage && !output.result.cellIsMatActive[cellIndex])
		{
			std::get<3>(values) += t.count;
		}
	}
	return coveragesPerGroupPerVariant;
}

std::vector<PseudobulkInfo> getVariantGroupPseudobulk(const EMOutput& output, const std::vector<CellMatch>& cellMatches, const double minConfidence, const std::unordered_map<std::string, std::string>& cellGrouping)
{
	std::vector<std::unordered_map<std::string, std::tuple<size_t, size_t, size_t, size_t>>> coveragesPerGroupPerVariant = getCellGroupedVariantPseudobulk(output, cellMatches, minConfidence, cellGrouping);
	std::vector<std::string> variantOrder = getVariantOrder(output.helpers.variantNameToIndex);
	std::vector<PseudobulkInfo> result;
	for (const std::string& name : variantOrder)
	{
		const size_t index = output.helpers.variantNameToIndex.at(name);
		std::vector<PseudobulkInfo> partsHere;
		for (const auto& pair : coveragesPerGroupPerVariant[index])
		{
			assert(pair.first != "");
			partsHere.emplace_back();
			partsHere.back().variantIndex = index;
			partsHere.back().name = name + "\t" + pair.first;
			partsHere.back().matXa = std::get<0>(pair.second);
			partsHere.back().matXi = std::get<1>(pair.second);
			partsHere.back().patXa = std::get<2>(pair.second);
			partsHere.back().patXi = std::get<3>(pair.second);
		}
		std::sort(partsHere.begin(), partsHere.end(), [](const auto& left, const auto& right) { return left.name < right.name; });
		result.insert(result.end(), partsHere.begin(), partsHere.end());
	}
	return result;
}

std::vector<PseudobulkInfo> getVariantPseudobulk(const EMOutput& output, const std::vector<CellMatch>& cellMatches, const double minConfidence)
{
	std::unordered_map<size_t, size_t> variantCoverageMatXa;
	std::unordered_map<size_t, size_t> variantCoverageMatXi;
	std::unordered_map<size_t, size_t> variantCoveragePatXa;
	std::unordered_map<size_t, size_t> variantCoveragePatXi;
	const std::vector<bool> includedVariants = getIncludedVariants(output, minConfidence);
	const std::vector<bool> includedCells = getIncludedCells(output, minConfidence);
	for (const auto& t : cellMatches)
	{
		const size_t variantIndex = output.helpers.variantNameToIndex.at(t.variant);
		const size_t cellIndex = output.helpers.cellNameToIndex.at(t.cell);
		if (!includedVariants[variantIndex]) continue;
		if (!includedCells[cellIndex]) continue;
		const bool activeCoverage = (output.result.variantIsMatRef[variantIndex] == output.result.cellIsMatActive[cellIndex]) == (!t.alt);
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
		const size_t index = output.helpers.variantNameToIndex.at(name);
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
	std::map<std::vector<std::pair<std::string, std::string>>, std::tuple<size_t, size_t, size_t, size_t>> counts;
	for (const auto& t : variantPseudobulk)
	{
		std::vector<std::pair<std::string, std::string>> key;
		for (const std::pair<std::string, std::string>& gene : variantGeneMatch[t.variantIndex])
		{
			key.emplace_back(gene.first, gene.second);
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
		std::vector<std::string> geneIds;
		std::vector<std::string> geneNames;
		for (const auto& pair2 : pair.first)
		{
			geneIds.push_back(pair2.first);
			geneNames.push_back(pair2.second);
		}
		result.back().name = join(',', geneIds) + "\t" + join(',', geneNames);
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

std::vector<PseudobulkInfo> getGeneGroupPseudobulk(const EMOutput& output, const std::vector<CellMatch>& cellMatches, const std::vector<std::vector<std::pair<std::string, std::string>>>& variantGeneMatch, const double minConfidence, const std::unordered_map<std::string, std::string>& cellGrouping)
{
	std::vector<std::unordered_map<std::string, std::tuple<size_t, size_t, size_t, size_t>>> coveragesPerGroupPerVariant = getCellGroupedVariantPseudobulk(output, cellMatches, minConfidence, cellGrouping);
	std::map<std::vector<std::pair<std::string, std::string>>, std::unordered_map<std::string, std::tuple<size_t, size_t, size_t, size_t>>> counts;
	for (size_t i = 0; i < output.helpers.numVariants(); i++)
	{
		std::vector<std::pair<std::string, std::string>> key;
		for (const std::pair<std::string, std::string>& gene : variantGeneMatch[i])
		{
			key.emplace_back(gene.first, gene.second);
		}
		if (key.size() == 0) continue;
		std::sort(key.begin(), key.end());
		for (const auto& t : coveragesPerGroupPerVariant[i])
		{
			std::tuple<size_t, size_t, size_t, size_t>& values = counts[key][t.first];
			std::get<0>(values) += std::get<0>(t.second);
			std::get<1>(values) += std::get<1>(t.second);
			std::get<2>(values) += std::get<2>(t.second);
			std::get<3>(values) += std::get<3>(t.second);
		}
	}
	std::vector<PseudobulkInfo> result;
	for (const auto& pair : counts)
	{
		assert(pair.first.size() != 0);
		std::vector<std::string> geneIds;
		std::vector<std::string> geneNames;
		for (const auto& pair2 : pair.first)
		{
			geneIds.push_back(pair2.first);
			geneNames.push_back(pair2.second);
		}
		for (const auto& pair2 : pair.second)
		{
			if (std::get<0>(pair2.second) == 0 && std::get<1>(pair2.second) == 0 && std::get<2>(pair2.second) == 0 && std::get<3>(pair2.second) == 0) continue;
			result.emplace_back();
			result.back().name = join(',', geneIds) + "\t" + join(',', geneNames) + "\t" + pair2.first;
			result.back().matXa = std::get<0>(pair2.second);
			result.back().matXi = std::get<1>(pair2.second);
			result.back().patXa = std::get<2>(pair2.second);
			result.back().patXi = std::get<3>(pair2.second);
		}
	}
	std::sort(result.begin(), result.end(), [](const auto& left, const auto& right) { return left.name < right.name; });
	return result;

}
