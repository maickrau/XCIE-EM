#ifndef OutputBuilder_h
#define OutputBuilder_h

#include <string>
#include <vector>
#include <unordered_map>
#include "../src/EM.h"
#include "InternalDeclarations.h"

namespace TestHelpers
{

// Builds an EMOutput by hand with the given phases and confidences, without running the EM.
// Variants and cells not mentioned in the maps get phase/activity false and confidence 0.
struct OutputBuilder
{
	std::vector<CellMatch> matches;
	std::unordered_map<std::string, bool> variantIsMatRef;
	std::unordered_map<std::string, double> variantConfidence;
	std::unordered_map<std::string, double> variantEscape;
	std::unordered_map<std::string, bool> cellIsMatActive;
	std::unordered_map<std::string, double> cellConfidence;
	std::unordered_map<std::string, double> cellConfidenceOnlyNonescape;
	std::unordered_map<std::string, double> cellEscape;
	EMOutput build() const
	{
		EMOutput output;
		output.helpers = getHelpers(matches, {});
		output.result = initializeResult(output.helpers);
		const size_t numVariants = output.helpers.numVariants();
		const size_t numCells = output.helpers.numCells();
		output.additions.variantPhaseConfidence.resize(numVariants, 0);
		output.additions.variantEscapeCILow.resize(numVariants, 0);
		output.additions.variantEscapeCIHigh.resize(numVariants, 0);
		output.additions.cellActiveConfidence.resize(numCells, 0);
		output.additions.cellActiveConfidenceOnlyNonescape.resize(numCells, 0);
		output.additions.cellEscapeCILow.resize(numCells, 0);
		output.additions.cellEscapeCIHigh.resize(numCells, 0);
		for (auto& escape : output.result.variantEscapeFraction) escape = 0.001;
		for (auto& escape : output.result.cellEscapeFraction) escape = 0.001;
		for (const auto& pair : output.helpers.variantNameToIndex)
		{
			const size_t i = pair.second;
			if (variantIsMatRef.count(pair.first)) output.result.variantIsMatRef[i] = variantIsMatRef.at(pair.first);
			if (variantConfidence.count(pair.first)) output.additions.variantPhaseConfidence[i] = variantConfidence.at(pair.first);
			if (variantEscape.count(pair.first)) output.result.variantEscapeFraction[i] = variantEscape.at(pair.first);
		}
		for (const auto& pair : output.helpers.cellNameToIndex)
		{
			const size_t i = pair.second;
			if (cellIsMatActive.count(pair.first)) output.result.cellIsMatActive[i] = cellIsMatActive.at(pair.first);
			if (cellConfidence.count(pair.first)) output.additions.cellActiveConfidence[i] = cellConfidence.at(pair.first);
			if (cellConfidenceOnlyNonescape.count(pair.first)) output.additions.cellActiveConfidenceOnlyNonescape[i] = cellConfidenceOnlyNonescape.at(pair.first);
			if (cellEscape.count(pair.first)) output.result.cellEscapeFraction[i] = cellEscape.at(pair.first);
		}
		return output;
	}
};

}

#endif
