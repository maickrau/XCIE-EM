#ifndef EM_h
#define EM_h

#include <string>
#include <vector>
#include <unordered_map>
#include "AlleleSpecificExpression.h"

struct EMResultAdditions
{
public:
	std::vector<double> variantPhaseConfidence;
	std::vector<double> variantEscapeCILow;
	std::vector<double> variantEscapeCIHigh;
	std::vector<double> cellActiveConfidence;
	std::vector<double> cellActiveConfidenceOnlyNonescape;
	std::vector<double> cellEscapeCILow;
	std::vector<double> cellEscapeCIHigh;
};

struct EMResult
{
public:
	std::vector<bool> cellIsMatActive;
	std::vector<double> cellEscapeFraction;
	std::vector<bool> variantIsMatRef;
	std::vector<double> variantEscapeFraction;
};

struct EMHelperVariables
{
public:
	std::unordered_map<std::string, size_t> variantNameToIndex;
	std::unordered_map<std::string, size_t> cellNameToIndex;
	std::vector<size_t> variantCoverage;
	std::vector<size_t> cellCoverage;
	std::vector<double> cellCoverageFraction;
	std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> activeCellsPerVariant;
	std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> activeVariantsPerCell;
	size_t numVariants() const;
	size_t numCells() const;
};

struct EMOutput
{
public:
	EMHelperVariables helpers;
	EMResult result;
	EMResultAdditions additions;
};

EMOutput runEM(const std::vector<CellMatch>& cellMatches, const std::unordered_map<std::string, bool>& forcedPhases, const size_t randomSeed, const double initialNoiseMagnitude, const double noiseDecay, const size_t numTries);

#endif
