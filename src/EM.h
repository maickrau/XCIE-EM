#ifndef EM_h
#define EM_h

#include <string>
#include <vector>
#include <unordered_map>

struct CellMatch
{
public:
	std::string cell;
	std::string variant;
	bool alt;
	size_t count;
};

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
	std::vector<std::unordered_map<size_t, size_t>> cellVariantRefCount;
	std::vector<std::unordered_map<size_t, size_t>> cellVariantAltCount;
	std::vector<std::vector<size_t>> activeCellsPerVariant;
	std::vector<std::vector<size_t>> activeVariantsPerCell;
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

struct PseudobulkVariantInfo
{
public:
	std::string variantName;
	size_t matXa;
	size_t matXi;
	size_t patXa;
	size_t patXi;
};

std::vector<CellMatch> filterOutHomozygousSites(const std::vector<CellMatch>& raw);
std::vector<CellMatch> excludeRegions(const std::vector<CellMatch>& raw, const std::vector<std::pair<size_t, size_t>>& excludedRegions);
EMOutput runEM(const std::vector<CellMatch>& cellMatches, const std::unordered_map<std::string, bool>& forcedPhases, const size_t randomSeed, const double initialNoiseMagnitude, const double noiseDecay, const size_t numTries);
std::vector<PseudobulkVariantInfo> getVariantPseudobulk(const EMOutput& output, const std::vector<CellMatch>& cellMatches, const double minConfidence);

#endif
