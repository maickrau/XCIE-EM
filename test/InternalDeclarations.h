#ifndef InternalDeclarations_h
#define InternalDeclarations_h

// Declarations of functions and types which are defined in the src/*.cpp files
// but not exposed in the public headers. They have external linkage so the unit
// tests can call them directly. If a signature in src/ changes, update it here too.

#include <string>
#include <vector>
#include <tuple>
#include <random>
#include <istream>
#include <unordered_map>
#include "../src/EM.h"
#include "../src/AlleleSpecificExpression.h"
#include "../src/PseudobulkAnalysis.h"

// Common.cpp
double logNChooseK(const size_t n, const size_t k);

// FileHandler.cpp
std::unordered_map<std::string, std::string> parseTagsGff3(const std::string& tagstr);
std::vector<std::tuple<size_t, size_t, std::string, std::string>> getGeneInfoFromStream(std::istream& stream, const bool onlyProteinCoding);

// SNPCounter.cpp
std::unordered_map<std::string, std::vector<std::tuple<size_t, char, char>>> readVariants(std::istream& file);
std::unordered_map<std::string, std::vector<std::tuple<size_t, char, char>>> readVariantsVcf(const std::string& vcfFile);
std::string parseFilenameNoDirectory(const std::string& filename);

// PseudobulkAnalysis.cpp
std::vector<bool> getIncludedVariants(const EMOutput& output, const double minConfidence);
std::vector<bool> getIncludedCells(const EMOutput& output, const double minConfidence);
std::vector<std::unordered_map<std::string, std::tuple<size_t, size_t, size_t, size_t>>> getCellGroupedVariantPseudobulk(const EMOutput& output, const std::vector<CellMatch>& cellMatches, const double minConfidence, const std::unordered_map<std::string, std::string>& cellGrouping);

// EM.cpp
// These two structs must be token-for-token identical to the definitions in src/EM.cpp
struct EMCache
{
public:
	std::vector<double> variantMatScores;
	std::vector<double> variantPatScores;
	std::vector<double> variantXeMat;
	std::vector<double> variantXePat;
	std::vector<bool> hasVariant;
	std::vector<double> cellMatScores;
	std::vector<double> cellPatScores;
	std::vector<double> cellCeMat;
	std::vector<double> cellCePat;
	std::vector<bool> hasCell;
};

struct NoiseMaker
{
public:
	double getNoise()
	{
		return magnitude * normalDistribution(generator);
	}
	void initializeSeed(const size_t seed)
	{
		generator = std::mt19937 { seed };
	}
	double magnitude;
private:
	std::mt19937 generator;
	std::normal_distribution<> normalDistribution;
};

void initializeRandomly(EMResult& result, const std::unordered_map<size_t, bool>& forcedPhases, const size_t randomSeed);
double logprob(const size_t n, const double cellCoverageFraction, const double variantCoverage, const double cellEscapeFraction, const double variantEscapeFraction, const double weight, const bool active);
double getCellLogProbDerivative(const EMResult& result, const EMHelperVariables& helpers, const size_t cell, const double Ce, const bool matActive, const std::vector<bool>& ignoreTheseVariantsForNow);
double getCellLogProb(const EMResult& result, const EMHelperVariables& helpers, const size_t cell, const double Ce, const bool matActive, const std::vector<bool>& ignoreTheseVariantsForNow);
double getVariantLogProbDerivative(const EMResult& result, const EMHelperVariables& helpers, const size_t variant, const double Xe, const bool matRef);
double getVariantLogProbs(const EMResult& result, const EMHelperVariables& helpers, const size_t variant, const double Xe, const bool matRef);
double binarySearchOptimalCe(const EMResult& result, const EMHelperVariables& helpers, const size_t cell, const bool mat, const std::vector<bool>& ignoreTheseVariantsForNow);
std::pair<double, double> getOptimalCellCe(const EMResult& result, const EMHelperVariables& helpers, const size_t cell, const std::vector<bool>& ignoreTheseVariantsForNow);
double binarySearchOptimalXe(const EMResult& result, const EMHelperVariables& helpers, const size_t variant, const bool mat);
std::pair<double, double> getOptimalVariantXe(const EMResult& result, const EMHelperVariables& helpers, const size_t variant);
void merge(std::vector<size_t>& variantParent, size_t v1, size_t v2);
size_t find(std::vector<size_t>& variantParent, size_t variant);
bool flipUnsatisfiedSubgraph(EMResult& result, EMCache& cache, const EMHelperVariables& helpers, const std::unordered_map<size_t, bool>& forcedPhases);
bool maximizeVariantStates(EMResult& result, EMCache& cache, const std::unordered_map<size_t, bool>& forcedPhases, const EMHelperVariables& helpers, const std::vector<bool>& ignoreTheseVariantsForNow, NoiseMaker& noise);
bool maximizeCellStates(EMResult& result, EMCache& cache, const EMHelperVariables& helpers, const std::vector<bool>& ignoreTheseVariantsForNow, NoiseMaker& noise);
double getNonnormalizedTotalLogProb(const EMResult& result, const EMHelperVariables& helpers, const std::vector<bool>& ignoreTheseVariantsForNow);
double getTotalLogProb(const EMResult& result, const EMHelperVariables& helpers, const std::vector<bool>& ignoreTheseVariantsForNow);
double getTotalLogProb(const EMResult& result, const EMHelperVariables& helpers);
std::pair<double, double> getVariantEscapeConfidenceInterval(const EMResult& result, const EMHelperVariables& helpers, const size_t variantIndex, const double escapePointEstimate, const bool matRef, const double confidenceIntervalWidth);
std::pair<double, double> getCellEscapeConfidenceInterval(const EMResult& result, const EMHelperVariables& helpers, const size_t cellIndex, const double escapePointEstimate, const bool matActive, const double confidenceIntervalWidth, const std::vector<bool>& ignoreTheseVariantsForNow);
EMHelperVariables getHelpers(const std::vector<CellMatch>& cellMatches, const std::vector<std::tuple<size_t, size_t, double>>& regionWeights);
void getMaximumLikelihoodEM(EMResult& result, const std::vector<CellMatch>& cellMatches, const std::unordered_map<size_t, bool>& forcedPhases, const EMHelperVariables& helpers, const size_t noiseSeed, const double initialNoiseMagnitude, const double noiseDecay);
EMResult initializeResult(const EMHelperVariables& helpers);
EMResultAdditions getEMAdditions(const std::vector<CellMatch>& cellMatches, const EMResult& result, const EMHelperVariables& helpers);
std::unordered_map<size_t, bool> parseForcedPhases(const std::unordered_map<std::string, bool>& forcedPhases, const EMHelperVariables& helpers);

#endif
