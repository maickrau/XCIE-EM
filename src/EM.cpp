#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <sstream>
#include <cassert>
#include <random>
#include <algorithm>

const double epsilon = 0.0001; // comparisons should be strict but add an epsilon because of float rounding issues
const double escapeBoundary = 0.001; // bounds X-chromosome inactivation escape to 0+this .. 1-this
const double variantEscapeBoundary = 0.5;
const double cellEscapeBoundary = 0.5;
const double inactiveExpressedPenalty = 3.0; // penalty for expressing inactive chrX more than active. only Xist should express inactive more than active. 3 corresponds to prior prob ~5%.

const size_t binarySearchRounds = 12;

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

struct CellMatch
{
public:
	std::string cell;
	std::string variant;
	bool alt;
	size_t count;
};

struct EMHelperVariables
{
public:
	std::unordered_map<std::string, size_t> variantNameToIndex;
	std::unordered_map<std::string, size_t> cellNameToIndex;
	std::vector<size_t> variantCoverage;
	std::vector<double> cellCoverageFraction;
	std::vector<std::unordered_map<size_t, size_t>> cellVariantRefCount;
	std::vector<std::unordered_map<size_t, size_t>> cellVariantAltCount;
	std::vector<std::vector<size_t>> activeCellsPerVariant;
	std::vector<std::vector<size_t>> activeVariantsPerCell;
	std::vector<std::vector<std::pair<size_t, bool>>> phaseBlocks;  // true: hap1 is alt. false: hap1 is ref
	std::vector<std::pair<size_t, bool>> variantPhaseBlock;         // true: hap1 is alt. false: hap1 is ref
};

struct EMResult
{
public:
	std::vector<bool> cellIsMatActive;
	std::vector<double> cellEscapeFraction;
	std::vector<double> variantEscapeFraction;
	std::vector<bool> phaseBlockIsHap1Mat;
	bool variantIsMatRef(const size_t variant, const EMHelperVariables& helpers) const
	{
		size_t block = helpers.variantPhaseBlock[variant].first;
		size_t allele = helpers.variantPhaseBlock[variant].second;
		return ~(phaseBlockIsHap1Mat[block] ^ allele);
	}
};

size_t getCount(const std::vector<std::unordered_map<size_t, size_t>>& cellVariantCount, const size_t cell, const size_t variant)
{
	if (cellVariantCount[cell].count(variant) == 0) return 0;
	return cellVariantCount[cell].at(variant);
}

std::vector<CellMatch> readMatchCounts(const std::string& matchTableFile)
{
	std::ifstream file { matchTableFile };
	std::vector<CellMatch> result;
	while (file.good())
	{
		std::string line;
		getline(file, line);
		std::stringstream sstr { line };
		if (line.size() < 7) continue;
		CellMatch match;
		std::string altstr;
		sstr >> match.cell >> match.variant >> altstr >> match.count;
		if (altstr == "ALT")
		{
			match.alt = true;
		}
		else
		{
			assert(altstr == "REF");
			match.alt = false;
		}
		result.emplace_back(match);
	}
	return result;
}

void initializeRandomly(EMResult& result, const std::unordered_map<size_t, bool>& forcedPhases, const size_t randomSeed)
{
	std::mt19937 mt(randomSeed);
	std::uniform_real_distribution<double> uniform(0, 1);
	for (size_t i = 0; i < result.cellIsMatActive.size(); i++)
	{
		if (uniform(mt) < 0.5)
		{
			result.cellIsMatActive[i] = true;
		}
	}
	for (size_t i = 0; i < result.phaseBlockIsHap1Mat.size(); i++)
	{
		if (forcedPhases.count(i) == 1)
		{
			result.phaseBlockIsHap1Mat[i] = forcedPhases.at(i);
		}
		else
		{
			if (uniform(mt) < 0.5)
			{
				result.phaseBlockIsHap1Mat[i] = true;
			}
		}
	}
	for (size_t i = 0; i < result.variantEscapeFraction.size(); i++)
	{
		result.variantEscapeFraction[i] = escapeBoundary + (2.0 - 2.0 * escapeBoundary) * uniform(mt);
	}
	for (size_t i = 0; i < result.cellEscapeFraction.size(); i++)
	{
		result.cellEscapeFraction[i] = escapeBoundary + (1.0 - 2.0 * escapeBoundary) * uniform(mt);
	}
}

double logprob(const size_t n, const double cellCoverageFraction, const double variantCoverage, const double cellEscapeFraction, const double variantEscapeFraction, const bool active)
{
	assert(cellCoverageFraction >= 0.0 - epsilon);
	assert(cellCoverageFraction <= 1.0 + epsilon);
	assert(cellEscapeFraction >= 0.0 - epsilon);
	assert(cellEscapeFraction <= 1.0 + epsilon);
	assert(variantEscapeFraction >= 0.0 - epsilon);
	assert(variantEscapeFraction <= 2.0 + epsilon);
	const double E = 1.0 - (1.0-cellEscapeFraction) * (1.0-variantEscapeFraction);
	double lambda;
	if (active)
	{
		double Ef = 1.0 - E * 0.5;
		lambda = cellCoverageFraction * variantCoverage * Ef;
	}
	else
	{
		double Ef = E * 0.5;
		lambda = cellCoverageFraction * variantCoverage * Ef;
	}
	if (n == 0) return -lambda;
	assert(lambda > 0);
	double result = n * log(lambda) - lambda;
	for (size_t i = 2; i <= n; i++)
	{
		result -= log(i);
	}
	assert(result < epsilon);
	return result;
}

// derivative by Ce
double logprobDerivativeCe(const size_t n, const double cellCoverageFraction, const double variantCoverage, const double cellEscapeFraction, const double variantEscapeFraction, const bool active)
{
	assert(cellCoverageFraction >= 0.0 - epsilon);
	assert(cellCoverageFraction <= 1.0 + epsilon);
	assert(cellEscapeFraction >= 0.0 - epsilon);
	assert(cellEscapeFraction <= 1.0 + epsilon);
	assert(variantEscapeFraction >= 0.0 - epsilon);
	assert(variantEscapeFraction <= 2.0 + epsilon);
	const double E = 1.0 - (1.0-cellEscapeFraction) * (1.0-variantEscapeFraction);
	double lambda;
	double lambdaDerivative;
	if (active)
	{
		double Ef = 1.0 - E * 0.5;
		lambda = cellCoverageFraction * variantCoverage * Ef;
		lambdaDerivative = cellCoverageFraction * variantCoverage * (-1) * (1.0 - variantEscapeFraction) / 2.0;
	}
	else
	{
		double Ef = E * 0.5;
		lambda = cellCoverageFraction * variantCoverage * Ef;
		lambdaDerivative = cellCoverageFraction * variantCoverage * (1.0 - variantEscapeFraction) / 2.0;
	}
	if (n == 0) return (-1) * lambdaDerivative;
	assert(lambda > 0);
	double result = ((double)n / lambda - 1.0) * lambdaDerivative;
	return result;
}

// derivative by Xe
double logprobDerivativeXe(const size_t n, const double cellCoverageFraction, const double variantCoverage, const double cellEscapeFraction, const double variantEscapeFraction, const bool active)
{
	assert(cellCoverageFraction >= 0.0 - epsilon);
	assert(cellCoverageFraction <= 1.0 + epsilon);
	assert(cellEscapeFraction >= 0.0 - epsilon);
	assert(cellEscapeFraction <= 1.0 + epsilon);
	assert(variantEscapeFraction >= 0.0 - epsilon);
	assert(variantEscapeFraction <= 2.0 + epsilon);
	const double E = 1.0 - (1.0-cellEscapeFraction) * (1.0-variantEscapeFraction);
	double lambda;
	double lambdaDerivative;
	if (active)
	{
		double Ef = 1.0 - E * 0.5;
		lambda = cellCoverageFraction * variantCoverage * Ef;
		lambdaDerivative = cellCoverageFraction * variantCoverage * (-1) * (1.0 - cellEscapeFraction) / 2.0;
	}
	else
	{
		double Ef = E * 0.5;
		lambda = cellCoverageFraction * variantCoverage * Ef;
		lambdaDerivative = cellCoverageFraction * variantCoverage * (1.0 - cellEscapeFraction) / 2.0;
	}
	if (n == 0) return (-1) * lambdaDerivative;
	assert(lambda > 0);
	double result = ((double)n / lambda - 1.0) * lambdaDerivative;
	return result;
}

double getCellLogProbDerivative(const EMResult& result, const EMHelperVariables& helpers, const size_t cell, const double Ce, const bool matActive, const std::vector<bool>& ignoreTheseVariantsForNow)
{
	double derivativeSum = 0;
	const double f_j = helpers.cellCoverageFraction[cell];
	assert(f_j >= 0.0 - epsilon);
	assert(f_j <= 1.0 + epsilon);
	for (const size_t variant : helpers.activeVariantsPerCell[cell])
	{
		if (ignoreTheseVariantsForNow[variant]) continue;
		const size_t c_i = helpers.variantCoverage[variant];
		const double Xe = result.variantEscapeFraction[variant];
		assert(Xe >= escapeBoundary - epsilon);
		assert(Xe <= 2.0 - escapeBoundary + epsilon);
		const size_t refCount = getCount(helpers.cellVariantRefCount, cell, variant);
		const size_t altCount = getCount(helpers.cellVariantAltCount, cell, variant);
		const bool activeMatchPhase = (result.variantIsMatRef(variant, helpers) == matActive);
		assert(refCount+altCount <= c_i);
		if (refCount > 0)
		{
			derivativeSum += logprobDerivativeCe(refCount, f_j, c_i, Ce, Xe, activeMatchPhase) - logprobDerivativeCe(0, f_j, c_i, Ce, Xe, activeMatchPhase);
		}
		if (altCount > 0)
		{
			derivativeSum += logprobDerivativeCe(altCount, f_j, c_i, Ce, Xe, !activeMatchPhase) - logprobDerivativeCe(0, f_j, c_i, Ce, Xe, !activeMatchPhase);
		}
	}
	return derivativeSum;
}

double getCellLogProb(const EMResult& result, const EMHelperVariables& helpers, const size_t cell, const double Ce, const bool matActive, const std::vector<bool>& ignoreTheseVariantsForNow)
{
	double logProbSum = 0;
	const double f_j = helpers.cellCoverageFraction[cell];
	assert(f_j >= 0.0 - epsilon);
	assert(f_j <= 1.0 + epsilon);
	for (const size_t variant : helpers.activeVariantsPerCell[cell])
	{
		if (ignoreTheseVariantsForNow[variant]) continue;
		const double Xe = result.variantEscapeFraction[variant];
		assert(Xe >= escapeBoundary - epsilon);
		assert(Xe <= 2.0 - escapeBoundary + epsilon);
		const size_t c_i = helpers.variantCoverage[variant];
		const size_t refCount = getCount(helpers.cellVariantRefCount, cell, variant);
		const size_t altCount = getCount(helpers.cellVariantAltCount, cell, variant);
		const bool activeMatchPhase = (result.variantIsMatRef(variant, helpers) == matActive);
		assert(refCount+altCount <= c_i);
		if (refCount > 0)
		{
			logProbSum += logprob(refCount, f_j, c_i, Ce, Xe, activeMatchPhase) - logprob(0, f_j, c_i, Ce, Xe, activeMatchPhase);
		}
		if (altCount > 0)
		{
			logProbSum += logprob(altCount, f_j, c_i, Ce, Xe, !activeMatchPhase) - logprob(0, f_j, c_i, Ce, Xe, !activeMatchPhase);
		}
	}
	return logProbSum;
}

double getVariantLogProbDerivative(const EMResult& result, const EMHelperVariables& helpers, const size_t variant, const double Xe, const bool matRef)
{
	const size_t c_i = helpers.variantCoverage[variant];
	double derivativeSum = 0;
	for (const size_t cell : helpers.activeCellsPerVariant[variant])
	{
		const double f_j = helpers.cellCoverageFraction[cell];
		const double Ce = result.cellEscapeFraction[cell];
		const size_t refCount = getCount(helpers.cellVariantRefCount, cell, variant);
		const size_t altCount = getCount(helpers.cellVariantAltCount, cell, variant);
		const bool activeMatchPhase = (result.cellIsMatActive[cell] == matRef);
		if (refCount > 0)
		{
			derivativeSum += logprobDerivativeXe(refCount, f_j, c_i, Ce, Xe, activeMatchPhase) - logprobDerivativeXe(0, f_j, c_i, Ce, Xe, activeMatchPhase);
		}
		if (altCount > 0)
		{
			derivativeSum += logprobDerivativeXe(altCount, f_j, c_i, Ce, Xe, !activeMatchPhase) - logprobDerivativeXe(0, f_j, c_i, Ce, Xe, !activeMatchPhase);
		}
	}
	return derivativeSum;
}

double getVariantLogProbs(const EMResult& result, const EMHelperVariables& helpers, const size_t variant, const double Xe, const bool matRef)
{
	const size_t c_i = helpers.variantCoverage.at(variant);
	double logProbSum = 0;
	for (const size_t cell : helpers.activeCellsPerVariant[variant])
	{
		const double f_j = helpers.cellCoverageFraction[cell];
		const double Ce = result.cellEscapeFraction[cell];
		const size_t refCount = getCount(helpers.cellVariantRefCount, cell, variant);
		const size_t altCount = getCount(helpers.cellVariantAltCount, cell, variant);
		const bool activeMatchPhase = (result.cellIsMatActive[cell] == matRef);
		if (refCount > 0)
		{
			logProbSum += logprob(refCount, f_j, c_i, Ce, Xe, activeMatchPhase) - logprob(0, f_j, c_i, Ce, Xe, activeMatchPhase);
		}
		if (altCount > 0)
		{
			logProbSum += logprob(altCount, f_j, c_i, Ce, Xe, !activeMatchPhase) - logprob(0, f_j, c_i, Ce, Xe, !activeMatchPhase);
		}
	}
	return logProbSum;
}

double binarySearchOptimalCe(const EMResult& result, const EMHelperVariables& helpers, const size_t cell, const bool mat, const std::vector<bool>& ignoreTheseVariantsForNow)
{
	double min = escapeBoundary;
	double max = 1.0 - escapeBoundary;
	// find derivative zero
	for (size_t i = 0; i < binarySearchRounds; i++)
	{
		double mid = (min+max) / 2.0;
		double derivativeHere = getCellLogProbDerivative(result, helpers, cell, mid, mat, ignoreTheseVariantsForNow);
		if (derivativeHere > 0)
		{
			min = mid;
		}
		else
		{
			max = mid;
		}
	}
	assert(max-min < 0.001);
	// round down to nearest 0.001
	double Ce = (int)(1000.0 * (max+min)/2.0) * 0.001;
	if (Ce < escapeBoundary) Ce = escapeBoundary;
	if (Ce > 1.0 - escapeBoundary) Ce = 1.0 - escapeBoundary;
	double biggerCe = Ce + 0.001;
	if (biggerCe > 1.0 - escapeBoundary) biggerCe = 1.0 - escapeBoundary;
	if (biggerCe != Ce)
	{
		double CeScore = getCellLogProb(result, helpers, cell, Ce, mat, ignoreTheseVariantsForNow);
		double biggerCeScore = getCellLogProb(result, helpers, cell, biggerCe, mat, ignoreTheseVariantsForNow);
		if (biggerCeScore > CeScore) return biggerCe;
	}
	return Ce;
}

std::pair<double, double> getOptimalCellCe(const EMResult& result, const EMHelperVariables& helpers, const size_t cell, const std::vector<bool>& ignoreTheseVariantsForNow)
{
	double matCe = binarySearchOptimalCe(result, helpers, cell, true, ignoreTheseVariantsForNow);
	double patCe = binarySearchOptimalCe(result, helpers, cell, false, ignoreTheseVariantsForNow);
	return std::make_pair(matCe, patCe);
}

double binarySearchOptimalXe(const EMResult& result, const EMHelperVariables& helpers, const size_t variant, const bool mat, const double maxEscape)
{
	double min = escapeBoundary;
	double max = maxEscape - escapeBoundary;
	// find derivative zero
	for (size_t i = 0; i < binarySearchRounds; i++)
	{
		double mid = (min+max) / 2.0;
		double derivativeHere = getVariantLogProbDerivative(result, helpers, variant, mid, mat);
		if (derivativeHere > 0)
		{
			min = mid;
		}
		else
		{
			max = mid;
		}
	}
	assert(max-min < 0.001);
	// round down to nearest 0.001
	double Xe = (int)(1000.0 * (max+min)/2.0) * 0.001;
	if (Xe < escapeBoundary) Xe = escapeBoundary;
	if (Xe > maxEscape - escapeBoundary) Xe = maxEscape - escapeBoundary;
	double biggerXe = Xe + 0.001;
	if (biggerXe > maxEscape - escapeBoundary) biggerXe = maxEscape - escapeBoundary;
	if (biggerXe != Xe)
	{
		double XeScore = getVariantLogProbs(result, helpers, variant, Xe, mat);
		double biggerXeScore = getVariantLogProbs(result, helpers, variant, biggerXe, mat);
		if (biggerXeScore >= XeScore) return biggerXe;
	}
	return Xe;
}

std::pair<double, double> getOptimalVariantXe(const EMResult& result, const EMHelperVariables& helpers, const size_t variant, const double maxEscape)
{
	double matXe = binarySearchOptimalXe(result, helpers, variant, true, maxEscape);
	double patXe = binarySearchOptimalXe(result, helpers, variant, false, maxEscape);
	return std::make_pair(matXe, patXe);
}

bool maximizeVariantStates(EMResult& result, const std::unordered_map<size_t, bool>& forcedPhases, const EMHelperVariables& helpers, const std::vector<bool>& ignoreTheseVariantsForNow, NoiseMaker& noise)
{
	bool changed = false;
	size_t phasesChanged = 0;
	size_t escapeChanged = 0;
	for (size_t block = 0; block < result.phaseBlockIsHap1Mat.size(); block++)
	{
		if (forcedPhases.count(block) == 0)
		{
			double blockHap1MatLogProbSum = 0;
			double blockHap1PatLogProbSum = 0;
			for (const auto& pair : helpers.phaseBlocks[block])
			{
				size_t variant = pair.first;
				if (ignoreTheseVariantsForNow[variant]) continue;
				bool varianthap1IsAlt = pair.second;
				double maxEscape = 1;
				if (helpers.phaseBlocks[block].size() > 1) maxEscape = 2;
				double matXe, patXe;
				std::tie(matXe, patXe) = getOptimalVariantXe(result, helpers, variant, maxEscape);
				double matRefLogProbSum = getVariantLogProbs(result, helpers, variant, matXe, true) + noise.getNoise();
				double patRefLogProbSum = getVariantLogProbs(result, helpers, variant, patXe, false) + noise.getNoise();
				if (!varianthap1IsAlt)
				{
					// variant hap1 is ref
					blockHap1MatLogProbSum += matRefLogProbSum;
					if (matXe > 1.0) blockHap1MatLogProbSum -= inactiveExpressedPenalty;
					blockHap1PatLogProbSum += patRefLogProbSum;
					if (patXe > 1.0) blockHap1PatLogProbSum -= inactiveExpressedPenalty;
				}
				else
				{
					blockHap1MatLogProbSum += patRefLogProbSum;
					if (patXe > 1.0) blockHap1MatLogProbSum -= inactiveExpressedPenalty;
					blockHap1PatLogProbSum += matRefLogProbSum;
					if (matXe > 1.0) blockHap1PatLogProbSum -= inactiveExpressedPenalty;
				}
			}
			if (blockHap1PatLogProbSum > blockHap1MatLogProbSum + epsilon)
			{
				if (result.phaseBlockIsHap1Mat[block])
				{
					result.phaseBlockIsHap1Mat[block] = false;
					changed = true;
					phasesChanged += 1;
				}
			}
			if (blockHap1MatLogProbSum > blockHap1PatLogProbSum + epsilon)
			{
				if (!result.phaseBlockIsHap1Mat[block])
				{
					result.phaseBlockIsHap1Mat[block] = true;
					changed = true;
					phasesChanged += 1;
				}
			}
		}
		else
		{
			assert(result.phaseBlockIsHap1Mat[block] == forcedPhases.at(block));
		}
	}
	for (size_t variant = 0; variant < helpers.variantNameToIndex.size(); variant++)
	{
		if (ignoreTheseVariantsForNow[variant]) continue;
		size_t block = helpers.variantPhaseBlock[variant].first;
		double maxEscape = 1;
		if (helpers.phaseBlocks[block].size() > 1) maxEscape = 2;
		double matXe, patXe;
		std::tie(matXe, patXe) = getOptimalVariantXe(result, helpers, variant, maxEscape);
		if (result.variantIsMatRef(variant, helpers))
		{
			if (matXe < result.variantEscapeFraction[variant] - epsilon || matXe > result.variantEscapeFraction[variant] + epsilon)
			{
				result.variantEscapeFraction[variant] = matXe;
				escapeChanged += 1;
			}
		}
		else
		{
			if (patXe < result.variantEscapeFraction[variant] - epsilon || patXe > result.variantEscapeFraction[variant] + epsilon)
			{
				result.variantEscapeFraction[variant] = patXe;
				escapeChanged += 1;
			}
		}
	}
//	std::cerr << phasesChanged << " variant phases changed, " << escapeChanged << " variant escapes changed" << std::endl;
	return changed;
}

bool maximizeCellStates(EMResult& result, const EMHelperVariables& helpers, const std::vector<bool>& ignoreTheseVariantsForNow, NoiseMaker& noise)
{
	bool changed = false;
	size_t cellsChanged = 0;
	size_t escapeChanged = 0;
	for (size_t cell = 0; cell < result.cellIsMatActive.size(); cell++)
	{
		double matCe, patCe;
		std::tie(matCe, patCe) = getOptimalCellCe(result, helpers, cell, ignoreTheseVariantsForNow);
		double matActiveLogProbSum = getCellLogProb(result, helpers, cell, matCe, true, ignoreTheseVariantsForNow) + noise.getNoise();
		double patActiveLogProbSum = getCellLogProb(result, helpers, cell, patCe, false, ignoreTheseVariantsForNow) + noise.getNoise();
		// should be strict comparison but add epsilon because of floating point rounding
		if (matActiveLogProbSum > patActiveLogProbSum + epsilon)
		{
			if (!result.cellIsMatActive[cell])
			{
				changed = true;
				result.cellIsMatActive[cell] = true;
				cellsChanged += 1;
			}
		}
		if (patActiveLogProbSum > matActiveLogProbSum + epsilon)
		{
			if (result.cellIsMatActive[cell])
			{
				changed = true;
				result.cellIsMatActive[cell] = false;
				cellsChanged += 1;
			}
		}
		if (result.cellIsMatActive[cell])
		{
			if (matCe > result.cellEscapeFraction[cell] + epsilon || matCe < result.cellEscapeFraction[cell] - epsilon)
			{
				changed = true;
				escapeChanged += 1;
			}
			result.cellEscapeFraction[cell] = matCe;
		}
		else
		{
			if (patCe > result.cellEscapeFraction[cell] + epsilon || patCe < result.cellEscapeFraction[cell] - epsilon)
			{
				changed = true;
				escapeChanged += 1;
			}
			result.cellEscapeFraction[cell] = patCe;
		}
	}
//	std::cerr << cellsChanged << " cell actives changed, " << escapeChanged << " cell escapes changed" << std::endl;
	return changed;
}

double getNonnormalizedTotalLogProb(const EMResult& result, const EMHelperVariables& helpers, const std::vector<bool>& ignoreTheseVariantsForNow)
{
	double total = 0;
	for (size_t variant = 0; variant < helpers.variantNameToIndex.size(); variant++)
	{
		if (ignoreTheseVariantsForNow[variant]) continue;
		const double Xe = result.variantEscapeFraction.at(variant);
		const double c_i = helpers.variantCoverage.at(variant);
		const bool variantIsMat = result.variantIsMatRef(variant, helpers);
		for (size_t cell : helpers.activeCellsPerVariant[variant])
		{
			const size_t refCount = getCount(helpers.cellVariantRefCount, cell, variant);
			const size_t altCount = getCount(helpers.cellVariantAltCount, cell, variant);
			const double Ce = result.cellEscapeFraction[cell];
			const bool cellIsMat = result.cellIsMatActive[cell];
			const double f_j = helpers.cellCoverageFraction[cell];
			if (refCount > 0) total += logprob(refCount, f_j, c_i, Ce, Xe, variantIsMat == cellIsMat) - logprob(0, f_j, c_i, Ce, Xe, variantIsMat == cellIsMat);
			if (altCount > 0) total += logprob(altCount, f_j, c_i, Ce, Xe, variantIsMat != cellIsMat) - logprob(0, f_j, c_i, Ce, Xe, variantIsMat != cellIsMat);
		}
	}
	return total;
}

double getTotalLogProb(const EMResult& result, const EMHelperVariables& helpers, const std::vector<bool>& ignoreTheseVariantsForNow)
{
	return getNonnormalizedTotalLogProb(result, helpers, ignoreTheseVariantsForNow);
}

double getTotalLogProb(const EMResult& result, const EMHelperVariables& helpers)
{
	std::vector<bool> empty;
	empty.resize(helpers.variantNameToIndex.size(), false);
	return getNonnormalizedTotalLogProb(result, helpers, empty);
}

std::vector<std::string> getVariantOrder(const std::unordered_map<std::string, size_t>& nameToIndex)
{
	std::vector<std::string> result;
	for (const auto& pair : nameToIndex)
	{
		result.emplace_back(pair.first);
	}
	std::sort(result.begin(), result.end(), [](const std::string& left, const std::string& right)
	{
		size_t leftPos = 0;
		for (size_t i = 0; i < left.size(); i++)
		{
			if (left[i] != ':') continue;
			for (size_t j = i+1; j < left.size(); j++)
			{
				if (left[j] != ':') continue;
				leftPos = std::stoull(left.substr(i+1, j-i-1));
				break;
			}
			break;
		}
		size_t rightPos = 0;
		for (size_t i = 0; i < right.size(); i++)
		{
			if (right[i] != ':') continue;
			for (size_t j = i+1; j < right.size(); j++)
			{
				if (right[j] != ':') continue;
				rightPos = std::stoull(right.substr(i+1, j-i-1));
				break;
			}
			break;
		}
		return leftPos < rightPos;
	});
	return result;
}

std::vector<std::string> getCellOrder(const std::unordered_map<std::string, size_t>& nameToIndex)
{
	std::vector<std::string> result;
	for (const auto& pair : nameToIndex)
	{
		result.emplace_back(pair.first);
	}
	std::sort(result.begin(), result.end());
	return result;
}

void writeResult(const EMResult& result, const std::vector<CellMatch>& cellMatches, const EMHelperVariables& helpers, std::ostream& stream)
{
	std::vector<std::string> variantOrder = getVariantOrder(helpers.variantNameToIndex);
	std::vector<std::string> cellOrder = getCellOrder(helpers.cellNameToIndex);
	std::unordered_map<std::string, size_t> variantCoverage;
	std::unordered_map<std::string, size_t> cellCoverage;
	std::vector<bool> ignoreNothing;
	ignoreNothing.resize(helpers.variantNameToIndex.size(), false);
	for (const auto& t : cellMatches)
	{
		variantCoverage[t.variant] += t.count;
		cellCoverage[t.cell] += t.count;
	}
	for (const std::string& variant : variantOrder)
	{
		const size_t variantIndex = helpers.variantNameToIndex.at(variant);
		const bool matRef = result.variantIsMatRef(variantIndex, helpers);
		double matXe, patXe;
		size_t block = helpers.variantPhaseBlock[variantIndex].first;
		double maxEscape = 1;
		if (helpers.phaseBlocks[block].size() > 1) maxEscape = 2;
		std::tie(matXe, patXe) = getOptimalVariantXe(result, helpers, variantIndex, maxEscape);
		double phaseScoreDifference = getVariantLogProbs(result, helpers, variantIndex, matRef ? matXe : patXe, matRef);
		phaseScoreDifference -= getVariantLogProbs(result, helpers, variantIndex, matRef ? patXe : matXe, !matRef);
		double escapeDifference = getVariantLogProbs(result, helpers, variantIndex, result.variantEscapeFraction.at(variantIndex), matRef);
		escapeDifference -= getVariantLogProbs(result, helpers, variantIndex, variantEscapeBoundary, matRef);
		stream << variant << "\t" << (matRef ? "mat" : "pat") << "\t" << result.variantEscapeFraction[variantIndex] << "\t" << variantCoverage.at(variant) << "\t" << phaseScoreDifference << "\t" << escapeDifference << std::endl;
	}
	for (const std::string& cell : cellOrder)
	{
		const size_t cellIndex = helpers.cellNameToIndex.at(cell);
		const bool matActive = result.cellIsMatActive[cellIndex];
		double matCe, patCe;
		std::tie(matCe, patCe) = getOptimalCellCe(result, helpers, cellIndex, ignoreNothing);
		double scoreDifference = getCellLogProb(result, helpers, cellIndex, matActive ? matCe : patCe, matActive, ignoreNothing);
		scoreDifference -= getCellLogProb(result, helpers, cellIndex, matActive ? patCe : matCe, !matActive, ignoreNothing);
		double escapeDifference = getCellLogProb(result, helpers, cellIndex, result.cellEscapeFraction[cellIndex], matActive, ignoreNothing);
		escapeDifference -= getCellLogProb(result, helpers, cellIndex, cellEscapeBoundary, matActive, ignoreNothing);
		stream << cell << "\t" << (matActive ? "mat" : "pat") << "\t" << result.cellEscapeFraction[cellIndex] << "\t" << cellCoverage.at(cell) << "\t" << scoreDifference << "\t" << escapeDifference << std::endl;
	}
}

EMHelperVariables getHelpers(const std::vector<CellMatch>& cellMatches)
{
	EMHelperVariables helpers;
	for (const auto& t : cellMatches)
	{
		if (helpers.variantNameToIndex.count(t.variant) == 0)
		{
			size_t index = helpers.variantNameToIndex.size();
			helpers.variantNameToIndex[t.variant] = index;
		}
		if (helpers.cellNameToIndex.count(t.cell) == 0)
		{
			size_t index = helpers.cellNameToIndex.size();
			helpers.cellNameToIndex[t.cell] = index;
		}
	}
	helpers.activeCellsPerVariant.resize(helpers.variantNameToIndex.size());
	helpers.activeVariantsPerCell.resize(helpers.cellNameToIndex.size());
	helpers.variantCoverage.resize(helpers.variantNameToIndex.size(), 0);
	helpers.cellVariantAltCount.resize(helpers.cellNameToIndex.size());
	helpers.cellVariantRefCount.resize(helpers.cellNameToIndex.size());
	helpers.cellCoverageFraction.resize(helpers.cellNameToIndex.size(), 0);
	for (const auto& t : cellMatches)
	{
		const size_t variantIndex = helpers.variantNameToIndex.at(t.variant);
		const size_t cellIndex = helpers.cellNameToIndex.at(t.cell);
		helpers.activeCellsPerVariant[variantIndex].emplace_back(cellIndex);
		helpers.activeVariantsPerCell[cellIndex].emplace_back(variantIndex);
		helpers.variantCoverage[variantIndex] += t.count;
		if (t.alt)
		{
			helpers.cellVariantAltCount[cellIndex][variantIndex] += t.count;
		}
		else
		{
			helpers.cellVariantRefCount[cellIndex][variantIndex] += t.count;
		}
	}
	for (size_t i = 0; i < helpers.activeCellsPerVariant.size(); i++)
	{
		std::unordered_set<size_t> uniques { helpers.activeCellsPerVariant[i].begin(), helpers.activeCellsPerVariant[i].end() };
		helpers.activeCellsPerVariant[i].clear();
		helpers.activeCellsPerVariant[i].insert(helpers.activeCellsPerVariant[i].end(), uniques.begin(), uniques.end());
	}
	for (size_t i = 0; i < helpers.activeVariantsPerCell.size(); i++)
	{
		std::unordered_set<size_t> uniques { helpers.activeVariantsPerCell[i].begin(), helpers.activeVariantsPerCell[i].end() };
		helpers.activeVariantsPerCell[i].clear();
		helpers.activeVariantsPerCell[i].insert(helpers.activeVariantsPerCell[i].end(), uniques.begin(), uniques.end());
	}
	size_t totalCoverage = 0;
	for (const auto& t : cellMatches)
	{
		const size_t cellIndex = helpers.cellNameToIndex.at(t.cell);
		totalCoverage += t.count;
		helpers.cellCoverageFraction[cellIndex] += t.count;
	}
	for (size_t cell = 0; cell < helpers.cellCoverageFraction.size(); cell++)
	{
		helpers.cellCoverageFraction[cell] /= (double)totalCoverage;
	}
	return helpers;
}

size_t parseVariantPosition(const std::string& name)
{
	size_t firstColon = name.size();
	for (size_t i = 0; i < name.size(); i++)
	{
		if (name[i] == ':')
		{
			if (firstColon == name.size())
			{
				firstColon = i;
			}
			else
			{
				return std::stoull(name.substr(firstColon+1, i-(firstColon+1)));
			}
		}
	}
	return 0;
}

std::vector<bool> getIgnoredSecondStepVariants(const EMHelperVariables& helpers, const std::vector<size_t>& secondStepVariants)
{
	std::vector<bool> ignored;
	ignored.resize(helpers.variantNameToIndex.size(), false);
	for (size_t variant : secondStepVariants)
	{
		ignored[variant] = true;
	}
	return ignored;
}

void getMaximumLikelihoodEM(EMResult& result, const std::vector<CellMatch>& cellMatches, const std::unordered_map<size_t, bool>& forcedPhases, const EMHelperVariables& helpers, const std::vector<size_t>& secondStepVariants, const size_t noiseSeed, const double initialNoiseMagnitude, const double noiseDecay)
{
//	assert(result.variantIsMatRef.size() == helpers.variantNameToIndex.size());
	assert(result.variantEscapeFraction.size() == helpers.variantNameToIndex.size());
	assert(result.cellIsMatActive.size() == helpers.cellNameToIndex.size());
	assert(result.cellEscapeFraction.size() == helpers.cellNameToIndex.size());
	size_t iteration = 0;
	double logprob = getTotalLogProb(result, helpers);
	std::vector<bool> ignoreTheseVariantsForNow = getIgnoredSecondStepVariants(helpers, secondStepVariants);
//	std::cerr << ignoreTheseVariantsForNow.size() << " variants ignored in first step" << std::endl;
//	std::cerr << "initial non-normalized log likelihood sum " << logprob << std::endl;
	NoiseMaker noise;
	noise.initializeSeed(noiseSeed);
	noise.magnitude = initialNoiseMagnitude;
	while (true)
	{
		bool variantChanged = maximizeVariantStates(result, forcedPhases, helpers, ignoreTheseVariantsForNow, noise);
		if (variantChanged)
		{
			logprob = getTotalLogProb(result, helpers, ignoreTheseVariantsForNow);
//			std::cerr << "iteration " << iteration << " non-normalized log likelihood sum " << logprob << std::endl;
		}
		bool cellChanged = maximizeCellStates(result, helpers, ignoreTheseVariantsForNow, noise);
		if (cellChanged)
		{
			logprob = getTotalLogProb(result, helpers, ignoreTheseVariantsForNow);
//			std::cerr << "iteration " << iteration << " non-normalized log likelihood sum " << logprob << std::endl;
		}
		iteration += 1;
		if (!cellChanged && !variantChanged)
		{
			if (ignoreTheseVariantsForNow.size() == 0)
			{
				if (noise.magnitude == 0)
				{
					break;
				}
				noise.magnitude = 0;
				continue;
			}
			ignoreTheseVariantsForNow.clear();

//			std::cerr << "switch to second step" << std::endl;
//			logprob = getTotalLogProb(result, helpers, ignoreTheseVariantsForNow);
//			std::cerr << "iteration " << iteration << " non-normalized log likelihood sum " << logprob << std::endl;
		}
		noise.magnitude *= noiseDecay;
	}
	logprob = getTotalLogProb(result, helpers);
	std::cerr << "final non-normalized log likelihood sum " << logprob << std::endl;
}

std::unordered_map<size_t, bool> readForcedVariantPhases(const std::string& filename, const EMHelperVariables& helpers)
{
	std::unordered_map<size_t, bool> forcedVariants;
	std::ifstream file { filename };
	size_t totalForcedVariants = 0;
	while (file.good())
	{
		std::string line;
		getline(file, line);
		if (line.size() < 3) continue;
		std::stringstream sstr { line };
		std::string variant;
		std::string origin;
		sstr >> variant >> origin;
		totalForcedVariants += 1;
		if (helpers.variantNameToIndex.count(variant) == 0) continue;
		if (origin == "mat")
		{
			forcedVariants[helpers.variantNameToIndex.at(variant)] = true;
		}
		else
		{
			assert(origin == "pat");
			forcedVariants[helpers.variantNameToIndex.at(variant)] = false;
		}
	}
//	std::cerr << "total forced variants " << totalForcedVariants << ", overlap with real variants " << forcedVariants.size() << std::endl;
	return forcedVariants;
}

EMResult initializeResult(const EMHelperVariables& helpers)
{
	EMResult result;
	result.phaseBlockIsHap1Mat.resize(helpers.phaseBlocks.size(), false);
//	result.variantIsMatRef.resize(helpers.variantNameToIndex.size(), false);
	result.variantEscapeFraction.resize(helpers.variantNameToIndex.size(), -1);
	result.cellIsMatActive.resize(helpers.cellNameToIndex.size(), false);
	result.cellEscapeFraction.resize(helpers.cellNameToIndex.size(), -1);
	return result;
}

std::vector<std::pair<std::string, std::string>> parseTags(const std::string& tags)
{
	std::vector<std::pair<std::string, std::string>> result;
	size_t start = 0;
	size_t equal = 0;
	for (size_t i = 0; i < tags.size(); i++)
	{
		if (tags[i] == ';')
		{
			result.emplace_back(tags.substr(start, equal-start), tags.substr(equal+1, i-(equal+1)));
			start = i+1;
			equal = i;
		}
		if (tags[i] == '=')
		{
			equal = i;
		}
	}
	result.emplace_back(tags.substr(start, equal-start), tags.substr(equal+1));
	return result;
}

std::vector<std::pair<size_t, size_t>> loadSecondStepRegions(const std::string& annotationFile, const std::string& secondStepGeneList)
{
	if (annotationFile == "/dev/null") return std::vector<std::pair<size_t, size_t>> {};
	if (secondStepGeneList == "/dev/null") return std::vector<std::pair<size_t, size_t>> {};
	std::unordered_set<std::string> secondStepGeneNames;
	{
		std::ifstream file { secondStepGeneList };
		while (file.good())
		{
			std::string gene;
			file >> gene;
			if (gene.size() > 0) secondStepGeneNames.insert(gene);
		}
	}
//	std::cerr << secondStepGeneNames.size() << " second step genes" << std::endl;
	std::vector<std::pair<size_t, size_t>> result;
	if (secondStepGeneNames.size() == 0) return std::vector<std::pair<size_t, size_t>> {};
	{
		std::ifstream file { annotationFile };
		while (file.good())
		{
			std::string line;
			std::getline(file, line);
			std::stringstream sstr { line };
			if (line.size() < 5) continue;
			if (line[0] == '#') continue;
			std::string chromosome, source, type;
			sstr >> chromosome >> source >> type;
			if (chromosome != "chrX" && chromosome != "23" && chromosome != "X") continue;
			if (type != "gene") continue;
			size_t startpos, endpos;
			std::string dummy, tags;
			sstr >> startpos >> endpos >> dummy >> dummy >> dummy >> tags;
			for (std::pair<std::string, std::string> tag : parseTags(tags))
			{
				if (tag.first != "gene_name") continue;
				if (secondStepGeneNames.count(tag.second) == 0) continue;
				result.emplace_back(startpos, endpos);
			}
		}
	}
//	std::cerr << result.size() << " second step regions" << std::endl;
	return result;
}

std::vector<size_t> loadSecondStepVariants(const std::string& annotationFile, const std::string& secondStepGeneList, const EMHelperVariables& helpers)
{
	auto regions = loadSecondStepRegions(annotationFile, secondStepGeneList);
	std::vector<size_t> result;
	if (regions.size() == 0) return result;
	for (const auto& pair : helpers.variantNameToIndex)
	{
		size_t parsedPosition = parseVariantPosition(pair.first);
		for (auto region : regions)
		{
			if (region.first <= parsedPosition && region.second >= parsedPosition)
			{
				result.emplace_back(pair.second);
				break;
			}
		}
	}
	return result;
}

void addUnphasedVariantsToIndividualPhaseBlocks(EMHelperVariables& helpers)
{
	std::vector<bool> presentInPhaseBlock;
	presentInPhaseBlock.resize(helpers.variantNameToIndex.size(), false);
	for (size_t i = 0; i < helpers.phaseBlocks.size(); i++)
	{
		for (const std::pair<size_t, bool>& pair : helpers.phaseBlocks[i])
		{
			size_t variant = pair.first;
			assert(!presentInPhaseBlock[variant]);
			presentInPhaseBlock[variant] = true;
		}
	}
	size_t countAdded = 0;
	for (size_t i = 0; i < presentInPhaseBlock.size(); i++)
	{
		if (presentInPhaseBlock[i]) continue;
		helpers.phaseBlocks.emplace_back();
		helpers.phaseBlocks.back().emplace_back(i, false);
		countAdded += 1;
	}
	std::cerr << countAdded << " variants without phase blocks" << std::endl;
}

std::vector<std::string> split(const std::string& str, const char separator)
{
	std::vector<std::string> result;
	size_t lastStart = 0;
	for (size_t i = 0; i < str.size(); i++)
	{
		if (str[i] != separator) continue;
		result.emplace_back(str.substr(lastStart, i-lastStart));
		lastStart = i+1;
	}
	result.emplace_back(str.substr(lastStart));
	return result;
}

std::vector<std::vector<std::pair<size_t, bool>>> readVCFPhaseBlocks(const EMHelperVariables& helpers, const std::string phasingVCF)
{
	std::ifstream file { phasingVCF };
	std::unordered_map<size_t, size_t> variantPhaseBlock;
	std::unordered_map<size_t, bool> variantPhaseAllele;
	while (file.good())
	{
		std::string line;
		std::getline(file, line);
		if (line.size() < 5) continue;
		if (line[0] == '@' || line[0] == '#') continue;
		auto parts = split(line, '\t');
		if (parts[0] != "23" && parts[0] != "chrX" && parts[0] != "X") continue;
		if (parts[6] != "PASS") continue;
		auto tags = split(parts[8], ':');
		size_t phaseBlockIndex = tags.size();
		size_t genotypeIndex = tags.size();
		for (size_t i = 0; i < tags.size(); i++)
		{
			if (tags[i] == "PL")
			{
				phaseBlockIndex = i;
			}
			if (tags[i] == "GT")
			{
				genotypeIndex = i;
			}
		}
		if (phaseBlockIndex == tags.size()) continue;
		if (genotypeIndex == tags.size()) continue;
		auto values = split(parts[9], ':');
		if (values.size() != tags.size()) continue;
		std::string genotype = values[genotypeIndex];
		if (genotype != "0|1" && genotype != "1|0") continue;
		size_t phaseBlock = std::numeric_limits<size_t>::max();
		phaseBlock = std::stoull(values[phaseBlockIndex]);
		bool phaseAllele = genotype[0] == '1';
		std::string variantName = "X:" + parts[1] + ":" + parts[3] + ":" + parts[4];
		if (helpers.variantNameToIndex.count(variantName) == 0) continue;
		size_t variantIndex = helpers.variantNameToIndex.at(variantName);
		variantPhaseBlock[variantIndex] = phaseBlock;
		variantPhaseAllele[variantIndex] = phaseAllele;
	}
	std::unordered_map<size_t, size_t> phaseBlockRenaming;
	for (const auto& pair : variantPhaseBlock)
	{
		if (phaseBlockRenaming.count(pair.second) == 1) continue;
		size_t newName = phaseBlockRenaming.size();
		phaseBlockRenaming[pair.second] = newName;
	}
	std::vector<std::vector<std::pair<size_t, bool>>> result;
	result.resize(phaseBlockRenaming.size());
	for (const auto& pair : variantPhaseBlock)
	{
		assert(variantPhaseAllele.count(pair.first) == 1);
		assert(phaseBlockRenaming.count(pair.second) == 1);
		result[phaseBlockRenaming.at(pair.second)].emplace_back(pair.first, variantPhaseAllele.at(pair.first));
	}
	std::cerr << result.size() << " phase blocks with " << variantPhaseBlock.size() << " variants" << std::endl;
	return result;
}

std::vector<std::pair<size_t, bool>> parsePhaseBlocks(const EMHelperVariables& helpers, const std::vector<std::vector<std::pair<size_t, bool>>>& parsePhaseBlocks)
{
	std::vector<std::pair<size_t, bool>> result;
	result.resize(helpers.variantNameToIndex.size(), std::make_pair(std::numeric_limits<size_t>::max(), true));
	for (size_t i = 0; i < parsePhaseBlocks.size(); i++)
	{
		for (const auto& pair : parsePhaseBlocks[i])
		{
			assert(pair.first < result.size());
			result[pair.first] = std::make_pair(i, pair.second);
		}
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		assert(result[i].first < parsePhaseBlocks.size());
	}
	return result;
}

void readReadPhasing(EMHelperVariables& helpers, const std::string phasingVCF)
{
	helpers.phaseBlocks = readVCFPhaseBlocks(helpers, phasingVCF);
	addUnphasedVariantsToIndividualPhaseBlocks(helpers);
	helpers.variantPhaseBlock = parsePhaseBlocks(helpers, helpers.phaseBlocks);
}

int main(int argc, char** argv)
{
	std::string matchTableFile { argv[1] };
	std::string forcedPhaseFile { argv[2] };
	size_t randomSeed = std::stoull(argv[3]);
	size_t numTries = std::stoull(argv[4]);
	std::string annotationFile { argv[5] };
	std::string secondStepGeneList { argv[6] };
	std::string readPhasingFile { argv[7] };
	double initialNoiseMagnitude = std::stod(argv[8]);
	double noiseDecay = std::stod(argv[9]);
//	std::cerr << "read match counts" << std::endl;
	std::vector<CellMatch> counts = readMatchCounts(matchTableFile);
//	std::cerr << "get helper variables" << std::endl;
	EMHelperVariables helpers = getHelpers(counts);
//	std::cerr << "read forced variant phases" << std::endl;
	auto forcedPhases = readForcedVariantPhases(forcedPhaseFile, helpers);
	readReadPhasing(helpers, readPhasingFile);
	std::vector<size_t> iterationsWithGoodScore;
	double bestScore = -100000.0;
	EMResult bestResult;
//	std::cerr << "read second step genes" << std::endl;
	std::vector<size_t> secondStepVariants = loadSecondStepVariants(annotationFile, secondStepGeneList, helpers);
//	std::cerr << secondStepVariants.size() << " second step variants" << std::endl;
	for (size_t iteration = 0; iteration < numTries; iteration++)
	{
		std::cerr << "bigiteration " << iteration << std::endl;
//		std::cerr << "initialize" << std::endl;
		EMResult result = initializeResult(helpers);
		size_t randomSeedHere = randomSeed + iteration;
//		std::cerr << "initialize with random seed " << randomSeedHere << std::endl;
		initializeRandomly(result, forcedPhases, randomSeedHere);
//		std::cerr << "run EM" << std::endl;
		getMaximumLikelihoodEM(result, counts, forcedPhases, helpers, secondStepVariants, randomSeedHere, initialNoiseMagnitude, noiseDecay);
//		std::cerr << "write results" << std::endl;
		double score = getTotalLogProb(result, helpers);
		if (iteration == 0 || score > bestScore)
		{
			bestScore = score;
			bestResult = result;
		}
	}
	std::cerr << "best score " << bestScore << std::endl;
	writeResult(bestResult, counts, helpers, std::cout);
}
