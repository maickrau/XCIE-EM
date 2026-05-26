#include <tuple>
#include <regex>
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
#include "Logger.h"
#include "EM.h"
#include "Common.h"
#include "FileHandler.h"

const double epsilon = 0.0001; // comparisons should be strict but add an epsilon because of float rounding issues
const double escapeBoundary = 0.001; // bounds X-chromosome inactivation escape to 0+this .. 1-this
const double maxEscape = 1.0;

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

size_t EMHelperVariables::numVariants() const
{
	return variantNameToIndex.size();
}

size_t EMHelperVariables::numCells() const
{
	return cellNameToIndex.size();
}

size_t getCount(const std::vector<std::unordered_map<size_t, size_t>>& cellVariantCount, const size_t cell, const size_t variant)
{
	if (cellVariantCount[cell].count(variant) == 0) return 0;
	return cellVariantCount[cell].at(variant);
}

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
	for (size_t i = 0; i < result.variantIsMatRef.size(); i++)
	{
		if (forcedPhases.count(i) == 1)
		{
			result.variantIsMatRef[i] = forcedPhases.at(i);
		}
		else
		{
			if (uniform(mt) < 0.5)
			{
				result.variantIsMatRef[i] = true;
			}
		}
	}
	for (size_t i = 0; i < result.variantEscapeFraction.size(); i++)
	{
		result.variantEscapeFraction[i] = escapeBoundary + (maxEscape - 2.0 * maxEscape * escapeBoundary) * uniform(mt);
	}
	for (size_t i = 0; i < result.cellEscapeFraction.size(); i++)
	{
		result.cellEscapeFraction[i] = escapeBoundary + (maxEscape - 2.0 * maxEscape * escapeBoundary) * uniform(mt);
	}
}

double logprob(const size_t n, const double cellCoverageFraction, const double variantCoverage, const double cellEscapeFraction, const double variantEscapeFraction, const bool active)
{
	assert(cellCoverageFraction >= 0.0 - epsilon);
	assert(cellCoverageFraction <= 1.0 + epsilon);
	assert(cellEscapeFraction >= 0.0 - epsilon);
	assert(cellEscapeFraction <= maxEscape + epsilon);
	assert(variantEscapeFraction >= 0.0 - epsilon);
	assert(variantEscapeFraction <= maxEscape + epsilon);
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
	assert(cellEscapeFraction <= maxEscape + epsilon);
	assert(variantEscapeFraction >= 0.0 - epsilon);
	assert(variantEscapeFraction <= maxEscape + epsilon);
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
	assert(cellEscapeFraction <= maxEscape + epsilon);
	assert(variantEscapeFraction >= 0.0 - epsilon);
	assert(variantEscapeFraction <= maxEscape + epsilon);
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
		assert(Xe <= maxEscape - escapeBoundary + epsilon);
		const size_t refCount = getCount(helpers.cellVariantRefCount, cell, variant);
		const size_t altCount = getCount(helpers.cellVariantAltCount, cell, variant);
		const bool activeMatchPhase = (result.variantIsMatRef[variant] == matActive);
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
		assert(Xe <= maxEscape - escapeBoundary + epsilon);
		const size_t c_i = helpers.variantCoverage[variant];
		const size_t refCount = getCount(helpers.cellVariantRefCount, cell, variant);
		const size_t altCount = getCount(helpers.cellVariantAltCount, cell, variant);
		const bool activeMatchPhase = (result.variantIsMatRef[variant] == matActive);
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
	double max = maxEscape - escapeBoundary;
	// find derivative zero
	for (size_t i = 0; i < 10; i++)
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
	if (Ce > maxEscape - escapeBoundary) Ce = maxEscape - escapeBoundary;
	double biggerCe = Ce + 0.001;
	if (biggerCe > maxEscape - escapeBoundary) biggerCe = maxEscape - escapeBoundary;
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

double binarySearchOptimalXe(const EMResult& result, const EMHelperVariables& helpers, const size_t variant, const bool mat)
{
	double min = escapeBoundary;
	double max = maxEscape - escapeBoundary;
	// find derivative zero
	for (size_t i = 0; i < 10; i++)
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

std::pair<double, double> getOptimalVariantXe(const EMResult& result, const EMHelperVariables& helpers, const size_t variant)
{
	double matXe = binarySearchOptimalXe(result, helpers, variant, true);
	double patXe = binarySearchOptimalXe(result, helpers, variant, false);
	return std::make_pair(matXe, patXe);
}

bool maximizeVariantStates(EMResult& result, const std::unordered_map<size_t, bool>& forcedPhases, const EMHelperVariables& helpers, const std::vector<bool>& ignoreTheseVariantsForNow, NoiseMaker& noise)
{
	bool changed = false;
	size_t phasesChanged = 0;
	size_t escapeChanged = 0;
	for (size_t variant = 0; variant < result.variantIsMatRef.size(); variant++)
	{
		if (ignoreTheseVariantsForNow[variant]) continue;
		double matXe = 0;
		double patXe = 0;
		std::tie(matXe, patXe) = getOptimalVariantXe(result, helpers, variant);
		if (forcedPhases.count(variant) == 0)
		{
			double matRefLogProbSum = getVariantLogProbs(result, helpers, variant, matXe, true) + noise.getNoise();
			double patRefLogProbSum = getVariantLogProbs(result, helpers, variant, patXe, false) + noise.getNoise();
			if (patRefLogProbSum > matRefLogProbSum + epsilon)
			{
				if (result.variantIsMatRef[variant])
				{
					result.variantIsMatRef[variant] = false;
					changed = true;
					phasesChanged += 1;
				}
			}
			if (matRefLogProbSum > patRefLogProbSum + epsilon)
			{
				if (!result.variantIsMatRef[variant])
				{
					result.variantIsMatRef[variant] = true;
					changed = true;
					phasesChanged += 1;
				}
			}
		}
		else
		{
			assert(result.variantIsMatRef[variant] == forcedPhases.at(variant));
		}
		if (result.variantIsMatRef[variant])
		{
			if (matXe > result.variantEscapeFraction[variant]+epsilon || matXe < result.variantEscapeFraction[variant]-epsilon)
			{
				changed = true;
				escapeChanged += 1;
			}
			result.variantEscapeFraction[variant] = matXe;
		}
		else
		{
			if (patXe > result.variantEscapeFraction[variant]+epsilon || patXe < result.variantEscapeFraction[variant]-epsilon)
			{
				changed = true;
				escapeChanged += 1;
			}
			result.variantEscapeFraction[variant] = patXe;
		}
	}
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << phasesChanged << " variant phases changed, " << escapeChanged << " variant escapes changed" << std::endl;
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
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << cellsChanged << " cell actives changed, " << escapeChanged << " cell escapes changed" << std::endl;
	return changed;
}

double getNonnormalizedTotalLogProb(const EMResult& result, const EMHelperVariables& helpers, const std::vector<bool>& ignoreTheseVariantsForNow)
{
	double total = 0;
	for (size_t variant = 0; variant < result.variantIsMatRef.size(); variant++)
	{
		if (ignoreTheseVariantsForNow[variant]) continue;
		const double Xe = result.variantEscapeFraction.at(variant);
		const double c_i = helpers.variantCoverage.at(variant);
		const bool variantIsMat = result.variantIsMatRef[variant];
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
	empty.resize(helpers.numVariants(), false);
	return getNonnormalizedTotalLogProb(result, helpers, empty);
}

std::pair<double, double> getVariantEscapeConfidenceInterval(const EMResult& result, const EMHelperVariables& helpers, const size_t variantIndex, const double escapePointEstimate, const bool matRef, const double confidenceIntervalWidth)
{
	double wantedScoreDifference = log((1.0-confidenceIntervalWidth)/2.0);
	double pointEstimateLogProb = getVariantLogProbs(result, helpers, variantIndex, escapePointEstimate, matRef);
	double minResult = escapeBoundary;
	double maxResult = maxEscape - escapeBoundary;
	for (int delta = 0; delta < 1000 && escapePointEstimate+(double)delta*0.001 < maxEscape-escapeBoundary; delta++)
	{
		double escape = escapePointEstimate + (double)delta * 0.001;
		double estimateHere = getVariantLogProbs(result, helpers, variantIndex, escape, matRef);
		if (estimateHere <= pointEstimateLogProb + wantedScoreDifference)
		{
			maxResult = escape;
			break;
		}
	}
	for (int delta = 0; delta < 1000 && escapePointEstimate-(double)delta*0.001 > escapeBoundary; delta++)
	{
		double escape = escapePointEstimate - (double)delta * 0.001;
		double estimateHere = getVariantLogProbs(result, helpers, variantIndex, escape, matRef);
		if (estimateHere <= pointEstimateLogProb + wantedScoreDifference)
		{
			minResult = escape;
			break;
		}
	}
	return std::make_pair(minResult, maxResult);
}

std::pair<double, double> getCellEscapeConfidenceInterval(const EMResult& result, const EMHelperVariables& helpers, const size_t cellIndex, const double escapePointEstimate, const bool matActive, const double confidenceIntervalWidth, const std::vector<bool>& ignoreTheseVariantsForNow)
{
	double wantedScoreDifference = log((1.0-confidenceIntervalWidth)/2.0);
	double pointEstimateLogProb = getCellLogProb(result, helpers, cellIndex, escapePointEstimate, matActive, ignoreTheseVariantsForNow);
	double minResult = escapeBoundary;
	double maxResult = maxEscape - escapeBoundary;
	for (int delta = 0; delta < 1000 && escapePointEstimate+(double)delta*0.001 < maxEscape-escapeBoundary; delta++)
	{
		double escape = escapePointEstimate + (double)delta * 0.001;
		double estimateHere = getCellLogProb(result, helpers, cellIndex, escape, matActive, ignoreTheseVariantsForNow);
		if (estimateHere <= pointEstimateLogProb + wantedScoreDifference)
		{
			maxResult = escape;
			break;
		}
	}
	for (int delta = 0; delta < 1000 && escapePointEstimate-(double)delta*0.001 > escapeBoundary; delta++)
	{
		double escape = escapePointEstimate - (double)delta * 0.001;
		double estimateHere = getCellLogProb(result, helpers, cellIndex, escape, matActive, ignoreTheseVariantsForNow);
		if (estimateHere <= pointEstimateLogProb + wantedScoreDifference)
		{
			minResult = escape;
			break;
		}
	}
	return std::make_pair(minResult, maxResult);
}

std::vector<PseudobulkVariantInfo> getVariantPseudobulk(const EMOutput& output, const std::vector<CellMatch>& cellMatches, const double minConfidence)
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
		const bool matRef = output.result.variantIsMatRef[variantIndex];
		double matXe, patXe;
		std::tie(matXe, patXe) = getOptimalVariantXe(output.result, output.helpers, variantIndex);
		double phaseScoreDifference = getVariantLogProbs(output.result, output.helpers, variantIndex, matRef ? matXe : patXe, matRef);
		phaseScoreDifference -= getVariantLogProbs(output.result, output.helpers, variantIndex, matRef ? patXe : matXe, !matRef);
		if (phaseScoreDifference < minConfidence) continue;
		includedVariants.insert(variantIndex);
	}
	std::unordered_set<size_t> includedCells;
	for (size_t cellIndex = 0; cellIndex < output.helpers.numCells(); cellIndex++)
	{
		const bool matActive = output.result.cellIsMatActive[cellIndex];
		double matCe, patCe;
		std::tie(matCe, patCe) = getOptimalCellCe(output.result, output.helpers, cellIndex, ignoreNothing);
		double scoreDifference = getCellLogProb(output.result, output.helpers, cellIndex, matActive ? matCe : patCe, matActive, ignoreNothing);
		scoreDifference -= getCellLogProb(output.result, output.helpers, cellIndex, matActive ? patCe : matCe, !matActive, ignoreNothing);
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
	std::vector<PseudobulkVariantInfo> result;
	result.reserve(output.helpers.numVariants());
	for (const std::string& name : variantOrder)
	{
		size_t index = output.helpers.variantNameToIndex.at(name);
		result.emplace_back();
		result.back().variantName = name;
		result.back().matXa = variantCoverageMatXa[index];
		result.back().matXi = variantCoverageMatXi[index];
		result.back().patXa = variantCoveragePatXa[index];
		result.back().patXi = variantCoveragePatXi[index];
	}
	return result;
}

EMHelperVariables getHelpers(const std::vector<CellMatch>& cellMatches)
{
	EMHelperVariables helpers;
	for (const auto& t : cellMatches)
	{
		if (helpers.variantNameToIndex.count(t.variant) == 0)
		{
			size_t index = helpers.numVariants();
			helpers.variantNameToIndex[t.variant] = index;
		}
		if (helpers.cellNameToIndex.count(t.cell) == 0)
		{
			size_t index = helpers.numCells();
			helpers.cellNameToIndex[t.cell] = index;
		}
	}
	helpers.activeCellsPerVariant.resize(helpers.numVariants());
	helpers.activeVariantsPerCell.resize(helpers.numCells());
	helpers.variantCoverage.resize(helpers.numVariants(), 0);
	helpers.cellCoverage.resize(helpers.numCells(), 0);
	helpers.cellVariantAltCount.resize(helpers.numCells());
	helpers.cellVariantRefCount.resize(helpers.numCells());
	helpers.cellCoverageFraction.resize(helpers.numCells(), 0);
	for (const auto& t : cellMatches)
	{
		const size_t variantIndex = helpers.variantNameToIndex.at(t.variant);
		const size_t cellIndex = helpers.cellNameToIndex.at(t.cell);
		helpers.activeCellsPerVariant[variantIndex].emplace_back(cellIndex);
		helpers.activeVariantsPerCell[cellIndex].emplace_back(variantIndex);
		helpers.variantCoverage[variantIndex] += t.count;
		helpers.cellCoverage[cellIndex] += t.count;
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

void getMaximumLikelihoodEM(EMResult& result, const std::vector<CellMatch>& cellMatches, const std::unordered_map<size_t, bool>& forcedPhases, const EMHelperVariables& helpers, const size_t noiseSeed, const double initialNoiseMagnitude, const double noiseDecay)
{
	assert(result.variantIsMatRef.size() == helpers.numVariants());
	assert(result.variantEscapeFraction.size() == helpers.numVariants());
	assert(result.cellIsMatActive.size() == helpers.numCells());
	assert(result.cellEscapeFraction.size() == helpers.numCells());
	size_t iteration = 0;
	double logprob = getTotalLogProb(result, helpers);
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "initial non-normalized log likelihood sum " << logprob << std::endl;
	NoiseMaker noise;
	noise.initializeSeed(noiseSeed);
	noise.magnitude = initialNoiseMagnitude;
	std::vector<bool> ignoreNothing;
	ignoreNothing.resize(helpers.numVariants(), false);
	while (true)
	{
		bool variantChanged = maximizeVariantStates(result, forcedPhases, helpers, ignoreNothing, noise);
		if (variantChanged)
		{
			logprob = getTotalLogProb(result, helpers, ignoreNothing);
			Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "iteration " << iteration << " non-normalized log likelihood sum " << logprob << std::endl;
		}
		bool cellChanged = maximizeCellStates(result, helpers, ignoreNothing, noise);
		if (cellChanged)
		{
			logprob = getTotalLogProb(result, helpers, ignoreNothing);
			Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "iteration " << iteration << " non-normalized log likelihood sum " << logprob << std::endl;
		}
		iteration += 1;
		if (!cellChanged && !variantChanged)
		{
			if (noise.magnitude == 0)
			{
				break;
			}
			noise.magnitude = 0;
			continue;
		}
		noise.magnitude *= noiseDecay;
	}
	logprob = getTotalLogProb(result, helpers);
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "final non-normalized log likelihood sum " << logprob << std::endl;
}

EMResult initializeResult(const EMHelperVariables& helpers)
{
	EMResult result;
	result.variantIsMatRef.resize(helpers.numVariants(), false);
	result.variantEscapeFraction.resize(helpers.numVariants(), -1);
	result.cellIsMatActive.resize(helpers.numCells(), false);
	result.cellEscapeFraction.resize(helpers.numCells(), -1);
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

EMResultAdditions getEMAdditions(const std::vector<CellMatch>& cellMatches, const EMResult& result, const EMHelperVariables& helpers)
{
	EMResultAdditions additions;
	additions.variantPhaseConfidence.resize(helpers.numVariants(), 0);
	additions.variantEscapeCILow.resize(helpers.numVariants(), 0);
	additions.variantEscapeCIHigh.resize(helpers.numVariants(), 0);
	additions.cellActiveConfidence.resize(helpers.numCells(), 0);
	additions.cellActiveConfidenceOnlyNonescape.resize(helpers.numCells(), 0);
	additions.cellEscapeCILow.resize(helpers.numCells(), 0);
	additions.cellEscapeCIHigh.resize(helpers.numCells(), 0);
	std::vector<bool> ignoreNothing;
	ignoreNothing.resize(helpers.numVariants(), false);
	std::vector<bool> ignoreEscapeVariants;
	ignoreEscapeVariants.resize(helpers.numVariants(), false);
	for (size_t i = 0 ; i < result.variantEscapeFraction.size(); i++)
	{
		if (result.variantEscapeFraction[i] <= escapeBoundary + epsilon) continue;
		ignoreEscapeVariants[i] = true;
	}
	for (size_t variantIndex = 0; variantIndex < helpers.numVariants(); variantIndex++)
	{
		const bool matRef = result.variantIsMatRef[variantIndex];
		double matXe, patXe;
		std::tie(matXe, patXe) = getOptimalVariantXe(result, helpers, variantIndex);
		double phaseScoreDifference = getVariantLogProbs(result, helpers, variantIndex, matRef ? matXe : patXe, matRef);
		phaseScoreDifference -= getVariantLogProbs(result, helpers, variantIndex, matRef ? patXe : matXe, !matRef);
		double escapeConfidenceIntervalMin, escapeConfidenceIntervalMax;
		std::tie(escapeConfidenceIntervalMin, escapeConfidenceIntervalMax) = getVariantEscapeConfidenceInterval(result, helpers, variantIndex, matRef ? matXe : patXe, matRef, 0.95);
		additions.variantPhaseConfidence[variantIndex] = phaseScoreDifference;
		additions.variantEscapeCILow[variantIndex] = escapeConfidenceIntervalMin;
		additions.variantEscapeCIHigh[variantIndex] = escapeConfidenceIntervalMax;
	}
	for (size_t cellIndex = 0; cellIndex < helpers.numCells(); cellIndex++)
	{
		const bool matActive = result.cellIsMatActive[cellIndex];
		double matCe, patCe;
		std::tie(matCe, patCe) = getOptimalCellCe(result, helpers, cellIndex, ignoreNothing);
		double scoreDifference = getCellLogProb(result, helpers, cellIndex, matActive ? matCe : patCe, matActive, ignoreNothing);
		scoreDifference -= getCellLogProb(result, helpers, cellIndex, matActive ? patCe : matCe, !matActive, ignoreNothing);
		double escapeConfidenceIntervalMin, escapeConfidenceIntervalMax;
		std::tie(escapeConfidenceIntervalMin, escapeConfidenceIntervalMax) = getCellEscapeConfidenceInterval(result, helpers, cellIndex, matActive ? matCe : patCe, matActive, 0.95, ignoreNothing);
		additions.cellActiveConfidence[cellIndex] = scoreDifference;
		additions.cellEscapeCILow[cellIndex] = escapeConfidenceIntervalMin;
		additions.cellEscapeCIHigh[cellIndex] = escapeConfidenceIntervalMax;
	}
	for (size_t cellIndex = 0; cellIndex < helpers.numCells(); cellIndex++)
	{
		const bool matActive = result.cellIsMatActive[cellIndex];
		double matCe, patCe;
		std::tie(matCe, patCe) = getOptimalCellCe(result, helpers, cellIndex, ignoreEscapeVariants);
		double scoreDifference = getCellLogProb(result, helpers, cellIndex, matActive ? matCe : patCe, matActive, ignoreEscapeVariants);
		scoreDifference -= getCellLogProb(result, helpers, cellIndex, matActive ? patCe : matCe, !matActive, ignoreEscapeVariants);
		additions.cellActiveConfidenceOnlyNonescape[cellIndex] = scoreDifference;
	}
	return additions;
}

std::unordered_map<size_t, bool> parseForcedPhases(const std::unordered_map<std::string, bool>& forcedPhases, const EMHelperVariables& helpers)
{
	std::unordered_map<size_t, bool> result;
	for (const auto& pair : forcedPhases)
	{
		if (helpers.variantNameToIndex.count(pair.first) == 0) continue;
		result[helpers.variantNameToIndex.at(pair.first)] = pair.second;
	}
	return result;
}

EMOutput runEM(const std::vector<CellMatch>& cellMatches, const std::unordered_map<std::string, bool>& forcedPhases, const size_t randomSeed, const double initialNoiseMagnitude, const double noiseDecay, const size_t numTries)
{
	double bestScore = -100000.0;
	Logger::Log.log(Logger::LogLevel::DebugInfo) << "get helper variables" << std::endl;
	EMHelperVariables helpers = getHelpers(cellMatches);
	EMResult bestResult;
	std::unordered_map<size_t, bool> parsedForcedPhases = parseForcedPhases(forcedPhases, helpers);
	for (size_t iteration = 0; iteration < numTries; iteration++)
	{
		Logger::Log.log(Logger::LogLevel::Always) << "EM run " << (iteration+1) << "/" << numTries << std::endl;
		Logger::Log.log(Logger::LogLevel::DebugInfo) << "initialize" << std::endl;
		EMResult result = initializeResult(helpers);
		size_t randomSeedHere = randomSeed + iteration;
		Logger::Log.log(Logger::LogLevel::DebugInfo) << "initialize with random seed " << randomSeedHere << std::endl;
		initializeRandomly(result, parsedForcedPhases, randomSeedHere);
		Logger::Log.log(Logger::LogLevel::DebugInfo) << "run EM" << std::endl;
		getMaximumLikelihoodEM(result, cellMatches, parsedForcedPhases, helpers, randomSeedHere, initialNoiseMagnitude, noiseDecay);
		double score = getTotalLogProb(result, helpers);
		Logger::Log.log(Logger::LogLevel::DebugInfo) << "got score " << score << std::endl;
		if (iteration == 0 || score > bestScore)
		{
			Logger::Log.log(Logger::LogLevel::DebugInfo) << "best so far" << std::endl;
			bestScore = score;
			bestResult = result;
		}
	}
	Logger::Log.log(Logger::LogLevel::Always) << "best score " << bestScore << std::endl;
	EMOutput output;
	std::swap(output.result, bestResult);
	std::swap(output.helpers, helpers);
	output.additions = getEMAdditions(cellMatches, output.result, output.helpers);
	return output;
}
