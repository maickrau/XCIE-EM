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
const double maxEscape = 1.0;

struct CellMatch
{
public:
	std::string cell;
	std::string variant;
	bool alt;
	size_t count;
};

struct EMResult
{
public:
	std::unordered_map<std::string, size_t> variantNameToIndex;
	std::unordered_map<std::string, size_t> cellNameToIndex;
	std::vector<bool> cellIsMatActive;
	std::vector<double> cellEscapeFraction;
	std::vector<bool> variantIsMatRef;
	std::vector<double> variantEscapeFraction;
};

struct EMHelperVariables
{
public:
	std::vector<size_t> variantCoverage;
	std::vector<double> cellCoverageFraction;
	std::vector<std::unordered_map<size_t, size_t>> cellVariantRefCount;
	std::vector<std::unordered_map<size_t, size_t>> cellVariantAltCount;
	std::vector<std::vector<size_t>> activeCellsPerVariant;
	std::vector<std::vector<size_t>> activeVariantsPerCell;
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

double getCellLogProbDerivative(const EMResult& result, const EMHelperVariables& helpers, const size_t cell, const double Ce, const bool matActive)
{
	double derivativeSum = 0;
	const double f_j = helpers.cellCoverageFraction[cell];
	assert(f_j >= 0.0 - epsilon);
	assert(f_j <= 1.0 + epsilon);
	for (const size_t variant : helpers.activeVariantsPerCell[cell])
	{
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

double getCellLogProb(const EMResult& result, const EMHelperVariables& helpers, const size_t cell, const double Ce, const bool matActive)
{
	double logProbSum = 0;
	const double f_j = helpers.cellCoverageFraction[cell];
	assert(f_j >= 0.0 - epsilon);
	assert(f_j <= 1.0 + epsilon);
	for (const size_t variant : helpers.activeVariantsPerCell[cell])
	{
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

double binarySearchOptimalCe(const EMResult& result, const EMHelperVariables& helpers, const size_t cell, const bool mat)
{
	double min = escapeBoundary;
	double max = maxEscape - escapeBoundary;
	// find derivative zero
	for (size_t i = 0; i < 10; i++)
	{
		double mid = (min+max) / 2.0;
		double derivativeHere = getCellLogProbDerivative(result, helpers, cell, mid, mat);
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
		double CeScore = getCellLogProb(result, helpers, cell, Ce, mat);
		double biggerCeScore = getCellLogProb(result, helpers, cell, biggerCe, mat);
		if (biggerCeScore > CeScore) return biggerCe;
	}
	return Ce;
}

std::pair<double, double> getOptimalCellCe(const EMResult& result, const EMHelperVariables& helpers, const size_t cell)
{
	double matCe = binarySearchOptimalCe(result, helpers, cell, true);
	double patCe = binarySearchOptimalCe(result, helpers, cell, false);
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

bool maximizeVariantStates(EMResult& result, const std::unordered_map<size_t, bool>& forcedPhases, const EMHelperVariables& helpers)
{
	bool changed = false;
	size_t phasesChanged = 0;
	size_t escapeChanged = 0;
	for (size_t variant = 0; variant < result.variantIsMatRef.size(); variant++)
	{
		double matXe = 0;
		double patXe = 0;
		std::tie(matXe, patXe) = getOptimalVariantXe(result, helpers, variant);
		if (forcedPhases.count(variant) == 0)
		{
			double matRefLogProbSum = getVariantLogProbs(result, helpers, variant, matXe, true);
			double patRefLogProbSum = getVariantLogProbs(result, helpers, variant, patXe, false);
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
	std::cerr << phasesChanged << " variant phases changed, " << escapeChanged << " variant escapes changed" << std::endl;
	return changed;
}

bool maximizeCellStates(EMResult& result, const EMHelperVariables& helpers)
{
	bool changed = false;
	size_t cellsChanged = 0;
	size_t escapeChanged = 0;
	for (size_t cell = 0; cell < result.cellIsMatActive.size(); cell++)
	{
		double matCe, patCe;
		std::tie(matCe, patCe) = getOptimalCellCe(result, helpers, cell);
		double matActiveLogProbSum = getCellLogProb(result, helpers, cell, matCe, true);
		double patActiveLogProbSum = getCellLogProb(result, helpers, cell, patCe, false);
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
	std::cerr << cellsChanged << " cell actives changed, " << escapeChanged << " cell escapes changed" << std::endl;
	return changed;
}

double getNonnormalizedTotalLogProb(const EMResult& result, const EMHelperVariables& helpers)
{
	double total = 0;
	for (size_t variant = 0; variant < result.variantIsMatRef.size(); variant++)
	{
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

double getTotalLogProb(const EMResult& result, const EMHelperVariables& helpers)
{
	return getNonnormalizedTotalLogProb(result, helpers);
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
	std::vector<std::string> variantOrder = getVariantOrder(result.variantNameToIndex);
	std::vector<std::string> cellOrder = getCellOrder(result.cellNameToIndex);
	std::unordered_map<std::string, size_t> variantCoverage;
	std::unordered_map<std::string, size_t> cellCoverage;
	for (const auto& t : cellMatches)
	{
		variantCoverage[t.variant] += t.count;
		cellCoverage[t.cell] += t.count;
	}
	for (const std::string& variant : variantOrder)
	{
		const size_t variantIndex = result.variantNameToIndex.at(variant);
		const bool matRef = result.variantIsMatRef[variantIndex];
		double matXe, patXe;
		std::tie(matXe, patXe) = getOptimalVariantXe(result, helpers, variantIndex);
		double phaseScoreDifference = getVariantLogProbs(result, helpers, variantIndex, matRef ? matXe : patXe, matRef);
		phaseScoreDifference -= getVariantLogProbs(result, helpers, variantIndex, matRef ? patXe : matXe, !matRef);
		double escapeDifference = getVariantLogProbs(result, helpers, variantIndex, result.variantEscapeFraction.at(variantIndex), matRef);
		escapeDifference -= getVariantLogProbs(result, helpers, variantIndex, 0.5, matRef);
		stream << variant << "\t" << (matRef ? "mat" : "pat") << "\t" << result.variantEscapeFraction[variantIndex] << "\t" << variantCoverage.at(variant) << "\t" << phaseScoreDifference << "\t" << escapeDifference << std::endl;
	}
	for (const std::string& cell : cellOrder)
	{
		const size_t cellIndex = result.cellNameToIndex.at(cell);
		const bool matActive = result.cellIsMatActive[cellIndex];
		double matCe, patCe;
		std::tie(matCe, patCe) = getOptimalCellCe(result, helpers, cellIndex);
		double scoreDifference = getCellLogProb(result, helpers, cellIndex, matActive ? matCe : patCe, matActive);
		scoreDifference -= getCellLogProb(result, helpers, cellIndex, matActive ? patCe : matCe, !matActive);
		double escapeDifference = getCellLogProb(result, helpers, cellIndex, result.cellEscapeFraction[cellIndex], matActive);
		escapeDifference -= getCellLogProb(result, helpers, cellIndex, 0.5, matActive);
		stream << cell << "\t" << (matActive ? "mat" : "pat") << "\t" << result.cellEscapeFraction[cellIndex] << "\t" << cellCoverage.at(cell) << "\t" << scoreDifference << "\t" << escapeDifference << std::endl;
	}
}

EMHelperVariables getHelpers(const std::vector<CellMatch>& cellMatches, const EMResult& result)
{
	EMHelperVariables helpers;
	helpers.activeCellsPerVariant.resize(result.variantIsMatRef.size());
	helpers.activeVariantsPerCell.resize(result.cellIsMatActive.size());
	helpers.variantCoverage.resize(result.variantIsMatRef.size(), 0);
	helpers.cellVariantAltCount.resize(result.cellIsMatActive.size());
	helpers.cellVariantRefCount.resize(result.cellIsMatActive.size());
	helpers.cellCoverageFraction.resize(result.cellIsMatActive.size(), 0);
	for (const auto& t : cellMatches)
	{
		const size_t variantIndex = result.variantNameToIndex.at(t.variant);
		const size_t cellIndex = result.cellNameToIndex.at(t.cell);
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
		const size_t cellIndex = result.cellNameToIndex.at(t.cell);
		totalCoverage += t.count;
		helpers.cellCoverageFraction[cellIndex] += t.count;
	}
	for (size_t cell = 0; cell < helpers.cellCoverageFraction.size(); cell++)
	{
		helpers.cellCoverageFraction[cell] /= (double)totalCoverage;
	}
	return helpers;
}

void getMaximumLikelihoodEM(EMResult& result, const std::vector<CellMatch>& cellMatches, const std::unordered_map<size_t, bool>& forcedPhases, const EMHelperVariables& helpers, const size_t randomSeed)
{
	std::cerr << "initialize with random seed " << randomSeed << std::endl;
	assert(result.variantIsMatRef.size() == result.variantNameToIndex.size());
	assert(result.variantEscapeFraction.size() == result.variantNameToIndex.size());
	assert(result.cellIsMatActive.size() == result.cellNameToIndex.size());
	assert(result.cellEscapeFraction.size() == result.cellNameToIndex.size());
	for (size_t variant = 0; variant < result.variantIsMatRef.size(); variant++)
	{
		result.variantIsMatRef[variant] = false;
		result.variantEscapeFraction[variant] = -1;
	}
	for (size_t cell = 0; cell < result.cellIsMatActive.size(); cell++)
	{
		result.cellIsMatActive[cell] = false;
		result.cellEscapeFraction[cell] = -1;
	}
	initializeRandomly(result, forcedPhases, randomSeed);
	size_t overlapBetweenForcedVariantsAndAllVariants = 0;
	for (size_t variant = 0; variant < result.variantIsMatRef.size(); variant++)
	{
		if (forcedPhases.count(variant) == 1)
		{
			overlapBetweenForcedVariantsAndAllVariants += 1;
		}
	}
	std::cerr << overlapBetweenForcedVariantsAndAllVariants << " overlap between forced variants and all variants" << std::endl;
	size_t iteration = 0;
	double logprob = getTotalLogProb(result, helpers);
	std::cerr << "initial non-normalized log likelihood sum " << logprob << std::endl;
	while (true)
	{
		bool cellChanged = maximizeCellStates(result, helpers);
		if (cellChanged)
		{
			logprob = getTotalLogProb(result, helpers);
			std::cerr << "iteration " << iteration << " non-normalized log likelihood sum " << logprob << std::endl;
		}
		bool variantChanged = maximizeVariantStates(result, forcedPhases, helpers);
		if (variantChanged)
		{
			logprob = getTotalLogProb(result, helpers);
			std::cerr << "iteration " << iteration << " non-normalized log likelihood sum " << logprob << std::endl;
		}
		iteration += 1;
		if (!cellChanged && !variantChanged) break;
	}
	logprob = getTotalLogProb(result, helpers);
	std::cerr << "final non-normalized log likelihood sum " << logprob << std::endl;
}

std::unordered_map<size_t, bool> readForcedVariantPhases(const std::string& filename, const EMResult& EMresult)
{
	std::unordered_map<size_t, bool> forcedVariants;
	std::ifstream file { filename };
	while (file.good())
	{
		std::string line;
		getline(file, line);
		if (line.size() < 3) continue;
		std::stringstream sstr { line };
		std::string variant;
		std::string origin;
		sstr >> variant >> origin;
		if (EMresult.variantNameToIndex.count(variant) == 0) continue;
		if (origin == "mat")
		{
			forcedVariants[EMresult.variantNameToIndex.at(variant)] = true;
		}
		else
		{
			assert(origin == "pat");
			forcedVariants[EMresult.variantNameToIndex.at(variant)] = false;
		}
	}
	return forcedVariants;
}

EMResult initializeResult(const std::vector<CellMatch>& counts)
{
	EMResult result;
	for (const auto& t : counts)
	{
		if (result.variantNameToIndex.count(t.variant) == 0)
		{
			size_t index = result.variantNameToIndex.size();
			result.variantNameToIndex[t.variant] = index;
		}
		if (result.cellNameToIndex.count(t.cell) == 0)
		{
			size_t index = result.cellNameToIndex.size();
			result.cellNameToIndex[t.cell] = index;
		}
	}
	result.variantIsMatRef.resize(result.variantNameToIndex.size(), false);
	result.variantEscapeFraction.resize(result.variantNameToIndex.size(), -1);
	result.cellIsMatActive.resize(result.cellNameToIndex.size(), false);
	result.cellEscapeFraction.resize(result.cellNameToIndex.size(), -1);
	return result;
}

int main(int argc, char** argv)
{
	std::string matchTableFile { argv[1] };
	std::string forcedPhaseFile { argv[2] };
	size_t randomSeed = std::stoull(argv[3]);
	std::cerr << "read match counts" << std::endl;
	std::vector<CellMatch> counts = readMatchCounts(matchTableFile);
	std::cerr << "initialize" << std::endl;
	EMResult result = initializeResult(counts);
	std::cerr << "read forced variant phases" << std::endl;
	auto forcedPhases = readForcedVariantPhases(forcedPhaseFile, result);
	std::cerr << forcedPhases.size() << " forced variant phases" << std::endl;
	std::cerr << "get helper variables" << std::endl;
	EMHelperVariables helpers = getHelpers(counts, result);
	std::cerr << "run EM" << std::endl;
	getMaximumLikelihoodEM(result, counts, forcedPhases, helpers, randomSeed);
	std::cerr << "write results" << std::endl;
	writeResult(result, counts, helpers, std::cout);
}
