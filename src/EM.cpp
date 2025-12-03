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
	std::unordered_map<std::string, bool> cellIsMatActive;
	std::unordered_map<std::string, double> cellEscapeFraction;
	std::unordered_map<std::string, bool> variantIsMatRef;
	std::unordered_map<std::string, double> variantEscapeFraction;
};

struct EMHelperVariables
{
public:
	std::unordered_map<std::string, size_t> variantCoverage;
	std::unordered_map<std::string, double> cellCoverageFraction;
	std::unordered_map<std::string, std::unordered_map<std::string, size_t>> cellVariantRefCount;
	std::unordered_map<std::string, std::unordered_map<std::string, size_t>> cellVariantAltCount;
	std::unordered_map<std::string, std::unordered_set<std::string>> activeCellsPerVariant;
	std::unordered_map<std::string, std::unordered_set<std::string>> activeVariantsPerCell;
};

size_t getCount(const std::unordered_map<std::string, std::unordered_map<std::string, size_t>>& cellVariantCount, const std::string& cell, const std::string& variant)
{
	if (cellVariantCount.count(cell) == 0) return 0;
	if (cellVariantCount.at(cell).count(variant) == 0) return 0;
	return cellVariantCount.at(cell).at(variant);
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

void initializeRandomly(EMResult& result, const std::unordered_map<std::string, bool>& forcedPhases, const size_t randomSeed)
{
	std::mt19937 mt(randomSeed);
	std::uniform_real_distribution<double> uniform(0, 1);
	for (auto& pair : result.cellIsMatActive)
	{
		if (uniform(mt) < 0.5)
		{
			pair.second = true;
		}
	}
	for (auto& pair : result.variantIsMatRef)
	{
		if (forcedPhases.count(pair.first) == 1)
		{
			pair.second = forcedPhases.at(pair.first);
		}
		else
		{
			if (uniform(mt) < 0.5)
			{
				pair.second = true;
			}
		}
	}
	for (auto& pair : result.variantEscapeFraction)
	{
		pair.second = escapeBoundary + (maxEscape - 2.0 * maxEscape * escapeBoundary) * uniform(mt);
	}
	for (auto& pair : result.cellEscapeFraction)
	{
		pair.second = escapeBoundary + (maxEscape - 2.0 * maxEscape * escapeBoundary) * uniform(mt);
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

double getCellLogProbDerivative(const EMResult& result, const std::vector<CellMatch>& cellMatches, const EMHelperVariables& helpers, const std::string& cell, const double Ce, const bool matActive)
{
	double derivativeSum = 0;
	const double f_j = helpers.cellCoverageFraction.at(cell);
	assert(f_j >= 0.0 - epsilon);
	assert(f_j <= 1.0 + epsilon);
	for (const std::string& variant : helpers.activeVariantsPerCell.at(cell))
	{
		const size_t c_i = helpers.variantCoverage.at(variant);
		assert(result.variantEscapeFraction.count(variant) == 1);
		const double Xe = result.variantEscapeFraction.at(variant);
		assert(Xe >= escapeBoundary - epsilon);
		assert(Xe <= maxEscape - escapeBoundary + epsilon);
		const size_t refCount = getCount(helpers.cellVariantRefCount, cell, variant);
		const size_t altCount = getCount(helpers.cellVariantAltCount, cell, variant);
		const bool activeMatchPhase = (result.variantIsMatRef.at(variant) == matActive);
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

double getCellLogProb(const EMResult& result, const std::vector<CellMatch>& cellMatches, const EMHelperVariables& helpers, const std::string& cell, const double Ce, const bool matActive)
{
	double logProbSum = 0;
	const double f_j = helpers.cellCoverageFraction.at(cell);
	assert(f_j >= 0.0 - epsilon);
	assert(f_j <= 1.0 + epsilon);
	for (const std::string& variant : helpers.activeVariantsPerCell.at(cell))
	{
		assert(result.variantEscapeFraction.count(variant) == 1);
		const double Xe = result.variantEscapeFraction.at(variant);
		assert(Xe >= escapeBoundary - epsilon);
		assert(Xe <= maxEscape - escapeBoundary + epsilon);
		const size_t c_i = helpers.variantCoverage.at(variant);
		const size_t refCount = getCount(helpers.cellVariantRefCount, cell, variant);
		const size_t altCount = getCount(helpers.cellVariantAltCount, cell, variant);
		const bool activeMatchPhase = (result.variantIsMatRef.at(variant) == matActive);
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

double getVariantLogProbDerivative(const EMResult& result, const std::vector<CellMatch>& cellMatches, const EMHelperVariables& helpers, const std::string& variant, const double Xe, const bool matRef)
{
	const size_t c_i = helpers.variantCoverage.at(variant);
	double derivativeSum = 0;
	for (const std::string& cell : helpers.activeCellsPerVariant.at(variant))
	{
		const double f_j = helpers.cellCoverageFraction.at(cell);
		const double Ce = result.cellEscapeFraction.at(cell);
		const size_t refCount = getCount(helpers.cellVariantRefCount, cell, variant);
		const size_t altCount = getCount(helpers.cellVariantAltCount, cell, variant);
		const bool activeMatchPhase = (result.cellIsMatActive.at(cell) == matRef);
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

double getVariantLogProbs(const EMResult& result, const std::vector<CellMatch>& cellMatches, const EMHelperVariables& helpers, const std::string& variant, const double Xe, const bool matRef)
{
	const size_t c_i = helpers.variantCoverage.at(variant);
	double logProbSum = 0;
	for (const std::string& cell : helpers.activeCellsPerVariant.at(variant))
	{
		const double f_j = helpers.cellCoverageFraction.at(cell);
		const double Ce = result.cellEscapeFraction.at(cell);
		const size_t refCount = getCount(helpers.cellVariantRefCount, cell, variant);
		const size_t altCount = getCount(helpers.cellVariantAltCount, cell, variant);
		const bool activeMatchPhase = (result.cellIsMatActive.at(cell) == matRef);
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

double binarySearchOptimalCe(const EMResult& result, const std::vector<CellMatch>& cellMatches, const EMHelperVariables& helpers, const std::string& cell, const bool mat)
{
	double min = escapeBoundary;
	double max = maxEscape - escapeBoundary;
	// find derivative zero
	for (size_t i = 0; i < 10; i++)
	{
		double mid = (min+max) / 2.0;
		double derivativeHere = getCellLogProbDerivative(result, cellMatches, helpers, cell, mid, mat);
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
		double CeScore = getCellLogProb(result, cellMatches, helpers, cell, Ce, mat);
		double biggerCeScore = getCellLogProb(result, cellMatches, helpers, cell, biggerCe, mat);
		if (biggerCeScore > CeScore) return biggerCe;
	}
	return Ce;
}

std::pair<double, double> getOptimalCellCe(const EMResult& result, const std::vector<CellMatch>& cellMatches, const EMHelperVariables& helpers, const std::string& cell)
{
	double matCe = binarySearchOptimalCe(result, cellMatches, helpers, cell, true);
	double patCe = binarySearchOptimalCe(result, cellMatches, helpers, cell, false);
	return std::make_pair(matCe, patCe);
}

double binarySearchOptimalXe(const EMResult& result, const std::vector<CellMatch>& cellMatches, const EMHelperVariables& helpers, const std::string& variant, const bool mat)
{
	double min = escapeBoundary;
	double max = maxEscape - escapeBoundary;
	// find derivative zero
	for (size_t i = 0; i < 10; i++)
	{
		double mid = (min+max) / 2.0;
		double derivativeHere = getVariantLogProbDerivative(result, cellMatches, helpers, variant, mid, mat);
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
		double XeScore = getVariantLogProbs(result, cellMatches, helpers, variant, Xe, mat);
		double biggerXeScore = getVariantLogProbs(result, cellMatches, helpers, variant, biggerXe, mat);
		if (biggerXeScore >= XeScore) return biggerXe;
	}
	return Xe;
}

std::pair<double, double> getOptimalVariantXe(const EMResult& result, const std::vector<CellMatch>& cellMatches, const EMHelperVariables& helpers, const std::string& variant)
{
	double matXe = binarySearchOptimalXe(result, cellMatches, helpers, variant, true);
	double patXe = binarySearchOptimalXe(result, cellMatches, helpers, variant, false);
	return std::make_pair(matXe, patXe);
}

bool maximizeVariantStates(EMResult& result, const std::unordered_map<std::string, bool>& forcedPhases, const std::vector<CellMatch>& cellMatches, const EMHelperVariables& helpers)
{
	bool changed = false;
	size_t phasesChanged = 0;
	size_t escapeChanged = 0;
	for (auto& pair : result.variantIsMatRef)
	{
		const std::string& variant = pair.first;
		double matXe = 0;
		double patXe = 0;
		std::tie(matXe, patXe) = getOptimalVariantXe(result, cellMatches, helpers, variant);
		if (forcedPhases.count(variant) == 0)
		{
			double matRefLogProbSum = getVariantLogProbs(result, cellMatches, helpers, variant, matXe, true);
			double patRefLogProbSum = getVariantLogProbs(result, cellMatches, helpers, variant, patXe, false);
			if (patRefLogProbSum > matRefLogProbSum + epsilon)
			{
				if (pair.second)
				{
					pair.second = false;
					changed = true;
					phasesChanged += 1;
				}
			}
			if (matRefLogProbSum > patRefLogProbSum + epsilon)
			{
				if (!pair.second)
				{
					pair.second = true;
					changed = true;
					phasesChanged += 1;
				}
			}
		}
		else
		{
			assert(pair.second == forcedPhases.at(variant));
		}
		if (pair.second)
		{
			if (matXe > result.variantEscapeFraction.at(pair.first)+epsilon || matXe < result.variantEscapeFraction.at(pair.first)-epsilon)
			{
				changed = true;
				escapeChanged += 1;
			}
			result.variantEscapeFraction[pair.first] = matXe;
		}
		else
		{
			if (patXe > result.variantEscapeFraction.at(pair.first)+epsilon || patXe < result.variantEscapeFraction.at(pair.first)-epsilon)
			{
				changed = true;
				escapeChanged += 1;
			}
			result.variantEscapeFraction[pair.first] = patXe;
		}
	}
	std::cerr << phasesChanged << " variant phases changed, " << escapeChanged << " variant escapes changed" << std::endl;
	return changed;
}

bool maximizeCellStates(EMResult& result, const std::vector<CellMatch>& cellMatches, const EMHelperVariables& helpers)
{
	bool changed = false;
	size_t cellsChanged = 0;
	size_t escapeChanged = 0;
	for (auto& pair : result.cellIsMatActive)
	{
		const std::string& cell = pair.first;
		double matCe, patCe;
		std::tie(matCe, patCe) = getOptimalCellCe(result, cellMatches, helpers, cell);
		double matActiveLogProbSum = getCellLogProb(result, cellMatches, helpers, cell, matCe, true);
		double patActiveLogProbSum = getCellLogProb(result, cellMatches, helpers, cell, patCe, false);
		// should be strict comparison but add epsilon because of floating point rounding
		if (matActiveLogProbSum > patActiveLogProbSum + epsilon)
		{
			if (!pair.second)
			{
				changed = true;
				pair.second = true;
				cellsChanged += 1;
			}
		}
		if (patActiveLogProbSum > matActiveLogProbSum + epsilon)
		{
			if (pair.second)
			{
				changed = true;
				pair.second = false;
				cellsChanged += 1;
			}
		}
		if (pair.second)
		{
			if (matCe > result.cellEscapeFraction.at(cell) + epsilon || matCe < result.cellEscapeFraction.at(cell) - epsilon)
			{
				changed = true;
				escapeChanged += 1;
			}
			result.cellEscapeFraction[cell] = matCe;
		}
		else
		{
			if (patCe > result.cellEscapeFraction.at(cell) + epsilon || patCe < result.cellEscapeFraction.at(cell) - epsilon)
			{
				changed = true;
				escapeChanged += 1;
			}
			result.cellEscapeFraction[cell] = patCe;
		}
	}
	std::cerr << cellsChanged << " cell statuses changed, " << escapeChanged << " cell escapes changed" << std::endl;
	return changed;
}

double getNonnormalizedTotalLogProb(const EMResult& result, const std::vector<CellMatch>& cellMatches, const EMHelperVariables& helpers)
{
	double total = 0;
	for (const auto& pair : result.variantIsMatRef)
	{
		const std::string& variant = pair.first;
		assert(result.variantEscapeFraction.count(variant) == 1);
		assert(helpers.variantCoverage.count(variant) == 1);
		const double Xe = result.variantEscapeFraction.at(variant);
		const double c_i = helpers.variantCoverage.at(variant);
		for (const std::string& cell : helpers.activeCellsPerVariant.at(variant))
		{
			assert(helpers.cellCoverageFraction.count(cell) == 1);
			const size_t refCount = getCount(helpers.cellVariantRefCount, cell, variant);
			const size_t altCount = getCount(helpers.cellVariantAltCount, cell, variant);
			const double Ce = result.cellEscapeFraction.at(cell);
			const bool cellIsMat = result.cellIsMatActive.at(cell);
			const double f_j = helpers.cellCoverageFraction.at(cell);
			if (refCount > 0) total += logprob(refCount, f_j, c_i, Ce, Xe, pair.second == cellIsMat) - logprob(0, f_j, c_i, Ce, Xe, pair.second == cellIsMat);
			if (altCount > 0) total += logprob(altCount, f_j, c_i, Ce, Xe, pair.second != cellIsMat) - logprob(0, f_j, c_i, Ce, Xe, pair.second != cellIsMat);
		}
	}
	return total;
}

double getTotalLogProb(const EMResult& result, const std::vector<CellMatch>& cellMatches, const EMHelperVariables& helpers)
{
	return getNonnormalizedTotalLogProb(result, cellMatches, helpers);
}

std::vector<std::string> getVariantOrder(const std::unordered_map<std::string, bool>& matRef)
{
	std::vector<std::string> result;
	for (const auto& pair : matRef)
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

std::vector<std::string> getCellOrder(const std::unordered_map<std::string, bool>& matActive)
{
	std::vector<std::string> result;
	for (const auto& pair : matActive)
	{
		result.emplace_back(pair.first);
	}
	std::sort(result.begin(), result.end());
	return result;
}

void writeResult(const EMResult& result, const std::vector<CellMatch>& cellMatches, const EMHelperVariables& helpers, std::ostream& stream)
{
	std::vector<std::string> variantOrder = getVariantOrder(result.variantIsMatRef);
	std::vector<std::string> cellOrder = getCellOrder(result.cellIsMatActive);
	std::unordered_map<std::string, size_t> variantCoverage;
	std::unordered_map<std::string, size_t> cellCoverage;
	for (const auto& t : cellMatches)
	{
		variantCoverage[t.variant] += t.count;
		cellCoverage[t.cell] += t.count;
	}
	for (const std::string& variant : variantOrder)
	{
		bool matRef = result.variantIsMatRef.at(variant);
		double matXe, patXe;
		std::tie(matXe, patXe) = getOptimalVariantXe(result, cellMatches, helpers, variant);
		double phaseScoreDifference = getVariantLogProbs(result, cellMatches, helpers, variant, matRef ? matXe : patXe, matRef);
		phaseScoreDifference -= getVariantLogProbs(result, cellMatches, helpers, variant, matRef ? patXe : matXe, !matRef);
		double escapeDifference = getVariantLogProbs(result, cellMatches, helpers, variant, result.variantEscapeFraction.at(variant), matRef);
		escapeDifference -= getVariantLogProbs(result, cellMatches, helpers, variant, 0.5, matRef);
		stream << variant << "\t" << (matRef ? "mat" : "pat") << "\t" << result.variantEscapeFraction.at(variant) << "\t" << variantCoverage.at(variant) << "\t" << phaseScoreDifference << "\t" << escapeDifference << std::endl;
	}
	for (const std::string& cell : cellOrder)
	{
		bool matActive = result.cellIsMatActive.at(cell);
		double matCe, patCe;
		std::tie(matCe, patCe) = getOptimalCellCe(result, cellMatches, helpers, cell);
		double scoreDifference = getCellLogProb(result, cellMatches, helpers, cell, matActive ? matCe : patCe, matActive);
		scoreDifference -= getCellLogProb(result, cellMatches, helpers, cell, matActive ? patCe : matCe, !matActive);
		double escapeDifference = getCellLogProb(result, cellMatches, helpers, cell, result.cellEscapeFraction.at(cell), matActive);
		escapeDifference -= getCellLogProb(result, cellMatches, helpers, cell, 0.5, matActive);
		stream << cell << "\t" << (matActive ? "mat" : "pat") << "\t" << result.cellEscapeFraction.at(cell) << "\t" << cellCoverage.at(cell) << "\t" << scoreDifference << "\t" << escapeDifference << std::endl;
	}
}

EMHelperVariables getHelpers(const std::vector<CellMatch>& cellMatches)
{
	EMHelperVariables helpers;
	for (const auto& t : cellMatches)
	{
		helpers.activeCellsPerVariant[t.variant].insert(t.cell);
		helpers.activeVariantsPerCell[t.cell].insert(t.variant);
		helpers.variantCoverage[t.variant] += t.count;
		if (t.alt)
		{
			helpers.cellVariantAltCount[t.cell][t.variant] += t.count;
		}
		else
		{
			helpers.cellVariantRefCount[t.cell][t.variant] += t.count;
		}
	}
	size_t totalCoverage = 0;
	for (const auto& t : cellMatches)
	{
		totalCoverage += t.count;
		helpers.cellCoverageFraction[t.cell] += t.count;
	}
	for (auto& pair : helpers.cellCoverageFraction)
	{
		pair.second /= (double)totalCoverage;
	}
	return helpers;
}

EMResult getMaximumLikelihoodEM(const std::vector<CellMatch>& cellMatches, const std::unordered_map<std::string, bool>& forcedPhases, const EMHelperVariables& helpers, const size_t randomSeed)
{
	EMResult result;
	std::cerr << "initialize with random seed " << randomSeed << std::endl;
	for (const auto& t : cellMatches)
	{
		result.cellIsMatActive[t.cell] = false;
		result.variantIsMatRef[t.variant] = false;
		result.variantEscapeFraction[t.variant] = -1;
		result.cellEscapeFraction[t.cell] = -1;
	}
	initializeRandomly(result, forcedPhases, randomSeed);
	{
		std::ofstream initial { "initial.txt" };
		writeResult(result, cellMatches, helpers, initial);
	}
	size_t overlapBetweenForcedVariantsAndAllVariants = 0;
	for (const auto& pair : result.variantIsMatRef)
	{
		if (forcedPhases.count(pair.first) == 1)
		{
			overlapBetweenForcedVariantsAndAllVariants += 1;
		}
	}
	std::cerr << overlapBetweenForcedVariantsAndAllVariants << " overlap between forced variants and all variants" << std::endl;
	size_t iteration = 0;
	std::cerr << "get initial log likelihood sum" << std::endl;
	double logprob = getTotalLogProb(result, cellMatches, helpers);
	std::cerr << "initial log likelihood sum " << logprob << std::endl;
	while (true)
	{
		bool cellChanged = maximizeCellStates(result, cellMatches, helpers);
		std::cerr << "iteration " << iteration << " cell change " << (cellChanged ? "yes" : "no") << std::endl;
		if (cellChanged)
		{
			logprob = getTotalLogProb(result, cellMatches, helpers);
			std::cerr << "iteration " << iteration << " log likelihood sum " << logprob << std::endl;
		}
		bool variantChanged = maximizeVariantStates(result, forcedPhases, cellMatches, helpers);
		std::cerr << "iteration " << iteration << " variant change " << (variantChanged ? "yes" : "no") << std::endl;
		if (variantChanged)
		{
			logprob = getTotalLogProb(result, cellMatches, helpers);
			std::cerr << "iteration " << iteration << " log likelihood sum " << logprob << std::endl;
		}
		iteration += 1;
		if (!cellChanged && !variantChanged) break;
	}
	return result;
}

std::unordered_map<std::string, bool> readForcedVariantPhases(const std::string& filename)
{
	std::unordered_map<std::string, bool> result;
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
		if (origin == "mat")
		{
			result[variant] = true;
		}
		else
		{
			assert(origin == "pat");
			result[variant] = false;
		}
	}
	return result;
}

int main(int argc, char** argv)
{
	std::string matchTableFile { argv[1] };
	std::string forcedPhaseFile { argv[2] };
	size_t randomSeed = std::stoull(argv[3]);
	std::cerr << "read forced variant phases" << std::endl;
	auto forcedPhases = readForcedVariantPhases(forcedPhaseFile);
	std::cerr << forcedPhases.size() << " forced variant phases" << std::endl;
	std::cerr << "read match counts" << std::endl;
	std::vector<CellMatch> counts = readMatchCounts(matchTableFile);
	std::cerr << "get helper variables" << std::endl;
	EMHelperVariables helpers = getHelpers(counts);
	std::cerr << "run EM" << std::endl;
	EMResult result = getMaximumLikelihoodEM(counts, forcedPhases, helpers, randomSeed);
	std::cerr << "write results" << std::endl;
	writeResult(result, counts, helpers, std::cout);
}
