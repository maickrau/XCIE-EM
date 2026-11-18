#include <cmath>
#include <random>
#include <algorithm>
#include <limits>
#include "../src/EM.h"
#include "../src/Common.h"
#include "../src/FileHandler.h"
#include "InternalDeclarations.h"
#include "TestHelpers.h"
#include "doctest.h"

using TestHelpers::makeMatch;

namespace
{

const double escapeBoundary = 0.001;
const double maxEscape = 0.999;

double poissonLogPmf(size_t n, double lambda)
{
	return n * log(lambda) - lambda - std::lgamma(n + 1.0);
}

std::string variantName(size_t i)
{
	return "X:" + std::to_string(i * 100 + 1) + ":A:G";
}

std::string cellName(size_t i)
{
	return "cell" + std::to_string(i);
}

struct SimulatedData
{
	std::vector<CellMatch> matches;
	std::vector<bool> variantIsMatRef;
	std::vector<bool> cellIsMatActive;
};

// Each cell covers each variant with the given probability. Reads come from the
// active chrX, except with probability inactiveFraction from the inactive chrX.
SimulatedData simulate(size_t numVariants, size_t numCells, double coverageProbability, double inactiveFraction, size_t seed)
{
	SimulatedData result;
	std::mt19937 rng(seed);
	std::uniform_real_distribution<double> uniform(0, 1);
	std::uniform_int_distribution<size_t> readCount(1, 5);
	for (size_t i = 0; i < numVariants; i++) result.variantIsMatRef.push_back(uniform(rng) < 0.5);
	for (size_t i = 0; i < numCells; i++) result.cellIsMatActive.push_back(uniform(rng) < 0.5);
	for (size_t cell = 0; cell < numCells; cell++)
	{
		bool anyCoverage = false;
		for (size_t variant = 0; variant < numVariants; variant++)
		{
			if (uniform(rng) >= coverageProbability && !(variant + 1 == numVariants && !anyCoverage)) continue;
			anyCoverage = true;
			size_t reads = readCount(rng);
			size_t refCount = 0;
			size_t altCount = 0;
			for (size_t r = 0; r < reads; r++)
			{
				bool fromActive = uniform(rng) >= inactiveFraction;
				bool refAlleleActive = (result.variantIsMatRef[variant] == result.cellIsMatActive[cell]);
				bool isRef = (fromActive == refAlleleActive);
				if (isRef) refCount += 1; else altCount += 1;
			}
			if (refCount > 0) result.matches.push_back(makeMatch(cellName(cell), variantName(variant), false, refCount));
			if (altCount > 0) result.matches.push_back(makeMatch(cellName(cell), variantName(variant), true, altCount));
		}
	}
	return result;
}

// fraction of variants (or cells) whose phase matches the truth, allowing a global flip
double accuracy(const std::vector<bool>& truth, const std::vector<bool>& predicted, bool flip)
{
	size_t correct = 0;
	for (size_t i = 0; i < truth.size(); i++)
	{
		if ((truth[i] != predicted[i]) == flip) correct += 1;
	}
	return (double)correct / (double)truth.size();
}

std::vector<bool> resultVariantPhasesInTruthOrder(const EMOutput& output, size_t numVariants)
{
	std::vector<bool> result;
	for (size_t i = 0; i < numVariants; i++) result.push_back(output.result.variantIsMatRef[output.helpers.variantNameToIndex.at(variantName(i))]);
	return result;
}

std::vector<bool> resultCellActivitiesInTruthOrder(const EMOutput& output, size_t numCells)
{
	std::vector<bool> result;
	for (size_t i = 0; i < numCells; i++) result.push_back(output.result.cellIsMatActive[output.helpers.cellNameToIndex.at(cellName(i))]);
	return result;
}

EMResult makeResult(const EMHelperVariables& helpers, bool variantMatRef, double Xe, bool cellMatActive, double Ce)
{
	EMResult result = initializeResult(helpers);
	for (size_t i = 0; i < helpers.numVariants(); i++)
	{
		result.variantIsMatRef[i] = variantMatRef;
		result.variantEscapeFraction[i] = Xe;
	}
	for (size_t i = 0; i < helpers.numCells(); i++)
	{
		result.cellIsMatActive[i] = cellMatActive;
		result.cellEscapeFraction[i] = Ce;
	}
	return result;
}

EMCache makeCache(const EMHelperVariables& helpers)
{
	EMCache cache;
	cache.variantPatScores.resize(helpers.numVariants(), 0);
	cache.variantMatScores.resize(helpers.numVariants(), 0);
	cache.variantXeMat.resize(helpers.numVariants(), 0);
	cache.variantXePat.resize(helpers.numVariants(), 0);
	cache.hasVariant.resize(helpers.numVariants(), false);
	cache.cellPatScores.resize(helpers.numCells(), 0);
	cache.cellMatScores.resize(helpers.numCells(), 0);
	cache.cellCeMat.resize(helpers.numCells(), 0);
	cache.cellCePat.resize(helpers.numCells(), 0);
	cache.hasCell.resize(helpers.numCells(), false);
	return cache;
}

NoiseMaker zeroNoise()
{
	NoiseMaker noise;
	noise.initializeSeed(0);
	noise.magnitude = 0;
	return noise;
}

// randomized EMResult with valid escapes
EMResult randomResult(const EMHelperVariables& helpers, size_t seed)
{
	EMResult result = initializeResult(helpers);
	initializeRandomly(result, {}, seed);
	return result;
}

double bestGridCellScore(const EMResult& result, const EMHelperVariables& helpers, size_t cell, bool mat, const std::vector<bool>& ignore)
{
	double best = -std::numeric_limits<double>::infinity();
	for (size_t k = 1; k <= 999; k++)
	{
		best = std::max(best, getCellLogProb(result, helpers, cell, k * 0.001, mat, ignore));
	}
	return best;
}

double bestGridVariantScore(const EMResult& result, const EMHelperVariables& helpers, size_t variant, bool mat)
{
	double best = -std::numeric_limits<double>::infinity();
	for (size_t k = 1; k <= 999; k++)
	{
		best = std::max(best, getVariantLogProbs(result, helpers, variant, k * 0.001, mat));
	}
	return best;
}

}

TEST_SUITE("EM")
{

TEST_CASE("getHelpers indexes, coverages and adjacency")
{
	std::vector<CellMatch> matches {
		makeMatch("c1", "X:100:A:G", false, 3),
		makeMatch("c2", "X:200:A:G", true, 4),
		makeMatch("c1", "X:100:A:G", true, 2),
		makeMatch("c1", "X:200:A:G", false, 1),
		makeMatch("c1", "X:100:A:G", false, 5),
	};
	EMHelperVariables helpers = getHelpers(matches, {});
	REQUIRE(helpers.numVariants() == 2);
	REQUIRE(helpers.numCells() == 2);
	// indices are assigned in order of first appearance
	CHECK(helpers.variantNameToIndex.at("X:100:A:G") == 0);
	CHECK(helpers.variantNameToIndex.at("X:200:A:G") == 1);
	CHECK(helpers.cellNameToIndex.at("c1") == 0);
	CHECK(helpers.cellNameToIndex.at("c2") == 1);
	CHECK(helpers.variantCoverage == std::vector<size_t> { 10, 5 });
	CHECK(helpers.cellCoverage == std::vector<size_t> { 11, 4 });
	REQUIRE(helpers.cellCoverageFraction.size() == 2);
	CHECK(helpers.cellCoverageFraction[0] == doctest::Approx(11.0 / 15.0));
	CHECK(helpers.cellCoverageFraction[1] == doctest::Approx(4.0 / 15.0));
	CHECK(helpers.variantWeight == std::vector<double> { 1, 1 });
	SUBCASE("ref and alt counts of the same cell-variant pair are combined and duplicates summed")
	{
		REQUIRE(helpers.activeVariantsPerCell[0].size() == 2);
		auto sortedEdges = helpers.activeVariantsPerCell[0];
		std::sort(sortedEdges.begin(), sortedEdges.end());
		CHECK(sortedEdges[0] == std::make_tuple<size_t, size_t, size_t>(0, 8, 2));
		CHECK(sortedEdges[1] == std::make_tuple<size_t, size_t, size_t>(1, 1, 0));
		REQUIRE(helpers.activeVariantsPerCell[1].size() == 1);
		CHECK(helpers.activeVariantsPerCell[1][0] == std::make_tuple<size_t, size_t, size_t>(1, 0, 4));
		REQUIRE(helpers.activeCellsPerVariant[0].size() == 1);
		CHECK(helpers.activeCellsPerVariant[0][0] == std::make_tuple<size_t, size_t, size_t>(0, 8, 2));
		REQUIRE(helpers.activeCellsPerVariant[1].size() == 2);
	}
	SUBCASE("adjacency lists are mirror images")
	{
		size_t edges = 0;
		for (size_t cell = 0; cell < helpers.numCells(); cell++)
		{
			for (const auto& t : helpers.activeVariantsPerCell[cell])
			{
				edges += 1;
				const auto& otherSide = helpers.activeCellsPerVariant[std::get<0>(t)];
				CHECK(std::find(otherSide.begin(), otherSide.end(), std::make_tuple(cell, std::get<1>(t), std::get<2>(t))) != otherSide.end());
			}
		}
		size_t edges2 = 0;
		for (const auto& list : helpers.activeCellsPerVariant) edges2 += list.size();
		CHECK(edges == edges2);
	}
}

TEST_CASE("getHelpers region weights")
{
	std::vector<CellMatch> matches {
		makeMatch("c1", "X:99:A:G", false, 1),
		makeMatch("c1", "X:100:A:G", false, 1),
		makeMatch("c1", "X:150:A:G", false, 1),
		makeMatch("c1", "X:200:A:G", false, 1),
		makeMatch("c1", "X:300:A:G", false, 1),
	};
	SUBCASE("regions are half-open [start, end)")
	{
		EMHelperVariables helpers = getHelpers(matches, { std::make_tuple<size_t, size_t, double>(100, 200, 0.1) });
		CHECK(helpers.variantWeight[helpers.variantNameToIndex.at("X:99:A:G")] == 1);
		CHECK(helpers.variantWeight[helpers.variantNameToIndex.at("X:100:A:G")] == 0.1);
		CHECK(helpers.variantWeight[helpers.variantNameToIndex.at("X:150:A:G")] == 0.1);
		CHECK(helpers.variantWeight[helpers.variantNameToIndex.at("X:200:A:G")] == 1);
		CHECK(helpers.variantWeight[helpers.variantNameToIndex.at("X:300:A:G")] == 1);
	}
	SUBCASE("later overlapping regions override earlier ones")
	{
		EMHelperVariables helpers = getHelpers(matches, { std::make_tuple<size_t, size_t, double>(0, 1000, 0.1), std::make_tuple<size_t, size_t, double>(140, 160, -1) });
		CHECK(helpers.variantWeight[helpers.variantNameToIndex.at("X:99:A:G")] == 0.1);
		CHECK(helpers.variantWeight[helpers.variantNameToIndex.at("X:150:A:G")] == -1);
	}
}

TEST_CASE("getHelpers with no matches")
{
	EMHelperVariables helpers = getHelpers({}, {});
	CHECK(helpers.numVariants() == 0);
	CHECK(helpers.numCells() == 0);
}

TEST_CASE("initializeResult")
{
	EMHelperVariables helpers = getHelpers({ makeMatch("c1", "X:1:A:G", false, 1), makeMatch("c2", "X:2:A:G", false, 1), makeMatch("c3", "X:2:A:G", false, 1) }, {});
	EMResult result = initializeResult(helpers);
	CHECK(result.variantIsMatRef == std::vector<bool> { false, false });
	CHECK(result.variantEscapeFraction == std::vector<double> { -1, -1 });
	CHECK(result.cellIsMatActive == std::vector<bool> { false, false, false });
	CHECK(result.cellEscapeFraction == std::vector<double> { -1, -1, -1 });
}

TEST_CASE("initializeRandomly")
{
	auto data = simulate(50, 100, 0.3, 0, 1);
	EMHelperVariables helpers = getHelpers(data.matches, {});
	std::unordered_map<size_t, bool> forced { { 0, true }, { 1, false }, { 2, true } };
	EMResult result = initializeResult(helpers);
	initializeRandomly(result, forced, 42);
	SUBCASE("forced phases are applied")
	{
		CHECK(result.variantIsMatRef[0] == true);
		CHECK(result.variantIsMatRef[1] == false);
		CHECK(result.variantIsMatRef[2] == true);
	}
	SUBCASE("escapes are within bounds")
	{
		for (double e : result.variantEscapeFraction)
		{
			CHECK(e >= escapeBoundary);
			CHECK(e <= maxEscape);
		}
		for (double e : result.cellEscapeFraction)
		{
			CHECK(e >= escapeBoundary);
			CHECK(e <= maxEscape);
		}
	}
	SUBCASE("both phases occur")
	{
		size_t matVariants = std::count(result.variantIsMatRef.begin(), result.variantIsMatRef.end(), true);
		size_t matCells = std::count(result.cellIsMatActive.begin(), result.cellIsMatActive.end(), true);
		CHECK(matVariants > 5);
		CHECK(matVariants < 45);
		CHECK(matCells > 10);
		CHECK(matCells < 90);
	}
	SUBCASE("deterministic for a seed and different between seeds")
	{
		EMResult same = initializeResult(helpers);
		initializeRandomly(same, forced, 42);
		CHECK(same.variantIsMatRef == result.variantIsMatRef);
		CHECK(same.cellIsMatActive == result.cellIsMatActive);
		CHECK(same.variantEscapeFraction == result.variantEscapeFraction);
		EMResult other = initializeResult(helpers);
		initializeRandomly(other, forced, 43);
		CHECK(other.cellIsMatActive != result.cellIsMatActive);
	}
}

TEST_CASE("parseForcedPhases")
{
	EMHelperVariables helpers = getHelpers({ makeMatch("c1", "X:1:A:G", false, 1), makeMatch("c1", "X:2:A:G", false, 1) }, {});
	auto parsed = parseForcedPhases({ { "X:2:A:G", true }, { "X:1:A:G", false }, { "X:3:A:G", true }, { "1:2:A:G", true } }, helpers);
	CHECK(parsed.size() == 2);
	CHECK(parsed.at(helpers.variantNameToIndex.at("X:1:A:G")) == false);
	CHECK(parsed.at(helpers.variantNameToIndex.at("X:2:A:G")) == true);
	CHECK(parseForcedPhases({}, helpers).empty());
}

TEST_CASE("logprob is a weighted Poisson log likelihood")
{
	const double f = 0.2;
	const double c = 50;
	for (double Ce : { 0.001, 0.1, 0.5, 0.999 })
	{
		for (double Xe : { 0.001, 0.3, 0.999 })
		{
			const double E = 1.0 - (1.0 - Ce) * (1.0 - Xe);
			const double lambdaActive = f * c * (1.0 - E / 2.0);
			const double lambdaInactive = f * c * E / 2.0;
			CHECK(lambdaActive + lambdaInactive == doctest::Approx(f * c));
			for (size_t n : { 0, 1, 2, 7, 30 })
			{
				CHECK(logprob(n, f, c, Ce, Xe, 1.0, true) == doctest::Approx(poissonLogPmf(n, lambdaActive)));
				CHECK(logprob(n, f, c, Ce, Xe, 1.0, false) == doctest::Approx(poissonLogPmf(n, lambdaInactive)));
				CHECK(logprob(n, f, c, Ce, Xe, 0.1, true) == doctest::Approx(0.1 * poissonLogPmf(n, lambdaActive)));
				CHECK(logprob(n, f, c, Ce, Xe, -1, false) == doctest::Approx(-poissonLogPmf(n, lambdaInactive)));
			}
		}
	}
}

TEST_CASE("logprob with no escape prefers active reads")
{
	CHECK(logprob(5, 0.5, 20, 0.001, 0.001, 1, true) > logprob(5, 0.5, 20, 0.001, 0.001, 1, false));
	// with full escape both alleles are expressed equally
	CHECK(logprob(5, 0.5, 20, 0.999, 0.999, 1, true) == doctest::Approx(logprob(5, 0.5, 20, 0.999, 0.999, 1, false)).epsilon(0.001));
}

TEST_CASE("union find merge and find")
{
	std::vector<size_t> parent { 0, 1, 2, 3, 4, 5 };
	for (size_t i = 0; i < parent.size(); i++) CHECK(find(parent, i) == i);
	merge(parent, 0, 1);
	merge(parent, 2, 3);
	CHECK(find(parent, 0) == find(parent, 1));
	CHECK(find(parent, 2) == find(parent, 3));
	CHECK(find(parent, 0) != find(parent, 2));
	CHECK(find(parent, 4) == 4);
	merge(parent, 1, 3);
	CHECK(find(parent, 0) == find(parent, 3));
	CHECK(find(parent, 2) == find(parent, 1));
	CHECK(find(parent, 5) == 5);
	merge(parent, 3, 0);
	CHECK(find(parent, 0) == find(parent, 2));
	SUBCASE("long chains")
	{
		std::vector<size_t> chain;
		for (size_t i = 0; i < 100; i++) chain.push_back(i);
		for (size_t i = 1; i < 100; i++) merge(chain, i, i-1);
		size_t root = find(chain, 0);
		for (size_t i = 0; i < 100; i++) CHECK(find(chain, i) == root);
	}
}

TEST_CASE("getCellLogProb sums logprob differences over the cell's variants")
{
	std::vector<CellMatch> matches {
		makeMatch("c1", "X:100:A:G", false, 3),
		makeMatch("c1", "X:100:A:G", true, 1),
		makeMatch("c1", "X:200:A:G", true, 2),
		makeMatch("c2", "X:200:A:G", false, 5),
	};
	EMHelperVariables helpers = getHelpers(matches, {});
	EMResult result = initializeResult(helpers);
	const size_t v1 = helpers.variantNameToIndex.at("X:100:A:G");
	const size_t v2 = helpers.variantNameToIndex.at("X:200:A:G");
	const size_t c1 = helpers.cellNameToIndex.at("c1");
	result.variantIsMatRef[v1] = true;
	result.variantIsMatRef[v2] = false;
	result.variantEscapeFraction[v1] = 0.1;
	result.variantEscapeFraction[v2] = 0.2;
	result.cellEscapeFraction[0] = 0.3;
	result.cellEscapeFraction[1] = 0.3;
	std::vector<bool> ignoreNothing { false, false };
	const double f = helpers.cellCoverageFraction[c1];
	const double Ce = 0.05;
	auto term = [&](size_t n, double c, double Xe, bool active) { return logprob(n, f, c, Ce, Xe, 1, active) - logprob(0, f, c, Ce, Xe, 1, active); };
	// mat active: v1 ref is on mat (active), v2 ref is on pat (inactive) so v2 alt is active
	double expectedMat = term(3, 4, 0.1, true) + term(1, 4, 0.1, false) + term(2, 7, 0.2, true);
	double expectedPat = term(3, 4, 0.1, false) + term(1, 4, 0.1, true) + term(2, 7, 0.2, false);
	CHECK(getCellLogProb(result, helpers, c1, Ce, true, ignoreNothing) == doctest::Approx(expectedMat));
	CHECK(getCellLogProb(result, helpers, c1, Ce, false, ignoreNothing) == doctest::Approx(expectedPat));
	CHECK(expectedMat > expectedPat);
	SUBCASE("ignored variants are skipped")
	{
		std::vector<bool> ignoreV2 { false, false };
		ignoreV2[v2] = true;
		CHECK(getCellLogProb(result, helpers, c1, Ce, true, ignoreV2) == doctest::Approx(term(3, 4, 0.1, true) + term(1, 4, 0.1, false)));
		std::vector<bool> ignoreAll { true, true };
		CHECK(getCellLogProb(result, helpers, c1, Ce, true, ignoreAll) == 0);
	}
	SUBCASE("variant weights scale the contribution")
	{
		EMHelperVariables weighted = getHelpers(matches, { std::make_tuple<size_t, size_t, double>(150, 250, 0.1) });
		double expectedWeighted = term(3, 4, 0.1, true) + term(1, 4, 0.1, false) + 0.1 * term(2, 7, 0.2, true);
		CHECK(getCellLogProb(result, weighted, c1, Ce, true, ignoreNothing) == doctest::Approx(expectedWeighted));
	}
}

TEST_CASE("getVariantLogProbs sums logprob differences over the variant's cells")
{
	std::vector<CellMatch> matches {
		makeMatch("c1", "X:100:A:G", false, 3),
		makeMatch("c1", "X:100:A:G", true, 1),
		makeMatch("c2", "X:100:A:G", true, 4),
		makeMatch("c2", "X:200:A:G", false, 2),
	};
	EMHelperVariables helpers = getHelpers(matches, {});
	EMResult result = initializeResult(helpers);
	const size_t v1 = helpers.variantNameToIndex.at("X:100:A:G");
	const size_t c1 = helpers.cellNameToIndex.at("c1");
	const size_t c2 = helpers.cellNameToIndex.at("c2");
	result.cellIsMatActive[c1] = true;
	result.cellIsMatActive[c2] = false;
	result.cellEscapeFraction[c1] = 0.1;
	result.cellEscapeFraction[c2] = 0.2;
	result.variantEscapeFraction[0] = 0.5;
	result.variantEscapeFraction[1] = 0.5;
	const double Xe = 0.05;
	auto term = [&](size_t n, size_t cell, bool active) { double f = helpers.cellCoverageFraction[cell]; double Ce = result.cellEscapeFraction[cell]; return logprob(n, f, 8, Ce, Xe, 1, active) - logprob(0, f, 8, Ce, Xe, 1, active); };
	double expectedMat = term(3, c1, true) + term(1, c1, false) + term(4, c2, true);
	double expectedPat = term(3, c1, false) + term(1, c1, true) + term(4, c2, false);
	CHECK(getVariantLogProbs(result, helpers, v1, Xe, true) == doctest::Approx(expectedMat));
	CHECK(getVariantLogProbs(result, helpers, v1, Xe, false) == doctest::Approx(expectedPat));
	CHECK(expectedMat > expectedPat);
}

TEST_CASE("analytic derivatives match numeric derivatives")
{
	auto data = simulate(20, 40, 0.4, 0.2, 7);
	EMHelperVariables helpers = getHelpers(data.matches, { std::make_tuple<size_t, size_t, double>(0, 500, 0.1) });
	EMResult result = randomResult(helpers, 3);
	std::vector<bool> ignoreNothing(helpers.numVariants(), false);
	std::vector<bool> ignoreSome(helpers.numVariants(), false);
	for (size_t i = 0; i < ignoreSome.size(); i += 3) ignoreSome[i] = true;
	const double h = 1e-6;
	for (double x : { 0.01, 0.2, 0.5, 0.8, 0.95 })
	{
		for (bool mat : { true, false })
		{
			for (size_t cell = 0; cell < helpers.numCells(); cell++)
			{
				for (const auto* ignore : { &ignoreNothing, &ignoreSome })
				{
					double numeric = (getCellLogProb(result, helpers, cell, x + h, mat, *ignore) - getCellLogProb(result, helpers, cell, x - h, mat, *ignore)) / (2 * h);
					double analytic = getCellLogProbDerivative(result, helpers, cell, x, mat, *ignore);
					CHECK(analytic == doctest::Approx(numeric).epsilon(1e-4).scale(1));
				}
			}
			for (size_t variant = 0; variant < helpers.numVariants(); variant++)
			{
				double numeric = (getVariantLogProbs(result, helpers, variant, x + h, mat) - getVariantLogProbs(result, helpers, variant, x - h, mat)) / (2 * h);
				double analytic = getVariantLogProbDerivative(result, helpers, variant, x, mat);
				CHECK(analytic == doctest::Approx(numeric).epsilon(1e-4).scale(1));
			}
		}
	}
}

TEST_CASE("cell and variant log likelihoods are concave in the escape fraction")
{
	auto data = simulate(10, 30, 0.5, 0.3, 11);
	EMHelperVariables helpers = getHelpers(data.matches, {});
	EMResult result = randomResult(helpers, 5);
	std::vector<bool> ignoreNothing(helpers.numVariants(), false);
	for (bool mat : { true, false })
	{
		for (size_t cell = 0; cell < helpers.numCells(); cell++)
		{
			double previous = std::numeric_limits<double>::infinity();
			for (size_t k = 1; k < 100; k++)
			{
				double derivative = getCellLogProbDerivative(result, helpers, cell, k * 0.01, mat, ignoreNothing);
				CHECK(derivative <= previous + 1e-9);
				previous = derivative;
			}
		}
		for (size_t variant = 0; variant < helpers.numVariants(); variant++)
		{
			double previous = std::numeric_limits<double>::infinity();
			for (size_t k = 1; k < 100; k++)
			{
				double derivative = getVariantLogProbDerivative(result, helpers, variant, k * 0.01, mat);
				CHECK(derivative <= previous + 1e-9);
				previous = derivative;
			}
		}
	}
}

TEST_CASE("binarySearchOptimalCe analytic optimum")
{
	// One cell, one variant, Xe at its minimum. With a active and b inactive reads
	// the likelihood a*log(1-E/2) + b*log(E/2) is maximized at E = 2b/(a+b),
	// and E = 1 - (1-Ce)(1-Xe) gives Ce = 1 - (1-E)/(1-Xe).
	auto check = [](size_t activeReads, size_t inactiveReads, double expectedCe)
	{
		std::vector<CellMatch> matches;
		if (activeReads > 0) matches.push_back(makeMatch("c", "X:1:A:G", false, activeReads));
		if (inactiveReads > 0) matches.push_back(makeMatch("c", "X:1:A:G", true, inactiveReads));
		EMHelperVariables helpers = getHelpers(matches, {});
		EMResult result = makeResult(helpers, true, escapeBoundary, true, 0.5);
		std::vector<bool> ignoreNothing { false };
		double Ce = binarySearchOptimalCe(result, helpers, 0, true, ignoreNothing);
		INFO("active ", activeReads, " inactive ", inactiveReads);
		CHECK(Ce == doctest::Approx(expectedCe).epsilon(0.0015).scale(1));
		CHECK(Ce >= escapeBoundary);
		CHECK(Ce <= maxEscape);
	};
	check(3, 1, 1.0 - 0.5 / (1.0 - escapeBoundary));
	check(9, 1, 1.0 - 0.8 / (1.0 - escapeBoundary));
	check(19, 1, 1.0 - 0.9 / (1.0 - escapeBoundary));
	// balanced expression: E=1 so Ce is at the upper boundary
	check(5, 5, maxEscape);
	// more inactive than active: optimum outside the range, clamped to the upper boundary
	check(1, 3, maxEscape);
	check(0, 3, maxEscape);
	// only active reads: Ce at the lower boundary
	check(10, 0, escapeBoundary);
}

TEST_CASE("binarySearchOptimalXe analytic optimum")
{
	auto check = [](size_t activeReads, size_t inactiveReads, double expectedXe)
	{
		std::vector<CellMatch> matches;
		if (activeReads > 0) matches.push_back(makeMatch("c", "X:1:A:G", false, activeReads));
		if (inactiveReads > 0) matches.push_back(makeMatch("c", "X:1:A:G", true, inactiveReads));
		EMHelperVariables helpers = getHelpers(matches, {});
		EMResult result = makeResult(helpers, true, 0.5, true, escapeBoundary);
		double Xe = binarySearchOptimalXe(result, helpers, 0, true);
		INFO("active ", activeReads, " inactive ", inactiveReads);
		CHECK(Xe == doctest::Approx(expectedXe).epsilon(0.0015).scale(1));
	};
	check(3, 1, 1.0 - 0.5 / (1.0 - escapeBoundary));
	check(9, 1, 1.0 - 0.8 / (1.0 - escapeBoundary));
	check(5, 5, maxEscape);
	check(10, 0, escapeBoundary);
	check(0, 10, maxEscape);
}

TEST_CASE("getOptimalCellCe and getOptimalVariantXe return mat and pat optima")
{
	std::vector<CellMatch> matches { makeMatch("c", "X:1:A:G", false, 3), makeMatch("c", "X:1:A:G", true, 1) };
	EMHelperVariables helpers = getHelpers(matches, {});
	EMResult result = makeResult(helpers, true, escapeBoundary, true, escapeBoundary);
	std::vector<bool> ignoreNothing { false };
	auto Ce = getOptimalCellCe(result, helpers, 0, ignoreNothing);
	CHECK(Ce.first == binarySearchOptimalCe(result, helpers, 0, true, ignoreNothing));
	CHECK(Ce.second == binarySearchOptimalCe(result, helpers, 0, false, ignoreNothing));
	// mat active: 3 active 1 inactive. pat active: 1 active 3 inactive -> upper boundary
	CHECK(Ce.first == doctest::Approx(0.4995).epsilon(0.0015).scale(1));
	CHECK(Ce.second == doctest::Approx(maxEscape));
	auto Xe = getOptimalVariantXe(result, helpers, 0);
	CHECK(Xe.first == binarySearchOptimalXe(result, helpers, 0, true));
	CHECK(Xe.second == binarySearchOptimalXe(result, helpers, 0, false));
	CHECK(Xe.second == doctest::Approx(maxEscape));
}

TEST_CASE("binary searched escapes are the best values on the 0.001 grid")
{
	for (size_t seed = 0; seed < 3; seed++)
	{
		auto data = simulate(15, 30, 0.4, 0.15, 100 + seed);
		EMHelperVariables helpers = getHelpers(data.matches, { std::make_tuple<size_t, size_t, double>(0, 300, 0.1) });
		EMResult result = randomResult(helpers, seed);
		std::vector<bool> ignoreNothing(helpers.numVariants(), false);
		for (size_t cell = 0; cell < helpers.numCells(); cell++)
		{
			for (bool mat : { true, false })
			{
				double Ce = binarySearchOptimalCe(result, helpers, cell, mat, ignoreNothing);
				CHECK(Ce >= escapeBoundary);
				CHECK(Ce <= maxEscape);
				double score = getCellLogProb(result, helpers, cell, Ce, mat, ignoreNothing);
				CHECK(score == doctest::Approx(bestGridCellScore(result, helpers, cell, mat, ignoreNothing)).epsilon(1e-9).scale(1));
			}
		}
		for (size_t variant = 0; variant < helpers.numVariants(); variant++)
		{
			for (bool mat : { true, false })
			{
				double Xe = binarySearchOptimalXe(result, helpers, variant, mat);
				CHECK(Xe >= escapeBoundary);
				CHECK(Xe <= maxEscape);
				double score = getVariantLogProbs(result, helpers, variant, Xe, mat);
				CHECK(score == doctest::Approx(bestGridVariantScore(result, helpers, variant, mat)).epsilon(1e-9).scale(1));
			}
		}
	}
}

TEST_CASE("binary searched escapes with a negatively weighted (XIST) variant" * doctest::should_fail())
{
	// BUG (minor): XIST regions get weight -1, which makes that variant's log likelihood term
	// convex instead of concave. binarySearchOptimalCe / binarySearchOptimalXe assume a concave
	// objective (monotonically decreasing derivative) and can then return an escape value
	// that is clearly worse than the best value on the grid.
	size_t worse = 0;
	size_t total = 0;
	double worstDifference = 0;
	for (size_t seed = 0; seed < 5; seed++)
	{
		auto data = simulate(15, 30, 0.5, 0.15, 200 + seed);
		EMHelperVariables helpers = getHelpers(data.matches, { std::make_tuple<size_t, size_t, double>(0, 200, -1) });
		EMResult result = randomResult(helpers, seed);
		std::vector<bool> ignoreNothing(helpers.numVariants(), false);
		for (size_t cell = 0; cell < helpers.numCells(); cell++)
		{
			for (bool mat : { true, false })
			{
				double Ce = binarySearchOptimalCe(result, helpers, cell, mat, ignoreNothing);
				double score = getCellLogProb(result, helpers, cell, Ce, mat, ignoreNothing);
				double best = bestGridCellScore(result, helpers, cell, mat, ignoreNothing);
				total += 1;
				if (score < best - 1e-6)
				{
					worse += 1;
					worstDifference = std::max(worstDifference, best - score);
				}
			}
		}
	}
	INFO(worse, " of ", total, " cell escape searches were worse than the grid optimum, by up to ", worstDifference);
	CHECK(worse == 0);
	size_t worseVariants = 0;
	size_t totalVariants = 0;
	double worstVariantDifference = 0;
	for (size_t seed = 0; seed < 5; seed++)
	{
		auto data = simulate(15, 30, 0.5, 0.15, 200 + seed);
		EMHelperVariables helpers = getHelpers(data.matches, { std::make_tuple<size_t, size_t, double>(0, 200, -1) });
		EMResult result = randomResult(helpers, seed);
		for (size_t variant = 0; variant < helpers.numVariants(); variant++)
		{
			for (bool mat : { true, false })
			{
				double Xe = binarySearchOptimalXe(result, helpers, variant, mat);
				double score = getVariantLogProbs(result, helpers, variant, Xe, mat);
				double best = bestGridVariantScore(result, helpers, variant, mat);
				totalVariants += 1;
				if (score < best - 1e-6)
				{
					worseVariants += 1;
					worstVariantDifference = std::max(worstVariantDifference, best - score);
				}
			}
		}
	}
	INFO(worseVariants, " of ", totalVariants, " variant escape searches were worse than the grid optimum, by up to ", worstVariantDifference);
	CHECK(worseVariants == 0);
}

TEST_CASE("total log probability equals the sum over variants and the sum over cells")
{
	auto data = simulate(20, 50, 0.3, 0.2, 21);
	EMHelperVariables helpers = getHelpers(data.matches, { std::make_tuple<size_t, size_t, double>(0, 700, 0.1) });
	EMResult result = randomResult(helpers, 9);
	std::vector<bool> ignoreNothing(helpers.numVariants(), false);
	double total = getTotalLogProb(result, helpers);
	CHECK(total == doctest::Approx(getTotalLogProb(result, helpers, ignoreNothing)));
	CHECK(total == doctest::Approx(getNonnormalizedTotalLogProb(result, helpers, ignoreNothing)));
	double variantSum = 0;
	for (size_t v = 0; v < helpers.numVariants(); v++) variantSum += getVariantLogProbs(result, helpers, v, result.variantEscapeFraction[v], result.variantIsMatRef[v]);
	double cellSum = 0;
	for (size_t c = 0; c < helpers.numCells(); c++) cellSum += getCellLogProb(result, helpers, c, result.cellEscapeFraction[c], result.cellIsMatActive[c], ignoreNothing);
	CHECK(total == doctest::Approx(variantSum));
	CHECK(total == doctest::Approx(cellSum));
	SUBCASE("ignoring variants")
	{
		std::vector<bool> ignore(helpers.numVariants(), false);
		ignore[0] = true;
		ignore[5] = true;
		double expected = total - getVariantLogProbs(result, helpers, 0, result.variantEscapeFraction[0], result.variantIsMatRef[0]) - getVariantLogProbs(result, helpers, 5, result.variantEscapeFraction[5], result.variantIsMatRef[5]);
		CHECK(getTotalLogProb(result, helpers, ignore) == doctest::Approx(expected));
	}
	SUBCASE("flipping every variant and cell keeps the total")
	{
		EMResult flipped = result;
		for (size_t i = 0; i < flipped.variantIsMatRef.size(); i++) flipped.variantIsMatRef[i] = !flipped.variantIsMatRef[i];
		for (size_t i = 0; i < flipped.cellIsMatActive.size(); i++) flipped.cellIsMatActive[i] = !flipped.cellIsMatActive[i];
		CHECK(getTotalLogProb(flipped, helpers) == doctest::Approx(total));
	}
	SUBCASE("the true phasing scores better than a random one")
	{
		EMResult truth = result;
		for (size_t v = 0; v < data.variantIsMatRef.size(); v++) truth.variantIsMatRef[helpers.variantNameToIndex.at(variantName(v))] = data.variantIsMatRef[v];
		for (size_t c = 0; c < data.cellIsMatActive.size(); c++) truth.cellIsMatActive[helpers.cellNameToIndex.at(cellName(c))] = data.cellIsMatActive[c];
		for (auto& e : truth.variantEscapeFraction) e = 0.01;
		for (auto& e : truth.cellEscapeFraction) e = 0.3;
		EMResult random = truth;
		initializeRandomly(random, {}, 77);
		for (auto& e : random.variantEscapeFraction) e = 0.01;
		for (auto& e : random.cellEscapeFraction) e = 0.3;
		CHECK(getTotalLogProb(truth, helpers) > getTotalLogProb(random, helpers));
	}
}

TEST_CASE("escape confidence intervals")
{
	auto check = [](size_t activeReads, size_t inactiveReads)
	{
		std::vector<CellMatch> matches;
		if (activeReads > 0) matches.push_back(makeMatch("c", "X:1:A:G", false, activeReads));
		if (inactiveReads > 0) matches.push_back(makeMatch("c", "X:1:A:G", true, inactiveReads));
		EMHelperVariables helpers = getHelpers(matches, {});
		EMResult result = makeResult(helpers, true, escapeBoundary, true, escapeBoundary);
		std::vector<bool> ignoreNothing { false };
		double Ce = binarySearchOptimalCe(result, helpers, 0, true, ignoreNothing);
		auto cellCI = getCellEscapeConfidenceInterval(result, helpers, 0, Ce, true, 0.95, ignoreNothing);
		double Xe = binarySearchOptimalXe(result, helpers, 0, true);
		auto variantCI = getVariantEscapeConfidenceInterval(result, helpers, 0, Xe, true, 0.95);
		return std::make_tuple(Ce, cellCI, Xe, variantCI);
	};
	SUBCASE("interval contains the point estimate and is within bounds")
	{
		for (auto counts : { std::make_pair(3, 1), std::make_pair(30, 10), std::make_pair(10, 0), std::make_pair(5, 5), std::make_pair(100, 3) })
		{
			auto [Ce, cellCI, Xe, variantCI] = check(counts.first, counts.second);
			INFO("counts ", counts.first, " ", counts.second);
			CHECK(cellCI.first <= Ce + 1e-9);
			CHECK(cellCI.second >= Ce - 1e-9);
			CHECK(cellCI.first >= escapeBoundary - 1e-9);
			CHECK(cellCI.second <= maxEscape + 1e-9);
			CHECK(variantCI.first <= Xe + 1e-9);
			CHECK(variantCI.second >= Xe - 1e-9);
			CHECK(variantCI.first >= escapeBoundary - 1e-9);
			CHECK(variantCI.second <= maxEscape + 1e-9);
		}
	}
	SUBCASE("more data gives a narrower interval")
	{
		auto [Ce1, ci1, Xe1, vci1] = check(3, 1);
		auto [Ce2, ci2, Xe2, vci2] = check(300, 100);
		CHECK(Ce1 == doctest::Approx(Ce2).epsilon(0.01));
		CHECK(ci2.second - ci2.first < ci1.second - ci1.first);
		CHECK(vci2.second - vci2.first < vci1.second - vci1.first);
		CHECK(ci2.first > escapeBoundary);
		CHECK(ci2.second < maxEscape);
	}
	SUBCASE("interval endpoints are where the log likelihood has dropped by log(0.025)")
	{
		std::vector<CellMatch> matches { makeMatch("c", "X:1:A:G", false, 300), makeMatch("c", "X:1:A:G", true, 100) };
		EMHelperVariables helpers = getHelpers(matches, {});
		EMResult result = makeResult(helpers, true, escapeBoundary, true, escapeBoundary);
		std::vector<bool> ignoreNothing { false };
		double Ce = binarySearchOptimalCe(result, helpers, 0, true, ignoreNothing);
		auto ci = getCellEscapeConfidenceInterval(result, helpers, 0, Ce, true, 0.95, ignoreNothing);
		double best = getCellLogProb(result, helpers, 0, Ce, true, ignoreNothing);
		double threshold = best + log(0.025);
		CHECK(getCellLogProb(result, helpers, 0, ci.first, true, ignoreNothing) <= threshold);
		CHECK(getCellLogProb(result, helpers, 0, ci.first + 0.001, true, ignoreNothing) > threshold);
		CHECK(getCellLogProb(result, helpers, 0, ci.second, true, ignoreNothing) <= threshold);
		CHECK(getCellLogProb(result, helpers, 0, ci.second - 0.001, true, ignoreNothing) > threshold);
	}
}

TEST_CASE("maximizeVariantStates picks the best phase for each variant")
{
	auto data = simulate(20, 60, 0.4, 0.05, 31);
	EMHelperVariables helpers = getHelpers(data.matches, {});
	EMResult result = randomResult(helpers, 4);
	// set cells to the truth so that the variant phases are determined by the data
	for (size_t c = 0; c < data.cellIsMatActive.size(); c++) result.cellIsMatActive[helpers.cellNameToIndex.at(cellName(c))] = data.cellIsMatActive[c];
	for (auto& e : result.cellEscapeFraction) e = escapeBoundary;
	EMCache cache = makeCache(helpers);
	NoiseMaker noise = zeroNoise();
	std::vector<bool> ignoreNothing(helpers.numVariants(), false);
	EMResult before = result;
	bool changed = maximizeVariantStates(result, cache, {}, helpers, ignoreNothing, noise);
	CHECK(changed);
	for (size_t v = 0; v < helpers.numVariants(); v++)
	{
		CHECK(cache.hasVariant[v]);
		auto Xe = getOptimalVariantXe(before, helpers, v);
		double matScore = getVariantLogProbs(before, helpers, v, Xe.first, true);
		double patScore = getVariantLogProbs(before, helpers, v, Xe.second, false);
		if (matScore > patScore + 0.0001) CHECK(result.variantIsMatRef[v]);
		if (patScore > matScore + 0.0001) CHECK(!result.variantIsMatRef[v]);
		CHECK(result.variantEscapeFraction[v] == (result.variantIsMatRef[v] ? Xe.first : Xe.second));
	}
	CHECK(accuracy(data.variantIsMatRef, resultVariantPhasesInTruthOrder(EMOutput { helpers, result, {} }, data.variantIsMatRef.size()), false) == 1.0);
	SUBCASE("cells of changed variants are invalidated")
	{
		for (size_t c = 0; c < helpers.numCells(); c++) CHECK(!cache.hasCell[c]);
	}
	SUBCASE("second call with valid cache changes nothing")
	{
		EMResult after = result;
		CHECK(!maximizeVariantStates(result, cache, {}, helpers, ignoreNothing, noise));
		CHECK(result.variantIsMatRef == after.variantIsMatRef);
		CHECK(result.variantEscapeFraction == after.variantEscapeFraction);
	}
	SUBCASE("total log probability does not decrease")
	{
		CHECK(getTotalLogProb(result, helpers) >= getTotalLogProb(before, helpers) - 1e-6);
	}
}

TEST_CASE("maximizeVariantStates respects forced phases")
{
	auto data = simulate(10, 40, 0.5, 0, 41);
	EMHelperVariables helpers = getHelpers(data.matches, {});
	EMResult result = randomResult(helpers, 4);
	for (size_t c = 0; c < data.cellIsMatActive.size(); c++) result.cellIsMatActive[helpers.cellNameToIndex.at(cellName(c))] = data.cellIsMatActive[c];
	// force variant 0 to the wrong phase
	const size_t forcedVariant = helpers.variantNameToIndex.at(variantName(0));
	std::unordered_map<size_t, bool> forced { { forcedVariant, !data.variantIsMatRef[0] } };
	result.variantIsMatRef[forcedVariant] = !data.variantIsMatRef[0];
	EMCache cache = makeCache(helpers);
	NoiseMaker noise = zeroNoise();
	std::vector<bool> ignoreNothing(helpers.numVariants(), false);
	maximizeVariantStates(result, cache, forced, helpers, ignoreNothing, noise);
	CHECK(result.variantIsMatRef[forcedVariant] == !data.variantIsMatRef[0]);
	for (size_t v = 1; v < data.variantIsMatRef.size(); v++)
	{
		CHECK(result.variantIsMatRef[helpers.variantNameToIndex.at(variantName(v))] == data.variantIsMatRef[v]);
	}
}

TEST_CASE("maximizeVariantStates skips ignored variants")
{
	auto data = simulate(10, 40, 0.5, 0, 43);
	EMHelperVariables helpers = getHelpers(data.matches, {});
	EMResult result = randomResult(helpers, 4);
	EMResult before = result;
	EMCache cache = makeCache(helpers);
	NoiseMaker noise = zeroNoise();
	std::vector<bool> ignoreAll(helpers.numVariants(), true);
	CHECK(!maximizeVariantStates(result, cache, {}, helpers, ignoreAll, noise));
	CHECK(result.variantIsMatRef == before.variantIsMatRef);
	CHECK(result.variantEscapeFraction == before.variantEscapeFraction);
}

TEST_CASE("maximizeCellStates picks the best activity for each cell")
{
	auto data = simulate(20, 60, 0.4, 0.05, 51);
	EMHelperVariables helpers = getHelpers(data.matches, {});
	EMResult result = randomResult(helpers, 4);
	for (size_t v = 0; v < data.variantIsMatRef.size(); v++) result.variantIsMatRef[helpers.variantNameToIndex.at(variantName(v))] = data.variantIsMatRef[v];
	for (auto& e : result.variantEscapeFraction) e = escapeBoundary;
	EMCache cache = makeCache(helpers);
	NoiseMaker noise = zeroNoise();
	std::vector<bool> ignoreNothing(helpers.numVariants(), false);
	EMResult before = result;
	bool changed = maximizeCellStates(result, cache, helpers, ignoreNothing, noise);
	CHECK(changed);
	for (size_t c = 0; c < helpers.numCells(); c++)
	{
		CHECK(cache.hasCell[c]);
		auto Ce = getOptimalCellCe(before, helpers, c, ignoreNothing);
		CHECK(cache.cellCeMat[c] == Ce.first);
		CHECK(cache.cellCePat[c] == Ce.second);
		CHECK(result.cellEscapeFraction[c] == (result.cellIsMatActive[c] ? Ce.first : Ce.second));
	}
	CHECK(accuracy(data.cellIsMatActive, resultCellActivitiesInTruthOrder(EMOutput { helpers, result, {} }, data.cellIsMatActive.size()), false) == 1.0);
	CHECK(getTotalLogProb(result, helpers) >= getTotalLogProb(before, helpers) - 1e-6);
	SUBCASE("variants of changed cells are invalidated and a second call changes nothing")
	{
		for (size_t v = 0; v < helpers.numVariants(); v++) CHECK(!cache.hasVariant[v]);
		CHECK(!maximizeCellStates(result, cache, helpers, ignoreNothing, noise));
	}
}

TEST_CASE("alternating maximization converges and never decreases the log probability")
{
	auto data = simulate(30, 80, 0.3, 0.1, 61);
	EMHelperVariables helpers = getHelpers(data.matches, {});
	EMResult result = randomResult(helpers, 8);
	EMCache cache = makeCache(helpers);
	NoiseMaker noise = zeroNoise();
	std::vector<bool> ignoreNothing(helpers.numVariants(), false);
	maximizeVariantStates(result, cache, {}, helpers, ignoreNothing, noise);
	double previous = getTotalLogProb(result, helpers);
	size_t iterations = 0;
	while (iterations < 1000)
	{
		bool changed = maximizeCellStates(result, cache, helpers, ignoreNothing, noise);
		double afterCells = getTotalLogProb(result, helpers);
		CHECK(afterCells >= previous - 1e-6);
		changed = maximizeVariantStates(result, cache, {}, helpers, ignoreNothing, noise) || changed;
		double afterVariants = getTotalLogProb(result, helpers);
		CHECK(afterVariants >= afterCells - 1e-6);
		previous = afterVariants;
		iterations += 1;
		if (!changed) break;
	}
	CHECK(iterations < 1000);
}

TEST_CASE("getMaximumLikelihoodEM reaches a local optimum")
{
	auto data = simulate(30, 100, 0.3, 0.05, 71);
	EMHelperVariables helpers = getHelpers(data.matches, {});
	EMResult result = initializeResult(helpers);
	initializeRandomly(result, {}, 1);
	double before = getTotalLogProb(result, helpers);
	getMaximumLikelihoodEM(result, data.matches, {}, helpers, 1, 20, 0.95);
	CHECK(getTotalLogProb(result, helpers) > before);
	// no single coordinate update can improve the result
	EMCache cache = makeCache(helpers);
	NoiseMaker noise = zeroNoise();
	std::vector<bool> ignoreNothing(helpers.numVariants(), false);
	CHECK(!maximizeVariantStates(result, cache, {}, helpers, ignoreNothing, noise));
	CHECK(!maximizeCellStates(result, cache, helpers, ignoreNothing, noise));
	SUBCASE("all escapes are within bounds")
	{
		for (double e : result.variantEscapeFraction)
		{
			CHECK(e >= escapeBoundary);
			CHECK(e <= maxEscape);
		}
		for (double e : result.cellEscapeFraction)
		{
			CHECK(e >= escapeBoundary);
			CHECK(e <= maxEscape);
		}
	}
}

TEST_CASE("getMaximumLikelihoodEM without noise")
{
	auto data = simulate(20, 60, 0.4, 0, 73);
	EMHelperVariables helpers = getHelpers(data.matches, {});
	EMResult result = initializeResult(helpers);
	initializeRandomly(result, {}, 1);
	getMaximumLikelihoodEM(result, data.matches, {}, helpers, 1, 0, 0.95);
	EMCache cache = makeCache(helpers);
	NoiseMaker noise = zeroNoise();
	std::vector<bool> ignoreNothing(helpers.numVariants(), false);
	CHECK(!maximizeVariantStates(result, cache, {}, helpers, ignoreNothing, noise));
	CHECK(!maximizeCellStates(result, cache, helpers, ignoreNothing, noise));
}

TEST_CASE("flipUnsatisfiedSubgraph")
{
	// Cells c1, c2 cover v1, v2 with only active expression. Truth: v1 mat, v2 pat; c1 mat, c2 pat.
	std::vector<CellMatch> matches {
		makeMatch("c1", "X:1:A:G", false, 3),
		makeMatch("c1", "X:2:A:G", true, 3),
		makeMatch("c2", "X:1:A:G", true, 3),
		makeMatch("c2", "X:2:A:G", false, 3),
	};
	EMHelperVariables helpers = getHelpers(matches, {});
	const size_t v1 = helpers.variantNameToIndex.at("X:1:A:G");
	const size_t v2 = helpers.variantNameToIndex.at("X:2:A:G");
	const size_t c1 = helpers.cellNameToIndex.at("c1");
	const size_t c2 = helpers.cellNameToIndex.at("c2");
	EMResult result = makeResult(helpers, true, escapeBoundary, true, escapeBoundary);
	EMCache cache = makeCache(helpers);
	for (size_t i = 0; i < cache.hasVariant.size(); i++) cache.hasVariant[i] = true;
	for (size_t i = 0; i < cache.hasCell.size(); i++) cache.hasCell[i] = true;
	SUBCASE("consistent phasing is not changed")
	{
		result.variantIsMatRef[v1] = true;
		result.variantIsMatRef[v2] = false;
		result.cellIsMatActive[c1] = true;
		result.cellIsMatActive[c2] = false;
		EMResult before = result;
		CHECK(!flipUnsatisfiedSubgraph(result, cache, helpers, {}));
		CHECK(result.variantIsMatRef == before.variantIsMatRef);
		CHECK(result.cellIsMatActive == before.cellIsMatActive);
	}
	SUBCASE("a cell with no satisfied edges is flipped")
	{
		result.variantIsMatRef[v1] = true;
		result.variantIsMatRef[v2] = false;
		result.cellIsMatActive[c1] = false;
		result.cellIsMatActive[c2] = false;
		CHECK(flipUnsatisfiedSubgraph(result, cache, helpers, {}));
		CHECK(result.cellIsMatActive[c1] == true);
		CHECK(!cache.hasCell[c1]);
		CHECK(!cache.hasVariant[v1]);
		CHECK(!cache.hasVariant[v2]);
		CHECK(!flipUnsatisfiedSubgraph(result, cache, helpers, {}));
	}
	SUBCASE("a conflicting cluster is flipped")
	{
		// c3 links v2 and v3. c4 covers v3 and v4. Phase v3/v4/c4 inconsistently with the rest.
		std::vector<CellMatch> more = matches;
		more.push_back(makeMatch("c3", "X:2:A:G", false, 3));
		more.push_back(makeMatch("c3", "X:3:A:G", false, 3));
		more.push_back(makeMatch("c4", "X:3:A:G", false, 3));
		more.push_back(makeMatch("c4", "X:4:A:G", false, 3));
		EMHelperVariables helpers2 = getHelpers(more, {});
		EMResult result2 = makeResult(helpers2, true, escapeBoundary, true, escapeBoundary);
		auto vi = [&](const std::string& name) { return helpers2.variantNameToIndex.at(name); };
		auto ci = [&](const std::string& name) { return helpers2.cellNameToIndex.at(name); };
		// consistent: v1 mat, v2 pat, c1 mat, c2 pat, c3 pat, v3 pat
		result2.variantIsMatRef[vi("X:1:A:G")] = true;
		result2.variantIsMatRef[vi("X:2:A:G")] = false;
		result2.cellIsMatActive[ci("c1")] = true;
		result2.cellIsMatActive[ci("c2")] = false;
		result2.cellIsMatActive[ci("c3")] = false;
		result2.variantIsMatRef[vi("X:3:A:G")] = true;
		result2.cellIsMatActive[ci("c4")] = true;
		result2.variantIsMatRef[vi("X:4:A:G")] = true;
		EMCache cache2 = makeCache(helpers2);
		size_t flips = 0;
		while (flipUnsatisfiedSubgraph(result2, cache2, helpers2, {}) && flips < 10) flips += 1;
		CHECK(flips >= 1);
		CHECK(flips < 10);
		// all edges are now satisfied: every observed allele is on the active chrX
		for (size_t cell = 0; cell < helpers2.numCells(); cell++)
		{
			for (const auto& t : helpers2.activeVariantsPerCell[cell])
			{
				bool refActive = result2.variantIsMatRef[std::get<0>(t)] == result2.cellIsMatActive[cell];
				CHECK((refActive ? std::get<2>(t) : std::get<1>(t)) == 0);
			}
		}
	}
}

TEST_CASE("getEMAdditions")
{
	auto data = simulate(20, 60, 0.4, 0.05, 81);
	EMHelperVariables helpers = getHelpers(data.matches, {});
	EMResult result = initializeResult(helpers);
	initializeRandomly(result, {}, 1);
	getMaximumLikelihoodEM(result, data.matches, {}, helpers, 1, 20, 0.95);
	EMResultAdditions additions = getEMAdditions(data.matches, result, helpers);
	REQUIRE(additions.variantPhaseConfidence.size() == helpers.numVariants());
	REQUIRE(additions.variantEscapeCILow.size() == helpers.numVariants());
	REQUIRE(additions.variantEscapeCIHigh.size() == helpers.numVariants());
	REQUIRE(additions.cellActiveConfidence.size() == helpers.numCells());
	REQUIRE(additions.cellActiveConfidenceOnlyNonescape.size() == helpers.numCells());
	REQUIRE(additions.cellEscapeCILow.size() == helpers.numCells());
	REQUIRE(additions.cellEscapeCIHigh.size() == helpers.numCells());
	std::vector<bool> ignoreNothing(helpers.numVariants(), false);
	SUBCASE("confidences are log likelihood differences to the flipped state")
	{
		for (size_t v = 0; v < helpers.numVariants(); v++)
		{
			auto Xe = getOptimalVariantXe(result, helpers, v);
			double mat = getVariantLogProbs(result, helpers, v, Xe.first, true);
			double pat = getVariantLogProbs(result, helpers, v, Xe.second, false);
			CHECK(additions.variantPhaseConfidence[v] == doctest::Approx(result.variantIsMatRef[v] ? mat - pat : pat - mat));
		}
		for (size_t c = 0; c < helpers.numCells(); c++)
		{
			auto Ce = getOptimalCellCe(result, helpers, c, ignoreNothing);
			double mat = getCellLogProb(result, helpers, c, Ce.first, true, ignoreNothing);
			double pat = getCellLogProb(result, helpers, c, Ce.second, false, ignoreNothing);
			CHECK(additions.cellActiveConfidence[c] == doctest::Approx(result.cellIsMatActive[c] ? mat - pat : pat - mat));
		}
	}
	SUBCASE("at a local optimum confidences are not negative")
	{
		for (double conf : additions.variantPhaseConfidence) CHECK(conf >= -0.0001);
		for (double conf : additions.cellActiveConfidence) CHECK(conf >= -0.0001);
	}
	SUBCASE("confidence intervals contain the escape estimates")
	{
		for (size_t v = 0; v < helpers.numVariants(); v++)
		{
			CHECK(additions.variantEscapeCILow[v] <= result.variantEscapeFraction[v] + 1e-9);
			CHECK(additions.variantEscapeCIHigh[v] >= result.variantEscapeFraction[v] - 1e-9);
		}
		for (size_t c = 0; c < helpers.numCells(); c++)
		{
			CHECK(additions.cellEscapeCILow[c] <= result.cellEscapeFraction[c] + 1e-9);
			CHECK(additions.cellEscapeCIHigh[c] >= result.cellEscapeFraction[c] - 1e-9);
		}
	}
	SUBCASE("only nonescape confidence ignores escaping variants")
	{
		EMResult modified = result;
		std::vector<bool> ignoreEscape(helpers.numVariants(), false);
		for (size_t v = 0; v < helpers.numVariants(); v += 2)
		{
			modified.variantEscapeFraction[v] = 0.5;
			ignoreEscape[v] = true;
		}
		for (size_t v = 1; v < helpers.numVariants(); v += 2) modified.variantEscapeFraction[v] = escapeBoundary;
		EMResultAdditions modifiedAdditions = getEMAdditions(data.matches, modified, helpers);
		for (size_t c = 0; c < helpers.numCells(); c++)
		{
			auto Ce = getOptimalCellCe(modified, helpers, c, ignoreEscape);
			double mat = getCellLogProb(modified, helpers, c, Ce.first, true, ignoreEscape);
			double pat = getCellLogProb(modified, helpers, c, Ce.second, false, ignoreEscape);
			CHECK(modifiedAdditions.cellActiveConfidenceOnlyNonescape[c] == doctest::Approx(modified.cellIsMatActive[c] ? mat - pat : pat - mat));
		}
	}
}

TEST_CASE("runEM phases perfectly consistent data correctly")
{
	auto data = simulate(25, 80, 0.3, 0, 91);
	EMOutput output = runEM(data.matches, {}, 0, 20, 0.95, 3, {});
	REQUIRE(output.helpers.numVariants() == 25);
	REQUIRE(output.helpers.numCells() == 80);
	auto variants = resultVariantPhasesInTruthOrder(output, 25);
	auto cells = resultCellActivitiesInTruthOrder(output, 80);
	bool flip = variants[0] != data.variantIsMatRef[0];
	CHECK(accuracy(data.variantIsMatRef, variants, flip) == 1.0);
	CHECK(accuracy(data.cellIsMatActive, cells, flip) == 1.0);
	SUBCASE("output vectors are sized consistently")
	{
		CHECK(output.result.variantEscapeFraction.size() == 25);
		CHECK(output.result.cellEscapeFraction.size() == 80);
		CHECK(output.additions.variantPhaseConfidence.size() == 25);
		CHECK(output.additions.cellActiveConfidence.size() == 80);
	}
	SUBCASE("no escape is estimated")
	{
		for (double e : output.result.variantEscapeFraction) CHECK(e <= 0.01);
		for (double e : output.result.cellEscapeFraction) CHECK(e <= 0.01);
	}
}

TEST_CASE("runEM with forced phases fixes the orientation")
{
	auto data = simulate(25, 80, 0.3, 0, 93);
	for (bool forcedValue : { true, false })
	{
		// force variant 0 to forcedValue, so the orientation is truth if it matches, flipped otherwise
		std::unordered_map<std::string, bool> forced { { variantName(0), forcedValue }, { "X:999999:A:G", true } };
		EMOutput output = runEM(data.matches, forced, 0, 20, 0.95, 2, {});
		auto variants = resultVariantPhasesInTruthOrder(output, 25);
		auto cells = resultCellActivitiesInTruthOrder(output, 80);
		CHECK(variants[0] == forcedValue);
		bool flip = forcedValue != data.variantIsMatRef[0];
		CHECK(accuracy(data.variantIsMatRef, variants, flip) == 1.0);
		CHECK(accuracy(data.cellIsMatActive, cells, flip) == 1.0);
	}
}

TEST_CASE("runEM with contradicting forced phases keeps both forced phases")
{
	auto data = simulate(25, 80, 0.3, 0, 95);
	// variant 0 forced to the truth, variant 1 forced against the truth: no orientation satisfies both
	std::unordered_map<std::string, bool> forced { { variantName(0), data.variantIsMatRef[0] }, { variantName(1), !data.variantIsMatRef[1] } };
	EMOutput output = runEM(data.matches, forced, 0, 20, 0.95, 2, {});
	auto variants = resultVariantPhasesInTruthOrder(output, 25);
	CHECK(variants[0] == data.variantIsMatRef[0]);
	CHECK(variants[1] == !data.variantIsMatRef[1]);
	// the free variants follow one of the two forced variants
	bool flip = variants[2] != data.variantIsMatRef[2];
	size_t consistent = 0;
	for (size_t i = 2; i < 25; i++) if ((variants[i] != data.variantIsMatRef[i]) == flip) consistent += 1;
	CHECK(consistent == 23);
	// the forced variant which disagrees with the rest has a negative phase confidence
	const size_t loser = flip ? 0 : 1;
	const size_t winner = flip ? 1 : 0;
	CHECK(output.additions.variantPhaseConfidence[output.helpers.variantNameToIndex.at(variantName(loser))] < 0);
	CHECK(output.additions.variantPhaseConfidence[output.helpers.variantNameToIndex.at(variantName(winner))] > 0);
}

TEST_CASE("runEM is deterministic for a given seed")
{
	auto data = simulate(20, 50, 0.3, 0.1, 97);
	EMOutput first = runEM(data.matches, {}, 5, 20, 0.95, 2, {});
	EMOutput second = runEM(data.matches, {}, 5, 20, 0.95, 2, {});
	CHECK(first.result.variantIsMatRef == second.result.variantIsMatRef);
	CHECK(first.result.cellIsMatActive == second.result.cellIsMatActive);
	CHECK(first.result.variantEscapeFraction == second.result.variantEscapeFraction);
	CHECK(first.result.cellEscapeFraction == second.result.cellEscapeFraction);
	CHECK(first.additions.variantPhaseConfidence == second.additions.variantPhaseConfidence);
}

TEST_CASE("runEM picks the best of multiple tries")
{
	auto data = simulate(30, 60, 0.15, 0.2, 99);
	EMHelperVariables helpers = getHelpers(data.matches, {});
	double bestSingle = -std::numeric_limits<double>::infinity();
	for (size_t seed = 0; seed < 4; seed++)
	{
		EMOutput single = runEM(data.matches, {}, seed, 20, 0.95, 1, {});
		bestSingle = std::max(bestSingle, getTotalLogProb(single.result, single.helpers));
	}
	EMOutput combined = runEM(data.matches, {}, 0, 20, 0.95, 4, {});
	CHECK(getTotalLogProb(combined.result, combined.helpers) == doctest::Approx(bestSingle));
}

TEST_CASE("runEM picks the best of multiple tries even when all scores are very low")
{
	// runEM starts with bestScore = -100000 but also always takes the first try,
	// so a dataset whose log likelihood is below -100000 must still work.
	// A uniform variant weight scales the log likelihood without changing the optimum.
	auto data = simulate(30, 60, 0.15, 0.2, 99);
	std::vector<std::tuple<size_t, size_t, double>> weights { std::make_tuple<size_t, size_t, double>(0, 1000000000, 10000) };
	EMOutput output = runEM(data.matches, {}, 0, 20, 0.95, 4, weights);
	double score = getTotalLogProb(output.result, output.helpers);
	CHECK(score < -100000);
	double bestSingle = -std::numeric_limits<double>::infinity();
	for (size_t seed = 0; seed < 4; seed++)
	{
		EMOutput single = runEM(data.matches, {}, seed, 20, 0.95, 1, weights);
		bestSingle = std::max(bestSingle, getTotalLogProb(single.result, single.helpers));
	}
	CHECK(score == doctest::Approx(bestSingle));
}

TEST_CASE("runEM estimates escape for an escaping variant")
{
	auto data = simulate(20, 150, 0.6, 0, 103);
	// make variant 0 escape: add inactive allele reads in every cell that covers it
	std::vector<CellMatch> matches;
	for (const auto& m : data.matches)
	{
		matches.push_back(m);
		if (m.variant == variantName(0)) matches.push_back(makeMatch(m.cell, m.variant, !m.alt, m.count));
	}
	EMOutput output = runEM(matches, {}, 0, 20, 0.95, 2, {});
	const size_t escaping = output.helpers.variantNameToIndex.at(variantName(0));
	CHECK(output.result.variantEscapeFraction[escaping] > 0.9);
	for (size_t v = 1; v < 20; v++)
	{
		CHECK(output.result.variantEscapeFraction[output.helpers.variantNameToIndex.at(variantName(v))] < 0.1);
	}
	auto variants = resultVariantPhasesInTruthOrder(output, 20);
	bool flip = variants[1] != data.variantIsMatRef[1];
	size_t correct = 0;
	for (size_t i = 1; i < 20; i++) if ((variants[i] != data.variantIsMatRef[i]) == flip) correct += 1;
	CHECK(correct == 19);
}

TEST_CASE("runEM with an XIST weighted variant expressed from the inactive chrX")
{
	auto data = simulate(20, 100, 0.4, 0, 105);
	// variant 0 is XIST: all its reads come from the inactive chrX
	std::vector<CellMatch> matches;
	for (const auto& m : data.matches)
	{
		if (m.variant == variantName(0))
		{
			matches.push_back(makeMatch(m.cell, m.variant, !m.alt, m.count));
		}
		else
		{
			matches.push_back(m);
		}
	}
	std::vector<std::tuple<size_t, size_t, double>> weights { std::make_tuple<size_t, size_t, double>(0, 50, -1) };
	EMOutput output = runEM(matches, {}, 0, 20, 0.95, 2, weights);
	CHECK(output.helpers.variantWeight[output.helpers.variantNameToIndex.at(variantName(0))] == -1);
	auto variants = resultVariantPhasesInTruthOrder(output, 20);
	auto cells = resultCellActivitiesInTruthOrder(output, 100);
	bool flip = variants[1] != data.variantIsMatRef[1];
	for (size_t i = 1; i < 20; i++) CHECK((variants[i] != data.variantIsMatRef[i]) == flip);
	CHECK(accuracy(data.cellIsMatActive, cells, flip) == 1.0);
	// the XIST variant is phased so that its expressed allele is on the inactive chrX,
	// which is the same phase as the simulated truth since the reads were swapped
	CHECK((variants[0] != data.variantIsMatRef[0]) == flip);
}

TEST_CASE("runEM on a single cell and variant")
{
	std::vector<CellMatch> matches { makeMatch("c", "X:1:A:G", false, 5) };
	EMOutput output = runEM(matches, {}, 0, 20, 0.95, 1, {});
	REQUIRE(output.helpers.numVariants() == 1);
	REQUIRE(output.helpers.numCells() == 1);
	CHECK(output.result.variantIsMatRef[0] == output.result.cellIsMatActive[0]);
	CHECK(output.additions.variantPhaseConfidence[0] > 0);
	CHECK(output.additions.cellActiveConfidence[0] > 0);
}

TEST_CASE("runEM on the simulated test data gives nonnegative confidences")
{
	auto matches = readMatchCounts("./test/testdata/table_lowXi.tsv");
	REQUIRE(matches.size() > 0);
	EMOutput output = runEM(matches, {}, 0, 20, 0.95, 1, {});
	for (double conf : output.additions.variantPhaseConfidence) CHECK(conf >= -0.0001);
	for (double conf : output.additions.cellActiveConfidence) CHECK(conf >= -0.0001);
}

}
