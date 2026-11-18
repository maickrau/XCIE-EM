#include <fstream>
#include "../src/Common.h"
#include "../src/EM.h"
#include "../src/FileHandler.h"
#include "doctest.h"

std::unordered_map<std::string, bool> readTruth(std::string filename)
{
	std::unordered_map<std::string, bool> result;
	std::ifstream file { filename };
	while (file.good())
	{
		std::string line;
		std::getline(file, line);
		auto parts = split(line, '\t');
		result[parts[0]] = (parts[1] == "pat");
	}
	return result;
}

std::pair<double, double> variantAndCellPhasingAccuracy(const std::unordered_map<std::string, bool>& variantTruth, const std::unordered_map<std::string, bool>& cellTruth, const EMOutput& output)
{
	std::vector<std::vector<size_t>> variantMatches;
	std::vector<std::vector<size_t>> cellMatches;
	variantMatches.emplace_back();
	variantMatches.emplace_back();
	variantMatches[0].emplace_back(0);
	variantMatches[0].emplace_back(0);
	variantMatches[1].emplace_back(0);
	variantMatches[1].emplace_back(0);
	cellMatches.emplace_back();
	cellMatches.emplace_back();
	cellMatches[0].emplace_back(0);
	cellMatches[0].emplace_back(0);
	cellMatches[1].emplace_back(0);
	cellMatches[1].emplace_back(0);
	for (const auto& pair : output.helpers.variantNameToIndex)
	{
		bool realHap = variantTruth.at(pair.first);
		bool ranHap = output.result.variantIsMatRef[pair.second];
		variantMatches[realHap ? 1 : 0][ranHap ? 1 : 0] += 1;
	}
	for (const auto& pair : output.helpers.cellNameToIndex)
	{
		bool realHap = cellTruth.at(pair.first);
		bool ranHap = output.result.cellIsMatActive[pair.second];
		cellMatches[realHap ? 1 : 0][ranHap ? 1 : 0] += 1;
	}
	double variantResult;
	double cellResult;
	if (variantMatches[0][0] + variantMatches[1][1] > variantMatches[1][0] + variantMatches[0][1])
	{
		variantResult = (double)(variantMatches[0][0] + variantMatches[1][1]) / (double)(variantMatches[0][0] + variantMatches[1][1] + variantMatches[1][0] + variantMatches[0][1]);
		cellResult = (double)(cellMatches[0][0] + cellMatches[1][1]) / (double)(cellMatches[0][0] + cellMatches[1][1] + cellMatches[1][0] + cellMatches[0][1]);
	}
	else
	{
		variantResult = (double)(variantMatches[0][1] + variantMatches[1][0]) / (double)(variantMatches[0][0] + variantMatches[1][1] + variantMatches[1][0] + variantMatches[0][1]);
		cellResult = (double)(cellMatches[0][1] + cellMatches[1][0]) / (double)(cellMatches[0][0] + cellMatches[1][1] + cellMatches[1][0] + cellMatches[0][1]);
	}
	return std::make_pair(variantResult, cellResult);
}

TEST_CASE("Simulated data with Xi=0 is correctly phased")
{
	auto cellMatches = readMatchCounts("./test/testdata/table_noXi.tsv");
	std::unordered_map<std::string, bool> forcedPhases;
	std::vector<std::tuple<size_t, size_t, double>> regionWeights;
	EMOutput output = runEM(cellMatches, forcedPhases, 0, 20, 0.95, 10, regionWeights);
	auto variantTruth = readTruth("./test/testdata/variants_noXi.tsv");
	auto cellTruth = readTruth("./test/testdata/cells_noXi.tsv");
	auto accuracy = variantAndCellPhasingAccuracy(variantTruth, cellTruth, output);
	CHECK(accuracy.first == 1.0);
	CHECK(accuracy.second == 1.0);
}

TEST_CASE("Simulated data with Xi=0.1 is decently phased")
{
	auto cellMatches = readMatchCounts("./test/testdata/table_lowXi.tsv");
	std::unordered_map<std::string, bool> forcedPhases;
	std::vector<std::tuple<size_t, size_t, double>> regionWeights;
	EMOutput output = runEM(cellMatches, forcedPhases, 0, 20, 0.95, 10, regionWeights);
	auto variantTruth = readTruth("./test/testdata/variants_lowXi.tsv");
	auto cellTruth = readTruth("./test/testdata/cells_lowXi.tsv");
	auto accuracy = variantAndCellPhasingAccuracy(variantTruth, cellTruth, output);
	CHECK(accuracy.first >= 0.95);
	CHECK(accuracy.second >= 0.95);
}

TEST_CASE("Simulated data with Xi=0.4 runs")
{
	auto cellMatches = readMatchCounts("./test/testdata/table_highXi.tsv");
	std::unordered_map<std::string, bool> forcedPhases;
	std::vector<std::tuple<size_t, size_t, double>> regionWeights;
	EMOutput output = runEM(cellMatches, forcedPhases, 0, 20, 0.95, 10, regionWeights);
	auto variantTruth = readTruth("./test/testdata/variants_highXi.tsv");
	auto cellTruth = readTruth("./test/testdata/cells_highXi.tsv");
	auto accuracy = variantAndCellPhasingAccuracy(variantTruth, cellTruth, output);
	CHECK(accuracy.first >= 0.45);
	CHECK(accuracy.second >= 0.45);
}
