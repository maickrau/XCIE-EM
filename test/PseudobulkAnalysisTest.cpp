#include <algorithm>
#include "../src/PseudobulkAnalysis.h"
#include "InternalDeclarations.h"
#include "TestHelpers.h"
#include "OutputBuilder.h"
#include "doctest.h"

using TestHelpers::makeMatch;
using TestHelpers::OutputBuilder;

namespace
{

// Variants: X:100 matRef (conf 3), X:200 patRef (conf 3), X:300 matRef (conf 1)
// Cells: c1 mat active (conf 3), c2 pat active (conf 3), c3 mat active (conf 1)
OutputBuilder makeBuilder()
{
	OutputBuilder builder;
	builder.matches = {
		makeMatch("c1", "X:100:A:G", false, 10),
		makeMatch("c1", "X:100:A:G", true, 1),
		makeMatch("c1", "X:200:C:T", false, 2),
		makeMatch("c1", "X:200:C:T", true, 6),
		makeMatch("c2", "X:100:A:G", false, 3),
		makeMatch("c2", "X:100:A:G", true, 7),
		makeMatch("c2", "X:200:C:T", false, 8),
		makeMatch("c2", "X:300:G:A", false, 4),
		makeMatch("c3", "X:100:A:G", false, 5),
		makeMatch("c3", "X:300:G:A", true, 9),
	};
	builder.variantIsMatRef = { { "X:100:A:G", true }, { "X:200:C:T", false }, { "X:300:G:A", true } };
	builder.variantConfidence = { { "X:100:A:G", 3 }, { "X:200:C:T", 3 }, { "X:300:G:A", 1 } };
	builder.cellIsMatActive = { { "c1", true }, { "c2", false }, { "c3", true } };
	builder.cellConfidence = { { "c1", 3 }, { "c2", 3 }, { "c3", 1 } };
	return builder;
}

void checkCounts(const PseudobulkInfo& info, const std::string& name, size_t matXa, size_t matXi, size_t patXa, size_t patXi)
{
	INFO("pseudobulk item ", info.name);
	CHECK(info.name == name);
	CHECK(info.matXa == matXa);
	CHECK(info.matXi == matXi);
	CHECK(info.patXa == patXa);
	CHECK(info.patXi == patXi);
}

std::vector<std::tuple<size_t, size_t, std::string, std::string>> makeGenes()
{
	return {
		std::make_tuple<size_t, size_t, std::string, std::string>(50, 150, "G1", "N1"),
		std::make_tuple<size_t, size_t, std::string, std::string>(100, 300, "G3", "N3"),
		std::make_tuple<size_t, size_t, std::string, std::string>(150, 250, "G2", "N2"),
	};
}

}

TEST_SUITE("PseudobulkAnalysis")
{

TEST_CASE("getIncludedVariants and getIncludedCells")
{
	EMOutput output = makeBuilder().build();
	auto variants2 = getIncludedVariants(output, 2);
	auto cells2 = getIncludedCells(output, 2);
	CHECK(variants2[output.helpers.variantNameToIndex.at("X:100:A:G")] == true);
	CHECK(variants2[output.helpers.variantNameToIndex.at("X:200:C:T")] == true);
	CHECK(variants2[output.helpers.variantNameToIndex.at("X:300:G:A")] == false);
	CHECK(cells2[output.helpers.cellNameToIndex.at("c1")] == true);
	CHECK(cells2[output.helpers.cellNameToIndex.at("c2")] == true);
	CHECK(cells2[output.helpers.cellNameToIndex.at("c3")] == false);
	SUBCASE("threshold is inclusive")
	{
		auto variants3 = getIncludedVariants(output, 3);
		CHECK(variants3[output.helpers.variantNameToIndex.at("X:100:A:G")] == true);
		auto cells1 = getIncludedCells(output, 1);
		CHECK(std::count(cells1.begin(), cells1.end(), true) == 3);
	}
	SUBCASE("confidence 0 includes everything with nonnegative confidence")
	{
		auto variants0 = getIncludedVariants(output, 0);
		CHECK(std::count(variants0.begin(), variants0.end(), true) == 3);
	}
	SUBCASE("negative confidences are excluded at threshold 0")
	{
		OutputBuilder builder = makeBuilder();
		builder.variantConfidence["X:300:G:A"] = -0.5;
		builder.cellConfidence["c3"] = -0.5;
		EMOutput output2 = builder.build();
		auto variants0 = getIncludedVariants(output2, 0);
		auto cells0 = getIncludedCells(output2, 0);
		CHECK(variants0[output2.helpers.variantNameToIndex.at("X:300:G:A")] == false);
		CHECK(cells0[output2.helpers.cellNameToIndex.at("c3")] == false);
	}
}

TEST_CASE("getVariantPseudobulk")
{
	OutputBuilder builder = makeBuilder();
	EMOutput output = builder.build();
	SUBCASE("confidence 2")
	{
		auto pseudobulk = getVariantPseudobulk(output, builder.matches, 2);
		REQUIRE(pseudobulk.size() == 3);
		checkCounts(pseudobulk[0], "X:100:A:G", 10, 1, 7, 3);
		checkCounts(pseudobulk[1], "X:200:C:T", 6, 2, 8, 0);
		checkCounts(pseudobulk[2], "X:300:G:A", 0, 0, 0, 0);
		CHECK(pseudobulk[0].variantIndex == output.helpers.variantNameToIndex.at("X:100:A:G"));
		CHECK(pseudobulk[2].variantIndex == output.helpers.variantNameToIndex.at("X:300:G:A"));
	}
	SUBCASE("confidence 0")
	{
		auto pseudobulk = getVariantPseudobulk(output, builder.matches, 0);
		REQUIRE(pseudobulk.size() == 3);
		checkCounts(pseudobulk[0], "X:100:A:G", 15, 1, 7, 3);
		checkCounts(pseudobulk[1], "X:200:C:T", 6, 2, 8, 0);
		checkCounts(pseudobulk[2], "X:300:G:A", 0, 9, 0, 4);
	}
	SUBCASE("total coverage is preserved at confidence 0")
	{
		auto pseudobulk = getVariantPseudobulk(output, builder.matches, 0);
		size_t total = 0;
		for (const auto& t : pseudobulk) total += t.matXa + t.matXi + t.patXa + t.patXi;
		size_t expected = 0;
		for (const auto& m : builder.matches) expected += m.count;
		CHECK(total == expected);
	}
	SUBCASE("flipping every phase and activity gives the same active/inactive split with haplotypes swapped")
	{
		OutputBuilder flipped = builder;
		for (auto& pair : flipped.variantIsMatRef) pair.second = !pair.second;
		for (auto& pair : flipped.cellIsMatActive) pair.second = !pair.second;
		auto original = getVariantPseudobulk(output, builder.matches, 0);
		auto pseudobulk = getVariantPseudobulk(flipped.build(), builder.matches, 0);
		REQUIRE(pseudobulk.size() == original.size());
		for (size_t i = 0; i < original.size(); i++)
		{
			checkCounts(pseudobulk[i], original[i].name, original[i].patXa, original[i].patXi, original[i].matXa, original[i].matXi);
		}
	}
}

TEST_CASE("getCellPseudobulk")
{
	OutputBuilder builder = makeBuilder();
	EMOutput output = builder.build();
	SUBCASE("confidence 2 excludes unconfident cells and variants")
	{
		auto pseudobulk = getCellPseudobulk(output, builder.matches, 2);
		REQUIRE(pseudobulk.size() == 2);
		checkCounts(pseudobulk[0], "c1", 16, 3, 0, 0);
		checkCounts(pseudobulk[1], "c2", 15, 3, 0, 0);
	}
	SUBCASE("confidence 0")
	{
		auto pseudobulk = getCellPseudobulk(output, builder.matches, 0);
		REQUIRE(pseudobulk.size() == 3);
		checkCounts(pseudobulk[0], "c1", 16, 3, 0, 0);
		checkCounts(pseudobulk[1], "c2", 15, 7, 0, 0);
		checkCounts(pseudobulk[2], "c3", 5, 9, 0, 0);
	}
	SUBCASE("no confident cells")
	{
		CHECK(getCellPseudobulk(output, builder.matches, 100).empty());
	}
	SUBCASE("empty output")
	{
		OutputBuilder empty;
		CHECK(getCellPseudobulk(empty.build(), {}, 0).empty());
	}
	SUBCASE("many cells with mixed confidence are filtered and sorted")
	{
		OutputBuilder many;
		for (size_t i = 0; i < 50; i++)
		{
			std::string cell = "cell" + std::to_string(i);
			many.matches.push_back(makeMatch(cell, "X:100:A:G", false, i+1));
			many.cellConfidence[cell] = (i % 3 == 0) ? 5 : 0;
			many.cellIsMatActive[cell] = true;
		}
		many.variantIsMatRef["X:100:A:G"] = true;
		many.variantConfidence["X:100:A:G"] = 5;
		auto pseudobulk = getCellPseudobulk(many.build(), many.matches, 2);
		REQUIRE(pseudobulk.size() == 17);
		for (size_t i = 1; i < pseudobulk.size(); i++) CHECK(pseudobulk[i-1].name < pseudobulk[i].name);
		for (const auto& t : pseudobulk)
		{
			size_t index = std::stoull(t.name.substr(4));
			CHECK(index % 3 == 0);
			CHECK(t.matXa == index+1);
			CHECK(t.matXi == 0);
		}
	}
}

TEST_CASE("getVariantGroupPseudobulk")
{
	OutputBuilder builder = makeBuilder();
	EMOutput output = builder.build();
	std::unordered_map<std::string, std::string> grouping { { "c1", "T" }, { "c2", "B" }, { "c3", "T" }, { "notacell", "Z" } };
	SUBCASE("confidence 2")
	{
		auto pseudobulk = getVariantGroupPseudobulk(output, builder.matches, 2, grouping);
		REQUIRE(pseudobulk.size() == 4);
		checkCounts(pseudobulk[0], "X:100:A:G\tB", 0, 0, 7, 3);
		checkCounts(pseudobulk[1], "X:100:A:G\tT", 10, 1, 0, 0);
		checkCounts(pseudobulk[2], "X:200:C:T\tB", 0, 0, 8, 0);
		checkCounts(pseudobulk[3], "X:200:C:T\tT", 6, 2, 0, 0);
	}
	SUBCASE("confidence 0")
	{
		auto pseudobulk = getVariantGroupPseudobulk(output, builder.matches, 0, grouping);
		REQUIRE(pseudobulk.size() == 6);
		checkCounts(pseudobulk[1], "X:100:A:G\tT", 15, 1, 0, 0);
		checkCounts(pseudobulk[4], "X:300:G:A\tB", 0, 0, 0, 4);
		checkCounts(pseudobulk[5], "X:300:G:A\tT", 0, 9, 0, 0);
	}
	SUBCASE("ungrouped cells are left out")
	{
		auto pseudobulk = getVariantGroupPseudobulk(output, builder.matches, 2, { { "c1", "T" } });
		REQUIRE(pseudobulk.size() == 2);
		checkCounts(pseudobulk[0], "X:100:A:G\tT", 10, 1, 0, 0);
		checkCounts(pseudobulk[1], "X:200:C:T\tT", 6, 2, 0, 0);
		CHECK(getVariantGroupPseudobulk(output, builder.matches, 2, {}).empty());
	}
	SUBCASE("group sums equal the ungrouped pseudobulk when every cell is grouped")
	{
		auto grouped = getVariantGroupPseudobulk(output, builder.matches, 0, grouping);
		auto ungrouped = getVariantPseudobulk(output, builder.matches, 0);
		for (const auto& variant : ungrouped)
		{
			size_t matXa = 0, matXi = 0, patXa = 0, patXi = 0;
			for (const auto& t : grouped)
			{
				if (t.name.substr(0, t.name.find('\t')) != variant.name) continue;
				CHECK(t.variantIndex == variant.variantIndex);
				matXa += t.matXa;
				matXi += t.matXi;
				patXa += t.patXa;
				patXi += t.patXi;
			}
			checkCounts(variant, variant.name, matXa, matXi, patXa, patXi);
		}
	}
}

TEST_CASE("getVariantGeneContainment")
{
	EMOutput output = makeBuilder().build();
	auto containment = getVariantGeneContainment(output, makeGenes());
	REQUIRE(containment.size() == 3);
	using GeneList = std::vector<std::pair<std::string, std::string>>;
	CHECK(containment[output.helpers.variantNameToIndex.at("X:100:A:G")] == GeneList { { "G1", "N1" }, { "G3", "N3" } });
	CHECK(containment[output.helpers.variantNameToIndex.at("X:200:C:T")] == GeneList { { "G2", "N2" }, { "G3", "N3" } });
	CHECK(containment[output.helpers.variantNameToIndex.at("X:300:G:A")] == GeneList { { "G3", "N3" } });
	SUBCASE("no genes")
	{
		auto empty = getVariantGeneContainment(output, {});
		REQUIRE(empty.size() == 3);
		for (const auto& t : empty) CHECK(t.empty());
	}
	SUBCASE("gene boundaries are inclusive")
	{
		auto boundary = getVariantGeneContainment(output, {
			std::make_tuple<size_t, size_t, std::string, std::string>(1, 99, "before", "before"),
			std::make_tuple<size_t, size_t, std::string, std::string>(1, 100, "endsAt", "endsAt"),
			std::make_tuple<size_t, size_t, std::string, std::string>(100, 100, "exact", "exact"),
			std::make_tuple<size_t, size_t, std::string, std::string>(101, 199, "between", "between"),
			std::make_tuple<size_t, size_t, std::string, std::string>(200, 1000, "startsAt", "startsAt"),
		});
		CHECK(boundary[output.helpers.variantNameToIndex.at("X:100:A:G")] == GeneList { { "endsAt", "endsAt" }, { "exact", "exact" } });
		CHECK(boundary[output.helpers.variantNameToIndex.at("X:200:C:T")] == GeneList { { "startsAt", "startsAt" } });
		CHECK(boundary[output.helpers.variantNameToIndex.at("X:300:G:A")] == GeneList { { "startsAt", "startsAt" } });
	}
}

TEST_CASE("getGenePseudobulk")
{
	OutputBuilder builder = makeBuilder();
	EMOutput output = builder.build();
	auto containment = getVariantGeneContainment(output, makeGenes());
	SUBCASE("confidence 2")
	{
		auto variantPseudobulk = getVariantPseudobulk(output, builder.matches, 2);
		auto genePseudobulk = getGenePseudobulk(variantPseudobulk, containment);
		// X:300 has zero counts at confidence 2 so its gene set {G3} is omitted
		REQUIRE(genePseudobulk.size() == 2);
		checkCounts(genePseudobulk[0], "G1,G3\tN1,N3", 10, 1, 7, 3);
		checkCounts(genePseudobulk[1], "G2,G3\tN2,N3", 6, 2, 8, 0);
	}
	SUBCASE("confidence 0")
	{
		auto variantPseudobulk = getVariantPseudobulk(output, builder.matches, 0);
		auto genePseudobulk = getGenePseudobulk(variantPseudobulk, containment);
		REQUIRE(genePseudobulk.size() == 3);
		checkCounts(genePseudobulk[2], "G3\tN3", 0, 9, 0, 4);
	}
	SUBCASE("variants with the same gene set are summed")
	{
		std::vector<std::vector<std::pair<std::string, std::string>>> sameGene;
		sameGene.resize(3, { { "G", "N" } });
		auto variantPseudobulk = getVariantPseudobulk(output, builder.matches, 0);
		auto genePseudobulk = getGenePseudobulk(variantPseudobulk, sameGene);
		REQUIRE(genePseudobulk.size() == 1);
		checkCounts(genePseudobulk[0], "G\tN", 21, 12, 15, 7);
	}
	SUBCASE("gene order within a variant doesn't matter")
	{
		std::vector<std::vector<std::pair<std::string, std::string>>> unsortedGenes;
		unsortedGenes.resize(3);
		unsortedGenes[output.helpers.variantNameToIndex.at("X:100:A:G")] = { { "B", "b" }, { "A", "a" } };
		unsortedGenes[output.helpers.variantNameToIndex.at("X:200:C:T")] = { { "A", "a" }, { "B", "b" } };
		auto genePseudobulk = getGenePseudobulk(getVariantPseudobulk(output, builder.matches, 2), unsortedGenes);
		REQUIRE(genePseudobulk.size() == 1);
		checkCounts(genePseudobulk[0], "A,B\ta,b", 16, 3, 15, 3);
	}
	SUBCASE("no genes")
	{
		std::vector<std::vector<std::pair<std::string, std::string>>> noGenes;
		noGenes.resize(3);
		CHECK(getGenePseudobulk(getVariantPseudobulk(output, builder.matches, 0), noGenes).empty());
	}
}

TEST_CASE("getGeneGroupPseudobulk")
{
	OutputBuilder builder = makeBuilder();
	EMOutput output = builder.build();
	auto containment = getVariantGeneContainment(output, makeGenes());
	std::unordered_map<std::string, std::string> grouping { { "c1", "T" }, { "c2", "B" }, { "c3", "T" } };
	SUBCASE("confidence 2")
	{
		auto pseudobulk = getGeneGroupPseudobulk(output, builder.matches, containment, 2, grouping);
		REQUIRE(pseudobulk.size() == 4);
		checkCounts(pseudobulk[0], "G1,G3\tN1,N3\tB", 0, 0, 7, 3);
		checkCounts(pseudobulk[1], "G1,G3\tN1,N3\tT", 10, 1, 0, 0);
		checkCounts(pseudobulk[2], "G2,G3\tN2,N3\tB", 0, 0, 8, 0);
		checkCounts(pseudobulk[3], "G2,G3\tN2,N3\tT", 6, 2, 0, 0);
	}
	SUBCASE("consistent with getGenePseudobulk when every cell is grouped")
	{
		auto grouped = getGeneGroupPseudobulk(output, builder.matches, containment, 0, grouping);
		auto ungrouped = getGenePseudobulk(getVariantPseudobulk(output, builder.matches, 0), containment);
		for (const auto& gene : ungrouped)
		{
			size_t matXa = 0, matXi = 0, patXa = 0, patXi = 0;
			for (const auto& t : grouped)
			{
				if (t.name.substr(0, t.name.rfind('\t')) != gene.name) continue;
				matXa += t.matXa;
				matXi += t.matXi;
				patXa += t.patXa;
				patXi += t.patXi;
			}
			checkCounts(gene, gene.name, matXa, matXi, patXa, patXi);
		}
	}
	SUBCASE("no grouping")
	{
		CHECK(getGeneGroupPseudobulk(output, builder.matches, containment, 0, {}).empty());
	}
}

}
