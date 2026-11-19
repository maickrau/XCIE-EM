#include <algorithm>
#include "../src/AlleleSpecificExpression.h"
#include "TestHelpers.h"
#include "doctest.h"

using TestHelpers::makeMatch;

namespace
{

size_t countForVariant(const std::vector<CellMatch>& matches, const std::string& variant)
{
	return std::count_if(matches.begin(), matches.end(), [&variant](const CellMatch& m) { return m.variant == variant; });
}

}

TEST_SUITE("AlleleSpecificExpression")
{

TEST_CASE("filterOutHomozygousSites")
{
	SUBCASE("empty input")
	{
		CHECK(filterOutHomozygousSites({}).empty());
	}
	SUBCASE("balanced heterozygous variant is kept")
	{
		std::vector<CellMatch> raw { makeMatch("c1", "X:1:A:G", false, 5), makeMatch("c1", "X:1:A:G", true, 5) };
		auto result = filterOutHomozygousSites(raw);
		REQUIRE(result.size() == 2);
		CHECK(result[0].alt == false);
		CHECK(result[1].alt == true);
	}
	SUBCASE("ref-only and alt-only variants are removed")
	{
		std::vector<CellMatch> raw {
			makeMatch("c1", "X:1:A:G", false, 5),
			makeMatch("c2", "X:1:A:G", false, 3),
			makeMatch("c1", "X:2:A:G", true, 5),
		};
		CHECK(filterOutHomozygousSites(raw).empty());
	}
	SUBCASE("counts are summed across cells before the ratio test")
	{
		// each cell is monoallelic but the variant as a whole is balanced
		std::vector<CellMatch> raw {
			makeMatch("c1", "X:1:A:G", false, 5),
			makeMatch("c2", "X:1:A:G", true, 5),
		};
		CHECK(filterOutHomozygousSites(raw).size() == 2);
	}
	SUBCASE("minor allele ratio boundary at 0.2")
	{
		std::vector<CellMatch> raw {
			// ref 2 alt 10: 2 >= 10*0.2 so kept
			makeMatch("c1", "X:1:A:G", false, 2),
			makeMatch("c1", "X:1:A:G", true, 10),
			// ref 1 alt 10: 1 < 2 so removed
			makeMatch("c1", "X:2:A:G", false, 1),
			makeMatch("c1", "X:2:A:G", true, 10),
			// symmetric with alt minor: ref 10 alt 2 kept
			makeMatch("c1", "X:3:A:G", false, 10),
			makeMatch("c1", "X:3:A:G", true, 2),
			// ref 10 alt 1 removed
			makeMatch("c1", "X:4:A:G", false, 10),
			makeMatch("c1", "X:4:A:G", true, 1),
		};
		auto result = filterOutHomozygousSites(raw);
		CHECK(countForVariant(result, "X:1:A:G") == 2);
		CHECK(countForVariant(result, "X:2:A:G") == 0);
		CHECK(countForVariant(result, "X:3:A:G") == 2);
		CHECK(countForVariant(result, "X:4:A:G") == 0);
	}
	SUBCASE("order and content of kept matches is preserved")
	{
		std::vector<CellMatch> raw {
			makeMatch("c3", "X:5:A:G", true, 4),
			makeMatch("c1", "X:9:C:T", true, 100),
			makeMatch("c2", "X:5:A:G", false, 3),
		};
		auto result = filterOutHomozygousSites(raw);
		REQUIRE(result.size() == 2);
		CHECK(result[0].cell == "c3");
		CHECK(result[0].count == 4);
		CHECK(result[1].cell == "c2");
		CHECK(result[1].count == 3);
	}
}

TEST_CASE("excludeRegions")
{
	std::vector<CellMatch> raw {
		makeMatch("c1", "X:5:A:G", false, 1),
		makeMatch("c1", "X:10:A:G", false, 1),
		makeMatch("c1", "X:15:A:G", false, 1),
		makeMatch("c1", "X:20:A:G", false, 1),
		makeMatch("c1", "X:25:A:G", false, 1),
		makeMatch("c1", "X:100:A:G", false, 1),
	};
	SUBCASE("no regions keeps everything")
	{
		CHECK(excludeRegions(raw, {}).size() == raw.size());
	}
	SUBCASE("region boundaries are inclusive")
	{
		auto result = excludeRegions(raw, { { 10, 20 } });
		REQUIRE(result.size() == 3);
		CHECK(result[0].variant == "X:5:A:G");
		CHECK(result[1].variant == "X:25:A:G");
		CHECK(result[2].variant == "X:100:A:G");
	}
	SUBCASE("multiple sorted regions")
	{
		auto result = excludeRegions(raw, { { 0, 5 }, { 15, 15 }, { 90, 1000 } });
		REQUIRE(result.size() == 3);
		CHECK(result[0].variant == "X:10:A:G");
		CHECK(result[1].variant == "X:20:A:G");
		CHECK(result[2].variant == "X:25:A:G");
	}
	SUBCASE("overlapping sorted regions")
	{
		auto result = excludeRegions(raw, { { 0, 12 }, { 8, 22 } });
		REQUIRE(result.size() == 2);
		CHECK(result[0].variant == "X:25:A:G");
		CHECK(result[1].variant == "X:100:A:G");
	}
	SUBCASE("region containing a smaller region")
	{
		auto result = excludeRegions(raw, { { 0, 50 }, { 10, 11 } });
		REQUIRE(result.size() == 1);
		CHECK(result[0].variant == "X:100:A:G");
	}
	SUBCASE("region with no variants")
	{
		CHECK(excludeRegions(raw, { { 30, 90 } }).size() == raw.size());
	}
}

TEST_CASE("filterToValidBarcodes")
{
	std::vector<CellMatch> raw {
		makeMatch("AAA", "X:5:A:G", false, 1),
		makeMatch("CCC", "X:5:A:G", true, 2),
		makeMatch("AAA", "X:6:A:G", true, 3),
		makeMatch("GGG", "X:6:A:G", true, 4),
	};
	SUBCASE("empty whitelist removes everything")
	{
		CHECK(filterToValidBarcodes(raw, {}).empty());
	}
	SUBCASE("only whitelisted barcodes are kept in order")
	{
		auto result = filterToValidBarcodes(raw, { "AAA", "GGG", "TTT" });
		REQUIRE(result.size() == 3);
		CHECK(result[0].count == 1);
		CHECK(result[1].count == 3);
		CHECK(result[2].count == 4);
	}
	SUBCASE("barcodes are matched exactly")
	{
		CHECK(filterToValidBarcodes(raw, { "aaa", "AA", "AAA-1" }).empty());
	}
}

TEST_CASE("hasVariantOverlap")
{
	std::vector<CellMatch> matches { makeMatch("c1", "X:5:A:G", false, 1), makeMatch("c2", "X:6:A:G", true, 2) };
	CHECK(!hasVariantOverlap(matches, {}));
	CHECK(!hasVariantOverlap({}, { { "X:5:A:G", true } }));
	CHECK(!hasVariantOverlap(matches, { { "X:7:A:G", true }, { "1:5:A:G", false } }));
	CHECK(hasVariantOverlap(matches, { { "X:7:A:G", true }, { "X:6:A:G", false } }));
	CHECK(hasVariantOverlap(matches, { { "X:5:A:G", false } }));
}

}
