#include <algorithm>
#include <unordered_set>
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
		CHECK(filterOutHomozygousSites({}, 0.15).empty());
		CHECK(filterOutHomozygousSites({}, 0).empty());
	}
	SUBCASE("balanced heterozygous variant is kept")
	{
		std::vector<CellMatch> raw { makeMatch("c1", "X:1:A:G", false, 5), makeMatch("c1", "X:1:A:G", true, 5) };
		for (double fraction : { 0.0, 0.1, 0.15, 0.25, 0.5 })
		{
			INFO("fraction ", fraction);
			auto result = filterOutHomozygousSites(raw, fraction);
			REQUIRE(result.size() == 2);
			CHECK(result[0].alt == false);
			CHECK(result[1].alt == true);
		}
	}
	SUBCASE("ref-only and alt-only variants are removed")
	{
		std::vector<CellMatch> raw {
			makeMatch("c1", "X:1:A:G", false, 5),
			makeMatch("c2", "X:1:A:G", false, 3),
			makeMatch("c1", "X:2:A:G", true, 5),
		};
		for (double fraction : { 0.001, 0.15, 0.5 })
		{
			INFO("fraction ", fraction);
			CHECK(filterOutHomozygousSites(raw, fraction).empty());
		}
	}
	SUBCASE("fraction 0 keeps everything, even homozygous sites")
	{
		std::vector<CellMatch> raw {
			makeMatch("c1", "X:1:A:G", false, 5),
			makeMatch("c1", "X:2:A:G", true, 5),
			makeMatch("c1", "X:3:A:G", false, 1),
			makeMatch("c1", "X:3:A:G", true, 1000),
		};
		CHECK(filterOutHomozygousSites(raw, 0).size() == raw.size());
	}
	SUBCASE("fraction 0.5 keeps only perfectly balanced sites")
	{
		std::vector<CellMatch> raw {
			makeMatch("c1", "X:1:A:G", false, 7),
			makeMatch("c2", "X:1:A:G", true, 7),
			makeMatch("c1", "X:2:A:G", false, 7),
			makeMatch("c2", "X:2:A:G", true, 8),
			makeMatch("c1", "X:3:A:G", false, 500),
			makeMatch("c2", "X:3:A:G", true, 499),
		};
		auto result = filterOutHomozygousSites(raw, 0.5);
		CHECK(countForVariant(result, "X:1:A:G") == 2);
		CHECK(countForVariant(result, "X:2:A:G") == 0);
		CHECK(countForVariant(result, "X:3:A:G") == 0);
	}
	SUBCASE("counts are summed across cells before the fraction test")
	{
		// each cell is monoallelic but the variant as a whole is balanced
		std::vector<CellMatch> raw {
			makeMatch("c1", "X:1:A:G", false, 5),
			makeMatch("c2", "X:1:A:G", true, 3),
			makeMatch("c3", "X:1:A:G", true, 2),
		};
		CHECK(filterOutHomozygousSites(raw, 0.5).size() == 3);
		// one cell with a lot of ref expression makes the site look homozygous overall
		raw.push_back(makeMatch("c4", "X:1:A:G", false, 100));
		CHECK(filterOutHomozygousSites(raw, 0.15).empty());
	}
	SUBCASE("minor allele fraction boundary is inclusive")
	{
		// these fractions are exactly representable, so the boundary is exact
		std::vector<CellMatch> raw {
			// ref 1 of 4: fraction exactly 0.25
			makeMatch("c1", "X:1:A:G", false, 1),
			makeMatch("c1", "X:1:A:G", true, 3),
			// alt 1 of 4: fraction exactly 0.25
			makeMatch("c1", "X:2:A:G", false, 3),
			makeMatch("c1", "X:2:A:G", true, 1),
			// ref 2 of 9: fraction 0.222 just under 0.25
			makeMatch("c1", "X:3:A:G", false, 2),
			makeMatch("c1", "X:3:A:G", true, 7),
			// alt 25 of 99: fraction 0.2525 just over 0.25
			makeMatch("c1", "X:4:A:G", false, 74),
			makeMatch("c1", "X:4:A:G", true, 25),
		};
		auto result = filterOutHomozygousSites(raw, 0.25);
		CHECK(countForVariant(result, "X:1:A:G") == 2);
		CHECK(countForVariant(result, "X:2:A:G") == 2);
		CHECK(countForVariant(result, "X:3:A:G") == 0);
		CHECK(countForVariant(result, "X:4:A:G") == 2);
		auto result2 = filterOutHomozygousSites(raw, 0.125);
		CHECK(countForVariant(result2, "X:3:A:G") == 2);
		auto result3 = filterOutHomozygousSites(raw, 0.375);
		CHECK(result3.empty());
	}
	SUBCASE("threshold is a fraction of total coverage, not a ratio of minor to major allele")
	{
		// ref 2 alt 10: ratio 0.2 but fraction 2/12 = 0.167
		std::vector<CellMatch> raw {
			makeMatch("c1", "X:1:A:G", false, 2),
			makeMatch("c1", "X:1:A:G", true, 10),
		};
		CHECK(filterOutHomozygousSites(raw, 0.2).empty());
		CHECK(filterOutHomozygousSites(raw, 0.15).size() == 2);
		CHECK(filterOutHomozygousSites(raw, 1.0 / 6.0 - 0.001).size() == 2);
		CHECK(filterOutHomozygousSites(raw, 1.0 / 6.0 + 0.001).empty());
	}
	SUBCASE("default fraction 0.15 is more lenient than the old minor/major ratio of 0.2")
	{
		// ratio 0.2 corresponds to fraction 1/6. Fractions between 0.15 and 1/6 are now kept.
		std::vector<CellMatch> raw {
			// 3 of 19: fraction 0.158, ratio 0.1875
			makeMatch("c1", "X:1:A:G", false, 3),
			makeMatch("c1", "X:1:A:G", true, 16),
			// 16 of 100: fraction 0.16, ratio 0.19
			makeMatch("c1", "X:2:A:G", false, 84),
			makeMatch("c1", "X:2:A:G", true, 16),
			// 14 of 100: fraction 0.14, removed
			makeMatch("c1", "X:3:A:G", false, 86),
			makeMatch("c1", "X:3:A:G", true, 14),
			// 1 of 10: fraction 0.1, removed
			makeMatch("c1", "X:4:A:G", false, 1),
			makeMatch("c1", "X:4:A:G", true, 9),
		};
		auto result = filterOutHomozygousSites(raw, 0.15);
		CHECK(countForVariant(result, "X:1:A:G") == 2);
		CHECK(countForVariant(result, "X:2:A:G") == 2);
		CHECK(countForVariant(result, "X:3:A:G") == 0);
		CHECK(countForVariant(result, "X:4:A:G") == 0);
		auto old = filterOutHomozygousSites(raw, 1.0 / 6.0);
		CHECK(countForVariant(old, "X:1:A:G") == 0);
		CHECK(countForVariant(old, "X:2:A:G") == 0);
	}
	SUBCASE("matches the exact integer rule for many counts and thresholds")
	{
		// with fraction k/8 a site is kept iff 8 * min(ref, alt) >= k * (ref + alt)
		std::vector<CellMatch> raw;
		for (size_t ref = 0; ref <= 30; ref++)
		{
			for (size_t alt = 0; alt <= 30; alt++)
			{
				if (ref + alt == 0) continue;
				std::string name = "X:" + std::to_string(ref * 100 + alt + 1) + ":A:G";
				if (ref > 0) raw.push_back(makeMatch("c1", name, false, ref));
				if (alt > 0) raw.push_back(makeMatch("c2", name, true, alt));
			}
		}
		for (size_t k = 0; k <= 4; k++)
		{
			double fraction = k / 8.0;
			auto result = filterOutHomozygousSites(raw, fraction);
			std::unordered_set<std::string> kept;
			for (const auto& m : result) kept.insert(m.variant);
			for (size_t ref = 0; ref <= 30; ref++)
			{
				for (size_t alt = 0; alt <= 30; alt++)
				{
					if (ref + alt == 0) continue;
					std::string name = "X:" + std::to_string(ref * 100 + alt + 1) + ":A:G";
					bool expected = 8 * std::min(ref, alt) >= k * (ref + alt);
					INFO("ref ", ref, " alt ", alt, " fraction ", fraction);
					CHECK((kept.count(name) == 1) == expected);
				}
			}
		}
	}
	SUBCASE("higher thresholds keep a subset of lower thresholds")
	{
		std::vector<CellMatch> raw;
		for (size_t i = 0; i < 40; i++)
		{
			std::string name = "X:" + std::to_string(i + 1) + ":A:G";
			raw.push_back(makeMatch("c1", name, false, 1 + i));
			raw.push_back(makeMatch("c2", name, true, 40 - i));
		}
		size_t previous = raw.size() + 1;
		for (double fraction : { 0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5 })
		{
			auto result = filterOutHomozygousSites(raw, fraction);
			CHECK(result.size() <= previous);
			previous = result.size();
		}
	}
	SUBCASE("order and content of kept matches is preserved")
	{
		std::vector<CellMatch> raw {
			makeMatch("c3", "X:5:A:G", true, 4),
			makeMatch("c1", "X:9:C:T", true, 100),
			makeMatch("c2", "X:5:A:G", false, 3),
		};
		auto result = filterOutHomozygousSites(raw, 0.15);
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
