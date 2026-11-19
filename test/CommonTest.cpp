#include <cmath>
#include <limits>
#include <sstream>
#include <unordered_set>
#include "../src/Common.h"
#include "../src/Logger.h"
#include "InternalDeclarations.h"
#include "doctest.h"

TEST_SUITE("Common")
{

TEST_CASE("matHapName and patHapName")
{
	CHECK(matHapName(true) == "mat");
	CHECK(matHapName(false) == "hap1");
	CHECK(patHapName(true) == "pat");
	CHECK(patHapName(false) == "hap2");
}

TEST_CASE("split")
{
	SUBCASE("basic")
	{
		auto parts = split("a\tb\tc", '\t');
		REQUIRE(parts.size() == 3);
		CHECK(parts[0] == "a");
		CHECK(parts[1] == "b");
		CHECK(parts[2] == "c");
	}
	SUBCASE("empty string gives one empty field")
	{
		auto parts = split("", ',');
		REQUIRE(parts.size() == 1);
		CHECK(parts[0] == "");
	}
	SUBCASE("no separator")
	{
		auto parts = split("abc", ',');
		REQUIRE(parts.size() == 1);
		CHECK(parts[0] == "abc");
	}
	SUBCASE("empty fields are kept")
	{
		auto parts = split(",a,,b,", ',');
		REQUIRE(parts.size() == 5);
		CHECK(parts[0] == "");
		CHECK(parts[1] == "a");
		CHECK(parts[2] == "");
		CHECK(parts[3] == "b");
		CHECK(parts[4] == "");
	}
	SUBCASE("only separators")
	{
		auto parts = split(":::", ':');
		REQUIRE(parts.size() == 4);
		for (const auto& part : parts) CHECK(part == "");
	}
	SUBCASE("variant name")
	{
		auto parts = split("X:12345:A:G", ':');
		REQUIRE(parts.size() == 4);
		CHECK(parts[0] == "X");
		CHECK(parts[1] == "12345");
		CHECK(parts[2] == "A");
		CHECK(parts[3] == "G");
	}
}

TEST_CASE("join")
{
	CHECK(join(',', {}) == "");
	CHECK(join(',', {"a"}) == "a");
	CHECK(join(',', {"a", "b", "c"}) == "a,b,c");
	CHECK(join('\t', {"", ""}) == "\t");
	CHECK(join(',', {"", "x", ""}) == ",x,");
}

TEST_CASE("split and join are inverses")
{
	std::vector<std::string> inputs { "", "a", "a,b", ",,", "abc,,def,", "x" };
	for (const auto& input : inputs)
	{
		CHECK(join(',', split(input, ',')) == input);
	}
}

TEST_CASE("lowercase")
{
	CHECK(lowercase("") == "");
	CHECK(lowercase("chrX") == "chrx");
	CHECK(lowercase("CHRX") == "chrx");
	CHECK(lowercase("already lower") == "already lower");
	CHECK(lowercase("MiXeD 123 _-:") == "mixed 123 _-:");
}

TEST_CASE("hasExtension")
{
	CHECK(hasExtension("file.vcf", ".vcf"));
	CHECK(hasExtension("file.vcf.gz", ".vcf.gz"));
	CHECK(hasExtension("file.vcf.gz", ".gz"));
	CHECK(!hasExtension("file.vcf.gz", ".vcf"));
	CHECK(!hasExtension("file.gff3", ".vcf"));
	CHECK(hasExtension("/some/dir/annotation.gff3", ".gff3"));
	SUBCASE("filename must be strictly longer than extension")
	{
		CHECK(!hasExtension(".vcf", ".vcf"));
		CHECK(!hasExtension("vcf", ".vcf"));
		CHECK(!hasExtension("", ".vcf"));
	}
	SUBCASE("case sensitive")
	{
		CHECK(!hasExtension("file.VCF", ".vcf"));
	}
}

TEST_CASE("round")
{
	CHECK(round(3.14159, 2) == "3.14");
	CHECK(round(3.14159, 0) == "3");
	CHECK(round(1.0, 3) == "1.000");
	CHECK(round(0.0, 2) == "0.00");
	CHECK(round(-1.25, 1).substr(0, 3) == "-1.");
	CHECK(round(2.675, 1) == "2.7");
	CHECK(round(99.999, 2) == "100.00");
	CHECK(round(12.5, 4) == "12.5000");
}

TEST_CASE("harmonizeVariantName")
{
	CHECK(harmonizeVariantName("X:100:A:G") == "X:100:A:G");
	CHECK(harmonizeVariantName("chrX:100:A:G") == "X:100:A:G");
	CHECK(harmonizeVariantName("ChrX:100:A:G") == "X:100:A:G");
	CHECK(harmonizeVariantName("x:100:A:G") == "X:100:A:G");
	CHECK(harmonizeVariantName("23:100:A:G") == "X:100:A:G");
	SUBCASE("other chromosomes are unchanged")
	{
		CHECK(harmonizeVariantName("chr1:100:A:G") == "chr1:100:A:G");
		CHECK(harmonizeVariantName("1:100:A:G") == "1:100:A:G");
		CHECK(harmonizeVariantName("chrY:100:A:G") == "chrY:100:A:G");
	}
	SUBCASE("names which are not in chrom:pos:ref:alt format are unchanged")
	{
		CHECK(harmonizeVariantName("chrX:100") == "chrX:100");
		CHECK(harmonizeVariantName("chrX:100:A:G:extra") == "chrX:100:A:G:extra");
		CHECK(harmonizeVariantName("") == "");
		CHECK(harmonizeVariantName("variant") == "variant");
	}
	SUBCASE("idempotent")
	{
		for (std::string name : { "chrX:1:A:G", "23:5:C:T", "chr2:3:G:A" }) CHECK(harmonizeVariantName(harmonizeVariantName(name)) == harmonizeVariantName(name));
	}
}

TEST_CASE("parseVariantPosition")
{
	CHECK(parseVariantPosition("X:12345:A:G") == 12345);
	CHECK(parseVariantPosition("chrX:1:A:G") == 1);
	CHECK(parseVariantPosition("X:0:A:T") == 0);
	CHECK(parseVariantPosition("23:155000000:C:T") == 155000000);
	CHECK(parseVariantPosition("X:12345:") == 12345);
	SUBCASE("names without two colons give position 0")
	{
		CHECK(parseVariantPosition("") == 0);
		CHECK(parseVariantPosition("X") == 0);
		CHECK(parseVariantPosition("X:12345") == 0);
	}
	SUBCASE("larger than 32 bits")
	{
		CHECK(parseVariantPosition("X:5000000000:A:G") == 5000000000ull);
	}
}

TEST_CASE("parseBedRegion valid regions")
{
	SUBCASE("start and end are parsed")
	{
		auto region = parseBedRegion("chrX:100-200");
		CHECK(std::get<1>(region) == 100);
		CHECK(std::get<2>(region) == 200);
	}
	SUBCASE("start zero")
	{
		auto region = parseBedRegion("X:0-3500000");
		CHECK(std::get<0>(region) == "X");
		CHECK(std::get<1>(region) == 0);
		CHECK(std::get<2>(region) == 3500000);
	}
	SUBCASE("chromosome is recognized as chrX by the logic used in main")
	{
		for (std::string name : { "chrX:1-2", "X:1-2", "x:1-2", "CHRX:1-2" })
		{
			auto region = parseBedRegion(name);
			std::string chromosome = std::get<0>(region);
			CHECK((chromosome == "23" || lowercase(chromosome) == "x" || lowercase(chromosome) == "chrx"));
		}
	}
}

TEST_CASE("parseBedRegion returns the full chromosome name")
{
	CHECK(std::get<0>(parseBedRegion("chrX:100-200")) == "chrX");
	CHECK(std::get<0>(parseBedRegion("23:100-200")) == "23");
	CHECK(std::get<0>(parseBedRegion("chr10:100-200")) == "chr10");
}

TEST_CASE("parseBedRegion invalid regions")
{
	const auto invalid = std::make_tuple<std::string, size_t, size_t>("", std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
	CHECK(parseBedRegion("") == invalid);
	CHECK(parseBedRegion("chrX") == invalid);
	CHECK(parseBedRegion("chrX:100") == invalid);
	CHECK(parseBedRegion("chrX:100-") == invalid);
	CHECK(parseBedRegion("chrX:-100") == invalid);
	CHECK(parseBedRegion("chrX:200-100") == invalid);
	CHECK(parseBedRegion("chrX:100-100") == invalid);
	CHECK(parseBedRegion("chrX:1,000-2,000") == invalid);
	CHECK(parseBedRegion("chrX:100-200 ") == invalid);
	CHECK(parseBedRegion(" chrX:100-200") == invalid);
	CHECK(parseBedRegion("chrX\t100\t200") == invalid);
	CHECK(parseBedRegion("chrZ:100-200") == invalid);
}

TEST_CASE("getCellOrder sorts lexicographically")
{
	std::unordered_map<std::string, size_t> nameToIndex { { "c", 0 }, { "a", 1 }, { "b", 2 }, { "aa", 3 } };
	auto order = getCellOrder(nameToIndex);
	REQUIRE(order.size() == 4);
	CHECK(order[0] == "a");
	CHECK(order[1] == "aa");
	CHECK(order[2] == "b");
	CHECK(order[3] == "c");
	CHECK(getCellOrder({}).empty());
}

TEST_CASE("getVariantOrder sorts numerically by position")
{
	std::unordered_map<std::string, size_t> nameToIndex { { "X:100:A:T", 0 }, { "X:20:A:T", 1 }, { "X:3:C:G", 2 }, { "X:1000000:G:A", 3 } };
	auto order = getVariantOrder(nameToIndex);
	REQUIRE(order.size() == 4);
	CHECK(order[0] == "X:3:C:G");
	CHECK(order[1] == "X:20:A:T");
	CHECK(order[2] == "X:100:A:T");
	CHECK(order[3] == "X:1000000:G:A");
	CHECK(getVariantOrder({}).empty());
}

TEST_CASE("getVariantOrder contains every name exactly once")
{
	std::unordered_map<std::string, size_t> nameToIndex;
	for (size_t i = 0; i < 100; i++)
	{
		nameToIndex["X:" + std::to_string((i * 7919) % 1000) + ":A:" + (i % 2 ? "T" : "G")] = i;
	}
	auto order = getVariantOrder(nameToIndex);
	CHECK(order.size() == nameToIndex.size());
	for (size_t i = 1; i < order.size(); i++)
	{
		CHECK(parseVariantPosition(order[i-1]) <= parseVariantPosition(order[i]));
	}
	std::unordered_set<std::string> seen { order.begin(), order.end() };
	CHECK(seen.size() == nameToIndex.size());
}

TEST_CASE("binomial p-values against known values")
{
	// P(X >= 3) for X ~ Bin(10, 0.1)
	CHECK(getBinomialPValueGreaterThan(0.1, 3, 10) == doctest::Approx(0.0701908264));
	// P(X <= 1) for X ~ Bin(10, 0.1)
	CHECK(getBinomialPValueLessThan(0.1, 1, 10) == doctest::Approx(0.7360989291));
	CHECK(getBinomialPValueGreaterThan(0.1, 1, 1) == doctest::Approx(0.1));
	CHECK(getBinomialPValueLessThan(0.1, 0, 2) == doctest::Approx(0.81));
	CHECK(getBinomialPValueGreaterThan(0.5, 2, 3) == doctest::Approx(0.5));
	CHECK(getBinomialPValueLessThan(0.5, 1, 3) == doctest::Approx(0.5));
}

TEST_CASE("binomial p-values edge cases")
{
	CHECK(getBinomialPValueGreaterThan(0.1, 0, 0) == doctest::Approx(1.0));
	CHECK(getBinomialPValueLessThan(0.1, 0, 0) == doctest::Approx(1.0));
	CHECK(getBinomialPValueGreaterThan(0.1, 0, 100) == doctest::Approx(1.0));
	CHECK(getBinomialPValueLessThan(0.1, 100, 100) == doctest::Approx(1.0));
	CHECK(getBinomialPValueGreaterThan(0.1, 100, 100) == doctest::Approx(1e-100));
	CHECK(getBinomialPValueLessThan(0.1, 0, 100) == doctest::Approx(pow(0.9, 100)));
}

TEST_CASE("binomial p-values are complementary")
{
	for (double p : { 0.1, 0.3, 0.5, 0.9 })
	{
		for (size_t n : { 1, 5, 20, 100 })
		{
			for (size_t k = 1; k <= n; k++)
			{
				CHECK(getBinomialPValueGreaterThan(p, k, n) + getBinomialPValueLessThan(p, k-1, n) == doctest::Approx(1.0));
			}
		}
	}
}

TEST_CASE("binomial p-values are monotonic")
{
	for (size_t k = 1; k <= 50; k++)
	{
		CHECK(getBinomialPValueGreaterThan(0.1, k, 50) <= getBinomialPValueGreaterThan(0.1, k-1, 50));
		CHECK(getBinomialPValueLessThan(0.1, k, 50) >= getBinomialPValueLessThan(0.1, k-1, 50));
	}
}

TEST_CASE("binomial p-values are fast for large numbers of trials")
{
	double over = getBinomialPValueGreaterThan(0.1, 3000, 30000);
	double under = getBinomialPValueLessThan(0.1, 3000, 30000);
	CHECK(over + under > 1.0);
}

TEST_CASE("binomial p-values with large number of trials")
{
	double over = getBinomialPValueGreaterThan(0.1, 1000, 10000);
	double under = getBinomialPValueLessThan(0.1, 1000, 10000);
	CHECK(std::isfinite(over));
	CHECK(std::isfinite(under));
	CHECK(over > 0.4);
	CHECK(over < 0.6);
	CHECK(under > 0.4);
	CHECK(under < 0.6);
	CHECK(getBinomialPValueGreaterThan(0.1, 2000, 10000) < 1e-100);
	CHECK(getBinomialPValueLessThan(0.1, 500, 10000) < 1e-50);
}

}

TEST_SUITE("Logger")
{

TEST_CASE("Logger verbosity filters messages")
{
	std::ostringstream stream;
	Logger logger { stream };
	SUBCASE("default verbosity only logs Always")
	{
		logger.log(Logger::LogLevel::Always) << "always" << 1;
		logger.log(Logger::LogLevel::DebugInfo) << "debug";
		logger.log(Logger::LogLevel::DetailedDebugInfo) << "detailed";
		CHECK(stream.str() == "always1");
	}
	SUBCASE("verbosity 1 logs debug info")
	{
		logger.setVerbosity(1);
		logger.log(Logger::LogLevel::Always) << "a";
		logger.log(Logger::LogLevel::DebugInfo) << "b";
		logger.log(Logger::LogLevel::DetailedDebugInfo) << "c";
		CHECK(stream.str() == "ab");
	}
	SUBCASE("verbosity 2 logs everything")
	{
		logger.setVerbosity(2);
		logger.log(Logger::LogLevel::Always) << "a";
		logger.log(Logger::LogLevel::DebugInfo) << "b";
		logger.log(Logger::LogLevel::DetailedDebugInfo) << "c";
		CHECK(stream.str() == "abc");
	}
	SUBCASE("verbosity can be lowered again")
	{
		logger.setVerbosity(2);
		logger.setVerbosity(0);
		logger.log(Logger::LogLevel::DebugInfo) << "b";
		CHECK(stream.str() == "");
	}
}

TEST_CASE("Logger handles std::endl and various types")
{
	std::ostringstream stream;
	Logger logger { stream };
	logger.log(Logger::LogLevel::Always) << "x " << 1.5 << " " << (size_t)7 << std::endl;
	CHECK(stream.str() == "x 1.5 7\n");
	logger.log(Logger::LogLevel::DebugInfo) << "hidden" << std::endl;
	CHECK(stream.str() == "x 1.5 7\n");
}

}
