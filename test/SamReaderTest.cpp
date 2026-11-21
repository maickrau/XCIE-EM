// Tests for the SAM reader (streamSNPsFromSam). These don't need bamtools, so they also
// run in builds with -DNOBAMTOOLSLIB. Equivalence with the bamtools reader is tested in
// SamBamEquivalenceTest.cpp.

#include <map>
#include <sstream>
#include <algorithm>
#include "../src/SNPCounter.h"
#include "InternalDeclarations.h"
#include "TestHelpers.h"
#include "AlignmentTestData.h"
#include "doctest.h"

using TestHelpers::TempFile;
using TestHelpers::abortsInChild;
using namespace AlignmentTestData;

namespace
{

StreamResult streamReads(const std::vector<ReadSpec>& reads, const std::string& vcf = threeVariantVcf, const std::vector<std::string>& references = defaultReferences)
{
	TempFile sam { ".sam" };
	writeSam(sam.path, reads, references);
	return streamSam(sam.path, vcf);
}

StreamResult streamText(const std::string& samContents, const std::string& vcf = threeVariantVcf)
{
	TempFile sam { ".sam", samContents };
	return streamSam(sam.path, vcf);
}

}

TEST_SUITE("SamReader")
{

TEST_CASE("parseCigar edge cases")
{
	CHECK(parseCigar("*").empty());
	CHECK(parseCigar("").empty());
	auto single = parseCigar("1M");
	REQUIRE(single.size() == 1);
	CHECK(single[0].Type == 'M');
	CHECK(single[0].Length == 1);
	auto large = parseCigar("1000000N");
	REQUIRE(large.size() == 1);
	CHECK(large[0].Length == 1000000);
	auto spliced = parseCigar("5S20M1000N30M2I10M3S");
	REQUIRE(spliced.size() == 7);
	CHECK(spliced[2].Type == 'N');
	CHECK(spliced[2].Length == 1000);
	CHECK(spliced[4].Type == 'I');
	CHECK(spliced[6].Type == 'S');
	CHECK(spliced[6].Length == 3);
}

TEST_CASE("readBarcode and readUMI from SAM fields")
{
	auto fields = [](std::vector<std::string> tags)
	{
		std::vector<std::string> parts { "read1", "0", "chrX", "1", "60", "4M", "*", "0", "0", "ACGT", "IIII" };
		parts.insert(parts.end(), tags.begin(), tags.end());
		return parts;
	};
	SUBCASE("primary tags")
	{
		auto parts = fields({ "CB:Z:AAA-1", "CR:Z:AAB", "UB:Z:UMI1", "UR:Z:UMI2", "XM:Z:UMI3" });
		CHECK(readBarcode(parts) == "AAA-1");
		CHECK(readUMI(parts) == "UMI1");
	}
	SUBCASE("tag order doesn't matter")
	{
		auto parts = fields({ "XM:Z:UMI3", "UR:Z:UMI2", "UB:Z:UMI1", "CR:Z:AAB", "CB:Z:AAA-1" });
		CHECK(readBarcode(parts) == "AAA-1");
		CHECK(readUMI(parts) == "UMI1");
	}
	SUBCASE("fallback to raw tags")
	{
		auto parts = fields({ "CR:Z:AAB", "UR:Z:UMI2" });
		CHECK(readBarcode(parts) == "AAB");
		CHECK(readUMI(parts) == "UMI2");
	}
	SUBCASE("dash is treated as missing")
	{
		auto parts = fields({ "CB:Z:-", "CR:Z:AAB", "UB:Z:-", "UR:Z:-", "XM:Z:UMI3" });
		CHECK(readBarcode(parts) == "AAB");
		CHECK(readUMI(parts) == "UMI3");
	}
	SUBCASE("only dashes or nothing")
	{
		CHECK(readBarcode(fields({ "CB:Z:-", "CR:Z:-" })) == "");
		CHECK(readUMI(fields({ "UB:Z:-", "UR:Z:-", "XM:Z:-" })) == "");
		CHECK(readBarcode(fields({})) == "");
		CHECK(readUMI(fields({})) == "");
	}
	SUBCASE("empty values are treated as missing")
	{
		CHECK(readBarcode(fields({ "CB:Z:", "CR:Z:AAB" })) == "AAB");
		CHECK(readUMI(fields({ "UB:Z:", "UR:Z:UMI2" })) == "UMI2");
	}
	SUBCASE("values can contain colons")
	{
		CHECK(readBarcode(fields({ "CB:Z:sample1:AAA-1" })) == "sample1:AAA-1");
		CHECK(readUMI(fields({ "UB:Z:a:b:c" })) == "a:b:c");
	}
	SUBCASE("tags of other types with the same name are ignored")
	{
		CHECK(readBarcode(fields({ "CB:i:5", "CR:Z:AAB" })) == "AAB");
		CHECK(readUMI(fields({ "UB:i:5" })) == "");
	}
	SUBCASE("similarly named tags are not confused")
	{
		CHECK(readBarcode(fields({ "CY:Z:FFFF", "XCB:Z:AAA" })) == "");
		CHECK(readUMI(fields({ "UY:Z:FFFF", "XUB:Z:AAA", "RG:Z:XM:Z:foo" })) == "");
	}
}

TEST_CASE("streamSNPsFromSam basic counting")
{
	auto result = streamReads({
		basicRead(5, { { 10, 'A' }, { 20, 'T' }, { 30, 'C' } }, "AAA", "U1"),
		basicRead(6, { { 10, 'A' }, { 20, 'T' } }, "AAA", "U1"),
		basicRead(7, { { 10, 'A' } }, "AAA", "U2"),
		basicRead(8, { { 30, 'A' } }, "CCC", "U1"),
	});
	CHECK(result.reads == 4);
	CHECK(result.readVariantMatches == 6);
	REQUIRE(result.matches.size() == 3);
	const SNPMatch* first = findMatch(result.matches, 11, "AAA");
	REQUIRE(first != nullptr);
	CHECK(first->chromosome == "X");
	CHECK(first->ref == 'A');
	CHECK(first->alt == 'G');
	CHECK(first->refFwCount == 2);
	CHECK(first->refBwCount == 0);
	CHECK(first->altFwCount == 0);
	CHECK(first->altBwCount == 0);
	const SNPMatch* second = findMatch(result.matches, 21, "AAA");
	REQUIRE(second != nullptr);
	CHECK(second->altFwCount == 1);
	const SNPMatch* third = findMatch(result.matches, 31, "CCC");
	REQUIRE(third != nullptr);
	CHECK(third->altFwCount == 1);
	CHECK(findMatch(result.matches, 31, "AAA") == nullptr);
	for (size_t i = 1; i < result.matches.size(); i++) CHECK(result.matches[i-1].position <= result.matches[i].position);
}

TEST_CASE("streamSNPsFromSam strand counting")
{
	ReadSpec reverseRef = basicRead(5, { { 10, 'A' } }, "AAA", "U1");
	reverseRef.reverse = true;
	ReadSpec reverseAlt = basicRead(5, { { 10, 'G' } }, "AAA", "U2");
	reverseAlt.reverse = true;
	ReadSpec forwardRef = basicRead(6, { { 10, 'A' } }, "AAA", "U1");
	auto result = streamReads({ reverseRef, reverseAlt, forwardRef });
	REQUIRE(result.matches.size() == 1);
	CHECK(result.matches[0].refFwCount == 1);
	CHECK(result.matches[0].refBwCount == 1);
	CHECK(result.matches[0].altFwCount == 0);
	CHECK(result.matches[0].altBwCount == 1);
}

TEST_CASE("streamSNPsFromSam filters reads")
{
	ReadSpec duplicate = basicRead(5, { { 10, 'A' } }, "D", "U1");
	duplicate.duplicate = true;
	ReadSpec secondary = basicRead(5, { { 10, 'A' } }, "S", "U1");
	secondary.secondary = true;
	ReadSpec failedQC = basicRead(5, { { 10, 'A' } }, "Q", "U1");
	failedQC.failedQC = true;
	ReadSpec noBarcode = basicRead(5, { { 10, 'A' } }, "", "U1");
	ReadSpec noUMI = basicRead(5, { { 10, 'A' } }, "N", "");
	ReadSpec good = basicRead(5, { { 10, 'A' } }, "G", "U1");
	ReadSpec unmapped = basicRead(5, { { 10, 'A' } }, "M", "U1");
	unmapped.unmapped = true;
	auto result = streamReads({ duplicate, secondary, failedQC, noBarcode, noUMI, good, unmapped });
	REQUIRE(result.matches.size() == 1);
	CHECK(result.matches[0].barcode == "G");
	// reads without barcode or UMI are counted as reads but not as matches
	CHECK(result.reads == 3);
	CHECK(result.readVariantMatches == 1);
}

TEST_CASE("streamSNPsFromSam counts supplementary alignments like the bamtools reader")
{
	// bamtools' IsPrimaryAlignment() only checks flag 0x100, so supplementary (0x800) alignments are counted
	ReadSpec supplementary = basicRead(5, { { 10, 'A' } }, "AAA", "U1");
	supplementary.supplementary = true;
	auto result = streamReads({ supplementary });
	REQUIRE(result.matches.size() == 1);
	CHECK(result.matches[0].refFwCount == 1);
}

TEST_CASE("streamSNPsFromSam handles CIGAR operations")
{
	auto single = [](int position, const std::string& cigar, const std::string& sequence)
	{
		ReadSpec read;
		read.position = position;
		read.cigar = cigar;
		read.sequence = sequence;
		read.CB = "AAA";
		read.UB = "U1";
		return streamReads({ read }).matches;
	};
	SUBCASE("soft clip")
	{
		auto matches = single(8, "3S10M", readSequence(13, { { 5, 'G' }, { 2, 'A' } }));
		REQUIRE(matches.size() == 1);
		CHECK(matches[0].position == 11);
		CHECK(matches[0].altFwCount == 1);
		CHECK(matches[0].refFwCount == 0);
	}
	SUBCASE("hard clip does not consume read bases")
	{
		auto matches = single(8, "5H2M5H", readSequence(2, {}) );
		CHECK(matches.empty());
		auto matches2 = single(10, "3H1M", "G");
		REQUIRE(matches2.size() == 1);
		CHECK(matches2[0].altFwCount == 1);
	}
	SUBCASE("deletion over the variant")
	{
		CHECK(single(8, "2M3D5M", "AAAAAAA").empty());
	}
	SUBCASE("insertion before the variant")
	{
		auto matches = single(8, "2M2I5M", readSequence(9, { { 2, 'A' }, { 3, 'A' }, { 4, 'G' } }));
		REQUIRE(matches.size() == 1);
		CHECK(matches[0].position == 11);
		CHECK(matches[0].altFwCount == 1);
	}
	SUBCASE("padding consumes neither")
	{
		auto matches = single(9, "1M2P5M", readSequence(6, { { 1, 'G' } }));
		REQUIRE(matches.size() == 1);
		CHECK(matches[0].altFwCount == 1);
	}
	SUBCASE("splice over a variant")
	{
		auto matches = single(9, "2M15N5M", readSequence(7, { { 1, 'A' } }));
		REQUIRE(matches.size() == 1);
		CHECK(matches[0].position == 11);
		CHECK(matches[0].refFwCount == 1);
	}
	SUBCASE("splice into a variant")
	{
		auto matches = single(5, "2M13N10M", readSequence(12, { { 2, 'C' } }));
		REQUIRE(matches.size() == 1);
		CHECK(matches[0].position == 21);
		CHECK(matches[0].refFwCount == 1);
	}
	SUBCASE("= and X operations")
	{
		auto matches = single(10, "1=10X", readSequence(11, { { 0, 'A' }, { 10, 'T' } }));
		REQUIRE(matches.size() == 2);
		CHECK(matches[0].refFwCount == 1);
		CHECK(matches[1].altFwCount == 1);
	}
	SUBCASE("variant at first and last base of an aligned block")
	{
		auto matches = single(10, "11M", readSequence(11, { { 0, 'G' }, { 10, 'T' } }));
		REQUIRE(matches.size() == 2);
		CHECK(matches[0].position == 11);
		CHECK(matches[1].position == 21);
	}
	SUBCASE("read ending just before a variant")
	{
		CHECK(single(0, "10M", readSequence(10, {})).empty());
	}
	SUBCASE("missing CIGAR counts nothing")
	{
		CHECK(single(8, "*", readSequence(10, { { 2, 'A' } })).empty());
	}
}

TEST_CASE("streamSNPsFromSam lowercase bases")
{
	// SAM can contain lowercase bases even though BAM can't
	auto result = streamReads({ basicRead(5, { { 10, 'a' }, { 20, 't' }, { 30, 'g' } }, "AAA", "U1") });
	REQUIRE(result.matches.size() == 3);
	CHECK(result.matches[0].refFwCount == 1);
	CHECK(result.matches[1].altFwCount == 1);
	CHECK(result.matches[2].refFwCount == 1);
}

TEST_CASE("streamSNPsFromSam with multiple chromosomes")
{
	std::string vcf = vcfHeader + vcfLine("chrX", 11, "A", "G") + vcfLine("chr1", 101, "C", "T") + vcfLine("chr2", 101, "C", "T");
	ReadSpec x = basicRead(5, { { 10, 'G' } }, "AAA", "U1");
	ReadSpec one = basicRead(90, { { 100, 'T' } }, "AAA", "U1");
	one.refId = 1;
	ReadSpec two = basicRead(90, { { 100, 'T' } }, "AAA", "U1");
	two.refId = 2;
	ReadSpec three = basicRead(90, { { 100, 'T' } }, "AAA", "U1");
	three.refId = 3;
	auto result = streamReads({ x, one, two, three }, vcf, { "chrX", "chr1", "chr2", "chr3" });
	REQUIRE(result.matches.size() == 3);
	CHECK(result.matches[0].chromosome == "X");
	CHECK(result.matches[1].chromosome == "1");
	CHECK(result.matches[2].chromosome == "2");
	// reads on chr3 have no variants and are not counted
	CHECK(result.reads == 3);
}

TEST_CASE("streamSNPsFromSam with many reads flushes all variants")
{
	std::string vcf = vcfHeader;
	std::vector<ReadSpec> reads;
	for (size_t i = 0; i < 50; i++) vcf += vcfLine("chrX", 1000 + i * 100 + 1, "A", "G");
	for (size_t i = 0; i < 50; i++)
	{
		size_t pos = 1000 + i * 100;
		reads.push_back(basicRead(pos - 5, { { pos, 'A' } }, "AAA", "U" + std::to_string(i)));
		reads.push_back(basicRead(pos - 3, { { pos, 'G' } }, "CCC", "U" + std::to_string(i)));
	}
	auto result = streamReads(reads, vcf);
	REQUIRE(result.matches.size() == 100);
	for (size_t i = 0; i < 50; i++)
	{
		const SNPMatch* ref = findMatch(result.matches, 1000 + i * 100 + 1, "AAA");
		const SNPMatch* alt = findMatch(result.matches, 1000 + i * 100 + 1, "CCC");
		REQUIRE(ref != nullptr);
		REQUIRE(alt != nullptr);
		CHECK(ref->refFwCount == 1);
		CHECK(alt->altFwCount == 1);
	}
}

TEST_CASE("streamSNPsFromSam text format edge cases")
{
	const std::string read = "r1\t0\tchrX\t6\t60\t40M\t*\t0\t0\t" + readSequence(40, { { 5, 'G' } }) + "\t" + std::string(40, 'I') + "\tCB:Z:AAA\tUB:Z:U1";
	SUBCASE("no header")
	{
		auto result = streamText(read + "\n");
		REQUIRE(result.matches.size() == 1);
		CHECK(result.matches[0].altFwCount == 1);
	}
	SUBCASE("no trailing newline")
	{
		auto result = streamText(samHeader(defaultReferences) + read);
		REQUIRE(result.matches.size() == 1);
	}
	SUBCASE("empty lines are skipped")
	{
		auto result = streamText(samHeader(defaultReferences) + "\n\n" + read + "\n\n");
		REQUIRE(result.matches.size() == 1);
	}
	SUBCASE("header only")
	{
		auto result = streamText(samHeader(defaultReferences));
		CHECK(result.matches.empty());
		CHECK(result.reads == 0);
	}
	SUBCASE("empty file")
	{
		auto result = streamText("");
		CHECK(result.matches.empty());
	}
	SUBCASE("optional fields of many types")
	{
		auto result = streamText(read + "\tNH:i:1\tAS:f:-1.5\tRE:A:E\tZZ:B:c,1,2,3\tHH:H:1AE301\n");
		REQUIRE(result.matches.size() == 1);
	}
}

TEST_CASE("streamSNPsFromSam handles Windows line endings" * doctest::should_fail())
{
	// BUG: lines are not stripped of '\r', so with CRLF line endings the last field of each
	// line keeps a trailing '\r'. When the CB tag is last the barcode becomes "AAA\r".
	auto reads = std::vector<ReadSpec> { basicRead(5, { { 10, 'G' } }, "", "U1"), basicRead(5, { { 10, 'A' } }, "", "U2") };
	// barcode tag last, so a leftover '\r' would end up in the barcode
	reads[0].extraTags.push_back("CB:Z:AAA");
	reads[1].extraTags.push_back("CB:Z:AAA");
	TempFile sam { ".sam", samText(reads, defaultReferences, "\r\n") };
	auto result = streamSam(sam.path, threeVariantVcf);
	REQUIRE(result.matches.size() == 1);
	CHECK(result.matches[0].barcode == "AAA");
	CHECK(result.matches[0].refFwCount == 1);
	CHECK(result.matches[0].altFwCount == 1);
}

TEST_CASE("streamSNPsFromSam aborts on unsorted input")
{
	TempFile unsortedPosition { ".sam" };
	writeSam(unsortedPosition.path, { basicRead(20, {}, "AAA", "U1"), basicRead(5, {}, "AAA", "U1") });
	TempFile unsortedChromosome { ".sam" };
	ReadSpec one = basicRead(5, {}, "AAA", "U1");
	one.refId = 1;
	writeSam(unsortedChromosome.path, { basicRead(5, {}, "AAA", "U1"), one, basicRead(5, {}, "AAA", "U1") });
	CHECK(abortsInChild([&]() { streamSam(unsortedPosition.path, threeVariantVcf); }));
	CHECK(abortsInChild([&]() { streamSam(unsortedChromosome.path, threeVariantVcf); }));
	SUBCASE("unsorted reads on a chromosome without variants are not checked, like the bamtools reader")
	{
		ReadSpec a = basicRead(20, {}, "AAA", "U1");
		a.refId = 1;
		ReadSpec b = basicRead(5, {}, "AAA", "U1");
		b.refId = 1;
		CHECK_NOTHROW(streamReads({ a, b }));
	}
}

TEST_CASE("streamSNPsFromSam missing file and BAM input")
{
	auto variants = parseVcf(threeVariantVcf);
	CHECK_THROWS_AS(streamSNPsFromSam("/nonexistent/file.sam", variants, [](const SNPMatch&) {}), std::runtime_error);
	CHECK(abortsInChild([&]() { streamSNPsFromSam("/nonexistent/file.bam", variants, [](const SNPMatch&) {}); }));
}

TEST_CASE("streamSNPsFromSam on random data is internally consistent")
{
	for (size_t seed = 0; seed < 20; seed++)
	{
		auto data = randomDataset(seed);
		TempFile sam { ".sam" };
		writeSam(sam.path, data.reads, data.references);
		auto result = streamSam(sam.path, data.vcf);
		size_t countSum = 0;
		for (size_t i = 0; i < result.matches.size(); i++)
		{
			const auto& m = result.matches[i];
			CHECK(m.refFwCount + m.refBwCount + m.altFwCount + m.altBwCount > 0);
			CHECK(m.barcode != "");
			CHECK(m.barcode != "-");
			CHECK((m.chromosome == "X" || m.chromosome == "1"));
			countSum += m.refFwCount + m.refBwCount + m.altFwCount + m.altBwCount;
			if (i > 0 && result.matches[i-1].chromosome == m.chromosome) CHECK(result.matches[i-1].position <= m.position);
		}
		// UMI deduplication can only reduce the number of counted read-variant matches
		CHECK(countSum <= result.readVariantMatches);
	}
}

#ifdef NOBAMTOOLSLIB

// In builds without bamtools the BAM entry points read SAM files

TEST_CASE("countSNPsFromBamVcf reads SAM without bamtools")
{
	TempFile sam { ".sam" };
	TempFile vcf { ".vcf", threeVariantVcf };
	writeSam(sam.path, { basicRead(5, { { 10, 'A' }, { 20, 'T' } }, "AAA", "U1") });
	auto matches = countSNPsFromBamVcf(vcf.path, sam.path);
	REQUIRE(matches.size() == 2);
	CHECK(matches[0].position == 11);
	CHECK(matches[0].refFwCount == 1);
	CHECK(matches[1].position == 21);
	CHECK(matches[1].altFwCount == 1);
}

TEST_CASE("getSNPMatchesFromBamVcf reads SAM without bamtools")
{
	TempFile vcf { ".vcf", threeVariantVcf + vcfLine("chr1", 101, "C", "T") };
	TempFile sam1 { ".sam" };
	TempFile sam2 { ".sam" };
	ReadSpec chr1Read = basicRead(90, { { 100, 'T' } }, "AAA", "U1");
	chr1Read.refId = 1;
	writeSam(sam1.path, {
		basicRead(5, { { 10, 'A' } }, "AAA", "U1"),
		basicRead(5, { { 10, 'G' } }, "AAA", "U2"),
		basicRead(5, { { 10, 'G' }, { 20, 'C' } }, "CCC", "U1"),
		chr1Read,
	});
	writeSam(sam2.path, { basicRead(5, { { 10, 'A' } }, "AAA", "U1") });
	auto single = getSNPMatchesFromBamVcf(sam1.path, vcf.path, {});
	REQUIRE(single.size() == 4);
	CHECK(single[0].variant == "X:11:A:G");
	CHECK(single[0].cell == "AAA");
	CHECK(single[3].variant == "X:21:C:T");
	auto unmerged = getSNPMatchesFromBamVcf(std::vector<std::string> { sam1.path, sam2.path }, vcf.path, {}, false);
	REQUIRE(unmerged.size() == 5);
	CHECK(unmerged[4].cell == parseFilenameNoDirectory(sam2.path) + ":1:AAA");
	CHECK(getSNPMatchesFromBamVcf(sam1.path, vcf.path, { "CCC" }).size() == 2);
}

TEST_CASE("BAM input aborts without bamtools")
{
	TempFile vcf { ".vcf", threeVariantVcf };
	CHECK(abortsInChild([&]() { getSNPMatchesFromBamVcf("input.bam", vcf.path, {}); }));
}

#endif

}
