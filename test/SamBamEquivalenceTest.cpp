// The SAM reader must give exactly the same results as the bamtools BAM reader.
// Each test writes the same reads both as SAM and as BAM, streams both, and
// compares the emitted SNPMatches (in emission order) and the returned counts.

#ifndef NOBAMTOOLSLIB

#include <sstream>
#include "../src/SNPCounter.h"
#include "InternalDeclarations.h"
#include "TestHelpers.h"
#include "AlignmentTestData.h"
#include "BamTestHelpers.h"
#include "doctest.h"

using TestHelpers::TempFile;
using namespace AlignmentTestData;

namespace
{

StreamResult streamBamtools(const std::string& bamPath, const std::string& vcfContents)
{
	auto variants = parseVcf(vcfContents);
	StreamResult result;
	std::tie(result.reads, result.readVariantMatches) = streamSNPsFromBamBamtools(bamPath, variants, [&result](const SNPMatch& match) { result.matches.push_back(match); });
	return result;
}

void checkEquivalent(const StreamResult& bam, const StreamResult& sam)
{
	CHECK(sam.reads == bam.reads);
	CHECK(sam.readVariantMatches == bam.readVariantMatches);
	CHECK(sam.matches.size() == bam.matches.size());
	for (size_t i = 0; i < std::min(sam.matches.size(), bam.matches.size()); i++)
	{
		INFO("match ", i, ": BAM ", describe(bam.matches[i]), " SAM ", describe(sam.matches[i]));
		CHECK(sameMatch(bam.matches[i], sam.matches[i]));
	}
}

// writes the reads as BAM and SAM, and checks that both readers agree
std::pair<StreamResult, StreamResult> compareReaders(const std::vector<ReadSpec>& reads, const std::string& vcf, const std::vector<std::string>& references = defaultReferences)
{
	TempFile bam { ".bam" };
	TempFile sam { ".sam" };
	BamTestHelpers::writeBam(bam.path, reads, references);
	writeSam(sam.path, reads, references);
	auto bamResult = streamBamtools(bam.path, vcf);
	auto samResult = streamSam(sam.path, vcf);
	checkEquivalent(bamResult, samResult);
	return std::make_pair(bamResult, samResult);
}

}

TEST_SUITE("SamBamEquivalence")
{

TEST_CASE("SAM and BAM readers agree on basic counting, UMIs and strands")
{
	ReadSpec reverseRef = basicRead(5, { { 10, 'A' } }, "AAA", "U3");
	reverseRef.reverse = true;
	ReadSpec reverseAlt = basicRead(6, { { 10, 'G' } }, "AAA", "U4");
	reverseAlt.reverse = true;
	auto results = compareReaders({
		basicRead(5, { { 10, 'A' }, { 20, 'T' }, { 30, 'C' } }, "AAA", "U1"),
		reverseRef,
		basicRead(6, { { 10, 'A' }, { 20, 'T' } }, "AAA", "U1"),
		reverseAlt,
		basicRead(7, { { 10, 'A' } }, "AAA", "U2"),
		basicRead(8, { { 30, 'A' } }, "CCC", "U1"),
	}, threeVariantVcf);
	CHECK(results.first.matches.size() == 3);
}

TEST_CASE("SAM and BAM readers agree on read filters")
{
	std::vector<ReadSpec> reads;
	for (size_t i = 0; i < 7; i++)
	{
		ReadSpec spec = basicRead(5, { { 10, 'A' } }, "C" + std::to_string(i), "U1");
		if (i == 0) spec.duplicate = true;
		if (i == 1) spec.secondary = true;
		if (i == 2) spec.failedQC = true;
		if (i == 3) spec.unmapped = true;
		if (i == 4) spec.supplementary = true;
		if (i == 5) { spec.reverse = true; spec.duplicate = true; spec.failedQC = true; }
		reads.push_back(spec);
	}
	reads.push_back(basicRead(5, { { 10, 'A' } }, "", "U1"));
	reads.push_back(basicRead(5, { { 10, 'A' } }, "NOUMI", ""));
	auto results = compareReaders(reads, threeVariantVcf);
	CHECK(results.first.matches.size() == 2);
}

TEST_CASE("SAM and BAM readers agree on barcode and UMI tag fallbacks")
{
	std::vector<ReadSpec> reads;
	auto add = [&reads](const std::string& CB, const std::string& CR, const std::string& UB, const std::string& UR, const std::vector<std::string>& extra)
	{
		ReadSpec spec = basicRead(5 + (int)reads.size() % 3, { { 10, 'A' }, { 20, 'T' } }, CB, UB);
		spec.CR = CR;
		spec.UR = UR;
		spec.extraTags = extra;
		reads.push_back(spec);
	};
	add("AAA", "", "U1", "", {});
	add("-", "CRB", "U1", "", {});
	add("-", "-", "U1", "", {});
	add("", "CRC", "-", "URU", {});
	add("AAA", "", "-", "-", { "XM:Z:XMU" });
	add("AAA", "", "", "", { "XM:Z:XMU2" });
	add("AAA", "", "", "", { "XM:Z:-" });
	add("AAA", "", "U9", "", { "NH:i:1", "RE:A:E", "GN:Z:XIST" });
	add("sample:AAA", "", "a:b", "", {});
	std::sort(reads.begin(), reads.end(), [](const ReadSpec& a, const ReadSpec& b) { return a.position < b.position; });
	compareReaders(reads, threeVariantVcf);
}

TEST_CASE("SAM and BAM readers agree on CIGAR operations")
{
	std::vector<std::pair<std::string, std::string>> cases {
		{ "3S10M", readSequence(13, { { 5, 'G' }, { 2, 'A' } }) },
		{ "2H3S10M4S1H", readSequence(17, { { 5, 'G' } }) },
		{ "2M3D5M", "AAAAAAA" },
		{ "2M2I5M", readSequence(9, { { 4, 'G' } }) },
		{ "1M2P5M", readSequence(6, { { 1, 'G' } }) },
		{ "2M15N5M", readSequence(7, { { 1, 'A' } }) },
		{ "1=10X", readSequence(11, { { 0, 'A' }, { 10, 'T' } }) },
		{ "11M", readSequence(11, { { 0, 'G' }, { 10, 'T' } }) },
		{ "1M1D1M1I1M1N20M", readSequence(24, { { 4, 'C' }, { 14, 'A' } }) },
	};
	for (const auto& t : cases)
	{
		INFO("CIGAR ", t.first);
		ReadSpec read;
		read.position = 10 - (t.first == "1M2P5M" || t.first == "2M15N5M" ? 1 : 0);
		read.cigar = t.first;
		read.sequence = t.second;
		read.CB = "AAA";
		read.UB = "U1";
		compareReaders({ read }, threeVariantVcf);
	}
}

TEST_CASE("SAM and BAM readers agree with multiple chromosomes and unplaced reads")
{
	std::string vcf = vcfHeader + vcfLine("chrX", 11, "A", "G") + vcfLine("chr1", 101, "C", "T") + vcfLine("chr2", 101, "C", "T") + vcfLine("chrY", 5, "C", "T");
	ReadSpec x = basicRead(5, { { 10, 'G' } }, "AAA", "U1");
	ReadSpec one = basicRead(90, { { 100, 'T' } }, "AAA", "U1");
	one.refId = 1;
	ReadSpec two = basicRead(90, { { 100, 'C' } }, "AAA", "U1");
	two.refId = 2;
	ReadSpec three = basicRead(90, { { 100, 'T' } }, "AAA", "U1");
	three.refId = 3;
	ReadSpec unplaced;
	unplaced.refId = -1;
	unplaced.position = -1;
	unplaced.unmapped = true;
	unplaced.sequence = "ACGT";
	unplaced.cigar = "*";
	unplaced.CB = "AAA";
	unplaced.UB = "U1";
	auto results = compareReaders({ x, one, two, three, unplaced, unplaced }, vcf, { "chrX", "chr1", "chr2", "chr3", "chrY" });
	CHECK(results.first.matches.size() == 3);
}

TEST_CASE("SAM and BAM readers agree on many variants")
{
	std::string vcf = vcfHeader;
	std::vector<ReadSpec> reads;
	for (size_t i = 0; i < 50; i++) vcf += vcfLine("chrX", 1000 + i * 100 + 1, "A", "G");
	for (size_t i = 0; i < 50; i++)
	{
		size_t pos = 1000 + i * 100;
		reads.push_back(basicRead(pos - 5, { { pos, 'A' } }, "AAA", "U" + std::to_string(i)));
		reads.push_back(basicRead(pos - 3, { { pos, 'G' } }, "CCC", "U" + std::to_string(i)));
		reads.push_back(basicRead(pos - 1, { { pos, 'G' }, { pos + 100, 'A' } }, "CCC", "V" + std::to_string(i), 110));
	}
	compareReaders(reads, vcf);
}

TEST_CASE("SAM and BAM readers agree on random data")
{
	// Random CIGARs (M,=,X,I,D,N,P,S,H), flags, tags, duplicate variant positions,
	// chromosomes with and without variants, and unplaced unmapped reads.
	for (size_t seed = 0; seed < 100; seed++)
	{
		INFO("seed ", seed);
		auto data = randomDataset(seed);
		auto results = compareReaders(data.reads, data.vcf, data.references);
		if (seed == 0) CHECK(results.first.matches.size() > 10);
	}
}

TEST_CASE("SAM and BAM readers agree on random data with dense coverage")
{
	for (size_t seed = 1000; seed < 1010; seed++)
	{
		INFO("seed ", seed);
		auto data = randomDataset(seed, 2000);
		compareReaders(data.reads, data.vcf, data.references);
	}
}

TEST_CASE("SAM and BAM readers agree on getSNPMatchesFromBamVcf level output")
{
	// getSNPMatchesFromBamVcf builds variant and cell names from the stream; check the full
	// CellMatch output of the SAM stream through the same lambda logic as the BAM path
	for (size_t seed = 0; seed < 10; seed++)
	{
		INFO("seed ", seed);
		auto data = randomDataset(seed + 500);
		TempFile bam { ".bam" };
		TempFile sam { ".sam" };
		TempFile vcf { ".vcf", data.vcf };
		BamTestHelpers::writeBam(bam.path, data.reads, data.references);
		writeSam(sam.path, data.reads, data.references);
		auto variants = readVariantsVcfParseFilename(vcf.path);
		std::vector<SNPMatch> fromSam;
		streamSNPsFromSam(sam.path, variants, [&fromSam](const SNPMatch& m) { fromSam.push_back(m); });
		auto fromBam = countSNPsFromBamVcf(vcf.path, bam.path);
		REQUIRE(fromSam.size() == fromBam.size());
		for (size_t i = 0; i < fromSam.size(); i++) CHECK(sameMatch(fromSam[i], fromBam[i]));
	}
}

TEST_CASE("SAM and BAM readers agree when a tag appears twice")
{
	ReadSpec read = basicRead(5, { { 10, 'A' } }, "FIRST", "U1");
	read.extraTags.push_back("CB:Z:SECOND");
	auto results = compareReaders({ read }, threeVariantVcf);
	REQUIRE(results.second.matches.size() == 1);
	CHECK(results.second.matches[0].barcode == "FIRST");
}

TEST_CASE("SAM and BAM readers agree when a non-tag field looks like a tag")
{
	ReadSpec read = basicRead(5, { { 10, 'A' } }, "", "U1");
	read.name = "CB:Z:FAKE";
	compareReaders({ read }, threeVariantVcf);
}

TEST_CASE("SAM reader with Windows line endings agrees with BAM")
{
	auto reads = std::vector<ReadSpec> {
		basicRead(5, { { 10, 'G' } }, "AAA", "U1"),
		basicRead(5, { { 10, 'A' } }, "AAA", "U2"),
	};
	// put the barcode tag last, so a leftover '\r' would end up in the barcode
	reads[0].CB = "";
	reads[0].extraTags.push_back("CB:Z:AAA");
	TempFile bam { ".bam" };
	TempFile sam { ".sam", samText(reads, defaultReferences, "\r\n") };
	BamTestHelpers::writeBam(bam.path, reads);
	checkEquivalent(streamBamtools(bam.path, threeVariantVcf), streamSam(sam.path, threeVariantVcf));
}

}

#endif
