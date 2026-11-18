#include <map>
#include <sstream>
#include <algorithm>
#include <api/BamWriter.h>
#include <api/BamAlignment.h>
#include <zstr.hpp>
#include "../src/SNPCounter.h"
#include "InternalDeclarations.h"
#include "TestHelpers.h"
#include "doctest.h"

using TestHelpers::TempFile;
using TestHelpers::abortsInChild;

namespace
{

std::vector<BamTools::CigarOp> parseCigar(const std::string& cigar)
{
	std::vector<BamTools::CigarOp> result;
	size_t length = 0;
	for (char c : cigar)
	{
		if (c >= '0' && c <= '9')
		{
			length = length * 10 + (c - '0');
			continue;
		}
		result.emplace_back(c, length);
		length = 0;
	}
	return result;
}

// read sequence of given length filled with N except the given bases at the given offsets
std::string readSequence(size_t length, const std::map<size_t, char>& bases)
{
	std::string result(length, 'N');
	for (const auto& pair : bases)
	{
		result[pair.first] = pair.second;
	}
	return result;
}

struct ReadSpec
{
	int refId = 0;
	int position = 0;
	std::string sequence;
	std::string cigar;
	std::string CB;
	std::string UB;
	std::string CR;
	std::string UR;
	bool reverse = false;
	bool duplicate = false;
	bool secondary = false;
	bool failedQC = false;
	bool unmapped = false;
};

BamTools::BamAlignment makeAlignment(const ReadSpec& spec, size_t readIndex)
{
	BamTools::BamAlignment aln;
	aln.Name = "read" + std::to_string(readIndex);
	aln.RefID = spec.refId;
	aln.Position = spec.position;
	aln.QueryBases = spec.sequence;
	aln.Qualities = std::string(spec.sequence.size(), 'I');
	aln.CigarData = parseCigar(spec.cigar.size() > 0 ? spec.cigar : std::to_string(spec.sequence.size()) + "M");
	aln.MapQuality = 60;
	aln.MateRefID = -1;
	aln.MatePosition = -1;
	aln.InsertSize = 0;
	aln.AlignmentFlag = 0;
	aln.SetIsReverseStrand(spec.reverse);
	aln.SetIsDuplicate(spec.duplicate);
	aln.SetIsPrimaryAlignment(!spec.secondary);
	aln.SetIsFailedQC(spec.failedQC);
	aln.SetIsMapped(!spec.unmapped);
	if (spec.CB != "") aln.AddTag<std::string>("CB", "Z", spec.CB);
	if (spec.CR != "") aln.AddTag<std::string>("CR", "Z", spec.CR);
	if (spec.UB != "") aln.AddTag<std::string>("UB", "Z", spec.UB);
	if (spec.UR != "") aln.AddTag<std::string>("UR", "Z", spec.UR);
	return aln;
}

const std::vector<std::string> defaultReferences { "chrX", "chr1" };

void writeBam(const std::string& path, const std::vector<ReadSpec>& reads, const std::vector<std::string>& references = defaultReferences)
{
	BamTools::RefVector refs;
	std::string header = "@HD\tVN:1.6\tSO:coordinate\n";
	for (const auto& name : references)
	{
		refs.emplace_back(name, 100000);
		header += "@SQ\tSN:" + name + "\tLN:100000\n";
	}
	BamTools::BamWriter writer;
	REQUIRE(writer.Open(path, header, refs));
	for (size_t i = 0; i < reads.size(); i++)
	{
		bool saved = writer.SaveAlignment(makeAlignment(reads[i], i));
		INFO(writer.GetErrorString());
		REQUIRE(saved);
	}
	writer.Close();
}

const std::string vcfHeader = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n";

std::string vcfLine(const std::string& chrom, size_t pos, const std::string& ref, const std::string& alt)
{
	return chrom + "\t" + std::to_string(pos) + "\t.\t" + ref + "\t" + alt + "\t50\tPASS\t.\tGT\t0/1\n";
}

std::vector<SNPMatch> streamAll(const std::string& bamPath, const std::string& vcfContents)
{
	std::istringstream stream { vcfContents };
	auto variants = readVariants(stream);
	for (auto& pair : variants) std::sort(pair.second.begin(), pair.second.end());
	std::vector<SNPMatch> result;
	streamSNPsFromBam(bamPath, variants, [&result](const SNPMatch& match) { result.push_back(match); });
	return result;
}

const SNPMatch* findMatch(const std::vector<SNPMatch>& matches, size_t position, const std::string& barcode)
{
	for (const auto& match : matches)
	{
		if (match.position == position && match.barcode == barcode) return &match;
	}
	return nullptr;
}

// Variants at 1-based positions 11 (A>G), 21 (C>T), 31 (G>A) on chrX, i.e. 0-based 10, 20, 30.
const std::string threeVariantVcf = vcfHeader + vcfLine("chrX", 11, "A", "G") + vcfLine("chrX", 21, "C", "T") + vcfLine("chrX", 31, "G", "A");

ReadSpec basicRead(int position, const std::map<size_t, char>& basesAtRefPositions, const std::string& CB, const std::string& UB, size_t length = 40)
{
	ReadSpec spec;
	spec.position = position;
	std::map<size_t, char> offsets;
	for (const auto& pair : basesAtRefPositions) offsets[pair.first - position] = pair.second;
	spec.sequence = readSequence(length, offsets);
	spec.CB = CB;
	spec.UB = UB;
	return spec;
}

}

TEST_SUITE("SNPCounter")
{

TEST_CASE("readVariants")
{
	std::string vcf = vcfHeader;
	vcf += vcfLine("chrX", 100, "A", "G");
	vcf += vcfLine("chrX", 50, "c", "t");
	vcf += vcfLine("chrX", 60, "AC", "G");
	vcf += vcfLine("chrX", 70, "A", "GT");
	vcf += vcfLine("chrX", 80, "A", "G,T");
	vcf += vcfLine("chrX", 90, ".", "G");
	vcf += vcfLine("chrX", 91, "A", ".");
	vcf += vcfLine("chrX", 92, "N", "A");
	vcf += vcfLine("chrX", 93, "A", "*");
	vcf += vcfLine("chr1", 10, "G", "C");
	vcf += "short\n";
	vcf += "chrX\t5\n";
	std::istringstream stream { vcf };
	auto variants = readVariants(stream);
	CHECK(variants.size() == 2);
	REQUIRE(variants.count("chrX") == 1);
	REQUIRE(variants.at("chrX").size() == 2);
	// positions are converted to 0-based and file order is kept
	CHECK(variants.at("chrX")[0] == std::make_tuple<size_t, char, char>(99, 'A', 'G'));
	CHECK(variants.at("chrX")[1] == std::make_tuple<size_t, char, char>(49, 'C', 'T'));
	REQUIRE(variants.count("chr1") == 1);
	CHECK(variants.at("chr1")[0] == std::make_tuple<size_t, char, char>(9, 'G', 'C'));
}

TEST_CASE("readVariants without a trailing newline reads the last line")
{
	std::istringstream stream { vcfHeader + "chrX\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1" };
	auto variants = readVariants(stream);
	REQUIRE(variants.count("chrX") == 1);
	CHECK(variants.at("chrX").size() == 1);
}

TEST_CASE("readVariants aborts on a malformed #CHROM header")
{
	CHECK(abortsInChild([]() { std::istringstream s { "#CHROM\tSTART\tID\tREF\tALT\n" }; readVariants(s); }));
	CHECK(abortsInChild([]() { std::istringstream s { "#CHROM\tPOS\tID\tALT\tREF\n" }; readVariants(s); }));
	CHECK(abortsInChild([]() { std::istringstream s { "#CHROM\tPOS\tID\tREF\tXXX\n" }; readVariants(s); }));
}

TEST_CASE("readVariantsVcfParseFilename reads .vcf and .vcf.gz and sorts")
{
	std::string vcf = vcfHeader + vcfLine("chrX", 300, "A", "G") + vcfLine("chrX", 100, "C", "T") + vcfLine("chrX", 200, "G", "A");
	TempFile plain { ".vcf", vcf };
	TempFile gz { ".vcf.gz" };
	{
		zstr::ofstream out { gz.path };
		out << vcf;
	}
	auto plainVariants = readVariantsVcfParseFilename(plain.path);
	auto gzVariants = readVariantsVcfParseFilename(gz.path);
	CHECK(plainVariants == gzVariants);
	REQUIRE(plainVariants.at("chrX").size() == 3);
	CHECK(std::get<0>(plainVariants.at("chrX")[0]) == 99);
	CHECK(std::get<0>(plainVariants.at("chrX")[1]) == 199);
	CHECK(std::get<0>(plainVariants.at("chrX")[2]) == 299);
	SUBCASE("unknown extension aborts")
	{
		TempFile wrong { ".bcf", vcf };
		CHECK(abortsInChild([&]() { readVariantsVcfParseFilename(wrong.path); }));
	}
}

TEST_CASE("parseFilenameNoDirectory")
{
	CHECK(parseFilenameNoDirectory("file.bam") == "file.bam");
	CHECK(parseFilenameNoDirectory("/a/b/file.bam") == "file.bam");
	CHECK(parseFilenameNoDirectory("relative/dir/file.bam") == "file.bam");
	CHECK(parseFilenameNoDirectory("C:\\data\\file.bam") == "file.bam");
	CHECK(parseFilenameNoDirectory("dir/") == "");
	CHECK(parseFilenameNoDirectory("") == "");
}

TEST_CASE("sortSNPMatches")
{
	auto make = [](const std::string& chrom, size_t pos, const std::string& barcode, char alt)
	{
		SNPMatch m;
		m.chromosome = chrom;
		m.position = pos;
		m.ref = 'A';
		m.alt = alt;
		m.barcode = barcode;
		m.altFwCount = m.altBwCount = m.refFwCount = m.refBwCount = 0;
		return m;
	};
	std::vector<SNPMatch> matches { make("X", 20, "b", 'G'), make("1", 30, "a", 'G'), make("X", 10, "c", 'G'), make("X", 20, "a", 'T'), make("X", 20, "a", 'C') };
	sortSNPMatches(matches);
	CHECK(matches[0].chromosome == "1");
	CHECK(matches[1].position == 10);
	CHECK(matches[2].barcode == "a");
	CHECK(matches[2].alt == 'C');
	CHECK(matches[3].barcode == "a");
	CHECK(matches[3].alt == 'T');
	CHECK(matches[4].barcode == "b");
	std::vector<SNPMatch> empty;
	sortSNPMatches(empty);
	CHECK(empty.empty());
}

TEST_CASE("parseSNPMatches counts UMIs per barcode")
{
	std::unordered_map<std::string, std::unordered_map<std::string, std::tuple<bool, bool, bool, bool>>> matches;
	matches["AAA"]["umi1"] = std::make_tuple(true, false, false, false);
	matches["AAA"]["umi2"] = std::make_tuple(true, true, false, false);
	matches["AAA"]["umi3"] = std::make_tuple(false, false, true, true);
	matches["CCC"]["umi1"] = std::make_tuple(false, false, false, true);
	auto parsed = parseSNPMatches("chrX", std::make_tuple<size_t, char, char>(99, 'A', 'G'), matches);
	sortSNPMatches(parsed);
	REQUIRE(parsed.size() == 2);
	CHECK(parsed[0].chromosome == "X");
	CHECK(parsed[0].position == 100);
	CHECK(parsed[0].ref == 'A');
	CHECK(parsed[0].alt == 'G');
	CHECK(parsed[0].barcode == "AAA");
	CHECK(parsed[0].refFwCount == 2);
	CHECK(parsed[0].refBwCount == 1);
	CHECK(parsed[0].altFwCount == 1);
	CHECK(parsed[0].altBwCount == 1);
	CHECK(parsed[1].barcode == "CCC");
	CHECK(parsed[1].refFwCount == 0);
	CHECK(parsed[1].refBwCount == 0);
	CHECK(parsed[1].altFwCount == 0);
	CHECK(parsed[1].altBwCount == 1);
	SUBCASE("chromosome names")
	{
		CHECK(parseSNPMatches("X", std::make_tuple<size_t, char, char>(0, 'A', 'G'), matches)[0].chromosome == "X");
		CHECK(parseSNPMatches("chr", std::make_tuple<size_t, char, char>(0, 'A', 'G'), matches)[0].chromosome == "chr");
		CHECK(parseSNPMatches("chr10", std::make_tuple<size_t, char, char>(0, 'A', 'G'), matches)[0].chromosome == "10");
		CHECK(parseSNPMatches("23", std::make_tuple<size_t, char, char>(0, 'A', 'G'), matches)[0].chromosome == "23");
	}
	SUBCASE("empty")
	{
		CHECK(parseSNPMatches("chrX", std::make_tuple<size_t, char, char>(0, 'A', 'G'), {}).empty());
	}
}

TEST_CASE("readBarcode and readUMI")
{
	ReadSpec spec;
	spec.sequence = "ACGT";
	SUBCASE("primary tags")
	{
		spec.CB = "AAA-1";
		spec.CR = "AAB";
		spec.UB = "UMI1";
		spec.UR = "UMI2";
		auto aln = makeAlignment(spec, 0);
		CHECK(readBarcode(aln) == "AAA-1");
		CHECK(readUMI(aln) == "UMI1");
	}
	SUBCASE("fallback to raw tags")
	{
		spec.CR = "AAB";
		spec.UR = "UMI2";
		auto aln = makeAlignment(spec, 0);
		CHECK(readBarcode(aln) == "AAB");
		CHECK(readUMI(aln) == "UMI2");
	}
	SUBCASE("dash is treated as missing")
	{
		spec.CB = "-";
		spec.CR = "AAB";
		spec.UB = "-";
		spec.UR = "UMI2";
		auto aln = makeAlignment(spec, 0);
		CHECK(readBarcode(aln) == "AAB");
		CHECK(readUMI(aln) == "UMI2");
	}
	SUBCASE("only dashes")
	{
		spec.CB = "-";
		spec.CR = "-";
		spec.UB = "-";
		spec.UR = "-";
		auto aln = makeAlignment(spec, 0);
		CHECK(readBarcode(aln) == "");
		CHECK(readUMI(aln) == "");
	}
	SUBCASE("no tags")
	{
		auto aln = makeAlignment(spec, 0);
		CHECK(readBarcode(aln) == "");
		CHECK(readUMI(aln) == "");
	}
	SUBCASE("XM fallback for UMI")
	{
		auto aln = makeAlignment(spec, 0);
		aln.AddTag<std::string>("XM", "Z", "UMI3");
		CHECK(readUMI(aln) == "UMI3");
	}
}

TEST_CASE("streamSNPsFromBam basic counting")
{
	TempFile bam { ".bam" };
	writeBam(bam.path, {
		// cell AAA: ref at 10, alt at 20, neither at 30
		basicRead(5, { { 10, 'A' }, { 20, 'T' }, { 30, 'C' } }, "AAA", "U1"),
		// same UMI again: not counted twice
		basicRead(6, { { 10, 'A' }, { 20, 'T' } }, "AAA", "U1"),
		// different UMI: counted
		basicRead(7, { { 10, 'A' } }, "AAA", "U2"),
		// cell CCC: alt at 30 (BAM can't store lowercase bases so those can't be tested here)
		basicRead(8, { { 30, 'A' } }, "CCC", "U1"),
	});
	auto matches = streamAll(bam.path, threeVariantVcf);
	REQUIRE(matches.size() == 3);
	const SNPMatch* first = findMatch(matches, 11, "AAA");
	REQUIRE(first != nullptr);
	CHECK(first->chromosome == "X");
	CHECK(first->ref == 'A');
	CHECK(first->alt == 'G');
	CHECK(first->refFwCount == 2);
	CHECK(first->refBwCount == 0);
	CHECK(first->altFwCount == 0);
	CHECK(first->altBwCount == 0);
	const SNPMatch* second = findMatch(matches, 21, "AAA");
	REQUIRE(second != nullptr);
	CHECK(second->refFwCount == 0);
	CHECK(second->altFwCount == 1);
	const SNPMatch* third = findMatch(matches, 31, "CCC");
	REQUIRE(third != nullptr);
	CHECK(third->altFwCount == 1);
	CHECK(findMatch(matches, 31, "AAA") == nullptr);
	SUBCASE("output is sorted by position")
	{
		for (size_t i = 1; i < matches.size(); i++)
		{
			CHECK(matches[i-1].position <= matches[i].position);
		}
	}
}

TEST_CASE("streamSNPsFromBam strand counting")
{
	TempFile bam { ".bam" };
	ReadSpec reverseRef = basicRead(5, { { 10, 'A' } }, "AAA", "U1");
	reverseRef.reverse = true;
	ReadSpec reverseAlt = basicRead(5, { { 10, 'G' } }, "AAA", "U2");
	reverseAlt.reverse = true;
	ReadSpec forwardRef = basicRead(6, { { 10, 'A' } }, "AAA", "U1");
	writeBam(bam.path, { reverseRef, reverseAlt, forwardRef });
	auto matches = streamAll(bam.path, threeVariantVcf);
	REQUIRE(matches.size() == 1);
	CHECK(matches[0].refFwCount == 1);
	CHECK(matches[0].refBwCount == 1);
	CHECK(matches[0].altFwCount == 0);
	CHECK(matches[0].altBwCount == 1);
}

TEST_CASE("streamSNPsFromBam filters reads")
{
	TempFile bam { ".bam" };
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
	writeBam(bam.path, { duplicate, secondary, failedQC, noBarcode, noUMI, good, unmapped });
	auto matches = streamAll(bam.path, threeVariantVcf);
	REQUIRE(matches.size() == 1);
	CHECK(matches[0].barcode == "G");
}

TEST_CASE("streamSNPsFromBam handles CIGAR operations")
{
	SUBCASE("soft clip")
	{
		TempFile bam { ".bam" };
		ReadSpec read;
		read.position = 8;
		read.cigar = "3S10M";
		// read index 3 is ref 8, so ref 10 is read index 5
		read.sequence = readSequence(13, { { 5, 'G' }, { 2, 'A' } });
		read.CB = "AAA";
		read.UB = "U1";
		writeBam(bam.path, { read });
		auto matches = streamAll(bam.path, threeVariantVcf);
		REQUIRE(matches.size() == 1);
		CHECK(matches[0].position == 11);
		CHECK(matches[0].altFwCount == 1);
		CHECK(matches[0].refFwCount == 0);
	}
	SUBCASE("deletion over the variant")
	{
		TempFile bam { ".bam" };
		ReadSpec read;
		read.position = 8;
		read.cigar = "2M3D5M";
		read.sequence = "AAAAAAA";
		read.CB = "AAA";
		read.UB = "U1";
		writeBam(bam.path, { read });
		CHECK(streamAll(bam.path, threeVariantVcf).empty());
	}
	SUBCASE("insertion before the variant")
	{
		TempFile bam { ".bam" };
		ReadSpec read;
		read.position = 8;
		read.cigar = "2M2I5M";
		// ref 8,9 = read 0,1; read 2,3 inserted; ref 10 = read 4
		read.sequence = readSequence(9, { { 2, 'A' }, { 3, 'A' }, { 4, 'G' } });
		read.CB = "AAA";
		read.UB = "U1";
		writeBam(bam.path, { read });
		auto matches = streamAll(bam.path, threeVariantVcf);
		REQUIRE(matches.size() == 1);
		CHECK(matches[0].position == 11);
		CHECK(matches[0].altFwCount == 1);
	}
	SUBCASE("splice over a variant")
	{
		TempFile bam { ".bam" };
		ReadSpec read;
		read.position = 9;
		read.cigar = "2M15N5M";
		// ref 9,10 = read 0,1; ref 11..25 skipped; ref 26.. = read 2..
		read.sequence = readSequence(7, { { 1, 'A' } });
		read.CB = "AAA";
		read.UB = "U1";
		writeBam(bam.path, { read });
		auto matches = streamAll(bam.path, threeVariantVcf);
		REQUIRE(matches.size() == 1);
		CHECK(matches[0].position == 11);
		CHECK(matches[0].refFwCount == 1);
	}
	SUBCASE("splice into a variant")
	{
		TempFile bam { ".bam" };
		ReadSpec read;
		read.position = 5;
		read.cigar = "2M13N10M";
		// ref 5,6 = read 0,1; ref 7..19 skipped; ref 20 = read 2
		read.sequence = readSequence(12, { { 2, 'C' } });
		read.CB = "AAA";
		read.UB = "U1";
		writeBam(bam.path, { read });
		auto matches = streamAll(bam.path, threeVariantVcf);
		REQUIRE(matches.size() == 1);
		CHECK(matches[0].position == 21);
		CHECK(matches[0].refFwCount == 1);
	}
	SUBCASE("variant at first and last base of an aligned block")
	{
		TempFile bam { ".bam" };
		ReadSpec read;
		read.position = 10;
		read.cigar = "11M";
		read.sequence = readSequence(11, { { 0, 'G' }, { 10, 'T' } });
		read.CB = "AAA";
		read.UB = "U1";
		writeBam(bam.path, { read });
		auto matches = streamAll(bam.path, threeVariantVcf);
		REQUIRE(matches.size() == 2);
		CHECK(matches[0].position == 11);
		CHECK(matches[0].altFwCount == 1);
		CHECK(matches[1].position == 21);
		CHECK(matches[1].altFwCount == 1);
	}
	SUBCASE("read ending just before a variant")
	{
		TempFile bam { ".bam" };
		ReadSpec read;
		read.position = 0;
		read.sequence = readSequence(10, {});
		read.CB = "AAA";
		read.UB = "U1";
		writeBam(bam.path, { read });
		CHECK(streamAll(bam.path, threeVariantVcf).empty());
	}
}

TEST_CASE("streamSNPsFromBam with multiple chromosomes")
{
	TempFile bam { ".bam" };
	std::string vcf = vcfHeader + vcfLine("chrX", 11, "A", "G") + vcfLine("chr1", 101, "C", "T") + vcfLine("chr2", 101, "C", "T");
	ReadSpec x = basicRead(5, { { 10, 'G' } }, "AAA", "U1");
	ReadSpec one = basicRead(90, { { 100, 'T' } }, "AAA", "U1");
	one.refId = 1;
	ReadSpec two = basicRead(90, { { 100, 'T' } }, "AAA", "U1");
	two.refId = 2;
	ReadSpec three = basicRead(90, { { 100, 'T' } }, "AAA", "U1");
	three.refId = 3;
	writeBam(bam.path, { x, one, two, three }, { "chrX", "chr1", "chr2", "chr3" });
	auto matches = streamAll(bam.path, vcf);
	REQUIRE(matches.size() == 3);
	CHECK(matches[0].chromosome == "X");
	CHECK(matches[0].position == 11);
	CHECK(matches[1].chromosome == "1");
	CHECK(matches[1].position == 101);
	CHECK(matches[2].chromosome == "2");
}

TEST_CASE("streamSNPsFromBam with many reads flushes all variants")
{
	TempFile bam { ".bam" };
	std::string vcf = vcfHeader;
	std::vector<ReadSpec> reads;
	for (size_t i = 0; i < 50; i++)
	{
		vcf += vcfLine("chrX", 1000 + i * 100 + 1, "A", "G");
	}
	for (size_t i = 0; i < 50; i++)
	{
		size_t pos = 1000 + i * 100;
		reads.push_back(basicRead(pos - 5, { { pos, 'A' } }, "AAA", "U" + std::to_string(i)));
		reads.push_back(basicRead(pos - 3, { { pos, 'G' } }, "CCC", "U" + std::to_string(i)));
	}
	writeBam(bam.path, reads);
	auto matches = streamAll(bam.path, vcf);
	REQUIRE(matches.size() == 100);
	for (size_t i = 0; i < 50; i++)
	{
		const SNPMatch* ref = findMatch(matches, 1000 + i * 100 + 1, "AAA");
		const SNPMatch* alt = findMatch(matches, 1000 + i * 100 + 1, "CCC");
		REQUIRE(ref != nullptr);
		REQUIRE(alt != nullptr);
		CHECK(ref->refFwCount == 1);
		CHECK(alt->altFwCount == 1);
	}
	for (size_t i = 1; i < matches.size(); i++)
	{
		CHECK(matches[i-1].position <= matches[i].position);
	}
}

TEST_CASE("streamSNPsFromBam aborts on unsorted BAM")
{
	TempFile unsortedPosition { ".bam" };
	writeBam(unsortedPosition.path, { basicRead(20, {}, "AAA", "U1"), basicRead(5, {}, "AAA", "U1") });
	TempFile unsortedChromosome { ".bam" };
	ReadSpec one = basicRead(5, {}, "AAA", "U1");
	one.refId = 1;
	writeBam(unsortedChromosome.path, { basicRead(5, {}, "AAA", "U1"), one, basicRead(5, {}, "AAA", "U1") });
	CHECK(abortsInChild([&]() { streamAll(unsortedPosition.path, threeVariantVcf); }));
	CHECK(abortsInChild([&]() { streamAll(unsortedChromosome.path, threeVariantVcf); }));
}

TEST_CASE("streamSNPsFromBam throws on missing BAM")
{
	std::unordered_map<std::string, std::vector<std::tuple<size_t, char, char>>> variants;
	CHECK_THROWS_AS(streamSNPsFromBam("/nonexistent/file.bam", variants, [](const SNPMatch&) {}), std::runtime_error);
}

TEST_CASE("countSNPsFromBamVcf")
{
	TempFile bam { ".bam" };
	TempFile vcf { ".vcf", threeVariantVcf };
	writeBam(bam.path, { basicRead(5, { { 10, 'A' }, { 20, 'T' } }, "AAA", "U1") });
	auto matches = countSNPsFromBamVcf(vcf.path, bam.path);
	REQUIRE(matches.size() == 2);
	CHECK(matches[0].position == 11);
	CHECK(matches[0].refFwCount == 1);
	CHECK(matches[1].position == 21);
	CHECK(matches[1].altFwCount == 1);
}

TEST_CASE("getSNPMatchesFromBamVcf")
{
	std::string vcfContents = threeVariantVcf + vcfLine("chr1", 101, "C", "T");
	TempFile vcf { ".vcf", vcfContents };
	TempFile bam1 { ".bam" };
	TempFile bam2 { ".bam" };
	ReadSpec chr1Read = basicRead(90, { { 100, 'T' } }, "AAA", "U1");
	chr1Read.refId = 1;
	writeBam(bam1.path, {
		basicRead(5, { { 10, 'A' } }, "AAA", "U1"),
		basicRead(5, { { 10, 'G' } }, "AAA", "U2"),
		basicRead(5, { { 10, 'G' }, { 20, 'C' } }, "CCC", "U1"),
		chr1Read,
	});
	writeBam(bam2.path, {
		basicRead(5, { { 10, 'A' } }, "AAA", "U1"),
	});
	SUBCASE("single file, only chrX, separate ref and alt entries")
	{
		auto matches = getSNPMatchesFromBamVcf(bam1.path, vcf.path, {});
		REQUIRE(matches.size() == 4);
		CHECK(matches[0].variant == "X:11:A:G");
		CHECK(matches[0].cell == "AAA");
		CHECK(matches[0].alt == false);
		CHECK(matches[0].count == 1);
		CHECK(matches[1].variant == "X:11:A:G");
		CHECK(matches[1].cell == "AAA");
		CHECK(matches[1].alt == true);
		CHECK(matches[1].count == 1);
		CHECK(matches[2].variant == "X:11:A:G");
		CHECK(matches[2].cell == "CCC");
		CHECK(matches[2].alt == true);
		CHECK(matches[3].variant == "X:21:C:T");
		CHECK(matches[3].cell == "CCC");
		CHECK(matches[3].alt == false);
	}
	SUBCASE("barcode whitelist")
	{
		auto matches = getSNPMatchesFromBamVcf(bam1.path, vcf.path, { "CCC" });
		REQUIRE(matches.size() == 2);
		CHECK(matches[0].cell == "CCC");
		CHECK(matches[1].cell == "CCC");
	}
	SUBCASE("multiple files without merging barcodes")
	{
		auto matches = getSNPMatchesFromBamVcf(std::vector<std::string> { bam1.path, bam2.path }, vcf.path, {}, false);
		REQUIRE(matches.size() == 5);
		const std::string name1 = parseFilenameNoDirectory(bam1.path);
		const std::string name2 = parseFilenameNoDirectory(bam2.path);
		CHECK(matches[0].cell == name1 + ":0:AAA");
		CHECK(matches[2].cell == name1 + ":0:CCC");
		CHECK(matches[4].cell == name2 + ":1:AAA");
		CHECK(matches[4].variant == "X:11:A:G");
		CHECK(matches[4].alt == false);
	}
	SUBCASE("multiple files with merging barcodes")
	{
		auto matches = getSNPMatchesFromBamVcf(std::vector<std::string> { bam1.path, bam2.path }, vcf.path, {}, true);
		REQUIRE(matches.size() == 5);
		CHECK(matches[0].cell == "AAA");
		CHECK(matches[4].cell == "AAA");
	}
	SUBCASE("whitelist applies to raw barcodes before file prefixing")
	{
		auto matches = getSNPMatchesFromBamVcf(std::vector<std::string> { bam1.path, bam2.path }, vcf.path, { "AAA" }, false);
		CHECK(matches.size() == 3);
	}
}

}
