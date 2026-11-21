#ifndef AlignmentTestData_h
#define AlignmentTestData_h

// Helpers for describing aligned reads in tests and writing them as SAM.
// No bamtools dependency, so this can be used in builds with -DNOBAMTOOLSLIB.
// BamTestHelpers.h writes the same ReadSpecs as BAM.

#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include "../src/SNPCounter.h"
#include "InternalDeclarations.h"
#include "TestHelpers.h"

namespace AlignmentTestData
{

struct ReadSpec
{
	int refId = 0; // -1 for unplaced unmapped reads
	int position = 0; // 0-based
	std::string sequence;
	std::string cigar; // empty means <length>M
	std::string CB;
	std::string UB;
	std::string CR;
	std::string UR;
	// extra tags in SAM format, eg "XM:Z:UMI" or "NH:i:1". Written after the standard tags.
	std::vector<std::string> extraTags;
	bool reverse = false;
	bool duplicate = false;
	bool secondary = false;
	bool failedQC = false;
	bool unmapped = false;
	bool supplementary = false;
	std::string name; // empty means read<index>
};

const std::vector<std::string> defaultReferences { "chrX", "chr1" };

// read sequence of given length filled with N except the given bases at the given offsets
inline std::string readSequence(size_t length, const std::map<size_t, char>& bases)
{
	std::string result(length, 'N');
	for (const auto& pair : bases)
	{
		result.at(pair.first) = pair.second;
	}
	return result;
}

inline std::string cigarOrDefault(const ReadSpec& spec)
{
	return spec.cigar.size() > 0 ? spec.cigar : std::to_string(spec.sequence.size()) + "M";
}

inline uint16_t samFlag(const ReadSpec& spec)
{
	uint16_t flag = 0;
	if (spec.unmapped) flag |= 0x4;
	if (spec.reverse) flag |= 0x10;
	if (spec.secondary) flag |= 0x100;
	if (spec.failedQC) flag |= 0x200;
	if (spec.duplicate) flag |= 0x400;
	if (spec.supplementary) flag |= 0x800;
	return flag;
}

inline std::vector<std::string> samTags(const ReadSpec& spec)
{
	std::vector<std::string> tags;
	if (spec.CB != "") tags.push_back("CB:Z:" + spec.CB);
	if (spec.CR != "") tags.push_back("CR:Z:" + spec.CR);
	if (spec.UB != "") tags.push_back("UB:Z:" + spec.UB);
	if (spec.UR != "") tags.push_back("UR:Z:" + spec.UR);
	tags.insert(tags.end(), spec.extraTags.begin(), spec.extraTags.end());
	return tags;
}

inline std::string samHeader(const std::vector<std::string>& references)
{
	std::string header = "@HD\tVN:1.6\tSO:coordinate\n";
	for (const auto& name : references)
	{
		header += "@SQ\tSN:" + name + "\tLN:100000\n";
	}
	header += "@PG\tID:unittest\tPN:unittest\n";
	return header;
}

inline std::string samLine(const ReadSpec& spec, size_t readIndex, const std::vector<std::string>& references)
{
	std::string line = (spec.name != "" ? spec.name : "read" + std::to_string(readIndex));
	line += "\t" + std::to_string(samFlag(spec));
	line += "\t" + (spec.refId >= 0 ? references.at(spec.refId) : std::string("*"));
	line += "\t" + std::to_string(spec.refId >= 0 ? spec.position + 1 : 0);
	line += "\t60";
	line += "\t" + cigarOrDefault(spec);
	line += "\t*\t0\t0";
	line += "\t" + (spec.sequence.size() > 0 ? spec.sequence : std::string("*"));
	line += "\t" + (spec.sequence.size() > 0 ? std::string(spec.sequence.size(), 'I') : std::string("*"));
	for (const auto& tag : samTags(spec)) line += "\t" + tag;
	return line;
}

inline std::string samText(const std::vector<ReadSpec>& reads, const std::vector<std::string>& references = defaultReferences, const std::string& newline = "\n")
{
	std::string result;
	for (const auto& line : TestHelpers::lines(samHeader(references))) result += line + newline;
	for (size_t i = 0; i < reads.size(); i++)
	{
		result += samLine(reads[i], i, references) + newline;
	}
	return result;
}

inline void writeSam(const std::string& path, const std::vector<ReadSpec>& reads, const std::vector<std::string>& references = defaultReferences)
{
	std::ofstream file { path, std::ios::binary };
	file << samText(reads, references);
}

const std::string vcfHeader = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n";

inline std::string vcfLine(const std::string& chrom, size_t pos, const std::string& ref, const std::string& alt)
{
	return chrom + "\t" + std::to_string(pos) + "\t.\t" + ref + "\t" + alt + "\t50\tPASS\t.\tGT\t0/1\n";
}

// Variants at 1-based positions 11 (A>G), 21 (C>T), 31 (G>A) on chrX, i.e. 0-based 10, 20, 30.
const std::string threeVariantVcf = vcfHeader + vcfLine("chrX", 11, "A", "G") + vcfLine("chrX", 21, "C", "T") + vcfLine("chrX", 31, "G", "A");

inline std::unordered_map<std::string, std::vector<std::tuple<size_t, char, char>>> parseVcf(const std::string& vcfContents)
{
	std::istringstream stream { vcfContents };
	auto variants = readVariants(stream);
	for (auto& pair : variants) std::sort(pair.second.begin(), pair.second.end());
	return variants;
}

struct StreamResult
{
	std::vector<SNPMatch> matches;
	size_t reads = 0;
	size_t readVariantMatches = 0;
};

inline StreamResult streamSam(const std::string& samPath, const std::string& vcfContents)
{
	auto variants = parseVcf(vcfContents);
	StreamResult result;
	std::tie(result.reads, result.readVariantMatches) = streamSNPsFromSam(samPath, variants, [&result](const SNPMatch& match) { result.matches.push_back(match); });
	return result;
}

inline const SNPMatch* findMatch(const std::vector<SNPMatch>& matches, size_t position, const std::string& barcode)
{
	for (const auto& match : matches)
	{
		if (match.position == position && match.barcode == barcode) return &match;
	}
	return nullptr;
}

inline ReadSpec basicRead(int position, const std::map<size_t, char>& basesAtRefPositions, const std::string& CB, const std::string& UB, size_t length = 40)
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

inline bool sameMatch(const SNPMatch& left, const SNPMatch& right)
{
	return left.chromosome == right.chromosome && left.position == right.position && left.ref == right.ref && left.alt == right.alt && left.barcode == right.barcode && left.refFwCount == right.refFwCount && left.refBwCount == right.refBwCount && left.altFwCount == right.altFwCount && left.altBwCount == right.altBwCount;
}

inline std::string describe(const SNPMatch& m)
{
	return m.chromosome + ":" + std::to_string(m.position) + ":" + m.ref + ":" + m.alt + " " + m.barcode + " refFw " + std::to_string(m.refFwCount) + " refBw " + std::to_string(m.refBwCount) + " altFw " + std::to_string(m.altFwCount) + " altBw " + std::to_string(m.altBwCount);
}

struct RandomDataset
{
	std::vector<std::string> references;
	std::vector<ReadSpec> reads;
	std::string vcf;
};

// Random coordinate sorted reads with random CIGARs, flags and tags over random variants.
// References: chrX and chr1 have variants, chr2 has none, chrY is in the header but has no reads.
inline RandomDataset randomDataset(size_t seed, size_t readsPerChromosome = 150)
{
	std::mt19937 rng(seed);
	auto randInt = [&rng](size_t low, size_t high) { return std::uniform_int_distribution<size_t>(low, high)(rng); };
	auto chance = [&rng](double p) { return std::uniform_real_distribution<double>(0, 1)(rng) < p; };
	const std::string bases = "ACGT";
	RandomDataset result;
	result.references = { "chrX", "chr1", "chr2", "chrY" };
	result.vcf = vcfHeader;
	const size_t regionSize = 400;
	for (const std::string chrom : { "chrX", "chr1" })
	{
		size_t numVariants = randInt(5, 25);
		for (size_t i = 0; i < numVariants; i++)
		{
			size_t pos = randInt(1, regionSize);
			char ref = bases[randInt(0, 3)];
			char alt = bases[(bases.find(ref) + randInt(1, 3)) % 4];
			result.vcf += vcfLine(chrom, pos, std::string(1, ref), std::string(1, alt));
			if (chance(0.15))
			{
				// second alt allele at the same position
				char alt2 = bases[(bases.find(alt) + 1) % 4];
				if (alt2 == ref) alt2 = bases[(bases.find(alt2) + 1) % 4];
				result.vcf += vcfLine(chrom, pos, std::string(1, ref), std::string(1, alt2));
			}
		}
	}
	const std::vector<std::string> barcodes { "AAAC-1", "CCCG-1", "GGGT-1", "TTTA-1" };
	const std::vector<std::string> umis { "U1", "U2", "U3", "U4", "U5" };
	for (int refId = 0; refId < 3; refId++)
	{
		std::vector<int> positions;
		for (size_t i = 0; i < readsPerChromosome; i++) positions.push_back(randInt(0, regionSize));
		std::sort(positions.begin(), positions.end());
		for (int position : positions)
		{
			ReadSpec spec;
			spec.refId = refId;
			spec.position = position;
			std::string cigar;
			size_t sequenceLength = 0;
			if (chance(0.1)) cigar += std::to_string(randInt(1, 5)) + "H";
			if (chance(0.3))
			{
				size_t len = randInt(1, 8);
				cigar += std::to_string(len) + "S";
				sequenceLength += len;
			}
			size_t numBlocks = randInt(1, 4);
			for (size_t block = 0; block < numBlocks; block++)
			{
				if (block > 0)
				{
					size_t gapType = randInt(0, 3);
					size_t len = randInt(1, 6);
					if (gapType == 0) { cigar += std::to_string(len) + "I"; sequenceLength += len; }
					if (gapType == 1) cigar += std::to_string(len) + "D";
					if (gapType == 2) cigar += std::to_string(randInt(10, 80)) + "N";
					if (gapType == 3) cigar += std::to_string(len) + "P";
				}
				size_t len = randInt(1, 40);
				const char matchTypes[] = { 'M', 'M', 'M', '=', 'X' };
				cigar += std::to_string(len) + matchTypes[randInt(0, 4)];
				sequenceLength += len;
			}
			if (chance(0.3))
			{
				size_t len = randInt(1, 8);
				cigar += std::to_string(len) + "S";
				sequenceLength += len;
			}
			if (chance(0.1)) cigar += std::to_string(randInt(1, 5)) + "H";
			spec.cigar = cigar;
			for (size_t i = 0; i < sequenceLength; i++) spec.sequence += "ACGTN"[randInt(0, chance(0.9) ? 3 : 4)];
			if (!chance(0.1)) spec.CB = chance(0.05) ? "-" : barcodes[randInt(0, barcodes.size()-1)];
			if (chance(0.15)) spec.CR = barcodes[randInt(0, barcodes.size()-1)];
			if (!chance(0.1)) spec.UB = chance(0.05) ? "-" : umis[randInt(0, umis.size()-1)];
			if (chance(0.15)) spec.UR = umis[randInt(0, umis.size()-1)];
			if (chance(0.1)) spec.extraTags.push_back("XM:Z:" + umis[randInt(0, umis.size()-1)]);
			if (chance(0.3)) spec.extraTags.push_back("NH:i:" + std::to_string(randInt(1, 3)));
			if (chance(0.2)) spec.extraTags.push_back("GN:Z:GENE" + std::to_string(randInt(1, 3)));
			if (chance(0.1)) spec.extraTags.push_back("RE:A:E");
			spec.reverse = chance(0.5);
			spec.duplicate = chance(0.05);
			spec.secondary = chance(0.05);
			spec.failedQC = chance(0.05);
			spec.unmapped = chance(0.03);
			spec.supplementary = chance(0.05);
			result.reads.push_back(spec);
		}
	}
	// unplaced unmapped reads at the end
	size_t unplaced = randInt(0, 3);
	for (size_t i = 0; i < unplaced; i++)
	{
		ReadSpec spec;
		spec.refId = -1;
		spec.position = -1;
		spec.unmapped = true;
		spec.sequence = "ACGTACGTAC";
		spec.cigar = "*";
		spec.CB = barcodes[0];
		spec.UB = umis[0];
		result.reads.push_back(spec);
	}
	return result;
}

}

#endif
