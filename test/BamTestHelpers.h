#ifndef BamTestHelpers_h
#define BamTestHelpers_h

#ifndef NOBAMTOOLSLIB

// Writes AlignmentTestData::ReadSpecs as BAM with bamtools.

#include <api/BamWriter.h>
#include <api/BamAlignment.h>
#include "AlignmentTestData.h"
#include "doctest.h"

namespace BamTestHelpers
{

using AlignmentTestData::ReadSpec;

inline std::vector<BamTools::CigarOp> parseCigarInternal(const std::string& cigar)
{
	std::vector<BamTools::CigarOp> result;
	if (cigar == "*") return result;
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

inline BamTools::BamAlignment makeAlignment(const ReadSpec& spec, size_t readIndex)
{
	BamTools::BamAlignment aln;
	aln.Name = spec.name != "" ? spec.name : "read" + std::to_string(readIndex);
	aln.RefID = spec.refId;
	aln.Position = spec.position;
	aln.QueryBases = spec.sequence;
	aln.Qualities = std::string(spec.sequence.size(), 'I');
	aln.CigarData = parseCigarInternal(AlignmentTestData::cigarOrDefault(spec));
	aln.MapQuality = 60;
	aln.MateRefID = -1;
	aln.MatePosition = -1;
	aln.InsertSize = 0;
	aln.AlignmentFlag = AlignmentTestData::samFlag(spec);
	for (const std::string& tag : AlignmentTestData::samTags(spec))
	{
		REQUIRE(tag.size() >= 5);
		std::string name = tag.substr(0, 2);
		char type = tag[3];
		std::string value = tag.substr(5);
		if (type == 'Z') aln.AddTag<std::string>(name, "Z", value);
		else if (type == 'i') aln.AddTag<int32_t>(name, "i", std::stoi(value));
		else if (type == 'A') aln.AddTag<int8_t>(name, "A", (int8_t)value[0]);
		else FAIL("unsupported tag type in test: " << tag);
	}
	return aln;
}

inline void writeBam(const std::string& path, const std::vector<ReadSpec>& reads, const std::vector<std::string>& references = AlignmentTestData::defaultReferences)
{
	BamTools::RefVector refs;
	for (const auto& name : references) refs.emplace_back(name, 100000);
	BamTools::BamWriter writer;
	REQUIRE(writer.Open(path, AlignmentTestData::samHeader(references), refs));
	for (size_t i = 0; i < reads.size(); i++)
	{
		bool saved = writer.SaveAlignment(makeAlignment(reads[i], i));
		INFO(writer.GetErrorString());
		REQUIRE(saved);
	}
	writer.Close();
}

}

#endif

#endif
