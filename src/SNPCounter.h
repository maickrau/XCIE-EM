#ifndef SNPCounter_h
#define SNPCounter_h

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <limits>
#include <vector>
#include <tuple>
#include <string>
#include <fstream>
#include <cassert>
#include <cstdint>
#include "AlleleSpecificExpression.h"
#include "Logger.h"
#include "Common.h"

struct SNPMatch
{
public:
	std::string chromosome;
	size_t position;
	char ref;
	char alt;
	std::string barcode;
	size_t altFwCount;
	size_t altBwCount;
	size_t refFwCount;
	size_t refBwCount;
};

struct CigarData
{
	char Type;
	size_t Length;
};

std::vector<SNPMatch> countSNPsFromBamVcf(const std::string& vcfFile, const std::string& bamFile);
std::vector<CellMatch> getSNPMatchesFromBamVcf(const std::string& bamFile, const std::string& vcfFile, const std::unordered_set<std::string>& barcodeWhitelist);
std::vector<CellMatch> getSNPMatchesFromBamVcf(const std::vector<std::string>& bamFiles, const std::string& vcfFile, const std::unordered_set<std::string>& barcodeWhitelist, const bool mergeBamBarcodes);
std::unordered_map<std::string, std::vector<std::tuple<size_t, char, char>>> readVariantsVcfParseFilename(const std::string& vcfFile);
std::unordered_map<std::string, std::vector<std::tuple<size_t, char, char>>> readVariantsVcfGz(const std::string& vcfGzFile);
void sortSNPMatches(std::vector<SNPMatch>& matches);
std::vector<SNPMatch> parseSNPMatches(const std::string& currentChromosome, const std::tuple<size_t, char, char>& currentVariant, const std::unordered_map<std::string, std::unordered_map<std::string, std::tuple<bool, bool, bool, bool>>>& matches);


std::string readBarcode(const std::vector<std::string>& parts);
std::string readUMI(const std::vector<std::string>& parts);
std::vector<CigarData> parseCigar(const std::string& str);

template <typename F>
std::tuple<size_t, size_t> streamSNPsFromSam(const std::string& bamFile, const std::unordered_map<std::string, std::vector<std::tuple<size_t, char, char>>>& refvariants, F callback)
{
	if (bamFile.size() >= 4 && bamFile.substr(bamFile.size()-4) == ".bam")
	{
		std::cerr << "BAM files are not supported. Recompile without -DNOBAMTOOLSLIB to enable BAM file support." << std::endl;
		std::abort();
	}
	Logger::Log.log(Logger::LogLevel::DebugInfo) << "read sam from " << bamFile << std::endl;
	std::ifstream file { bamFile };
	if (!file.good()) throw std::runtime_error { "Cannot open sam" };
	std::string line;
	std::string currentChromosome = "";
	std::unordered_map<size_t, std::unordered_map<std::string, std::unordered_map<std::string, std::tuple<bool, bool, bool, bool>>>> matches;
	std::unordered_set<std::string> processedChromosomes;
	size_t currentChromosomeSNPIndex = 0;
	size_t reads = 0;
	size_t readsWithTag = 0;
	size_t matchCount = 0;
	size_t lastReadStartPosition = 0;
	while (file.good())
	{
		std::getline(file, line);
		if (line.size() == 0) continue;
		if (line[0] == '@') continue;
		auto parts = split(line, '\t');
		uint64_t flags = std::stoull(parts[1]);
		bool secondary = (flags & 256) != 0;
		if (secondary) continue;
		bool duplicate = (flags & 1024) != 0;
		if (duplicate) continue;
		bool isFailedQC = (flags & 512) != 0;
		if (isFailedQC) continue;
		bool isUnMapped = (flags & 4) != 0;
		if (isUnMapped) continue;
		std::string read = parts[0];
		std::string chromosome = parts[2];
		std::string queryBases = parts[9];
		bool isReverseStrand = (flags & 16) != 0;
		if (chromosome != currentChromosome)
		{
			std::vector<SNPMatch> parsed;
			for (const auto& pospair : matches)
			{
				auto parts = parseSNPMatches(currentChromosome, refvariants.at(currentChromosome)[pospair.first], pospair.second);
				parsed.insert(parsed.end(), parts.begin(), parts.end());
			}
			sortSNPMatches(parsed);
			for (size_t i = 0; i < parsed.size(); i++)
			{
				callback(parsed[i]);
			}
			matches.clear();
			currentChromosome = chromosome;
			Logger::Log.log(Logger::LogLevel::DebugInfo) << "begin chromosome " << currentChromosome << std::endl;
			if (processedChromosomes.count(currentChromosome) == 1)
			{
				std::cerr << "Input BAM is not sorted! BAM needs to be sorted." << std::endl;
				std::abort();
			}
			processedChromosomes.emplace(currentChromosome);
			currentChromosomeSNPIndex = 0;
			lastReadStartPosition = 0;
		}
		if (refvariants.count(currentChromosome) == 0) continue;
		reads += 1;
		const auto& variants = refvariants.at(currentChromosome);
		size_t position = std::stoull(parts[3])-1;
		if (position < lastReadStartPosition)
		{
			std::cerr << "Input BAM is not sorted! BAM needs to be sorted." << std::endl;
			std::abort();
		}
		lastReadStartPosition = position;
		while (currentChromosomeSNPIndex < variants.size() && position > std::get<0>(variants[currentChromosomeSNPIndex]))
		{
			if (matches.count(currentChromosomeSNPIndex) == 1)
			{
				std::vector<SNPMatch> parsed = parseSNPMatches(currentChromosome, refvariants.at(currentChromosome)[currentChromosomeSNPIndex], matches.at(currentChromosomeSNPIndex));
				sortSNPMatches(parsed);
				for (size_t i = 0; i < parsed.size(); i++)
				{
					callback(parsed[i]);
				}
				matches.erase(currentChromosomeSNPIndex);
			}
			currentChromosomeSNPIndex += 1;
		}
		size_t j = currentChromosomeSNPIndex;
		size_t refpos = position;
		size_t readpos = 0;
		std::string barcode = readBarcode(parts);
		std::string umi = readUMI(parts);
		if (barcode == "") continue;
		if (umi == "") continue;
		readsWithTag += 1;
		auto cigarData = parseCigar(parts[5]);
		for (const auto& cigar : cigarData)
		{
			if (cigar.Type != 'M' && cigar.Type != '=' && cigar.Type != 'X')
			{
				if (cigar.Type == 'I' || cigar.Type == 'S')
				{
					readpos += cigar.Length;
				}
				if (cigar.Type == 'D' || cigar.Type == 'N')
				{
					refpos += cigar.Length;
				}
				continue;
			}
			assert(cigar.Type == 'M' || cigar.Type == '=' || cigar.Type == 'X');
			while (j < variants.size() && refpos > std::get<0>(variants[j])) j += 1;
			while (j < variants.size() && refpos+cigar.Length > std::get<0>(variants[j]))
			{
				std::string chromosome = currentChromosome;
				size_t variantPositionInRead = readpos+(std::get<0>(variants[j])-refpos);
				assert(variantPositionInRead < queryBases.size());
				char readNucleotide = queryBases[variantPositionInRead];
				if (readNucleotide == 'a') readNucleotide = 'A';
				if (readNucleotide == 'c') readNucleotide = 'C';
				if (readNucleotide == 'g') readNucleotide = 'G';
				if (readNucleotide == 't') readNucleotide = 'T';
				if (readNucleotide == std::get<1>(variants[j]))
				{
					matchCount += 1;
					if (isReverseStrand)
					{
						std::get<1>(matches[j][barcode][umi]) = true;
					}
					else
					{
						std::get<0>(matches[j][barcode][umi]) = true;
					}
				}
				if (readNucleotide == std::get<2>(variants[j]))
				{
					matchCount += 1;
					if (isReverseStrand)
					{
						std::get<3>(matches[j][barcode][umi]) = true;
					}
					else
					{
						std::get<2>(matches[j][barcode][umi]) = true;
					}
				}
				j += 1;
			}
			readpos += cigar.Length;
			refpos += cigar.Length;
		}
	}
	std::vector<SNPMatch> parsed;
	for (const auto& pospair : matches)
	{
		auto parts = parseSNPMatches(currentChromosome, refvariants.at(currentChromosome)[pospair.first], pospair.second);
		parsed.insert(parsed.end(), parts.begin(), parts.end());
	}
	sortSNPMatches(parsed);
	for (size_t i = 0; i < parsed.size(); i++)
	{
		callback(parsed[i]);
	}
	return std::make_tuple(reads, matchCount);
}

#ifdef NOBAMTOOLSLIB

template <typename F>
std::tuple<size_t, size_t> streamSNPsFromBam(const std::string& bamFile, const std::unordered_map<std::string, std::vector<std::tuple<size_t, char, char>>>& refvariants, F callback)
{
	return streamSNPsFromSam(bamFile, refvariants, callback);
}

#else

#include <api/BamReader.h>

std::string readBarcode(const BamTools::BamAlignment& aln);
std::string readUMI(const BamTools::BamAlignment& aln);

template <typename F>
std::tuple<size_t, size_t> streamSNPsFromBamBamtools(const std::string& bamFile, const std::unordered_map<std::string, std::vector<std::tuple<size_t, char, char>>>& refvariants, F callback)
{
	Logger::Log.log(Logger::LogLevel::DebugInfo) << "read bam from " << bamFile << std::endl;
	BamTools::BamReader reader;
	if (!reader.Open(bamFile))
	{
		throw std::runtime_error { "Cannot open bam" };
	}
	if (!reader.IsOpen())
	{
		throw std::runtime_error { "Cannot open bam" };
	}
	auto references = reader.GetReferenceData();
	BamTools::BamAlignment aln;
	size_t currentChromosomeSNPIndex = 0;
	size_t refId = std::numeric_limits<size_t>::max();
	std::string currentChromosome = "";
	std::unordered_map<size_t, std::unordered_map<std::string, std::unordered_map<std::string, std::tuple<bool, bool, bool, bool>>>> matches;
	std::unordered_set<std::string> processedChromosomes;
	size_t reads = 0;
	size_t readsWithTag = 0;
	size_t matchCount = 0;
	size_t lastReadStartPosition = 0;
	while (reader.GetNextAlignment(aln))
	{
		if (!aln.IsPrimaryAlignment()) continue;
		if (aln.IsDuplicate()) continue;
		if (aln.IsFailedQC()) continue;
		if (!aln.IsMapped()) continue;
//		if (aln.MapQuality < 20) continue;
		assert(aln.RefID >= 0);
		if ((size_t)aln.RefID != refId)
		{
			std::vector<SNPMatch> parsed;
			for (const auto& pospair : matches)
			{
				auto parts = parseSNPMatches(currentChromosome, refvariants.at(currentChromosome)[pospair.first], pospair.second);
				parsed.insert(parsed.end(), parts.begin(), parts.end());
			}
			sortSNPMatches(parsed);
			for (size_t i = 0; i < parsed.size(); i++)
			{
				callback(parsed[i]);
			}
			matches.clear();
			refId = aln.RefID;
			currentChromosome = references[refId].RefName;
			Logger::Log.log(Logger::LogLevel::DebugInfo) << "begin chromosome " << currentChromosome << std::endl;
			if (processedChromosomes.count(currentChromosome) == 1)
			{
				std::cerr << "Input BAM is not sorted! BAM needs to be sorted." << std::endl;
				std::abort();
			}
			processedChromosomes.emplace(currentChromosome);
			currentChromosomeSNPIndex = 0;
			lastReadStartPosition = 0;
		}
		if (refvariants.count(currentChromosome) == 0) continue;
		reads += 1;
		const auto& variants = refvariants.at(currentChromosome);
		assert(aln.Position >= 0);
		if ((size_t)aln.Position < lastReadStartPosition)
		{
			std::cerr << "Input BAM is not sorted! BAM needs to be sorted." << std::endl;
			std::abort();
		}
		lastReadStartPosition = (size_t)aln.Position;
		while (currentChromosomeSNPIndex < variants.size() && (size_t)aln.Position > std::get<0>(variants[currentChromosomeSNPIndex]))
		{
			if (matches.count(currentChromosomeSNPIndex) == 1)
			{
				std::vector<SNPMatch> parsed = parseSNPMatches(currentChromosome, refvariants.at(currentChromosome)[currentChromosomeSNPIndex], matches.at(currentChromosomeSNPIndex));
				sortSNPMatches(parsed);
				for (size_t i = 0; i < parsed.size(); i++)
				{
					callback(parsed[i]);
				}
				matches.erase(currentChromosomeSNPIndex);
			}
			currentChromosomeSNPIndex += 1;
		}
		size_t j = currentChromosomeSNPIndex;
		size_t refpos = aln.Position;
		size_t readpos = 0;
		std::string barcode = readBarcode(aln);
		std::string umi = readUMI(aln);
		if (barcode == "") continue;
		if (umi == "") continue;
		readsWithTag += 1;
		for (const auto& cigar : aln.CigarData)
		{
			if (cigar.Type != 'M' && cigar.Type != '=' && cigar.Type != 'X')
			{
				if (cigar.Type == 'I' || cigar.Type == 'S')
				{
					readpos += cigar.Length;
				}
				if (cigar.Type == 'D' || cigar.Type == 'N')
				{
					refpos += cigar.Length;
				}
				continue;
			}
			assert(cigar.Type == 'M' || cigar.Type == '=' || cigar.Type == 'X');
			while (j < variants.size() && refpos > std::get<0>(variants[j])) j += 1;
			while (j < variants.size() && refpos+cigar.Length > std::get<0>(variants[j]))
			{
				std::string chromosome = currentChromosome;
				size_t variantPositionInRead = readpos+(std::get<0>(variants[j])-refpos);
				assert(variantPositionInRead < aln.QueryBases.size());
				char readNucleotide = aln.QueryBases[variantPositionInRead];
				if (readNucleotide == 'a') readNucleotide = 'A';
				if (readNucleotide == 'c') readNucleotide = 'C';
				if (readNucleotide == 'g') readNucleotide = 'G';
				if (readNucleotide == 't') readNucleotide = 'T';
				if (readNucleotide == std::get<1>(variants[j]))
				{
					matchCount += 1;
					if (aln.IsReverseStrand())
					{
						std::get<1>(matches[j][barcode][umi]) = true;
					}
					else
					{
						std::get<0>(matches[j][barcode][umi]) = true;
					}
				}
				if (readNucleotide == std::get<2>(variants[j]))
				{
					matchCount += 1;
					if (aln.IsReverseStrand())
					{
						std::get<3>(matches[j][barcode][umi]) = true;
					}
					else
					{
						std::get<2>(matches[j][barcode][umi]) = true;
					}
				}
				j += 1;
			}
			readpos += cigar.Length;
			refpos += cigar.Length;
		}
	}
	std::vector<SNPMatch> parsed;
	for (const auto& pospair : matches)
	{
		auto parts = parseSNPMatches(currentChromosome, refvariants.at(currentChromosome)[pospair.first], pospair.second);
		parsed.insert(parsed.end(), parts.begin(), parts.end());
	}
	sortSNPMatches(parsed);
	for (size_t i = 0; i < parsed.size(); i++)
	{
		callback(parsed[i]);
	}
	return std::make_tuple(reads, matchCount);
}

template <typename F>
std::tuple<size_t, size_t> streamSNPsFromBam(const std::string& bamFile, const std::unordered_map<std::string, std::vector<std::tuple<size_t, char, char>>>& refvariants, F callback)
{
	return streamSNPsFromBamBamtools(bamFile, refvariants, callback);
}

#endif

#endif
