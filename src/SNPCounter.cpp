#include <limits>
#include <vector>
#include <tuple>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <api/BamReader.h>
#include <zstr.hpp> //https://github.com/mateidavid/zstr
#include "Common.h"
#include "SNPCounter.h"
#include "Logger.h"

std::unordered_map<std::string, std::vector<std::tuple<size_t, char, char>>> readVariants(std::istream& file)
{
	size_t chromColumnIndex = 0;
	size_t posColumnIndex = 1;
	size_t refColumnIndex = 3;
	size_t altColumnIndex = 4;
	std::unordered_map<std::string, std::vector<std::tuple<size_t, char, char>>> result;
	while (file.good())
	{
		std::string line;
		std::getline(file, line);
		if (line.size() < 5) continue;
		auto parts = split(line, '\t');
		if (parts[0] == "#CHROM")
		{
			if (parts[1] != "POS")
			{
				std::cerr << "VCF format wrong. Expected column 2 to be position." << std::endl;
				std::abort();
			}
			if (parts[3] != "REF")
			{
				std::cerr << "VCF format wrong. Expected column 4 to be ref allele." << std::endl;
				std::abort();
			}
			if (parts[4] != "ALT")
			{
				std::cerr << "VCF format wrong. Expected column 5 to be alt allele." << std::endl;
				std::abort();
			}
		}
		if (line[0] == '#') continue;
		if (parts.size() < 5) continue;
		std::string chromosome = parts[chromColumnIndex];
		size_t pos = std::stoull(parts[posColumnIndex])-1;
		std::string ref = parts[refColumnIndex];
		std::string alt = parts[altColumnIndex];
		// skip non-SNPs
		if (ref.size() != 1) continue;
		if (alt.size() != 1) continue;
		result[chromosome].emplace_back(pos, ref[0], alt[0]);
	}
	return result;
}

std::unordered_map<std::string, std::vector<std::tuple<size_t, char, char>>> readVariantsVcf(const std::string& vcfFile)
{
	std::ifstream file { vcfFile };
	return readVariants(file);
}

std::unordered_map<std::string, std::vector<std::tuple<size_t, char, char>>> readVariantsVcfGz(const std::string& vcfGzFile)
{
	zstr::ifstream file { vcfGzFile };
	return readVariants(file);
}

std::unordered_map<std::string, std::vector<std::tuple<size_t, char, char>>> readVariantsVcfParseFilename(const std::string& vcfFile)
{
	Logger::Log.log(Logger::LogLevel::DebugInfo) << "read variants from " << vcfFile << std::endl;
	std::unordered_map<std::string, std::vector<std::tuple<size_t, char, char>>> refVariants;
	if (hasExtension(vcfFile, ".vcf"))
	{
		refVariants = readVariantsVcf(vcfFile);
	}
	else if (hasExtension(vcfFile, ".vcf.gz"))
	{
		refVariants = readVariantsVcfGz(vcfFile);
	}
	else
	{
		std::cerr << "Unknown file extension for vcf file. Valid extensions are .vcf and .vcf.gz" << std::endl;
		std::abort();
	}
	Logger::Log.log(Logger::LogLevel::DebugInfo) << "variants in " << refVariants.size() << " chromosomes" << std::endl;
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "sort variants" << std::endl;
	size_t totalVariants = 0;
	for (auto& pair : refVariants)
	{
		std::sort(pair.second.begin(), pair.second.end());
		totalVariants += pair.second.size();
	}
	Logger::Log.log(Logger::LogLevel::DebugInfo) << totalVariants << " variants" << std::endl;
	return refVariants;
}

std::string readBarcode(const BamTools::BamAlignment& aln)
{
	std::string result;
	aln.GetTag("CB", result);
	if (result == "" || result == "-") aln.GetTag("CR", result);
	if (result == "-") result = ""
	return result;
}

std::vector<SNPMatch> countSNPsFromBamVcf(const std::string& vcfFile, const std::string& bamFile)
{
	std::vector<SNPMatch> result;
	std::tuple<size_t, size_t> counts;
	counts = streamSNPsFromBam(bamFile, vcfFile, [&result](SNPMatch& match)
	{
		result.emplace_back();
		std::swap(result.back(), match);
	});
	Logger::Log.log(Logger::LogLevel::DebugInfo) << std::get<0>(counts) << " reads" << std::endl;
	Logger::Log.log(Logger::LogLevel::DebugInfo) << std::get<1>(counts) << " read-variant matches" << std::endl;
	return result;
}

std::vector<CellMatch> getSNPMatchesFromBamVcf(const std::string& bamFile, const std::string& vcfFile)
{
	std::vector<CellMatch> result;
	streamSNPsFromBam(bamFile, vcfFile, [&result](const SNPMatch& item)
	{
		if (item.chromosome != "23" && lowercase(item.chromosome) != "x" && lowercase(item.chromosome) != "chrx")
		{
			return;
		}
		if (item.refFwCount + item.refBwCount > 0)
		{
			result.emplace_back();
			result.back().variant = item.chromosome + ":" + std::to_string(item.position) + ":" + item.ref + ":" + item.alt;
			result.back().cell = item.barcode;
			result.back().alt = false;
			result.back().count = item.refFwCount + item.refBwCount;
		}
		if (item.altFwCount + item.altBwCount > 0)
		{
			result.emplace_back();
			result.back().variant = item.chromosome + ":" + std::to_string(item.position) + ":" + item.ref + ":" + item.alt;
			result.back().cell = item.barcode;
			result.back().alt = true;
			result.back().count = item.altFwCount + item.altBwCount;
		}
	});
	return result;
}

void sortSNPMatches(std::vector<SNPMatch>& parsed)
{
	std::sort(parsed.begin(), parsed.end(), [](const auto& left, const auto& right)
	{
		if (left.chromosome < right.chromosome) return true;
		if (left.chromosome > right.chromosome) return false;
		if (left.position < right.position) return true;
		if (left.position > right.position) return false;
		if (left.barcode < right.barcode) return true;
		if (left.barcode > right.barcode) return false;
		assert(left.ref == right.ref);
		if (left.alt < right.alt) return true;
		if (left.alt > right.alt) return false;
		return false;
	});
}

std::vector<SNPMatch> parseSNPMatches(const std::string& currentChromosome, const std::tuple<size_t, char, char>& currentVariant, const std::unordered_map<std::string, std::tuple<size_t, size_t, size_t, size_t>>& matches)
{
	std::vector<SNPMatch> parsed;
	for (const auto& cellpair : matches)
	{
		parsed.emplace_back();
		parsed.back().chromosome = currentChromosome;
		if (currentChromosome.size() >= 4 && currentChromosome.substr(0, 3) == "chr") parsed.back().chromosome = currentChromosome.substr(3);
		parsed.back().position = std::get<0>(currentVariant)+1;
		parsed.back().ref = std::get<1>(currentVariant);
		parsed.back().alt = std::get<2>(currentVariant);
		parsed.back().barcode = cellpair.first;
		parsed.back().refFwCount = std::get<0>(cellpair.second);
		parsed.back().refBwCount = std::get<1>(cellpair.second);
		parsed.back().altFwCount = std::get<2>(cellpair.second);
		parsed.back().altBwCount = std::get<3>(cellpair.second);
	}
	return parsed;
}
