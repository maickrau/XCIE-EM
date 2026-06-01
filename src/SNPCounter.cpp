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

std::string readBarcode(const BamTools::BamAlignment& aln)
{
	std::string result;
	aln.GetTag("CB", result);
	return result;
}

std::vector<SNPMatch> countSNPsFromBam(const std::string& bamFile, const std::unordered_map<std::string, std::vector<std::tuple<size_t, char, char>>>& refvariants)
{
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
	std::unordered_map<std::string, std::unordered_map<size_t, std::unordered_map<std::string, std::tuple<size_t, size_t, size_t, size_t>>>> matches;
	size_t bamAlns = 0;
	size_t reads = 0;
	size_t readsWithTag = 0;
	size_t matchCount = 0;
	size_t lastReadStartPosition = 0;
	while (reader.GetNextAlignment(aln))
	{
		bamAlns += 1;
		if (!aln.IsPrimaryAlignment()) continue;
		if (aln.IsDuplicate()) continue;
		if (aln.IsFailedQC()) continue;
		if (!aln.IsMapped()) continue;
//		if (aln.MapQuality < 20) continue;
		if ((size_t)aln.RefID != refId)
		{
			refId = aln.RefID;
			currentChromosome = references[refId].RefName;
			Logger::Log.log(Logger::LogLevel::DebugInfo) << "begin chromosome " << currentChromosome << std::endl;
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
		while (currentChromosomeSNPIndex < variants.size() && (size_t)aln.Position > std::get<0>(variants[currentChromosomeSNPIndex])) currentChromosomeSNPIndex += 1;
		size_t j = currentChromosomeSNPIndex;
		size_t refpos = aln.Position;
		size_t readpos = 0;
		std::string barcode = readBarcode(aln);
		if (barcode == "") continue;
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
				if (readNucleotide == std::get<1>(variants[j]))
				{
					matchCount += 1;
					if (aln.IsReverseStrand())
					{
						std::get<1>(matches[chromosome][j][barcode]) += 1;
					}
					else
					{
						std::get<0>(matches[chromosome][j][barcode]) += 1;
					}
				}
				if (readNucleotide == std::get<2>(variants[j]))
				{
					matchCount += 1;
					if (aln.IsReverseStrand())
					{
						std::get<3>(matches[chromosome][j][barcode]) += 1;
					}
					else
					{
						std::get<2>(matches[chromosome][j][barcode]) += 1;
					}
				}
				j += 1;
			}
			readpos += cigar.Length;
			refpos += cigar.Length;
		}
	}
	Logger::Log.log(Logger::LogLevel::DebugInfo) << reads << " reads" << std::endl;
	Logger::Log.log(Logger::LogLevel::DebugInfo) << readsWithTag << " reads with barcode" << std::endl;
	Logger::Log.log(Logger::LogLevel::DebugInfo) << matchCount << " read-variant matches" << std::endl;
	std::vector<SNPMatch> parsed;
	for (const auto& variantpair : matches)
	{
		for (const auto& pospair : variantpair.second)
		{
			for (const auto& cellpair : pospair.second)
			{
				parsed.emplace_back();
				parsed.back().chromosome = variantpair.first;
				if (variantpair.first.size() >= 4 && variantpair.first.substr(0, 3) == "chr") parsed.back().chromosome = variantpair.first.substr(3);
				parsed.back().position = std::get<0>(refvariants.at(variantpair.first)[pospair.first])+1;
				parsed.back().ref = std::get<1>(refvariants.at(variantpair.first)[pospair.first]);
				parsed.back().alt = std::get<2>(refvariants.at(variantpair.first)[pospair.first]);
				parsed.back().barcode = cellpair.first;
				parsed.back().refFwCount = std::get<0>(cellpair.second);
				parsed.back().refBwCount = std::get<1>(cellpair.second);
				parsed.back().altFwCount = std::get<2>(cellpair.second);
				parsed.back().altBwCount = std::get<3>(cellpair.second);
			}
		}
	}
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
	return parsed;
}

std::vector<SNPMatch> countSNPsFromBamVcf(const std::string& vcfFile, const std::string& bamFile)
{
	Logger::Log.log(Logger::LogLevel::DebugInfo) << "read variants from " << vcfFile << std::endl;
	std::unordered_map<std::string, std::vector<std::tuple<size_t, char, char>>> variants;
	if (hasExtension(vcfFile, ".vcf"))
	{
		variants = readVariantsVcf(vcfFile);
	}
	else if (hasExtension(vcfFile, ".vcf.gz"))
	{
		variants = readVariantsVcfGz(vcfFile);
	}
	else
	{
		std::cerr << "Unknown file extension for vcf file. Valid extensions are .vcf and .vcf.gz" << std::endl;
		std::abort();
	}
	Logger::Log.log(Logger::LogLevel::DebugInfo) << "variants in " << variants.size() << " chromosomes" << std::endl;
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "sort variants" << std::endl;
	size_t totalVariants = 0;
	for (auto& pair : variants)
	{
		std::sort(pair.second.begin(), pair.second.end());
		totalVariants += pair.second.size();
	}
	Logger::Log.log(Logger::LogLevel::DebugInfo) << totalVariants << " variants" << std::endl;
	Logger::Log.log(Logger::LogLevel::DebugInfo) << "read bam from " << bamFile << std::endl;
	auto result = countSNPsFromBam(bamFile, variants);
	return result;
}

std::vector<CellMatch> convertSNPMatchesToCellMatches(const std::vector<SNPMatch>& snpMatches)
{
	std::vector<CellMatch> result;
	for (const auto& item : snpMatches)
	{
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
	}
	return result;
}
