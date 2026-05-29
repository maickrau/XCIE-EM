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

std::vector<CellMatch> countSNPsFromBam(const std::string& bamFile, const std::unordered_map<std::string, std::vector<std::tuple<size_t, char, char>>>& refvariants)
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
	std::unordered_map<std::string, std::unordered_map<std::string, std::pair<size_t, size_t>>> matches;
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
				std::string variantName = currentChromosome.substr(3) + ":" + std::to_string(std::get<0>(variants[j])+1) + ":" + std::get<1>(variants[j]) + ":" + std::get<2>(variants[j]);
				size_t variantPositionInRead = readpos+(std::get<0>(variants[j])-refpos);
				assert(variantPositionInRead < aln.QueryBases.size());
				char readNucleotide = aln.QueryBases[variantPositionInRead];
				if (readNucleotide == std::get<1>(variants[j]))
				{
					matchCount += 1;
					matches[variantName][barcode].first += 1;
				}
				if (readNucleotide == std::get<2>(variants[j]))
				{
					matchCount += 1;
					matches[variantName][barcode].second += 1;
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
	std::vector<CellMatch> parsed;
	for (const auto& variantpair : matches)
	{
		for (const auto& cellpair : variantpair.second)
		{
			if (cellpair.second.first > 0)
			{
				parsed.emplace_back();
				parsed.back().cell = cellpair.first;
				parsed.back().variant = variantpair.first;
				parsed.back().alt = false;
				parsed.back().count = cellpair.second.first;
			}
			if (cellpair.second.second > 0)
			{
				parsed.emplace_back();
				parsed.back().cell = cellpair.first;
				parsed.back().variant = variantpair.first;
				parsed.back().alt = true;
				parsed.back().count = cellpair.second.second;
			}
		}
	}
	return parsed;
}

std::vector<CellMatch> countSNPsFromBamVcf(const std::string& vcfFile, const std::string& bamFile)
{
	Logger::Log.log(Logger::LogLevel::DebugInfo) << "read variants from " << vcfFile << std::endl;
	std::unordered_map<std::string, std::vector<std::tuple<size_t, char, char>>> variants;
	if (vcfFile.size() >= 4 && vcfFile.substr(vcfFile.size() - 4) == ".vcf")
	{
		variants = readVariantsVcf(vcfFile);
	}
	else if (vcfFile.size() >= 7 && vcfFile.substr(vcfFile.size() - 7) == ".vcf.gz")
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
	std::sort(result.begin(), result.end(), [](const auto& left, const auto& right)
	{
		size_t leftpos = parseVariantPosition(left.variant);
		size_t rightpos = parseVariantPosition(right.variant);
		if (leftpos < rightpos) return true;
		if (rightpos < leftpos) return false;
		if (left.cell < right.cell) return true;
		if (right.cell < left.cell) return false;
		if (!left.alt && right.alt) return true;
		return false;
	});
	return result;
}

