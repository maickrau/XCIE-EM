#include <algorithm>
#include <sstream>
#include <tuple>
#include <cassert>
#include <fstream>
#include <zstr.hpp>
#include "Common.h"
#include "Logger.h"
#include "FileHandler.h"

std::vector<CellMatch> readScReadCountsMatchCounts(const std::string& scReadCountsFile)
{
	std::ifstream file { scReadCountsFile };
	if (!file.good())
	{
		std::cerr << "Input scReadCounts file can't be read!" << std::endl;
		std::abort();
	}
	std::string header;
	std::getline(file, header);
	if (!file.good())
	{
		std::cerr << "Input scReadCounts file can't be read!" << std::endl;
		std::abort();
	}
	auto parts = split(header, '\t');
	if (parts.size() != 14 || parts[0] != "CHROM" || parts[1] != "POS" || parts[2] != "REF" || parts[3] != "ALT" || parts[4] != "ReadGroup" || parts[9] != "SNVCount" || parts[10] != "RefCount")
	{
		std::cerr << "Input scReadCounts file format is wrong. Double check that you are using the sparse matrix tsv file (output.tsv) and not the dense matrix file (output.cnt.matrix.tsv or output.vaf-m1.matrix.tsv)" << std::endl;
		std::abort();
	}
	std::vector<CellMatch> matches;
	while (file.good())
	{
		std::string line;
		std::getline(file, line);
		if (!file.good()) break;
		parts = split(line, '\t');
		if (parts[0] != "23" && lowercase(parts[0]) != "x" && lowercase(parts[0]) != "chrx")
		{
			continue;
		}
		std::string variantName = parts[0] + ":" + parts[1] + ":" + parts[2] + ":" + parts[3];
		std::string barcode = parts[4];
		size_t refMatch = std::stoull(parts[10]);
		size_t altMatch = std::stoull(parts[9]);
		if (refMatch > 0)
		{
			CellMatch match;
			match.variant = variantName;
			match.cell = barcode;
			match.alt = false;
			match.count = refMatch;
			matches.emplace_back(match);
		}
		if (altMatch > 0)
		{
			CellMatch match;
			match.variant = variantName;
			match.cell = barcode;
			match.alt = true;
			match.count = altMatch;
			matches.emplace_back(match);
		}
	}
	return matches;
}

std::vector<CellMatch> readMatchCounts(const std::string& matchTableFile)
{
	std::ifstream file { matchTableFile };
	std::vector<CellMatch> result;
	while (file.good())
	{
		std::string line;
		getline(file, line);
		if (!file.good()) break;
		if (line.size() == 0) continue;
		auto parts = split(line, '\t');
		if (parts.size() != 4 || (parts[2] != "REF" && parts[2] != "ALT"))
		{
			std::cerr << "Input table has invalid format" << std::endl;
			std::abort();
		}
		CellMatch match;
		match.cell = parts[0];
		match.variant = parts[1];
		match.count = std::stoull(parts[3]);
		if (parts[2] == "ALT")
		{
			match.alt = true;
		}
		else
		{
			assert(parts[2] == "REF");
			match.alt = false;
		}
		result.emplace_back(match);
	}
	return result;
}

void writeResultCellOnlyNonescapeVariants(const EMOutput& result, const bool phasesAreMatPat, std::ostream& stream)
{
	const std::string matName = matHapName(phasesAreMatPat);
	const std::string patName = patHapName(phasesAreMatPat);
	std::vector<std::string> cellOrder = getCellOrder(result.helpers.cellNameToIndex);
	stream << "cell\tcoverage\tactive_chrX\tactive_chrX_confidence\tescape_estimate\tescape_ci_low\tescape_ci_high" << "\n";
	for (const std::string& cell : cellOrder)
	{
		const size_t cellIndex = result.helpers.cellNameToIndex.at(cell);
		const bool matActive = result.result.cellIsMatActive[cellIndex];
		const double escape = result.result.cellEscapeFraction[cellIndex];
		const double scoreDifference = result.additions.cellActiveConfidenceOnlyNonescape[cellIndex];
		const double escapeConfidenceIntervalMin = result.additions.cellEscapeCILow[cellIndex];
		const double escapeConfidenceIntervalMax = result.additions.cellEscapeCIHigh[cellIndex];
		stream << cell << "\t" << result.helpers.cellCoverage[cellIndex] << "\t" << (matActive ? matName : patName) << "\t" << scoreDifference << "\t" << escape << "\t" << escapeConfidenceIntervalMin << "\t" << escapeConfidenceIntervalMax << "\n";
	}
}

void writeResultCells(const EMOutput& result, const bool phasesAreMatPat, std::ostream& stream)
{
	const std::string matName = matHapName(phasesAreMatPat);
	const std::string patName = patHapName(phasesAreMatPat);
	std::vector<std::string> cellOrder = getCellOrder(result.helpers.cellNameToIndex);
	stream << "cell\tcoverage\tactive_chrX\tactive_chrX_confidence\tescape_estimate\tescape_ci_low\tescape_ci_high" << "\n";
	for (const std::string& cell : cellOrder)
	{
		const size_t cellIndex = result.helpers.cellNameToIndex.at(cell);
		const bool matActive = result.result.cellIsMatActive[cellIndex];
		const double escape = result.result.cellEscapeFraction[cellIndex];
		const double scoreDifference = result.additions.cellActiveConfidence[cellIndex];
		const double escapeConfidenceIntervalMin = result.additions.cellEscapeCILow[cellIndex];
		const double escapeConfidenceIntervalMax = result.additions.cellEscapeCIHigh[cellIndex];
		stream << cell << "\t" << result.helpers.cellCoverage[cellIndex] << "\t" << (matActive ? matName : patName) << "\t" << scoreDifference << "\t" << escape << "\t" << escapeConfidenceIntervalMin << "\t" << escapeConfidenceIntervalMax << "\n";
	}
}

void writeResultVariants(const EMOutput& result, const bool phasesAreMatPat, std::ostream& stream)
{
	const std::string matName = matHapName(phasesAreMatPat);
	const std::string patName = patHapName(phasesAreMatPat);
	std::vector<std::string> variantOrder = getVariantOrder(result.helpers.variantNameToIndex);
	stream << "variant\tcoverage\tphase\tphase_confidence\tescape_estimate\tescape_ci_low\tescape_ci_high" << "\n";
	for (const std::string& variant : variantOrder)
	{
		const size_t variantIndex = result.helpers.variantNameToIndex.at(variant);
		const bool matRef = result.result.variantIsMatRef[variantIndex];
		const double escape = result.result.variantEscapeFraction[variantIndex];
		const double phaseScoreDifference = result.additions.variantPhaseConfidence[variantIndex];
		const double escapeConfidenceIntervalMin = result.additions.variantEscapeCILow[variantIndex];
		const double escapeConfidenceIntervalMax = result.additions.variantEscapeCIHigh[variantIndex];
		stream << variant << "\t" << result.helpers.variantCoverage[variantIndex] << "\t" <<  (matRef ? matName : patName) << "\t" << phaseScoreDifference << "\t" << escape << "\t" << escapeConfidenceIntervalMin << "\t" << escapeConfidenceIntervalMax << "\n";
	}
}

void writePseudobulkResults(const std::vector<PseudobulkInfo>& pseudobulk, const bool phasesAreMatPat, const std::string& firstColumnName, const std::string filename)
{
	const std::string matName = matHapName(phasesAreMatPat);
	const std::string patName = patHapName(phasesAreMatPat);
	std::ofstream file { filename };
	file << firstColumnName << "\t" << matName << "_active_expression" << "\t" << matName << "_inactive_expression" << "\t" << patName << "_active_expression" << "\t" << patName << "_inactive_expression" << "\t" << "Xi" << "\t" << "pvalue_Xiover10_unadjusted" << "\t" << "pvalue_Xiunder10_unadjusted" << "\n";
	for (const auto& t : pseudobulk)
	{
		double pValueOver = getBinomialPValueGreaterThan(0.1, t.matXi+t.patXi, t.matXa+t.matXi+t.patXa+t.patXi);
		double pValueUnder = getBinomialPValueLessThan(0.1, t.matXi+t.patXi, t.matXa+t.matXi+t.patXa+t.patXi);
		double Xi = (double)(t.matXi+t.patXi) / (double)(t.matXi+t.matXa+t.patXi+t.patXa);
		if (t.matXi == 0 && t.matXa == 0 && t.patXi == 0 && t.patXa == 0) Xi = 0;
		file << t.name << "\t" << t.matXa << "\t" << t.matXi << "\t" << t.patXa << "\t" << t.patXi << "\t" << Xi << "\t" << pValueOver << "\t" << pValueUnder << "\n";
	}
}

std::pair<std::unordered_map<std::string, bool>, bool> readForcedVariantPhases(const std::string& filename)
{
	std::unordered_map<std::string, bool> forcedVariants;
	bool forcedPhasesAreMatPat = false;
	bool foundAny = false;
	if (filename == "") return std::make_pair(forcedVariants, forcedPhasesAreMatPat);
	std::ifstream file { filename };
	size_t totalForcedVariants = 0;
	while (file.good())
	{
		std::string line;
		getline(file, line);
		if (!file.good()) break;
		if (line.size() < 3) continue;
		std::stringstream sstr { line };
		std::string variant;
		std::string origin;
		sstr >> variant >> origin;
		totalForcedVariants += 1;
		if (origin == "hap1" || origin == "mat")
		{
			forcedVariants[variant] = true;
		}
		if (origin == "hap2" || origin == "pat")
		{
			forcedVariants[variant] = false;
		}
		if (origin == "hap1" || origin == "hap2")
		{
			if (foundAny && forcedPhasesAreMatPat)
			{
				std::cerr << "Forced phase file has invalid format: mixing hap1/hap2 and mat/pat. Phases should be either hap1/hap2 or mat/pat but not both." << "\n";
				std::abort();
			}
			foundAny = true;
			forcedPhasesAreMatPat = false;
		}
		if (origin == "mat" || origin == "pat")
		{
			if (foundAny && !forcedPhasesAreMatPat)
			{
				std::cerr << "Forced phase file has invalid format: mixing hap1/hap2 and mat/pat. Phases should be either hap1/hap2 or mat/pat but not both." << "\n";
				std::abort();
			}
			foundAny = true;
			forcedPhasesAreMatPat = true;
		}
		if (origin != "hap1" && origin != "hap2" && origin != "mat" && origin != "pat")
		{
			std::cerr << "Invalid phase in forced phase file: \"" << origin << "\"" << "\n";
			std::abort();
		}
	}
	Logger::Log.log(Logger::LogLevel::DebugInfo) << "total forced variants " << totalForcedVariants << "\n";
	return std::make_pair(forcedVariants, forcedPhasesAreMatPat);
}

void writeCellMatchCounts(const std::vector<CellMatch>& cellMatches, const std::string& filename)
{
	std::ofstream file { filename };
	for (const auto& match : cellMatches)
	{
		file << match.cell << "\t" << match.variant << "\t" << (match.alt ? "ALT" : "REF") << "\t" << match.count << "\n";
	}
}

std::unordered_map<std::string, std::string> parseTagsGff3(const std::string& tagstr)
{
	auto parts = split(tagstr, ';');
	std::unordered_map<std::string, std::string> result;
	for (const std::string& tag : parts)
	{
		auto parts2 = split(tag, '=');
		if (parts2.size() != 2) continue;
		result[parts2[0]] = parts2[1];
	}
	return result;
}

std::vector<std::tuple<size_t, size_t, std::string, std::string>> getGeneInfoFromStream(std::istream& stream, const bool onlyProteinCoding)
{
	std::vector<std::tuple<size_t, size_t, std::string, std::string>> result;
	while (stream.good())
	{
		std::string line;
		std::getline(stream, line);
		if (!stream.good()) break;
		if (line.size() < 10) continue;
		if (line[0] == '#') continue;
		auto parts = split(line, '\t');
		std::string chromosome = parts[0];
		if (chromosome != "23" && lowercase(chromosome) != "x" && lowercase(chromosome) != "chrx") continue;
		std::string type = parts[2];
		if (type != "gene") continue;
		auto tags = parseTagsGff3(parts[8]);
		if (onlyProteinCoding)
		{
			if (tags.count("gene_type") == 0 || tags.at("gene_type") != "protein_coding") continue;
		}
		std::string geneId = tags.count("gene_id") == 1 ? tags.at("gene_id") : "";
		std::string geneName = tags.count("gene_id") == 1 ? tags.at("gene_name") : "";
		size_t startPos = std::stoull(parts[3]);
		size_t endPos = std::stoull(parts[4]);
		if (geneName == "") geneName = geneId;
		result.emplace_back(startPos, endPos, geneId, geneName);
	}
	std::sort(result.begin(), result.end());
	return result;
}

std::vector<std::tuple<size_t, size_t, std::string, std::string>> getGeneInfoGff3(const std::string gff3Path, const bool onlyProteinCoding)
{
	std::ifstream file { gff3Path };
	return getGeneInfoFromStream(file, onlyProteinCoding);
}

std::vector<std::tuple<size_t, size_t, std::string, std::string>> getGeneInfoGff3Gz(const std::string gff3Path, const bool onlyProteinCoding)
{
	zstr::ifstream file { gff3Path };
	return getGeneInfoFromStream(file, onlyProteinCoding);
}

std::vector<std::tuple<size_t, size_t, std::string, std::string>> getGeneInfo(const std::string gff3Path, const bool onlyProteinCoding)
{
	if (hasExtension(gff3Path, ".gff3"))
	{
		return getGeneInfoGff3(gff3Path, onlyProteinCoding);
	}
	else if (hasExtension(gff3Path, ".gff3.gz"))
	{
		return getGeneInfoGff3Gz(gff3Path, onlyProteinCoding);
	}
	std::cerr << "Unknown annotation extension. Valid extensions are .gff3 and .gff3.gz" << std::endl;
	std::abort();
}

void writeGenesPerVariant(const EMOutput& output, const std::vector<std::vector<std::pair<std::string, std::string>>>& variantToGeneMatch, const std::string& filename)
{
	std::ofstream file { filename };
	file << "variant\tcount_genes\tgene_ids\tgene_names" << "\n";
	std::vector<std::string> variantOrder = getVariantOrder(output.helpers.variantNameToIndex);
	for (const std::string& variantName : variantOrder)
	{
		const size_t index = output.helpers.variantNameToIndex.at(variantName);
		std::vector<std::string> geneIds;
		std::vector<std::string> geneNames;
		for (const auto& pair : variantToGeneMatch[index])
		{
			geneIds.emplace_back(pair.first);
			geneNames.emplace_back(pair.second);
		}
		file << variantName << "\t" << variantToGeneMatch[index].size() << "\t" << join(',', geneIds) << "\t" << join(',', geneNames) << "\n";
	}
}

std::unordered_map<std::string, std::string> readCellGrouping(const std::string& filename)
{
	std::ifstream file { filename };
	std::unordered_map<std::string, std::string> result;
	{
		std::string header;
		std::getline(file, header);
		auto parts = split(header, '\t');
		if (parts.size() != 2 || (lowercase(parts[0]) != "barcode" && lowercase(parts[0]) != "cell"))
		{
			std::cerr << "Invalid cell group file format. Format should have two columns with barcode on first column." << std::endl;
			std::abort();
		}
	}
	while (file.good())
	{
		std::string line;
		std::getline(file, line);
		if (!file.good()) break;
		auto parts = split(line, '\t');
		if (result.count(parts[0]) == 1)
		{
			std::cerr << "Cell group file has duplicates. Duplicate barcode: " << parts[0] << std::endl;
			std::abort();
		}
		result[parts[0]] = parts[1];
	}
	return result;
}

void writeCellGroupStatistics(const std::vector<CellMatch>& cellMatches, const std::unordered_map<std::string, std::string>& cellGrouping, const std::string& filename)
{
	std::unordered_map<std::string, size_t> ASEperGroup;
	size_t ASEtotal = 0;
	std::unordered_map<std::string, std::unordered_set<std::string>> cellsPerGroup;
	std::unordered_set<std::string> cellsTotal;
	size_t ASENoGroup = 0;
	std::unordered_set<std::string> cellsNoGroup;
	for (const auto& t : cellMatches)
	{
		ASEtotal += t.count;
		cellsTotal.emplace(t.cell);
		if (cellGrouping.count(t.cell) == 1)
		{
			std::string group = cellGrouping.at(t.cell);
			ASEperGroup[group] += t.count;
			cellsPerGroup[group].emplace(t.cell);
		}
		else
		{
			ASENoGroup += t.count;
			cellsNoGroup.emplace(t.cell);
		}
	}
	std::vector<std::string> groupNames;
	for (const auto& pair : cellsPerGroup)
	{
		groupNames.emplace_back(pair.first);
	}
	std::sort(groupNames.begin(), groupNames.end());
	std::ofstream file { filename };
	file << "Group\tASE\tASE_cells\tASE_percent\tASE_cells_percent\n";
	for (const std::string& group : groupNames)
	{
		assert(ASEperGroup.count(group) == 1);
		file << group << "\t" << ASEperGroup.at(group) << "\t" << cellsPerGroup.at(group).size() << "\t" << (double)ASEperGroup.at(group)/ASEtotal*100 << "\t" << (double)cellsPerGroup.at(group).size()/cellsTotal.size()*100 << "\n";
	}
	file << "No group\t" << ASENoGroup << "\t" << cellsNoGroup.size() << "\t" << (double)ASENoGroup/ASEtotal*100 << "\t" << (double)cellsNoGroup.size()/cellsTotal.size()*100 << "\n";
	file << "Total\t" << ASEtotal << "\t" << cellsTotal.size() << "\t100\t100" << "\n";
}
