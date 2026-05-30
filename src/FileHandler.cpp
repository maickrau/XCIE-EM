#include <sstream>
#include <tuple>
#include <cassert>
#include <fstream>
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

void writePseudobulkResults(const std::vector<PseudobulkInfo>& pseudobulk, const bool phasesAreMatPat, const std::string filename)
{
	const std::string matName = matHapName(phasesAreMatPat);
	const std::string patName = patHapName(phasesAreMatPat);
	std::ofstream file { filename };
	file << "variant\t" << matName << "_active_expression" << "\t" << matName << "_inactive_expression" << "\t" << patName << "_active_expression" << "\t" << patName << "_inactive_expression" << "\t" << "Xi" << "\t" << "pvalue_Xiover10" << "\t" << "pvalue_Xiunder10" << "\n";
	for (const auto& t : pseudobulk)
	{
		double pValueOver = getBinomialPValueGreaterThan(0.1, t.matXi+t.patXi, t.matXa+t.matXi+t.patXa+t.patXi);
		double pValueUnder = getBinomialPValueLessThan(0.1, t.matXi+t.patXi, t.matXa+t.matXi+t.patXa+t.patXi);
		double Xi = (double)(t.matXi+t.patXi) / (double)(t.matXi+t.matXa+t.patXi+t.patXa);
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
