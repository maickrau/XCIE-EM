#include <fstream>
#include <iostream>
#include <cxxopts.hpp>
#include "FileHandler.h"
#include "Logger.h"
#include "EM.h"
#include "Common.h"

int main(int argc, char** argv)
{
	cxxopts::Options options { "XCIE-EM" };
	options.add_options()
		("h,help", "Print help")
		("v,version", "Print version")
		("input-screadcounts", "Input scReadCounts cell/variant match table", cxxopts::value<std::string>())
		("input-preprocessed-table", "Input prerocessed table of cell/variant matches", cxxopts::value<std::string>())
		("o,output-prefix", "Output prefix", cxxopts::value<std::string>()->default_value("./result"))
		("force-phase", "File with pre-phased trio variants", cxxopts::value<std::string>())
		("EM-noise-magnitude", "initial EM noise magnitude", cxxopts::value<double>()->default_value("20"))
		("EM-noise-decay", "EM noise decay", cxxopts::value<double>()->default_value("0.95"))
		("EM-random-seed", "Random seed for EM initialization", cxxopts::value<size_t>()->default_value("1"))
		("EM-num-runs", "Number of runs for EM", cxxopts::value<size_t>()->default_value("10"))
		("exclude-region", "Exclude regions", cxxopts::value<std::vector<std::string>>())
		("exclude-PAR", "Exclude the PAR region. Equivalent to \"--exclude-region chrX:0-3500000\"")
		("exclude-XIST-grch38", "Exclude the XIST and TSIX genes in grch38 coordinates. Equivalent to \"--exclude-region chrX:73792204-73852753\"")
		("exclude-XIST-grch37", "Exclude the XIST and TSIX genes in grch37 coordinates. Equivalent to \"--exclude-region chrX:73012040-73072588\"")
		("exclude-XIST-chm13", "Exclude the XIST and TSIX genes in chm13 coordinates. Equivalent to \"--exclude-region chrX:72225527-72286069\"")
		("verbose", "Print more information while running")
	;
	cxxopts::ParseResult params;
	try
	{
		params = options.parse(argc, argv);
	}
	catch (const cxxopts::exceptions::no_such_option& e)
	{
		std::cerr << e.what() << std::endl;;
		std::abort();
	}
	if (params.count("v") == 1)
	{
		std::cerr << "Version: " << VERSION << std::endl;
		std::exit(0);
	}
	if (params.count("h") == 1)
	{
		std::cerr << options.help() << std::endl;
		std::exit(0);
	}
	bool paramError = false;
	if (params.count("input-preprocessed-table") == 0 && params.count("input-screadcounts") == 0)
	{
		std::cerr << "Input is required" << std::endl;
		paramError = true;
	}
	if (params.count("input-preprocessed-table") + params.count("input-screadcounts") > 1)
	{
		std::cerr << "Use only one input" << std::endl;
		paramError = true;
	}
	if (params["EM-noise-decay"].as<double>() >= 1.0)
	{
		std::cerr << "Noise decay must be less than 1" << std::endl;
		paramError = true;
	}
	if (params["EM-noise-magnitude"].as<double>() < 0)
	{
		std::cerr << "Noise magnitude must be 0 or positive" << std::endl;
		paramError = true;
	}
	if (params.count("exclude-region") > 0)
	{
		for (std::string value : params["exclude-region"].as<std::vector<std::string>>())
		{
			std::tuple<std::string, size_t, size_t> region = parseBedRegion(value);
			if (std::get<0>(region) == "" && std::get<1>(region) == std::numeric_limits<size_t>::max() && std::get<2>(region) == std::numeric_limits<size_t>::max())
			{
				std::cerr << "Could not parse region \"" << value << "\". Regions should be in format chrX:1-3000000" << std::endl;
				paramError = true;
			}
		}
	}
	if (paramError)
	{
		std::abort();
	}
	Logger::Log.setVerbosity(params.count("verbose"));
	std::string outputPrefix = params["o"].as<std::string>();
	std::string forcedPhaseFile = "";
	if (params.count("force-phase") > 0)
	{
		forcedPhaseFile = params["force-phase"].as<std::string>();
	}
	size_t randomSeed = params["EM-random-seed"].as<size_t>();
	size_t numTries = params["EM-num-runs"].as<size_t>();
//	std::string annotationFile { argv[5] };
//	std::string secondStepGeneList { argv[6] };
	double initialNoiseMagnitude = params["EM-noise-magnitude"].as<double>();
	double noiseDecay = params["EM-noise-decay"].as<double>();
	std::vector<CellMatch> cellMatches;
	std::vector<std::pair<size_t, size_t>> excludedRegions;
	if (params.count("exclude-region") > 0)
	{
		for (std::string value : params["exclude-region"].as<std::vector<std::string>>())
		{
			std::tuple<std::string, size_t, size_t> region = parseBedRegion(value);
			std::string chromosome = std::get<0>(region);
			if (chromosome != "23" && lowercase(chromosome) != "x" && lowercase(chromosome) != "chrx") continue;
			excludedRegions.emplace_back(std::get<1>(region), std::get<2>(region));
		}
	}
	if (params.count("exclude-PAR"))
	{
		excludedRegions.emplace_back(0, 3500000);
	}
	if (params.count("exclude-XIST-grch38"))
	{
		excludedRegions.emplace_back(73792204, 73852753);
	}
	if (params.count("exclude-XIST-grch37"))
	{
		excludedRegions.emplace_back(73012040, 73072588);
	}
	if (params.count("exclude-XIST-chm13"))
	{
		excludedRegions.emplace_back(72225527, 72286069);
	}
	std::sort(excludedRegions.begin(), excludedRegions.end());
	if (params.count("input-preprocessed-table") > 0)
	{
		std::string matchTableFile = params["input-preprocessed-table"].as<std::string>();
		Logger::Log.log(Logger::LogLevel::DebugInfo) << "read preprocessed counts from file " << matchTableFile << std::endl;
		cellMatches = readMatchCounts(matchTableFile);
		Logger::Log.log(Logger::LogLevel::DebugInfo) << cellMatches.size() << " count items" << std::endl;
	}
	if (params.count("input-screadcounts") > 0)
	{
		std::string scReadCountsFile = params["input-screadcounts"].as<std::string>();
		Logger::Log.log(Logger::LogLevel::DebugInfo) << "read screadcounts from file " << scReadCountsFile << std::endl;
		cellMatches = readScReadCountsMatchCounts(scReadCountsFile);
		Logger::Log.log(Logger::LogLevel::DebugInfo) << cellMatches.size() << " count items" << std::endl;
		Logger::Log.log(Logger::LogLevel::DebugInfo) << "filter out homozygous sites" << std::endl;
		cellMatches = filterOutHomozygousSites(cellMatches);
		Logger::Log.log(Logger::LogLevel::DebugInfo) << cellMatches.size() << " count items" << std::endl;
	}
	if (excludedRegions.size() > 0)
	{
		Logger::Log.log(Logger::LogLevel::Always) << "excluding regions:";
		for (auto t : excludedRegions)
		{
			Logger::Log.log(Logger::LogLevel::Always) << " " << t.first << "-" << t.second;
		}
		Logger::Log.log(Logger::LogLevel::Always) << std::endl;
		cellMatches = excludeRegions(cellMatches, excludedRegions);
	}
	writeCellMatchCounts(cellMatches, outputPrefix + ".preprocessed_matches.tsv");
	Logger::Log.log(Logger::LogLevel::DebugInfo) << "read forced variant phases" << std::endl;
	std::unordered_map<std::string, bool> forcedPhases;
	bool forcedPhasesAreMatPat;
	std::tie(forcedPhases, forcedPhasesAreMatPat) = readForcedVariantPhases(forcedPhaseFile);
	bool phasesAreMatPat = (forcedPhases.size() > 0) && forcedPhasesAreMatPat;
	EMOutput output = runEM(cellMatches, forcedPhases, randomSeed, initialNoiseMagnitude, noiseDecay, numTries);
	{
		std::ofstream variantResult { outputPrefix + ".variants.tsv" };
		writeResultVariants(output, phasesAreMatPat, variantResult);
	}
	{
		std::ofstream cellResultWithEscapeVariants { outputPrefix + ".cells.withescapevariants.tsv" };
		writeResultCells(output, phasesAreMatPat, cellResultWithEscapeVariants);
	}
	{
		std::ofstream cellResult { outputPrefix + ".cells.tsv" };
		writeResultCellOnlyNonescapeVariants(output, phasesAreMatPat, cellResult);
	}
	auto pseudobulk2 = getVariantPseudobulk(output, cellMatches, 2);
	writePseudobulkVariantResults(pseudobulk2, phasesAreMatPat, outputPrefix + ".pseudobulk.variants.confidence2.tsv");
	auto pseudobulk0 = getVariantPseudobulk(output, cellMatches, 0);
	writePseudobulkVariantResults(pseudobulk0, phasesAreMatPat, outputPrefix + ".pseudobulk.variants.confidence0.tsv");
}
