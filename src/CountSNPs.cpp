#include <iostream>
#include "SNPCounter.h"
#include "Logger.h"

int main(int argc, char** argv)
{
	std::string vcfFile { argv[1] };
	std::string bamFile { argv[2] };
	Logger::Log.setVerbosity(1);
	auto matches = countSNPsFromBamVcf(vcfFile, bamFile);
	Logger::Log.log(Logger::LogLevel::DebugInfo) << "write results" << std::endl;
	std::cout << "CHROM\tPOS\tREF\tALT\tReadGroup\tSNVCountForward\tSNVCountReverse\tRefCountForward\tRefCountReverse\tSNVCount\tRefCount" << std::endl;
	for (const auto& match : matches)
	{
		printf("%s\t%lu\t%c\t%c\t%s\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\n", match.chromosome.c_str(), match.position, match.ref, match.alt, match.barcode.c_str(), match.altFwCount, match.altBwCount, match.refFwCount, match.refBwCount, match.altFwCount+match.altBwCount, match.refFwCount+match.refBwCount);
	}
	Logger::Log.log(Logger::LogLevel::DebugInfo) << "done" << std::endl;
}
