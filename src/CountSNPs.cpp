#include <iostream>
#include "SNPCounter.h"
#include "Logger.h"

int main(int argc, char** argv)
{
	std::string vcfFile { argv[1] };
	std::string bamFile { argv[2] };
	Logger::Log.setVerbosity(1);
	std::cout << "CHROM\tPOS\tREF\tALT\tReadGroup\tSNVCountForward\tSNVCountReverse\tRefCountForward\tRefCountReverse\tSNVCount\tRefCount" << std::endl;
	auto matches = streamSNPsFromBam(bamFile, vcfFile, [](const auto& match)
	{
		printf("%s\t%lu\t%c\t%c\t%s\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\n", match.chromosome.c_str(), match.position, match.ref, match.alt, match.barcode.c_str(), match.altFwCount, match.altBwCount, match.refFwCount, match.refBwCount, match.altFwCount+match.altBwCount, match.refFwCount+match.refBwCount);
	});
	Logger::Log.log(Logger::LogLevel::DebugInfo) << std::get<0>(matches) << " reads" << std::endl;
	Logger::Log.log(Logger::LogLevel::DebugInfo) << std::get<1>(matches) << " read-variant matches" << std::endl;
	Logger::Log.log(Logger::LogLevel::DebugInfo) << "done" << std::endl;
}
