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
		std::cout << match.chromosome << "\t" << match.position << "\t" << match.ref << "\t" << match.alt << "\t" << match.barcode << "\t" << match.altFwCount << "\t" << match.altBwCount << "\t" << match.refFwCount << "\t" << match.refBwCount << "\t" << (match.altFwCount + match.altBwCount) << "\t" << (match.refFwCount + match.refBwCount) << std::endl;
	}
	Logger::Log.log(Logger::LogLevel::DebugInfo) << "done" << std::endl;
}
