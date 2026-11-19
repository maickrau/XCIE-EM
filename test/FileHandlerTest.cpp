#include <sstream>
#include <algorithm>
#include <zstr.hpp>
#include "../src/FileHandler.h"
#include "../src/Common.h"
#include "InternalDeclarations.h"
#include "TestHelpers.h"
#include "OutputBuilder.h"
#include "doctest.h"

using TestHelpers::makeMatch;
using TestHelpers::TempFile;
using TestHelpers::abortsInChild;
using TestHelpers::OutputBuilder;

namespace
{

const std::string scReadCountsHeader = "CHROM\tPOS\tREF\tALT\tReadGroup\tSNVCountForward\tSNVCountReverse\tRefCountForward\tRefCountReverse\tSNVCount\tRefCount\tGoodReads\t%BadRead\tVAF\n";

std::string scReadCountsLine(const std::string& chrom, size_t pos, char ref, char alt, const std::string& barcode, size_t snvCount, size_t refCount)
{
	return chrom + "\t" + std::to_string(pos) + "\t" + ref + "\t" + alt + "\t" + barcode + "\t0\t0\t0\t0\t" + std::to_string(snvCount) + "\t" + std::to_string(refCount) + "\t" + std::to_string(snvCount+refCount) + "\t0\t0\n";
}

bool sameMatch(const CellMatch& left, const CellMatch& right)
{
	return left.cell == right.cell && left.variant == right.variant && left.alt == right.alt && left.count == right.count;
}

std::string gffLine(const std::string& chrom, const std::string& type, size_t start, size_t end, const std::string& tags)
{
	return chrom + "\tHAVANA\t" + type + "\t" + std::to_string(start) + "\t" + std::to_string(end) + "\t.\t+\t.\t" + tags + "\n";
}

}

TEST_SUITE("FileHandler")
{

TEST_CASE("readMatchCounts parses the preprocessed table")
{
	TempFile file { ".tsv", "cell1\tX:100:A:G\tREF\t5\ncell2\tX:200:C:T\tALT\t3\n\ncell1\tX:200:C:T\tREF\t1\n" };
	auto matches = readMatchCounts(file.path);
	REQUIRE(matches.size() == 3);
	CHECK(sameMatch(matches[0], makeMatch("cell1", "X:100:A:G", false, 5)));
	CHECK(sameMatch(matches[1], makeMatch("cell2", "X:200:C:T", true, 3)));
	CHECK(sameMatch(matches[2], makeMatch("cell1", "X:200:C:T", false, 1)));
}

TEST_CASE("readMatchCounts on an empty or missing file returns nothing")
{
	TempFile file { ".tsv", "" };
	CHECK(readMatchCounts(file.path).empty());
	CHECK(readMatchCounts("/nonexistent/path/to/file.tsv").empty());
}

TEST_CASE("readMatchCounts aborts on malformed input")
{
	TempFile wrongColumns { ".tsv", "cell1\tX:100:A:G\tREF\n" };
	TempFile wrongAllele { ".tsv", "cell1\tX:100:A:G\tref\t5\n" };
	TempFile extraColumn { ".tsv", "cell1\tX:100:A:G\tREF\t5\t1\n" };
	CHECK(abortsInChild([&]() { readMatchCounts(wrongColumns.path); }));
	CHECK(abortsInChild([&]() { readMatchCounts(wrongAllele.path); }));
	CHECK(abortsInChild([&]() { readMatchCounts(extraColumn.path); }));
}

TEST_CASE("readMatchCounts reads the last line even without a trailing newline")
{
	TempFile file { ".tsv", "cell1\tX:100:A:G\tREF\t5\ncell2\tX:200:C:T\tALT\t3" };
	CHECK(readMatchCounts(file.path).size() == 2);
}

TEST_CASE("writeCellMatchCounts and readMatchCounts roundtrip")
{
	std::vector<CellMatch> matches {
		makeMatch("AAAC-1", "X:100:A:G", false, 5),
		makeMatch("AAAC-1", "X:100:A:G", true, 1),
		makeMatch("bam.bam:0:GGGT-1", "X:155000000:C:T", true, 1000000),
	};
	TempFile file { ".tsv" };
	writeCellMatchCounts(matches, file.path);
	CHECK(file.read() == "AAAC-1\tX:100:A:G\tREF\t5\nAAAC-1\tX:100:A:G\tALT\t1\nbam.bam:0:GGGT-1\tX:155000000:C:T\tALT\t1000000\n");
	auto readBack = readMatchCounts(file.path);
	REQUIRE(readBack.size() == matches.size());
	for (size_t i = 0; i < matches.size(); i++)
	{
		CHECK(sameMatch(readBack[i], matches[i]));
	}
}

TEST_CASE("readScReadCountsMatchCounts")
{
	std::string contents = scReadCountsHeader;
	contents += scReadCountsLine("X", 100, 'A', 'G', "cell1", 2, 3);
	contents += scReadCountsLine("chrX", 200, 'C', 'T', "cell2", 0, 4);
	contents += scReadCountsLine("23", 300, 'G', 'A', "cell3", 6, 0);
	contents += scReadCountsLine("x", 400, 'T', 'C', "cell4", 1, 1);
	contents += scReadCountsLine("ChrX", 500, 'T', 'C', "cell4", 1, 0);
	contents += scReadCountsLine("1", 100, 'A', 'G', "cell1", 2, 3);
	contents += scReadCountsLine("chrY", 100, 'A', 'G', "cell1", 2, 3);
	contents += scReadCountsLine("X", 600, 'T', 'C', "cell5", 0, 0);
	TempFile file { ".tsv", contents };
	auto matches = readScReadCountsMatchCounts(file.path);
	REQUIRE(matches.size() == 7);
	CHECK(sameMatch(matches[0], makeMatch("cell1", "X:100:A:G", false, 3)));
	CHECK(sameMatch(matches[1], makeMatch("cell1", "X:100:A:G", true, 2)));
	// variant names for chrX spelled differently are checked in a separate test below
	CHECK(matches[2].cell == "cell2");
	CHECK(matches[2].alt == false);
	CHECK(matches[2].count == 4);
	CHECK(parseVariantPosition(matches[2].variant) == 200);
	CHECK(sameMatch(matches[3], makeMatch("cell3", "23:300:G:A", true, 6)));
	CHECK(matches[4].alt == false);
	CHECK(matches[4].count == 1);
	CHECK(matches[5].alt == true);
	CHECK(matches[5].variant == matches[4].variant);
	CHECK(parseVariantPosition(matches[5].variant) == 400);
	CHECK(matches[6].alt == true);
	CHECK(parseVariantPosition(matches[6].variant) == 500);
}

TEST_CASE("readScReadCountsMatchCounts names chrX variants like the BAM input and the forced phase format" * doctest::should_fail())
{
	// BUG: the variant name is built from the raw CHROM column, so scReadCounts input with
	// CHROM "chrX" gives names like "chrX:200:C:T". The BAM input strips "chr" ("X:200:C:T"),
	// and the README specifies "X:(position):(ref):(alt)" for --force-phase. With chrX-named
	// scReadCounts input every trio-phased variant is silently ignored.
	std::string contents = scReadCountsHeader;
	contents += scReadCountsLine("chrX", 200, 'C', 'T', "cell2", 1, 4);
	contents += scReadCountsLine("ChrX", 300, 'C', 'T', "cell2", 1, 4);
	TempFile file { ".tsv", contents };
	auto matches = readScReadCountsMatchCounts(file.path);
	REQUIRE(matches.size() == 4);
	CHECK(matches[0].variant == "X:200:C:T");
	CHECK(matches[2].variant == "X:300:C:T");
}

TEST_CASE("readScReadCountsMatchCounts with only a header")
{
	TempFile file { ".tsv", scReadCountsHeader };
	CHECK(readScReadCountsMatchCounts(file.path).empty());
}

TEST_CASE("readScReadCountsMatchCounts aborts on invalid input")
{
	TempFile empty { ".tsv", "" };
	TempFile wrongHeader { ".tsv", "CHROM\tPOS\tREF\tALT\n" };
	std::string denseHeader = scReadCountsHeader;
	denseHeader.replace(denseHeader.find("ReadGroup"), 9, "cell1");
	TempFile dense { ".tsv", denseHeader };
	CHECK(abortsInChild([&]() { readScReadCountsMatchCounts("/nonexistent/file.tsv"); }));
	CHECK(abortsInChild([&]() { readScReadCountsMatchCounts(empty.path); }));
	CHECK(abortsInChild([&]() { readScReadCountsMatchCounts(wrongHeader.path); }));
	CHECK(abortsInChild([&]() { readScReadCountsMatchCounts(dense.path); }));
}

TEST_CASE("readForcedVariantPhases")
{
	SUBCASE("empty filename")
	{
		auto result = readForcedVariantPhases("");
		CHECK(result.first.empty());
		CHECK(result.second == false);
	}
	SUBCASE("mat/pat phases")
	{
		TempFile file { ".tsv", "X:100:A:G\tmat\nX:200:C:T\tpat\n\nX:300:G:A mat\n" };
		auto result = readForcedVariantPhases(file.path);
		CHECK(result.second == true);
		REQUIRE(result.first.size() == 3);
		CHECK(result.first.at("X:100:A:G") == true);
		CHECK(result.first.at("X:200:C:T") == false);
		CHECK(result.first.at("X:300:G:A") == true);
	}
	SUBCASE("hap1/hap2 phases")
	{
		TempFile file { ".tsv", "X:100:A:G\thap2\nX:200:C:T\thap1\n" };
		auto result = readForcedVariantPhases(file.path);
		CHECK(result.second == false);
		REQUIRE(result.first.size() == 2);
		CHECK(result.first.at("X:100:A:G") == false);
		CHECK(result.first.at("X:200:C:T") == true);
	}
	SUBCASE("later lines override earlier lines for the same variant")
	{
		TempFile file { ".tsv", "X:100:A:G\tmat\nX:100:A:G\tpat\n" };
		auto result = readForcedVariantPhases(file.path);
		REQUIRE(result.first.size() == 1);
		CHECK(result.first.at("X:100:A:G") == false);
	}
	SUBCASE("non-chrX variants are kept in the map")
	{
		TempFile file { ".tsv", "1:100:A:G\tmat\n" };
		auto result = readForcedVariantPhases(file.path);
		CHECK(result.first.size() == 1);
	}
}

TEST_CASE("readForcedVariantPhases aborts on invalid input")
{
	TempFile mixed1 { ".tsv", "X:100:A:G\tmat\nX:200:C:T\thap1\n" };
	TempFile mixed2 { ".tsv", "X:100:A:G\thap2\nX:200:C:T\tpat\n" };
	TempFile invalidPhase { ".tsv", "X:100:A:G\tMAT\n" };
	TempFile missingPhase { ".tsv", "X:100:A:G\n" };
	CHECK(abortsInChild([&]() { readForcedVariantPhases(mixed1.path); }));
	CHECK(abortsInChild([&]() { readForcedVariantPhases(mixed2.path); }));
	CHECK(abortsInChild([&]() { readForcedVariantPhases(invalidPhase.path); }));
	CHECK(abortsInChild([&]() { readForcedVariantPhases(missingPhase.path); }));
}

TEST_CASE("parseTagsGff3")
{
	auto tags = parseTagsGff3("ID=gene1;gene_id=ENSG0001;gene_name=XIST;gene_type=lncRNA");
	CHECK(tags.size() == 4);
	CHECK(tags.at("ID") == "gene1");
	CHECK(tags.at("gene_id") == "ENSG0001");
	CHECK(tags.at("gene_name") == "XIST");
	CHECK(tags.at("gene_type") == "lncRNA");
	SUBCASE("malformed tags are skipped")
	{
		auto tags2 = parseTagsGff3("a=1;noequals;b=2=3;;c=");
		CHECK(tags2.size() == 2);
		CHECK(tags2.at("a") == "1");
		CHECK(tags2.at("c") == "");
	}
	SUBCASE("empty")
	{
		CHECK(parseTagsGff3("").empty());
	}
}

TEST_CASE("getGeneInfoFromStream")
{
	std::string gff = "##gff-version 3\n";
	gff += gffLine("chrX", "gene", 500, 900, "ID=g2;gene_id=G2;gene_type=protein_coding;gene_name=GENE2");
	gff += gffLine("chrX", "transcript", 500, 900, "ID=t2;gene_id=G2;gene_type=protein_coding;gene_name=GENE2");
	gff += gffLine("chrX", "gene", 100, 400, "ID=g1;gene_id=G1;gene_type=protein_coding;gene_name=GENE1");
	gff += gffLine("chrX", "gene", 100, 300, "ID=g3;gene_id=G3;gene_type=lncRNA;gene_name=LNC3");
	gff += gffLine("chr1", "gene", 100, 300, "ID=g4;gene_id=G4;gene_type=protein_coding;gene_name=AUTO4");
	gff += gffLine("X", "gene", 1000, 2000, "ID=g5;gene_id=G5;gene_type=protein_coding;gene_name=GENE5");
	gff += gffLine("23", "gene", 3000, 4000, "ID=g6;gene_id=G6;gene_type=protein_coding;gene_name=GENE6");
	gff += "#" + gffLine("chrX", "gene", 1, 2, "ID=g7;gene_id=G7;gene_type=protein_coding;gene_name=COMMENTED");
	gff += "short\n";
	SUBCASE("only protein coding")
	{
		std::istringstream stream { gff };
		auto genes = getGeneInfoFromStream(stream, true);
		REQUIRE(genes.size() == 4);
		CHECK(genes[0] == std::make_tuple<size_t, size_t, std::string, std::string>(100, 400, "G1", "GENE1"));
		CHECK(genes[1] == std::make_tuple<size_t, size_t, std::string, std::string>(500, 900, "G2", "GENE2"));
		CHECK(genes[2] == std::make_tuple<size_t, size_t, std::string, std::string>(1000, 2000, "G5", "GENE5"));
		CHECK(genes[3] == std::make_tuple<size_t, size_t, std::string, std::string>(3000, 4000, "G6", "GENE6"));
	}
	SUBCASE("all genes")
	{
		std::istringstream stream { gff };
		auto genes = getGeneInfoFromStream(stream, false);
		REQUIRE(genes.size() == 5);
		CHECK(genes[0] == std::make_tuple<size_t, size_t, std::string, std::string>(100, 300, "G3", "LNC3"));
		CHECK(genes[1] == std::make_tuple<size_t, size_t, std::string, std::string>(100, 400, "G1", "GENE1"));
	}
	SUBCASE("gene_type missing is excluded when only protein coding")
	{
		std::istringstream stream { gffLine("chrX", "gene", 1, 2, "ID=g;gene_id=G;gene_name=N") };
		CHECK(getGeneInfoFromStream(stream, true).empty());
	}
	SUBCASE("empty stream")
	{
		std::istringstream stream { "" };
		CHECK(getGeneInfoFromStream(stream, true).empty());
	}
}

TEST_CASE("getGeneInfoFromStream handles a missing gene_name")
{
	std::istringstream stream1 { gffLine("chrX", "gene", 1, 20, "ID=g;gene_id=G;gene_type=protein_coding") };
	std::vector<std::tuple<size_t, size_t, std::string, std::string>> genes;
	CHECK_NOTHROW(genes = getGeneInfoFromStream(stream1, true));
	REQUIRE(genes.size() == 1);
	CHECK(std::get<3>(genes[0]) == "G");
	std::istringstream stream2 { gffLine("chrX", "gene", 1, 20, "ID=g;gene_name=NAME;gene_type=protein_coding") };
	genes = getGeneInfoFromStream(stream2, true);
	REQUIRE(genes.size() == 1);
	CHECK(std::get<3>(genes[0]) == "NAME");
}

TEST_CASE("getGeneInfo reads .gff3 and .gff3.gz")
{
	std::string gff = "##gff-version 3\n";
	gff += gffLine("chrX", "gene", 500, 900, "ID=g2;gene_id=G2;gene_type=protein_coding;gene_name=GENE2");
	gff += gffLine("chrX", "gene", 100, 400, "ID=g1;gene_id=G1;gene_type=protein_coding;gene_name=GENE1");
	TempFile plain { ".gff3", gff };
	TempFile gz { ".gff3.gz" };
	{
		zstr::ofstream out { gz.path };
		out << gff;
	}
	auto genesPlain = getGeneInfo(plain.path, true);
	auto genesGz = getGeneInfo(gz.path, true);
	REQUIRE(genesPlain.size() == 2);
	CHECK(genesPlain == genesGz);
	CHECK(std::get<2>(genesPlain[0]) == "G1");
	SUBCASE("unknown extension aborts")
	{
		TempFile wrong { ".gtf", gff };
		CHECK(abortsInChild([&]() { getGeneInfo(wrong.path, true); }));
	}
}

TEST_CASE("readCellGrouping")
{
	SUBCASE("barcode header")
	{
		TempFile file { ".tsv", "barcode\tcelltype\nAAA\tTcell\nCCC\tBcell\nGGG\tTcell\n" };
		auto grouping = readCellGrouping(file.path);
		REQUIRE(grouping.size() == 3);
		CHECK(grouping.at("AAA") == "Tcell");
		CHECK(grouping.at("CCC") == "Bcell");
		CHECK(grouping.at("GGG") == "Tcell");
	}
	SUBCASE("cell header, case insensitive")
	{
		TempFile file { ".tsv", "Cell\tgroup\nAAA\tTcell\n" };
		auto grouping = readCellGrouping(file.path);
		CHECK(grouping.size() == 1);
		TempFile file2 { ".tsv", "BARCODE\tgroup\nAAA\tTcell\n" };
		CHECK(readCellGrouping(file2.path).size() == 1);
	}
	SUBCASE("invalid header aborts")
	{
		TempFile wrongName { ".tsv", "id\tgroup\nAAA\tTcell\n" };
		TempFile wrongColumns { ".tsv", "barcode\tgroup\textra\nAAA\tTcell\tx\n" };
		CHECK(abortsInChild([&]() { readCellGrouping(wrongName.path); }));
		CHECK(abortsInChild([&]() { readCellGrouping(wrongColumns.path); }));
		CHECK(abortsInChild([&]() { readCellGrouping("/nonexistent/file.tsv"); }));
	}
	SUBCASE("duplicate barcodes abort")
	{
		TempFile file { ".tsv", "barcode\tgroup\nAAA\tTcell\nAAA\tBcell\n" };
		CHECK(abortsInChild([&]() { readCellGrouping(file.path); }));
	}
}

TEST_CASE("readCellGrouping reads the last line even without a trailing newline")
{
	TempFile file { ".tsv", "barcode\tgroup\nAAA\tTcell\nCCC\tBcell" };
	CHECK(readCellGrouping(file.path).size() == 2);
}

TEST_CASE("readForcedVariantPhases reads the last line even without a trailing newline")
{
	TempFile file { ".tsv", "X:100:A:G\tmat\nX:200:C:T\tpat" };
	CHECK(readForcedVariantPhases(file.path).first.size() == 2);
}

TEST_CASE("readScReadCountsMatchCounts reads the last line even without a trailing newline")
{
	std::string contents = scReadCountsHeader + scReadCountsLine("X", 100, 'A', 'G', "cell1", 2, 3) + scReadCountsLine("X", 200, 'A', 'G', "cell1", 2, 3);
	contents.pop_back();
	TempFile file { ".tsv", contents };
	CHECK(readScReadCountsMatchCounts(file.path).size() == 4);
}

TEST_CASE("getGeneInfoFromStream reads the last line even without a trailing newline")
{
	std::string gff = gffLine("chrX", "gene", 1, 20, "ID=g;gene_id=G1;gene_type=protein_coding;gene_name=N1") + gffLine("chrX", "gene", 30, 40, "ID=g;gene_id=G2;gene_type=protein_coding;gene_name=N2");
	gff.pop_back();
	std::istringstream stream { gff };
	CHECK(getGeneInfoFromStream(stream, true).size() == 2);
}

TEST_CASE("readCellGrouping handles Windows line endings" * doctest::should_fail())
{
	// BUG (minor): no reader strips '\r', so a cell group file saved with CRLF line endings
	// gives group names like "Tcell\r", which then end up inside the output tables.
	TempFile file { ".tsv", "barcode\tgroup\r\nAAA\tTcell\r\n" };
	auto grouping = readCellGrouping(file.path);
	REQUIRE(grouping.count("AAA") == 1);
	CHECK(grouping.at("AAA") == "Tcell");
}

TEST_CASE("getGeneInfoFromStream handles Windows line endings" * doctest::should_fail())
{
	// BUG (minor): the last GFF3 attribute keeps the '\r', eg gene_name "GENE1\r"
	std::istringstream stream { "chrX\tHAVANA\tgene\t1\t20\t.\t+\t.\tgene_id=G1;gene_type=protein_coding;gene_name=GENE1\r\n" };
	auto genes = getGeneInfoFromStream(stream, true);
	REQUIRE(genes.size() == 1);
	CHECK(std::get<3>(genes[0]) == "GENE1");
}

TEST_CASE("readBarcodeWhitelist")
{
	TempFile file { ".txt", "AAA-1\nCCC-1\n\n  GGG-1  \nAAA-1\nTTT-1" };
	auto whitelist = readBarcodeWhitelist(file.path);
	CHECK(whitelist.size() == 4);
	CHECK(whitelist.count("AAA-1") == 1);
	CHECK(whitelist.count("CCC-1") == 1);
	CHECK(whitelist.count("GGG-1") == 1);
	CHECK(whitelist.count("TTT-1") == 1);
	CHECK(readBarcodeWhitelist("/nonexistent/file").empty());
}

TEST_CASE("writeCellGroupStatistics")
{
	std::vector<CellMatch> matches {
		makeMatch("c1", "X:1:A:G", false, 10),
		makeMatch("c1", "X:2:A:G", true, 10),
		makeMatch("c2", "X:1:A:G", false, 20),
		makeMatch("c3", "X:1:A:G", false, 40),
		makeMatch("c4", "X:1:A:G", false, 20),
	};
	std::unordered_map<std::string, std::string> grouping { { "c1", "B" }, { "c2", "A" }, { "c3", "B" }, { "unused", "C" } };
	TempFile file { ".tsv" };
	writeCellGroupStatistics(matches, grouping, file.path);
	auto lines = file.readLines();
	REQUIRE(lines.size() == 5);
	CHECK(lines[0] == "Group\tASE\tASE_cells\tASE_percent\tASE_cells_percent");
	CHECK(lines[1] == "A\t20\t1\t20\t25");
	CHECK(lines[2] == "B\t60\t2\t60\t50");
	CHECK(lines[3] == "No group\t20\t1\t20\t25");
	CHECK(lines[4] == "Total\t100\t4\t100\t100");
}

TEST_CASE("writeResultVariants, writeResultCells and writeResultCellOnlyNonescapeVariants")
{
	OutputBuilder builder;
	builder.matches = {
		makeMatch("cellB", "X:100:A:G", false, 3),
		makeMatch("cellA", "X:20:C:T", true, 2),
		makeMatch("cellA", "X:100:A:G", true, 1),
	};
	builder.variantIsMatRef = { { "X:100:A:G", true }, { "X:20:C:T", false } };
	builder.variantConfidence = { { "X:100:A:G", 2.5 }, { "X:20:C:T", 0.5 } };
	builder.variantEscape = { { "X:100:A:G", 0.25 } };
	builder.cellIsMatActive = { { "cellA", false }, { "cellB", true } };
	builder.cellConfidence = { { "cellA", 1.5 }, { "cellB", 3 } };
	builder.cellConfidenceOnlyNonescape = { { "cellA", 1 }, { "cellB", 4 } };
	builder.cellEscape = { { "cellB", 0.5 } };
	EMOutput output = builder.build();
	SUBCASE("variants, mat/pat")
	{
		std::ostringstream stream;
		writeResultVariants(output, true, stream);
		auto lines = TestHelpers::lines(stream.str());
		REQUIRE(lines.size() == 3);
		CHECK(lines[0] == "variant\tcoverage\tphase\tphase_confidence\tescape_estimate\tescape_ci_low\tescape_ci_high");
		CHECK(lines[1] == "X:20:C:T\t2\tpat\t0.5\t0.001\t0\t0");
		CHECK(lines[2] == "X:100:A:G\t4\tmat\t2.5\t0.25\t0\t0");
	}
	SUBCASE("variants, hap1/hap2")
	{
		std::ostringstream stream;
		writeResultVariants(output, false, stream);
		auto lines = TestHelpers::lines(stream.str());
		REQUIRE(lines.size() == 3);
		CHECK(lines[1] == "X:20:C:T\t2\thap2\t0.5\t0.001\t0\t0");
		CHECK(lines[2] == "X:100:A:G\t4\thap1\t2.5\t0.25\t0\t0");
	}
	SUBCASE("cells")
	{
		std::ostringstream stream;
		writeResultCells(output, true, stream);
		auto lines = TestHelpers::lines(stream.str());
		REQUIRE(lines.size() == 3);
		CHECK(lines[0] == "cell\tcoverage\tactive_chrX\tactive_chrX_confidence\tescape_estimate\tescape_ci_low\tescape_ci_high");
		CHECK(lines[1] == "cellA\t3\tpat\t1.5\t0.001\t0\t0");
		CHECK(lines[2] == "cellB\t3\tmat\t3\t0.5\t0\t0");
	}
	SUBCASE("cells only nonescape")
	{
		std::ostringstream stream;
		writeResultCellOnlyNonescapeVariants(output, false, stream);
		auto lines = TestHelpers::lines(stream.str());
		REQUIRE(lines.size() == 3);
		CHECK(lines[1] == "cellA\t3\thap2\t1\t0.001\t0\t0");
		CHECK(lines[2] == "cellB\t3\thap1\t4\t0.5\t0\t0");
	}
}

TEST_CASE("writePseudobulkResults")
{
	std::vector<PseudobulkInfo> pseudobulk;
	pseudobulk.push_back(PseudobulkInfo { 0, "X:1:A:G", 9, 1, 0, 0 });
	pseudobulk.push_back(PseudobulkInfo { 1, "X:2:A:G", 0, 0, 0, 0 });
	pseudobulk.push_back(PseudobulkInfo { 2, "X:3:A:G", 5, 0, 3, 2 });
	TempFile file { ".tsv" };
	writePseudobulkResults(pseudobulk, true, "variant", file.path);
	auto lines = file.readLines();
	REQUIRE(lines.size() == 4);
	CHECK(lines[0] == "variant\tmat_active_expression\tmat_inactive_expression\tpat_active_expression\tpat_inactive_expression\tXi\tpvalue_Xiover10_unadjusted\tpvalue_Xiunder10_unadjusted");
	auto parts = split(lines[1], '\t');
	REQUIRE(parts.size() == 8);
	CHECK(parts[0] == "X:1:A:G");
	CHECK(parts[1] == "9");
	CHECK(parts[2] == "1");
	CHECK(parts[3] == "0");
	CHECK(parts[4] == "0");
	CHECK(std::stod(parts[5]) == doctest::Approx(0.1));
	CHECK(std::stod(parts[6]) == doctest::Approx(getBinomialPValueGreaterThan(0.1, 1, 10)));
	CHECK(std::stod(parts[7]) == doctest::Approx(getBinomialPValueLessThan(0.1, 1, 10)));
	parts = split(lines[2], '\t');
	REQUIRE(parts.size() == 8);
	CHECK(parts[5] == "0");
	CHECK(std::stod(parts[6]) == doctest::Approx(1));
	CHECK(std::stod(parts[7]) == doctest::Approx(1));
	parts = split(lines[3], '\t');
	REQUIRE(parts.size() == 8);
	CHECK(std::stod(parts[5]) == doctest::Approx(0.2));
	SUBCASE("hap names and multi-column first column")
	{
		TempFile file2 { ".tsv" };
		writePseudobulkResults({ PseudobulkInfo { 0, "G1\tGENE1", 1, 1, 1, 1 } }, false, "gene_id\tgene_name", file2.path);
		auto lines2 = file2.readLines();
		REQUIRE(lines2.size() == 2);
		CHECK(lines2[0].substr(0, 50) == std::string("gene_id\tgene_name\thap1_active_expression\thap1_inactive_expression").substr(0, 50));
		CHECK(split(lines2[0], '\t').size() == 9);
		CHECK(split(lines2[1], '\t').size() == 9);
		CHECK(split(lines2[1], '\t')[6] == "0.5");
	}
}

TEST_CASE("writePseudobulkCellResults")
{
	OutputBuilder builder;
	builder.matches = {
		makeMatch("cellA", "X:1:A:G", false, 1),
		makeMatch("cellB", "X:1:A:G", false, 1),
	};
	builder.cellIsMatActive = { { "cellA", true }, { "cellB", false } };
	EMOutput output = builder.build();
	std::vector<PseudobulkInfo> pseudobulk {
		PseudobulkInfo { 0, "cellA", 3, 1, 0, 0 },
		PseudobulkInfo { 1, "cellB", 0, 0, 0, 0 },
	};
	TempFile file { ".tsv" };
	writePseudobulkCellResults(pseudobulk, output, true, { { "cellA", "Tcell" } }, file.path);
	auto lines = file.readLines();
	REQUIRE(lines.size() == 3);
	CHECK(lines[0] == "cell\tgroup\tactive_chrX\tactive_expression\tinactive_expression\tXi");
	CHECK(lines[1] == "cellA\tTcell\tmat\t3\t1\t0.25");
	CHECK(lines[2] == "cellB\t-\tpat\t0\t0\t0");
}

TEST_CASE("writeGenesPerVariant")
{
	OutputBuilder builder;
	builder.matches = {
		makeMatch("c", "X:300:A:G", false, 1),
		makeMatch("c", "X:100:A:G", false, 1),
		makeMatch("c", "X:200:A:G", false, 1),
	};
	EMOutput output = builder.build();
	std::vector<std::vector<std::pair<std::string, std::string>>> genes;
	genes.resize(3);
	genes[output.helpers.variantNameToIndex.at("X:100:A:G")] = { { "G1", "GENE1" }, { "G2", "GENE2" } };
	genes[output.helpers.variantNameToIndex.at("X:300:A:G")] = { { "G3", "GENE3" } };
	TempFile file { ".tsv" };
	writeGenesPerVariant(output, genes, file.path);
	auto lines = file.readLines();
	REQUIRE(lines.size() == 4);
	CHECK(lines[0] == "variant\tcount_genes\tgene_ids\tgene_names");
	CHECK(lines[1] == "X:100:A:G\t2\tG1,G2\tGENE1,GENE2");
	CHECK(lines[2] == "X:200:A:G\t0\t\t");
	CHECK(lines[3] == "X:300:A:G\t1\tG3\tGENE3");
}

TEST_CASE("writeResultsSummary")
{
	OutputBuilder builder;
	builder.matches = {
		makeMatch("c1", "X:10:A:G", false, 9),
		makeMatch("c1", "X:10:A:G", true, 1),
		makeMatch("c1", "X:20:C:T", true, 4),
		makeMatch("c2", "X:10:A:G", true, 5),
		makeMatch("c2", "X:30:G:A", false, 3),
		makeMatch("c3", "X:20:C:T", false, 7),
	};
	builder.variantIsMatRef = { { "X:10:A:G", true }, { "X:20:C:T", false }, { "X:30:G:A", true } };
	builder.variantConfidence = { { "X:10:A:G", 5 }, { "X:20:C:T", 2 }, { "X:30:G:A", 1.99 } };
	builder.cellIsMatActive = { { "c1", true }, { "c2", false }, { "c3", true } };
	builder.cellConfidence = { { "c1", 5 }, { "c2", 2 }, { "c3", 0 } };
	EMOutput output = builder.build();
	TempFile file { ".txt" };
	writeResultsSummary(file.path, true, output);
	auto lines = file.readLines();
	REQUIRE(lines.size() == 9);
	CHECK(lines[0] == "Total allele specific expression: 29");
	CHECK(lines[1] == "Total heterozygously expressed variants: 3");
	CHECK(lines[2] == "Confidently phased variants: 2");
	CHECK(lines[3] == "Total cells with allele specific expression: 3");
	CHECK(lines[4] == "Confidently assigned cells: 2");
	CHECK(lines[5] == "Cells confidently assigned to mat: 1");
	CHECK(lines[6] == "Cells confidently assigned to pat: 1");
	CHECK(lines[7] == "Skew: 50/50");
	// c1: X:10 ref 9 active, alt 1 inactive; X:20 (pat ref) alt 4 active; c2: X:10 alt 5 active
	CHECK(lines[8] == "Average inactive expression over chrX (Xi/(Xi+Xa)): 5.26%");
	SUBCASE("hap names and skew rounding")
	{
		OutputBuilder builder2 = builder;
		builder2.cellConfidence = { { "c1", 5 }, { "c2", 2 }, { "c3", 2 } };
		TempFile file2 { ".txt" };
		writeResultsSummary(file2.path, false, builder2.build());
		auto lines2 = file2.readLines();
		REQUIRE(lines2.size() == 9);
		CHECK(lines2[5] == "Cells confidently assigned to hap1: 2");
		CHECK(lines2[6] == "Cells confidently assigned to hap2: 1");
		CHECK(lines2[7] == "Skew: 67/33");
	}
}

TEST_CASE("writeResultsSummary when no cells are confidently assigned")
{
	OutputBuilder builder;
	builder.matches = { makeMatch("c1", "X:10:A:G", false, 1), makeMatch("c1", "X:10:A:G", true, 1) };
	TempFile file { ".txt" };
	writeResultsSummary(file.path, true, builder.build());
	auto lines = file.readLines();
	REQUIRE(lines.size() == 9);
	CHECK(lines[7].find("Skew: ") == 0);
	auto skew = split(lines[7].substr(6), '/');
	REQUIRE(skew.size() == 2);
	CHECK(std::stoi(skew[0]) >= 0);
	CHECK(std::stoi(skew[0]) <= 100);
	CHECK(lines[8].find("nan") == std::string::npos);
}

}
