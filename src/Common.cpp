#include <algorithm>
#include <limits>
#include <regex>
#include "Common.h"

std::string matHapName(const bool phasesAreMatPat)
{
	return phasesAreMatPat ? "mat" : "hap1";
}

std::string patHapName(const bool phasesAreMatPat)
{
	return phasesAreMatPat ? "pat" : "hap2";
}

std::vector<std::string> split(const std::string& raw, const char separator)
{
	std::vector<std::string> result;
	size_t lastSplit = 0;
	for (size_t i = 0; i < raw.size(); i++)
	{
		if (raw[i] != separator) continue;
		result.emplace_back(raw.begin()+lastSplit, raw.begin()+i);
		lastSplit = i+1;
	}
	result.emplace_back(raw.begin()+lastSplit, raw.end());
	return result;
}

std::string lowercase(std::string raw)
{
	for (size_t i = 0; i < raw.size(); i++)
	{
		raw[i] = std::tolower(raw[i]);
	}
	return raw;
}

std::vector<std::string> getVariantOrder(const std::unordered_map<std::string, size_t>& nameToIndex)
{
	std::vector<std::string> result;
	for (const auto& pair : nameToIndex)
	{
		result.emplace_back(pair.first);
	}
	std::sort(result.begin(), result.end(), [](const std::string& left, const std::string& right)
	{
		size_t leftPos = 0;
		for (size_t i = 0; i < left.size(); i++)
		{
			if (left[i] != ':') continue;
			for (size_t j = i+1; j < left.size(); j++)
			{
				if (left[j] != ':') continue;
				leftPos = std::stoull(left.substr(i+1, j-i-1));
				break;
			}
			break;
		}
		size_t rightPos = 0;
		for (size_t i = 0; i < right.size(); i++)
		{
			if (right[i] != ':') continue;
			for (size_t j = i+1; j < right.size(); j++)
			{
				if (right[j] != ':') continue;
				rightPos = std::stoull(right.substr(i+1, j-i-1));
				break;
			}
			break;
		}
		return leftPos < rightPos;
	});
	return result;
}

std::vector<std::string> getCellOrder(const std::unordered_map<std::string, size_t>& nameToIndex)
{
	std::vector<std::string> result;
	for (const auto& pair : nameToIndex)
	{
		result.emplace_back(pair.first);
	}
	std::sort(result.begin(), result.end());
	return result;
}

size_t parseVariantPosition(const std::string& name)
{
	size_t firstColon = name.size();
	for (size_t i = 0; i < name.size(); i++)
	{
		if (name[i] == ':')
		{
			if (firstColon == name.size())
			{
				firstColon = i;
			}
			else
			{
				return std::stoull(name.substr(firstColon+1, i-(firstColon+1)));
			}
		}
	}
	return 0;
}

std::tuple<std::string, size_t, size_t> parseBedRegion(const std::string& region)
{
	std::regex regex { "^([chrCHR0123456789XYMxymt])+:([0-9]+)-([0-9]+)$" };
	std::smatch matches;
	if (std::regex_match(region, matches, regex))
	{
		std::string chromosome = matches[1];
		size_t start = std::stoull(matches[2]);
		size_t end = std::stoull(matches[3]);
		if (end > start)
		{
			return std::make_tuple(chromosome, start, end);
		}
	}
	return std::make_tuple<std::string, size_t, size_t>("", std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
}
